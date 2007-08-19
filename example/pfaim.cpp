//include variables set by configure
#include "parfe_ConfigDefs.h"
#include <string>
#include <fstream>
#include <vector>
#include <Epetra_Time.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <AztecOO_StatusTestResNorm.h>
#include <time.h>
#include "DistMesh.h"
#include <AztecOO.h>
//IO
#include "MEDIT_MeshWriter.hpp"
#include "PMVIS_MeshWriter.hpp"
#include "MEDIT_DispWriter.hpp"
#include "IBT_ProblemReader.hpp"
#include "IBT_ProblemWriter.hpp"
#include "MATLAB_MatrixWriter.hpp"
#include "MATLAB_VectorWriter.hpp"
#include "IBT_DispWriter.hpp"
#include "IBT_SolutionWriter.hpp"
#include "EbeElasticityProblem.h"
#include "VbrElasticityProblem.h"
#include "Jacobi.h"
#include <unistd.h>
#include <getopt.h>
//header file for memory information
#if defined(HAVE_MALLINFO)
#include <malloc.h>
#elif defined(HAVE_HEAP_INFO)
#include <catamount/catmalloc.h>
#endif
#include "FEParameterList.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>
#ifdef HAVE_MLLIB
#define HAVE_ML_MATRIX_FREE //comment this line if don't have ml_MatrixFree
#include <ml_MultiLevelPreconditioner.h>
#ifdef HAVE_ML_MATRIX_FREE
#include <ml_MatrixFreePreconditioner.h>
#endif
#endif
#ifdef HAVE_IFPACKLIB
#include <Ifpack_IC.h>
#include <Ifpack_AdditiveSchwarz.h>
#endif

using namespace ML_Epetra;
using namespace Teuchos;

//Helper function to compute memory usage
unsigned meminfo()
{
#if defined(HAVE_MALLINFO)
    long int fragments, total_free, largest_free, total_used;
    /* int percent; */

    /* use system call to get memory used information */

    struct mallinfo M = mallinfo();
    fragments = M.ordblks + M.smblks + M.hblks;
    total_free = M.fsmblks + M.fordblks;
    total_used = M.hblkhd + M.usmblks + M.uordblks;
    /*  total_free = ml_total_mem - total_used;  */
    largest_free = -1024;

    /* convert to Mbytes */

    return( (unsigned)(total_used/(1024*1024)) );

#elif defined(HAVE_HEAP_INFO)
    size_t fragments;
    unsigned long total_free;
    unsigned long largest_free;
    unsigned long total_used;

    heap_info(&fragments, &total_free, &largest_free, &total_used);
   
    return( (unsigned)(total_used/(1024*1024)) );

#else
    return(0);
#endif
}

static const char* help_message = "Usage: pfaim.exe [options] filename\n" \
"  -h, --help                This small usage guide\n" \
"  -i, --maxiters=n          Set the maximum number of solver iterations\n" \
"  -t, --tolerance=r         Set the solver tolerance\n" \
"  -p, --precond=name        Specify the preconditioner to be used by the\n" \
"                            solver (One of: ml, matrixfree,\n"
"                            ic[k] (with optional k = fill-in), Jacobi)\n" \
"      --memory=quantity     Specify memory consumption (One of: low/normal/high)\n" \
"  -a, --ascii               Read/Write ASCII encoded files, instead of HDF5\n" \
"      --norepart            Don't apply load balancing\n" \
"      --nosolve             Do not solve the problem\n" \
"      --printmedit          Print mesh and displacements to be viewed with MEDIT\n" \
"                            http://www.ann.jussieu.fr/~frey/logiciels/medit.html\n" \
"      --printpmvis          Print mesh partition to be viewed with PMVIS\n" \
"                            http://www-users.cs.umn.edu/~oztekin/pmvis/\n" \
"      --printmatrix         Print out the matrix\n" \
"      --printrhs            Print out the right hand side vector\n" \
"      --printlhs            Print out the left hand side vector\n" \
"      --toASCII             Convert to ascii format and exit\n" \
"      --toHDF5              Convert to hdf5 format and exit\n" \
"  -v, --verbose             Get more output.\n";

#ifdef HAVE_MPI
#define EXIT(exit_code) \
    MPI_Finalize(); \
    exit(exit_code);
#elif
#define EXIT(exit_code) exit(exit_code);
#endif

int main( int argc, char* argv[] )
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    // define an Epetra communicator
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
  
    // get the proc ID of this process
    int MyPID = comm.MyPID();
  
    // get the total number of processes
    int NumProc = comm.NumProc();

    Epetra_Time timer(comm);
    double start_time = timer.WallTime();
    double temp_time = start_time;

//Macros to print information
#define INFO(mesg) \
    if (MyPID == 0) cout << "*** INFO ***: " << mesg << endl;

#define WARN(mesg) \
    if (MyPID == 0) cout << "*** WARNING ***: " << mesg << endl;

#define ERROR(mesg, exit_code)\
    if (MyPID == 0) cerr << "*** ERROR ***: " << mesg << endl; \
    EXIT(exit_code);		


    FEParameterList fe_param(NumProc);

    //Command line parsing
    static struct option long_options[] = {
	{"ebe", 0, 0, 'e'},
	{"maxiters", 1, 0, 'i'},
	{"tolerance", 1, 0, 't'},
	{"precond", 1, 0, 'p'},
	{"help", 0, 0, 'h'},
	{"norepart",0, 0, 'N'},
	{"ascii",0,0,'a'},
	{"nosolve",0,0,'S'},
	{"printmedit",0,0,'M'},
	{"printpmvis",0,0,'P'},
	{"printmatrix",0,0,'A'},
	{"printrhs",0,0,'B'},
	{"printlhs",0,0,'X'},
	{"toASCII", 0,0,'C'},
	{"toHDF5", 0,0,'H'},
	{"verbose",0,0,'v'},
	{"memory",1,0, 'U'},
	{0, 0, 0, 0}
    };

    std::string tmp;
    int option_index = 0;
    char optchar;
    while ((optchar = getopt_long(argc, argv, "hi:t:p:av", long_options, &option_index)) != -1) {
	switch (optchar)
	{
	    case 'i':
		fe_param.set("iteration limit", atoi(optarg));
		break;
	    case 't':
		fe_param.set("tolerance", atof(optarg));
		break;
	    case 'p':
		tmp = optarg;
		if (tmp.find( "ic", 0 ) == 0) {
		    tmp.resize(2);
		    fe_param.set("preconditioner", tmp);
		    fe_param.set("ic: fill-in", atoi(optarg+2));
		} else {
		    fe_param.set("preconditioner", std::string(optarg));
		}
		break;
	    case 'a':
		fe_param.set("hdf5 I/O", false);
		break;
	    case 'v':
		fe_param.set("verbose", true);
		break;
	    case 'N':
		fe_param.set("load balancing", false);
		break;
	    case 'S':
		fe_param.set("solve", false);
		break;
	    case 'M':
		fe_param.set("print for MEDIT", true);
		break;
	    case 'P':
		fe_param.set("print for PMVIS", true);
		break;
	    case 'A':
		fe_param.set("print matrix", true);
		break;
	    case 'X':
		fe_param.set("print lhs", true);
		break;
	    case 'B':
		fe_param.set("print rhs", true);
		break;
	    case 'C':
		fe_param.set("convert to ascii", true);
		fe_param.set("load balancing", false);
		break;
	    case 'H':
		fe_param.set("convert to hdf5", true);
		fe_param.set("hdf5 I/O", false);
		fe_param.set("load balancing", false);
		break;
	    case 'U':
		fe_param.set("memory usage", std::string(optarg));
		break;
	    case 'h':
		ERROR(help_message,0);
	    default:
		ERROR(help_message,1);
	}
    }
  
    Epetra_MultiVector* BCs = 0;

    if (optind != argc-1) {
	ERROR("missing file argument" << endl << help_message, 1)
    }
  
    std::string ifile_base(argv[argc-1]);
    int ext_pos =  ifile_base.rfind(".mesh");
    ifile_base = ifile_base.substr(0, ext_pos);
  
    temp_time = timer.WallTime();
  
    bool hdf5 = fe_param.get<bool>("hdf5 I/O");
    ProblemReader* pr = 0;
    try {
	if (hdf5) {
	    fe_param.set("input", ifile_base + ".mesh.h5");
	    pr = new IBT_ProblemReader<HDF5_GReader>(fe_param.get<std::string>("input"), comm);
	} else {
	    fe_param.set("input", ifile_base + ".mesh");
	    pr = new IBT_ProblemReader<C_ASCII_GReader>(fe_param.get<std::string>("input"), comm);
	}
    } catch (std::string err) { 
	ERROR(err,1);
    }

    fe_param.Scan(pr);

    //If these parameters are still unset, provide some default values
    fe_param.get("iteration limit", 1000);
    fe_param.get("tolerance", 1e-6);

    //Check options
    if (MyPID == 0) fe_param.print(cout);


    DistMesh* mesh = new DistMesh(fe_param.get<int>("#nodes"), fe_param.get<int>("#dimensions"), fe_param.get<int>("#elements"), fe_param.get<int>("#nodes per element"), fe_param.get<IntMap>("material ids"), comm);
  
    mesh->Scan(pr);
    BoundaryCondition* bcond = new BoundaryCondition(fe_param.get<int>("#dofs per node"), comm);
    bcond->Scan(pr);
  
    //TODO: Fix destructor call
    if (hdf5)
	delete dynamic_cast<IBT_ProblemReader<HDF5_GReader>*>(pr);
    else 
	delete dynamic_cast<IBT_ProblemReader<C_ASCII_GReader>*>(pr);
    //--------------------------end of input--------------------------------
    INFO("Input Time: " << timer.WallTime() - temp_time);
    temp_time = timer.WallTime();
    //cerr << "Proc " << MyPID << ": timer::WallTime() called." << endl;
    mesh->Redistribute(fe_param.get<bool>("load balancing"), fe_param.get<bool>("verbose"));
    if (fe_param.get<bool>("load balancing")) {
	INFO("Time used to partition mesh: " << timer.WallTime() - temp_time);
    } else {
	WARN("No load balancing used");
    }
  
    //redistribute the boundary condition data according to the node map
    bcond->Redistribute(mesh->NodeMap()->NumMyElements(), mesh->NodeMap()->MyGlobalElements());

    if (fe_param.get<bool>("print for MEDIT")) {
	std::string meshfile = ifile_base + ".medit.mesh";
	MEDIT_MeshWriter* meshw = new MEDIT_MeshWriter(meshfile, fe_param.get<int>("#dimensions"), comm);
	mesh->Print(meshw);
	delete meshw;
    }

    if (fe_param.get<bool>("print for PMVIS")) {
	std::string con = ifile_base + ".con";
	std::string xyz = ifile_base + ".xyz";
	std::string part = ifile_base + ".part";
	PMVIS_MeshWriter* meshw = new PMVIS_MeshWriter(xyz, con, part, comm);
	mesh->Print(meshw);
	delete meshw;
    }

    int memused = meminfo();
    INFO("Mem used to hold mesh information: " << memused <<  "MB");
  
    INFO("Total Mesh setup time: " << timer.WallTime() - start_time);

    ProblemWriter* pw = 0;
    if (fe_param.get<bool>("convert to hdf5")) {
	pw = new IBT_ProblemWriter<HDF5_GWriter>(ifile_base+".mesh.h5", comm);
    } else if (fe_param.get<bool>("convert to ascii")) {
	pw = new IBT_ProblemWriter<C_ASCII_GWriter>(ifile_base+".mesh", comm);
    }

    if (pw) {
	temp_time = timer.WallTime();
	fe_param.Print(pw);
	mesh->Print(pw);
	bcond->Print(pw);
	delete pw;
       
	if (fe_param.get<bool>("convert to hdf5")) {
	    INFO("Time for hdf5 output: " << timer.WallTime() -temp_time);
	}
	else if (fe_param.get<bool>("convert to ascii")) {
	    INFO("Time for ascii output: " << timer.WallTime() -temp_time);
	}

#ifdef HAVE_MPI
	MPI_Finalize();
#endif
	exit( EXIT_SUCCESS );
    }

    //-------------------------beginning of assembly------------------------
    temp_time = timer.WallTime();
    if (fe_param.get<std::string>("preconditioner") == "Jacobi" ||
	fe_param.get<std::string>("preconditioner") == "matrixfree") {
	fe_param.set("element by element", true);
    }

    ElasticityProblem* ep;
    if (fe_param.get<bool>("element by element"))
      ep = new EbeElasticityProblem(*mesh, fe_param);
    else
      ep = new VbrElasticityProblem(*mesh, fe_param);

    ep->Impose(*bcond);

    BCs = new Epetra_MultiVector(*ep->GetRHS());

    INFO("Time used to assemble matrix: " << timer.WallTime() - temp_time);
  
    memused = meminfo() - memused;
    INFO("Mem used to hold stiffness matrix operator: " << memused << "MB");
  
    Epetra_Vector* X = dynamic_cast<Epetra_Vector*>(ep->GetLHS());
    Epetra_Vector* B = dynamic_cast<Epetra_Vector*>(ep->GetRHS());
    if (X == 0 || B == 0) {
	ERROR("Use Epetra_Vector in ElastictyProblem", 1);
    }

    if (fe_param.get<bool>("print matrix")) {
	MATLAB_MatrixWriter<HDF5_GWriter>* mw = new MATLAB_MatrixWriter<HDF5_GWriter>(ifile_base + ".matrix.h5", comm);
	ep->PrintMatrix(mw);
	delete mw;
    }
  
    if (fe_param.get<bool>("print lhs")) {
	MATLAB_VectorWriter<HDF5_GWriter>* vw = new MATLAB_VectorWriter<HDF5_GWriter>(ifile_base + ".lhs.h5", comm);
	ep->PrintLHS(vw);
	delete vw;
    }

    if (fe_param.get<bool>("print rhs")) {
	MATLAB_VectorWriter<HDF5_GWriter>* vw = new MATLAB_VectorWriter<HDF5_GWriter>(ifile_base + ".rhs.h5", comm);
	ep->PrintRHS(vw);
	delete vw;
    }

    vector<double> null_space;

    if  (fe_param.get<bool>("solve")) {
	//Compute null space for ml preconditioner
	if (fe_param.get<std::string>("preconditioner") == "ml" ||
	    fe_param.get<std::string>("preconditioner") == "matrixfree") {
	    // MS // the following seems to be wrong for ebe

	    int num_rows = 3 * mesh->NodeMap()->NumMyPoints();
	    null_space.resize(6 * num_rows);
      
	    double* row_it = &null_space[0];
	    //clear vector
	    for (double* it = &null_space[0];it < &null_space[0]+6*num_rows; ++it)
		*it = 0;
	    //first coordinate is constant
	    for (double* it = row_it;it < row_it+num_rows; it+=3)
		*it = 1;
	    //second coordinate is constant
	    row_it += num_rows;
	    for (double* it = row_it+1;it < row_it+num_rows; it+=3)
		*it = 1;
	    //third coordinate is constant
	    row_it += num_rows;
	    for (double* it = row_it+2;it < row_it+num_rows; it+=3)
		*it = 1;
	    // u = y and v = -x
	    double* values = mesh->Coordinates()->Values();
	    row_it += num_rows;
	    for (double* it = row_it;it < row_it+num_rows; it+=3) {
		*it = *(values+1);
		*(it+1) = -*values;
		values+=3;
	    }
	    // v = z and w = -y
	    values = mesh->Coordinates()->Values();
	    row_it += num_rows;
	    for (double* it = row_it;it < row_it+num_rows; it+=3) {
		*(it+1) = *(values+2);
		*(it+2) = -*(values+1);
		values+=3;
	    }
	    // u = z and w = -x
	    values = mesh->Coordinates()->Values();
	    row_it += num_rows;
	    for (double* it = row_it;it < row_it+num_rows; it+=3) {
		*it = *(values+2);
		*(it+2) = -*values;
		values+=3;
	    }
	}
    
	RefCountPtr<Epetra_CrsGraph> graph;

	if (fe_param.get<std::string>("preconditioner") == "matrixfree") 
	{
	    ElementByElementMatrix* ebe = dynamic_cast<ElementByElementMatrix*>(ep->GetOperator());
	    if (ebe != 0)
	    {
		graph = rcp(new Epetra_CrsGraph(Copy, ebe->Map(), 0));
		mesh->MatrixGraph(*(graph.get()));
	    }
	}
	//free unused data
	//delete mesh;
	//delete bcond;
    
	//------------------------beginning of solution------------------------ 

	RefCountPtr<Epetra_Operator> Prec;
	RefCountPtr<Epetra_Vector> diagonal;

	AztecOO solver(*ep);

	solver.SetAztecOption(AZ_solver, AZ_cg);
	solver.SetAztecOption(AZ_output, 16);
	//AztecOO_StatusTestResNorm restest(*ep->Operator(), *X, *B, fe_param.get<double>("tolerance"));
    
	temp_time = timer.WallTime();

	if (fe_param.get<std::string>("preconditioner") == "Jacobi") {
	    ElementByElementMatrix* ebe = dynamic_cast<ElementByElementMatrix*>(ep->GetOperator());

	    // impose the boundary conditions for nodes with fixed displacement.
	    // This is done by applying the operator without the fix for BC,
	    // then resetting the nodes corresponding to fixed nodes in RHS.
	    ebe->Apply_NoResetBoundaries(*B, *X);
	    *B = *X;
	    ebe->ResetBoundaries(B);
	    X->PutScalar(0.0);

	    //Use user defined jacobi preconditioner
	    diagonal = rcp(new Epetra_Vector(ebe->Map()));
	    ebe->ExtractDiagonalCopy(*diagonal);
	    Prec = rcp(new Jacobi(*diagonal));  
	    solver.SetPrecOperator(Prec.get());

	    // FIXME

	    for (int i = 0; i < B->MyLength(); ++i)
		if ((*B)[i] != 0.0)
		    (*X)[i] = (*B)[i];
	}
	else if (fe_param.get<std::string>("preconditioner") == "matrixfree") 
	{
	    ElementByElementMatrix* ebe = dynamic_cast<ElementByElementMatrix*>(ep->GetOperator());
	    
	    // impose the boundary conditions for nodes with fixed displacement.
	    // This is done by applying the operator without the fix for BC,
	    // then resetting the nodes corresponding to fixed nodes in RHS.
	    ebe->Apply_NoResetBoundaries(*B, *X);
	    *B = *X;
	
	    ebe->ResetBoundaries(B);
	    X->PutScalar(0.0);
	    
	    ParameterList MLList;
	    MLList.set("prec: type", "hybrid");
	    MLList.set("smoother: type", "Chebyshev");
	    MLList.set("smoother: degree", 3);
	    MLList.set("low memory", true);
	    MLList.set("aggregation: type", "METIS");
	    MLList.set("eigen-analysis: max iters", 20);
	    if (!fe_param.get<bool>("verbose"))
		MLList.set("output", 0);
	    MLList.sublist("ML list").set("max levels", 10);

	    SetDefaults("SA", MLList.sublist("ML list"));
	    MLList.sublist("ML list").set("aggregation: damping factor", 1.333);
	    MLList.sublist("ML list").set("coarse: max size", 1024);
	    MLList.sublist("ML list").set("smoother: type", "symmetric Gauss-Seidel");
	    MLList.sublist("ML list").set("aggregation: type", "Uncoupled-MIS");
	    MLList.sublist("ML list").set("coarse: type", "Amesos-KLU");
	    MLList.sublist("ML list").set("low memory usage", true);

	    if (fe_param.get<std::string>("memory usage") == "very-low")
	    {
		MLList.set("smoother: degree", 8);
		MLList.set("AP allocation factor", 0.3);
		MLList.set("aggregation: nodes per aggregate", 512);
		MLList.sublist("ML list").set("cycle applications", 10);
	    }
	    else if (fe_param.get<std::string>("memory usage") == "low")
	    {
		MLList.set("smoother: degree", 5);
		MLList.set("AP allocation factor", 0.5);
		MLList.set("aggregation: nodes per aggregate", 256);
		MLList.sublist("ML list").set("cycle applications", 5);
	    }
	    else
	    {
		MLList.set("smoother: degree", 3);
		MLList.set("AP allocation factor", 0.8);
		MLList.set("aggregation: nodes per aggregate", 64);
		MLList.sublist("ML list").set("cycle applications", 3);
	    }

	    Epetra_MultiVector NullSpace(Copy, graph->Map(),
					 &null_space[0], graph->Map().NumMyPoints(), 6);

	    diagonal = rcp(new Epetra_Vector(ebe->Map()));
	    ebe->ExtractDiagonalCopy(*diagonal);
#ifdef HAVE_ML_MATRIX_FREE
	    try
	    {
		Prec = rcp(new MatrixFreePreconditioner(*ebe, *graph, NullSpace,
							*diagonal, MLList));
	    }
	    catch (...)
	    {
		ERROR("Caught generic exception in the construction" << endl \
                   << "of the matrix-free multilevel preconditioner," << endl \
		   << "maybe this is a memory allocation error.", 1);
	    }
	    graph = null;
	    null_space.resize(0);
	    solver.SetPrecOperator(Prec.get());
#endif
	    // FIXME

	    for (int i = 0; i < B->MyLength(); ++i)
		if ((*B)[i] != 0.0)
		    (*X)[i] = (*B)[i];
	}
      
#ifdef HAVE_MLLIB
	else if (fe_param.get<std::string>("preconditioner") == "ml")
	{
	    // allocate the ML preconditioner
	    ParameterList MLList;
	    SetDefaults("SA", MLList);
	    MLList.set("smoother: type", "MLS");
	    MLList.set("coarse: max size", 1024);
	    MLList.set("null space: dimension", 6);
	    MLList.set("null space: type", "pre-computed");
	    MLList.set("null space: vectors", &null_space[0]);
	    MLList.set("aggregation: type (level 0)", "METIS");
	    MLList.set("aggregation: type (level 1)", "Uncoupled");
	    MLList.set("aggregation: type (level 2)", "MIS");
	    MLList.set("aggregation: nodes per aggregate", 50);
	    MLList.set("low memory usage", true);
	    if (!fe_param.get<bool>("verbose"))
		MLList.set("output", 0);
	    Prec = rcp(new MultiLevelPreconditioner(*ep->GetMatrix(), MLList, true));
	    solver.SetPrecOperator(Prec.get());
	    null_space.resize(0);
	} 
#endif
#ifdef HAVE_IFPACKLIB
	else if (fe_param.get<std::string>("preconditioner") == "ic") {
	    ParameterList SchwarzList;
	    SchwarzList.set("schwarz: reordering type", "metis");
	    SchwarzList.set("fact: level-of-fill", fe_param.get<int>("ic: fill-in"));
	    Ifpack_Preconditioner* IFPACKPrec = new Ifpack_AdditiveSchwarz<Ifpack_IC>(ep->GetMatrix());
	    IFPACKPrec->SetParameters(SchwarzList); 
	    IFPACKPrec->Initialize();
	    IFPACKPrec->Compute();
	    Prec = rcp(IFPACKPrec);
	    solver.SetPrecOperator(Prec.get());
	} 
#endif
	else {
	    solver.SetAztecOption(AZ_precond, AZ_none);
	    WARN("No preconditioner used");
	}

	INFO("Time used to build preconditioner: " << timer.WallTime() - temp_time);
	memused = meminfo() - memused;
	INFO("Mem used to hold preconditioner: " << memused << "MB");
      
        ep->GetLHS()->PutScalar(0.0);

        solver.Iterate(fe_param.get<int>("iteration limit"), fe_param.get<double>("tolerance"));
      
        // need to re-impose the boundary conditions for restrained nodes
        if (fe_param.get<std::string>("preconditioner") == "matrixfree" ||
            fe_param.get<std::string>("preconditioner") == "Jacobi")
        {
        for (int i = 0; i < BCs->MyLength(); ++i)
          if ((*BCs)[0][i] != 0.0)
            (*ep->GetLHS())[0][i] = (*BCs)[0][i];
        }

	diagonal = null;
	Prec = null;
    
	//------------------------------Print solution-------------------------
	temp_time = timer.WallTime();
	
	// VectorWriter* vw = 0;
// 	if (hdf5) {
// 	    //For hdf5, write into the same file
// 	    vw = new IBT_DispWriter<HDF5_GWriter>(fe_param.get<std::string>("input"), comm);
// 	} else {
// 	    vw = new IBT_DispWriter<C_ASCII_GWriter>(ifile_base+".disp", comm);
// 	}
	
// 	ep->PrintLHS(vw);
// 	delete vw;
	
	SolutionWriter* sw = new IBT_SolutionWriter<HDF5_GWriter>(fe_param.get<std::string>("input"), comm);
	ep->PrintSolution(sw, fe_param.get<Epetra_SerialDenseMatrix>("material properties"));
	delete sw;

	if (fe_param.get<bool>("print for MEDIT")) {
	    std::string bbfile = ifile_base + ".medit.bb";
	    MEDIT_DispWriter* dispw = new MEDIT_DispWriter(bbfile, comm, fe_param.get<int>("#dimensions") );
	    ep->PrintLHS(dispw);
	    delete dispw;
	}
	
	INFO("output time: " << timer.WallTime() -temp_time);
	
    }
  
    comm.Barrier();
    INFO("Total time used: " << timer.WallTime() -start_time);
    
// ======================= //
// Finalize MPI and exit //
// ----------------------- //
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    
    exit( EXIT_SUCCESS );
    
}
