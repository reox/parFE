/*
 * ParFE: a micro-FE solver for trabecular bone modeling
 * Copyright (C) 2006, Uche Mennel and Marzio Sala
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
 * 02110-1301, USA.
 */

#include "GReader.hpp"
#include <cstring>


//C functions
C_ASCII_GReader::C_ASCII_GReader(MPI_Comm a_comm)
    : file(NULL), comm(a_comm), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
}

C_ASCII_GReader::C_ASCII_GReader(const std::string& a_filename, MPI_Comm a_comm)
  : filename(a_filename), file(NULL), comm(a_comm), mpi_rank(0), mpi_size(1)   //save a copy of the filename
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
  Open();
}

C_ASCII_GReader::~C_ASCII_GReader()
{
  Close();
}

bool C_ASCII_GReader::Open()
{
  file = fopen(filename.c_str(), "r");
  return file;
}

bool C_ASCII_GReader::Open(const std::string& a_filename)
{
  filename = a_filename;
  return Open();
}

void C_ASCII_GReader::Close()
{
  fclose(file);
}

int C_ASCII_GReader::Read(int& i)
{
  char buffer[64];
  char c;
  int count=0, diff;
  char *last = buffer;
  while ( isspace( c = fgetc(file)) );
  buffer[count++] = c;
  while ( !isspace(buffer[count++] = fgetc(file)) && count < 64);
  buffer[count] = '\0';
  i = strtol(buffer, &last, 10);
  diff = last - buffer;
  fseek(file, diff - count, SEEK_CUR);
  return diff;
}

int C_ASCII_GReader::Read(double& d)
{
  char buffer[64];
  int count=0, diff;
  char *last = buffer;
  char c;
  while ( isspace( c = fgetc(file)) );
  buffer[count++] = c;
  while ( !isspace(buffer[count++] = fgetc(file)) && count < 64);
  buffer[count] = '\0';
  d = strtod(buffer, &last);
  diff = last - buffer;
  fseek(file, diff - count, SEEK_CUR);
  return diff;
}

int C_ASCII_GReader::Read(std::string& s)
{
  char buffer[256];
  int count=0;
  char c;
  while ( isspace( c = fgetc(file)) );
  buffer[count++] = c;
  while ( !isspace(buffer[count++] = fgetc(file)) && count < 256);
  ungetc(buffer[--count], file);
  buffer[count] = '\0';
  s = buffer;
  return count;
}

bool C_ASCII_GReader::operator!()
{
  return !file;
}

int C_ASCII_GReader::Skip(char comment)
{
  char c;
  char buffer[1024];
  int count = 0;
  while (isspace( c = fgetc(file)) );
  while (c == comment) {
      fgets(buffer, 1024, file);
      while (isspace( c = fgetc(file)) );
      ++count;
  }
  ungetc(c, file);
  return count;
}

bool C_ASCII_GReader::Select(const std::string& s)
{
  //nothing to do here
  return true;
}

bool C_ASCII_GReader::Seek(const std::string& s)
{
  if (marks.count(s) == 0)
    fgetpos(file, &marks[s]);
  return (fsetpos(file, &marks[s]) == 0);
}

//CPP functions
CPP_ASCII_GReader::CPP_ASCII_GReader(MPI_Comm a_comm)
    :comm(a_comm), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
}

CPP_ASCII_GReader::CPP_ASCII_GReader(const std::string& a_filename, MPI_Comm a_comm)
    : filename(a_filename), comm(a_comm), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
  Open();
}

CPP_ASCII_GReader::~CPP_ASCII_GReader()
{
  Close();
}

bool CPP_ASCII_GReader::Open()
{
  file.open(filename.c_str());
  return file.good();
}

bool CPP_ASCII_GReader::Open(const std::string& a_filename)
{
  filename = a_filename;
  return Open();
}

void CPP_ASCII_GReader::Close()
{
  file.close();
}

bool CPP_ASCII_GReader::operator!()
{
  return (!file);
}


int CPP_ASCII_GReader::Skip(char comment)
{
  std::string s;
  unsigned char c;
  int count = 0;

  file >> c;
  while (c == comment) { 
    std::getline( file, s );
    file >> c;
    ++count;
  }
  file.putback(c);

  return count;
}


bool CPP_ASCII_GReader::Select(const std::string& s)
{
  //nothing to do here
  return true;
}


bool CPP_ASCII_GReader::Seek(const std::string& s)
{
  if (marks.count(s) == 0)
    marks[s] = file.tellg();
  return file.seekg(marks[s]).good();
}


HDF5_GReader::HDF5_GReader(MPI_Comm a_comm)
    : file(0), plist(H5Pcreate(H5P_FILE_ACCESS)), index(0),  comm(a_comm), file_info(MPI_INFO_NULL), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Info_create(&file_info);
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif 
//Initialize driver specific properties
#ifdef H5_HAVE_PARALLEL
  H5Pset_fapl_mpio(plist, comm, file_info);
#endif
  offset[0] = 0;
  offset[1] = 0;
  stride[0] = 1;
  stride[1] = 1;
}


HDF5_GReader::HDF5_GReader(const std::string& a_filename, MPI_Comm a_comm)
    : file(0), plist(H5Pcreate(H5P_FILE_ACCESS)), filename(a_filename), index(0), comm(a_comm), file_info(MPI_INFO_NULL), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Info_create(&file_info);
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#ifdef __QK_USER__
  //Optimization parameters needed if running on the CRAY XT3
  MPI_Info_set(file_info, "cb_config_list", "*:*");
  MPI_Info_set(file_info, "cb_buffer_size", "1048576");
  MPI_Info_set(file_info, "romio_cb_read", "enable");
#endif
#endif
 //Initialize driver specific properties
#ifdef H5_HAVE_PARALLEL
  H5Pset_fapl_mpio(plist, comm, file_info);
#endif
  stride[0] = 1;
  stride[1] = 1;
  offset[0] = 0;
  offset[1] = 0;
  Open();
}


HDF5_GReader::~HDF5_GReader()
{
  Close();
#ifdef HAVE_MPI
  MPI_Info_free(&file_info);
#endif
}


bool HDF5_GReader::Open(const std::string& a_filename)
{
  filename = a_filename;
  return Open();
}


bool HDF5_GReader::Open()
{
  //Suppress error messges.
  //H5Eset_auto( NULL, NULL );
  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist);
  H5Pclose(plist);
  group = H5Gopen(file, "/");
  plist = H5P_DEFAULT;
  H5Eset_auto( NULL, NULL );
  return (group >= 0);
}


void HDF5_GReader::Close()
{

  H5Gclose(group);
  H5Fclose(file);
}


bool HDF5_GReader::operator!()
{
  return (file<0);
}


int HDF5_GReader::Skip(char comment)
{
  //nothing to do here
  return 0;
}


bool HDF5_GReader::Select(std::string s)
{
  hid_t temp;
  //Test group
  temp = group;
  group = H5Gopen(temp, s.c_str());
  if (group < 0) {
    return false;
  }
  H5Gclose(temp);
  return true;
}


int HDF5_GReader::Read(const std::string& name, hid_t type, void* data, int num_global_elems, int num_my_elems, int elem_size, int my_offset)
{
  hsize_t len = num_my_elems*elem_size;
  hsize_t my_dims[2];
  my_dims[0] = num_my_elems;
  my_dims[1] = elem_size;
  
  hid_t dataset = H5Dopen(group, name.c_str());
  if (dataset < 0)
    return -1;
  
  //Get dataspace of the dataset.
  hid_t dataspace = H5Dget_space(dataset);
  offset[0] = my_offset;
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, stride,  my_dims, NULL );

  hid_t memspace = H5Screate_simple( 1, &len, NULL );
  if (memspace <= 0) {
    hsize_t one=1;
    memspace = H5Screate_simple( 1, &one, NULL );
    H5Sselect_none(memspace);
  }
  
  H5Dread(dataset, type, memspace, dataspace, plist, data);
  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Dclose(dataset);
  return len;
}


int HDF5_GReader::Read(const std::string& name, hid_t type, void* data)
{
  hid_t dataset = H5Dopen(group, name.c_str());
  if (dataset < 0)
    return -1;
  H5Dread(dataset, type, H5S_ALL, H5S_ALL, plist, data);
  H5Dclose(dataset);
  return 1;
}


int HDF5_GReader::Read(const std::string& name, std::string& data)
{
  const int buf_size = 80;
  int res;
  char c_str[buf_size];
  if (mpi_rank == 0)
    res = Read(name, getNativeType(data), c_str);
  MPI_Bcast(c_str, buf_size, MPI_CHAR, 0, comm);
  data = c_str;
  return res;
}


hid_t HDF5_GReader::getNativeType(int, hsize_t i) {
  if (i == 1)
    return H5T_NATIVE_INT;
  return H5Tarray_create(H5T_NATIVE_INT, 1, &i, NULL);
}

hid_t HDF5_GReader::getNativeType(double, hsize_t i) {
  if (i == 1)
    return H5T_NATIVE_DOUBLE;
  return H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &i, NULL);
}

hid_t HDF5_GReader::getNativeType(float, hsize_t i) {
  if (i==1)
    return H5T_NATIVE_FLOAT; 
  return H5Tarray_create(H5T_NATIVE_FLOAT, 1, &i, NULL);
}

hid_t HDF5_GReader::getNativeType(const std::string&, hsize_t i) {
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, 32);
  if (i==1)
    return type; 
  return H5Tarray_create(type, 1, &i, NULL);
}

void HDF5_GReader::copy(void* dest, const void* src, int nblocks, int block_size, int stride)
{
  for(int i=0; i<nblocks; ++i) {
    memcpy(dest, src, block_size);
    dest = (void*)((unsigned long)(dest)+block_size);
    src = (void*)((unsigned long)(src)+stride);
  }
}
