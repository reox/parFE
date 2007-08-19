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

#ifndef PARFE_CONFIGDEFS_H
#define PARFE_CONFIGDEFS_H

#ifndef __cplusplus
#define __cplusplus
#endif

#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include <parfe_config.h>

/* #ifdef HAVE_CSTDLIB */
/* #include <cstdlib> */
/* #else */
/* #include <stdlib.h> */
/* #endif */

/* #ifdef HAVE_CSTDIO */
/* #include <cstdio> */
/* #else */
/* #include <stdio.h> */
/* #endif */

/* #ifdef HAVE_CASSERT */
/* #include <cassert> */
/* #else */
/* #include <assert.h> */
/* #endif */

/* #ifdef HAVE_STRING */
/* #include <string> */
/* #else */
/* #include <string.h> */
/* #endif */

/* #ifdef HAVE_IOSTREAM */
/* #include <iostream> */
/* #else */
/* #include <iostream.h> */
/* #endif */

/* #ifdef HAVE_CMATH */
/* #include <cmath> */
/* #else */
/* #include <math.h> */
/* #endif */
using namespace std;

#endif /* PARFE_CONFIGDEFS_H */
