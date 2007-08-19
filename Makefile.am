#
# ParFE: a micro-FE solver for trabecular bone modeling
# Copyright (C) 2006, Uche Mennel and Marzio Sala
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
# 02110-1301, USA.
#

AUTOMAKE_OPTIONS = foreign

ACLOCAL	= aclocal -I config

include config/aminclude.am


MOSTLYCLEANFILES = $(DX_CLEANFILES)


EXTRA_DIST = \
  README \
  ChangeLog \
  LICENCE \
  config/depcomp \
  config/install-sh \
  config/missing \
  config/mkinstalldirs \
  config/config.guess \
  config/config.sub \
  bootstrap \
  mesh/simple_test.mesh \
  mesh/cube_s01_fea.mesh.h5 \
  doc/index.doc \
  doc/invoke-configure-xt3.sh \
  doc/invoke-configure-mac.sh \
  $(DX_CONFIG)

AUX_DIST = config/install-sh config/missing config/mkinstalldirs

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	src/common/config-h.in src/common/stamp-h.in \
	$(AUX_DIST) 

all-local:
	@echo "ParFE built successfully"

SUBDIRS = src example test
