#!/bin/bash
#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the SHiELD Build System.
#*
#* The SHiELD Build System free software: you can redistribute it
#* and/or modify it under the terms of the
#* GNU Lesser General Public License as published by the
#* Free Software Foundation, either version 3 of the License, or
#* (at your option) any later version.
#*
#* The SHiELD Build System distributed in the hope that it will be
#* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
#* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#* See the GNU General Public License for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with theSHiELD Build System
#* If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************
#
#  DISCLAIMER: This script is provided as-is and as such is unsupported.
#

      echo " Stellar environment "
      module list
      
     # . ${MODULESHOME}/init/sh
      echo $SHELL
      module load intel/2021.1.2
      module load openmpi/intel-2021.1/4.1.2
      module load netcdf/intel-2021.1/hdf5-1.10.6/4.7.4
      module load hdf5/intel-2021.1/1.10.6
      

      set FC=mpif90
      set CC=mpicc
      set CXX=mpicxx
      set LD=mpif90
      set TEMPLATE=site/intel.mk

# highest level of AVX support
      set AVX_LEVEL=-march=core-avx2




#      export NETCDF_DIR=${NETCDFDIR}
      
      # make your compiler selections here
#      export FC=mpif90
#      export CC=mpicc
#      export CXX=mpicxx
#      export LD=mpif90
#      export TEMPLATE=site/intel.mk
#
## highest level of AVX support
#      export AVX_LEVEL=-march=core-avx2


      module list
      echo " DONE "

