!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @file
!> @ingroup mpp_mod

#if defined(use_libMPI)
#define SYSTEM_CLOCK system_clock_mpi

!#######################################################################
subroutine system_clock_mpi( count, count_rate, count_max )
! There can be one ONE baseline count0 and this routine is
! included in multiple places.
!mimics F90 SYSTEM_CLOCK intrinsic
      integer(i8_kind), intent(out), optional :: count, count_rate, count_max
!count must return a number between 0 and count_max
      integer(i8_kind), parameter :: maxtick=HUGE(count_max)
      if(first_call_system_clock_mpi)then
         first_call_system_clock_mpi=.false.
         mpi_count0 = MPI_WTime()
         mpi_tick_rate = real(1.d0/MPI_WTick(), r8_kind)
      endif
      if( PRESENT(count) )then
          count = int((MPI_WTime()-mpi_count0)*mpi_tick_rate, i8_kind)
      end if
      if( PRESENT(count_rate) )then
          count_rate = int(mpi_tick_rate, i8_kind)
      end if
      if( PRESENT(count_max) )then
          count_max = maxtick-1
      end if
      return
    end subroutine system_clock_mpi

#else
#define SYSTEM_CLOCK system_clock_default
subroutine system_clock_default( count, count_rate, count_max )
!mimics F90 SYSTEM_CLOCK intrinsic
      integer(i8_kind), optional :: count, count_rate, count_max
!count must return a number between 0 and count_max
      integer                      :: count_int, count_rate_int, count_max_int
      call system_clock( count_int, count_rate_int, count_max_int)
      if( PRESENT(count) )      count      = count_int
      if( PRESENT(count_rate) ) count_rate = count_rate_int
      if( PRESENT(count_max) )  count_max  = count_max_int
      return
    end subroutine system_clock_default
#endif
