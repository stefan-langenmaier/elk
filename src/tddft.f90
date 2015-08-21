
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddft
use modmain
use modtddft
implicit none
! initialise global variables
call init0
call init1
! read time-dependent A-field from file
call readafieldt


return
end subroutine

