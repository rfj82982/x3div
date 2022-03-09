!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

module case

  use param
  use decomp_2d, only : mytype
  use variables

  implicit none

  private ! All functions/subroutines private by default
  public :: init

contains
  !##################################################################
  subroutine init (ux1, uy1, uz1, dux1, duy1, duz1, &
       pp3, px1, py1, pz1)

    use mom, only : vel
    use decomp_2d, only : xsize, ph1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzm, npress) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1, py1, pz1

    integer :: it, is
    integer :: i, j, k 

    !! Zero out the pressure field
    do concurrent (k=1:nzm, j=ph1%zst(2):ph1%zen(2), i=ph1%zst(1):ph1%zen(1))
      pp3(i,j,k,1) = zero
    enddo
    
    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1)) 
      px1(i,j,k) = zero
      py1(i,j,k) = zero
      pz1(i,j,k) = zero
    enddo

    call vel(ux1, uy1, uz1)
    
    !! Setup old arrays
    do it = 1, ntime
      do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1)) 
        dux1(i,j,k,it)=ux1(i,j,k)
        duy1(i,j,k,it)=uy1(i,j,k)
        duz1(i,j,k,it)=uz1(i,j,k)
      enddo
    enddo


  end subroutine init
end module case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! case.f90 ends here
