!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module bc_tgv

  use MPI

  use decomp_2d_constants, only : mytype, real_type
  use decomp_2d_mpi, only : nrank, nproc, decomp_2d_warning, decomp_2d_abort
  use decomp_2d, only : xsize, xstart
  use x3d_precision, only : pi, twopi
  use variables, only : nx, ny, nz
  use param, only : dx, dy, dz, xlx, yly, zlz, dt, xnu
  use param, only : zero, one, two

  implicit none

  integer, save :: tgv_iounit
  integer, parameter :: tgv_verbose = 0
  character(len=*), parameter :: tgv_file = "tgv.dat"

  ! Make everything private unless declared public
  private
  public :: tgv_boot, &
            tgv_listing, &
            tgv_init, &
            tgv_postprocess, &
            tgv_finalize

contains

  !
  ! Initialize case-specific parameters and IO
  !
  subroutine tgv_boot()

    use variables, only : nxm, nym, nzm
    use param, only : dt, dx, dy, dz, dx2, dy2, dz2, &
                      xlx, yly, zlz, re, xnu, one

    implicit none

    ! Local variable
    logical :: fexists

    ! Default domain size is 2 pi x 2 pi x 2 pi
    ! This should be inside input.i3d
    xlx = twopi
    yly = twopi
    zlz = twopi
    dx = xlx / real(nxm, mytype)
    dy = yly / real(nym, mytype)
    dz = zlz / real(nzm, mytype)
    dx2 = dx * dx
    dy2 = dy * dy
    dz2 = dz * dz

    ! Default time step : CFL = 0.2 and U = 1
    ! This should be inside input.i3d ?
    dt = 0.1_mytype * dx

    ! Default Re is 1600 when Re = 0
    ! This should be inside input.i3d ?
    if (abs(re) <= epsilon(re)) then
      re = 1600._mytype
      xnu = one / re
    endif

    ! Check if the file is present and get IO unit
    if (nrank == 0) then

      inquire(file=tgv_file, exist=fexists)
      if (fexists) then
        if (tgv_verbose > 0) call decomp_2d_warning(1, "TGV: file "//tgv_file//" replaced.")
        open(newunit=tgv_iounit, file=tgv_file, action='write', status='replace')
      else
        open(newunit=tgv_iounit, file=tgv_file, action='write', status='new')
      endif

    endif

  end subroutine tgv_boot

  !
  ! Case-specific parameters in the listing
  !
  subroutine tgv_listing()

    implicit none

    if (nrank == 0) then
      write(*,*)' 3D TGV test case'
      write(*,*)'==========================================================='
    endif

  end subroutine tgv_listing

  !
  ! Initialize 3D fields
  !
  subroutine tgv_init(ux1, uy1, uz1)

    implicit none

    ! Arguments
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

    ! Local variables
    real(mytype) :: x, y, z
    integer :: i, j, k, ip, it
    integer :: xe1, xe2, xe3
    integer :: xs1, xs2, xs3

    xe1 = xsize(1)
    xe2 = xsize(2)
    xe3 = xsize(3)
    xs1 = xstart(1)
    xs2 = xstart(2)
    xs3 = xstart(3)

    ! Initial velocity
    !$acc kernels default(present)
    do concurrent (k=1:xe3, j=1:xe2, i=1:xe1)
       z = (k + xs3 - 2) * dz
       y = (j + xs2 - 2) * dy
       x = (i + xs1 - 2) * dx
       ux1(i, j, k) = sin(twopi * (x / xlx)) * cos(twopi * (y / yly)) * cos(twopi * (z / zlz))
       uy1(i, j, k) = -cos(twopi * (x / xlx)) * sin(twopi * (y / yly)) * cos(twopi * (z / zlz))
       uz1(i, j, k) = zero
    enddo
    !$acc end kernels

    ! Check initial TKE, dissipation and enstrophy
    call tgv_postprocess(ux1, uy1, uz1, 1)

  endsubroutine tgv_init

  !
  ! Compute and log various statistics
  !
  subroutine tgv_postprocess(ux1, uy1, uz1, ndt)

    use decomp_2d, only : xsize, ysize, zsize
    use param, only : half, two, xnu, dt
    use variables, only : nx, ny, nz, ppy
    use var, only : ux2,ux3
    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1
    use var, only : ta2
    use var, only : ta3
    use x3d_derive
    use x3d_transpose
    use x3d_operator_1d

    implicit none

    ! Arguments
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    integer, intent(in) :: ndt

    ! Local variables
    real(mytype) :: x, y, z, tke, eps, eps2, enst
    integer :: code
    integer :: i, j, k , n1, n2, n3

    n1 = xsize(1)
    n2 = xsize(2)
    n3 = xsize(3)

    ! Collect statistics at each time step currently
    if (mod(ndt, 1) /= 0) return

    ! This is needed to compute derivatives
    !$acc data create(ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,ux2,ta2,ux3,ta3) present(ux1,uy1,uz1)
    ! Compute X derivative
    call derx (ta1,ux1,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    call derx (tb1,uy1,x3d_op_derxp,xsize(1),xsize(2),xsize(3))
    call derx (tc1,uz1,x3d_op_derxp,xsize(1),xsize(2),xsize(3))
    
    ! For transpose we use only dummy array ux => BE CAREFUL TO HAVE THE CORRECT VELOCITY
    ! Compute Y derivative and transpose back to X
    call x3d_transpose_x_to_y(ux1, ux2)
    call dery (ta2,ux2,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))
    call x3d_transpose_y_to_x(ta2, td1)
    
    call x3d_transpose_x_to_y(uy1, ux2)
    call dery (ta2,ux2,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    call x3d_transpose_y_to_x(ta2, te1)
    
    call x3d_transpose_x_to_y(uz1, ux2)
    call dery (ta2,ux2,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))
    call x3d_transpose_y_to_x(ta2, tf1)

    ! Compute Z derivative and transpose back to X
    call x3d_transpose_x_to_y(ux1, ux2)
    call x3d_transpose_y_to_z(ux2, ux3)
    call derz (ta3,ux3,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    call x3d_transpose_z_to_y(ta3, ta2)
    call x3d_transpose_y_to_x(ta2, tg1)

    call x3d_transpose_x_to_y(uy1, ux2)
    call x3d_transpose_y_to_z(ux2, ux3)
    call derz (ta3,ux3,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    call x3d_transpose_z_to_y(ta3, ta2)
    call x3d_transpose_y_to_x(ta2, th1)

    call x3d_transpose_x_to_y(uz1, ux2)
    call x3d_transpose_y_to_z(ux2, ux3)
    call derz (ta3,ux3,x3d_op_derz ,zsize(1),zsize(2),zsize(3))
    call x3d_transpose_z_to_y(ta3, ta2)
    call x3d_transpose_y_to_x(ta2, ti1)
    !du/dx=ta1   du/dy=td1   du/dz=tg1
    !dv/dx=tb1   dv/dy=te1   dv/dz=th1
    !dw/dx=tc1   dw/dy=tf1   dw/dz=ti1

    ! Space-average of enstrophy
    enst = zero
    !$acc parallel loop default(present) collapse(3) reduction(+:enst)
    do k=1,n3
        do j=1,n2
            do i=1,n1
                enst = enst+  ((tf1(i,j,k)-th1(i,j,k))**2 + &
                               (tg1(i,j,k)-tc1(i,j,k))**2 + &
                               (tb1(i,j,k)-td1(i,j,k))**2)
            enddo
        enddo
    enddo
    !$acc end parallel

    call MPI_ALLREDUCE(MPI_IN_PLACE,enst,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, &
                                        "MPI_ALLREDUCE")
    enst = half * enst / (nx*ny*nz)

    ! Space-average of energy dissipation
    eps = zero
    !$acc parallel loop default(present) collapse(3) reduction(+:eps)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             eps = eps + ((two*ta1(i,j,k))**2 + (two*te1(i,j,k))**2 + (two*ti1(i,j,k))**2 + &
                           two*(td1(i,j,k)+tb1(i,j,k))**2 + &
                           two*(tg1(i,j,k)+tc1(i,j,k))**2 + &
                           two*(th1(i,j,k)+tf1(i,j,k))**2)
          enddo
       enddo
    enddo
    !$acc end parallel

    call MPI_ALLREDUCE(MPI_IN_PLACE,eps,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, &
                                        "MPI_ALLREDUCE")
    eps = half * xnu * eps / (nx*ny*nz)

    ! Space-average of TKE
    tke = zero
    !$acc parallel loop default(present) collapse(3) reduction(+:tke)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             tke = tke + (ux1(i,j,k)**2 + uy1(i,j,k)**2 + uz1(i,j,k)**2)
          enddo
       enddo
    enddo
    !$acc end parallel
    call MPI_ALLREDUCE(MPI_IN_PLACE,tke,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, &
                                        "MPI_ALLREDUCE")
    tke = half * tke / (nx*ny*nz)

    ! Compute X second derivative
    call derxx(ta1,ux1,x3d_op_derxx, xsize(1),xsize(2),xsize(3))
    call derxx(tb1,uy1,x3d_op_derxxp,xsize(1),xsize(2),xsize(3))
    call derxx(tc1,uz1,x3d_op_derxxp,xsize(1),xsize(2),xsize(3))

    ! Compute Y second derivative and transpose back to X
    call x3d_transpose_x_to_y(ux1, ux2)
    call deryy(ta2,ux2,x3d_op_deryyp,ysize(1),ysize(2),ysize(3))
    call x3d_transpose_y_to_x(ta2, td1)
    
    call x3d_transpose_x_to_y(uy1, ux2)
    call deryy(ta2,ux2,x3d_op_deryy ,ysize(1),ysize(2),ysize(3))
    call x3d_transpose_y_to_x(ta2, te1)
    
    call x3d_transpose_x_to_y(uz1, ux2)
    call deryy(ta2,ux2,x3d_op_deryyp,ysize(1),ysize(2),ysize(3))
    call x3d_transpose_y_to_x(ta2, tf1)

    ! Compute Z second derivative and transpose back to X
    call x3d_transpose_x_to_y(ux1, ux2)
    call x3d_transpose_y_to_z(ux2, ux3)
    call derzz(ta3,ux3,x3d_op_derzzp,zsize(1),zsize(2),zsize(3))
    call x3d_transpose_z_to_y(ta3, ta2)
    call x3d_transpose_y_to_x(ta2, tg1)
    
    call x3d_transpose_x_to_y(uy1, ux2)
    call x3d_transpose_y_to_z(ux2, ux3)
    call derzz(ta3,ux3,x3d_op_derzzp,zsize(1),zsize(2),zsize(3))
    call x3d_transpose_z_to_y(ta3, ta2)
    call x3d_transpose_y_to_x(ta2, th1)
    
    call x3d_transpose_x_to_y(uz1, ux2)
    call x3d_transpose_y_to_z(ux2, ux3)
    call derzz(ta3,ux3,x3d_op_derzz ,zsize(1),zsize(2),zsize(3))
    call x3d_transpose_z_to_y(ta3, ta2)
    call x3d_transpose_y_to_x(ta2, ti1)
    !!$acc wait
    !d²u/dxx=ta1   d²u/dyy=td1   d²u/dzz=tg1
    !d²v/dxx=tb1   d²v/dyy=te1   d²v/dzz=th1
    !d²w/dxx=tc1   d²w/dyy=tf1   d²w/dzz=ti1

    ! Space average of energy dissipation with second derivatives
    eps2 = zero
    !$acc parallel loop default(present) collapse(3) reduction(+:eps2)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             eps2 = eps2 + (ux1(i,j,k)*(ta1(i,j,k)+td1(i,j,k)+tg1(i,j,k)) + &
                            uy1(i,j,k)*(tb1(i,j,k)+te1(i,j,k)+th1(i,j,k)) + &
                            uz1(i,j,k)*(tc1(i,j,k)+tf1(i,j,k)+ti1(i,j,k)))
          enddo
       enddo
    enddo
    !$acc end parallel
    call MPI_ALLREDUCE(MPI_IN_PLACE,eps2,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, &
                                        "MPI_ALLREDUCE")
    eps2 = -xnu * eps2 / (nx*ny*nz)
    !$acc end data

    ! Log
    if (nrank == 0) then
      write(tgv_iounit,'(20e20.12)') (ndt-1)*dt, tke, eps, eps2, enst
      call flush(tgv_iounit)
    endif

  end subroutine tgv_postprocess

  !
  ! Finalize case-specific IO
  !
  subroutine tgv_finalize()

    implicit none

    if (nrank == 0) close(tgv_iounit)

  end subroutine tgv_finalize

end module bc_tgv
