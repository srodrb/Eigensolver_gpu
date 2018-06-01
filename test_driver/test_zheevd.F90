!
! Copyright (c) 2016, NVIDIA CORPORATION. All rights reserved.
! 
! 
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the "Software"),
! to deal in the Software without restriction, including without limitation
! the rights to use, copy, modify, merge, publish, distribute, sublicense,
! and/or sell copies of the Software, and to permit persons to whom the
! Software is furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.
!

module funcs
  contains

  ! Creates pseudo-random positive-definite hermetian matrix
  subroutine create_random_hermetian_pd(A, N)
    use cudafor
    use cublas
    complex(8), allocatable, dimension(:,:)         :: A, temp
    complex(8), allocatable, dimension(:,:), device :: A_d, temp_d
    complex(8)                                      :: val
    real(8)                                         :: rv, iv
    integer                                         :: i, j, N

    allocate(A(N,N))
    allocate(temp(N,N))

    ! Create general hermetian temp
    do j = 1, N
      do i = 1, N
        if (i > j) then
          call random_number(rv)
          call random_number(iv)
          temp(i,j) = cmplx(rv, iv, 8)
          temp(j,i) = conjg(temp(i,j))
        else if (i == j) then
          call random_number(rv)
          temp(i,j) = rv
        end if
      end do
    end do

    allocate(A_d, source = A)
    allocate(temp_d, source = temp)

    ! Multiply temp by conjugate transpose of temp to get positive definite hermetian A
    call cublaszgemm('N', 'C', N, N, N, cmplx(1.0, 0.0, 8), temp_d, N, temp_d, N, cmplx(0.0, 0.0, 8), A_d, N)

    A = A_d
    deallocate(temp)
    deallocate(A_d)
    deallocate(temp_d)
        
  end subroutine
end module funcs

program main
  use cudafor
  use cublas
  use cusolverDn
  use eigsolve_vars, ONLY: init_eigsolve_gpu
  use zheevd_gpu
  use nvtx_inters
  use funcs
  use compare_utils
  implicit none
  
  integer                                         :: N, M, i, j, info, lda, istat
  integer                                         :: lwork_d, lrwork_d, lwork, lrwork, liwork, il, iu
  character(len=20)                               :: arg
  real(8)                                         :: ts, te, wallclock
  complex(8), dimension(:,:), allocatable         :: A1, A2, Aref
  complex(8), dimension(:,:), allocatable, pinned :: Z1, Z2
  complex(8), dimension(:,:), allocatable, device :: A2_d, Z2_d
  complex(8), dimension(:), allocatable, pinned   :: work !, rwork !rwork is double precision for zheevd
  real(8), dimension(:), allocatable, pinned      :: w1, w2, rwork
  integer, dimension(:), allocatable, pinned      :: iwork
  complex(8), dimension(:), allocatable, device   :: work_d
  real(8), dimension(:), allocatable, device      :: w2_d, rwork_d
  integer, device                                 :: devInfo_d
  type(cusolverDnHandle)                          :: h

  ! Parse command line arguments
  i = command_argument_count()

  if (i >= 1) then
    ! If N is provided, generate random symmetric matrices for A
    print*, "Using randomly-generated matrices..."
    call get_command_argument(1, arg)
    read(arg, *)  N
    lda = N

    ! Create random positive-definite hermetian matrices on host
    call create_random_hermetian_pd(Aref, N)

  else
    print*, "Usage:\n\t ./main [N]"
    call exit
  endif

  print*, "Running with N = ", N

  ! Allocate/Copy matrices to device
  allocate(A1, source = Aref)
  allocate(A2, source = Aref)
  allocate(A2_d, source = Aref)
  allocate(Z1, source = Aref)
  allocate(Z2, source = Aref)
  allocate(Z2_d, source = Aref)

  allocate(w1(N), w2(N))
  allocate(w2_d, source = w2)

  ! Initialize solvers
  call init_eigsolve_gpu()

  istat = cublasInit
  if (istat /= CUBLAS_STATUS_SUCCESS) write(*,*) 'cublas intialization failed'

  istat = cusolverDnCreate(h)
  if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'handle creation failed'

#ifdef HAVE_MAGMA
  call magmaf_init
#endif


  !! Solving generalized eigenproblem using DSYGVD
  ! CASE 1: CPU _____________________________________________
  !print*
  !print*, "CPU_____________________"
  lwork = 2*N + N*N
  lrwork = 1 + 5*N + 2*N*N
  liwork = 3 + 5*N
  !allocate(iwork(liwork))
  !allocate(rwork(lrwork))
  !allocate(work(lwork))
  !call zheevd('V', 'U', N, A1, lda, w1, work, -1, rwork, -1, iwork, -1, istat)
  !if (istat /= 0) write(*,*) 'CPU zheevd worksize failed'
  !lwork = work(1);; liwork = iwork(1);; lrwork = rwork(1);
  !deallocate(iwork, rwork, work)
  allocate(iwork(liwork), rwork(lrwork), work(lwork) )

  A1 = Aref
  ! Run once before timing
  !call zheevd('V', 'U', N, A1, lda, w1, work, lwork, iwork, liwork, istat)
  ! call zheevd('V', 'U', N, A1, lda, w1, work, lwork, rwork, lrwork, iwork, liwork, istat)
  !if (istat /= 0) write(*,*) 'CPU zheevd failed. istat = ', istat

  A1 = Aref
  !ts = wallclock()
  call nvtxStartRange("CPU ZHEEVD",1)
  !call zheevd('V', 'U', N, A1, lda, w1, work, lwork, rwork, lrwork, iwork, liwork, istat)
  call nvtxEndRange
  !te = wallclock()
  !if (istat /= 0) write(*,*) 'CPU zheevd failed. istat = ', istat

  !print*, "\tTime for CPU zheevd = ", (te - ts)*1000.0
  !print*

  ! CASE 4: using CUSTOM ____________________________________________________________________
  print*
  print*, "CUSTOM_____________________"
  A2 = Aref
  w2 = 0
  A2_d = A2
  w2_d = w2
  il = 1
  iu = N

  deallocate(work, rwork, iwork)
  lwork  = N
  lrwork = 1+5*N+2*N*N
  liwork = 3+5*N
  allocate(work(lwork), iwork(liwork), rwork(lrwork))

  deallocate(work_d)
  lwork_d = 2*64*64 + 66 * N
  lrwork_d = N
  allocate(work_d(1*lwork_d))
  allocate(rwork_d(1*lwork_d))

  ts = wallclock()
  call nvtxStartRange("Custom",0)

  call zheevd_gpu( 'V',        & !jobz
                   'U',        & !uplo
                   il,         & !il
                   iu,         & !ui
                   N,          & !N
                   A2_d,       & !A
                   lda,        & !lda
                   Z2_d,       & !Z
                   lda,        & !ldz
                   w2_d,       & !w
                   work_d,     & !work
                   lwork_d,    & !lwork
                   rwork_d,    & !rwork
                   lrwork_d,   & !lrwork
                   work,       & !work_h
                   lwork,      & !lwork_h
                   rwork,      & !rwork_h
                   lrwork,     & !lrwork_h
                   iwork,      & !iwork_h
                   liwork,     & !liwork_h,
                   Z2,         & !Z_h
                   lda,        & !ldz_h
                   w2,         & !w_h
                   istat)        !info

  call nvtxEndRange
  te = wallclock()

  ! print*, "evalues/evector accuracy: (compared to CPU results)"
  ! call compare(w1, w2, iu)
  ! call compare(A1, Z2, N, iu)
  ! print*

  print*, "Time for CUSTOM zheevd/x = ", (te - ts)*1000.0
  if (istat /= 0) write(*,*) 'zheevd_gpu failed'

end program
