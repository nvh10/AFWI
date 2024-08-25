!    subroutine pardiso_cluster(nbc2,DSN22,ia,ja,P2,Dis2,nmatrl,msglvl,phase)
!      use mkl_cluster_sparse_solver
!implicit none
!include 'mpif.h'
!    double complex,allocatable :: DSN22(:)
!    integer,allocatable :: ia(:), ja(:)
!    double complex,allocatable :: Dis2(:,:),P2(:,:)
!    TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), allocatable  :: pt(:)
!    integer maxfct, mnum, mtype, phase, nrhs, error, msglvl,nbc2,nmatrl,i,idum(1)
!    integer*4 mpi_stat, rank
!    integer(4) MKL_COMM
!
!    integer, allocatable :: iparm( : )
!    allocate (ia(nbc2+1))
!    allocate (ja(ia(nbc2+1)-1))
!    allocate( Dis2(nbc2,nmatrl), P2(nbc2,nmatrl), DSN22(ia(nbc2+1)-1))
!    MKL_COMM=MPI_COMM_WORLD
!    call mpi_init(mpi_stat)
!    call mpi_comm_rank(MKL_COMM, rank, mpi_stat)
!    allocate( pt ( 64 ) )
!    allocate( iparm ( 64 ) )
!
!    do i = 1, 64
!        iparm(i) = 0
!    end do
!    do i = 1, 64
!        pt(i)%dummy = 0
!    end do
!
!    iparm(1) = 1 ! no solver default
!    iparm(2) = 3 ! fill-in reordering from METIS
!    iparm(6) = 0 ! =0 solution on the first n compoments of x
!    iparm(8) = 2 ! numbers of iterative refinement steps
!    iparm(10) = 13 ! perturbe the pivot elements with 1E-13
!    iparm(11) = 0 ! use nonsymmetric permutation and scaling MPS
!    iparm(13) = 1 ! maximum weighted matching algorithm is switched-off
!    iparm(40) = 0 ! Input: matrix/rhs/solution stored on master
!    error = 0  ! initialize error flag
!    mtype = 13  ! complex and nonsymmetric
!    mnum=1
!    !phase = 13 !Analysis, numerical factorization, solve, iterative refinement
!    !pardiso(pt, 1, 1, 13, 13, nbc(2), DSN22, ia, ja, idum, 1, iparm, 1, P2(:,1), Dis2, ierr)
!    call cluster_sparse_solver (pt, maxfct, mnum, mtype, phase, nbc2, DSN22, ia, ja,idum, nmatrl, iparm, msglvl, P2, DIS2, MKL_COMM, error )
!
!
!    if ( allocated( ia ) )      deallocate( ia )
!    if ( allocated( ja ) )      deallocate( ja )
!    if ( allocated( DSN22 ) )       deallocate( DSN22 )
!    if ( allocated( P2 ) )       deallocate( P2 )
!    if ( allocated( Dis2 ) )       deallocate( Dis2 )
!
!    if ( allocated( pt ) )      deallocate( pt )
!    if ( allocated( iparm ) )   deallocate( iparm )
!
!
!    call mpi_finalize(mpi_stat)
!
!    endsubroutine pardiso_cluster
