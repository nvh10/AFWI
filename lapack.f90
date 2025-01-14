    subroutine InverseMatrixD(N, A, invA)
    implicit none
    integer N, IPIV(N), INFO,ix,iz,LWORK
    double precision A(N,N), invA(N,N)
    double precision,allocatable ::WORK(:)

    LWORK=128*N
    allocate(WORK(LWORK))

    !invA(:,:) = A(:,:)
    !call DGETRF (N, N, invA(:,:), N, IPIV(:), INFO)
    !call DGETRI (N, invA(:,:), N, IPIV(:), WORK(:), N, INFO)
    invA = A
    call DSYTRF('U',N,invA(:,:),N,IPIV,WORK,LWORK,INFO)
    call DSYTRI('U',N,INVA(:,:),N,IPIV,WORK,INFO)
    do ix=1,N-1
        INVA(ix+1:N,ix)=INVA(ix,ix+1:N)
    enddo
    end subroutine InverseMatrixD

    subroutine solveAxB(N,M, A, B,X)
    implicit none
    integer N, IPIV(N), INFO,ix,iz,LWORK,M
    double precision A(N,N), invA(N,N),B(N,M),X(N,M)
    double precision,allocatable ::WORK(:)

    LWORK=128*N
    allocate(WORK(LWORK))
    invA = A
    X=B
    call DSYTRF('U',N,invA(:,:),N,IPIV,WORK,LWORK,INFO)

    Call DSYTRS('U', N, M, invA, N, IPIV, X, N,INFO)



    end subroutine solveAxB

    subroutine solveAxB_Z(N,M, A, B,X)
    implicit none
    integer N, IPIV(N), INFO,ix,iz,LWORK,M
    double complex A(N,N), invA(N,N),B(N,M),X(N,M)

    double complex,allocatable ::WORK(:)

    LWORK=128*N
    allocate(WORK(LWORK))
    invA = A
    X=B
    
    call ZPOTRF('U',N,invA,N,info)
    call ZPOTRS('U', N, 1, invA, N, X, N, info )

    
    !call ZSYTRF('U',N,invA(:,:),N,IPIV,WORK,LWORK,INFO)
    !
    !Call ZSYTRS('U', N, M, invA, N, IPIV, X, N,INFO)
    !


    end subroutine solveAxB_Z




    subroutine InverseMatrixZ(N, A, invA)
    implicit none
    integer N, IPIV(N), INFO,LWORK,ix,iz
    double complex A(N,N), invA(N,N),cholA(N,N)
    double complex,allocatable ::WORK(:)
    
LWORK=128*N
    allocate(WORK(LWORK))
 
    invA(:,:) = A(:,:)
    !$omp parallel do
    do iz=1,N
        do ix=iz+1,N
            invA(ix,iz)=0d0
        enddo
    enddo
    !$omp end parallel do
    
     
    call ZPOTRF('U',N,invA,N,info)
    call ZPOTRI('U', N, invA, N, info )
    do iz=1,N
        do ix=iz+1,N
            invA(ix,iz)=dconjg(invA(iz,ix))
        enddo
    enddo

    
    !   allocate(WORK(N))
    !call ZGETRF (N, N, invA(:,:), N, IPIV(:), INFO)
    !call ZGETRI (N, invA(:,:), N, IPIV(:), WORK(:), N, INFO)


    !LWORK=128*N
    !allocate(WORK(LWORK))
    !invA = A
    !call ZSYTRF('U',N,invA(:,:),N,IPIV,WORK,LWORK,INFO)
    !call ZSYTRI('U',N,INVA(:,:),N,IPIV,WORK,INFO)
    !do ix=1,N-1
    !    INVA(ix+1:N,ix)=INVA(ix,ix+1:N)
    !enddo
    !invA=A
    !
    !!$omp parallel do
    !do iz=1,N
    !    do ix=iz+1,N
    !        INVA(ix,iz)=dcmplx(1d0,0d0)
    !    enddo
    !enddo
    !!$omp end parallel do
    !call ZPOTRF('U',N,INVA,N,info)
    !call ZPOTRI( 'U', N, INVA, N, info )
    !!$omp parallel do
    !do ix=1,N-1
    !    INVA(ix+1:N,ix)=INVA(ix,ix+1:N)
    !enddo
    !!$omp end parallel do


    end subroutine InverseMatrixZ

    subroutine eigenD(N, A, B, ALPHA, BETA, EVEC)
    implicit none
    integer N, N1
    double precision A(N,N), B(N,N), A1(N,N), B1(N,N)
    double precision ALPHAR(N), ALPHAI(N), BETA(N)
    double precision VL(N,N), EVEC(N,N), WORK(8*N)
    double complex ALPHA(N)

    integer INFO, itmp, jtmp
    double precision dswap(N)
    double complex zswap(N)

    A1(:,:) = A(:,:)
    B1(:,:) = B(:,:)

    call DGGEV('N', 'V', N, A1(:,:), N, B1(:,:), N, ALPHAR(:), ALPHAI(:), BETA(:), VL(:,:), N, EVEC(:,:), N, WORK(:), 8*N, INFO)

    ALPHA(:) = dcmplx(ALPHAR(:), ALPHAI(:))

    jtmp = N

    do itmp = N, 1, -1

        if (dabs(BETA(itmp)) == 0.d0) then

            zswap(1) = ALPHA(itmp)
            ALPHA(itmp) = ALPHA(jtmp)
            ALPHA(jtmp) = zswap(1)

            dswap(1) = BETA(itmp)
            BETA(itmp) = BETA(jtmp)
            BETA(jtmp) = dswap(1)

            dswap(:) = EVEC(:,itmp)
            EVEC(:,itmp) = EVEC(:,jtmp)
            EVEC(:,jtmp) = dswap(:)

            jtmp = jtmp - 1

        endif

    enddo

    N1 = jtmp

    do itmp = 1 , N1-1

        do jtmp = itmp+1 , N1

            if (cdabs(ALPHA(itmp)/BETA(itmp)) > cdabs(ALPHA(jtmp)/BETA(jtmp))) then

                zswap(1) = ALPHA(itmp)
                ALPHA(itmp) = ALPHA(jtmp)
                ALPHA(jtmp) = zswap(1)

                dswap(1) = BETA(itmp)
                BETA(itmp) = BETA(jtmp)
                BETA(jtmp) = dswap(1)

                dswap(:) = EVEC(:,itmp)
                EVEC(:,itmp) = EVEC(:,jtmp)
                EVEC(:,jtmp) = dswap(:)

            endif

        enddo

    enddo

    end subroutine eigenD

    subroutine eigenZ(N, A, B, ALPHA, BETA, EVEC)
    implicit none
    integer N, N1
    double complex A(N,N), B(N,N), A1(N,N), B1(N,N)
    double complex ALPHA(N), BETA(N), EVEC(N,N)
    double complex VL(N,N), WORK(4*N)
    double precision RWORK(8*N)

    integer INFO, itmp, jtmp
    double complex zswap(N)

    A1(:,:) = A(:,:)
    B1(:,:) = B(:,:)

    call ZGGEV('N', 'V', N, A1(:,:), N, B1(:,:), N, ALPHA(:), BETA(:), VL(:,:), N, EVEC(:,:), N, WORK(:), 4*N, RWORK(:), INFO)

    jtmp = N

    do itmp = N, 1, -1

        if (cdabs(BETA(itmp)) == 0.d0) then

            zswap(1) = ALPHA(itmp)
            ALPHA(itmp) = ALPHA(jtmp)
            ALPHA(jtmp) = zswap(1)

            zswap(1) = BETA(itmp)
            BETA(itmp) = BETA(jtmp)
            BETA(jtmp) = zswap(1)

            zswap(:) = EVEC(:,itmp)
            EVEC(:,itmp) = EVEC(:,jtmp)
            EVEC(:,jtmp) = zswap(:)

            jtmp = jtmp - 1

        endif

    enddo

    N1 = jtmp

    do itmp = 1 , N1-1

        do jtmp = itmp+1 , N1

            if (cdabs(ALPHA(itmp)/BETA(itmp)) > cdabs(ALPHA(jtmp)/BETA(jtmp))) then

                zswap(1) = ALPHA(itmp)
                ALPHA(itmp) = ALPHA(jtmp)
                ALPHA(jtmp) = zswap(1)

                zswap(1) = BETA(itmp)
                BETA(itmp) = BETA(jtmp)
                BETA(jtmp) = zswap(1)

                zswap(:) = EVEC(:,itmp)
                EVEC(:,itmp) = EVEC(:,jtmp)
                EVEC(:,jtmp) = zswap(:)

            endif

        enddo

    enddo

    end subroutine eigenZ

    subroutine DetterminantD(N, A,detA)
    implicit none
    integer N, IPIV(N), INFO, ii
    double precision A(N,N),  loA(N,N)
    real(kind=16)::detA,detU,detL, detP
    loA(:,:) = A(:,:)
    call dgetrf(N,N,loA,N,IPIV,info)
    detU = 1d0
    do ii=1,N
        detU =  detU*loA(ii,ii)
    enddo
    detL = 1d0
    detP = 1d0

    do ii=1,N
        if (ipiv(ii)/=ii)then
            detP = -detP
        endif
    enddo

    detA = detP*detL*detU
    end subroutine DetterminantD

    subroutine cholesky_lower(N, A,cholA)
    implicit none
    integer N, INFO, ix,iz
    double precision A(N,N),  cholA(N,N)
    cholA=A
    call DPOTRF('L',N,cholA,N,info)
    !$omp parallel do
    do iz=1,N
        do ix=iz+1,N
            cholA(iz,ix)=0d0
        enddo
    enddo
    !$omp end parallel do
    end subroutine cholesky_lower

    subroutine cholesky_lower_Z(N, A,cholA)
    implicit none
    integer N, INFO, ix,iz
    double complex A(N,N),  cholA(N,N)
    cholA=A
    !$omp parallel do
    do iz=1,N
        do ix=iz+1,N
            cholA(ix,iz)=0d0
        enddo
    enddo
    !$omp end parallel do
    call ZPOTRF('U',N,cholA,N,info)

    end subroutine cholesky_lower_Z