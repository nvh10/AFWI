    subroutine Cmd_Cdd_ensemble(GA,CS0,nmatrl,nobs,Cmd,Ldd,invcovUz,sigma_M2,da_cal,da_obs,CS)
    implicit none
    integer nobs,nmatrl,iz,ix,itmp,INFO,LWORK

    double complex Ga(nobs,nmatrl),Cmd(nmatrl,nobs),Cdm(nobs,nmatrl),Ldd(nobs,nobs),K(nobs,nobs),GaCOV(nobs,nmatrl)
    double precision CS0(nmatrl),invcovUz(nobs),sigma_M2(nmatrl,nmatrl)
    double complex,allocatable ::WORK(:)
    double complex dtheta_c(nmatrl),da_cal(nobs),da_obs(nobs)
    double precision CS(nmatrl),dtheta(nmatrl),CStemp(nmatrl)


    !call ZGEMM('N','C',nmatrl,nobs,nmatrl,dcmplx(1d0,0d0),dcmplx(sigma_M2,0d0),nmatrl,GA,nobs,dcmplx(0d0,0d0),Cmd,nmatrl)
    !$omp parallel do
    do itmp=1, nmatrl
        Cmd(itmp,:)=sigma_M2(itmp,itmp)*dconjg(Ga(:,itmp))
    enddo
    !$omp end parallel do
    call ZGEMM('N','N',nobs,nobs,nmatrl,dcmplx(1d0,0d0),GA,nobs,Cmd,nmatrl,dcmplx(0d0,0d0),K,nobs)


    !call ZGEMM('N','N',nobs,nmatrl,nmatrl,dcmplx(1d0,0d0),GA,nobs,dcmplx(sigma_M2,0d0),nmatrl,dcmplx(0d0,0d0),GaCOV,nobs)
    !call ZGEMM('N','C',nobs,nobs,nmatrl,dcmplx(1d0,0d0),GaCOV,nobs,GA,nobs,dcmplx(0d0,0d0),K,nobs)
    !$omp parallel do
    do itmp=1,nobs
        K(itmp,itmp)=K(itmp,itmp)+1/dcmplx(invcovUz(1),0d0)
    enddo
    !$omp end parallel do
    Ldd=K
    !call cholesky_lower_Z(nobs, K,Ldd)

    !call solution_ensemble(GA,Cmd,Lk,da_cal,da_obs,CS,CS0,nmatrl,nobs,dtheta)

    endsubroutine Cmd_Cdd_ensemble




    subroutine solution_ensemble(GA,Cmd,Ldd_1,CA_Uz_cal,CA_Uz_obs,CS,CS0,nmatrl,nmeas,dCS)
    implicit none
    integer nmatrl,nmeas,itmp
    double complex GA(nmeas,nmatrl),Cmd(nmatrl,nmeas),Ldd_1(nmeas,nmeas),CA_Uz_cal(nmeas),CA_Uz_obs(nmeas),D(nmatrl,nmeas),b(nmeas),gg(nmeas),dCSc(nmatrl),b1(nmeas),b2(nmeas),inv_Ldd_1(nmeas,nmeas)
    double precision CS(nmatrl),CS0(nmatrl),dCS(nmatrl)

    call ZGEMM('N','N',nmeas,1,nmatrl,dcmplx(1d0,0d0),GA,nmeas,dcmplx(CS-CS0,0d0),nmatrl,dcmplx(0d0,0d0),b,nmeas)
    b=CA_Uz_cal-CA_Uz_obs-b
    call solveAxB_Z(nmeas,1, Ldd_1, b,b1)
    call ZGEMM('N','N',nmatrl,1,nmeas,dcmplx(1d0,0d0),Cmd,nmatrl,b1,nmeas,dcmplx(0d0,0d0),dCsc,nmatrl)
    dCs=dreal(dCsc)

    !D=Cmd
    !call ZTRSM('R','U','N','N',nmatrl,nmeas,dcmplx(1d0,0d0),Ldd_1,nmeas,D,nmatrl)
    !call ZGEMM('N','N',nmeas,1,nmatrl,dcmplx(1d0,0d0),GA,nmeas,dcmplx(CS-CS0,0d0),nmatrl,dcmplx(0d0,0d0),b,nmeas)
    !b=CA_Uz_cal-CA_Uz_obs-b
    !call ZTRSM('L','U','C','N',nmeas,1,dcmplx(1d0,0d0),Ldd_1,nmeas,b,nmeas)
    !call ZGEMM('N','N',nmatrl,1,nmeas,dcmplx(1d0,0d0),D,nmatrl,b,nmeas,dcmplx(0d0,0d0),dCsc,nmatrl)
    !dCs=dreal(dCsc)
    endsubroutine solution_ensemble

    subroutine hessgrad(dUz_estimated,Uz_estimated,Uz_measured,invcovUz, inv_sigma_M2,CSest,CS0,nobs,nmatrl,hess,grad)
    implicit none
    integer nobs,nmatrl,iz,ix,itmp
    double complex Uz_measured(nobs),Uz_estimated(nobs), dUz_estimated(nobs,nmatrl),hess(nmatrl,nmatrl),grad(nmatrl)
    double complex,allocatable ::grad1(:)
    double precision inv_sigma_M2(nmatrl,nmatrl),invcovUz(nobs),CS0(nmatrl),CSest(nmatrl)
    double complex,allocatable ::dUzTinvCovUz(:,:),dUzconj(:,:),hess1(:,:),Uz_temp(:),dUzconjT(:,:)
    double precision,allocatable ::CStemp(:),grad2(:)
    integer(kind=8) :: tclock11, tclock22, clock_rate2
    real(kind=8) :: elapsed_time2
    integer(kind=8) :: tclock1, tclock2
    !allocate(dUzconj(nobs,nmatrl))
    !write (*,*) '2'
    !dUzconj=dconjg(dUz_estimated)
    !write (*,*) '2.1'
    !allocate(dUzconjT(nmatrl,nobs))
    !dUzconjT=transpose(dUzconj)
    !deallocate(dUzconj)
    allocate(dUzTinvCovUz(nmatrl,nobs))

    !dUzTinvCovUz=dUzconjT*dcmplx(invcovUz(1),0d0)
    !call system_clock(tclock2)
    !dUzTinvCovUz=transpose(dconjg(dUz_estimated))*dcmplx(invcovUz(1),0d0)
    dUzTinvCovUz=(dconjg(dUz_estimated))
    !call system_clock(tclock22, clock_rate2)
    !elapsed_time2 = float(tclock22 - tclock2) / float(clock_rate2)
    !write(*,*) "---transpose conjugate", elapsed_time2
    !deallocate(dUzconjT)

    allocate(hess1(nmatrl,nmatrl),Uz_temp(nobs),grad1(nmatrl))
    !call zsymm('R', 'U', nmatrl, nobs, dcmplx(1d0,0d0),dcmplx(invcovUz,0d0),nobs, dUzconjT, nmatrl, dcmplx(0d0,0d0), dUzTinvCovUz, nmatrl)
    !!call zgemm('T', 'N', nmatrl, nobs, nobs, dcmplx(1d0,0d0), dUzconj, nobs, dcmplx(invcovUz,0d0),nobs, dcmplx(0d0,0d0), dUzTinvCovUz, nmatrl)
    !call system_clock(tclock2)
    call zgemm('T', 'N', nmatrl, nmatrl, nobs, dcmplx(invcovUz(1),0d0), dUzTinvCovUz, nobs, dUz_estimated,nobs, dcmplx(0d0,0d0), hess1, nmatrl)
    !
    !call zgemm('N', 'N', nmatrl, nmatrl, nobs, dcmplx(1d0,0d0), dUzTinvCovUz, nmatrl, dUz_estimated,nobs, dcmplx(0d0,0d0), hess1, nmatrl)
    hess=hess1

    hess=hess+inv_sigma_M2

    !call system_clock(tclock22, clock_rate2)
    !    elapsed_time2 = float(tclock22 - tclock2) / float(clock_rate2)
    !    write(*,*) "---hess", elapsed_time2
    !hess=hess+inv_sigma_M2
    !call zsyrk('L','T',nmatrl,nobs,dcmplx(1d0,0d0),dUzconj,nobs,dcmplx(0d0,0d0),hess1, nmatrl)
    !hess=0d0
    !!$omp parallel do
    !    do iz=1,nmatrl
    !
    !        hess(iz,iz)=dreal(hess1(iz,iz))
    !
    !    enddo
    !    !$omp end parallel do
    !    hess=hess+invcovM
    !hess=dreal(hess1)*invcovUz(1,1)+invcovM

    !call system_clock(tclock2)
    Uz_temp=Uz_estimated- Uz_measured
    call zgemm('T', 'N', nmatrl, 1,nobs,  dcmplx(invcovUz(1),0d0), dUzTinvCovUz, nobs, Uz_temp,nobs, dcmplx(0d0,0d0), grad1, nmatrl)
    !call zgemm('N', 'N', nmatrl, 1,nobs,  dcmplx(1d0,0d0), dUzTinvCovUz, nmatrl, Uz_temp,nobs, dcmplx(0d0,0d0), grad1, nmatrl)
    !call system_clock(tclock22, clock_rate2)
    !    elapsed_time2 = float(tclock22 - tclock2) / float(clock_rate2)
    !    write(*,*) "---grad1", elapsed_time2

    !call system_clock(tclock2)
    allocate(CStemp(nmatrl),grad2(nmatrl))
    CStemp=CSest - CS0
    call dsymm('L', 'U', nmatrl, 1, 1d0,inv_sigma_M2,nmatrl, CStemp, nmatrl, 0d0, grad2, nmatrl)
    !do itmp=1,nmatrl
    !    grad2( itmp)=inv_sigma_M2( itmp)*CStemp(itmp)
    !enddo

    !call system_clock(tclock22, clock_rate2)
    !elapsed_time2 = float(tclock22 - tclock2) / float(clock_rate2)
    !write(*,*) "---grad2", elapsed_time2
    !grad2=matmul(inv_sigma_M2,CStemp)

    !call system_clock(tclock2)
    !call dgemm('N', 'N', nmatrl,1, nmatrl, 1d0, invcovM, nmatrl, CStemp,nmatrl, 0d0, grad2, nmatrl)
    grad=grad1+grad2
    !call system_clock(tclock22, clock_rate2)
    !   elapsed_time2 = float(tclock22 - tclock2) / float(clock_rate2)
    !   write(*,*) "---grad", elapsed_time2

    deallocate(dUzTinvCovUz,hess1,Uz_temp)
    deallocate(CStemp,grad2)

    endsubroutine hessgrad

    subroutine postprop(hess,CStemp,coeff,nmatrl,logprop)
    implicit none
    integer nmatrl
    double precision hess(nmatrl,nmatrl),CStemp(nmatrl),CStemp1(nmatrl)
    double precision logprop,coeff

    call dsymm('L', 'U', nmatrl, 1, 1d0,hess,nmatrl, CStemp, nmatrl, 0d0, CStemp1, nmatrl)

    call dgemm('N', 'N', nmatrl,1, nmatrl, 1d0, hess, nmatrl, CStemp,nmatrl, 0d0, CStemp1, nmatrl)

    logprop = -0.5d0*coeff*dot_product(CStemp1,CStemp)

    endsubroutine postprop

    subroutine random_stduniform(u)
    implicit none
    double precision,intent(out) :: u
    double precision :: r
    call random_number(r)
    u = 1 - r
    endsubroutine random_stduniform

    subroutine random_stdnormal(x,nx)
    implicit none
    double precision,intent(out) :: x(nx)
    double precision :: pi
    double precision :: u1,u2
    integer nx,ix
    pi = 4.d0 * datan(1.d0)
    !$omp parallel do
    do ix=1,nx
        call random_stduniform(u1)
        call random_stduniform(u2)
        x(ix) = sqrt(-2*log(u1))*cos(2*pi*u2)
    enddo
    !$omp end parallel do
    endsubroutine random_stdnormal





    subroutine create_CSR_formartZ(A,xA,yA,ia,ja,vA,nnzA)
    implicit none


    integer xA,yA,nnzA,ia(xA+1),ja(nnzA),t,ix,iz,nnz
    double complex A(xA,yA)
    double complex vA(nnzA)

    t=1
    do ix=1,xA
        do iz=1,yA
            if (A(ix,iz)/=dcmplx(0d0,0d0)) then
                ja(t)=iz
                vA(t)=A(ix,iz)
                t=t+1;
            endif
        enddo
    enddo
    ia(1)=1;

    do ix=1,xA
        nnz=count(A(ix,:)/=dcmplx(0d0,0d0))
        ia(ix+1)=ia(ix)+nnz
    enddo
    end subroutine create_CSR_formartZ

    subroutine As_m_Bd(v_A,c_A,pB_A,pE_A,nnz_A,ndof_u,nelmt,Ne,B,C)
    implicit none
    integer ndof_u,nelmt,i,j,k,nnz_A,Ne
    integer pB_A(ndof_u),pE_A(ndof_u),c_A(nnz_A)
    double complex B(nelmt,Ne),v_A(nnz_A),C(ndof_u,Ne)

    C=dcmplx(0d0,0d0)
    !$omp parallel do
    do j = 1, ndof_u
        do k = pB_A(j), pE_A(j)-1
            do i = 1, Ne
                C(j,i) = C(j,i) + B(c_A(k),i) * v_A(k)
            end do
        end do
    end do
    !$omp end parallel do
    end subroutine As_m_Bd


    subroutine AdT_m_Bs(v_B,c_B,pB_B,pE_B,nnz_B,ndof_u,nmatrl,nout,A,C)
    implicit none
    integer ndof_u,nmatrl,i,j,k,nnz_B,nout
    integer pB_B(ndof_u),pE_B(ndof_u),c_B(nnz_B)
    double complex A(ndof_u,nout),v_B(nnz_B),C(nout,nmatrl)

    !C=dcmplx(0d0,0d0)
    !!$omp parallel do
    do j = 1, nout
        do k = pB_B(j), pE_B(j)-1
            do i = 1, nmatrl
                C(j,i) = C(j,i) + A(c_B(k),j) * v_B(k)
            end do
        end do
    end do
    !!$omp end parallel do
    end subroutine AdT_m_Bs


    subroutine sparse_get_valueZ(nrow,ncol,c,pB,pE,v,nnz,A)
    integer i,j,k, nrow,ncol,pB(nrow),pE(nrow),c(nnz)
    double complex A(nrow,ncol),v(nnz)
    !!$omp parallel do
    do j = 1, nrow
        do k = pB(j), pE(j)-1
            !do i = 1, ncol
            v(k)=A(j,c(k))
            !end do
        end do
    end do
    !!$omp end parallel do

    end subroutine sparse_get_valueZ


    subroutine sparse_get_valueTZ(nshot,nrow,ncol,c,pB,pE,v,nnz,A)
    implicit none
    integer i,j,k,nnz, nrow,ncol,pB(nrow),pE(nrow),c(nnz),nshot
    double complex A(ncol,nshot,nrow),v(nnz,nshot)
    !!$omp parallel do
    do j = 1, nrow
        do k = pB(j), pE(j)-1
            !do i = 1, ncol
            v(k,:)=A(c(k),:,j)
            !end do
        end do
    end do
    !!$omp end parallel do

    end subroutine sparse_get_valueTZ

    subroutine PDML_parameters(nz,nx,npml,nmatrl,pi,width,depth,nnode,nelmt,CS_equi,CP_equi,enode,matrl,etype,epara,x,z,dx,dz,nuy,CS)
    implicit none
    double precision cs_left,cs_right,cs_bottom,dx,dz,nuy, pi,width,depth,c1
    integer nz,nx,npml,imatrl,nmatrl,nelmt,itmp,jtmp,ielmt,inode,nnode

    integer enode(nelmt,0:9), matrl(nelmt), etype(nelmt)
    double complex epara(nelmt,2)
    double precision x(nnode), z(nnode),CS(nmatrl),CSsr(nmatrl)
    double complex CS_equi(nelmt),CP_equi(nelmt)
    do imatrl=1,nmatrl
        CSsr(imatrl)=dsqrt(CS(imatrl))
    enddo


    c1=0d0
    do imatrl=1,(nx-2*npml)
        c1=c1+CSsr(imatrl)
    enddo
    do imatrl=1,((nz-npml)-1)
        c1=c1+CSsr(imatrl*(nx-2*npml)+1)
    enddo
    do imatrl=1,((nz-npml)-1)
        c1=c1+CSsr(imatrl*(nx-2*npml)+(nx-2*npml))
    enddo
    c1=c1/((nx-2*npml)+2*((nz-npml)-1))
    cs_bottom=c1
    cs_left=c1
    cs_right=c1



    !cs_bottom=0d0
    !do imatrl=1,(nx-2*npml)
    !    cs_bottom=cs_bottom+CSsr(imatrl)
    !enddo
    !cs_bottom=cs_bottom/(nx-2*npml)
    !
    !cs_left=0d0
    !do imatrl=1,((nz-npml)-1)
    !    cs_left=cs_left+CSsr(imatrl*(nx-2*npml)+1)
    !enddo
    !cs_left=cs_left/((nz-npml)-1)
    !
    !cs_right=0d0
    !do imatrl=1,((nz-npml)-1)
    !    cs_right=cs_right+CSsr(imatrl*(nx-2*npml)+(nx-2*npml))
    !enddo
    !cs_right=cs_right/((nz-npml)-1)
    !
    CS_equi=dcmplx(0d0,0d0)
    CP_equi=dcmplx(0d0,0d0)

    write (*,*)

    do jtmp = 1 , nz + 1

        do itmp	= 1 , nx + 1

            inode = (nx + 1) * (jtmp - 1) + itmp
            x(inode) = (itmp - 1) * dx
            z(inode) = (jtmp - 1) * dz - npml*2 * dz

        enddo

    enddo

    do jtmp = 1 , nz

        do itmp = 1 , nx

            ielmt = nx * (jtmp - 1) + itmp

            enode(ielmt,0) = 4
            enode(ielmt,1) = (nx + 1) * (jtmp - 1) + itmp
            enode(ielmt,2) = (nx + 1) * (jtmp - 1) + itmp + 1
            enode(ielmt,3) = (nx + 1) *  jtmp      + itmp + 1
            enode(ielmt,4) = (nx + 1) *  jtmp      + itmp

            matrl(ielmt) = (nx - 2*npml) * (jtmp - npml - 1) + itmp-npml

            etype(ielmt) = 0
            epara(ielmt,1) = 0.d0
            epara(ielmt,2) = 0.d0

            if ((itmp >= npml+1) .and. (itmp <= nx-npml) .and. (jtmp <= npml-3)) then

                matrl(ielmt) = 0
                epara(ielmt,1) = 0.d0
                epara(ielmt,2) = (2*(npml-3+1-jtmp)-1) * pi /  width
                etype(ielmt) = 22
                CS_equi(ielmt)=dcmplx(cs_bottom,0d0)
            else if ((itmp >= nx-npml+1) .and. (itmp <= nx-3-1) .and. (jtmp <= npml-3)) then

                matrl(ielmt) = 0
                if (itmp==nx-npml+1) then
                    epara(ielmt,1) = dcmplx(cs_right,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_right,0d0)

                else if (itmp==nx-npml+2) then
                    epara(ielmt,1) = dcmplx(cs_right*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_right*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (itmp==nx-npml+3) then
                    epara(ielmt,1) = 1000d0
                endif
                CS_equi(ielmt)=dcmplx(cs_right,0d0)
                epara(ielmt,2) = (2*(npml-3+1-jtmp)-1) * pi /  width
                etype(ielmt) = 32

            else if ((itmp >= npml-3+1) .and. (itmp <= npml) .and. (jtmp <= npml-3)) then

                matrl(ielmt) = 0
                if (itmp==npml-3+3) then
                    epara(ielmt,1) = dcmplx(cs_left,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_left,0d0)
                else if (itmp==npml-3+2) then
                    epara(ielmt,1) = dcmplx(cs_left*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_left*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)

                else if (itmp==npml-3+1) then
                    epara(ielmt,1) = 1000d0
                endif
                CS_equi(ielmt)=dcmplx(cs_left,0d0)
                epara(ielmt,2) = (2*(npml-3+1-jtmp)-1) * pi /  width
                etype(ielmt) = 32

            else if ((itmp >=1) .and. (itmp <= npml-3)  .and. (jtmp <= npml-3)) then

                matrl(ielmt) = 0
                epara(ielmt,1) = (2*(npml-3+1-itmp)-1) * pi / (2*depth)
                epara(ielmt,2) = (2*(npml-3+1-jtmp)-1) * pi /  width
                etype(ielmt) = 34
                CS_equi(ielmt)=0.5d0*(dcmplx(cs_left,0d0)+dcmplx(cs_bottom,0d0))
            else if ((itmp >=nx-npml+3) .and. (jtmp <= npml-3))then

                matrl(ielmt) = 0
                epara(ielmt,1) = (2*(itmp-nx+npml-3)-1) * pi / (2*depth)
                epara(ielmt,2) = (2*(npml-3+1-jtmp)-1) * pi /  width
                etype(ielmt) = 34
                CS_equi(ielmt)=0.5d0*(dcmplx(cs_right,0d0)+dcmplx(cs_bottom,0d0))
            else if ((itmp >= npml+1) .and. (itmp <= nx-npml) .and. (jtmp > npml-3) .and. (jtmp <= npml)) then

                matrl(ielmt) = 0
                epara(ielmt,1) = 0d0
                if (jtmp==npml-3+3) then
                    epara(ielmt,2) =dcmplx(cs_bottom,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_bottom,0d0)
                else if (jtmp==npml-3+2) then
                    epara(ielmt,2) =  dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (jtmp==npml-3+1) then
                    epara(ielmt,2) = 1000d0
                endif
                CS_equi(ielmt)=dcmplx(cs_bottom,0d0)
                etype(ielmt) = 21

            else if ((itmp >= nx-npml+1) .and. (itmp <= nx-3-1).and. (jtmp > npml-3) .and. (jtmp <= npml)) then

                matrl(ielmt) =0
                if (itmp==nx-npml+1) then
                    epara(ielmt,1) = dcmplx(cs_right,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_right,0d0)
                else if (itmp==nx-npml+2) then
                    epara(ielmt,1) =dcmplx(cs_right*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_right*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (itmp==nx-npml+3) then
                    epara(ielmt,1) = 1000d0
                endif
                if (jtmp==npml-3+3) then
                    epara(ielmt,2) =dcmplx(cs_bottom,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_bottom,0d0)
                else if (jtmp==npml-3+2) then
                    epara(ielmt,2) =  dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (jtmp==npml-3+1) then
                    epara(ielmt,2) = 1000d0
                endif
                CS_equi(ielmt)=0.5d0*(dcmplx(cs_right,0d0)+dcmplx(cs_bottom,0d0))
                etype(ielmt) = 31

            else if ((itmp >= npml-3+1) .and. (itmp <= npml) .and. (jtmp > npml-3) .and. (jtmp <= npml)) then

                matrl(ielmt) = 0
                if (itmp==npml-3+3) then
                    epara(ielmt,1) = dcmplx(cs_left,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_left,0d0)
                else if (itmp==npml-3+2) then
                    epara(ielmt,1) = dcmplx(cs_left*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_left*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (itmp==npml-3+1) then
                    epara(ielmt,1) = 1000d0
                endif
                if (jtmp==npml-3+3) then
                    epara(ielmt,2) =dcmplx(cs_bottom,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_bottom,0d0)
                else if (jtmp==npml-3+2) then
                    epara(ielmt,2) =  dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (jtmp==npml-3+1) then
                    epara(ielmt,2) = 1000d0
                endif
                CS_equi(ielmt)=0.5d0*(dcmplx(cs_left,0d0)+dcmplx(cs_bottom,0d0))
                etype(ielmt) = 31

            else if ((itmp >=1) .and. (itmp <= npml-3) .and. (jtmp > npml-3) .and. (jtmp <= npml)) then

                matrl(ielmt) = 0
                epara(ielmt,1) = (2*(npml-3+1-itmp)-1) * pi / (2*depth)
                if (jtmp==npml-3+3) then
                    epara(ielmt,2) =dcmplx(cs_bottom,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_bottom,0d0)
                else if (jtmp==npml-3+2) then
                    epara(ielmt,2) =  dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (jtmp==npml-3+1) then
                    epara(ielmt,2) = 1000d0
                endif
                CS_equi(ielmt)=dcmplx(cs_bottom,0d0)
                etype(ielmt) = 33

            else if ((itmp >=nx-npml+3) .and. (jtmp > npml-3) .and. (jtmp <= npml)) then

                matrl(ielmt) =0
                epara(ielmt,1) = (2*(itmp-nx+npml-3)-1) * pi / (2*depth)
                if (jtmp==npml-3+3) then
                    epara(ielmt,2) =dcmplx(cs_bottom,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_bottom,0d0)
                else if (jtmp==npml-3+2) then
                    epara(ielmt,2) =  dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_bottom*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (jtmp==npml-3+1) then
                    epara(ielmt,2) = 1000d0
                endif
                etype(ielmt) = 33
                CS_equi(ielmt)=dcmplx(cs_bottom,0d0)
            else if ((itmp >= npml-3+1) .and. (itmp <= npml)  .and. (jtmp > npml) ) then

                matrl(ielmt) = 0
                if (itmp==npml-3+3) then
                    epara(ielmt,1) = dcmplx(cs_left,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_left,0d0)
                else if (itmp==npml-3+2) then
                    epara(ielmt,1) = dcmplx(cs_left*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_left*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (itmp==npml-3+1) then
                    epara(ielmt,1) = 1000d0
                endif
                CS_equi(ielmt)=dcmplx(cs_left,0d0)
                epara(ielmt,2) = 0.d0
                etype(ielmt) = 11

            else if ((itmp >= nx-npml+1) .and. (itmp <= nx-3-1).and. (jtmp > npml-3) .and. (jtmp > npml)) then

                matrl(ielmt) =0
                if (itmp==nx-npml+1) then
                    epara(ielmt,1)  = dcmplx(cs_right,0d0)
                    !CS_equi(ielmt)=dcmplx(cs_right,0d0)
                else if (itmp==nx-npml+2) then
                    epara(ielmt,1) = dcmplx(cs_right*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                    !CP_equi(ielmt)=dcmplx(cs_right*dsqrt((1d0-nuy)/(0.5d0-nuy)),0d0)
                else if (itmp==nx-npml+3) then
                    epara(ielmt,1) = 1000d0
                endif
                CS_equi(ielmt)=dcmplx(cs_right,0d0)
                epara(ielmt,2) = 0.d0
                etype(ielmt) = 11

            else if ((itmp >=nx-npml+3) .and. (jtmp > npml) ) then

                matrl(ielmt) = 0
                epara(ielmt,1) = (2*(itmp-nx+npml-3)-1) * pi / (2*depth)
                epara(ielmt,2) = 0.d0
                etype(ielmt) = 12
                CS_equi(ielmt)=dcmplx(cs_right,0d0)
            else if ((itmp >=1) .and. (itmp <= npml-3) .and. (jtmp > npml)) then

                matrl(ielmt) = 0
                epara(ielmt,1) = (2*(npml-3+1-itmp)-1) * pi / (2*depth)
                epara(ielmt,2) = 0.d0
                etype(ielmt) = 12
                CS_equi(ielmt)=dcmplx(cs_left,0d0)
            endif

        enddo

    enddo
    end subroutine PDML_parameters

    subroutine update_inv_sigma_M2(dUz_estimated,invcovUz, inv_sigma_M2,nobs,nmatrl)
    implicit none
    integer nobs,nmatrl,iz,ix,itmp,imatrl
    double complex Uz_measured(nobs),Uz_estimated(nobs), dUz_estimated(nobs,nmatrl)
    double precision inv_sigma_M2(nmatrl,nmatrl),invcovUz(nobs),CSest(nmatrl),hess22(nmatrl,nmatrl)
    double complex,allocatable ::dUzTinvCovUz(:,:),dUzconj(:,:),hess1(:,:),hess2(:,:),hess(:,:),Uz_temp(:),dUzconjT(:,:),inv_hess(:,:),inv_hess_hess1(:,:)


    !allocate(dUzTinvCovUz(nmatrl,nobs))
    !dUzTinvCovUz=dconjg(transpose(dUz_estimated))*dcmplx(invcovUz(1),0d0)


    allocate(hess1(nmatrl,nmatrl),hess(nmatrl,nmatrl),Uz_temp(nobs),inv_hess_hess1(nmatrl,nmatrl),hess2(nmatrl,nmatrl),inv_hess(nmatrl,nmatrl))
    !call zgemm('T', 'N', nmatrl, nmatrl, nobs, dcmplx(invcovUz(1),0d0), dUzTinvCovUz, nobs, dUz_estimated,nobs, dcmplx(0d0,0d0), hess1, nmatrl)
    !
    !call zgemm('N', 'N', nmatrl, nmatrl, nobs, dcmplx(1d0,0d0), dUzTinvCovUz, nmatrl, dUz_estimated,nobs, dcmplx(0d0,0d0), hess1, nmatrl)
    call zgemm('C', 'N', nmatrl, nmatrl, nobs, dcmplx(1d0,0d0), dUz_estimated, nobs, dUz_estimated,nobs, dcmplx(0d0,0d0), hess1, nmatrl)
    hess1=hess1*dcmplx(invcovUz(1),0d0)
    hess=hess1+inv_sigma_M2


    write(*,*) '1'
    call inversematrixz(nmatrl, hess, inv_hess)
    write(*,*) '2'

    inv_hess_hess1=matmul(inv_hess,hess1)
    
    inv_hess=transpose((inv_hess))
    hess2=matmul(inv_hess_hess1,inv_hess)
    inv_hess_hess1=-inv_hess_hess1
!    call ZGEMM('N','N',nmatrl,nmatrl,nmatrl,dcmplx(1d0,0d0),inv_hess,nmatrl,hess1,nmatrl,dcmplx(0d0,0d0),inv_hess_hess1,nmatrl)
!inv_hess=transpose((inv_hess))
!    call ZGEMM('N','N',nmatrl,nmatrl,nmatrl,dcmplx(1d0,0d0),inv_hess_hess1,nmatrl,inv_hess,nmatrl,dcmplx(0d0,0d0),hess2,nmatrl)
!inv_hess_hess1=-inv_hess_hess1

    write(*,*) '3'
    do itmp=1,nmatrl
        inv_hess_hess1(itmp,itmp)=1d0+inv_hess_hess1(itmp,itmp)
    enddo
    do itmp=1,nmatrl
        inv_sigma_m2(itmp,itmp)=2d0/inv_sigma_m2(itmp,itmp)
    enddo

    hess=matmul(inv_hess_hess1,inv_sigma_m2)
        !call DGEMM('N','N',nmatrl,nmatrl,nmatrl,1d0,dreal(inv_hess_hess1),nmatrl,inv_sigma_m2,nmatrl,0d0,hess,nmatrl)

    write(*,*) '4'
    inv_hess_hess1=transpose((inv_hess_hess1))
    hess=matmul(hess,inv_hess_hess1)
    !hess22=dreal(inv_hess_hess1)
    ! call DGEMM('N','N',nmatrl,nmatrl,nmatrl,1d0,hess,nmatrl,hess22,nmatrl,0d0,hess,nmatrl)

    write(*,*) '5'
    hess=hess+hess2
    

    inv_sigma_M2=dreal(hess)

    !!$omp parallel do
    !do imatrl=1,nmatrl
    !    call zgemm('N', 'N', 1, 1, nobs, dcmplx(1d0,0d0), dUzTinvCovUz(imatrl,:), 1, dUz_estimated(:,imatrl),nobs, dcmplx(0d0,0d0), hess1(imatrl),1)
    !enddo
    !!$omp end parallel do
    !inv_sigma_M2=dreal(hess1)+inv_sigma_M2


    deallocate(hess1,Uz_temp)


    endsubroutine update_inv_sigma_M2


    subroutine update_inv_sigma_M2_rev(dUz_estimated,invcovUz, sigma_M2,nobs,nmatrl,inv_sigma_M2)
    implicit none
    integer nobs,nmatrl,iz,ix,itmp,imatrl
    double complex Uz_measured(nobs),Uz_estimated(nobs), dUz_estimated(nobs,nmatrl),Cmd(nmatrl,nobs),F(nmatrl,nobs),K(nobs,nobs),inv_K(nobs,nobs)
    double precision inv_sigma_M2(nmatrl,nmatrl),sigma_M2(nmatrl,nmatrl),invcovUz(nobs),CSest(nmatrl)
    double precision,allocatable :: h1(:,:),h2(:,:),FFr(:,:),FFsigma(:,:)
    double complex,allocatable ::dUzTinvCovUz(:,:),dUzconj(:,:),FF(:,:),h11(:,:),h22(:,:)
    allocate(FF(nmatrl,nmatrl),FFr(nmatrl,nmatrl),FFsigma(nmatrl,nmatrl),h1(nmatrl,nmatrl),h11(nmatrl,nmatrl),h2(nmatrl,nmatrl),h22(nmatrl,nmatrl))


    write(*,*) '1'
    !$omp parallel do
    do itmp=1, nmatrl
        Cmd(itmp,:)=sigma_M2(itmp,itmp)*dconjg(dUz_estimated(:,itmp))
    enddo
    !$omp end parallel do
    call ZGEMM('N','N',nobs,nobs,nmatrl,dcmplx(1d0,0d0),dUz_estimated,nobs,Cmd,nmatrl,dcmplx(0d0,0d0),K,nobs)

    write(*,*) '2'
    !$omp parallel do
    do itmp=1,nobs
        K(itmp,itmp)=K(itmp,itmp)+1/dcmplx(invcovUz(1),0d0)
    enddo
    !$omp end parallel do

    write(*,*) '3'
    call inversematrixz(nobs, K, inv_K)
    !F=matmul(Cmd,inv_K)
    call ZGEMM('N','N',nmatrl,nobs,nobs,dcmplx(1d0,0d0),Cmd,nmatrl,inv_K,nobs,dcmplx(0d0,0d0),F,nmatrl)
    write(*,*) '4'
    !FF=matmul(F,dUz_estimated)
    call ZGEMM('N','N',nmatrl,nmatrl,nobs,dcmplx(1d0,0d0),F,nmatrl,dUz_estimated,nobs,dcmplx(0d0,0d0),FF,nmatrl)

    write(*,*) '5'
    do itmp=1,nmatrl
        FF(itmp,itmp)=FF(itmp,itmp)-dcmplx(1d0,0d0)
    enddo
    FF=-FF
    write(*,*) '4'
    FFr=dreal(FF)
    
    !h1=matmul(FFr,sigma_M2)
    !FFr=transpose(FFr)
    !h1=2d0*matmul(h1,FFr)

    call DGEMM('N','N',nmatrl,nmatrl,nmatrl,1d0,FFr,nmatrl,sigma_M2,nmatrl,0d0,FFsigma,nmatrl)
    call DGEMM('N','T',nmatrl,nmatrl,nmatrl,1d0,FFsigma,nmatrl,FFr,nmatrl,0d0,h1,nmatrl)    
    h1=2d0*h1
    
    write(*,*) '5'
    !h2=dreal(1/dcmplx(invcovUz(1),0d0)*matmul(F,transpose(dconjg(F))))

    call ZGEMM('N','C',nmatrl,nmatrl,nobs,dcmplx(1d0,0d0),F,nmatrl,F,nmatrl,dcmplx(0d0,0d0),h22,nmatrl)
    h2=(1/invcovUz(1))*dreal(h22)

    write(*,*) '6'

    inv_sigma_M2=(h1+h2)

    endsubroutine update_inv_sigma_M2_rev