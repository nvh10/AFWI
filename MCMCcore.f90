    subroutine MCMCcore(group_frequency,CS_equi,nshot,nx, nz, npml, dx, dz, nnode, nelmt, nmatrl, CS, rhos,nuy,beta, x, z, enode, matrl, etype, epara, nbc, bc, id, &
        nfreq_total, freq_total,nfreq,ngroup,width,depth, Pext, inpar, inobs,outgama, outpost1, outpost0,outprop01,outprop10, outest, outrsp, outrsl,outcovM, nobs, onode, CSest, sigma_M2, covUz, nSample,alpha_coef, beta_coef)
    implicit none
    integer outgama, outpost1, outpost0,outprop01,outprop10,nshot,igroup

    integer  nfreq_total,itype,ishot,nfreq,ngroup,dd
    integer nnode, nelmt, nmatrl, nSample, norder, nobs, onode(nobs), npml
    double precision alpha_coef, beta_coef, lb(nmatrl),ub(nmatrl),u_rand,dtmp2,dtmp1,width,depth
    real(kind=16)::detHa,detH,detHb
    double precision rhos, CS(nmatrl), CSest(nmatrl,2),  CS0(nmatrl),nuy,beta,dCs(nmatrl),CStemp(nmatrl)

    integer*8 pt(64)
    double precision,allocatable:: cholhess(:,:),invHess(:,:),hess1(:,:),grad1(:),inv_sigma_M2(:,:),hess_real(:,:),grad_real(:)
    integer INFO
    double precision sigma_M2(nmatrl,nmatrl),covUz(2*nobs*nfreq*nshot), invcovUz(2*nobs*nfreq*nshot)
    double precision prop01,prop10,post0,post1,gama
    double precision Uz_measuredt(2*nobs*nfreq_total*nshot),  derr(nobs*nfreq*nshot), randCs(nmatrl),derr_CS(nmatrl)
    double complex Uz_measured(nobs*nfreq_total*nshot),a1,a2,a4,CS_equi(nelmt),CP_equi(nelmt),dStiffness_all(8,8,nmatrl),dCs_c(nmatrl)
    double complex,allocatable ::hess(:,:),grad(:),Hmat(:,:),dprd(:),dobs(:)


    ! node data
    double precision x(nnode), z(nnode)

    ! element data
    integer enode(nelmt,0:9), matrl(nelmt), etype(nelmt)
    double complex epara(nelmt,2)

    ! boundary condition data
    integer nbc(3), bc(nnode), id(nnode,2), ndof, ndof2

    double precision Pext(nnode*2,nshot) !
    ! frequency data
    double precision  freq_total(nfreq_total),freq(nfreq),group_frequency(nfreq,ngroup)
    integer outsgb, outrhb, outest, inexa, inpar, inrfa, outerr, inobs, outrsp, outrsl,outcovM,nnz_P2
    double precision pi, ww    !
    double complex,allocatable :: MSN22(:), DPN22(:), STN22(:), RRN22(:),DSN22(:)
    double complex,allocatable :: eK(:,:),v_P2(:,:),C_md(:,:),Ldd(:,:)
    integer,allocatable :: ia_P2(:),ja_p2(:),pB_P2(:),pE_P2(:)
    double complex,allocatable :: Uz_estimated(:,:), dUz_estimated(:,:,:)
    integer nn, n1, n2, nfix, ne, nm, ifreq, itmp, jtmp, inode, jnode, idof, jdof, ielmt, iSample, imatrl, ktmp
    integer,allocatable :: ia(:), ja(:), idcheck(:)
    integer node,  iparm(64), idum, ierr, iia
    integer nx, nz, ltmp, nx1, otime(10), iflag,iflag_group
    double precision dx, dz, x1, z1, outdata(5,nnode)
    integer(kind=8) :: tclock1, tclock2, clock_rate
    real(kind=8) :: elapsed_time
    integer(kind=8) :: tclock11, tclock22, clock_rate2
    real(kind=8) :: elapsed_time2
    double complex,allocatable :: invDSN22T(:,:),Dis2(:,:),GA(:,:),CA_Uz_cal(:),CA_Uz_obs(:)
    double complex,allocatable ::  dDis2dm(:,:), P2(:,:,:),dDis2dmT(:,:), P2T(:,:)
    double precision,allocatable ::  inv_Ca_D(:),Ca_D(:)
    pi = 4.d0 * datan(1.d0)
    allocate(dprd(nobs*nfreq*nshot),dobs(nobs*nfreq*nshot))
    allocate(Uz_estimated(nfreq*nobs*nshot,2),dUz_estimated (nfreq*nobs*nshot,nmatrl,1))
    allocate(hess(nmatrl,nmatrl),grad(nmatrl),inv_sigma_M2(nmatrl,nmatrl))
    allocate (GA(2*nfreq*nobs*nshot,nmatrl),CA_Uz_cal(2*nfreq*nobs*nshot),CA_Uz_obs(2*nfreq*nobs*nshot),inv_Ca_D(nfreq*nobs*nshot*2),Ca_D(nfreq*nobs*nshot*2))
    allocate(C_md(nmatrl,2*nfreq*nobs*nshot),Ldd(nfreq*nobs*nshot*2,nfreq*nobs*nshot*2))
allocate(hess_real(nmatrl,nmatrl),grad_real(nmatrl))
    !allocate(cholhess(nmatrl,nmatrl),hess_ap(nmatrl,nmatrl))
    !allocate(hess1(nmatrl,nmatrl),grad1(nmatrl))



    !call InverseMatrixD(nmatrl, covM, invcovM)
    inv_sigma_M2=0d0
    !call InverseMatrixD(nmatrl, sigma_M2, inv_sigma_M2)
    do imatrl=1,nmatrl
        inv_sigma_M2(imatrl,imatrl)=1d0/sigma_M2(imatrl,imatrl)
    enddo

    !call InverseMatrixD(nfreq*nobs*nshot, covUz, invcovUz)
    invcovUz=0d0
    do itmp=1, 2*nfreq*nobs*nshot
        invcovUz(itmp)=1d0/covUz(itmp)
    enddo
    CS0=CSest(:,1)
    allocate (ia(nbc(2)+1))
    inv_Ca_D=1d0/covUz(1)
    Ca_D=1d0*covUz(1)
    do itmp = 1 , nbc(2)+1
        read (inpar,*) ia(itmp)
    enddo
    read (inpar,*)

    allocate (ja(ia(nbc(2)+1)-1))
    allocate ( MSN22(ia(nbc(2)+1)-1), DPN22(ia(nbc(2)+1)-1), STN22(ia(nbc(2)+1)-1), RRN22(ia(nbc(2)+1)-1))
    allocate (DSN22(ia(nbc(2)+1)-1))



    do itmp = 1 , ia(nbc(2)+1)-1
        read (inpar,*) ja(itmp)
    enddo

    call pardisoinit(pt, 13, iparm)

    !write (outest,'(10x,<nmatrl>i15)') (itmp, itmp = 1 , nmatrl)
    write (outest,'(i10,i10,<nmatrl>f15.5)') 0,0, (CSest(1:nmatrl,1))
    write (outrsp,'(10x,<2*nobs>i16)') (onode(inode), onode(inode), inode = 1 , nobs)

    Uz_measuredt = 0.d0
    Uz_measured = 0.d0
    read (inobs,*)

    do ifreq = 1 , nfreq_total*nshot

        read (inobs,*) dtmp1, (Uz_measuredt(2*nobs*(ifreq-1)+2*itmp-1), Uz_measuredt(2*nobs*(ifreq-1)+2*itmp), itmp = 1 , nobs)

    enddo
    do ifreq = 1 , nfreq_total*nshot
        do itmp = 1 , nobs
            Uz_measured(nobs*(ifreq-1)+itmp)=dcmplx(Uz_measuredt(2*nobs*(ifreq-1)+2*itmp-1),Uz_measuredt(2*nobs*(ifreq-1)+2*itmp))
        enddo
    enddo

    iSample = 0
    iflag = 1

    allocate(eK(nbc(2),nobs))
    eK=dcmplx(0d0,0d0)

    do itmp=1,nobs

        eK(id(onode(itmp),2),itmp)=dcmplx(1d0,0d0)
    enddo
    allocate( Dis2(nbc(2),1))
    allocate(P2(nbc(2),1,nmatrl))
    !$omp parallel do
    do imatrl = 1 , nmatrl
        call dSdmU_elastic_pre(nx,nz,npml,1,onode,nobs,nnode, nelmt, nmatrl, CSest(:,1), rhos,nuy,beta, x, z, enode, matrl, etype, epara, nbc, id, ia, ja, imatrl,ww, Dis2, P2(:,:,imatrl),dStiffness_all)
    enddo
    !$omp end parallel do


    nnz_P2=count(P2(:,1,:)/=dcmplx(0d0,0d0))
    !allocate(v_P2(nnz_P2),ja_P2(nnz_P2),ia_P2(nbc(2)+1),pB_P2(nbc(2)),pE_P2(nbc(2)))
    allocate(v_P2(nnz_P2,nshot),ja_P2(nnz_P2),ia_P2(nmatrl+1),pB_P2(nmatrl),pE_P2(nmatrl))
    allocate(P2T(nmatrl,nbc(2)))
    P2T=transpose(P2(:,1,:))
    call create_CSR_formartZ(P2T,nmatrl,nbc(2),ia_P2,ja_P2,v_P2(:,1),nnz_P2)
    pB_P2=ia_P2(1:nmatrl)
    pE_P2=ia_P2(2:nmatrl+1)
    ! pB_P2=ia_P2(1:nbc(2))
    !pE_P2=ia_P2(2:nbc(2)+1)
    deallocate(P2T)
    deallocate( Dis2,P2)




    dd=(nfreq_total-nfreq)/(ngroup-2)
    igroup=1
    freq=group_frequency(:,igroup)
    !dobs=Uz_measured(1:nfreq*nshot*nobs)


    call FiniteElement_PlaneStrain(CS_equi,nnode, nelmt, nmatrl, CS, rhos,nuy,beta, MSN22, DPN22,  STN22, RRN22,x, z, enode, matrl, etype, epara, nbc, id, ia, ja)


    do ifreq = 1 , nfreq
        allocate(P2(nbc(2),nshot,1))
        ww = 2.d0*pi * freq(ifreq)
        DSN22 = -ww**2.d0 * MSN22 + dcmplx(0.d0,ww) * DPN22 + STN22 + 1.d0/dcmplx(0.d0,ww) * RRN22
        !do itmp=1,nshot
        P2(:,:,1) = dcmplx(Pext(1:nbc(2),:),0d0)
        !enddo
        !iparm(3) = 8	! number of processors
        iparm(8) = 0
        !iparm(59) = 2
        allocate( Dis2(nbc(2),nshot))
        if ((iSample == 0) .and. (ifreq == 1)) then
            call pardiso(pt, 1, 1, 13, 13, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 1, P2(:,1:nshot,1), Dis2, ierr)
        else
            call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 0, P2(:,1:nshot,1), Dis2, ierr)
        endif

        do ishot=1,nshot
            do itmp = 1 , nobs
                dobs(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp)=Dis2(id(onode(itmp),2),ishot)
            enddo
        enddo
        deallocate(P2)
        deallocate( Dis2)

    enddo


    call FiniteElement_PlaneStrain(CS_equi,nnode, nelmt, nmatrl, CSest(:,1), rhos,nuy,beta, MSN22, DPN22,  STN22, RRN22,x, z, enode, matrl, etype, epara, nbc, id, ia, ja)


    do ifreq = 1 , nfreq
        allocate(P2(nbc(2),nshot,1))
        ww = 2.d0*pi * freq(ifreq)
        DSN22 = -ww**2.d0 * MSN22 + dcmplx(0.d0,ww) * DPN22 + STN22 + 1.d0/dcmplx(0.d0,ww) * RRN22
        !do itmp=1,nshot
        P2(:,:,1) = dcmplx(Pext(1:nbc(2),:),0d0)
        !enddo
        !iparm(3) = 8	! number of processors
        iparm(8) = 0
        !iparm(59) = 2
        allocate( Dis2(nbc(2),nshot))
        if ((iSample == 0) .and. (ifreq == 1)) then
            call pardiso(pt, 1, 1, 13, 13, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 1, P2(:,1:nshot,1), Dis2, ierr)
        else
            call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 0, P2(:,1:nshot,1), Dis2, ierr)
        endif

        do ishot=1,nshot
            do itmp = 1 , nobs
                dprd(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp)=Dis2(id(onode(itmp),2),ishot)
            enddo
        enddo
        deallocate(P2)
        !allocate(P2(nbc(2),nshot,1))
        !P2 = dcmplx(0d0,0d0)
        !dDis2dm = dcmplx(0d0,0d0)
        v_P2=dcmplx(0d0,0d0)
        !$omp parallel do
        do imatrl = 1 , nmatrl
            call dSdmU_elastic(nshot,onode,nobs,nnode, nelmt, nmatrl, CSest(:,1), rhos,nuy,beta, x, z, enode, matrl, etype, epara, nbc, id, ia, ja, imatrl,ww, Dis2, v_P2,nnz_P2,pB_P2,pE_P2,ja_P2)
            !call dSdmU_elastic_assemble(nnode, nelmt, nmatrl, CSest(:,2), enode, matrl, nbc, id, ia, ja, imatrl, ww, Dis2, P2(:,imatrl),dStiffness,his_n,his_jtmp,his_yy)
        enddo
        !$omp end parallel do
        deallocate( Dis2)


        !deallocate (P2)

        !call sparse_get_valueTZ(nshot,nmatrl,nbc(2),ja_P2,pB_P2,pE_P2,v_P2,nnz_p2,P2)
        !deallocate(P2)
        allocate (invDSN22T(nbc(2),nobs))
        call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nobs, iparm, 0, eK, invDSN22T, ierr)

        allocate(dDis2dm(nobs,nmatrl),dDis2dmT(nmatrl,nobs))
        do ishot=1,nshot
            call As_m_Bd(v_P2(:,ishot),ja_P2,pB_P2,pE_P2,nnz_P2,nmatrl,nbc(2),nobs,invDSN22T,dDis2dmT)
            dDis2dm=transpose(dDis2dmT)
            !call zgemm('T', 'N', nobs, nmatrl, nbc(2), dcmplx(1d0,0d0), invDSN22T, nbc(2), P2,nbc(2), dcmplx(0d0,0d0), dDis2dm, nobs)

            do itmp = 1 , nobs
                !Hmat(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp,:) = dDis2dm(itmp,:)
                dUz_estimated(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp,:,1)= dDis2dm(itmp,:)
            enddo
        enddo
        deallocate(dDis2dm,dDis2dmT)
        deallocate (invDSN22T)

    enddo

    Uz_estimated(:,1)=dprd
    GA(1:nfreq*nobs*nshot,:)=dUz_estimated(:,:,1)
    GA(1+nfreq*nobs*nshot:2*nfreq*nobs*nshot,:)=dconjg(dUz_estimated(:,:,1))
    CA_Uz_cal(1:nfreq*nobs*nshot)=Uz_estimated(:,1)
    CA_Uz_cal(nfreq*nobs*nshot+1:nfreq*nobs*nshot*2)=dconjg(Uz_estimated(:,1))
    CA_Uz_obs(1:nfreq*nobs*nshot)=dobs
    CA_Uz_obs(nfreq*nobs*nshot+1:nfreq*nobs*nshot*2)=dconjg(dobs)





    do while (iflag)
        iflag_group=1
        call system_clock(tclock1)
        write (*,*) '- Group :',igroup,' / ',ngroup
        iSample = iSample + 1
        write (*,*) '       + Iteration :',iSample,' / ',nSample
        call hessgrad(GA,CA_Uz_cal,CA_Uz_obs,inv_Ca_D,inv_sigma_M2,CSest(:,1),CS0,2*nfreq*nobs*nshot,nmatrl,hess,grad)

        !call Cmd_Cdd_ensemble(GA,CS0,nmatrl,2*nfreq*nobs*nshot,C_md,Ldd,invcovUz,sigma_M2,CA_Uz_cal,CA_Uz_obs,Csest(:,1))
        !call solution_ensemble(GA,C_md,Ldd,CA_Uz_cal,CA_Uz_obs,CSest(:,1),CS0,nmatrl,2*nfreq*nobs*nshot,dCS)


        !hess=hess1
        !grad=grad1
        !call random_stdnormal(randCs,nmatrl)

        hess_real=dreal(hess)
        grad_real=dreal(grad)
        call system_clock(tclock2)
        call solveAxB(nmatrl,1, hess_real,  grad_real,dCs)

        call system_clock(tclock22, clock_rate2)
        elapsed_time2 = float(tclock22 - tclock2) / float(clock_rate2)
        write(*,*) "dCs", elapsed_time2


        !call cholesky_lower(nmatrl, hess,cholhess)
        CStemp =randCs
        randCs=0d0
        !call system_clock(tclock2)
        !call DTRTRS('L','T','N',nmatrl,1,cholhess,nmatrl,CStemp,nmatrl,INFO)
        CSest(:,2) =CSest(:,1) -alpha_coef*dCs+beta_coef*CStemp
        !CSest(:,2) = CS0-alpha_coef*dCs+beta_coef*CStemp
        !call system_clock(tclock22, clock_rate2)
        !elapsed_time2 = float(tclock22 - tclock2) / float(clock_rate2)
        !write(*,*) "DTRTRS", elapsed_time2


        !CStemp=CSest(:,1) - (CSest(:,1)+dCs)
        !call system_clock(tclock2)
        !call postprop(hess,CStemp,1d0,nmatrl,post0)
        !call system_clock(tclock22, clock_rate2)
        !elapsed_time2 = float(tclock22 - tclock2) / float(clock_rate2)
        !write(*,*) "postprop", elapsed_time2*4
        !
        !CStemp = CSest(:,2) - (CSest(:,1)+dCs)
        !call postprop(hess,CStemp,1d0,nmatrl,post1)
        !CStemp = CSest(:,2) - (CSest(:,1)+alpha_coef*dCs)
        !call postprop(hess,CStemp,1d0/beta_coef**2d0,nmatrl,prop01)

        call FiniteElement_PlaneStrain(CS_equi,nnode, nelmt, nmatrl, CSest(:,2), rhos,nuy,beta, MSN22, DPN22,  STN22, RRN22,x, z, enode, matrl, etype, epara, nbc, id, ia, ja)

        do ifreq = 1 , nfreq
            allocate(P2(nbc(2),nshot,1))
            ww = 2.d0*pi * freq(ifreq)
            DSN22 = -ww**2.d0 * MSN22 + dcmplx(0.d0,ww) * DPN22 + STN22 + 1.d0/dcmplx(0.d0,ww) * RRN22
            P2(:,:,1) = dcmplx(Pext(1:nbc(2),:),0d0)
            iparm(8) = 0
            allocate( Dis2(nbc(2),nshot))
            iparm(8) = 0
            call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 0, P2(:,1:nshot,1), Dis2, ierr)


            do ishot=1,nshot
                do itmp = 1 , nobs
                    dprd(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp)=Dis2(id(onode(itmp),2),ishot)
                enddo
            enddo
            deallocate(P2)
            v_P2=dcmplx(0d0,0d0)
            !$omp parallel do
            do imatrl = 1 , nmatrl
                call dSdmU_elastic(nshot,onode,nobs,nnode, nelmt, nmatrl, CSest(:,2), rhos,nuy,beta, x, z, enode, matrl, etype, epara, nbc, id, ia, ja, imatrl,ww, Dis2, v_P2,nnz_P2,pB_P2,pE_P2,ja_P2)
            enddo
            !$omp end parallel do
            deallocate( Dis2)
            allocate (invDSN22T(nbc(2),nobs))
            call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nobs, iparm, 0, eK, invDSN22T, ierr)

            allocate(dDis2dm(nobs,nmatrl),dDis2dmT(nmatrl,nobs))
            do ishot=1,nshot
                call As_m_Bd(v_P2(:,ishot),ja_P2,pB_P2,pE_P2,nnz_P2,nmatrl,nbc(2),nobs,invDSN22T,dDis2dmT)
                dDis2dm=transpose(dDis2dmT)
                do itmp = 1 , nobs
                    dUz_estimated(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp,:,1)= dDis2dm(itmp,:)
                enddo
            enddo
            deallocate(dDis2dm,dDis2dmT)
            deallocate (invDSN22T)

        enddo


        Uz_estimated(:,2)=dprd

        GA(1:nfreq*nobs*nshot,:)=dUz_estimated(:,:,1)
        GA(1+nfreq*nobs*nshot:2*nfreq*nobs*nshot,:)=dconjg(dUz_estimated(:,:,1))
        CA_Uz_cal(1:nfreq*nobs*nshot)=Uz_estimated(:,2)
        CA_Uz_cal(nfreq*nobs*nshot+1:nfreq*nobs*nshot*2)=dconjg(Uz_estimated(:,2))
        !call hessgrad(GA,CA_Uz_cal,CA_Uz_obs,inv_Ca_D,inv_sigma_M2,CSest(:,2),CS0,2*nfreq*nobs*nshot,nmatrl,hess,grad)
        !call Cmd_Cdd_ensemble(GA,CS0,nmatrl,2*nfreq*nobs*nshot,C_md,Ldd,invcovUz,sigma_M2,CA_Uz_cal,CA_Uz_obs,Csest(:,2))
        !call solution_ensemble(GA,C_md,Ldd,CA_Uz_cal,CA_Uz_obs,CSest(:,2),CS0,nmatrl,2*nfreq*nobs*nshot,dCS)

        derr_CS = dabs( CSest(:,2)- CSest(:,1))
        Uz_estimated(:,1) = Uz_estimated(:,2)
        CSest(:,1) = CSest(:,2)


        write (outest,'(i10,i10,<nmatrl>f15.5)') igroup,iSample, (CSest(1:nmatrl,1))
        !write (outgama,'(i10,<nmatrl>f15.5)') iSample, gama



        dtmp2 = 0.d0
        dtmp1 = 0.d0
        do imatrl=1,nmatrl

            dtmp2 = dtmp2 + derr_CS(imatrl)**2d0
            dtmp1 = dtmp1 + CSest(imatrl,1)**2d0
        enddo

        !dtmp1 = dsqrt(dtmp1)
        dtmp2 = dsqrt(dtmp2/dtmp1)
        write(*,*) '       + rate of change : ', dtmp2*100d0
        if ((dtmp2 < 0.001d0) .and. (iSample>5) ) iflag_group = 0
        if (iSample == nSample) iflag_group = 0
        call system_clock(tclock2, clock_rate)
        elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
        write(*,*) "       + time for 1 iterations : ", elapsed_time

        if (iflag_group == 0) then
            igroup=igroup+1
            if (igroup <= ngroup) then
                iSample=0
                freq=group_frequency(:,igroup)
                !call update_inv_sigma_M2(GA,CA_Uz_cal,CA_Uz_obs,inv_Ca_D, inv_sigma_M2,CSest(:,1),2*nfreq*nobs*nshot,nmatrl)

                call FiniteElement_PlaneStrain(CS_equi,nnode, nelmt, nmatrl, CS, rhos,nuy,beta, MSN22, DPN22,  STN22, RRN22,x, z, enode, matrl, etype, epara, nbc, id, ia, ja)

                do ifreq = 1 , nfreq
                    allocate(P2(nbc(2),nshot,1))
                    ww = 2.d0*pi * freq(ifreq)
                    DSN22 = -ww**2.d0 * MSN22 + dcmplx(0.d0,ww) * DPN22 + STN22 + 1.d0/dcmplx(0.d0,ww) * RRN22
                    !do itmp=1,nshot
                    P2(:,:,1) = dcmplx(Pext(1:nbc(2),:),0d0)
                    !enddo
                    !iparm(3) = 8	! number of processors
                    iparm(8) = 0
                    !iparm(59) = 2
                    allocate( Dis2(nbc(2),nshot))

                    call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 0, P2(:,1:nshot,1), Dis2, ierr)
                    do ishot=1,nshot
                        do itmp = 1 , nobs
                            dobs(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp)=Dis2(id(onode(itmp),2),ishot)
                        enddo
                    enddo
                    deallocate(P2)
                    deallocate( Dis2)

                enddo
                CA_Uz_obs(1:nfreq*nobs*nshot)=dobs
                CA_Uz_obs(nfreq*nobs*nshot+1:nfreq*nobs*nshot*2)=dconjg(dobs)
                call FiniteElement_PlaneStrain(CS_equi,nnode, nelmt, nmatrl, CSest(:,2), rhos,nuy,beta, MSN22, DPN22,  STN22, RRN22,x, z, enode, matrl, etype, epara, nbc, id, ia, ja)

                do ifreq = 1 , nfreq
                    allocate(P2(nbc(2),nshot,1))
                    ww = 2.d0*pi * freq(ifreq)
                    DSN22 = -ww**2.d0 * MSN22 + dcmplx(0.d0,ww) * DPN22 + STN22 + 1.d0/dcmplx(0.d0,ww) * RRN22
                    P2(:,:,1) = dcmplx(Pext(1:nbc(2),:),0d0)
                    iparm(8) = 0
                    allocate( Dis2(nbc(2),nshot))
                    iparm(8) = 0
                    call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 0, P2(:,1:nshot,1), Dis2, ierr)


                    do ishot=1,nshot
                        do itmp = 1 , nobs
                            dprd(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp)=Dis2(id(onode(itmp),2),ishot)
                        enddo
                    enddo
                    deallocate(P2)
                    v_P2=dcmplx(0d0,0d0)
                    !$omp parallel do
                    do imatrl = 1 , nmatrl
                        call dSdmU_elastic(nshot,onode,nobs,nnode, nelmt, nmatrl, CSest(:,2), rhos,nuy,beta, x, z, enode, matrl, etype, epara, nbc, id, ia, ja, imatrl,ww, Dis2, v_P2,nnz_P2,pB_P2,pE_P2,ja_P2)
                    enddo
                    !$omp end parallel do
                    deallocate( Dis2)
                    allocate (invDSN22T(nbc(2),nobs))
                    call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nobs, iparm, 0, eK, invDSN22T, ierr)

                    allocate(dDis2dm(nobs,nmatrl),dDis2dmT(nmatrl,nobs))
                    do ishot=1,nshot
                        call As_m_Bd(v_P2(:,ishot),ja_P2,pB_P2,pE_P2,nnz_P2,nmatrl,nbc(2),nobs,invDSN22T,dDis2dmT)
                        dDis2dm=transpose(dDis2dmT)
                        do itmp = 1 , nobs
                            dUz_estimated(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp,:,1)= dDis2dm(itmp,:)
                        enddo
                    enddo
                    deallocate(dDis2dm,dDis2dmT)
                    deallocate (invDSN22T)

                enddo


                Uz_estimated(:,2)=dprd

                GA(1:nfreq*nobs*nshot,:)=dUz_estimated(:,:,1)
                GA(1+nfreq*nobs*nshot:2*nfreq*nobs*nshot,:)=dconjg(dUz_estimated(:,:,1))
                CA_Uz_cal(1:nfreq*nobs*nshot)=Uz_estimated(:,2)
                CA_Uz_cal(nfreq*nobs*nshot+1:nfreq*nobs*nshot*2)=dconjg(Uz_estimated(:,2))
                call hessgrad(GA,CA_Uz_cal,CA_Uz_obs,inv_Ca_D,inv_sigma_M2,CSest(:,2),CS0,2*nfreq*nobs*nshot,nmatrl,hess,grad)
                !call Cmd_Cdd_ensemble(GA,CS0,nmatrl,2*nfreq*nobs*nshot,C_md,Ldd,invcovUz,sigma_M2,CA_Uz_cal,CA_Uz_obs,Csest(:,2))
                !call solution_ensemble(GA,C_md,Ldd,CA_Uz_cal,CA_Uz_obs,CSest(:,2),CS0,nmatrl,2*nfreq*nobs*nshot,dCS)



                CS0=CSest(1:nmatrl,2)
            endif
        endif

        if (iflag_group == 0) then
            derr_CS = dabs( (CSest(:,2))- (CS))

            dtmp2 = 0.d0
            dtmp1 = 0.d0
            do imatrl=1,nmatrl

                dtmp2 = dtmp2 + derr_CS(imatrl)**2d0
                dtmp1 = dtmp1 + CS(imatrl)**2d0
            enddo

            !dtmp1 = dsqrt(dtmp1)
            dtmp2 = dsqrt(dtmp2/dtmp1)
            write(*,*) '       + accuracy : ', dtmp2*100d0,' %'
        endif

        if (igroup== ngroup+1) iflag = 0
    enddo

    deallocate (GA,CA_Uz_cal,CA_Uz_obs,inv_Ca_D,Uz_estimated,dUz_estimated,dprd,dobs )
    allocate(Uz_estimated(nfreq_total*nobs*nshot,2),dUz_estimated (nfreq_total*nobs*nshot,nmatrl,1))
    allocate(dprd(nobs*nfreq_total*nshot),dobs(nobs*nfreq_total*nshot))

    allocate (GA(2*nfreq_total*nobs*nshot,nmatrl),CA_Uz_cal(2*nfreq_total*nobs*nshot),CA_Uz_obs(2*nfreq_total*nobs*nshot),inv_Ca_D(nfreq_total*nobs*nshot*2))
    inv_Ca_D=1d0/covUz(1)

    call FiniteElement_PlaneStrain(CS_equi,nnode, nelmt, nmatrl, CSest(1:nmatrl,2), rhos,nuy,beta, MSN22, DPN22,  STN22, RRN22,x, z, enode, matrl, etype, epara, nbc, id, ia, ja)

    do ifreq = 1 , nfreq_total
        allocate(P2(nbc(2),nshot,1))
        ww = 2.d0*pi * freq_total(ifreq)
        DSN22 = -ww**2.d0 * MSN22 + dcmplx(0.d0,ww) * DPN22 + STN22 + 1.d0/dcmplx(0.d0,ww) * RRN22
        P2(:,:,1) = dcmplx(Pext(1:nbc(2),:),0d0)
        iparm(8) = 0
        allocate( Dis2(nbc(2),nshot))
        iparm(8) = 0
        call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 0, P2(:,1:nshot,1), Dis2, ierr)

        deallocate(P2)
        v_P2=dcmplx(0d0,0d0)
        !$omp parallel do
        do imatrl = 1 , nmatrl
            call dSdmU_elastic(nshot,onode,nobs,nnode, nelmt, nmatrl, CSest(:,1), rhos,nuy,beta, x, z, enode, matrl, etype, epara, nbc, id, ia, ja, imatrl,ww, Dis2, v_P2,nnz_P2,pB_P2,pE_P2,ja_P2)
        enddo
        !$omp end parallel do
        deallocate( Dis2)
        allocate (invDSN22T(nbc(2),nobs))
        call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nobs, iparm, 0, eK, invDSN22T, ierr)

        allocate(dDis2dm(nobs,nmatrl),dDis2dmT(nmatrl,nobs))
        do ishot=1,nshot
            call As_m_Bd(v_P2(:,ishot),ja_P2,pB_P2,pE_P2,nnz_P2,nmatrl,nbc(2),nobs,invDSN22T,dDis2dmT)
            dDis2dm=transpose(dDis2dmT)
            do itmp = 1 , nobs
                dUz_estimated(nobs*nshot*(ifreq-1)+(ishot-1)*nobs+itmp,:,1)= dDis2dm(itmp,:)
            enddo
        enddo
        deallocate(dDis2dm,dDis2dmT)
        deallocate (invDSN22T)

    enddo

    GA(1:nfreq_total*nobs*nshot,:)=dUz_estimated(:,:,1)
    deallocate(dUz_estimated)
    do imatrl=1,nmatrl
        GA(1+nfreq_total*nobs*nshot:2*nfreq_total*nobs*nshot,imatrl)=dconjg(GA(1:nfreq_total*nobs*nshot,imatrl))
    enddo

    !call update_inv_sigma_M2(GA,CA_Uz_cal,CA_Uz_obs,inv_Ca_D, inv_sigma_M2,CSest,2*nfreq*nobs*nshot,nmatrl)

    !call update_inv_sigma_M2(GA,inv_Ca_D, inv_sigma_M2,2*nfreq_total*nobs*nshot,nmatrl)
    call update_inv_sigma_M2_rev(GA,inv_Ca_D, sigma_M2,2*nfreq_total*nobs*nshot,nmatrl,inv_sigma_M2)
    
    !write (outcovM,'(<nmatrl>f15.5)')  1d0/dsqrt(inv_sigma_M2(1:nmatrl))
    !call update_inv_sigma_M2(GA,CA_Uz_cal,CA_Uz_obs,inv_Ca_D, inv_sigma_M2,CSest(:,1),2*nfreq*nobs*nshot,nmatrl)

    do imatrl =1,nmatrl
        write (outcovM,'(<nmatrl>f15.5)')  inv_sigma_M2(imatrl,1:nmatrl)
    enddo
    deallocate (ia, ja)

    end subroutine MCMCcore

