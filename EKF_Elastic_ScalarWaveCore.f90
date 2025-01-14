    subroutine EKF_Elastic_ScalarWaveCore(nx, nz, npml, dx, dz, nnode, nelmt, nmatrl, CS, rhos,nuy,beta, x, z, enode, matrl, etype, epara, nbc, bc, id, &
        nfreq, freq, Pext, inpar, inobs, outest, outrsp, outrsl,outPst, nout, onode, CSest, Pst, iQmat, Robs, niter, elmtM_all, elmtC_all, elmtK_all, elmtR_all, elmtM_all_n, elmtC_all_n, elmtK_all_n, elmtR_all_n,dobs)
    !use mkl_cluster_sparse_solver
    implicit none
    !include 'mpif.h'
    integer nnode, nelmt, nmatrl, niter, norder, nout, onode(nout), npml
    double precision rhos(nmatrl), CS(nmatrl), CSest(nmatrl), CS0(nmatrl),nuy(nmatrl),beta(nmatrl)
    !TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), allocatable  :: pt(:)
    integer*8 pt(64)

    integer INFO
    double precision Pst(nmatrl,nmatrl), Robs(2*nout*nfreq,2*nout*nfreq)
    double precision dobs(2*nout*nfreq), dprd(2*nout*nfreq), &
        derr(2*nout*nfreq)
    double precision,allocatable ::Hmat(:,:), Dmat(:,:), Smat(:,:),        Kgain(:,:)
    double precision Cmat(nmatrl,nmatrl), &
        Emat(nmatrl,nmatrl), dcon(nmatrl), &
        iQmat(nmatrl,nmatrl)
    double complex elmtM_all(18,18,nelmt), elmtC_all(18,18,nelmt), elmtK_all(18,18,nelmt), elmtR_all(18,18,nelmt)
    double complex elmtM_all_n(18,18,nelmt), elmtC_all_n(18,18,nelmt), elmtK_all_n(18,18,nelmt), elmtR_all_n(18,18,nelmt)

    ! node data
    double precision x(nnode), z(nnode)

    ! element data
    integer enode(nelmt,0:9), matrl(nelmt), etype(nelmt)
    double precision epara(nelmt,2)

    ! boundary condition data
    integer nbc(3), bc(nnode), id(nnode,2), ndof, ndof2

    double precision Pext(nnode*2) !, Phis(0:niter)

    ! frequency data
    integer nfreq
    double precision freq(nfreq)

    integer outsgb, outrhb, outest, inexa, inpar, inrfa, outerr, inobs, outrsp, outrsl,outPst

    double precision pi, ww


    double complex,allocatable :: MSN22(:), DPN22(:), STN22(:), RRN22(:),DSN22(:)
    !double complex,allocatable :: MSN22(:,:), DPN22(:,:), STN22(:,:), RRN22(:,:),DSN22(:,:)
    double complex,allocatable :: invDSN22(:,:),eK(:,:),invDSN22T(:,:)
    double complex,allocatable ::  dDis2dm(:,:), P2(:,:)
    !double complex,allocatable ::  dDis2dm(:,:), P2(:,:)
    double complex,allocatable :: Dis2(:)

    integer nn, n1, n2, nfix, ne, nm, ifreq, itmp, jtmp, inode, jnode, idof, jdof, ielmt, iCapp, iregion, itime, imatrl, ktmp


    integer,allocatable :: ia(:), ja(:), idcheck(:)
    integer equation, node, position, iparm(64), idum, ierr, iia,i


    integer nx, nz, ltmp, nx1, jtime, otime(10), iflag, ndiv,idiv,ninterval
    double precision dx, dz, x1, z1, dtmp0, dtmp1, dtmp2, err1, err2, outdata(5,nnode)
    integer(kind=8) :: tclock1, tclock2, clock_rate
    real(kind=8) :: elapsed_time
    integer(kind=8) :: tclock11, tclock22, clock_rate2
    real(kind=8) :: elapsed_time2

    ndiv=1
    pi = 4.d0 * datan(1.d0)
    ninterval=nmatrl/ndiv
    allocate(Dis2(nbc(2)))
    !allocate(dDis2dm(nbc(2),ninterval), P2(nbc(2),ninterval))
    allocate (ia(nbc(2)+1))

    do itmp = 1 , nbc(2)+1

        read (inpar,*) ia(itmp)

    enddo

    read (inpar,*)

    allocate (ja(ia(nbc(2)+1)-1))
    !allocate ( MSN22(ia(nbc(2)+1)-1), DPN22(ia(nbc(2)+1)-1), STN22(ia(nbc(2)+1)-1), RRN22(ia(nbc(2)+1)-1), DSN22(ia(nbc(2)+1)-1))
    do itmp = 1 , ia(nbc(2)+1)-1

        read (inpar,*) ja(itmp)

    enddo

    call pardisoinit(pt, 13, iparm)

    write (outest,'(10x,<nmatrl>i15)') (itmp, itmp = 1 , nmatrl)
    write (outest,'(i10,<nmatrl>f15.5)') 0, CSest(1:nmatrl)

    write (outrsp,'(10x,<2*nout>i16)') (onode(inode), onode(inode), inode = 1 , nout)

    dobs = 0.d0

    read (inobs,*)

    do ifreq = 1 , nfreq

        read (inobs,*) dtmp1, (dobs(2*nout*(ifreq-1)+2*itmp-1), dobs(2*nout*(ifreq-1)+2*itmp), itmp = 1 , nout)

    enddo

    itime = 0
    CS0 = CSest



    iflag = 1

    allocate(eK(nbc(2),nout))
    eK=0d0

    do itmp=1,nout
        eK(id(onode(itmp),2),itmp)=1d0
    enddo

    do while (iflag)
        call system_clock(tclock1)
        allocate( P2(nbc(2),nmatrl))
        allocate(dDis2dm(nout,nmatrl))
        !allocate(dDis2dm(nbc(2),nmatrl))
        allocate ( MSN22(ia(nbc(2)+1)-1), DPN22(ia(nbc(2)+1)-1), STN22(ia(nbc(2)+1)-1), RRN22(ia(nbc(2)+1)-1), DSN22(ia(nbc(2)+1)-1))
        !allocate ( MSN22(nbc(2),nbc(2)), DPN22(nbc(2),nbc(2)), STN22(nbc(2),nbc(2)), RRN22(nbc(2),nbc(2)), DSN22(nbc(2),nbc(2)))

        allocate (invDSN22T(nbc(2),nout))

        allocate(Hmat(2*nout*nfreq,nmatrl))
        itime = itime + 1

        write (*,'(1h+,a14,i10,a3,i10)') '   Iteration :',itime,' / ',niter

        MSN22 = 0.d0
        DPN22 = 0.d0
        STN22 = 0.d0
        RRN22 = 0.d0
        call FiniteElement_PlaneStrain(nnode, nelmt, nmatrl, CSest, rhos,nuy,beta, MSN22, DPN22,  STN22, RRN22,x, z, enode, matrl, etype, epara, nbc, id, ia, ja)
        !call FEM_assemble(nnode, nelmt, nmatrl,CSest, rhos,nuy,beta,x, z, enode, matrl, etype, epara, nbc, id,ia,ja, elmtM_all, elmtC_all, elmtK_all, elmtR_all, elmtM_all_n, elmtC_all_n, elmtK_all_n, elmtR_all_n,MSN22,DPN22,STN22,RRN22)

        call system_clock(tclock11)
        do ifreq = 1 , nfreq
            !write (*,'(1h+,a14,i10,a3,i10)') '   frequency :',ifreq,' / ',nfreq

            ww = 2.d0*pi * freq(ifreq)

            DSN22 = -ww**2.d0 * MSN22 + dcmplx(0.d0,ww) * DPN22 + STN22 + 1.d0/dcmplx(0.d0,ww) * RRN22

            P2(:,1) = Pext(1:nbc(2))

            iparm(3) = 8	! number of processors
            iparm(8) = 0

            Dis2 = 0.d0


            if ((itime == 1) .and. (ifreq == 1)) then

                call pardiso(pt, 1, 1, 13, 13, nbc(2), DSN22, ia, ja, idum, 1, iparm, 1, P2(:,1), Dis2, ierr)


            else

                call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, 1, iparm, 0, P2(:,1), Dis2, ierr)


            endif
            invDSN22T=0d0
            call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nout, iparm, 0, eK, invDSN22T, ierr)
            !invDSN22=transpose(invDSN22T)

            do itmp = 1 , nout

                dprd(2*nout*(ifreq-1)+2*itmp-1) = dreal(Dis2(id(onode(itmp),2)))
                dprd(2*nout*(ifreq-1)+2*itmp  ) = dimag(Dis2(id(onode(itmp),2)))

            enddo

            !
            P2 = 0.d0

            dDis2dm = 0.d0
            do imatrl = 1 , nmatrl

                call dSdmU_elastic(nnode, nelmt, nmatrl, CSest, rhos,nuy,beta, x, z, enode, matrl, etype, epara, nbc, id, ia, ja, imatrl,ww, Dis2, P2(:,imatrl))

            enddo
            call zgemm('T', 'N', nout, nmatrl, nbc(2), dcmplx(1d0,0d0), invDSN22T, nbc(2), P2,nbc(2), dcmplx(0d0,0d0), dDis2dm, nout)

            !call zgemm('N', 'N', nout, nmatrl, nbc(2), dcmplx(1d0,0d0), invDSN22, nout, P2,nbc(2), dcmplx(0d0,0d0), dDis2dm, nout)

            do itmp = 1 , nout
                Hmat(2*nout*(ifreq-1)+2*itmp-1,:) = dreal(dDis2dm(itmp,:))
                Hmat(2*nout*(ifreq-1)+2*itmp ,:) = dimag(dDis2dm(itmp,:))
            enddo

            !call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nmatrl, iparm, 0, P2, dDis2dm, ierr)
            !!
            !do itmp = 1 , nout
            !    Hmat(2*nout*(ifreq-1)+2*itmp-1,:) = dreal(dDis2dm(id(onode(itmp),2),:))
            !    Hmat(2*nout*(ifreq-1)+2*itmp ,:) = dimag(dDis2dm(id(onode(itmp),2),:))
            !enddo


        enddo
        deallocate(dDis2dm, P2)
        deallocate ( MSN22, DPN22, STN22, RRN22, DSN22)
        deallocate (invDSN22T)
        allocate( Dmat(2*nout*nfreq,nmatrl), Smat(2*nout*nfreq,2*nout*nfreq),  Kgain(nmatrl,2*nout*nfreq))

        call system_clock(tclock22, clock_rate2)
        elapsed_time2 = float(tclock22 - tclock11) / float(clock_rate2)
        write(*,*) "Linearization", elapsed_time2

        Dmat = 0.d0
        call dgemm('N', 'N', 2*nout*nfreq, nmatrl, nmatrl, 1.d0, Hmat, 2*nout*nfreq, Pst, nmatrl, 0.d0, Dmat(1:2*nout*nfreq,:), 2*nout*nfreq)

        Smat = 0.d0
        call dgemm('N', 'T', 2*nout*nfreq, 2*nout*nfreq, nmatrl, 1.d0, Dmat(1:2*nout*nfreq,:), 2*nout*nfreq, Hmat, 2*nout*nfreq, 0.d0, Smat(1:2*nout*nfreq,1:2*nout*nfreq), 2*nout*nfreq)

        Smat = Smat + Robs


        call dpotrf('L', 2*nout*nfreq, Smat, 2*nout*nfreq, INFO)

        Kgain = transpose(Dmat)

        call dtrsm('R', 'L', 'T', 'N', nmatrl, 2*nout*nfreq, 1.d0, Smat, 2*nout*nfreq, Kgain, nmatrl)
        call dtrsm('R', 'L', 'N', 'N', nmatrl, 2*nout*nfreq, 1.d0, Smat, 2*nout*nfreq, Kgain, nmatrl)

        derr(1:2*nout*nfreq) = dobs(1:2*nout*nfreq) - dprd(1:2*nout*nfreq)

        CSest = CSest + matmul(Kgain, derr)

        call dgemm('N', 'N', nmatrl, nmatrl, 2*nout*nfreq, -1.d0, Kgain, nmatrl, Dmat, 2*nout*nfreq, 1.d0, Pst, nmatrl)

        write (outest,'(i10,<nmatrl>f15.5)') itime, CSest(1:nmatrl)

        dtmp1 = 0.d0
        dtmp2 = 0.d0

        do ifreq = 1 , nfreq
            do itmp = 1 , nout

                dtmp1 = dtmp1 + derr(2*nout*(ifreq-1)+2*itmp-1)**2 + derr(2*nout*(ifreq-1)+2*itmp)**2
                dtmp2 = dtmp2 + dobs(2*nout*(ifreq-1)+2*itmp-1)**2 + dobs(2*nout*(ifreq-1)+2*itmp)**2

            enddo
        enddo

        dtmp1 = dsqrt(dtmp1)
        dtmp2 = dsqrt(dtmp2)

        if (dtmp1/dtmp2 < 0.001d0) iflag = 0
        if (itime == niter) iflag = 0
        deallocate(Dmat, Smat, Kgain,Hmat)
        call system_clock(tclock2, clock_rate)
        elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
        write(*,*) "inversion time =", elapsed_time
        write(*,*) "CSest(1)=", CSest(1),"CSest(end)=", CSest((nx-2*npml)*(nz-npml))
        if (iflag==200) then
            do itmp = 1 , nz - npml

                write (outPst,'(i10,<nx-2*npml>e16.5E3)') itmp, (Pst((nx-2*npml)*(itmp-1)+jtmp,(nx-2*npml)*(itmp-1)+jtmp), jtmp = 1 , nx-2*npml)

            enddo

        endif
    enddo
    do itmp = 1 , nfreq

        write (outrsp,'(f10.5,<2*nout>e16.5E3)') freq(itmp), (dprd(2*nout*(itmp-1)+2*inode-1), dprd(2*nout*(itmp-1)+2*inode), inode = 1 , nout)

    enddo

    write (outrsl) CSest

    do itmp = 1 , nz - npml

        write (outPst,'(i10,<nx-2*npml>e16.5E3)') itmp, (Pst((nx-2*npml)*(itmp-1)+jtmp,(nx-2*npml)*(itmp-1)+jtmp), jtmp = 1 , nx-2*npml)

    enddo
    call pardiso(pt, 1, 1, 13, -1, nbc(2), DSN22, ia, ja, idum, 1, iparm, 0, P2, Dis2, ierr)

    deallocate (ia, ja)

    end subroutine EKF_Elastic_ScalarWaveCore