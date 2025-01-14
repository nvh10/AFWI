    subroutine MCMCReference(CS_equi_true,nshot,nx, nz, npml, dx, dz, nnode, nelmt, nmatrl, CS, rhos,nuy,beta, x, z, enode, matrl, etype, epara, nbc, bc, id, &
        nfreq_total, freq, Pext, inpar, outobs, nout, onode,dobs, iflag, Pfac)
    implicit none
      integer nfreq_total
          integer nbc(3), bc(nnode), id(nnode,2), ndof, ndof2,nshot
    integer nx, nz, npml, nnode, nelmt, nmatrl, nout, onode(nout),iflag
    double precision rhos,nuy,beta, CS(nmatrl),  dobs(2*nout*nfreq_total),Pfac
double complex CS_equi_true(nelmt)
    ! node data
    double precision x(nnode), z(nnode)

    ! element data
    integer enode(nelmt,0:9), matrl(nelmt), etype(nelmt)
    double complex epara(nelmt,2)

    ! boundary condition data


    double precision Pext(nnode*2,nshot)
    !double precision Pext(nnode) !anisotropic
    ! frequency data
  
    double precision dfreq, freq(nfreq_total)

    integer outsgb, outrhb, outobs, inexa, inpar, inrfa, outerr

    double precision pi, ww
   

    double complex,allocatable :: MSN22(:), DPN22(:), STN22(:), RRN22(:),DSN22(:)
    double complex  Dis2(nbc(2),nshot), P2(nbc(2),nshot)

    integer nn, n1, n2, nfix, ne, nm, ifreq, itmp, jtmp, inode, jnode, idof, jdof, ielmt, iCapp, iregion, itime


    integer,allocatable :: ia(:), ja(:), idcheck(:)
    integer equation, node, position, iparm(64), idum, ierr, iia
    integer*8 pt(64)

    double precision delta, alph1, a0, a1, a2, a3, a4, a5, a6, a7

    integer ltmp, nx1, jtime, otime(10)
    double precision dx, dz, x1, z1, dtmp0, dtmp1, err1, err2, outdata(5,nnode*2)

    pi = 4.d0 * datan(1.d0)

    allocate (ia(nbc(2)+1))

    do itmp = 1 , nbc(2)+1

        read (inpar,*) ia(itmp)

    enddo

    read (inpar,*)
    allocate (ja(ia(nbc(2)+1)-1))
    allocate ( MSN22(ia(nbc(2)+1)-1), DPN22(ia(nbc(2)+1)-1), STN22(ia(nbc(2)+1)-1), RRN22(ia(nbc(2)+1)-1), DSN22(ia(nbc(2)+1)-1))

    do itmp = 1 , ia(nbc(2)+1)-1

        read (inpar,*) ja(itmp)

    enddo

    MSN22 = 0.d0
    DPN22 = 0.d0
    STN22 = 0.d0
    RRN22 = 0.d0


    !call FiniteElement(nnode, nelmt, nmatrl, CS, rhos, MSN11, MSN12, MSN22, DPN11, DPN12, DPN22, STN11, STN12, STN22, RRN11, RRN12, RRN22, &
    !    x, z, enode, matrl, etype, epara, nbc, id, ia, ja)
    !call FEM_assemble(nnode, nelmt, nmatrl,CS, rhos,nuy,beta,x, z, enode, matrl, etype, epara, nbc, id,ia,ja, elmtM_all, elmtC_all, elmtK_all, elmtR_all, elmtM_all_n, elmtC_all_n, elmtK_all_n, elmtR_all_n,MSN22,DPN22,STN22,RRN22)
    call FiniteElement_PlaneStrain(CS_equi_true,nnode, nelmt, nmatrl, CS, rhos,nuy,beta, MSN22, DPN22,  STN22, RRN22,x, z, enode, matrl, etype, epara, nbc, id, ia, ja)





    call pardisoinit(pt, 13, iparm)
    if (iflag == 1) write (outobs,'(10x,<2*nout>i16)') (onode(inode), onode(inode), inode = 1 , nout)

    Pfac = 0.d0

    do ifreq = 1 , nfreq_total

        write (*,'(1h+,a14,i10,a3,i10)') '   Freq step :',ifreq,' / ',nfreq_total

        ww = 2.d0*pi * freq(ifreq)

        DSN22 = -ww**2.d0 * MSN22 + dcmplx(0.d0,ww) * DPN22 + STN22 + 1.d0/dcmplx(0.d0,ww) * RRN22
       do itmp=1,nshot
        P2(:,itmp) = dcmplx(Pext(1:nbc(2),itmp),0d0)
        enddo


        !iparm(3) = 8	! number of processors
        iparm(8) = 0

        Dis2 = dcmplx(0d0,0d0)

        if (ifreq == 1) then

            call pardiso(pt, 1, 1, 13, 13, nbc(2), DSN22, ia, ja, idum,nshot, iparm, 1, P2, Dis2, ierr)

        else

            call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 0, P2, Dis2, ierr)

        endif

        do itmp=1,nshot
        if (iflag == 1) write (outobs,'(f10.5,<2*nout>e16.5E3)') freq(ifreq), (dreal(Dis2(id(onode(inode),2),itmp)), dimag(Dis2(id(onode(inode),2),itmp)), inode = 1 , nout)

        if (maxval(cdabs(Dis2(id(onode(1:nout),2),itmp))) > Pfac) Pfac = maxval(cdabs(Dis2(id(onode(1:nout),2),itmp)))
         enddo
        
        !do itmp = 1 , nout
        !
        !    dobs(2*nout*(ifreq-1)+2*itmp-1) = dreal(Dis2(id(onode(itmp),2)))
        !    dobs(2*nout*(ifreq-1)+2*itmp  ) = dimag(Dis2(id(onode(itmp),2)))
        !
        !enddo
        !
    enddo

    Pfac = 1.d0 / Pfac


    call pardiso(pt, 1, 1, 13, -1, nbc(2), DSN22, ia, ja, idum, nshot, iparm, 0, P2, Dis2, ierr)

    deallocate (ia, ja, MSN22, DPN22, STN22, RRN22, DSN22)

    end subroutine MCMCReference