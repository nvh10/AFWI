    program MCMC
    !use mkl_cluster_sparse_solver

    implicit none
    !include 'mpif.h'
    ! input & output
    integer ininp, infar, outsgb, outrhb, outout,outcovM, inexa, inpar, outerr, outdat, outhis, nobs, inobs, outobs, outrsp, outest, outrsl, nshot, nfreq_total, ngroup
    integer outgama, outpost1, outpost0,outprop01,outprop10,ifreq,igroup
    character*30 flname
    double precision bx,bz,dd
    double precision alpha_coef, beta_coef,depth,width,c1,c2
    ! control data
    integer atype, nnode, nelmt, nmatrl, norder, neig, iopt
    double complex,allocatable :: CS_equi(:),CP_equi(:),CS_equi_true(:),CP_equi_true(:)

    double precision,allocatable :: CS(:), CSest(:,:), covM(:,:), covUz(:),lb(:),ub(:),group_frequency(:,:)
    double precision rhos,nuy,beta
    ! node data
    double precision,allocatable :: x(:), z(:)

    ! element data
    integer,allocatable :: enode(:,:), matrl(:), etype(:)
    double precision,allocatable :: freq_total(:)
    double complex,allocatable:: epara(:,:),epara_true(:,:)
    double precision,allocatable :: Pext(:,:), Phis(:),dobs(:), sigma_M2(:,:)
    double precision Pfac
    ! boundary condition data
    integer nbc(3), ndof
    integer,allocatable :: bc(:), id(:,:), onode(:)

    ! earthquake data
    integer nfreq, nSample, jtmp, ktmp, dof,ltmp
    double precision dfreq, dt, d_sigma

    integer itmp, nx, nz, npml, inode, ielmt, nx1, jdof, idof, nxout, nzout, itime, imatrl
    double precision pi, dx, dz, x1, z1, dtmp


    pi = 4.d0 * datan(1.d0)

    print *, 'Input File Name'
    read *, flname

    !	ininp  = 11
    !	open (ininp , file=trim(flname)//'.inp', status='old')

    outdat = 101
    open (outdat, file=trim(flname)//'.dat')

    outhis = 102
    open (outhis, file=trim(flname)//'.his')

    print *, 'Input job type'
    print *, ' [1] Preprocessing, Reference calculation & Main job'
    print *, ' [2] Reference calculation & Main job'
    print *, ' [3] Main job'
    print *, ' [4] Main job with pevious results'
    read  *, iopt

    npml = 7
    !nx = 245+ 2 * npml
    !nz = 79 + npml
    !nx = 322 + 2 * npml
    !nz = 105 + npml

    nx = 100+ 2 * npml
    nz = 100+ npml
    !!
    ngroup =15

    dx=0.005d0
    dz=0.005d0
    !dx = 9.2d0 / (nx - 2 * npml)
    !dz = 3.0d0 / (nz - npml)
    nfreq = 2
    nfreq_total = 30
    dfreq = (30d0-1d0)/(nfreq_total-1)
    dof=2
    nSample = 50

    depth=(nz-npml)*dz
    width=(nx-2*npml)*dx
    allocate(freq_total(nfreq_total))
    allocate(group_frequency(nfreq,ngroup))
    

    do itmp = 1, nfreq_total

        !freq_total(itmp) = 33d0-(itmp-1)*dfreq
        !freq(itmp) = 1d0*dfreq
        freq_total(itmp) = 1d0+(itmp-1)*dfreq
    enddo
    !dd=6d0!(nfreq_total-nfreq)/(ngroup-2) 
    dd=(nfreq_total-nfreq)/(ngroup-2)-1
     
    do igroup=1,ngroup
        group_frequency(:,igroup)=freq_total(1+dd*(igroup-1):nfreq+dd*(igroup-1))
    enddo
 !do igroup=1,ngroup
 !   do ifreq=1,nfreq
 !   
 !       group_frequency(ifreq,ngroup-igroup+1)=freq_total(igroup)+dd*(ifreq-1)
 !   enddo
 !   enddo
    
    nelmt = nx * nz
    nnode = (nx + 1) * (nz + 1)
    nmatrl = (nx - npml*2) * (nz - npml)

    allocate (x(nnode), z(nnode))
    allocate (enode(nelmt,0:9), matrl(nelmt), etype(nelmt), epara(nelmt,2),epara_true(nelmt,2))
    allocate (CS(nmatrl))

    allocate (bc(nnode), id(nnode,2))

    !CS=1.d0
    !do itmp = 1 , 2
    !
    !    jtmp = (nx - npml*2) * (itmp - 1)
    !    CS(jtmp+1:jtmp+nx-npml*2) = 1.5d0
    !
    !enddo
    !do itmp = 1 , 2
    !
    !    jtmp = (nx - npml*2) * (2 + itmp - 1)
    !    CS(jtmp+1:jtmp+nx-npml*2) = 1.4d0
    !
    !enddo
    !do itmp = 1 , 2
    !
    !    jtmp = (nx - npml*2) * (4 + itmp - 1)
    !    CS(jtmp+1:jtmp+nx-npml*2) = 1.3d0
    !
    !enddo
    !do itmp = 1 , 2
    !
    !    jtmp = (nx - npml*2) * (8 + itmp - 1)
    !    CS(jtmp+1:jtmp+nx-npml*2) = 1.2d0
    !
    !enddo
    !
    !CS(:) = 1.d0
    !
    !do itmp = 1 , 4
    !
    !    jtmp = (nx - npml*2) * (itmp - 1)
    !    CS(jtmp+1:jtmp+nx-npml*2) = 1.5d0
    !
    !enddo
    !
    !do itmp = 1 , 16
    !
    !    jtmp = (nx - npml*2) * (4 + itmp - 1)
    !    CS(jtmp+1:jtmp+nx-npml*2) = 1.3d0
    !
    !enddo
    !
    !do itmp = 1 , 8
    !
    !    jtmp = (nx - npml*2) * (20 + itmp - 1)
    !    CS(jtmp+1:jtmp+nx-npml*2) = 1.0d0
    !
    !enddo
    !
    !do itmp = 5 , 8
    !
    !    do ktmp	= 11 , 20
    !
    !        jtmp = (nx - npml*2) * (20 + itmp - 1)
    !        CS(jtmp+ktmp) = 1.5d0
    !
    !    enddo
    !
    !enddo
    !
    !
    !do itmp = 1 , 12
    !
    !    jtmp = (nx - npml*2) * (28 + itmp - 1)
    !    CS(jtmp+1:jtmp+nx-npml*2) = 1.15d0
    !
    !enddo
    !
    !do itmp = 1 , 4
    !
    !    do ktmp	= 11 , 20
    !
    !        jtmp = (nx - npml*2) * (28 + itmp - 1)
    !        CS(jtmp+ktmp) = 1.5d0
    !
    !    enddo
    !
    !enddo
    !
    CS(:) = 1.d0

    do itmp = 1 , 12

        jtmp = (nx - 2 * npml) * (itmp - 1)
        CS(jtmp+1:jtmp+nx-2*npml) = 1.5d0

    enddo

    do itmp = 1 , 28

        jtmp = (nx - 2 * npml) * (12 + itmp - 1)
        CS(jtmp+1:jtmp+nx-2*npml) = 1.35d0

    enddo

    do itmp = 1 , 20

        jtmp = (nx - 2 * npml) * (40 + itmp - 1)
        CS(jtmp+1:jtmp+nx-2*npml) = 1.15d0

    enddo

    do itmp = 11 , 20

        do ktmp	= 21 , 40

            jtmp = (nx - 2 * npml) * (40 + itmp - 1)
            CS(jtmp+ktmp) = 1.5d0

        enddo

    enddo

    do itmp = 1 , 24

        jtmp = (nx - 2 * npml) * (60 + itmp - 1)
        CS(jtmp+1:jtmp+nx-2*npml) = 1.25d0

    enddo

    do itmp = 1 , 10

        do ktmp	= 21 , 40

            jtmp = (nx - 2 * npml) * (60 + itmp - 1)
            CS(jtmp+ktmp) = 1.5d0

        enddo

    enddo

    do itmp = 1 , 16

        jtmp = (nx - 2 * npml) * (84 + itmp - 1)
        CS(jtmp+1:jtmp+nx-2*npml) = 1.00d0

    enddo

    !open (200, file='marmousi_105x322_Vs.txt')
    !open (200, file='marmousi_79x245_Vs.txt')
    !
    !read (200,*) CS
    !close (200)
    !open (200, file='Marmousi_target.txt')
    !
    !do itmp = 1 , nz - npml
    !
    !    read (200,*) CS((nx-2*npml)*(itmp-1)+1:(nx-2*npml)*itmp)
    !
    !enddo
    !
    !close (200)
    !CS=1d0
    allocate(CSest(nmatrl,2))

    CSest=0d0
    !
    open (200, file='CSguess.txt')
    !open (200, file='CSguess_marmousi.txt')

    do itmp = 1 , nz - npml

        read (200,*) CSest((nx-2*npml)*(itmp-1)+1:(nx-2*npml)*itmp,1)

    enddo
    rhos = 2d0
    nuy = 0.3d0
    beta = 0.05d0


    !do imatrl=1,nmatrl
    !    CS(imatrl)=CS(imatrl)**2d0
    !    CSest(imatrl,1)=CSest(imatrl,1)**2d0
    !enddo
    allocate (CS_equi(nelmt),CP_equi(nelmt),CS_equi_true(nelmt),CP_equi_true(nelmt))
    call PDML_parameters(nz,nx,npml,nmatrl,pi,width,depth,nnode,nelmt,CS_equi,CP_equi,enode,matrl,etype,epara,x,z,dx,dz,nuy,CSest(:,1))
    call PDML_parameters(nz,nx,npml,nmatrl,pi,width,depth,nnode,nelmt,CS_equi_true,CP_equi_true,enode,matrl,etype,epara_true,x,z,dx,dz,nuy,CS)


    bc = 2

    do inode = 1 , nx + 1
        bc(inode) = 3
    enddo

    do inode = 1 , nz + 1

        bc((nx+1)* inode) = 3
        bc((nx+1)* inode-nx) = 3

    enddo

    itmp = 0

    nbc(1) = 0

    do inode = 1 , nnode

        if (bc(inode) == 1) then			! ground nodes on rigid foundation

            itmp = itmp + 2
            nbc(1) = nbc(1) + 2
            id(inode,1) = itmp - 1
            id(inode,2) = itmp

        endif

    enddo

    nbc(2) = 0

    do inode = 1 , nnode

        if (bc(inode) == 2) then			! other ground nodes

            do idof = 1 , dof

                itmp = itmp + 1
                nbc(2) = nbc(2) + 1

                id(inode,idof) = itmp

            enddo

        endif

    enddo

    nbc(3) = 0

    do inode = 1 , nnode

        if (bc(inode) == 3) then			! fixed nodes

            itmp = itmp + 2
            nbc(3) = nbc(3) + 2

            id(inode,1) = itmp - 1
            id(inode,2) = itmp

        endif

    enddo

    nobs = 33!nx - 2*npml + 1-2
    !nobs = 246!marmousi

    allocate ( onode(nobs))
    allocate (lb(nmatrl),ub(nmatrl))
    !
    do inode = 1 , nobs

        itmp = inode

        !onode(itmp) =  (nx + 1) * nz +npml+1+ 1*(itmp-1)
        onode(itmp) =  (nx + 1) * nz +npml+1+ 3*(itmp-1)+2


    enddo
    !nshot=1
    !nshot=50!marmousi
    nshot = 17

    allocate (Pext(nnode*2,nshot))
    allocate( covUz(2*nobs*nfreq*nshot))
    Pext = 0.d0


    do itmp = 1 ,nshot
        !inode=12363
        inode =  (nx + 1) * nz +npml+1+ 6*(itmp-1)+2
        !inode = (nx + 1) * nz +npml+1+ 5*(itmp-1) !marmousi

        write(*,*)  inode


        x1 = (itmp - 1) * dx - nx/2 * dx
        z1 = (jtmp - 1) * dz - npml * dz
        if (bc(inode) == 2) Pext(id(inode,2),itmp) = 1.d0
        !if ((bc(inode) == 1) .and. (itmp == 1)) Pext(id(inode)) = 0.5d0 * Pext(id(inode))

    enddo


    !Pext=Pext/40d0
    write (outdat,'(a)') '/ nx, nz, dx, dz'
    write (outdat,'(2i5,2f10.5,i5,f10.5)') nx, nz, dx, dz

    write (outdat,'(a)') '/ # of node, x, z'

    do inode = 1 , nnode

        write (outdat,'(i7,2f10.5)') inode, x(inode), z(inode)

    enddo

    !write (outdat,'(a)') '/ # of element, enode, matrl, etype, epara'
    !
    !do ielmt = 1 , nelmt
    !
    !    write (outdat,'(8i7,2f10.5)') ielmt, enode(ielmt,0), enode(ielmt,1:enode(ielmt,0)), matrl(ielmt), etype(ielmt), epara(ielmt,1:2)
    !
    !enddo

    !write (outdat,'(a)') '/ # of material, CS, rhos'

    !do imatrl = 1 , nmatrl
    !
    !    write (outdat,'(i5,2f10.5)') imatrl, CS(imatrl)
    !
    !enddo

    write (outdat,'(a)') '/ boundary condition'
    write (outdat,'(i5)') nx+1

    do inode = 1 , nx + 1

        write (outdat,'(3i5)') inode, bc(inode)

    enddo

    write (outdat,'(a)') '/ output node'
    write (outdat,'(i5)') nobs

    do inode = 1 , nobs

        write (outdat,'(3i5)') inode, onode(inode)

    enddo

    write (outdat,'(a)') '/ exciting force'

    if (iopt == 1) then

        print *, 'Preprocessing'

        outout = 24
        open (outout, file=trim(flname)//'.par')

        call MCMCPardiso(nnode, nelmt,  enode,nbc, id, outout)
        !call FiniteElement_PlaneStrain_pardiso(nnode, nelmt, nmatrl, CS, rhos,nuy,beta,x, z, enode, matrl, etype, epara, nbc, id, outout)
        close (outout)
    endif

    !
    !
    !

    if ((iopt == 1) .or. (iopt == 2) .or. (iopt == 4)) then

        print *, 'Reference calculation'

        inpar  = 31
        outobs = 24

        open (inpar , file=trim(flname)//'.par')
        open (outobs, file=trim(flname)//'.obs')

        call MCMCReference(CS_equi_true,nshot,nx, nz, npml, dx, dz, nnode, nelmt, nmatrl, CS, rhos,nuy,beta, x, z, enode, matrl, etype, epara_true, &
            nbc, bc, id, nfreq_total, freq_total, Pext, inpar, outobs, nobs, onode, dobs, 0, Pfac)

        close (inpar)


        Pext(:,:) = Pext(:,:) * Pfac


        inpar  = 31
        open (inpar , file=trim(flname)//'.par')
        call MCMCReference(CS_equi_true,nshot,nx, nz, npml, dx, dz, nnode, nelmt, nmatrl, CS, rhos,nuy,beta, x, z, enode, matrl, etype, epara_true, &
            nbc, bc, id, nfreq_total, freq_total, Pext, inpar, outobs, nobs, onode, dobs, 1, Pfac)


        close (inpar)

        close (outobs)

    endif




    if ((iopt == 1) .or. (iopt == 2) .or. (iopt == 3) .or. (iopt == 4)) then





        allocate( sigma_M2(nmatrl,nmatrl))
        d_sigma = (0.15d0-0.05d0)/(nz-npml-1)

        
        do itmp = 1, nz-npml
            do jtmp=(itmp-1)*(nx-2*npml)+1,itmp*(nx-2*npml)
            sigma_M2( jtmp,jtmp) = ((0.15d0-(itmp-1)*d_sigma))
            enddo
        enddo
        !do imatrl=1,nmatrl
        !    sigma_M2(imatrl,imatrl)=(sigma_M2(imatrl,imatrl)*1.5d0)**2d0
        !enddo
        do imatrl=1,nmatrl
            sigma_M2(imatrl,imatrl)=(0.1*CSest(imatrl,1))**2d0
        enddo
          open (1001 , file='inclusion.priC')
        do imatrl =1,nmatrl
        write (1001 ,'(<nmatrl>f15.5)')  sigma_M2(imatrl,1:nmatrl)
        enddo
        
        sigma_M2=0.5d0*sigma_M2
        !do itmp = 1 , nz - npml
        !    do jtmp = 1 , nx - 2 * npml
        !
        !        do ktmp = 1 , nz - npml
        !            do ltmp = 1 , nx - 2 * npml
        !
        !                sigma_M2((nx-2*npml)*(itmp-1)+jtmp,(nx-2*npml)*(ktmp-1)+ltmp) = (0.3d0)**2d0 * dexp(-iabs(itmp-ktmp)/5d0) * dexp(-iabs(jtmp-ltmp)/5d0)
        !
        !            enddo
        !        enddo
        !
        !    enddo
        !enddo
        
        covUz = 0.d0


        do itmp = 1 , nobs*nfreq*nshot

            covUz(itmp) = (0.002d0*1.d0)**2
            !covUz(itmp) = 0.06

        enddo
        alpha_coef=1d0
        beta_coef=1d0
        lb=0.5d0
        ub=2.5d0
    endif

    if ((iopt == 1) .or. (iopt == 2) .or. (iopt == 3)) then

        print *, 'Main job'

        inpar  = 31
        inobs  = 32
        outest = 24
        outrsp = 25
        outrsl = 30
        outcovM = 40
        outgama=50
        outpost1=60
        outpost0=70
        outprop01=80
        outprop10=90

        open (inpar , file=trim(flname)//'.par')
        open (inobs , file=trim(flname)//'.obs')
        open (outest, file=trim(flname)//'.est')
        open (outrsp, file=trim(flname)//'.prd')
        open (outrsl, file=trim(flname)//'.rsl', form='unformatted')
        open (outcovM, file=trim(flname)//'.covM')
        open (outgama, file=trim(flname)//'.oug')
        !open (outpost1, file=trim(flname)//'.op1')
        !open (outpost0, file=trim(flname)//'.op2')
        !open (outprop01, file=trim(flname)//'.op3')
        !open (outprop10, file=trim(flname)//'.op4')
        !
        call MCMCcore(group_frequency,CS_equi,nshot,nx, nz, npml, dx, dz, nnode, nelmt, nmatrl, CS, rhos,nuy,beta, x, z, enode, matrl, etype, epara, &
            nbc, bc, id, nfreq_total, freq_total,nfreq,ngroup,width,depth, Pext, inpar, inobs,outgama, outpost1, outpost0,outprop01,outprop10, outest, outrsp, outrsl,outcovM, nobs, onode, CSest, sigma_M2,  covUz, nSample,alpha_coef, beta_coef)

        close (inpar)
        close (inobs)
        close (outest)
        close (outrsp)
        close (outrsl)
        close (outcovM)
        !close (outgama)
        !close (outpost1)
        !close (outpost0)
        !close (outprop01)
        !close (outprop10)

    endif

    !
    !
    !
    if (iopt == 4) then

        print *, 'Main job with previous results'

        inpar  = 31
        inobs  = 32
        outest = 24
        outrsp = 25
        outrsl = 30
        outcovM = 40

        open (inpar , file=trim(flname)//'.par')
        open (inobs , file=trim(flname)//'.obs')
        open (outest, file=trim(flname)//'.est')
        open (outrsp, file=trim(flname)//'.prd')
        open (outrsl, file=trim(flname)//'.rsl', form='unformatted')

        read (outrsl) CSest
        read (outrsl) covM

        close (outrsl)

        open (outrsl, file=trim(flname)//'.rsl', form='unformatted')

        call MCMCcore(group_frequency,CS_equi,nshot,nx, nz, npml, dx, dz, nnode, nelmt, nmatrl, CS, rhos,nuy,beta, x, z, enode, matrl, etype, epara, &
            nbc, bc, id, nfreq_total, freq_total,nfreq,ngroup,width,depth, Pext, inpar, inobs,outgama, outpost1, outpost0,outprop01,outprop10, outest, outrsp, outrsl,outcovM, nobs, onode, CSest, sigma_M2, covUz, nSample, alpha_coef, beta_coef)

        close (inpar)
        close (inobs)
        close (outest)
        close (outrsp)
        close (outrsl)
        close (outcovM)
        close (outgama)
        close (outpost1)
        close (outpost0)
        close (outprop01)
        close (outprop10)

    endif

    deallocate (x, z)
    deallocate (enode, matrl, etype, epara)
    deallocate (CS, CSest, covUz, Pext)
    deallocate (bc, id, onode)

    close (outdat)
    close (outhis)

    end program MCMC
