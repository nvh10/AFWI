    subroutine dSdmU_elastic(nshot,onode,nobs,nnode, nelmt, nmatrl, CS, rhos,nuy,beta, x, y, enode, matrl, etype, epara, nbc, id, ia, ja, imatrl, ww, Dis2, v_P2,nnz_P2,pB_P2,pE_P2,ja_P2)
    implicit none
    integer nnz_P2,nshot,nobs,nnode, nelmt, nmatrl, nbc(3), id(nnode,2), ia(nbc(2)+1), ja(ia(nbc(2)+1)-1), imatrl,ixi1,ixi2,iobs
    double precision ww, CS(nmatrl), rhos,nuy,beta
    double precision x(nnode), y(nnode), tx, tz, nx, nz, length, delta, alpha, dt
    double complex Dis2(nbc(2),nshot),xtemp
    integer enode(nelmt,0:9), matrl(nelmt), etype(nelmt),onode(nobs)
    double complex epara(nelmt,2)
    double complex v_P2(nnz_P2,nshot)
    double complex,allocatable:: P2(:,:)
    !double complex,allocatable:: P2(:,:)
    integer pB_P2(nmatrl),pE_P2(nmatrl),ja_P2(nnz_P2)
    integer(kind=8) :: tclock1, tclock2, clock_rate
    real(kind=8) :: elapsed_time
    integer(kind=8) :: tclock11, tclock22, clock_rate2
    real(kind=8) :: elapsed_time2
    ! double precision MSN11(ia(nbc(1)+1)-1), MSN12(nbc(1),nbc(2)), MSN22(nbc(2),nbc(2)), &
    ! DPN11(ia(nbc(1)+1)-1), DPN12(nbc(1),nbc(2)), DPN22(nbc(2),nbc(2)), &
    ! STN11(ia(nbc(1)+1)-1), STN12(nbc(1),nbc(2)), STN22(nbc(2),nbc(2))
    integer nenode, ngp1, ngp2, ielmt, ixi, ieta, inode, idof, jnode, jdof, itmp, jtmp, iia
    double precision pi, N(2,18), shp(9), dshp(9,2), J(2,2), invJ(2,2), detJ, xi1(10), et1(10), wt1(10), xi2(10), et2(10), wt2(10), r
    double precision xlength, ylength, Gxx, Gyy, Gxy
    double complex D(3,3), Drr(2,2), Dss(2,2), Drs(2,2), Dsr(2,2)
    double precision B(3,36), Bp(1,36), Brr(2,36), Bss(2,36), Brs(2,36), Bsr(2,36), xlambda, ylambda
    double complex dMdm(18,18), dCdm(18,18), dKdm(18,18), dRdm(18,18),dStiffness(18,18),a1,a2,a4


    pi = 4.d0 * datan(1.d0)
    !P2=0d0

    !call system_clock(tclock2)
    allocate(P2(nbc(2),nshot))
    do itmp = pB_P2(imatrl), pE_P2(imatrl)-1
        !do i = 1, ncol
        P2(ja_P2(itmp),:)=dcmplx(0.d0,0d0)
        !end do
    end do
    !dStiffness=dStiffness_all(:,:,imatrl)
    nenode = 4
    do ielmt = 1,nelmt
        if (matrl(ielmt) == imatrl) then


            ylength = y(enode(ielmt,4)) - y(enode(ielmt,1))
            xlength = x(enode(ielmt,2)) - x(enode(ielmt,1))
            
            dMdm(:,:) = 0.d0
            dCdm(:,:) = 0.d0
            dKdm(:,:) = 0.d0
            dRdm(:,:) = 0.d0
            
            if (nenode == 4) ngp1 = 2
            if (nenode == 4) ngp2 = 2
            if (nenode == 9) ngp1 = 3
            if (nenode == 9) ngp2 = 3
            
            if (etype(ielmt) == 11) then ! edge CFABC for horizontally propagating waves
            
                if (nenode == 4) ngp1 = 1
                if (nenode == 9) ngp1 = 2
            
            else if (etype(ielmt) == 12) then ! edge CFABC for horizontally evanescent waves
            
                if (nenode == 4) ngp1 = 1
                if (nenode == 9) ngp1 = 2
            
            else if (etype(ielmt) == 21) then ! edge CFABC for vertically propagating waves
            
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp2 = 2
            
            else if (etype(ielmt) == 22) then ! edge CFABC for vertically evanescent waves
            
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp2 = 2
            
            else if (etype(ielmt) == 31) then ! corner CFABC: 11 * 21
            
                if (nenode == 4) ngp1 = 1
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp1 = 2
                if (nenode == 9) ngp2 = 2
            
            else if (etype(ielmt) == 32) then ! corner CFABC: 11 * 22
            
                if (nenode == 4) ngp1 = 1
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp1 = 2
                if (nenode == 9) ngp2 = 2
            
            else if (etype(ielmt) == 33) then ! corner CFABC: 12 * 21
            
                if (nenode == 4) ngp1 = 1
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp1 = 2
                if (nenode == 9) ngp2 = 2
            
            else if (etype(ielmt) == 34) then ! corner CFABC: 12 * 22
            
                if (nenode == 4) ngp1 = 1
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp1 = 2
                if (nenode == 9) ngp2 = 2
            
            endif
            
            call GaussPt(ngp1, xi1(1:ngp1), wt1(1:ngp1))
            call GaussPt(ngp2, xi2(1:ngp2), wt2(1:ngp2))
            
            
            
            do ixi1 = 1 , ngp1
            
                do ixi2 = 1 , ngp2
            
                    call shpfnc(nenode, shp(1:nenode), dshp(1:nenode,1:2), xi1(ixi1), xi2(ixi2))
                    call CalcJ(nenode, J(:,:), xi1(ixi1), xi2(ixi2), x(enode(ielmt,1:nenode)), y(enode(ielmt,1:nenode)))
                    detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1)
                    invJ(1,1) = J(2,2)
                    invJ(2,2) = J(1,1)
                    invJ(1,2) = -J(1,2)
                    invJ(2,1) = -J(2,1)
                    invJ(:,:) = invJ(:,:) / detJ
            
                    B(:,:) = 0.d0
            
                    Brr(:,:) = 0.d0
                    Bss(:,:) = 0.d0
                    Brs(:,:) = 0.d0
                    Bsr(:,:) = 0.d0
            
                    do itmp = 1 , nenode
            
                        B(1,2*itmp-1) = dshp(itmp,1) * invJ(1,1) + dshp(itmp,2) * invJ(1,2) ! dN / dx
                        B(2,2*itmp ) = dshp(itmp,1) * invJ(2,1) + dshp(itmp,2) * invJ(2,2) ! dN / dz
                        B(3,2*itmp-1) = dshp(itmp,1) * invJ(2,1) + dshp(itmp,2) * invJ(2,2)
                        B(3,2*itmp ) = dshp(itmp,1) * invJ(1,1) + dshp(itmp,2) * invJ(1,2)
            
                        Brr(1,2*itmp-1) = dshp(itmp,1) ! dN / dr
                        Brr(2,2*itmp ) = dshp(itmp,1)
            
                        Bss(1,2*itmp-1) = dshp(itmp,2) ! dN / ds
                        Bss(2,2*itmp ) = dshp(itmp,2)
            
                        Brs(1,2*itmp-1) = dshp(itmp,1)
                        Brs(2,2*itmp ) = dshp(itmp,2)
            
                        Bsr(1,2*itmp-1) = dshp(itmp,2)
                        Bsr(2,2*itmp ) = dshp(itmp,1)
            
                    enddo
            
                    D=0d0
                D(1,1)=(2d0*CS(matrl(ielmt)))*(2*rhos*(nuy - 1))/(2*nuy - 1)
                D(1,2)= (2d0*CS(matrl(ielmt)))*(-(2*nuy*rhos)/(2*nuy - 1))
                D(2,1)=D(1,2)
                D(2,2)=(2d0*CS(matrl(ielmt)))*(2*rhos*(nuy - 1))/(2*nuy - 1)
                D(3,3)=(2d0*CS(matrl(ielmt)))*rhos
                D=D*dcmplx(1.d0,2.d0*beta)
            
                    !
                    !Drr(:,:) = 0.d0
                    !Drr(1,1) = (CS(imatrl)*2d0)*(2*rhos*(nuy - 1))/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    !Drr(2,2) = (CS(imatrl)*2d0)*rhos * dcmplx(1.d0,2.d0*beta)
                    !
                    !Dss(:,:) = 0.d0
                    !Dss(1,1) =(CS(imatrl)*2d0)*rhos * dcmplx(1.d0,2.d0*beta)
                    !Dss(2,2) = (CS(imatrl)*2d0)*(2*rhos*(nuy - 1))/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    !
                    !Drs(:,:) = 0.d0
                    !Drs(1,2) = -(8d0*CS(imatrl)*nuy*rhos)/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    !Drs(2,1) = -(8d0*CS(imatrl)*nuy*rhos)/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    !
                    !Dsr(:,:) = 0.d0
                    !Dsr(1,2) = (CS(imatrl)*2d0)*rhos * dcmplx(1.d0,2.d0*beta)
                    !Dsr(2,1) = (CS(imatrl)*2d0)*rhos* dcmplx(1.d0,2.d0*beta)
            
                    if (etype(ielmt) == 0) then
            
                        dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                            + matmul(matmul(transpose(B(1:3,1:2*nenode)), D(1:3,1:3)), B(1:3,1:2*nenode)) * detJ * wt1(ixi1) * wt2(ixi2)
            
                        !else if (etype(ielmt) == 11) then
                        !
                        !    dCdm(1:2*nenode,1:2*nenode) = dCdm(1:2*nenode,1:2*nenode) &
                        !        + ylength/2.d0 / epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dRdm(1:2*nenode,1:2*nenode) = dRdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,1) / (ylength/2.d0) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 12) then
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + ylength/2.d0 * epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + 1.d0 / (ylength/2.d0) / epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 21) then
                        !
                        !    dCdm(1:2*nenode,1:2*nenode) = dCdm(1:2*nenode,1:2*nenode) &
                        !        + xlength/2.d0 / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dRdm(1:2*nenode,1:2*nenode) = dRdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,2) / (xlength/2.d0) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 22) then
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + 1.d0 / (xlength/2.d0) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + xlength/2.d0 * epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 31) then
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 32) then
                        !
                        !    dCdm(1:2*nenode,1:2*nenode) = dCdm(1:2*nenode,1:2*nenode) &
                        !        + 1.d0 / epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dRdm(1:2*nenode,1:2*nenode) = dRdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 33) then
                        !
                        !    dCdm(1:2*nenode,1:2*nenode) = dCdm(1:2*nenode,1:2*nenode) &
                        !        + 1.d0 / epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dRdm(1:2*nenode,1:2*nenode) = dRdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 34) then
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
            
                    endif
            
                enddo
            
            enddo
            
            dStiffness= dKdm


          
            !call system_clock(tclock2)
            !!$omp parallel do
            do inode = 1 , nenode!enode(ielmt,0)
                do idof=1,2
                    itmp = id(enode(ielmt,inode),idof)
                    if ((itmp <= nbc(2))) then
                        do jnode = 1 , nenode!enode(ielmt,0)
                            do jdof=1,2
                                jtmp = id(enode(ielmt,jnode),jdof)
                                if ((jtmp <= nbc(2) )) then

                                    !P2(itmp,:) =  P2(itmp,:)-dStiffness(2*(inode-1)+idof,2*(jnode-1)+jdof) * Dis2(jtmp,:)
                                    P2(itmp,:) =  P2(itmp,:)-dStiffness(2*(inode-1)+idof,2*(jnode-1)+jdof) * Dis2(jtmp,:)
                                endif
                            enddo
                        enddo
                    endif
                enddo
            enddo
        endif
    enddo
    do itmp = pB_P2(imatrl), pE_P2(imatrl)-1
        !do i = 1, ncol
        v_P2(itmp,:)=P2(ja_P2(itmp),:)
        !end do
    end do
    deallocate(P2)
    end subroutine dSdmU_elastic


    subroutine dSdmU_elastic_pre(nx,nz,npml,nshot,onode,nobs,nnode, nelmt, nmatrl, CS, rhos,nuy,beta, x, y, enode, matrl, etype, epara, nbc, id, ia, ja, imatrl, ww, Dis2, P2,dStiffness_all)
    implicit none
    integer nshot,nobs,nnode, nelmt, nmatrl, nbc(3), id(nnode,2), ia(nbc(2)+1), ja(ia(nbc(2)+1)-1), imatrl,ixi1,ixi2,iobs,nx,nz,npml
    double precision ww, CS(nmatrl), rhos,nuy,beta
    double precision x(nnode), y(nnode), tx, tz, length, delta, alpha, dt
    double complex Dis2(nbc(2),nshot), P2(nbc(2),nshot),xtemp
    integer enode(nelmt,0:9), matrl(nelmt), etype(nelmt),onode(nobs)
    double complex epara(nelmt,2)

    integer(kind=8) :: tclock1, tclock2, clock_rate
    real(kind=8) :: elapsed_time
    integer(kind=8) :: tclock11, tclock22, clock_rate2
    real(kind=8) :: elapsed_time2
    ! double precision MSN11(ia(nbc(1)+1)-1), MSN12(nbc(1),nbc(2)), MSN22(nbc(2),nbc(2)), &
    ! DPN11(ia(nbc(1)+1)-1), DPN12(nbc(1),nbc(2)), DPN22(nbc(2),nbc(2)), &
    ! STN11(ia(nbc(1)+1)-1), STN12(nbc(1),nbc(2)), STN22(nbc(2),nbc(2))
    integer nenode, ngp1, ngp2, ielmt, ixi, ieta, inode, idof, jnode, jdof, itmp, jtmp, iia
    double precision pi, N(2,18), shp(9), dshp(9,2), J(2,2), invJ(2,2), detJ, xi1(10), et1(10), wt1(10), xi2(10), et2(10), wt2(10), r
    double precision xlength, ylength, Gxx, Gyy, Gxy
    double complex D(3,3), Drr(2,2), Dss(2,2), Drs(2,2), Dsr(2,2)
    double precision B(3,36), Bp(1,36), Brr(2,36), Bss(2,36), Brs(2,36), Bsr(2,36), xlambda, ylambda,CS_surf((nx - npml*2) * (nz - npml))
    double complex dMdm(18,18), dCdm(18,18), dKdm(18,18), dRdm(18,18),dStiffness(18,18),a1,a2,a4,dStiffness_all(8,8,nmatrl)


    pi = 4.d0 * datan(1.d0)
    !P2=0d0
    a1=dcmplx(ww**2.d0,0d0)
    a2= dcmplx(0.d0,ww)
    a4=dcmplx(1d0,0d0)/dcmplx(0.d0,ww)
    do itmp=1,nbc(2)
        Dis2(itmp,1)=dcmplx(itmp**1.5,itmp**-1.5)
    enddo
    !call system_clock(tclock2)
    do ielmt = 1,nelmt
        if (matrl(ielmt) == imatrl) then
            !call system_clock(tclock2)
            nenode = enode(ielmt,0)
            ylength = y(enode(ielmt,4)) - y(enode(ielmt,1))
            xlength = x(enode(ielmt,2)) - x(enode(ielmt,1))

            dMdm(:,:) = 0.d0
            dCdm(:,:) = 0.d0
            dKdm(:,:) = 0.d0
            dRdm(:,:) = 0.d0

            if (nenode == 4) ngp1 = 2
            if (nenode == 4) ngp2 = 2
            if (nenode == 9) ngp1 = 3
            if (nenode == 9) ngp2 = 3

            if (etype(ielmt) == 11) then ! edge CFABC for horizontally propagating waves

                if (nenode == 4) ngp1 = 1
                if (nenode == 9) ngp1 = 2

            else if (etype(ielmt) == 12) then ! edge CFABC for horizontally evanescent waves

                if (nenode == 4) ngp1 = 1
                if (nenode == 9) ngp1 = 2

            else if (etype(ielmt) == 21) then ! edge CFABC for vertically propagating waves

                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp2 = 2

            else if (etype(ielmt) == 22) then ! edge CFABC for vertically evanescent waves

                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp2 = 2

            else if (etype(ielmt) == 31) then ! corner CFABC: 11 * 21

                if (nenode == 4) ngp1 = 1
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp1 = 2
                if (nenode == 9) ngp2 = 2

            else if (etype(ielmt) == 32) then ! corner CFABC: 11 * 22

                if (nenode == 4) ngp1 = 1
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp1 = 2
                if (nenode == 9) ngp2 = 2

            else if (etype(ielmt) == 33) then ! corner CFABC: 12 * 21

                if (nenode == 4) ngp1 = 1
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp1 = 2
                if (nenode == 9) ngp2 = 2

            else if (etype(ielmt) == 34) then ! corner CFABC: 12 * 22

                if (nenode == 4) ngp1 = 1
                if (nenode == 4) ngp2 = 1
                if (nenode == 9) ngp1 = 2
                if (nenode == 9) ngp2 = 2

            endif

            call GaussPt(ngp1, xi1(1:ngp1), wt1(1:ngp1))
            call GaussPt(ngp2, xi2(1:ngp2), wt2(1:ngp2))



            do ixi1 = 1 , ngp1

                do ixi2 = 1 , ngp2

                    call shpfnc(nenode, shp(1:nenode), dshp(1:nenode,1:2), xi1(ixi1), xi2(ixi2))
                    call CalcJ(nenode, J(:,:), xi1(ixi1), xi2(ixi2), x(enode(ielmt,1:nenode)), y(enode(ielmt,1:nenode)))
                    detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1)
                    invJ(1,1) = J(2,2)
                    invJ(2,2) = J(1,1)
                    invJ(1,2) = -J(1,2)
                    invJ(2,1) = -J(2,1)
                    invJ(:,:) = invJ(:,:) / detJ

                    B(:,:) = 0.d0

                    Brr(:,:) = 0.d0
                    Bss(:,:) = 0.d0
                    Brs(:,:) = 0.d0
                    Bsr(:,:) = 0.d0

                    do itmp = 1 , nenode

                        B(1,2*itmp-1) = dshp(itmp,1) * invJ(1,1) + dshp(itmp,2) * invJ(1,2) ! dN / dx
                        B(2,2*itmp ) = dshp(itmp,1) * invJ(2,1) + dshp(itmp,2) * invJ(2,2) ! dN / dz
                        B(3,2*itmp-1) = dshp(itmp,1) * invJ(2,1) + dshp(itmp,2) * invJ(2,2)
                        B(3,2*itmp ) = dshp(itmp,1) * invJ(1,1) + dshp(itmp,2) * invJ(1,2)

                        Brr(1,2*itmp-1) = dshp(itmp,1) ! dN / dr
                        Brr(2,2*itmp ) = dshp(itmp,1)

                        Bss(1,2*itmp-1) = dshp(itmp,2) ! dN / ds
                        Bss(2,2*itmp ) = dshp(itmp,2)

                        Brs(1,2*itmp-1) = dshp(itmp,1)
                        Brs(2,2*itmp ) = dshp(itmp,2)

                        Bsr(1,2*itmp-1) = dshp(itmp,2)
                        Bsr(2,2*itmp ) = dshp(itmp,1)

                    enddo

                    D=0d0
                    D(1,1)=(2*rhos*(nuy - 1))/(2*nuy - 1)
                    D(1,2)= (-(2*nuy*rhos)/(2*nuy - 1))
                    D(2,1)=D(1,2)
                    D(2,2)=(2*rhos*(nuy - 1))/(2*nuy - 1)
                    D(3,3)=rhos
                    D=D*dcmplx(1.d0,2.d0*beta)

                    !
                    !Drr(:,:) = 0.d0
                    !Drr(1,1) = (CS_equi*2d0)*(2*rhos*(nuy - 1))/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    !Drr(2,2) = (CS(imatrl)*2d0)*rhos * dcmplx(1.d0,2.d0*beta)
                    !
                    !Dss(:,:) = 0.d0
                    !Dss(1,1) =(CS(imatrl)*2d0)*rhos * dcmplx(1.d0,2.d0*beta)
                    !Dss(2,2) = (CS(imatrl)*2d0)*(2*rhos*(nuy - 1))/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    !
                    !Drs(:,:) = 0.d0
                    !Drs(1,2) = -(8d0*CS(imatrl)*nuy*rhos)/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    !Drs(2,1) = -(8d0*CS(imatrl)*nuy*rhos)/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    !
                    !Dsr(:,:) = 0.d0
                    !Dsr(1,2) = (CS(imatrl)*2d0)*rhos * dcmplx(1.d0,2.d0*beta)
                    !Dsr(2,1) = (CS(imatrl)*2d0)*rhos* dcmplx(1.d0,2.d0*beta)

                    if (etype(ielmt) == 0) then

                        dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                            + matmul(matmul(transpose(B(1:3,1:2*nenode)), D(1:3,1:3)), B(1:3,1:2*nenode)) * detJ * wt1(ixi1) * wt2(ixi2)

                        !else if (etype(ielmt) == 11) then
                        !
                        !    dCdm(1:2*nenode,1:2*nenode) = dCdm(1:2*nenode,1:2*nenode) &
                        !        + ylength/2.d0 / epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dRdm(1:2*nenode,1:2*nenode) = dRdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,1) / (ylength/2.d0) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 12) then
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + ylength/2.d0 * epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + 1.d0 / (ylength/2.d0) / epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 21) then
                        !
                        !    dCdm(1:2*nenode,1:2*nenode) = dCdm(1:2*nenode,1:2*nenode) &
                        !        + xlength/2.d0 / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dRdm(1:2*nenode,1:2*nenode) = dRdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,2) / (xlength/2.d0) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 22) then
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + 1.d0 / (xlength/2.d0) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + xlength/2.d0 * epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 31) then
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 32) then
                        !
                        !    dCdm(1:2*nenode,1:2*nenode) = dCdm(1:2*nenode,1:2*nenode) &
                        !        + 1.d0 / epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dRdm(1:2*nenode,1:2*nenode) = dRdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 33) then
                        !
                        !    dCdm(1:2*nenode,1:2*nenode) = dCdm(1:2*nenode,1:2*nenode) &
                        !        + 1.d0 / epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !    dRdm(1:2*nenode,1:2*nenode) = dRdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)
                        !
                        !else if (etype(ielmt) == 34) then
                        !
                        !    dKdm(1:2*nenode,1:2*nenode) = dKdm(1:2*nenode,1:2*nenode) &
                        !        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        !        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    endif

                enddo

            enddo
            !     delmtM_all(:,:,ielmt)=dMdm(:,:)
            !delmtK_all(:,:,ielmt)=dKdm(:,:)
            !delmtC_all(:,:,ielmt)=dCdm(:,:)
            !delmtR_all(:,:,ielmt)=dRdm(:,:)
            dStiffness=dKdm
            dStiffness_all(1:8,1:8,imatrl)=dStiffness(1:8,1:8)
        
            !call system_clock(tclock22, clock_rate2)
            !elapsed_time2 = float(tclock22 - tclock2) / float(clock_rate2)
            !write(*,*) "dStiffness", elapsed_time2



            !call system_clock(tclock2)
            !!$omp parallel do
            do inode = 1 , nenode!enode(ielmt,0)
                do idof=1,2
                    itmp = id(enode(ielmt,inode),idof)
                    if ((itmp <= nbc(2))) then
                        do jnode = 1 , nenode!enode(ielmt,0)
                            do jdof=1,2
                                jtmp = id(enode(ielmt,jnode),jdof)
                                if ((jtmp <= nbc(2) )) then

                                    P2(itmp,:) =  P2(itmp,:)-dStiffness(2*(inode-1)+idof,2*(jnode-1)+jdof) * Dis2(jtmp,:)
                                    !P2(itmp) = P2(itmp) - (dcmplx(ww**2.d0,0d0)*dMdm(2*(inode-1)+idof,2*(jnode-1)+jdof) + dcmplx(0.d0,ww) * dCdm(2*(inode-1)+idof,2*(jnode-1)+jdof) + dKdm(2*(inode-1)+idof,2*(jnode-1)+jdof) + dcmplx(1d0,0d0)/dcmplx(0.d0,ww)*dRdm(2*(inode-1)+idof,2*(jnode-1)+jdof)) * Dis2(jtmp)
                                endif
                            enddo
                        enddo
                    endif
                enddo
            enddo
        endif
    enddo
    
    end subroutine dSdmU_elastic_pre