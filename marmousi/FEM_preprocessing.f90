    subroutine FEM_preprocessing(nnode, nelmt, nmatrl, rhos,nuy,beta,x, y, enode, matrl, etype, epara, nbc, id, elmtM_all, elmtC_all, elmtK_all, elmtR_all, elmtM_all_n, elmtC_all_n, elmtK_all_n, elmtR_all_n)
    implicit none
    integer nnode, nelmt, nmatrl, nbc(3), id(nnode,2)
    double precision rhos(nmatrl), nuy(nmatrl),beta(nmatrl),rho(2,2)
    double precision x(nnode), y(nnode), tx, tz, nx, nz, length
    integer enode(nelmt,0:9), matrl(nelmt), etype(nelmt)
    double precision epara(nelmt,2)



    !double complex,allocatable :: MSN22(:), DPN22(:), STN22(:), RRN22(:),DSN22(:)
    integer nenode, ngp1, ngp2, ielmt, ixi, ieta, inode, idof, jnode, jdof, itmp, jtmp, iia,ixi1, ixi2
    double precision pi, N(2,18), shp(9), dshp(9,2), J(2,2), invJ(2,2), detJ, xi1(10), et1(10), wt1(10), xi2(10), et2(10), wt2(10), r
    !double precision xlength, ylength, Gxx, Gyy, Gxy, D(3,3), Dr(2,2), Ds(2,2), B(3,36), Br(2,36), Bs(2,36)
    double precision xlength, ylength, Gxx, Gyy, Gxy
    double complex D(3,3), Drr(2,2), Dss(2,2), Drs(2,2), Dsr(2,2)
    double precision B(3,36), Bp(1,36), Brr(2,36), Bss(2,36), Brs(2,36), Bsr(2,36), xlambda, ylambda
    double complex elmtM_all(18,18,nelmt), elmtC_all(18,18,nelmt), elmtK_all(18,18,nelmt), elmtR_all(18,18,nelmt)
    double complex elmtM_all_n(18,18,nelmt), elmtC_all_n(18,18,nelmt), elmtK_all_n(18,18,nelmt), elmtR_all_n(18,18,nelmt)
    elmtM_all=0d0
    elmtC_all=0d0
    elmtK_all=0d0
    elmtR_all=0d0
    elmtM_all_n=0d0
    elmtC_all_n=0d0
    elmtK_all_n=0d0
    elmtR_all_n=0d0

    pi = 4.d0 * datan(1.d0)

    do ielmt = 1,nelmt

        nenode = enode(ielmt,0)
        ylength = y(enode(ielmt,4)) - y(enode(ielmt,1))
        xlength = x(enode(ielmt,2)) - x(enode(ielmt,1))

        elmtM_all(:,:,ielmt) = 0.d0
        elmtC_all(:,:,ielmt) = 0.d0
        elmtK_all(:,:,ielmt) = 0.d0
        elmtR_all(:,:,ielmt) = 0.d0
        elmtM_all_n(:,:,ielmt) = 0.d0
        elmtC_all_n(:,:,ielmt) = 0.d0
        elmtK_all_n(:,:,ielmt) = 0.d0
        elmtR_all_n(:,:,ielmt) = 0.d0

        if (nenode == 4) ngp1 = 2
        if (nenode == 4) ngp2 = 2
        if (nenode == 9) ngp1 = 3
        if (nenode == 9) ngp2 = 3

        if (etype(ielmt) == 11) then		! edge CFABC for horizontally propagating waves

            if (nenode == 4) ngp1 = 1
            if (nenode == 9) ngp1 = 2

        else if (etype(ielmt) == 12) then	! edge CFABC for horizontally evanescent waves

            if (nenode == 4) ngp1 = 1
            if (nenode == 9) ngp1 = 2

        else if (etype(ielmt) == 21) then	! edge CFABC for vertically propagating waves

            if (nenode == 4) ngp2 = 1
            if (nenode == 9) ngp2 = 2

        else if (etype(ielmt) == 22) then	! edge CFABC for vertically evanescent waves

            if (nenode == 4) ngp2 = 1
            if (nenode == 9) ngp2 = 2

        else if (etype(ielmt) == 31) then	! corner CFABC: 11 * 21

            if (nenode == 4) ngp1 = 1
            if (nenode == 4) ngp2 = 1
            if (nenode == 9) ngp1 = 2
            if (nenode == 9) ngp2 = 2

        else if (etype(ielmt) == 32) then	! corner CFABC: 11 * 22

            if (nenode == 4) ngp1 = 1
            if (nenode == 4) ngp2 = 1
            if (nenode == 9) ngp1 = 2
            if (nenode == 9) ngp2 = 2

        else if (etype(ielmt) == 33) then	! corner CFABC: 12 * 21

            if (nenode == 4) ngp1 = 1
            if (nenode == 4) ngp2 = 1
            if (nenode == 9) ngp1 = 2
            if (nenode == 9) ngp2 = 2

        else if (etype(ielmt) == 34) then	! corner CFABC: 12 * 22

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
                invJ(1,1) =  J(2,2)
                invJ(2,2) =  J(1,1)
                invJ(1,2) = -J(1,2)
                invJ(2,1) = -J(2,1)
                invJ(:,:) = invJ(:,:) / detJ

                N(:,:) = 0.d0

                do itmp = 1 , nenode

                    N(1,2*itmp-1) = shp(itmp)
                    N(2,2*itmp  ) = shp(itmp)

                enddo

                rho(:,:) = 0.d0
                rho(1,1) = rhos(matrl(ielmt))
                rho(2,2) = rhos(matrl(ielmt))

                if (etype(ielmt) == 0) then

                    elmtM_all_n(1:2*nenode,1:2*nenode,ielmt) = elmtM_all_n(1:2*nenode,1:2*nenode,ielmt) &
                        + matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * detJ * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 11) then

                    elmtC_all_n(1:2*nenode,1:2*nenode,ielmt) = elmtC_all_n(1:2*nenode,1:2*nenode,ielmt) &
                        + ylength/2.d0 * epara(ielmt,1) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 12) then

                    elmtM_all_n(1:2*nenode,1:2*nenode,ielmt) = elmtM_all_n(1:2*nenode,1:2*nenode,ielmt) &
                        + ylength/2.d0 / epara(ielmt,1) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 21) then

                    elmtC_all_n(1:2*nenode,1:2*nenode,ielmt) = elmtC_all_n(1:2*nenode,1:2*nenode,ielmt) &
                        + xlength/2.d0 * epara(ielmt,1) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 22) then

                    elmtM_all_n(1:2*nenode,1:2*nenode,ielmt) = elmtM_all_n(1:2*nenode,1:2*nenode,ielmt) &
                        + xlength/2.d0 / epara(ielmt,1) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 31) then

                    elmtK_all_n(1:2*nenode,1:2*nenode,ielmt) = elmtK_all_n(1:2*nenode,1:2*nenode,ielmt) &
                        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 32) then

                    elmtC_all_n(1:2*nenode,1:2*nenode,ielmt) = elmtC_all_n(1:2*nenode,1:2*nenode,ielmt) &
                        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 33) then

                    elmtC_all_n(1:2*nenode,1:2*nenode,ielmt) = elmtC_all_n(1:2*nenode,1:2*nenode,ielmt) &
                        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 34) then

                    elmtM_all_n(1:2*nenode,1:2*nenode,ielmt) = elmtM_all_n(1:2*nenode,1:2*nenode,ielmt) &
                        + 1.d0 / epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                endif

            enddo

        enddo

        do ixi1 = 1 , ngp1

            do ixi2 = 1 , ngp2

                call shpfnc(nenode, shp(1:nenode), dshp(1:nenode,1:2), xi1(ixi1), xi2(ixi2))
                call CalcJ(nenode, J(:,:), xi1(ixi1), xi2(ixi2), x(enode(ielmt,1:nenode)), y(enode(ielmt,1:nenode)))
                detJ = J(1,1) * J(2,2) - J(1,2) * J(2,1)
                invJ(1,1) =  J(2,2)
                invJ(2,2) =  J(1,1)
                invJ(1,2) = -J(1,2)
                invJ(2,1) = -J(2,1)
                invJ(:,:) = invJ(:,:) / detJ

                B(:,:) = 0.d0

                Brr(:,:) = 0.d0
                Bss(:,:) = 0.d0
                Brs(:,:) = 0.d0
                Bsr(:,:) = 0.d0

                do itmp = 1 , nenode

                    B(1,2*itmp-1) = dshp(itmp,1) * invJ(1,1) + dshp(itmp,2) * invJ(1,2)		! dN / dx
                    B(2,2*itmp  ) = dshp(itmp,1) * invJ(2,1) + dshp(itmp,2) * invJ(2,2)		! dN / dz
                    B(3,2*itmp-1) = dshp(itmp,1) * invJ(2,1) + dshp(itmp,2) * invJ(2,2)
                    B(3,2*itmp  ) = dshp(itmp,1) * invJ(1,1) + dshp(itmp,2) * invJ(1,2)

                    Brr(1,2*itmp-1) = dshp(itmp,1)		! dN / dr
                    Brr(2,2*itmp  ) = dshp(itmp,1)

                    Bss(1,2*itmp-1) = dshp(itmp,2)		! dN / ds
                    Bss(2,2*itmp  ) = dshp(itmp,2)

                    Brs(1,2*itmp-1) = dshp(itmp,1)
                    Brs(2,2*itmp  ) = dshp(itmp,2)

                    Bsr(1,2*itmp-1) = dshp(itmp,2)
                    Bsr(2,2*itmp  ) = dshp(itmp,1)

                enddo

                !D(:,:) = 0.d0
                !D(1,1) = (la(ielmt) + 2.d0*mu(ielmt)) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                !D(1,2) = la(ielmt) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                !D(2,1) = la(ielmt) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                !D(2,2) = (la(ielmt) + 2.d0*mu(ielmt)) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                !D(3,3) = mu(ielmt) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                D=0d0
                D(1,1)=(2*rhos(matrl(ielmt))*(nuy(matrl(ielmt)) - 1))/(2*nuy(matrl(ielmt)) - 1)
                D(1,2)= (-(2*nuy(matrl(ielmt))*rhos(matrl(ielmt)))/(2*nuy(matrl(ielmt)) - 1))
                D(2,1)=D(1,2)
                D(2,2)=(2*rhos(matrl(ielmt))*(nuy(matrl(ielmt)) - 1))/(2*nuy(matrl(ielmt)) - 1)
                D(3,3)=rhos(matrl(ielmt))
                D=D*dcmplx(1.d0,2.d0*beta(matrl(ielmt)))

                !Drr(:,:) = 0.d0
                !Drr(1,1) = (la(ielmt) + 2.d0*mu(ielmt)) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                !Drr(2,2) = mu(ielmt) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                !
                !Dss(:,:) = 0.d0
                !Dss(1,1) = mu(ielmt) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                !Dss(2,2) = (la(ielmt) + 2.d0*mu(ielmt)) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                Drr(:,:) = 0.d0
                Drr(1,1) = (2*rhos(matrl(ielmt))*(nuy(matrl(ielmt)) - 1))/(2*nuy(matrl(ielmt)) - 1) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                Drr(2,2) = rhos(matrl(ielmt)) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))

                Dss(:,:) = 0.d0
                Dss(1,1) =rhos(matrl(ielmt)) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                Dss(2,2) =(2*rhos(matrl(ielmt))*(nuy(matrl(ielmt)) - 1))/(2*nuy(matrl(ielmt)) - 1) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))

                Drs(:,:) = 0.d0
                Drs(1,2) = -(2d0*nuy(matrl(ielmt))*rhos(matrl(ielmt)))/(2*nuy(matrl(ielmt)) - 1) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                Drs(2,1) = -(2d0*nuy(matrl(ielmt))*rhos(matrl(ielmt)))/(2*nuy(matrl(ielmt)) - 1) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))

                Dsr(:,:) = 0.d0
                Dsr(1,2) = rhos(matrl(ielmt)) * dcmplx(1.d0,2.d0*beta(matrl(ielmt)))
                Dsr(2,1) = rhos(matrl(ielmt))* dcmplx(1.d0,2.d0*beta(matrl(ielmt)))

                if (etype(ielmt) == 0) then

                    elmtK_all(1:2*nenode,1:2*nenode,ielmt) = elmtK_all(1:2*nenode,1:2*nenode,ielmt) &
                        + matmul(matmul(transpose(B(1:3,1:2*nenode)), D(1:3,1:3)), B(1:3,1:2*nenode)) * detJ * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 11) then

                    elmtC_all(1:2*nenode,1:2*nenode,ielmt) = elmtC_all(1:2*nenode,1:2*nenode,ielmt) &
                        + ylength/2.d0 / epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtK_all(1:2*nenode,1:2*nenode,ielmt) = elmtK_all(1:2*nenode,1:2*nenode,ielmt) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtR_all(1:2*nenode,1:2*nenode,ielmt) = elmtR_all(1:2*nenode,1:2*nenode,ielmt) &
                        + epara(ielmt,1) / (ylength/2.d0) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 12) then

                    elmtK_all(1:2*nenode,1:2*nenode,ielmt) = elmtK_all(1:2*nenode,1:2*nenode,ielmt) &
                        + ylength/2.d0 * epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + 1.d0 / (ylength/2.d0) / epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 21) then

                    elmtC_all(1:2*nenode,1:2*nenode,ielmt) = elmtC_all(1:2*nenode,1:2*nenode,ielmt) &
                        + xlength/2.d0 / epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtK_all(1:2*nenode,1:2*nenode,ielmt) = elmtK_all(1:2*nenode,1:2*nenode,ielmt) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtR_all(1:2*nenode,1:2*nenode,ielmt) = elmtR_all(1:2*nenode,1:2*nenode,ielmt) &
                        + epara(ielmt,1) / (xlength/2.d0) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 22) then

                    elmtK_all(1:2*nenode,1:2*nenode,ielmt) = elmtK_all(1:2*nenode,1:2*nenode,ielmt) &
                        + 1.d0 / (xlength/2.d0) / epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + xlength/2.d0 * epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 31) then

                    elmtK_all(1:2*nenode,1:2*nenode,ielmt) = elmtK_all(1:2*nenode,1:2*nenode,ielmt) &
                        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 32) then

                    elmtC_all(1:2*nenode,1:2*nenode,ielmt) = elmtC_all(1:2*nenode,1:2*nenode,ielmt) &
                        + 1.d0 / epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtK_all(1:2*nenode,1:2*nenode,ielmt) = elmtK_all(1:2*nenode,1:2*nenode,ielmt) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtR_all(1:2*nenode,1:2*nenode,ielmt) = elmtR_all(1:2*nenode,1:2*nenode,ielmt) &
                        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 33) then

                    elmtC_all(1:2*nenode,1:2*nenode,ielmt) = elmtC_all(1:2*nenode,1:2*nenode,ielmt) &
                        + 1.d0 / epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtK_all(1:2*nenode,1:2*nenode,ielmt) = elmtK_all(1:2*nenode,1:2*nenode,ielmt) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtR_all(1:2*nenode,1:2*nenode,ielmt) = elmtR_all(1:2*nenode,1:2*nenode,ielmt) &
                        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 34) then

                    elmtK_all(1:2*nenode,1:2*nenode,ielmt) = elmtK_all(1:2*nenode,1:2*nenode,ielmt) &
                        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                endif

            enddo

        enddo

    enddo

    end subroutine FEM_preprocessing

    !!!!!!!!!!!!!!!!

    subroutine FEM_assemble(nnode, nelmt, nmatrl,CS, rhos,nuy,beta,x, y, enode, matrl, etype, epara, nbc, id,ia,ja, elmtM_all, elmtC_all, elmtK_all, elmtR_all, elmtM_all_n, elmtC_all_n, elmtK_all_n, elmtR_all_n,MSN22,DPN22,STN22,RRN22)
    implicit none
    integer nnode, nelmt, nmatrl, nbc(3), id(nnode,2), ia(nbc(2)+1),ja(ia(nbc(2)+1)-1)
    double precision CS(nmatrl),rhos(nmatrl), nuy(nmatrl),beta(nmatrl),rho(2,2)
    double precision x(nnode), y(nnode), tx, tz, nx, nz, length
    integer enode(nelmt,0:9), matrl(nelmt), etype(nelmt)
    double precision epara(nelmt,2)
    double complex MSN22(ia(nbc(2)+1)-1), DPN22(ia(nbc(2)+1)-1), STN22(ia(nbc(2)+1)-1), RRN22(ia(nbc(2)+1)-1), DSN22(ia(nbc(2)+1)-1)

    integer nenode, ngp1, ngp2, ielmt, ixi, ieta, inode, idof, jnode, jdof, itmp, jtmp, iia,ixi1, ixi2
    double precision pi
    double precision xlength, ylength

    double complex elmtM_all(18,18,nelmt), elmtC_all(18,18,nelmt), elmtK_all(18,18,nelmt), elmtR_all(18,18,nelmt)
    double complex elmtM_all_n(18,18,nelmt), elmtC_all_n(18,18,nelmt), elmtK_all_n(18,18,nelmt), elmtR_all_n(18,18,nelmt)

    double complex elmtM(18,18), elmtC(18,18), elmtK(18,18), elmtR(18,18)

    pi = 4.d0 * datan(1.d0)
    do ielmt = 1,nelmt
        elmtM=elmtM_all_n(:,:,ielmt)+(CS(matrl(ielmt))**2d0)*elmtM_all(:,:,ielmt)
        elmtK=elmtK_all_n(:,:,ielmt)+(CS(matrl(ielmt))**2d0)*elmtK_all(:,:,ielmt)
        elmtC=elmtC_all_n(:,:,ielmt)+(CS(matrl(ielmt))**2d0)*elmtC_all(:,:,ielmt)
        elmtR=elmtR_all_n(:,:,ielmt)+(CS(matrl(ielmt))**2d0)*elmtR_all(:,:,ielmt)


        do inode = 1 , enode(ielmt,0)
            do idof=1,2

                itmp = id(enode(ielmt,inode),idof)

                do jnode = 1 , enode(ielmt,0)
                    do jdof=1,2

                        jtmp = id(enode(ielmt,jnode),jdof)

                        if ((itmp > nbc(1)) .and. (itmp <= nbc(1)+nbc(2)) .and. &
                            (jtmp > nbc(1)) .and. (jtmp <= nbc(1)+nbc(2))) then



                            do iia = ia(itmp-nbc(1)) , ia(itmp+1-nbc(1))-1

                                if (ja(iia) == jtmp-nbc(1)) then

                                    MSN22(iia) = MSN22(iia) + elmtM(2*(inode-1)+idof,2*(jnode-1)+jdof)
                                    DPN22(iia) = DPN22(iia) + elmtC(2*(inode-1)+idof,2*(jnode-1)+jdof)
                                    STN22(iia) = STN22(iia) + elmtK(2*(inode-1)+idof,2*(jnode-1)+jdof)
                                    RRN22(iia) = RRN22(iia) + elmtR(2*(inode-1)+idof,2*(jnode-1)+jdof)
                                endif

                            enddo


                        else if ((itmp <= nbc(2)) .and. (jtmp > nbc(2))) then


                        else if ((itmp > nbc(2)) .and. (jtmp > nbc(2))) then

                        endif
                    enddo
                enddo
            enddo
        enddo

    enddo

    end subroutine FEM_assemble


    subroutine dSdmU_elastic_assemble(nnode, nelmt, nmatrl, CS, enode, matrl, nbc, id, ia, ja, imatrl, ww, Dis2, P2,dStiffness,his_n,his_jtmp,his_yy)
    implicit none
    integer nnode, nelmt, nmatrl, nbc(3), id(nnode,2), ia(nbc(2)+1), ja(ia(nbc(2)+1)-1), imatrl,ixi1,ixi2,NN
    double precision ww, CS(nmatrl)
    double complex Dis2(nbc(2)), P2(nbc(2))
    integer enode(nelmt,0:9),matrl(nelmt),his_itmp(64,nelmt),his_jtmp(64,nelmt),his_xx(64,nelmt),his_yy(64,nelmt),his_n(nelmt)
    integer nenode, ngp1, ngp2, ielmt, ixi, ieta, inode, idof, jnode, jdof, itmp, jtmp, iia
    double precision pi
    double complex dStiffness(18,18,nelmt)
    integer(kind=8) :: tclock1, tclock2, clock_rate
    real(kind=8) :: elapsed_time
    integer(kind=8) :: tclock11, tclock22, clock_rate2
    real(kind=8) :: elapsed_time2
    double complex, allocatable :: A(:,:),B(:),C(:)
    !
    do ielmt = 1,nelmt
        if (matrl(ielmt) == imatrl) then
            do inode = 1 , enode(ielmt,0)


                do idof=1,2
                    itmp = id(enode(ielmt,inode),idof)
                    if ((itmp <= nbc(2))) then
                        do jnode = 1 , enode(ielmt,0)
                            do jdof=1,2
                                jtmp = id(enode(ielmt,jnode),jdof)
                                if ((jtmp <= nbc(2) )) then

                                    P2(itmp) =  P2(itmp)-dStiffness(2*(inode-1)+idof,2*(jnode-1)+jdof,ielmt) * Dis2(jtmp)
                                    !P2(itmp) = P2(itmp) - (dcmplx(ww**2.d0,0d0)*dMdm(2*(inode-1)+idof,2*(jnode-1)+jdof) + dcmplx(0.d0,ww) * dCdm(2*(inode-1)+idof,2*(jnode-1)+jdof) + dKdm(2*(inode-1)+idof,2*(jnode-1)+jdof) + dcmplx(1d0,0d0)/dcmplx(0.d0,ww)*dRdm(2*(inode-1)+idof,2*(jnode-1)+jdof)) * Dis2(jtmp)

                                endif

                            enddo
                        enddo

                    endif
                enddo

            enddo
                  endif
            enddo
            
            !do ielmt = 1,nelmt
            !    if (matrl(ielmt) == imatrl) then
            !        NN=his_n(ielmt)**0.5
            !        allocate(A(NN,NN),B(NN),C(NN))
            !        A=dStiffness(his_yy(1:NN,ielmt),his_yy(1:NN,ielmt),ielmt)
            !        B=Dis2(his_jtmp(1:NN,ielmt))
            !        C=P2(his_jtmp(1:NN,ielmt))
            !        call zsymm('L','L',NN,1,dcmplx(-1d0,0d0),A,NN,B,NN,dcmplx(1d0,0d0), C,NN)
            !        P2(his_jtmp(1:NN,ielmt))=C
            !    deallocate(A,B,C)
            !    endif
            !enddo
        end subroutine dSdmU_elastic_assemble
