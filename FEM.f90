    subroutine FiniteElement_PlaneStrain(CS_equi,nnode, nelmt, nmatrl, CS, rhos,nuy,beta, MSN22, DPN22, STN22,RRN22,x, y, enode, matrl, etype, epara, nbc, id, ia, ja)
    implicit none
    integer nnode, nelmt, nmatrl, nbc(3), id(nnode,2), ia(nbc(2)+1),ja(ia(nbc(2)+1)-1)

    double precision ww, CS(nmatrl), rhos, nuy,beta,rho(2,2)
    double precision x(nnode), y(nnode), tx, tz, nx, nz, length
    integer enode(nelmt,0:9), matrl(nelmt), etype(nelmt)
    double complex epara(nelmt,2),CS_equi(nelmt)

    !double complex  MSN23(nbc(2),nbc(3)), MSN33(nbc(3),nbc(3)), &
    !    DPN23(nbc(2),nbc(3)), DPN33(nbc(3),nbc(3)), &
    !    STN23(nbc(2),nbc(3)), STN33(nbc(3),nbc(3)), &
    !    RRN23(nbc(2),nbc(3)), RRN33(nbc(3),nbc(3))
    !double complex MSN22(nbc(2),nbc(2)), DPN22(nbc(2),nbc(2)), STN22(nbc(2),nbc(2)), RRN22(nbc(2),nbc(2)), DSN22(nbc(2),nbc(2))
    double complex MSN22(ia(nbc(2)+1)-1), DPN22(ia(nbc(2)+1)-1), STN22(ia(nbc(2)+1)-1), RRN22(ia(nbc(2)+1)-1), DSN22(ia(nbc(2)+1)-1)
    !double complex,allocatable :: MSN22(:), DPN22(:), STN22(:), RRN22(:),DSN22(:)
    integer nenode, ngp1, ngp2, ielmt, ixi, ieta, inode, idof, jnode, jdof, itmp, jtmp, iia,ixi1, ixi2
    double precision pi, N(2,18), shp(9), dshp(9,2), J(2,2), invJ(2,2), detJ, xi1(10), et1(10), wt1(10), xi2(10), et2(10), wt2(10), r
    !double precision xlength, ylength, Gxx, Gyy, Gxy, D(3,3), Dr(2,2), Ds(2,2), B(3,36), Br(2,36), Bs(2,36)
    double precision xlength, ylength, Gxx, Gyy, Gxy
    double complex D(3,3), Drr(2,2), Dss(2,2), Drs(2,2), Dsr(2,2)
    double precision B(3,36), Bp(1,36), Brr(2,36), Bss(2,36), Brs(2,36), Bsr(2,36), xlambda, ylambda
    double complex elmtM(18,18), elmtC(18,18), elmtK(18,18), elmtR(18,18)
    double complex elmtM_all(18,18,nelmt), elmtC_all(18,18,nelmt), elmtK_all(18,18,nelmt), elmtR_all(18,18,nelmt)
    !double precision elmtM(9,9), elmtC(9,9), elmtK(9,9), elmtR(9,9), elmtq(9,2)
    !allocate ( MSN22(nbc(2),nbc(2)), DPN22(nbc(2),nbc(2)), STN22(nbc(2),nbc(2)), RRN22(nbc(2),nbc(2)), DSN22(nbc(2),nbc(2)))

    pi = 4.d0 * datan(1.d0)
    MSN22=dcmplx(0d0,0d0)
    STN22=dcmplx(0d0,0d0)
    DPN22=dcmplx(0d0,0d0)
    RRN22=dcmplx(0d0,0d0)




    do ielmt = 1,nelmt

        nenode = enode(ielmt,0)
        ylength = y(enode(ielmt,4)) - y(enode(ielmt,1))
        xlength = x(enode(ielmt,2)) - x(enode(ielmt,1))

        elmtM(:,:) = 0.d0
        elmtC(:,:) = 0.d0
        elmtK(:,:) = 0.d0
        elmtR(:,:) = 0.d0

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
                rho(1,1) = rhos
                rho(2,2) = rhos

                if (etype(ielmt) == 0) then

                    elmtM(1:2*nenode,1:2*nenode) = elmtM(1:2*nenode,1:2*nenode) &
                        + matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * detJ * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 11) then

                    elmtC(1:2*nenode,1:2*nenode) = elmtC(1:2*nenode,1:2*nenode) &
                        + ylength/2.d0 * epara(ielmt,1) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 12) then

                    elmtM(1:2*nenode,1:2*nenode) = elmtM(1:2*nenode,1:2*nenode) &
                        + ylength/2.d0 / epara(ielmt,1) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 21) then

                    elmtC(1:2*nenode,1:2*nenode) = elmtC(1:2*nenode,1:2*nenode) &
                        + xlength/2.d0 * epara(ielmt,2) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 22) then

                    elmtM(1:2*nenode,1:2*nenode) = elmtM(1:2*nenode,1:2*nenode) &
                        + xlength/2.d0 / epara(ielmt,2) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 31) then

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 32) then

                    elmtC(1:2*nenode,1:2*nenode) = elmtC(1:2*nenode,1:2*nenode) &
                        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 33) then

                    elmtC(1:2*nenode,1:2*nenode) = elmtC(1:2*nenode,1:2*nenode) &
                        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(N(1:2,1:2*nenode)), rho(1:2,1:2)), N(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 34) then

                    elmtM(1:2*nenode,1:2*nenode) = elmtM(1:2*nenode,1:2*nenode) &
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
                if (matrl(ielmt) /= 0) then
                    D=0d0
                    D(1,1)=(CS(matrl(ielmt))**2)*(2*rhos*(nuy - 1))/(2*nuy - 1)
                    D(1,2)= (CS(matrl(ielmt))**2)*(-(2*nuy*rhos)/(2*nuy - 1))
                    D(2,1)=D(1,2)
                    D(2,2)=(CS(matrl(ielmt))**2)*(2*rhos*(nuy - 1))/(2*nuy - 1)
                    D(3,3)=(CS(matrl(ielmt))**2)*rhos
                    D=D*dcmplx(1.d0,2.d0*beta)
                else
                    Drr(:,:) = 0.d0
                    Drr(1,1) = (CS_equi(ielmt)**2d0)*(2*rhos*(nuy - 1))/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    Drr(2,2) = (CS_equi(ielmt)**2d0)*rhos * dcmplx(1.d0,2.d0*beta)

                    Dss(:,:) = 0.d0
                    Dss(1,1) =(CS_equi(ielmt)**2d0)*rhos * dcmplx(1.d0,2.d0*beta)
                    Dss(2,2) = (CS_equi(ielmt)**2d0)*(2*rhos*(nuy - 1))/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)

                    Drs(:,:) = 0.d0
                    Drs(1,2) = -(2d0*CS_equi(ielmt)**2d0*nuy*rhos)/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)
                    Drs(2,1) = -(2d0*CS_equi(ielmt)**2d0*nuy*rhos)/(2*nuy - 1) * dcmplx(1.d0,2.d0*beta)

                    Dsr(:,:) = 0.d0
                    Dsr(1,2) = (CS_equi(ielmt)**2d0)*rhos * dcmplx(1.d0,2.d0*beta)
                    Dsr(2,1) = (CS_equi(ielmt)**2d0)*rhos* dcmplx(1.d0,2.d0*beta)
                endif
                if (etype(ielmt) == 0) then

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + matmul(matmul(transpose(B(1:3,1:2*nenode)), D(1:3,1:3)), B(1:3,1:2*nenode)) * detJ * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 11) then

                    elmtC(1:2*nenode,1:2*nenode) = elmtC(1:2*nenode,1:2*nenode) &
                        + ylength/2.d0 / epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtR(1:2*nenode,1:2*nenode) = elmtR(1:2*nenode,1:2*nenode) &
                        + epara(ielmt,1) / (ylength/2.d0) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 12) then

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + ylength/2.d0 * epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + 1.d0 / (ylength/2.d0) / epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 21) then

                    elmtC(1:2*nenode,1:2*nenode) = elmtC(1:2*nenode,1:2*nenode) &
                        + xlength/2.d0 / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtR(1:2*nenode,1:2*nenode) = elmtR(1:2*nenode,1:2*nenode) &
                        + epara(ielmt,2) / (xlength/2.d0) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 22) then

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + 1.d0 / (xlength/2.d0) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + xlength/2.d0 * epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 31) then

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 32) then

                    elmtC(1:2*nenode,1:2*nenode) = elmtC(1:2*nenode,1:2*nenode) &
                        + 1.d0 / epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtR(1:2*nenode,1:2*nenode) = elmtR(1:2*nenode,1:2*nenode) &
                        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 33) then

                    elmtC(1:2*nenode,1:2*nenode) = elmtC(1:2*nenode,1:2*nenode) &
                        + 1.d0 / epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                    elmtR(1:2*nenode,1:2*nenode) = elmtR(1:2*nenode,1:2*nenode) &
                        + epara(ielmt,1) * epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                else if (etype(ielmt) == 34) then

                    elmtK(1:2*nenode,1:2*nenode) = elmtK(1:2*nenode,1:2*nenode) &
                        + epara(ielmt,1) / epara(ielmt,2) * matmul(matmul(transpose(Brr(1:2,1:2*nenode)), Drr(1:2,1:2)), Brr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + epara(ielmt,2) / epara(ielmt,1) * matmul(matmul(transpose(Bss(1:2,1:2*nenode)), Dss(1:2,1:2)), Bss(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Brs(1:2,1:2*nenode)), Drs(1:2,1:2)), Brs(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2) &
                        + matmul(matmul(transpose(Bsr(1:2,1:2*nenode)), Dsr(1:2,1:2)), Bsr(1:2,1:2*nenode)) * wt1(ixi1) * wt2(ixi2)

                endif

            enddo

        enddo

        !elmtM_all(:,:,ielmt)=elmtM(:,:)
        !elmtK_all(:,:,ielmt)=elmtK(:,:)
        !elmtC_all(:,:,ielmt)=elmtC(:,:)
        !elmtR_all(:,:,ielmt)=elmtR(:,:)

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
                            !MSN22(itmp-nbc(1),jtmp-nbc(1))  = MSN22(itmp-nbc(1),jtmp-nbc(1))  + elmtM(2*(inode-1)+idof,2*(jnode-1)+jdof)
                            !DPN22(itmp-nbc(1),jtmp-nbc(1))  = DPN22(itmp-nbc(1),jtmp-nbc(1))  + elmtC(2*(inode-1)+idof,2*(jnode-1)+jdof)
                            !STN22(itmp-nbc(1),jtmp-nbc(1))  = STN22(itmp-nbc(1),jtmp-nbc(1))  + elmtK(2*(inode-1)+idof,2*(jnode-1)+jdof)
                            !RRN22(itmp-nbc(1),jtmp-nbc(1))  = RRN22(itmp-nbc(1),jtmp-nbc(1))  + elmtR(2*(inode-1)+idof,2*(jnode-1)+jdof)


                        else if ((itmp <= nbc(2)) .and. (jtmp > nbc(2))) then

                            !MSN23(itmp,jtmp-nbc(2)) = MSN23(itmp,jtmp-nbc(2)) + elmtM(2*(inode-1)+idof,2*(jnode-1)+jdof)
                            !DPN23(itmp,jtmp-nbc(2)) = DPN23(itmp,jtmp-nbc(2)) + elmtC(2*(inode-1)+idof,2*(jnode-1)+jdof)
                            !STN23(itmp,jtmp-nbc(2)) = STN23(itmp,jtmp-nbc(2)) + elmtK(2*(inode-1)+idof,2*(jnode-1)+jdof)
                            !RRN23(itmp,jtmp-nbc(2)) = RRN23(itmp,jtmp-nbc(2)) + elmtR(2*(inode-1)+idof,2*(jnode-1)+jdof)

                        else if ((itmp > nbc(2)) .and. (jtmp > nbc(2))) then

                            !MSN33(itmp-nbc(2),jtmp-nbc(2)) = MSN33(itmp-nbc(2),jtmp-nbc(2)) + elmtM(2*(inode-1)+idof,2*(jnode-1)+jdof)
                            !DPN33(itmp-nbc(2),jtmp-nbc(2)) = DPN33(itmp-nbc(2),jtmp-nbc(2)) + elmtC(2*(inode-1)+idof,2*(jnode-1)+jdof)
                            !STN33(itmp-nbc(2),jtmp-nbc(2)) = STN33(itmp-nbc(2),jtmp-nbc(2)) + elmtK(2*(inode-1)+idof,2*(jnode-1)+jdof)
                            !RRN33(itmp-nbc(2),jtmp-nbc(2)) = RRN33(itmp-nbc(2),jtmp-nbc(2)) + elmtR(2*(inode-1)+idof,2*(jnode-1)+jdof)

                        endif
                    enddo
                enddo
            enddo
        enddo

    enddo

    end subroutine FiniteElement_PlaneStrain


    subroutine CalcJ(nenode, J, xi, eta, x, y)
    implicit none
    integer nenode
    double precision J(2,2), xi, eta, x(nenode), y(nenode)

    integer itmp
    double precision shp(nenode), dshp(nenode,2)

    call shpfnc(nenode, shp(:), dshp(:,:), xi, eta)

    J(:,:) = 0.d0

    do itmp = 1 , nenode

        J(1,1) = J(1,1) + dshp(itmp,1) * x(itmp)
        J(1,2) = J(1,2) + dshp(itmp,1) * y(itmp)
        J(2,1) = J(2,1) + dshp(itmp,2) * x(itmp)
        J(2,2) = J(2,2) + dshp(itmp,2) * y(itmp)

    enddo

    end subroutine CalcJ

    subroutine shpfnc(nenode, shp, dshp, xi, eta)
    implicit none
    integer nenode
    double precision shp(nenode), dshp(nenode,2), xi, eta

    if (nenode == 3) then

        shp(1) = 1.d0 - xi - eta
        shp(2) = xi
        shp(3) = eta

        dshp(1,1) = -1.d0
        dshp(2,1) =  1.d0
        dshp(3,1) =  0.d0

        dshp(1,2) = -1.d0
        dshp(2,2) =  0.d0
        dshp(3,2) =  1.d0

    else if (nenode == 6) then

        shp(4) = 4.d0 *  xi * (1.d0 - xi - eta)
        shp(5) = 4.d0 *  xi * eta
        shp(6) = 4.d0 * eta * (1.d0 - xi - eta)
        shp(1) = 1.d0 - xi - eta - 0.5d0 * shp(4) - 0.5d0 * shp(6)
        shp(2) = xi              - 0.5d0 * shp(4) - 0.5d0 * shp(5)
        shp(3) = eta             - 0.5d0 * shp(5) - 0.5d0 * shp(6)

        dshp(4,1) =  4.d0 * (1.d0 - 2.d0 * xi - eta)
        dshp(5,1) =  4.d0 * eta
        dshp(6,1) = -4.d0 * eta
        dshp(1,1) = -1.d0 - 0.5d0 * dshp(4,1) - 0.5d0 * dshp(6,1)
        dshp(2,1) =  1.d0 - 0.5d0 * dshp(4,1) - 0.5d0 * dshp(5,1)
        dshp(3,1) =       - 0.5d0 * dshp(5,1) - 0.5d0 * dshp(6,1)

        dshp(4,2) = -4.d0 * xi
        dshp(5,2) =  4.d0 * xi
        dshp(6,2) =  4.d0 * (1.d0 - xi - 2.d0 * eta)
        dshp(1,2) = -1.d0 - 0.5d0 * dshp(4,2) - 0.5d0 * dshp(6,2)
        dshp(2,2) =       - 0.5d0 * dshp(4,2) - 0.5d0 * dshp(5,2)
        dshp(3,2) =  1.d0 - 0.5d0 * dshp(5,2) - 0.5d0 * dshp(6,2)

    else if (nenode == 4) then

        shp(1) = (1.d0-xi) * (1.d0-eta) / 4.d0
        shp(2) = (1.d0+xi) * (1.d0-eta) / 4.d0
        shp(3) = (1.d0+xi) * (1.d0+eta) / 4.d0
        shp(4) = (1.d0-xi) * (1.d0+eta) / 4.d0

        dshp(1,1) = -(1.d0-eta) / 4.d0
        dshp(2,1) =  (1.d0-eta) / 4.d0
        dshp(3,1) =  (1.d0+eta) / 4.d0
        dshp(4,1) = -(1.d0+eta) / 4.d0

        dshp(1,2) = -(1.d0-xi) / 4.d0
        dshp(2,2) = -(1.d0+xi) / 4.d0
        dshp(3,2) =  (1.d0+xi) / 4.d0
        dshp(4,2) =  (1.d0-xi) / 4.d0

    else if (nenode == 8) then

        shp(5) = (1.d0-xi**2) * (1.d0-eta   ) / 2.d0
        shp(6) = (1.d0+xi   ) * (1.d0-eta**2) / 2.d0
        shp(7) = (1.d0-xi**2) * (1.d0+eta   ) / 2.d0
        shp(8) = (1.d0-xi   ) * (1.d0-eta**2) / 2.d0
        shp(1) = (1.d0-xi   ) * (1.d0-eta   ) / 4.d0 - (shp(8) + shp(5)) / 2.d0
        shp(2) = (1.d0+xi   ) * (1.d0-eta   ) / 4.d0 - (shp(5) + shp(6)) / 2.d0
        shp(3) = (1.d0+xi   ) * (1.d0+eta   ) / 4.d0 - (shp(6) + shp(7)) / 2.d0
        shp(4) = (1.d0-xi   ) * (1.d0+eta   ) / 4.d0 - (shp(7) + shp(8)) / 2.d0

        dshp(5,1) = -2.d0*xi * (1.d0-eta) / 2.d0
        dshp(6,1) =  (1.d0-eta**2) / 2.d0
        dshp(7,1) = -2.d0*xi * (1.d0+eta) / 2.d0
        dshp(8,1) = -(1.d0-eta**2) / 2.d0
        dshp(1,1) = -(1.d0-eta) / 4.d0 - (dshp(8,1) + dshp(5,1)) / 2.d0
        dshp(2,1) =  (1.d0-eta) / 4.d0 - (dshp(5,1) + dshp(6,1)) / 2.d0
        dshp(3,1) =  (1.d0+eta) / 4.d0 - (dshp(6,1) + dshp(7,1)) / 2.d0
        dshp(4,1) = -(1.d0+eta) / 4.d0 - (dshp(7,1) + dshp(8,1)) / 2.d0

        dshp(5,2) = -(1.d0-xi**2) / 2.d0
        dshp(6,2) = -2.d0*eta * (1.d0+xi) / 2.d0
        dshp(7,2) =  (1.d0-xi**2) / 2.d0
        dshp(8,2) = -2.d0*eta * (1.d0-xi) / 2.d0
        dshp(1,2) = -(1.d0-xi) / 4.d0 - (dshp(8,2) + dshp(5,2)) / 2.d0
        dshp(2,2) = -(1.d0+xi) / 4.d0 - (dshp(5,2) + dshp(6,2)) / 2.d0
        dshp(3,2) =  (1.d0+xi) / 4.d0 - (dshp(6,2) + dshp(7,2)) / 2.d0
        dshp(4,2) =  (1.d0-xi) / 4.d0 - (dshp(7,2) + dshp(8,2)) / 2.d0

    else if (nenode == 9) then

        shp(9) = (1.d0-xi**2) * (1.d0-eta**2)
        shp(5) = (1.d0-xi**2) * (1.d0-eta   ) / 2.d0 - shp(9) / 2.d0
        shp(6) = (1.d0+xi   ) * (1.d0-eta**2) / 2.d0 - shp(9) / 2.d0
        shp(7) = (1.d0-xi**2) * (1.d0+eta   ) / 2.d0 - shp(9) / 2.d0
        shp(8) = (1.d0-xi   ) * (1.d0-eta**2) / 2.d0 - shp(9) / 2.d0
        shp(1) = (1.d0-xi   ) * (1.d0-eta   ) / 4.d0 - (shp(8) + shp(5)) / 2.d0 - shp(9) / 4.d0
        shp(2) = (1.d0+xi   ) * (1.d0-eta   ) / 4.d0 - (shp(5) + shp(6)) / 2.d0 - shp(9) / 4.d0
        shp(3) = (1.d0+xi   ) * (1.d0+eta   ) / 4.d0 - (shp(6) + shp(7)) / 2.d0 - shp(9) / 4.d0
        shp(4) = (1.d0-xi   ) * (1.d0+eta   ) / 4.d0 - (shp(7) + shp(8)) / 2.d0 - shp(9) / 4.d0

        dshp(9,1) = -2.d0*xi * (1.d0-eta**2)
        dshp(5,1) = -2.d0*xi * (1.d0-eta) / 2.d0 - dshp(9,1) / 2.d0
        dshp(6,1) =  (1.d0-eta**2) / 2.d0 - dshp(9,1) / 2.d0
        dshp(7,1) = -2.d0*xi * (1.d0+eta) / 2.d0 - dshp(9,1) / 2.d0
        dshp(8,1) = -(1.d0-eta**2) / 2.d0 - dshp(9,1) / 2.d0
        dshp(1,1) = -(1.d0-eta) / 4.d0 - (dshp(8,1) + dshp(5,1)) / 2.d0 - dshp(9,1) / 4.d0
        dshp(2,1) =  (1.d0-eta) / 4.d0 - (dshp(5,1) + dshp(6,1)) / 2.d0 - dshp(9,1) / 4.d0
        dshp(3,1) =  (1.d0+eta) / 4.d0 - (dshp(6,1) + dshp(7,1)) / 2.d0 - dshp(9,1) / 4.d0
        dshp(4,1) = -(1.d0+eta) / 4.d0 - (dshp(7,1) + dshp(8,1)) / 2.d0 - dshp(9,1) / 4.d0

        dshp(9,2) = -2.d0*eta * (1.d0-xi**2)
        dshp(5,2) = -(1.d0-xi**2) / 2.d0 - dshp(9,2) / 2.d0
        dshp(6,2) = -2.d0*eta * (1.d0+xi) / 2.d0 - dshp(9,2) / 2.d0
        dshp(7,2) =  (1.d0-xi**2) / 2.d0 - dshp(9,2) / 2.d0
        dshp(8,2) = -2.d0*eta * (1.d0-xi) / 2.d0 - dshp(9,2) / 2.d0
        dshp(1,2) = -(1.d0-xi) / 4.d0 - (dshp(8,2) + dshp(5,2)) / 2.d0 - dshp(9,2) / 4.d0
        dshp(2,2) = -(1.d0+xi) / 4.d0 - (dshp(5,2) + dshp(6,2)) / 2.d0 - dshp(9,2) / 4.d0
        dshp(3,2) =  (1.d0+xi) / 4.d0 - (dshp(6,2) + dshp(7,2)) / 2.d0 - dshp(9,2) / 4.d0
        dshp(4,2) =  (1.d0-xi) / 4.d0 - (dshp(7,2) + dshp(8,2)) / 2.d0 - dshp(9,2) / 4.d0

    endif

    end subroutine shpfnc


    !----------------------------------------------------------------------c
    subroutine GaussPt(ngpr,gp,wt)
    !----------------------------------------------------------------------c
    !
    !   gauss points and weights for 9-node rectangular element
    !
    implicit none
    integer ngpr

    double precision gp(ngpr),wt(ngpr)

    !
    !   initialize
    !
    gp=0.d0
    wt=0.d0

    !
    !   ngpr routine decision
    !
    if(ngpr.eq. 1) goto  100
    if(ngpr.eq. 2) goto  200
    if(ngpr.eq. 3) goto  300
    if(ngpr.eq.10) goto 1000

    !
    !   if number of gausspoint(ngpr) = 1
    !
100 gp(1) = 0.d0
    wt(1) = 2.d0
    return

    !
    !   if number of gausspoint(ngpr) = 2
    !
200 gp(1) =-1.d0/dsqrt(3.d0)
    gp(2) =-gp(1)
    wt(1) = 1.d0
    wt(2) = 1.d0
    return

    !
    !   if number of gausspoint(ngpr) = 3
    !
300 gp(1) =-dsqrt(0.6d0)
    gp(2) = 0.d0
    gp(3) =-gp(1)
    wt(1) = 5.d0/9.d0
    wt(2) = 8.d0/9.d0
    wt(3) = wt(1)
    return

1000 gp(1) =-0.9739065285d0
    gp(2) =-0.8650633666d0
    gp(3) =-0.6794095683d0
    gp(4) =-0.4333953941d0
    gp(5) =-0.1488743389d0
    gp(6) =-gp(5)
    gp(7) =-gp(4)
    gp(8) =-gp(3)
    gp(9) =-gp(2)
    gp(10)=-gp(1)
    wt(1) = 0.0666713443d0
    wt(2) = 0.1494513491d0
    wt(3) = 0.2190863625d0
    wt(4) = 0.2692667193d0
    wt(5) = 0.2955242247d0
    wt(6) = wt(5)
    wt(7) = wt(4)
    wt(8) = wt(3)
    wt(9) = wt(2)
    wt(10)= wt(1)
    return

    end	subroutine GaussPt
