subroutine TB0(nnode, nelmt, ww, Cw, rhow, alphr, DSTB, PETB, z, Vfr)
	implicit none
	integer nnode, nelmt
	double precision pi, ww, rhow, alphr, r1, z(nnode)
	double complex Cw
	double complex DSTB(nnode,nnode), PETB(nnode), invDS(nnode,nnode), Phfr(nnode), Vfr(nnode)
	
	double precision elmtG(2,2), elmtH(2,2)
	double complex cReservoir, GF(nnode,nnode), HF(nnode,nnode), CF(nnode,nnode), DF(nnode-1,nnode-1), BF(nnode-1,nnode)
	
	double complex ALPHA(nnode-1), BETA(nnode-1), Psi(nnode-1,nnode-1), la(nnode-1,nnode-1)
	double complex k, H0, H1, dH, c
	
	integer ielmt, itmp, ieig, inode
	double precision edep
	
	pi = 4.d0 * datan(1.d0)
	
	cReservoir = 1.d0/Cw * (1.d0 - alphr) / (1.d0 + alphr);

	GF = 0.d0
	HF = 0.d0
	CF = 0.d0

	do ielmt = 1 , nelmt
    
		edep = z(ielmt+1) - z(ielmt)
    
	    elmtG = 0.d0
	    elmtG(1,1) = edep / 3.d0
		elmtG(1,2) = edep / 6.d0
	    elmtG(2,1) = edep / 6.d0
		elmtG(2,2) = edep / 3.d0
            
	    elmtH = 0.d0
		elmtH(1,1) =  1.d0 / edep
	    elmtH(1,2) = -1.d0 / edep
		elmtH(2,1) = -1.d0 / edep
	    elmtH(2,2) =  1.d0 / edep
    
		GF(ielmt:ielmt+1,ielmt:ielmt+1) = GF(ielmt:ielmt+1,ielmt:ielmt+1) + 1/Cw**2 * elmtG
	    HF(ielmt:ielmt+1,ielmt:ielmt+1) = HF(ielmt:ielmt+1,ielmt:ielmt+1) + elmtH
    
	enddo

	CF(1,1) = cReservoir
	
    DSTB = -ww**2 * GF + dcmplx(0.d0,ww) * CF + HF
    call InverseMatrixZ(nnode-1, DSTB(1:nnode-1,1:nnode-1), invDS(1:nnode-1,1:nnode-1))
    
    PETB = 0.d0
    PETB(1:nnode-1) = rhow * invDS(1:nnode-1,1)
    
    Phfr = 0.d0
    Phfr = 1.d0/(-rhow) * PETB
    
    Vfr(1) = (Phfr(2) - Phfr(1)) / (z(2) - z(1))
    
    do inode = 2 , nnode - 1
    
        Vfr(inode) = 0.5d0 * (Phfr(inode+1) - Phfr(inode)) / (z(inode+1) - z(inode)) + 0.5d0 * (Phfr(inode) - Phfr(inode-1)) / (z(inode) - z(inode-1))
    
    enddo
    
    Vfr(nnode) = (Phfr(nnode) - Phfr(nnode-1)) / (z(nnode) - z(nnode-1))
    
    Psi = 0.d0
    la = 0.d0
    
    call eigenZ(nnode-1, dcmplx(0.d0,ww)*CF(1:nnode-1,1:nnode-1)+HF(1:nnode-1,1:nnode-1), Cw**2*GF(1:nnode-1,1:nnode-1), &
		ALPHA, BETA, Psi)
		
	la = matmul(matmul(transpose(Psi),Cw**2*GF(1:nnode-1,1:nnode-1)),Psi)
	
	do itmp = 1 , nnode - 1

		Psi(:,itmp) = Psi(:,itmp) / cdsqrt(la(itmp,itmp))
		
		la(itmp,itmp) = cdsqrt(ALPHA(itmp) / BETA(itmp))
		
	enddo
	
    DF = 0.d0

	do ieig = 1 , nnode-1

        k = cdsqrt(ww**2/Cw**2 - la(ieig,ieig)**2)
        if (dimag(k) > 0.d0) k = -k

        DF(ieig,ieig) = dcmplx(0.d0,1.d0) * k
        
	enddo
    
    BF = 0.d0
    BF = Cw**2 * matmul(transpose(Psi), GF(1:nnode-1,:))
    
    DSTB = matmul(matmul(transpose(BF), DF), BF)
    
    PETB = matmul(DSTB, PETB)
	
end subroutine TB0