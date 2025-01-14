subroutine Fourier1(nnode, nelmt, Cw, rhow, alphr, r, z, enode, ngelmt, gnode, &
	nbc, bc, id, cnode, r0, z0, nfreq, dfreq, outsgb, outrhb, outout)
	implicit none
	integer nnode(2), nelmt(2), ngelmt
	double precision rhow, alphr
	double complex Cw

	! node data
	double precision r(2,nnode(1)), z(2,nnode(1))

	! element data
	integer enode(nelmt(1),0:9), gnode(ngelmt,0:3)

	! boundary condition data
	integer nbc(3), bc(nnode(1)), id(nnode(1)), cnode(nnode(2))
	double precision r0, z0

	! frequency data
	integer nfreq
	double precision dfreq

	integer outsgb, outrhb, outout

	double precision pi, ww
	integer cdof(nnode(1))

	double complex DSN11(nbc(1),nbc(1)+1), DSN21(nbc(2),nbc(1)+1), &
		DSN13(nbc(1),nbc(3)), DSN23(nbc(2),nbc(3)), DSNtmp(nbc(2),nbc(1)+1), invDSN11(nbc(1),nbc(1)), P1(nbc(1))
	double complex,allocatable :: DSN22(:)
	
	double complex DSTB(nnode(2),nnode(2)), PETB(nnode(2))

	integer nn, n1, n2, nfix, ne, nm, ifreq, itmp, jtmp, inode, jnode, idof, jdof, ielmt, iCapp 
	double precision rtmi, rtmj, ztmi, ztmj

	integer,allocatable :: ia(:), ja(:), idcheck(:)
	integer equation, node, position, iparm(64), idum, ierr, iia
	integer*8 pt(64)

	double complex argJ, J0, J1, Jc

	pi = 4.d0 * datan(1.d0)

	allocate (ia(nbc(2)+1), idcheck(nbc(2)))

	ia(1) = 1

	do equation = 1 , nbc(2)

		node = 0

		do inode = 1 , nnode(1)

			if (id(inode) == equation+nbc(1)) then

				node = inode

			endif

		enddo

		idcheck(:) = 0

		do ielmt = 1 , nelmt(1)
		do inode = 1 , enode(ielmt,0)

			if (enode(ielmt,inode) == node) then

				do jnode = 1 , enode(ielmt,0)

					if ((id(enode(ielmt,jnode)) >= nbc(1)+1) .and. &
						(id(enode(ielmt,jnode)) <= nbc(1)+nbc(2))) then

						idcheck(id(enode(ielmt,jnode))-nbc(1)) = 1

					endif

				enddo				

			endif

		enddo
		enddo

		do inode = 1 , nnode(2)

			if (cnode(inode) == node) then

				do jnode = 1 , nnode(2)

					if ((id(cnode(jnode)) >= nbc(1)+1) .and. &
						(id(cnode(jnode)) <= nbc(1)+nbc(2))) then

						idcheck(id(cnode(jnode))-nbc(1)) = 1

					endif

				enddo

			endif

		enddo

		ia(equation+1) = ia(equation)

		do idof = 1 , nbc(2)

			if (idcheck(idof) == 1) then

				ia(equation+1) = ia(equation+1) + 1

			endif

		enddo

	enddo

	allocate (ja(ia(nbc(2)+1)-1), DSN22(ia(nbc(2)+1)-1))
		
	position = 0

	do equation = 1 , nbc(2)

		node = 0

		do inode = 1 , nnode(1)

			if (id(inode) == equation+nbc(1)) then

				node = inode

			endif

		enddo

		idcheck(:) = 0

		do ielmt = 1 , nelmt(1)
		do inode = 1 , enode(ielmt,0)

			if (enode(ielmt,inode) == node) then

				do jnode = 1 , enode(ielmt,0)

					if ((id(enode(ielmt,jnode)) >= nbc(1)+1) .and. &
						(id(enode(ielmt,jnode)) <= nbc(1)+nbc(2))) then

						idcheck(id(enode(ielmt,jnode))-nbc(1)) = 1

					endif

				enddo				

			endif

		enddo
		enddo

		do inode = 1 , nnode(2)

			if (cnode(inode) == node) then

				do jnode = 1 , nnode(2)

					if ((id(cnode(jnode)) >= nbc(1)+1) .and. &
						(id(cnode(jnode)) <= nbc(1)+nbc(2))) then

						idcheck(id(cnode(jnode))-nbc(1)) = 1

					endif

				enddo

			endif

		enddo

		do idof = 1 , nbc(2)

			if (idcheck(idof) == 1) then

				position = position + 1

				ja(position) = idof

			endif

		enddo

	enddo

	call pardisoinit(pt, 13, iparm)

	do ifreq = 1 , nfreq
			   
		write (*,'(1h+,a14,i5,a3,i5)') '   Freq step :',ifreq,' / ',nfreq

		ww = 2.d0 * pi * ifreq * dfreq

		! FEM for near field
		nn = nnode(1)
		ne = nelmt(1)

		DSN11 = 0.d0
		DSN21 = 0.d0
		DSN22 = 0.d0
		DSN13 = 0.d0
		DSN23 = 0.d0

		call FEM1(nn, ne, ww, Cw, rhow, alphr, DSN11, DSN21, DSN22, DSN13, DSN23, &
			r(1,1:nn), z(1,1:nn), enode, nbc, bc, id, ia, ja, ngelmt, gnode)
			
		nn = nnode(2)
		ne = nelmt(2)

		DSTB = 0.d0
		PETB = 0.d0

		call TB1(nn, ne, ww, Cw, rhow, alphr, DSTB, PETB, r(2,1), z(2,1:nn))

		cdof(:) = 0

		do itmp = 1 , nnode(2)

			cdof(itmp) = id(cnode(itmp))

		enddo

		do inode = 1 , nnode(2)
			
			itmp = cdof(inode)
			
			do jnode = 1 , nnode(2)
				
				jtmp = cdof(jnode)
				
				if ((itmp <= nbc(1)) .and. (jtmp <= nbc(1))) then

					DSN11(itmp,jtmp) = DSN11(itmp,jtmp) + DSTB(inode,jnode)

				else if ((itmp > nbc(1)) .and. (itmp <= nbc(1)+nbc(2)) .and. (jtmp <= nbc(1))) then

					DSN21(itmp-nbc(1),jtmp) = DSN21(itmp-nbc(1),jtmp) + DSTB(inode,jnode)

				else if ((itmp > nbc(1)) .and. (itmp <= nbc(1)+nbc(2)) .and. &
					  (jtmp > nbc(1)) .and. (jtmp <= nbc(1)+nbc(2))) then

					do iia = ia(itmp-nbc(1)) , ia(itmp+1-nbc(1))-1

						if (ja(iia) == jtmp-nbc(1)) then

							DSN22(iia) = DSN22(iia) + DSTB(inode,jnode)

						endif

					enddo

				endif

			enddo				
			
!			if (itmp <= nbc(1)) then

!				DSN11(itmp,nbc(1)+1) = DSN11(itmp,nbc(1)+1) + PETB(inode)

!			else if ((itmp > nbc(1)) .and. (itmp <= nbc(1)+nbc(2))) then

!				DSN21(itmp-nbc(1),nbc(1)+1) = DSN21(itmp-nbc(1),nbc(1)+1) + PETB(inode)
				
!			endif

		enddo

		iparm(3) = 8	! number of processors
		iparm(8) = 0

		if (ifreq == 1) then

			call pardiso(pt, 1, 1, 13, 13, nbc(2), DSN22, ia, ja, idum, nbc(1)+1, iparm, 1, DSN21, DSNtmp, ierr)
			
		else

			call pardiso(pt, 1, 1, 13, 23, nbc(2), DSN22, ia, ja, idum, nbc(1)+1, iparm, 0, DSN21, DSNtmp, ierr)

		endif

		DSN11 = DSN11 - matmul(transpose(DSN21(:,1:nbc(1))), DSNtmp)
		
		call InverseMatrixZ(nbc(1), DSN11(1:nbc(1),1:nbc(1)), invDSN11)
		
!		P1 = - matmul(matmul(invDSN11, transpose(DSN21(:,1:nbc(1)))), DSNtmp(:,nbc(1)+1))
		
		write (outsgb) invDSN11
!		write (outsgb) P1

	enddo
	
	call pardiso(pt, 1, 1, 13, -1, nbc(2), DSN22, ia, ja, idum, nbc(1)+1, iparm, 0, DSN21, DSNtmp, ierr)

	deallocate (ia, ja, idcheck, DSN22)

end subroutine Fourier1