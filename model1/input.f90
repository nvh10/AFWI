subroutine inputHyd(inhyd, nx, nz, dx, dz, nnode, nelmt, x, z, enode, matrl, nbc, bc, id)
	implicit none
	! input & output
	integer inhyd

	! cotrol data
	integer nnode, nelmt

	! node data
	double precision x(nnode), z(nnode)

	! element data
	integer enode(nelmt,0:9), matrl(nelmt)

	! boundary condition data
	integer nbc(2), bc(nnode), id(nnode)
	
	integer inode, ielmt, imatrl, ibc, idof, itmp, jtmp, helmt, nx, nz
	double precision dx, dz

	! read from HX file
!	read (inhyd,*)
	
!	do inode = 1 , nnode

!		read (inhyd,*) itmp, x(itmp), z(itmp)

!	enddo
	
	do jtmp = 1 , nz + 1
		
		do itmp	= 1 , nx + 1
			
			inode = (nx + 1) * (jtmp - 1) + itmp
			x(inode) = (itmp - 1) * dx
			z(inode) = (jtmp - 1) * dz			
			
		enddo
		
	enddo
	
!	read (inhyd,*)

!	do ielmt = 1 , nelmt

!		read (inhyd,*) itmp, enode(itmp,0), (enode(itmp,jtmp),jtmp=1,enode(itmp,0))

!	enddo
	
	do jtmp = 1 , nz
		
		do itmp = 1 , nx
			
			ielmt = nx * (jtmp - 1) + itmp
			enode(ielmt,0) = 4
			enode(ielmt,1) = (nx + 1) * (jtmp - 1) + itmp
			enode(ielmt,2) = (nx + 1) * (jtmp - 1) + itmp + 1
			enode(ielmt,3) = (nx + 1) *  jtmp      + itmp + 1
			enode(ielmt,4) = (nx + 1) *  jtmp      + itmp
			
			matrl(ielmt) = 1
			if (jtmp >= 31) matrl(ielmt) = 2
			
		enddo
		
	enddo
	
	bc = 1

!	read (inhyd,*)
!	read (inhyd,*) itmp

!	do ibc = 1 , itmp

!		read (inhyd,*) jtmp, bc(jtmp)

!	enddo

	nbc(1) = 0		! free DOFs that will be condensed out

	do inode = 1 , nnode

		if (bc(inode) == 1) then			! ground nodes on rigid foundation

			nbc(1) = nbc(1) + 1
				
			id(inode) = nbc(1)

		endif

	enddo

	nbc(2) = 0		! free DOFs that will be condensed out

	do inode = 1 , nnode

		if (bc(inode) == 2) then			! ground nodes on rigid foundation

			nbc(2) = nbc(2) + 1
				
			id(inode) = nbc(1) + nbc(2)

		endif

	enddo

end subroutine inputHyd