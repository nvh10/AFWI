subroutine scaling(nn, AA, BB)
	implicit none
	integer nn
	double complex AA(nn,nn), BB(nn,nn)

	integer itmp, jtmp
	double precision mAA

	!	Scaling
	do itmp = 1 , nn

		mAA = 0.d0

		do jtmp = 1 , nn

			if (cdabs(AA(jtmp,itmp)) > mAA) mAA = cdabs(AA(jtmp,itmp))

		enddo

		AA(1:nn,itmp) = AA(1:nn,itmp) / mAA
		BB(1:nn,itmp) = BB(1:nn,itmp) / mAA

	enddo

end subroutine scaling

subroutine condense(n1, AA1, AA2, n2, BB1, BB2, cdof)
	implicit none
	integer n1, n2, cdof(n1)
	double complex AA1(n1,n1), AA2(n1), BB1(n2,n2), BB2(n2)

	integer nC11, nC22, itmp, jtmp, ktmp, ltmp
	double complex C11(n1,n1), C22(n1,n1), C12(n1,n1), D1(n1), D2(n1), invC22(n1,n1)
		
	C11(:,:) = 0.d0
	C22(:,:) = 0.d0
	C12(:,:) = 0.d0

	ktmp = 1

	do itmp = 1 , n1

		if (cdof(itmp) == 1) then

			ltmp = 1

			do jtmp = 1 , n1

				if (cdof(jtmp) == 1) then
			
					C11(ktmp,ltmp) = AA1(itmp,jtmp)
					D1 (ktmp)      = AA2(itmp)

					ltmp = ltmp + 1

				endif

			enddo

			ktmp = ktmp + 1

		endif

	enddo

	nC11 = ktmp - 1
	
	ktmp = 1

	do itmp = 1 , n1

		if (cdof(itmp) == -1) then

			ltmp = 1

			do jtmp = 1 , n1

				if (cdof(jtmp) == -1) then
			
					C22(ktmp,ltmp) = AA1(itmp,jtmp)
					D2 (ktmp)      = AA2(itmp)

					ltmp = ltmp + 1

				endif

			enddo

			ktmp = ktmp + 1

		endif

	enddo

	nC22 = ktmp - 1
	
	ktmp = 1

	do itmp = 1 , n1

		if (cdof(itmp) == 1) then

			ltmp = 1

			do jtmp = 1 , n1

				if (cdof(jtmp) == -1) then
			
					C12(ktmp,ltmp) = AA1(itmp,jtmp)

					ltmp = ltmp + 1

				endif

			enddo

			ktmp = ktmp + 1

		endif

	enddo

	call InverseMatrixZ(nC22, C22(1:nC22,1:nC22), invC22(1:nC22,1:nC22))

	BB1(:,:) = C11(1:nC11,1:nC11) - matmul(C12(1:nC11,1:nC22),matmul(invC22(1:nC22,1:nC22),transpose(C12(1:nC11,1:nC22))))
	BB2(:)   = D1 (1:nC11)        - matmul(C12(1:nC11,1:nC22),matmul(invC22(1:nC22,1:nC22),          D2 (1:nC22)))

end subroutine condense

subroutine SubstituteMatrix(n1, AA1, n2, BB1, cdof)
	implicit none
	integer n1, n2, cdof(n1)
	double complex AA1(n1,n1), BB1(n2,n2)

	integer itmp, jtmp

	do itmp = 1 , n1

		do jtmp = 1 , n1

			if ((cdof(itmp) /= 0) .and. (cdof(jtmp) /= 0)) &
				BB1(cdof(itmp),cdof(jtmp)) = BB1(cdof(itmp),cdof(jtmp)) + AA1(itmp,jtmp)

		enddo

	enddo

end subroutine SubstituteMatrix

subroutine SubstituteMatrix2(n1, AA1, AA2, n2, BB1, BB2, cdof)
	implicit none
	integer n1, n2, cdof(n1)
	double complex AA1(n1,n1), AA2(n1), BB1(n2,n2), BB2(n2)

	integer itmp, jtmp

	do itmp = 1 , n1

		if (cdof(itmp) /= 0) BB2(cdof(itmp)) = BB2(cdof(itmp)) + AA2(itmp)

		do jtmp = 1 , n1

			if ((cdof(itmp) /= 0) .and. (cdof(jtmp) /= 0)) &
				BB1(cdof(itmp),cdof(jtmp)) = BB1(cdof(itmp),cdof(jtmp)) + AA1(itmp,jtmp)

		enddo

	enddo

end subroutine SubstituteMatrix2
