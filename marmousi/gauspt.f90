!----------------------------------------------------------------------c
subroutine gausrect(ngpr,gp1,gp2,wt)
!----------------------------------------------------------------------c
!
!   gauss points and weights for 9-node rectangular element
!
	implicit none
	integer ngpr
	
	double precision gp1(ngpr),gp2(ngpr),wt(ngpr)

!
!   initialize
!
	gp1=0.d0
	gp2=0.d0
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
100	gp1(1) = 0.d0
	gp2(:) = gp1(:)
	wt(1) = 2.d0
	return

!
!   if number of gausspoint(ngpr) = 2
!
200	gp1(1) =-1.d0/dsqrt(3.d0)
	gp1(2) =-gp1(1)
	gp2(:) = gp1(:)
	wt(1) = 1.d0
	wt(2) = 1.d0
	return

!
!   if number of gausspoint(ngpr) = 3
!
300	gp1(1) =-dsqrt(0.6d0)
	gp1(2) = 0.d0
	gp1(3) =-gp1(1)
	gp2(:) = gp1(:)
	wt(1) = 5.d0/9.d0
	wt(2) = 8.d0/9.d0
	wt(3) = wt(1)
	return

1000	gp1(1) =-0.9739065285d0
		gp1(2) =-0.8650633666d0
		gp1(3) =-0.6794095683d0
		gp1(4) =-0.4333953941d0
		gp1(5) =-0.1488743389d0
		gp1(6) =-gp1(5)
		gp1(7) =-gp1(4)
		gp1(8) =-gp1(3)
		gp1(9) =-gp1(2)
		gp1(10)=-gp1(1)
		gp2(:) = gp1(:)
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

end	subroutine gausrect

!----------------------------------------------------------------------c
subroutine gaustri(ngpt,gp1,gp2,wt)
!----------------------------------------------------------------------c
!
!   gauss points and weights for 6-node triangular element
!
	implicit none
	integer ngpt

	double precision gp1(ngpt),gp2(ngpt),wt(ngpt)

!
!   initialize
!
	gp1(:)=0.d0
	gp2(:)=0.d0
	wt (:)=0.d0

!
!   ngpt routine decision
!
	if(ngpt.eq. 7) goto  700
!	if(ngpt.eq.12) goto 1200
!	if(ngpt.eq.13) goto 1300
!
!   if number of gausspoint(ngpt) = 7
!
700	gp1(1) = 0.101286507323456d0
	gp1(2) = 0.797426985353087d0
	gp1(3) = gp1(1)
	gp1(4) = 0.470142064105115d0
	gp1(5) = gp1(4)
	gp1(6) = 0.059715871789770d0
	gp1(7) = 1.d0/3.d0

	gp2(1) = gp1(1)
	gp2(2) = gp1(1)
	gp2(3) = gp1(2)
	gp2(4) = gp1(6)
	gp2(5) = gp1(4)
	gp2(6) = gp1(4)
	gp2(7) = gp1(7)

	wt(1)  = 0.125939180544827d0
	wt(2)  = wt(1)
	wt(3)  = wt(1)
	wt(4)  = 0.132394152788506d0
	wt(5)  = wt(4)
	wt(6)  = wt(4)
	wt(7)  = 0.225d0

end subroutine gaustri