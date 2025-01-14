      subroutine bessel (z,j0,j1,c)
c
c     computation of the bessel functions of zero and first order
c     z=complex argument=(x,y)   -180. lt arg(z) lt 180. degrees
c     c=j0/j1
c
c     written by e.kausel
c
      implicit double precision (a-h,o-z)
      complex*16      z,c,p,a,a0,a1,b0,b1,j0,j1

      x=dble(z)
      y=dimag(z)
      r=dsqrt(x*x+y*y)
      if (r.le.10.d0) go to 30
      phi=datan2(y,x)
      l=0
      if (x.ge.0.d0) go to 10
      l=1
      z=-z
      phi=phi-3.14159265358793d0*dsign(1.d0,phi)
      x=-x
      y=-y
   10 j=2.d0*r
      c=(0.d0,-0.125d0)/z
      k=2
      p=c*c
      a=4.5d0*p
      p=7.5d0*p
      a0=c
      b0=1.d0+a
      a1=-3.d0*c
      b1=1.d0-p
   20 i=4*k
      k=k+1
      di=i
      dk=k
      a=a*c*(di+1.d0/dk)
      p=p*c*(di-3.d0/dk)
      a0=a0+a
      a1=a1-p
      i=4*k
      k=k+1
      di=i
      dk=k
      a=a*c*(di+1.d0/dk)
      p=p*c*(di-3.d0/dk)
      b0=b0+a
      b1=b1-p
      if(dabs(dble(p))+dabs(dimag(p)).gt.1.e-16.and.k.lt.j) go to 20
      cc=dcos(x)
      s=dsin(x)
      ss=s+cc
      cc=s-cc
      th=dtanh(y)
      ch=0.56418958354776d0*dcosh(y)/dsqrt(r)
      s=phi/2.d0
      a=ch*dcmplx(dcos(s),-dsin(s))
      j0=b0*dcmplx(ss,-cc*th)
      j0=j0+a0*dcmplx(-ss*th,cc)
      j1=b1*dcmplx(cc,ss*th)
      j1=j1-a1*dcmplx(cc*th,ss)
      if (l.eq.1) j1=-j1
      c=j0/j1
      j0=j0*a
      j1=j1*a
      return
   30 a=z/2.d0
      c=-a*a
      j0=1.d0
      j1=1.d0
      p=1.d0
      k=1
   40 p=p*c
      j0=j0+p
      k=k+1
      dk=k
      p=p/dk
      j1=j1+p
      p=p/dk
      if (dabs(dble(p))+dabs(dimag(p)).gt.1.e-16) go to 40
      j1=j1*a
      c=j0/j1

      return
      end

c/////
      subroutine hankel (zz,h0,h1,c,ind)
c
c     computation of hankel funtion
c     z  : complex argument,-3.1415...le.arg(z).le. 3.1415...
c     h0 : hankel funtion of ind'th kind and zero order
c     h1 : hankel funtion of ind'th kind and first order
c     ind: 1 for first kind Hankel function
c          2 for second kind Hankel function
c     c  : h0/h1
c
      implicit double precision (a-h,o-z)
      complex*16      z,h0,h1,c,a,e,e2,zh,p, g0,g1,zz

      eulernum=0.577215664901533d0
      z=zz
      sg=1.0d0
      if (ind.eq.1) sg=-1.0d0
      x=dble(z)
      y=dimag(z)*sg
      r=dsqrt(x*x+y*y)
      phi=datan2(y,x)
      if (r.le.10.0d0) go to 60
      j=2.0d0*r
      ii=3
      h0=0.0d0
      h1=0.0d0
      if (phi.lt.1.57d0) go to 20
      ii=0
      x=-x
      z=-z
      sg=-sg
      phi=3.14159265358979d0-phi
   10 y=-y
      phi=-phi
      sg=-sg
      ii=ii+1
   20 g0=2.0d0*h0
      g1=2.0d0*h1
      c=dcmplx(0.0d0,sg*0.125d0)/z
      k=2
      p=c*c
      a=4.5d0*p
      p=7.5d0*p
      h0=1.d0+c+a
      h1=1.d0-3.d0*c-p
   30 i=4*k
      k=k+1
      di=i
      dk=k
      a=a*c*(di+1.d0/dk)
      p=p*c*(di-3.d0/dk)
      h0=h0+a
      h1=h1-p
      if(dabs(dble(p))+dabs(dimag(p)).gt.1.0e-16.and.k.lt.j) go to 30
      c=h0/h1*dcmplx(0.0d0,-sg)
      ar=(0.785398163397448d0-x-phi/2.d0)*sg
      e=0.797884560802865d0/dsqrt(r)*dexp(y)*dcmplx(dcos(ar),dsin(ar))
      if (x.ne.0.d0.or.phi.gt.0.d0) go to 40
      e=dcmplx(0.0d0,dimag(e))
   40 h0=h0*e
      h1=h1*e*dcmplx(0.0d0,sg)
      h0=h0+g0
      h1=h1+g1
      go to (10,50,90), ii
   50 h1=-h1
      go to 80
   60 zh=z/2.d0
      c=-zh*zh
      e=(0.d0,.318309886183791d0)*sg
      e2=e*2.d0
      a=1.d0-e2*( eulernum +dlog(r/2.d0))+phi*0.636619772367582d0
      p=1.d0
      k=1
      h0=a
      h1=a+e*(1.d0-1.d0/c)
   70 a=a+e2/k
      p=p*c
      h0=h0+a*p
      k=k+1
      p=p/(k*k)
      h1=h1+(a*k+e)*p
      if (dabs(dble(p))+dabs(dimag(p)).gt.1.e-16) go to 70
c
      h1=h1*zh
      if (x.ne.0.d0.or.phi.gt.0.d0) go to 80
      h0=dcmplx(0.0d0,dimag(h0))
      h1=dble(h1)
   80 c=h0/h1
   90 return
      end
