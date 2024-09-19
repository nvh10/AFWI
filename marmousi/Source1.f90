   ! program test
   !
   !
   ! implicit none
   !
   ! double complex A(6,5)
   ! double complex, allocatable :: vA(:)
   ! integer, allocatable :: ja(:)
   ! integer ia(7),nnz
   ! 
   ! !A(1,:)=(/0, 7, 0, 6, -2/)
   ! !A(2,:)=(/1, -1, 0, 0, 2/)
   ! !A(3,:)=(/3, 4, 0, 0, 3/)
   ! !A(4,:)=(/2, 5, 0, -2, 2/)
   ! !A(5,:)=(/0, 3, 0, 1, 4/)
   ! !A(6,:)=(/0, 3, 0, 1, 4/)
   !  A(1,:)=(/dcmplx(0d0,0d0), dcmplx(7d0,0d0), dcmplx(0d0,0d0), dcmplx(6d0,0d0), dcmplx(-2d0,0d0)/)
   ! A(2,:)=(/dcmplx(1d0,0d0),dcmplx(-1d0,0d0), dcmplx(0d0,0d0), dcmplx(0d0,0d0), dcmplx(2d0,0d0)/)
   ! A(3,:)=(/dcmplx(3d0,0d0), dcmplx(4d0,0d0),dcmplx(0d0,0d0),dcmplx(0d0,0d0), dcmplx(3d0,0d0)/)
   ! A(4,:)=(/dcmplx(2d0,0d0), dcmplx(5d0,0d0), dcmplx(0d0,0d0), dcmplx(-2d0,0d0), dcmplx(2d0,0d0)/)
   ! A(5,:)=(/dcmplx(0d0,0d0), dcmplx(3d0,0d0), dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(4d0,0d0)/)
   ! A(6,:)=(/dcmplx(0d0,0d0), dcmplx(3d0,0d0),dcmplx(0d0,0d0), dcmplx(1d0,0d0), dcmplx(4d0,0d0)/)
   !
   !  nnz=count(A/=dcmplx(0d0,0d0))
   !  allocate(ja(nnz),vA(nnz))
   !  
   !  
   !  call create_CSR_formartZ(A,6,5,ia,ja,vA,nnz)
   !end program
