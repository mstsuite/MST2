module WriteMatrixModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ten2m8, ten2m12
!
public :: writeMatrix
   interface writeMatrix
      module procedure WriteMatrixr1, WriteMatrixr2, WriteMatrixr3
      module procedure WriteMatrixc1, WriteMatrixc2, WriteMatrixc3
   end interface
!
private
!
   real (kind=RealKind), parameter :: tol0 = ten2m12
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeMatrixr1(a,x,n,tol_in)
!  ===================================================================
   implicit   none
!
   character (len=*), intent(in) :: a
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
! 
   real (kind=RealKind), intent(in) :: x(n)
   real (kind=RealKind), intent(in), optional :: tol_in
   real (kind=RealKind) :: tol
! 
!  *******************************************************************
!  * writes out the non-zero elements (> 10**-8) of a N complex 
!  * matrix
!  *******************************************************************
   if (present(tol_in)) then
      tol = tol_in
   else
      tol = tol0
   endif
! 
   write(6,'(/,80(''-''))')
   write(6,'(/,27x,a)') '***************************'
   write(6,'(  27x,a )')'* Output from writeMatrix *'
   write(6,'( 27x,a,/)')'***************************'
!
   write(6,'('' TITLE: '',a)') a
   write(6,'('' Matrix(1:n). n = '',i5)') n
   write(6,'(a)')'   i         Matrix(i)'
   do i=1,n
      if (abs(x(i)) > tol) then
         write(6,'(i4,2x,d18.8)')i,x(i)
      endif
   enddo
   write(6,'(80(''=''))')
   end subroutine writeMatrixr1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeMatrixr2(a,x,n1,n2,tol_in)
!  ===================================================================
   implicit   none
!
   character (len=*), intent(in) :: a
!
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind) :: i, j
! 
   real (kind=RealKind), intent(in) :: x(n1, n2)
   real (kind=RealKind), intent(in), optional :: tol_in
   real (kind=RealKind) :: tol
! 
!  *******************************************************************
!  * writes out the non-zero elements (> 10**-8) of a N1*N2 complex 
!  * matrix
!  *******************************************************************
   if (present(tol_in)) then
      tol = tol_in
   else
      tol = tol0
   endif
! 
   write(6,'(/,80(''-''))')
   write(6,'(/,27x,a)') '***************************'
   write(6,'(  27x,a )')'* Output from writeMatrix *'
   write(6,'( 27x,a,/)')'***************************'
!
   write(6,'('' TITLE: '',a)') a
   write(6,'('' Matrix(1:n1, 1:n2). n1, n2 = '',2i5)') n1, n2
   write(6,'(a)')'   i   j       Matrix(i,j)'
   do j=1,n2
      do i=1,n1
         if (abs(x(i,j)) > tol) then
            write(6,'(2i4,2x,d18.8)')i,j,x(i,j)
         endif
      enddo
   enddo
   write(6,'(80(''=''))')
   end subroutine writeMatrixr2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeMatrixr3(a,x,n1,n2,tol_in)
!  ===================================================================
   implicit   none
!
   character (len=*), intent(in) :: a
! 
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind) :: i, j, k
! 
   real (kind=RealKind), intent(in) :: x(n1*n2)
   real (kind=RealKind), intent(in), optional :: tol_in
   real (kind=RealKind) :: tol
! 
!  *******************************************************************
!  * writes out the non-zero elements (> 10**-8) of a N1*N2 complex 
!  * matrix
!  *******************************************************************
   if (present(tol_in)) then
      tol = tol_in
   else
      tol = tol0
   endif
! 
   write(6,'(/,80(''-''))')
   write(6,'(/,27x,a)') '***************************'
   write(6,'(  27x,a )')'* Output from writeMatrix *'
   write(6,'( 27x,a,/)')'***************************'
!
   write(6,'('' TITLE: '',a)') a
   write(6,'('' Matrix(1:n1, 1:n2). n1, n2 = '',2i5)') n1, n2
   write(6,'(a)')'   i   j       Matrix(i,j)'
   do j=1,n2
      k=(j-1)*n1
      do i=1,n1
         if (abs(x(i+k)) > tol) then
            write(6,'(2i4,2x,d18.8)')i,j,x(i+k)
         endif
      enddo
   enddo
   write(6,'(80(''=''))')
   end subroutine writeMatrixr3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeMatrixc1(a,x,n,tol_in)
!  ===================================================================
   implicit   none
!
   character (len=*), intent(in) :: a
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
! 
   complex (kind=CmplxKind), intent(in) :: x(n)
   real (kind=RealKind), intent(in), optional :: tol_in
   real (kind=RealKind) :: tol
! 
!  *******************************************************************
!  * writes out the non-zero elements (> 10**-8) of a N complex 
!  * matrix
!  *******************************************************************
   if (present(tol_in)) then
      tol = tol_in
   else
      tol = tol0
   endif
! 
   write(6,'(/,80(''-''))')
   write(6,'(/,27x,a)') '***************************'
   write(6,'(  27x,a )')'* Output from writeMatrix *'
   write(6,'( 27x,a,/)')'***************************'
!
   write(6,'('' TITLE: '',a)') a
   write(6,'('' Matrix(1:n). n = '',i5)') n
   write(6,'(a)')'   i                  Matrix(i)'
   do i=1,n
      if (abs(x(i)) > tol) then
         write(6,'(i4,2x,2d18.8)')i,x(i)
      endif
   enddo
   write(6,'(80(''=''))')
   end subroutine writeMatrixc1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeMatrixc2(a,x,n1,n2,tol_in)
!  ===================================================================
   implicit   none
!
   character (len=*), intent(in) :: a
!
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind) :: i, j
! 
   complex (kind=CmplxKind), intent(in) :: x(n1, n2)
   real (kind=RealKind), intent(in), optional :: tol_in
   real (kind=RealKind) :: tol
! 
!  *******************************************************************
!  * writes out the non-zero elements (> 10**-8) of a N1*N2 complex 
!  * matrix
!  *******************************************************************
   if (present(tol_in)) then
      tol = tol_in
   else
      tol = tol0
   endif
! 
   write(6,'(/,80(''-''))')
   write(6,'(/,27x,a)') '***************************'
   write(6,'(  27x,a )')'* Output from writeMatrix *'
   write(6,'( 27x,a,/)')'***************************'
!
   write(6,'('' TITLE: '',a)') a
   write(6,'('' Matrix(1:n1, 1:n2). n1, n2 = '',2i5)') n1, n2
   write(6,'(a)')'   i   j                 Matrix(i,j)'
   do j=1,n2
      do i=1,n1
         if (abs(x(i,j)) > tol) then
            write(6,'(2i4,2x,2d18.8)')i,j,x(i,j)
         endif
      enddo
   enddo
   write(6,'(80(''=''))')
   end subroutine writeMatrixc2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeMatrixc3(a,x,n1,n2,tol_in)
!  ===================================================================
   implicit   none
!
   character (len=*), intent(in) :: a
! 
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind) :: i, j, k
! 
   complex (kind=CmplxKind), intent(in) :: x(n1*n2)
   real (kind=RealKind), intent(in), optional :: tol_in
   real (kind=RealKind) :: tol
! 
!  *******************************************************************
!  * writes out the non-zero elements (> 10**-8) of a N1*N2 complex 
!  * matrix
!  *******************************************************************
   if (present(tol_in)) then
      tol = tol_in
   else
      tol = tol0
   endif
! 
   write(6,'(/,80(''-''))')
   write(6,'(/,27x,a)') '***************************'
   write(6,'(  27x,a )')'* Output from writeMatrix *'
   write(6,'( 27x,a,/)')'***************************'
!
   write(6,'('' TITLE: '',a)') a
   write(6,'('' Matrix(1:n1, 1:n2). n1, n2 = '',2i5)') n1, n2
   write(6,'(a)')'   i   j                 Matrix(i,j)'
   do j=1,n2
      k=(j-1)*n1
      do i=1,n1
         if (abs(x(i+k)) > tol) then
            write(6,'(2i4,2x,2d18.8)')i,j,x(i+k)
         endif
      enddo
   enddo
   write(6,'(80(''=''))')
   end subroutine writeMatrixc3
!  ===================================================================
end module WriteMatrixModule
