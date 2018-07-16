module WriteFunctionModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind   
   use ErrorHandlerModule, only : ErrorHandler
!
public :: writeFunction
!
   interface writeFunction
      module procedure writeFunction_r1, writeFunction_r2,            &
                       writeFunction_c1, writeFunction_c2
   end interface writeFunction
!
private
   integer (kind=IntKind), parameter :: funit = 101
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_r1(file_name,nrs,r,f,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: f(:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (present(rpow)) then
      if (rpow >= 0 .and. rpow < 10) then
         write(funit,'(a,i1)') ' ir        r              f(r)*r^',rpow
      else
         write(funit,'(a,i2)') ' ir        r             f(r)*r^',rpow
      endif
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,5x,d15.8)')ir,r(ir),f(ir)*r(ir)**rpow
      enddo
   else
      write(funit,'(a)') ' ir        r               f(r)'
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,5x,d15.8)')ir,r(ir),f(ir)
      enddo
   endif
!
   close(funit)
!
   end subroutine writeFunction_r1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_r2(file_name,nrs,ns,r,f,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs, ns
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: f(:,:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (ns < 1 .or. ns > 2) then
      call ErrorHandler('writeFunction_r2','Not implemented for a real function with ns other than 1 and 2',ns)
   endif
!
   if (present(rpow)) then
      if (ns == 1) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1)')' ir        r              f(r)*r^',rpow
         else
            write(funit,'(a,i2)')' ir        r             f(r)*r^',rpow
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,5x,d15.8)')ir,r(ir),f(ir,1)*r(ir)**rpow
         enddo
      else if (ns == 2) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1,a,i1,a)')' ir        r            f(r)*r^',rpow,'-up             f(r)*r^',rpow,'-down'
         else
            write(funit,'(a,i2,a,i2,a)')' ir        r           f(r)*r^',rpow,'-up            f(r)*r^',rpow,'-down'
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,5x,d15.8,6x,d15.8)')ir,r(ir),f(ir,1)*r(ir)**rpow,f(ir,2)*r(ir)**rpow
         enddo
      endif
   else
      if (ns == 1) then
         write(funit,'(a)')' ir        r               f(r)'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,5x,d15.8)')ir,r(ir),f(ir,1)
         enddo
      else if (ns == 2) then
         write(funit,'(a)')' ir        r             f(r)-up              f(r)-down'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,5x,d15.8,6x,d15.8)')ir,r(ir),f(ir,1),f(ir,2)
         enddo
      endif
   endif
!
   close(funit)
!
   end subroutine writeFunction_r2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_c1(file_name,nrs,r,f,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   complex (kind=CmplxKind), intent(in) :: f(:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (present(rpow)) then
      if (rpow >= 0 .and. rpow < 10) then
         write(funit,'(a,i1,a,i1)') ' ir        r           Re[f(r)]*r^',rpow, &
                                                '           Im[f(r)]r^]',rpow
      else
         write(funit,'(a,i2,a,i2)') ' ir        r          Re[f(r)]*r^',rpow,  &
                                                '          Im[f(r)]*r^',rpow
      endif
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f(ir)*r(ir)**rpow
      enddo
   else
      write(funit,'(a)') ' ir        r            Re[f(r)]            Im[f(r)]'
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f(ir)
      enddo
   endif
!
   close(funit)
!
   end subroutine writeFunction_c1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_c2(file_name,nrs,ns,r,f,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs, ns
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   complex (kind=CmplxKind), intent(in) :: f(:,:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (ns < 1 .or. ns > 2) then
      call ErrorHandler('writeFunction_r2','Not implemented for a real function with ns other than 1 and 2',ns)
   endif
!
   if (present(rpow)) then
      if (ns == 1) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1,a,i1)')' ir        r             Re[f(r)]*r^',rpow, &
                                                  '             Im[f(r)]*r^',rpow
         else
            write(funit,'(a,i2,a,i2)')' ir        r             Re[f(r)]*r^',rpow, &
                                                  '             Im[f(r)]*r^',rpow
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f(ir,1)*r(ir)**rpow
         enddo
      else if (ns == 2) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1,a,i1,a,i1,a,i1,a)')' ir        r            Re[f(r)]*r^',rpow,    &
                                                              '-up            Im[f(r)]*r^',rpow, &
                                                              '-up            Re[f(r)]*r^',rpow,  &
                                                              '-down          Im[f(r)]*r^',rpow,'-down'
         else
            write(funit,'(a,i2,a,i2,a)')' ir        r            Re[f(r)]*r^',rpow,  &
                                                    '-up          Re[f(r)]*r^',rpow, &
                                                    '-up          Im[f(r)]*r^',rpow, &
                                                    '-down        Im[f(r)]*r^',rpow,'-down'
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8,6x,d15.8))')ir,r(ir),f(ir,1)*r(ir)**rpow,f(ir,2)*r(ir)**rpow
         enddo
      endif
   else
      if (ns == 1) then
         write(funit,'(a)')' ir        r            Re[f(r)]            Im[f(r)]'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f(ir,1)
         enddo
      else if (ns == 2) then
         write(funit,'(a)')' ir        r           Re[f(r)]-up         Im[f(r)]-up      Re[f(r)]-down      Im[f(r)]-down'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8,6x,d15.8))')ir,r(ir),f(ir,1),f(ir,2)
         enddo
      endif
   endif
!
   close(funit)
!
   end subroutine writeFunction_c2
!  ===================================================================
end module WriteFunctionModule
