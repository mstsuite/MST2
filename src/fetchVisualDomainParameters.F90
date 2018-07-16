!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine fetchVisualDomainParameters(na,nb,nc,vcell,v0)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO
!
   use ScfDataModule, only : TableID
!
   use InputModule, only : getKeyValue, getNumKeyValues
!
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   character (len=80) :: svalue
!
   integer (kind=IntKind), intent(out) :: na, nb, nc
   integer (kind=IntKind) :: GridDim
!
   real (kind=RealKind), intent(out) :: vcell(3,3), v0(3)
   real (kind=RealKind) :: a0_vgrid
!
   if ( getKeyValue(TableID,'Visual Grid Type (0<D<4)', svalue) == 0 ) then
      read(svalue,*)GridDim
      if (GridDim < 1 .or. GridDim > 3) then
         call ErrorHandler('fetchVisualDomainParameters','Invalid visual grid type input',GridDim)
      endif
   else
      call ErrorHandler('fetchVisualDomainParameters','Visual Grid Type (0<D<4)','Not exist')
   endif
!
   if ( getKeyValue(TableID,'Origin Grid Vector', svalue) == 0 ) then
      read(svalue,*)v0
   else
      call ErrorHandler('fetchVisualDomainParameters','Origin Grid Vector','Not exist')
   endif
!
   if ( getNumKeyValues(TableID,'Grid Vector') < GridDim ) then
      call ErrorHandler('fetchVisualDomainParameters','number of grid vectors < dimension')
   else if (GridDim > 3) then
      call ErrorHandler('fetchVisualDomainParameters','dimension > 3',GridDim)
   endif
!
   vcell = ZERO
   if ( getKeyValue(TableID,'Grid Vector', 3, vcell, GridDim) /= 0 ) then
      call ErrorHandler('fetchVisualDomainParameters','Grid Vector','Data is flawed')
   endif
!
   if ( getKeyValue(TableID,'Grid Points', svalue) == 0 ) then
      read(svalue,*)na, nb, nc
   else
      call ErrorHandler('fetchVisualDomainParameters','Visual Grid Type (0<D<4)','Not exist')
   endif
!
   if ( getKeyValue(TableID,'Grid Scale', a0_vgrid) == 0 ) then
      if ( a0_vgrid > ZERO ) then
         v0  = a0_vgrid*v0
         vcell = a0_vgrid*vcell
      endif
   endif
!
   end subroutine fetchVisualDomainParameters
!  ===================================================================
