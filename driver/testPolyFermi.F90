program testPolyFermi
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use SystemModule, only : getNumAtoms
!
   use ScfDataModule, only : n_spin_pola, n_spin_cant, Temperature, NumEs
   use ScfDataModule, only : istop, EvBottom, isNonRelativisticCore
!
   use AtomModule, only : getPhiLmax, getStepFuncLmax, getPotLmax
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, CZERO, SQRTm1, TEN2m8, PI4, PI
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : getVolume, getInscrSphVolume
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use OutputModule, only : getStandardOutputLevel
!
   use Atom2ProcModule, only : getLocalNumAtoms
!
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential
!
   implicit   none
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: NumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: lmax_max
   integer (kind=IntKind) :: id, ig, ne, ie
!
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:), lmax_step(:)
!
   real (kind=RealKind) :: Temp
!
   complex (kind=CmplxKind), allocatable :: xg(:), wg(:)
!
   interface
      subroutine polyfermi(T,xg,wg,ne,npoles)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: ne
         integer (kind=IntKind), intent(in), optional :: npoles
         real (kind=RealKind), intent(in) :: T
         complex (kind=CmplxKind), intent(out) :: xg(ne),wg(ne)
      end subroutine polyfermi
   end interface
!
!  -------------------------------------------------------------------
   call startProcess()
   NumAtoms = getNumAtoms()
   LocalNumAtoms=getLocalNumAtoms()
!  -------------------------------------------------------------------
!
   allocate(atom_print_level(1:LocalNumAtoms))
   allocate(lmax_step(LocalNumAtoms),lmax_pot(LocalNumAtoms))
!
   lmax_max = 0
   do id = 1, LocalNumAtoms
      lmax_pot(id) = getPotLmax(id)
      lmax_step(id)  = getStepFuncLmax(id)
      lmax_max = max(lmax_max,2*getPhiLmax(id),lmax_step(id))
      atom_print_level(id) = getStandardOutputLevel(id)
   enddo
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(lmax_max,istop,iprint)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call setupRadGridAndCell(LocalNumAtoms,lmax_max)
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot,lmax_step,               &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
!  ===================================================================
!  read potential data
!  -------------------------------------------------------------------
   call readPotential()
!  -------------------------------------------------------------------
!
   Temp = Temperature
   ne = NumEs
   write(6,'(a,f10.5)')'Temperature = ',Temp
   write(6,'(a,i5)')'Number of Gaussian points on the contour = ',ne
   allocate(xg(ne),wg(ne))
!  -------------------------------------------------------------------
   call polyfermi(Temp,xg,wg,ne)
!  -------------------------------------------------------------------
   do ie = 1, ne
      write(6,'(a,i5,2x,2d15.7,2x,2d15.7)')'ie, xg, wg = ',ie,xg(ie),wg(ie)
   enddo
!
   deallocate(atom_print_level,lmax_step,lmax_pot)
   deallocate(xg,wg)
!
!  -------------------------------------------------------------------
   call endPotential()
   call endGauntFactors()
   call endSphericalHarmonics()
   call finishProcess()
!  -------------------------------------------------------------------
   stop 'Ok'
!
end program testPolyFermi
