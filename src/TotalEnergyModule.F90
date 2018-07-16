Module TotalEnergyModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler
   use MathParamModule, only : ZERO, CZERO, THIRD, HALF, ONE, TWO, THREE, &
                               FOUR, FIVE, SIX, PI4, TEN2m10, SQRT_PI
   use IntegerFactorsModule, only : lofk, mofk, jofk, lofj, mofj, kofj, m1m

!
   implicit none
!
public :: initTotalEnergy,          &
          computeEnergyFunctional,  &
          printTotalEnergy,         &
          getEnergyPerAtom,         &
          getPressurePerAtom,       &
          endTotalEnergy
!
private
   character (len=50) :: StopRoutine
!
   logical :: Initialized = .false.
   logical :: InitFactors = .false.
   logical :: Computed = .false.
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: NumVacancies
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: PotentialType
   integer (kind=IntKind), allocatable :: Print_Level(:)
!
   real (kind=RealKind) :: total_energy
   real (kind=RealKind) :: energy_xc
   real (kind=RealKind) :: energy_coul
   real (kind=RealKind) :: energy_kinetic
   real (kind=RealKind) :: pressure
#ifdef DEBUG_EPRINT
   real (kind=RealKind) :: e_array(12)
#endif

   real (kind=RealKind), allocatable :: SiteEnPres(:,:,:)
!
   real (kind=RealKind), allocatable :: GlobalIndex(:)
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
!
   interface janake
      module procedure janake1, janake2
   end interface
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initTotalEnergy(nlocal,num_atoms,num_vacancies,         &
                              npola,istop,iprint)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use Atom2ProcModule, only : getGlobalIndex
   use AtomModule, only : getLocalnumSpecies
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: nlocal,num_atoms,num_vacancies
   integer (kind=IntKind), intent(in) :: npola
   integer (kind=IntKind), intent(in) :: iprint(nlocal)
!
   integer (kind=IntKind) :: id, ia, max_ns
!
   LocalNumAtoms = nlocal
   GlobalNumAtoms = num_atoms
   NumVacancies = num_vacancies
   n_spin_pola = npola
!
   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   allocate(Print_Level(nlocal))
   Print_Level(1:nlocal) = iprint(1:nlocal)
!
   total_energy = ZERO
   energy_xc = ZERO
   energy_coul = ZERO
   energy_kinetic = ZERO
   pressure = ZERO
!
   allocate( GlobalIndex(LocalNumAtoms) )
   do id = 1, LocalNumAtoms
      GlobalIndex(id) = getGlobalIndex(id)
   enddo
!
   max_ns = 0
   do id=1,nlocal
      max_ns=max(max_ns,getLocalnumSpecies(id))
   enddo
   allocate(SiteEnPres(2,max_ns,nlocal))
   InitFactors = .false.
   Initialized = .true.
   Computed = .false.
!
   end subroutine initTotalEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endTotalEnergy()
!  ===================================================================
   use IntegerFactorsModule, only : endIntegerFactors
!
   total_energy = ZERO
   energy_xc = ZERO
   energy_coul = ZERO
   energy_kinetic = ZERO
   pressure = ZERO
!
   deallocate( GlobalIndex, Print_Level, SiteEnPres )
!
   if ( InitFactors ) then
      call endIntegerFactors()
      InitFactors = .false.
   endif
   Initialized = .false.
   Computed = .false.
!
   end subroutine endTotalEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEnergyPerAtom() result (e)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: e
!
   if (.not.Initialized) then
      call ErrorHandler('getEnergyPerAtom',                      &
                        'Need to initialize TotalEnergyModule first')
   else if (.not.Computed) then
      call ErrorHandler('getEnergyPerAtom',                      &
                        'Need to call computeEnergyFunctional first')
   endif
!
   if (GlobalNumAtoms == NumVacancies) then
      e = total_energy/real(GlobalNumAtoms,RealKind)
   else
      e = total_energy/real(GlobalNumAtoms-NumVacancies,RealKind)
   endif
!
   end function getEnergyPerAtom
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPressurePerAtom() result(p)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: p
!
   if (.not.Initialized) then
      call ErrorHandler('getPressurePerAtom',                  &
                        'Need to initialize TotalEnergyModule first')
   else if (.not.Computed) then
      call ErrorHandler('getPressurePerAtom',                  &
                        'Need to call computeEnergyFunctional first')
   endif
!
   if (GlobalNumAtoms == NumVacancies) then
      p = pressure/real(GlobalNumAtoms,RealKind)
   else
      p = pressure/real(GlobalNumAtoms-NumVacancies,RealKind)
   endif
!
   end function getPressurePerAtom
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeEnergyFunctional(isMT)
!  ===================================================================
   use PotentialTypeModule, only : isFullPotential
   implicit none
!
   logical, optional :: isMT
!
   logical :: isMTon
!
   if (.not.Initialized) then
      call ErrorHandler('computeEnergyFunctional',                    &
                        'Need to initialize TotalEnergyModule first')
   endif
!
   isMTon =.false.
   if ( present(isMT) ) then
      if (isMT) then
          isMTon=.true.
      endif
   endif
!
   if ( .not.isFullPotential() .or. isMTon ) then
      if (Print_Level(1) >= 0) then
         write(6,'(/,a,/)') "     Total Energy:  -  SphPot mode "
      endif
!     ----------------------------------------------------------------
      call computeSphericalTotalEnergy()
!     ----------------------------------------------------------------
   else
      if (Print_Level(1) >= 0) then
         write(6,'(/,a,/)') "     Total Energy:  -  FullPot mode "
      endif
!     ----------------------------------------------------------------
      call computeFullTotalEnergy()
!     ----------------------------------------------------------------
   endif
!
   end subroutine computeEnergyFunctional
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeFullTotalEnergy()
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
   use GroupCommModule, only : GlobalSumInGroup
!
   use PublicTypeDefinitionsModule, only : GridStruct
   use RadialGridModule, only : getGrid
!
   use Atom2ProcModule, only : getGlobalIndex
   use AtomModule, only : getLocalAtomicNumber, getLocalnumSpecies,   &
                          getLocalSpeciesContent
!
   use PolyhedraModule, only : getVolume
!
   use ChargeDensityModule, only : getMultipoleMoment, getRhoLmax,     &
                                   getChargeDensity, getMomentDensity
!
   use PotentialModule, only : getOldPotential => getPotential
   use PotentialModule, only : getOldPotLmax => getPotLmax, getOldTruncPotLmax => getTruncPotLmax
   use PotentialModule, only : getOldTruncPotential => getTruncatedPotential
!
   use PotentialGenerationModule, only : getNewPotential => getPotential, &
                                         getPotLmax, getVshift, getVcoulomb_R0
!
!  use CoreStatesModule, only : getDeepCoreEnergy, getSemiCoreEnergy
   use CoreStatesModule, only : getDeepCoreKineticEnergy, getSemiCoreKineticEnergy
!
   use LdaCorrectionModule, only : checkLdaCorrection, getEnergyCorrection
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use InterpolationModule, only : FitInterp
!
   use SystemModule, only : setAtomEnergy
!
   use ScfDataModule, only : getSingleSiteSolverMethod
!
   use ValenceDensityModule, only : getValenceKineticEnergy
!
   implicit   none
!
   character (len=23), parameter :: sname='computeFullTotalEnergy'
!
   type (GridStruct), pointer :: Grid
!
   integer (kind=IntKind) :: jend, jend_max, jmax_max, lmax_max
   integer (kind=IntKind) :: jmax_pot, jmax_rho, kmax_pot, kmax_rho
   integer (kind=IntKind) :: lmax_rho, lmax_pot
   integer (kind=IntKind) :: lmax_pot_old, jmax_pot_old
   integer (kind=IntKind) :: lmax_prod, jmax_prod, kmax_prod
   integer (kind=IntKind) :: is, ir, na, jl, lmax, ia, nr
!
   integer (kind=IntKind), parameter :: formula_form = 1
!
   real (kind=RealKind) :: msgbuf(4)
   real (kind=RealKind) :: Zi, content
   real (kind=RealKind) :: fact1, vint_mt
   real (kind=RealKind) :: ecorr, e_coul
   real (kind=RealKind) :: e_tot(2), press(2), e_v(2), e_term1(2)
   real (kind=RealKind) :: e_kinetic(2), e_exchg(2), emu_exchg(2)
   real (kind=RealKind) :: e_terms(4,2)
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: rho_tot(:,:)
   complex (kind=CmplxKind), pointer :: mom_tot(:,:)
!
   complex (kind=CmplxKind), pointer :: rho_tmp(:,:)
   complex (kind=CmplxKind), pointer :: v_tmp(:,:)
   complex (kind=CmplxKind), pointer :: v_old(:,:)
!  complex (kind=CmplxKind), pointer :: v_new(:,:)
   complex (kind=CmplxKind), pointer :: v_coulomb(:,:)
   complex (kind=CmplxKind), pointer :: v_xc(:,:)
   complex (kind=CmplxKind), pointer :: e_xc(:,:)
   complex (kind=CmplxKind), pointer :: prod(:,:)
!
   complex (kind=CmplxKind), allocatable :: ws_rho(:)
   complex (kind=CmplxKind), allocatable :: ws_pot(:)
   complex (kind=CmplxKind), allocatable :: ws_prod(:)
!
   interface  
      subroutine computeProdExpan(n,lf,f,lg,g,lh,h)
        use KindParamModule, only : IntKind, CmplxKind
        integer (kind=IntKind), intent(in) :: n, lf, lg, lh
        complex (kind=CmplxKind), intent(in) :: f(:,:), g(:,:)
        complex (kind=CmplxKind), intent(out) :: h(:,:)
      end subroutine computeProdExpan
   end interface
!
   jend_max = 0
   jmax_max = 0
   lmax_max = 0
   do na = 1, LocalNumAtoms
       lmax = getRhoLmax(na)
       lmax_max = max(lmax_max,lmax)
       lmax = getPotLmax(na)
       lmax_max = max(lmax_max,lmax)
!      if (getSingleSiteSolverMethod() == 2) then
       if (getSingleSiteSolverMethod() >= -2) then
          lmax_max = max(lmax_max,getOldTruncPotLmax(na))
       endif
       Grid => getGrid(na)
       jend_max = max(jend_max,Grid%jend)
   enddo
   jmax_max = (lmax_max+1)*(lmax_max+2)/2
!
   if ( .not.InitFactors ) then
      call initIntegerFactors(lmax_max)
      InitFactors = .true.
   endif
   allocate(ws_rho(jend_max*jmax_max))
   allocate(ws_pot(jend_max*jmax_max))
   allocate(ws_prod(jend_max*(2*lmax_max+1)*(lmax_max+1)))
!
   SiteEnPres  = ZERO
   e_tot       = ZERO
   press       = ZERO
   do na = 1, LocalNumAtoms
      if (Print_Level(na) >= 0) then
         write(6,'(/,80(''-''))')
         write(6,'(/,20x,a)')'***************************************'
         write(6,'( 20x,a )')'*    Some Output Data from Routine    *'
         write(6,'( 20x,a )')'*  computeFullTotalEnergy and janake  *'
         write(6,'(20x,a,/)')'***************************************'
         write(6,'(80(''=''))')
         write(6,'(/,a,i5)')'Local Atom Index::',na
      endif
      Grid => getGrid(na)
      jend = Grid%jend
      r_mesh => Grid%r_mesh(1:jend)
      lmax_rho = getRhoLmax(na)
      jmax_rho = (lmax_rho+1)*(lmax_rho+2)/2
      kmax_rho = (lmax_rho+1)*(lmax_rho+1)
      lmax_pot = getPotLmax(na)
      jmax_pot = (lmax_pot+1)*(lmax_pot+2)/2
      kmax_pot = (lmax_pot+1)*(lmax_pot+1)
!
      e_terms   = ZERO
      e_term1   = ZERO
      e_coul    = ZERO
      e_kinetic = ZERO
      e_v       = ZERO
      e_exchg   = ZERO
      emu_exchg = ZERO
!
      lmax_pot_old = getOldPotLmax(na)
!     if (getSingleSiteSolverMethod() == 2) then
      if (getSingleSiteSolverMethod() >= -2) then
         lmax_pot_old = max(lmax_pot_old, getOldTruncPotLmax(na))
      endif
      jmax_pot_old = (lmax_pot_old+1)*(lmax_pot_old+2)/2
      v_tmp => aliasArray2_c(ws_pot,jend,jmax_pot_old)
      rho_tmp => aliasArray2_c(ws_rho,jend,jmax_rho)
!
      lmax_prod = lmax_rho+max(lmax_pot,lmax_pot_old)
      jmax_prod = (lmax_prod+1)*(lmax_prod+2)/2
      kmax_prod = (lmax_prod+1)**2
      prod => aliasArray2_c(ws_prod,jend,jmax_prod)
!
      do ia = 1, getLocalNumSpecies(na)
         Zi = getLocalAtomicNumber(na,ia)
         content = getLocalSpeciesContent(na,ia)
         rho_tot => getChargeDensity("TotalNew",na,ia)
         if (n_spin_pola == 2) then
            mom_tot => getMomentDensity('TotalNew',na,ia)
         endif
!
         do is=1,n_spin_pola
            v_tmp = CZERO
!
            e_v(is) = getValenceKineticEnergy(na,ia,is)  +            &
                      getSemiCoreKineticEnergy(na,ia,is) +            &
                      getDeepCoreKineticEnergy(na,ia,is)
            e_kinetic(is) = e_kinetic(is) + e_v(is)*content
            e_terms(1,is) = e_terms(1,is) + e_v(is)*content
            e_terms(4,is) = e_terms(4,is) + e_v(is)*content
            e_term1(is) = e_term1(is) + e_v(is)*content
!
            e_xc  => getNewPotential("En_Exchg",na,ia,is)
            v_xc  => getNewPotential("Exchg",na,ia,is)
!
            rho_tmp = CZERO
            if ( n_spin_pola==2 ) then
               do jl = 1,jmax_rho
                  do ir = 1,jend
                     rho_tmp(ir,jl) = (rho_tot(ir,jl) +               &
                                   (3-2*is)*mom_tot(ir,jl))*HALF
                  enddo
               enddo
            else
               do jl = 1,jmax_rho
                  do ir = 1,jend
                     rho_tmp(ir,jl) = rho_tot(ir,jl)
                  enddo
               enddo
            endif
!
            do jl = 1,jmax_pot
               do ir = 1,jend
                  v_tmp(ir,jl) =  e_xc(ir,jl)
               enddo
            enddo
            lmax_prod = lmax_rho+lmax_pot
            jmax_prod = (lmax_prod+1)*(lmax_prod+2)/2
            kmax_prod = (lmax_prod+1)**2
!           ----------------------------------------------------------
            call computeProdExpan(jend,lmax_rho,rho_tmp,lmax_pot,v_tmp,lmax_prod,prod)
            fact1 = getVolumeIntegration( na, jend, r_mesh(1:jend),      &
                                          kmax_prod, jmax_prod, 0, prod, vint_mt)
!           ----------------------------------------------------------
            e_exchg(is) = e_exchg(is) + fact1*content
            e_term1(is) = e_term1(is) + fact1*content
            e_terms(3,is) = e_terms(3,is) + fact1*content
            e_terms(4,is) = e_terms(4,is) + fact1*content
!
            do jl = 1,jmax_pot
               do ir = 1,jend
                  v_tmp(ir,jl) =  e_xc(ir,jl)-v_xc(ir,jl)
               enddo
            enddo
!           ----------------------------------------------------------
            call computeProdExpan(jend,lmax_rho,rho_tmp,lmax_pot,v_tmp,lmax_prod,prod)
            fact1 = getVolumeIntegration( na, jend, r_mesh(1:jend),      &
                                          kmax_prod, jmax_prod, 0, prod, vint_mt)
!           ----------------------------------------------------------
            emu_exchg(is) = emu_exchg(is) + fact1*content
!
         enddo
!
         v_coulomb => getNewPotential("Coulomb",na,ia,1)
!
         do jl = 1,jmax_pot
            do ir = 1,jend
               v_tmp(ir,jl) = v_coulomb(ir,jl)
            enddo
         enddo
         lmax_prod = lmax_rho+lmax_pot
         jmax_prod = (lmax_prod+1)*(lmax_prod+2)/2
         kmax_prod = (lmax_prod+1)**2
!        -------------------------------------------------------------
         call computeProdExpan(jend,lmax_rho,rho_tot,lmax_pot,v_tmp,lmax_prod,prod)
         fact1 = HALF*getVolumeIntegration( na, jend, r_mesh(1:jend),    &
                                            kmax_prod, jmax_prod, 0, prod, vint_mt)
!        -------------------------------------------------------------
         e_coul = e_coul + fact1*content
         e_terms(2,1) = e_terms(2,1) + fact1*content
         e_terms(4,1) = e_terms(4,1) + fact1*content
!
         fact1 = - HALF*Zi*getVcoulomb_R0(na)
         e_coul = e_coul + fact1*content
         e_terms(2,1) = e_terms(2,1) + fact1*content
         e_terms(4,1) = e_terms(4,1) + fact1*content
!
         do is = 1,n_spin_pola
!
            SiteEnPres(1,ia,na) = SiteEnPres(1,ia,na) + e_term1(is)
            SiteEnPres(2,ia,na) = SiteEnPres(2,ia,na) + &
                                  TWO*e_kinetic(is) - THREE*emu_exchg(is)
            e_tot(is) = e_tot(is) + e_term1(is)
            press(is) = press(is) + TWO*e_kinetic(is) - THREE*emu_exchg(is)
!
            if (Print_Level(na) >= 0) then
               write(6,'(a,i5)') "TotalEnergy:: Spin   = ", is
               write(6,'(a,f18.8)') "TotalEnergy:: e_kinetic+e_exchg = ", e_term1(is)
               write(6,'(a,f18.8)') "TotalEnergy:: e_rho             = ", e_kinetic(is)-e_v(is)
               write(6,'(a,f18.8)') "TotalEnergy:: e_vsum            = ", e_v(is)
               write(6,'(a,f18.8)') "TotalEnergy:: emu_exchg         = ", emu_exchg(is)
               write(6,'(a,f18.8)') "TotalEnergy:: pressure          = ", press(is)
               write(6,'(a,f18.8)') "TotalEnergy:: e_kinetic         = ", e_kinetic(is)
               write(6,'(a,f18.8)') "TotalEnergy:: e_exchg           = ", e_exchg(is)
            endif
         enddo
         if (Print_Level(na) >= 0) then
            write(6,'(a,f18.8)')    "TotalEnergy:: e_coul            = ", e_coul
            if (n_spin_pola == 1) then
               write(6,'(a,f18.8)') 'e_kinetic+e_exchg+e_coul        = ', &
                 e_kinetic(1)+e_exchg(1)+e_coul
            else
               write(6,'(a,f18.8)') 'e_kinetic+e_exchg+e_coul        = ', &
                 e_kinetic(1)+e_kinetic(2)+e_exchg(1)+e_exchg(2)+e_coul
            endif
         endif
!
         SiteEnPres(1,ia,na) = SiteEnPres(1,ia,na) + e_coul
         SiteEnPres(2,ia,na) = SiteEnPres(2,ia,na) + e_coul
         e_tot(1) = e_tot(1) + e_coul
         press(1) = press(1) + e_coul
!
!        ============================================================
!        Check if energy correction (e.g. LDA+U) is needed
!        ============================================================
         if ( checkLdaCorrection(na,ia) ) then ! Temporary fix ........
            ecorr = getEnergyCorrection(na,ia)
            SiteEnPres(1,ia,na) = SiteEnPres(1,ia,na) + ecorr
            e_tot(1) = e_tot(1) + ecorr
            e_terms(4,1) = e_terms(4,1) + ecorr
         endif
      enddo
   enddo
!
   nullify(rho_tot, mom_tot, rho_tmp, v_tmp, prod)
   deallocate(ws_rho, ws_pot, ws_prod)
!  ===================================================================
!  Perform global sums for the Energy and pressure
!  ===================================================================
   msgbuf = ZERO
   do is = 1,n_spin_pola
      msgbuf(1) = msgbuf(1) + e_tot(is)
      msgbuf(2) = msgbuf(2) + press(is)
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,msgbuf(1:2),2)
!  -------------------------------------------------------------------
   total_energy=msgbuf(1)
   pressure=msgbuf(2)
!
   do na =1, LocalNumAtoms
      call setAtomEnergy(getGlobalIndex(na),SiteEnPres(1:2,1,na))
   enddo
   Computed = .true.
!
!  ===================================================================
   if( StopRoutine .eq. sname ) then
       call StopHandler(sname)
   else
       return
   endif
!
   end subroutine computeFullTotalEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeSphericalTotalEnergy()
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use SystemVolumeModule, only : getTotalInterstitialVolume
!
   use ScfDataModule, only : isInterstitialElectronPolarized
!
   use PublicTypeDefinitionsModule, only : GridStruct
   use RadialGridModule, only : getGrid
!
   use Atom2ProcModule, only : getGlobalIndex
   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies,   &
                          getLocalSpeciesContent
!
   use PolyhedraModule, only : getVolume
!
   use PotentialTypeModule, only : isMuffintinASAPotential,           &
                                   isMuffintinPotential, &
                                   isASAPotential
!
   use ChargeDistributionModule, only : getInterstitialElectronDensity, &
                                        getInterstitialMomentDensity,   &
                                        getGlobalMTSphereElectronTable, &
                                        getGlobalVPCellElectronTable
!
   use PotentialModule, only : getOldSphPotr => getSphPotr
!
   use PotentialGenerationModule, only : getNewSphPotr => getSphPotr
!
   use CoreStatesModule, only : getDeepCoreEnergy, getSemiCoreEnergy,   &
                                getDeepCoreDensity, getSemiCoreDensity
!
   use ExchCorrFunctionalModule, only : getExchCorrPot, getExchCorrEnDen
!
!  use DataServiceCenterModule, only : getDataStorage, RealMark
   use ChargeDensityModule, only : getSphChargeDensity, getSphMomentDensity
!
   use LdaCorrectionModule, only : checkLdaCorrection, getEnergyCorrection
!
   use SystemModule, only : setAtomEnergy
!
   use ValenceDensityModule, only : getBandEnergy, getSphRho
!
   implicit   none
!
   character (len=27), parameter :: sname='computeSphericalTotalEnergy'
!
   type (GridStruct), pointer :: Grid
!
   integer (kind=IntKind) :: jmt
   integer (kind=IntKind) :: is, ir, na, j, ia
!
   integer (kind=IntKind), parameter :: janak_form = 2
!
   real (kind=RealKind), pointer :: rho_tot(:)
   real (kind=RealKind), pointer :: mom_tot(:)
!
   real (kind=RealKind), allocatable, target :: mom_tmp(:)
   real (kind=RealKind), allocatable :: rho_spin(:)
   real (kind=RealKind), pointer :: vrold(:)
   real (kind=RealKind), pointer :: vrnew(:)
   real (kind=RealKind), pointer :: valden_rho(:)
   real (kind=RealKind), pointer :: valden_mom(:)
   real (kind=RealKind), pointer :: deepcore(:)
   real (kind=RealKind), pointer :: semicore(:)
   real (kind=RealKind), pointer :: vx(:)
   real (kind=RealKind), pointer :: enxc(:)
   real (kind=RealKind), pointer :: rr(:)
   real (kind=RealKind), pointer :: Qmt_Table(:)
   real (kind=RealKind), pointer :: Qvp_Table(:)
!
   real (kind=RealKind) :: msgbuf(4)
   real (kind=RealKind) :: evalsum
   real (kind=RealKind) :: ecorv
   real (kind=RealKind) :: esemv
   real (kind=RealKind) :: rmt
   real (kind=RealKind) :: omegmt
   real (kind=RealKind) :: content
   real (kind=RealKind) :: ztotss
   real (kind=RealKind) :: dq
   real (kind=RealKind) :: omega_vp
   real (kind=RealKind) :: emad, dummy
   real (kind=RealKind) :: emadp
   real (kind=RealKind) :: rhoint, mdenint
   real (kind=RealKind) :: qint, mint
   real (kind=RealKind) :: u0
   real (kind=RealKind) :: sfac
   real (kind=RealKind) :: fac
   real (kind=RealKind) :: etot_is, press_is
   real (kind=RealKind) :: etot, press, ecorr
!
   Qmt_Table => getGlobalMTSphereElectronTable()
   Qvp_Table => getGlobalVPCellElectronTable()
!
   jmt = 0
   do na = 1, LocalNumAtoms
       Grid => getGrid(na)
       jmt = max(jmt,Grid%jmt)
   enddo
   allocate(rho_spin(jmt))
   if (n_spin_pola == 1) then
      allocate( mom_tmp(jmt) )
      do ir = 1, jmt
         mom_tmp(ir) = ZERO
      enddo
   endif
!
#ifdef DEBUG_EPRINT
   e_array = ZERO
#endif
   SiteEnPres = ZERO
   etot=ZERO
   press=ZERO
   do na = 1,LocalNumAtoms
      if (Print_Level(na) >= 0) then
         write(6,'(/,80(''-''))')
         write(6,'(/,20x,a)')'***************************************'
         write(6,'( 20x,a )')'*    Some Output Data from Routine    *'
         write(6,'( 20x,a )')'*  computeEnergyFunctiona and janake  *'
         write(6,'(20x,a,/)')'***************************************'
         write(6,'(80(''=''))')
         write(6,'(/,a,i5)')'  Local Atom Index:',na
      endif
!     rho_tot => getDataStorage(na,'NewSphericalElectronDensity',jmt,RealMark)
      Grid => getGrid(na)
      jmt = Grid%jmt
      rr => Grid%r_mesh(1:jmt)
      rmt = Grid%rmt
      omega_vp = getVolume(na)
      do ia = 1, getLocalNumSpecies(na)
         ztotss = getLocalAtomicNumber(na,ia)
         content = getLocalSpeciesContent(na,ia)
         rho_tot => getSphChargeDensity('TotalNew',na,ia)
         if (n_spin_pola == 1) then
            mom_tot => mom_tmp(1:jmt)
         else
!           mom_tot => getDataStorage(na,'NewSphericalMomentDensity',jmt,RealMark)
            mom_tot => getSphMomentDensity('TotalNew',na,ia)
         endif
         do is = 1,n_spin_pola
            if ( Print_Level(na) >= 0 ) then
               write(6,'(a,i2,1x,39(''-''))')'     Spin Index :',is
            endif
            vrold => getOldSphPotr(na,ia,is)
            vrnew => getNewSphPotr(na,ia,is)
!
            valden_rho => getSphRho(na,ia,is)
!
            vx => getExchCorrPot(jmt,rho_tot,mom_tot,is)
            enxc => getExchCorrEnDen(jmt,rho_tot,mom_tot,is)
!
            deepcore => getDeepCoreDensity(na,ia,is)
            semicore => getSemiCoreDensity(na,ia,is)
            do ir = 1,jmt
               rho_spin(ir) = deepcore(ir) + semicore(ir) + valden_rho(ir)
            enddo
!
            evalsum = getBandEnergy(na,ia,is)
            ecorv = getDeepCoreEnergy(na,ia,is)
            esemv = getSemiCoreEnergy(na,ia,is)
            if (janak_form == 1) then
!              -------------------------------------------------------
               call janake(vrold,vrnew,rho_spin,deepcore,                &
                           rr,jmt,ztotss,omega_vp,                       &
                           vx,enxc,                                      &
                           evalsum,ecorv,esemv,etot_is,press_is,Print_Level(na))
!              -------------------------------------------------------
            else
!              -------------------------------------------------------
               call janake(vrold,rho_tot,rho_spin,                       &
                           rr,rmt,jmt,ztotss,omega_vp,                   &
                           vx,enxc,                                      &
                           evalsum,ecorv,esemv,etot_is,press_is,Print_Level(na))
!              -------------------------------------------------------
            endif
            dummy = etot_is*content
            SiteEnPres(1,ia,na) = SiteEnPres(1,ia,na) + dummy
            etot = etot + dummy
            dummy = press_is*content
            SiteEnPres(2,ia,na) = SiteEnPres(2,ia,na) + dummy
            press = press + dummy
         enddo
      enddo
!
!     ================================================================
!      Check if energy correction (e.g. LDA+U) is needed
!     ================================================================
      do ia = 1,getLocalNumSpecies(na)
         if ( checkLdaCorrection(na,ia) ) then
            ecorr = getEnergyCorrection(na,ia)
#ifdef DEBUG_EPRINT
            e_array(10) = e_array(10) + ecorr
#endif
            SiteEnPres(1,ia,na) = SiteEnPres(1,ia,na) + ecorr
         else
            ecorr = ZERO
         endif
      enddo
      etot = etot + ecorr
   enddo
   deallocate( rho_spin )
   nullify( rho_tot, mom_tot )
   if (n_spin_pola == 1) then
      deallocate( mom_tmp )
   endif
!  ===================================================================
!  Perform global sums for the Energy and pressure
!  ===================================================================
   msgbuf(1)=etot
   msgbuf(2)=press
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,msgbuf(1:2),2)
!  -------------------------------------------------------------------
   total_energy=msgbuf(1)
   pressure=msgbuf(2)
#ifdef DEBUG_EPRINT
   call GlobalSumInGroup(GroupID,e_array(1:10),10)
#endif
!
   rhoint = getInterstitialElectronDensity()
   if (isInterstitialElectronPolarized()) then
      mdenint = getInterstitialMomentDensity()
   else
      mdenint = ZERO
   endif
!
   if ( isASAPotential() ) then
      qint = ZERO
      mint = ZERO
   else
      qint = rhoint*getTotalInterstitialVolume()
      mint = mdenint*getTotalInterstitialVolume()
   endif
!
   if ( n_spin_pola == 1 ) then
      sfac = ONE
   else
      sfac = HALF
   endif
!
   emad = ZERO
   emadp = ZERO
   if ( isMuffintinPotential() ) then
      do is = 1,n_spin_pola
         fac = sfac*(qint+(3-2*is)*mint)
         emad = emad + fac*getExchCorrEnDen(rhoint,mdenint,is)
         emadp= emadp + fac*THREE*(getExchCorrPot(rhoint,mdenint,is) &
                                  -getExchCorrEnDen(rhoint,mdenint,is))
      enddo
      do na = 1, LocalNumAtoms
         do ia = 1,getLocalNumSpecies(na)
         SiteEnPres(1,ia,na) = SiteEnPres(1,ia,na) + emad/GlobalNumAtoms
         SiteEnPres(2,ia,na) = SiteEnPres(2,ia,na) + emadp/GlobalNumAtoms
         enddo
      enddo
   else if (isMuffintinASAPotential()) then
!     ================================================================
!     emadp is not implemented???
!     ================================================================
      do na = 1, LocalNumAtoms
         Grid => getGrid(na)
         rmt = Grid%rmt
         omegmt = PI4*THIRD*rmt**3
         j = GlobalIndex(na)
         dq = Qvp_Table(j) - Qmt_Table(j) - rhoint*(getVolume(na)-omegmt)
         do is=1,n_spin_pola
            dummy = getExchCorrPot(rhoint,mdenint,is)*dq
            do ia = 1,getLocalNumSpecies(na)
               SiteEnPres(1,ia,na) = SiteEnPres(1,ia,na) + sfac*dummy
            enddo
            emad=emad+dummy
         enddo
      enddo
      emad = sfac*emad
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,emad)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   u0 = getu0()
!  -------------------------------------------------------------------
   do na = 1, LocalNumAtoms
      do ia = 1,getLocalNumSpecies(na)
         if ( NumVacancies/=GlobalNumAtoms ) then
         SiteEnPres(1,ia,na) = SiteEnPres(1,ia,na) + u0/(GlobalNumAtoms-NumVacancies)
         SiteEnPres(2,ia,na) = SiteEnPres(2,ia,na) + u0/(GlobalNumAtoms-NumVacancies)
         else
         SiteEnPres(1,ia,na) = SiteEnPres(1,ia,na) + u0/GlobalNumAtoms
         SiteEnPres(2,ia,na) = SiteEnPres(2,ia,na) + u0/GlobalNumAtoms
         endif
      enddo
      call setAtomEnergy(getGlobalIndex(na),SiteEnPres(1:2,1,na))
   enddo
   total_energy=total_energy+u0+emad
   pressure=pressure+u0+emadp
!
#ifdef DEBUG_EPRINT
   e_array(10) = e_array(10)+u0+emad
   e_array(11) = u0
   e_array(12) = emad
#endif
!
   do na = 1, LocalNumAtoms
      if(Print_Level(na) >= 0) then
      write(6,'(5x,a)')'******************************************************'
      write(6,'(10x,''emad'',t30,''='',f22.11)')emad
      write(6,'(10x,''emadp'',t30,''='',f22.11)') emadp
      write(6,'(10x,''u0'',t30,''='',f22.11)') u0
      exit
      endif
   enddo
!
   Computed = .true.
!
!  ===================================================================
   if ( StopRoutine .eq. sname ) then
      call StopHandler(sname)
   else
      return
   endif
!
   end subroutine computeSphericalTotalEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine janake2(vrold,rhotot,rho,                               &
                      rr,rmt,jmt,ztotss,omega_vp,                     &
                      vx,enxc,                                        &
                      evalsum,ecorv,esemv,etot,press,iprint)
!  ===================================================================
   use MathParamModule, only : EIGHT
   use InterpolationModule, only : FitInterp
   use IntegrationModule, only : calIntegration
   implicit   none
!
   character (len=7), parameter :: sname = 'janake2'
!
   integer (kind=IntKind), intent(in) :: jmt, iprint
!
   integer (kind=IntKind) :: ir
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: rmt
   real (kind=RealKind), intent(in) :: vrold(jmt)
   real (kind=RealKind), intent(in) :: rho(jmt)
   real (kind=RealKind), intent(in) :: rhotot(jmt)
   real (kind=RealKind), intent(in) :: rr(jmt)
   real (kind=RealKind), intent(in) :: ztotss
   real (kind=RealKind), intent(in) :: omega_vp
   real (kind=RealKind), intent(in) :: vx(jmt)
   real (kind=RealKind), intent(in) :: enxc(jmt)
   real (kind=RealKind), intent(in) :: evalsum
   real (kind=RealKind), intent(in) :: ecorv
   real (kind=RealKind), intent(in) :: esemv
   real (kind=RealKind), intent(out) :: etot
   real (kind=RealKind), intent(out) :: press
!
   real (kind=RealKind) :: qtmp
   real (kind=RealKind) :: rspin
   real (kind=RealKind) :: ekinetic
   real (kind=RealKind) :: erho
   real (kind=RealKind) :: ezrho
   real (kind=RealKind) :: ecoulomb
!
   real (kind=RealKind), allocatable :: sqrt_r(:)
   real (kind=RealKind), allocatable :: bndint(:)
   real (kind=RealKind), allocatable :: bnd(:)
   real (kind=RealKind) :: dummy
   real (kind=RealKind) :: exchen
   real (kind=RealKind) :: pterm1
   real (kind=RealKind) :: ezpt
   real (kind=RealKind) :: tpzpt
   real (kind=RealKind) :: evssum
!
   real (kind=RealKind), parameter :: PI8 = PI4*TWO
!
   allocate( bndint(0:jmt), bnd(0:jmt), sqrt_r(0:jmt) )
!
   rspin = n_spin_pola
!
   sqrt_r(0)=ZERO
   do ir=1,jmt
       sqrt_r(ir)=sqrt(rr(ir))
   enddo
!  ==================================================================
!  calculate the zeropoint energy....................................
!  ------------------------------------------------------------------
   call zeropt(ezpt,tpzpt,omega_vp,ztotss)
!  ------------------------------------------------------------------
!
!  ==================================================================
!  calculate the energy eigenvalue sum for both valence and
!  sem-core electrons, and store it in evssum..................
!  ==================================================================
   evssum = evalsum+esemv
!
!  ==================================================================
!  Look at some terms of interest
!  ==================================================================
   do i=1,jmt
       bndint(i)=rho(i)*vrold(i)
   enddo
!  ------------------------------------------------------------------
   call FitInterp(4,sqrt_r(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,sqrt_r(0:jmt),bndint(0:jmt),bnd(0:jmt),3)
   call FitInterp(jmt+1,sqrt_r(0:jmt),bnd(0:jmt),sqrt(rmt),ekinetic,dummy)
!  ------------------------------------------------------------------
   ekinetic=PI8*ekinetic
   ekinetic=ecorv+evssum-ekinetic
!
!  ==================================================================
#ifdef No_BLAS
   do i=1,jmt
       bndint(i)=rhotot(i)
   enddo
#else
!  ------------------------------------------------------------------
   call dcopy(jmt,rhotot(1:jmt),1,bndint(1:jmt),1)
!  ------------------------------------------------------------------
#endif
!  ------------------------------------------------------------------
   call FitInterp(4,sqrt_r(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,sqrt_r(0:jmt),bndint(0:jmt),bnd(0:jmt),5)
!  ------------------------------------------------------------------
   do i=1,jmt
       bndint(i)=rho(i)*bnd(i)/rr(i)
   enddo
!  ------------------------------------------------------------------
   call FitInterp(jmt+1,sqrt_r(0:jmt),bndint(0:jmt),sqrt(rmt),qtmp,dummy)
!  ------------------------------------------------------------------
   qtmp=TWO*FOUR*PI4*PI4*qtmp-ztotss
!
!  ------------------------------------------------------------------
   call FitInterp(4,sqrt_r(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,sqrt_r(0:jmt),bndint(0:jmt),bnd(0:jmt),5)
   call FitInterp(jmt+1,sqrt_r(0:jmt),bnd(0:jmt),sqrt(rmt),erho,dummy)
!  ------------------------------------------------------------------
   erho = TWO*FOUR*PI4*PI4*erho
!
!  ==================================================================
#ifdef No_BLAS
   do i=1,jmt
       bndint(i)=rho(i)
   enddo
#else
!  ------------------------------------------------------------------
   call dcopy(jmt,rho(1:jmt),1,bndint(1:jmt),1)
!  ------------------------------------------------------------------
#endif
!  ------------------------------------------------------------------
   call FitInterp(4,sqrt_r(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,sqrt_r(0:jmt),bndint(0:jmt),bnd(0:jmt),3)
   call FitInterp(jmt+1,sqrt_r(0:jmt),bnd(0:jmt),sqrt(rmt),ezrho,dummy)
!  ------------------------------------------------------------------
!  ezrho=-PI4*FOUR*ztotss*bnd(jmt)  !-rho(1)*ztotss
   ezrho = -ezrho*PI4*FOUR*ztotss
   ecoulomb=erho+ezrho
   do i=1,jmt
       bndint(i)=rho(i)*enxc(i)
   enddo
!  ------------------------------------------------------------------
   call FitInterp(4,sqrt_r(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,sqrt_r(0:jmt),bndint(0:jmt),bnd(0:jmt),5)
   call FitInterp(jmt+1,sqrt_r(0:jmt),bnd(0:jmt),sqrt(rmt),exchen,dummy)
!  ------------------------------------------------------------------
!  exchen=TWO*PI4*bnd(jmt)
   exchen=TWO*PI4*exchen
   do i=1,jmt
       bndint(i)=rho(i)*(enxc(i)-vx(i))
   enddo
!  ------------------------------------------------------------------
   call FitInterp(4,sqrt_r(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,sqrt_r(0:jmt),bndint(0:jmt),bnd(0:jmt),5)
   call FitInterp(jmt+1,sqrt_r(0:jmt),bnd(0:jmt),sqrt(rmt),pterm1,dummy)
!  ------------------------------------------------------------------
!  pterm1=TWO*PI4*bnd(jmt)
   pterm1=TWO*PI4*pterm1
!
!  ******************************************************************
   etot = ekinetic+ecoulomb+exchen+ezpt/rspin
   press= TWO*ekinetic+ecoulomb-THREE*pterm1+tpzpt/rspin
   if(iprint >= 0) then
       write(6,'(10x,''esemv'',t30,''='',f22.11)')esemv
       write(6,'(10x,''ecorev'',t30,''='',f22.11)')ecorv
       write(6,'(10x,''evssum'',t30,''='',f22.11)')evssum
       write(6,'(10x,''Kinetic E'',t30,''='',f22.11)') ekinetic
       write(6,'(10x,''Coulomb E(rho)'',t30,''='',f22.11)') erho
       write(6,'(10x,''Coulomb E(Z_rho)'',t30,''='',f22.11)') ezrho
       write(6,'(10x,''Coulomb E'',t30,''='',f22.11)') ecoulomb
       write(6,'(10x,''Exch E'',t30,''='',f22.11)') exchen
       write(6,'(10x,''pterm1'',t30,''='',f22.11)')pterm1
       write(6,'(10x,''ezpt/spin'',t30,''='',f22.11)')ezpt/rspin
       write(6,'(10x,''tpzpt'',t30,''='',f22.11)')tpzpt/rspin
   endif
#ifdef DEBUG_EPRINT
   e_array(1) = e_array(1) + esemv 
   e_array(2) = e_array(2) + ecorv 
   e_array(3) = e_array(3) + evssum 
   e_array(4) = e_array(4) + ekinetic 
   e_array(5) = e_array(5) + erho 
   e_array(6) = e_array(6) + ezrho 
   e_array(7) = e_array(7) + ecoulomb 
   e_array(8) = e_array(8) + exchen 
   e_array(9) = e_array(9) +  ezpt/rspin
   e_array(10) = e_array(10) +  etot
#endif
!
   deallocate( bndint, bnd, sqrt_r )
!
   if( StopRoutine .eq. sname ) then
       call StopHandler(sname)
   else
       return
   endif
!
   end subroutine janake2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine janake1(vrold,vrnew,rho,corden,                         &
                      rr,jmt,ztotss,omega_vp,                         &
                      vx,enxc,                                        &
                      evalsum,ecorv,esemv,etot,press,iprint)
!  ===================================================================
   use InterpolationModule, only : FitInterp
   use IntegrationModule, only : calIntegration
   implicit   none
!
   character (len=7), parameter :: sname = 'janake1'
!
   integer (kind=IntKind), intent(in) :: jmt,iprint
!
   integer (kind=IntKind) :: ir
   integer (kind=IntKind) :: i
   integer (kind=IntKind), parameter :: iformu = 1
!
   real (kind=RealKind), intent(in) :: vrold(jmt)
   real (kind=RealKind), intent(in) :: vrnew(jmt)
   real (kind=RealKind), intent(in) :: rho(jmt)
   real (kind=RealKind), intent(in) :: corden(jmt)
   real (kind=RealKind), intent(in) :: rr(jmt)
   real (kind=RealKind), intent(in) :: ztotss
   real (kind=RealKind), intent(in) :: omega_vp
   real (kind=RealKind), intent(in) :: vx(jmt)
   real (kind=RealKind), intent(in) :: enxc(jmt)
   real (kind=RealKind), intent(in) :: evalsum
   real (kind=RealKind), intent(in) :: ecorv
   real (kind=RealKind), intent(in) :: esemv
   real (kind=RealKind), intent(out) :: etot
   real (kind=RealKind), intent(out) :: press
!
   real (kind=RealKind), allocatable :: rtmp(:)
   real (kind=RealKind), allocatable :: vmvold(:)
   real (kind=RealKind), allocatable :: rhov(:)
   real (kind=RealKind), allocatable :: derv(:)
   real (kind=RealKind), allocatable :: bndint(:)
   real (kind=RealKind), allocatable :: bnd(:)
!
   real (kind=RealKind) :: rspin
   real (kind=RealKind) :: dummy
   real (kind=RealKind) :: exchen
   real (kind=RealKind) :: coren
   real (kind=RealKind) :: valen
   real (kind=RealKind) :: correc
   real (kind=RealKind) :: cor3pv
   real (kind=RealKind) :: pterm1
   real (kind=RealKind) :: ezpt
   real (kind=RealKind) :: tpzpt
   real (kind=RealKind) :: evssum
   real (kind=RealKind) :: xcmt
!
   allocate( vmvold(1:jmt), rhov(1:jmt), derv(1:jmt),   &
             bndint(0:jmt), bnd(0:jmt), rtmp(0:jmt) )
!
   rspin = n_spin_pola
!
   rtmp(0)=ZERO
   do ir=1,jmt
       rtmp(ir)=sqrt(rr(ir))
   enddo
!  ==================================================================
!  calculate the zeropoint energy....................................
!  ------------------------------------------------------------------
   call zeropt(ezpt,tpzpt,omega_vp,ztotss)
!  ------------------------------------------------------------------
!
!  ==================================================================
!  calculate the energy eigenvalue sum for both valence and
!  sem-core electrons, and store it in evssum..................
!  ==================================================================
   evssum = evalsum+esemv
!
!  ==================================================================
!  note: rhov contains both valence and semi-core charge
!        density...............................................
!  ==================================================================
   do ir=1,jmt
       rhov(ir)=rho(ir)-corden(ir)
       vmvold(ir)=vrold(ir)-vrnew(ir)
   enddo
!
!  ==================================================================
!  calculate coren for each component where one piece of 
!  dv/dr is equal to ecorv. add this piece to derv term and 
!  calculate: coren =-int r3 dr 2pi*corden*dv/dr         iformu
!                   = int r2 dr 2pi*corden*(v-d(rv)/dr)       1
!                   = sum ec - int r2 dr 4pi*corden*d(rv)/dr  0
!  ==================================================================
!  ------------------------------------------------------------------
   call newder(vrold(1),derv,rtmp(1),jmt)
!  ------------------------------------------------------------------
   if(iformu.eq.0) then
       do i=1,jmt
           bndint(i)=-rr(i)*corden(i)*derv(i)*rtmp(i)
       enddo
!      --------------------------------------------------------------
       call FitInterp(4,rtmp(1:4),bndint(1:4),ZERO,bndint(0),dummy)
       call calIntegration(jmt+1,rtmp,bndint,bnd,1)
!      --------------------------------------------------------------
       coren=ecorv+bnd(jmt)*PI4
   else
       do i=1,jmt
           bndint(i)=rr(i)*corden(i)*(vrold(i)-HALF*rtmp(i)*derv(i))
       enddo
!      --------------------------------------------------------------
       call FitInterp(4,rtmp(1:4),bndint(1:4),ZERO,bndint(0),dummy)
       call calIntegration(jmt+1,rtmp,bndint,bnd,1)
!      --------------------------------------------------------------
       coren = bnd(jmt)*PI4
   endif
!
!  ==================================================================
!  calculate: pterm1 =-int r1 dr 4pi*valden d(r2v)/dr
!                    =-int r2 dr 4pi*valden (d(rv)/dr+v)
!  ==================================================================
   do i=1,jmt
       bndint(i)=rhov(i)*(HALF*derv(i)*rtmp(i)+vrold(i))
   enddo
!  ------------------------------------------------------------------
   call FitInterp(4,rtmp(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,rtmp,bndint,bnd,3)
!  ------------------------------------------------------------------
   pterm1=-PI4*TWO*bnd(jmt)
!
!  ==================================================================
!  calculate: valen =-int r2 dr 4pi rhov*d(rv)/dr..............
!  ==================================================================
!  ------------------------------------------------------------------
   call newder(vrnew(1),derv,rtmp(1),jmt)
!  ------------------------------------------------------------------
   do i=1,jmt
       bndint(i)=-rhov(i)*derv(i)*rtmp(i)
   enddo
!  ------------------------------------------------------------------
   call FitInterp(4,rtmp(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,rtmp,bndint,bnd,3)
!  ------------------------------------------------------------------
   valen = bnd(jmt)*PI4
!
!  ==================================================================
!  calculate: exchen = int r2 dr 4pi rho (4exc-3vxc)
!  ==================================================================
   do i=1,jmt
       bndint(i)=rr(i)*rho(i)*(FOUR*enxc(i)-THREE*vx(i))
   enddo
!  ------------------------------------------------------------------
   call FitInterp(4,rtmp(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,rtmp,bndint,bnd,3)
!  ------------------------------------------------------------------
   exchen=PI4*TWO*bnd(jmt)
!
!  ==================================================================
!  calculate terms that make the energy variational: 
!            correc = int 4pi dr r3*corden*d[vold-vnew]/dr
!                    -int 4pi dr r2*rhov*[vold-vnew]
!                   = int 4pi dr r2*corden*d[r*(vold-vnew)]/dr
!                    -int 4pi dr r2*rho*[vold-vnew]
!  calculate terms that are corrections to the pressure:
!            cor3pv = int 8pi dr r3*corden*d[vold-vnew]/dr
!                   = int 8pi dr r2*corden*d[r*(vold-vnew)]/dr
!                    -int 8pi dr r2*corden*[vold-vnew]
!  ==================================================================
!  ------------------------------------------------------------------
   call newder(vmvold,derv,rtmp(1),jmt)
!  ------------------------------------------------------------------
   do i=1,jmt
       bndint(i)=corden(i)*derv(i)*rtmp(i)
   enddo
!  ------------------------------------------------------------------
   call FitInterp(4,rtmp(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,rtmp,bndint,bnd,3)
!  ------------------------------------------------------------------
   correc=bnd(jmt)*PI4
   cor3pv=TWO*correc
   do i=1,jmt
       bndint(i)=rho(i)*vmvold(i)
   enddo
!  ------------------------------------------------------------------
   call FitInterp(4,rtmp(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,rtmp,bndint,bnd,3)
!  ------------------------------------------------------------------
   correc=correc-TWO*bnd(jmt)*PI4
   do i=1,jmt
       bndint(i)=corden(i)*vmvold(i)
   enddo
!  ------------------------------------------------------------------
   call FitInterp(4,rtmp(1:4),bndint(1:4),ZERO,bndint(0),dummy)
   call calIntegration(jmt+1,rtmp,bndint,bnd,3)
!  ------------------------------------------------------------------
   cor3pv=cor3pv-FOUR*bnd(jmt)*PI4
!
   xcmt=PI4*rr(jmt)**3*rho(jmt)*( enxc(jmt)-vx(jmt) ) 
!
!  ==================================================================
!  sum over terms to get the total energy and 3PV
!  ==================================================================
   etot = coren+valen+exchen+correc+ezpt/rspin+evssum-xcmt
   press= pterm1+cor3pv+tpzpt/rspin+TWO*evssum-xcmt
!
!  *******************************************************************
!  Major print out for current sublattice
!  *******************************************************************
   if(iprint >= 0) then
!      write(6,'(10x,''ecorv'',t30,''='',f22.11)')ecorv
!      write(6,'(10x,''esemv'',t30,''='',f22.11)')esemv
!      write(6,'(10x,''evalsum'',t30,''='',f22.11)')evalsum
       write(6,'(10x,''coren'',t30,''='',f22.11)')coren
       write(6,'(10x,''valen'',t30,''='',f22.11)')valen
       write(6,'(10x,''exchen'',t30,''='',f22.11)')exchen
       write(6,'(10x,''correc'',t30,''='',f22.11)') correc
       write(6,'(10x,''evssum'',t30,''='',f22.11)')evssum
       write(6,'(10x,''-xcmt'',t30,''='',f22.11)')-xcmt
       write(6,'(10x,''c3pv'',t30,''='',f22.11)')cor3pv
       write(6,'(10x,''pterm1'',t30,''='',f22.11)')pterm1
       write(6,'(10x,''ezpt/spin'',t30,''='',f22.11)')ezpt/rspin
       write(6,'(10x,''tpzpt'',t30,''='',f22.11)')tpzpt/rspin
   endif
!
   deallocate( vmvold, rhov, derv, bndint, bnd, rtmp )
!
   if( StopRoutine .eq. sname ) then
       call StopHandler(sname)
   else
       return
   endif
!
   end subroutine janake1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine zeropt(ezpt,tpzpt,omegws,ztotss)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: idebye(49)
   integer (kind=IntKind) :: iz
!
   real (kind=RealKind), intent(out) :: ezpt
   real (kind=RealKind), intent(out) :: tpzpt
   real (kind=RealKind), intent(in) :: omegws
   real (kind=RealKind), intent(in) :: ztotss
!
   real (kind=RealKind) :: bolts
   real (kind=RealKind) :: grune(49)
   real (kind=RealKind) :: expvol(49)
   real (kind=RealKind) :: ezero
   real (kind=RealKind) :: tpvzer
!
   data  bolts/6.33870d-06/
   data  idebye/0,0,344,1440,0,0,0,0,0,0,158,400,428,0,0,0,0,0,91,  &
                230,360,420,380,630,410,467,445,450,343,            &
                327,320,0,0,0,0,0,56,147,280,                       &
                291,275,450,350,600,480,274,225,209,108/
   data  grune/0.0d0,0.0d0,1.18d0,1.18d0,0.00d0,0.00d0,0.00d0,      &
               0.0d0,0.0d0,0.00d0,1.31d0,1.48d0,2.19d0,0.00d0,      &
               0.0d0,0.0d0,0.00d0,0.00d0,1.37d0,1.16d0,1.17d0,      &
               1.18d0,1.05d0,1.30d0,2.07d0,1.66d0,                  &
               1.93d0,1.88d0,2.00d0,2.01d0,2.00d0,                  &
               0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,1.67d0,1.00d0,    &
               0.89d0,0.83d0,1.58d0,1.60d0,2.60d0,3.20d0,2.23d0,    &
               2.28d0,2.36d0,2.23d0,2.37d0/
   data  expvol/0.0d0,143.7d0,0.0d0,54.54d0,0.0d0,0.0d0,0.0d0,      &
                0.0d0,0.0d0,0.0d0,254.5d0,151.4d0,109.9d0,0.0d0,    &
                0.0d0,0.0d0,0.0d0,  0.0d0,481.3d0,291.1d0,          &
                168.7d0,120.3d0,93.48d0,80.63d0,82.84d0,78.95d0,    &
                74.72d0,73.42d0,78.92d0,99.35d0,132.4d0,            &
                0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,598.9d0,373.6d0,      &
                194.7d0,139.9d0,119.2d0,102.7d0,97.25d0,92.54d0,    &
                93.70d0,99.67d0,111.9d0,142.9d0,179.2d0/
!
!  ===================================================================
   iz=ztotss+.1
!  write(6,'('' zeropt:: iz,deb,grun,expv,omeg:'',2i5,3d12.4)')     &
!                        iz,idebye(iz),grune(iz),expvol(iz),omegws
   if(iz.lt.1 .or. iz.gt.49) then
       ezero=ZERO
       tpvzer=ZERO
   else
       if(expvol(iz) > TEN2m10) then
           ezero=1.125d0*bolts*idebye(iz)*(expvol(iz)/omegws)**grune(iz)
       else
           ezero = ZERO
       endif
       tpvzer =3.0d0*grune(iz)*ezero
   endif
   ezpt=ezero
   tpzpt=tpvzer
!
   end subroutine zeropt
!  ===================================================================
!
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getu0() result(u0)
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use PolyhedraModule, only : getVolume
!
   use SystemModule, only : getAtomicNumber
!
   use SystemVolumeModule, only : getAtomicVPVolume
!
   use MadelungModule, only : getMadelungMatrix
!
   use PotentialTypeModule, only : isMuffintinASAPotential,           &
                                   isMuffintinPotential
!
   use ChargeDistributionModule, only : getInterstitialElectronDensity, &
                                        getGlobalOnSiteElectronTable,   &
                                        getGlobalMTSphereElectronTable, &
                                        getGlobalVPCellElectronTable
!
   use RadialGridModule, only : getGrid
!
   implicit   none
!
   integer (kind=IntKind) :: i, j, k
!
   type (GridStruct), pointer :: Grid
!
   real (kind=RealKind), pointer :: madmat(:)
   real (kind=RealKind), pointer :: Q_Table(:)
   real (kind=RealKind), pointer :: Qmt_Table(:)
   real (kind=RealKind), pointer :: Qvp_Table(:)
!
   real (kind=RealKind) :: rhoint
   real (kind=RealKind) :: rmt
   real (kind=RealKind) :: omegmt
   real (kind=RealKind) :: surfamt
   real (kind=RealKind) :: u0
   real (kind=RealKind) :: dq_mt,dq
   real (kind=RealKind) :: qsub_i, qsub_j
!
   real (kind=RealKind), parameter :: fifth=ONE/FIVE
   real (kind=RealKind), parameter :: sixfifth=SIX/FIVE
!
!  calculate muffin-tin zero potential and its contribution to the
!  Coulomb energy..................................................
!  note: lattice constant factor is included in madmat.............
!  ===================================================================
   rhoint = getInterstitialElectronDensity()
!
   Q_Table => getGlobalOnSiteElectronTable()
   Qmt_Table => getGlobalMTSphereElectronTable()
   Qvp_Table => getGlobalVPCellElectronTable()
!
   u0 = ZERO
   if(isMuffintinPotential() .or. isMuffintinASAPotential()) then
       do i=1, LocalNumAtoms
           Grid => getGrid(i)
           rmt = Grid%rmt
           j = GlobalIndex(i)
           surfamt=PI4*rmt*rmt
           omegmt=surfamt*rmt*THIRD
           qsub_j = getAtomicNumber(j)-Q_Table(j)+rhoint*getAtomicVPVolume(j)
           u0=u0+rhoint*omegmt*(-sixfifth*rhoint*omegmt+THREE*qsub_j)/rmt
           if(isMuffintinASAPotential()) then
               dq_mt = Qmt_Table(j) - getAtomicNumber(j)
               dq = Qvp_Table(j) - Qmt_Table(j) - rhoint*(getVolume(i)-omegmt)
               u0=u0+dq*(TWO*dq_mt+dq)/rmt
           endif
       enddo
   endif
!
   do k = 1, LocalNumAtoms
       madmat => getMadelungMatrix(k)
       j = GlobalIndex(k)
       qsub_j = getAtomicNumber(j)-Q_Table(j) + rhoint*getAtomicVPVolume(j)
       do i=1,GlobalNumAtoms
           qsub_i = getAtomicNumber(i)-Q_Table(i) + rhoint*getAtomicVPVolume(i)
           u0=u0+madmat(i)*qsub_j*qsub_i
       enddo
   enddo
!
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,u0)
!  -------------------------------------------------------------------
!
   end function getu0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printTotalEnergy(isMT)
!  ===================================================================
#ifdef DEBUG_EPRINT
   use PotentialTypeModule, only : isFullPotential
#endif
   implicit none
!
   logical, optional :: isMT
   logical :: isMTon
!
   if (.not.Initialized) then
      call ErrorHandler('printTotalEnergy',                    &
                        'Need to initialize TotalEnergyModule first')
   else if (.not.Computed) then
      call ErrorHandler('printTotalEnergy',                    &
                        'Need to call computeEnergyFunctional first')
   endif
!
   isMTon=.false.
   if ( present(isMT) ) then
      if (isMT) then
         isMTon=.true.
      endif
   endif
!  *******************************************************************
!  Major printout: Total Energy.......................................
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')'********************************'
   write(6,'( 24x,a )')'* Output from printTotalEnergy *'
   write(6,'(24x,a,/)')'********************************'
   write(6,'(80(''=''))')
!  write(6,'(8x,''Total Energy      ='',f17.8)') total_energy
   if (GlobalNumAtoms == NumVacancies) then
      write(6,'(8x,''Total Energy Per Atom ='',f17.8)')                        &
                     total_energy/real(GlobalNumAtoms,RealKind)
!     write(6,'(8x,''Pressure          ='',f17.8)') pressure
      write(6,'(8x,''Pressure Per Atom     ='',f17.8)')                        &
                     pressure/real(GlobalNumAtoms,RealKind)
   else
      write(6,'(8x,''Total Energy Per Atom ='',f17.8)')                        &
                     total_energy/real(GlobalNumAtoms-NumVacancies,RealKind)
!     write(6,'(8x,''Pressure          ='',f17.8)') pressure
      write(6,'(8x,''Pressure Per Atom     ='',f17.8)')                        &
                     pressure/real(GlobalNumAtoms-NumVacancies,RealKind)
   endif
   write(6,'(80(''=''))')
#ifdef DEBUG_EPRINT
   if (GlobalNumAtoms == NumVacancies) then
      e_array=e_array/real(GlobalNumAtoms,RealKind)
   else
      e_array=e_array/real(GlobalNumAtoms-NumVacancies,RealKind)
   endif
   if (.not.isFullPotential() .or. isMTon) then
      write(6,'(a)') "#E_PRINT:: Tc, E_vs, E_rho, E_zrho, E_Coul, E_XC, E_Vmad, U, Etot"
      write(6,'(a,9(1x,f17.8))') "#E_PRINT::", e_array(4),e_array(3),   &
            e_array(5), e_array(6), e_array(7), e_array(8), e_array(12),&
            e_array(11), e_array(10)
   else
      write(6,'(a)') "#E_PRINT:: Tc, E_rho, E_vs, E_XC, E_Coul, E_Vmad, U, Etot"
      write(6,'(a,8(1x,f17.8))') "#E_PRINT::", e_array(1:8)
   endif
#endif

!  *******************************************************************
!
   end subroutine printTotalEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calFactors(lmax)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jl, kl, l, m, jmax, kmax, n
!
   if ( InitFactors ) then
      if ( (lmax+2)*(lmax+1)/2 > size(lofj) ) then
         deallocate( lofj, mofj, kofj, lofk, mofk, jofk, m1m )
      else
         return
      endif
   endif
!
   jmax = (lmax+2)*(lmax+1)/2
   kmax = (lmax+1)*(lmax+1)
   allocate( lofj(1:jmax), mofj(1:jmax), kofj(1:jmax) )
   allocate( lofk(1:kmax), mofk(1:kmax), jofk(1:kmax), m1m(-lmax:lmax) )
!
   jl=0
   kl=0
   m1m(0) = 1
   do l = 0, lmax
      if (l/=0) then
         m1m(l) = m1m(l-1)*(-1)
         m1m(-l) = m1m(l)
      endif
      n=(l+1)*(l+2)/2-l
      do m = -l, l
         kl = kl + 1
         lofk(kl) = l
         mofk(kl) = m
         jofk(kl) = n+abs(m)
         if (m >= 0) then
            jl = jl + 1
            lofj(jl) = l
            mofj(jl) = m
            kofj(jl) = kl
         endif
      enddo
   enddo
   InitFactors = .true.
!
   end subroutine calFactors
!  ===================================================================
end Module TotalEnergyModule
