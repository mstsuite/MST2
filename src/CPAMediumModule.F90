module CPAMediumModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initCPAMedium,            &
          endCPAMedium,             &
          computeCPAMedium,         &
          getNumCPAMediums,         &
          getNumSpeciesInCPAMedium, &
          getCPAMatrix
!
private
   integer (kind=IntKind) :: MaxIterations = 50
   integer (kind=IntKind) :: iteration = 0
!
   integer (kind=IntKind) :: GlobalNumSites, LocalNumSites
   integer (kind=IntKind) :: nSpinCant, NumCPAMediums = 0
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind), allocatable :: MediumIndex(:)
   integer (kind=IntKind) :: kmax_kkr_max
   integer (kind=IntKind) :: ndim_Tmat
   integer (kind=IntKind) :: print_instruction
!
   type CPAMatrixStruct
      real (kind=RealKind) :: content
      complex (kind=CmplxKind), pointer :: tmat_a(:,:)
      complex (kind=CmplxKind), pointer :: tau_a(:,:,:)
   end type CPAMatrixStruct
!
   type CPAMediumStruct
      integer (kind=IntKind) :: dsize
      integer (kind=IntKind) :: local_index
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: num_species
      complex (kind=CmplxKind), pointer :: tau_c(:,:,:)
      complex (kind=CmplxKind), pointer :: Tcpa(:,:)
      complex (kind=CmplxKind), pointer :: TcpaInv(:,:)
      complex (kind=CmplxKind), pointer :: Tcpa_old(:,:)
      complex (kind=CmplxKind), pointer :: TcpaInv_old(:,:)
      type (CPAMatrixStruct), pointer :: CPAMatrix(:)
   end type CPAMediumStruct
!
   type (CPAMediumStruct), allocatable :: CPAMedium(:)
!
   logical :: isRelativistic = .false.
!
   character (len=50) :: stop_routine
!
   complex (kind=CmplxKind), allocatable, target :: Tcpa(:)
   complex (kind=CmplxKind), allocatable, target :: TcpaInv(:)
   complex (kind=CmplxKind), allocatable, target :: Tcpa_old(:)
   complex (kind=CmplxKind), allocatable, target :: TcpaInv_old(:)
!  complex (kind=CmplxKind), allocatable, target :: Tau(:)
!  complex (kind=CmplxKind), allocatable, target :: TauA(:)
   complex (kind=CmplxKind), allocatable, target :: Tmat_global(:,:)
   complex (kind=CmplxKind), allocatable, target :: WORK0(:), WORK1(:), WORK2(:)
!
   real (kind=RealKind) :: CPA_tolerance = TEN2m8
   real (kind=RealKind) :: CPA_alpha = 0.15d0
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initCPAMedium(cant, lmax_kkr, rel, mix_type, max_iter,  &
                            cpa_mix, cpa_tol, istop, iprint)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : getNumGroups, getGroupLabel, isGroupExisting
!
   use CrystalMatrixModule, only : initCrystalMatrix
   use MediumHostModule, only  : getNumSites, getLocalNumSites, getNumSpecies, &
                                 getLmaxKKR, getGlobalSiteIndex, getSpeciesContent
!
   use AccelerateCPAModule, only : initAccelerateCPA
!
   use EmbeddedClusterModule, only : initEmbeddedCluster
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: cant, rel, mix_type, max_iter
   integer (kind=IntKind), intent(in) :: iprint(:)
   integer (kind=IntKind), intent(in) :: lmax_kkr(:)
   integer (kind=IntKind) :: i, ig, kmaxi, ic, n, NumImpurities
   integer (kind=IntKind) :: aid, num, dsize, lmaxi, nsize
!
   real (kind=RealKind), intent(in) :: cpa_mix, cpa_tol
!
   character (len=14) :: sname = "initCPAMedium"
!
   GlobalNumSites = getNumSites()
   LocalNumSites = getLocalNumSites()
   nSpinCant = cant
!
   if (isGroupExisting('Medium Cell')) then
      GroupID = getGroupID('Medium Cell')
      NumPEsInGroup = getNumPEsInGroup(GroupID)
      MyPEinGroup = getMyPEinGroup(GroupID)
   else
      GroupID = -1
      NumPEsInGroup = 1
      MyPEinGroup = 0
   endif
!
   if (rel>1) then
      isRelativistic = .true.
   else
      isRelativistic = .false.
   endif
!
   MaxIterations = max_iter
   CPA_tolerance = cpa_tol
   CPA_alpha = cpa_mix
!
   kmax_kkr_max = 0
   do ig = 1, GlobalNumSites
      lmaxi = getLmaxKKR(ig)
      kmaxi = (lmaxi+1)**2
      kmax_kkr_max = max( kmax_kkr_max,kmaxi )
   enddo
!
   allocate(WORK0(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant))
   allocate(WORK1(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant))
   allocate(WORK2(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant))
   allocate(MediumIndex(LocalNumSites))
!
   NumCPAMediums = 0
   NumImpurities = 0
   do i = 1, LocalNumSites
      ig = getGlobalSiteIndex(i)
      if (getNumSpecies(ig) > 1) then
         NumCPAMediums = NumCPAMediums + 1
         NumImpurities = NumImpurities + getNumSpecies(ig)
      endif
   enddo
!
   if (NumCPAMediums > 0) then
      allocate(CPAMedium(NumCPAMediums))
      allocate(Tcpa(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(TcpaInv(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(Tcpa_old(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(TcpaInv_old(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
!     allocate(Tau(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
!     allocate(TauA(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumImpurities))
      if (nSpinCant == 2) then
         allocate(Tmat_global(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant,NumImpurities))
      endif
      Tcpa = CZERO; TcpaInv = CZERO
!     Tau = CZERO; TauA=CZERO
   endif
!
   n = 0
   NumImpurities = 0
   MediumIndex = 0
   ndim_Tmat = 0; aid = 0
   do i = 1, LocalNumSites
      ig = getGlobalSiteIndex(i)
      num = getNumSpecies(ig)
      if (num > 1) then
         n = n + 1
         MediumIndex(i) = n
         CPAMedium(n)%num_species = num
         CPAMedium(n)%local_index = i
         CPAMedium(n)%global_index = ig
         lmaxi = getLmaxKKR(ig)
         dsize = (lmaxi+1)**2
         CPAMedium(n)%dsize = dsize
         nsize = dsize*nSpinCant
         CPAMedium(n)%Tcpa => aliasArray2_c(Tcpa(ndim_Tmat+1:),nsize,nsize)
         CPAMedium(n)%TcpaInv => aliasArray2_c(TcpaInv(ndim_Tmat+1:),nsize,nsize)
         CPAMedium(n)%Tcpa_old => aliasArray2_c(Tcpa_old(ndim_Tmat+1:),nsize,nsize)
         CPAMedium(n)%TcpaInv_old => aliasArray2_c(TcpaInv_old(ndim_Tmat+1:),nsize,nsize)
!        CPAMedium(n)%tau_c => aliasArray2_c(Tau(ndim_Tmat+1:),dsize,dsize,nSpinCant*nSpinCant)
         nullify(CPAMedium(n)%tau_c)
         ndim_Tmat = ndim_Tmat + nsize*nsize
         allocate(CPAMedium(n)%CPAMatrix(num))
         do ic = 1, num
            NumImpurities = NumImpurities + 1
            CPAMedium(n)%CPAMatrix(ic)%content = getSpeciesContent(ic,ig)
            if (nSpinCant == 2) then ! Spin-canted case, use the locally allocated space
               CPAMedium(n)%CPAMatrix(ic)%tmat_a => aliasArray2_c(Tmat_global(:,NumImpurities),nsize,nsize)
            else ! In (non-)spin-polarized case, use the Tmat space in SSSolverModule
               nullify(CPAMedium(n)%CPAMatrix(ic)%tmat_a)
            endif
!           CPAMedium(n)%CPAMatrix(ic)%tau_a => aliasArray2_c(TauA(aid+1:),dsize,dsize,nSpinCant*nSpinCant)
            nullify(CPAMedium(n)%CPAMatrix(ic)%tau_a)
            aid = aid + nsize*nsize
         enddo
      endif
   enddo
!
   iteration = 0
   print_instruction = maxval(iprint)
!
!  -------------------------------------------------------------------
   call initCrystalMatrix(LocalNumSites, cant, lmax_kkr, rel, istop, iprint)
!  -------------------------------------------------------------------
   call initEmbeddedCluster(cant, kmax_kkr_max, istop, print_instruction)
!  -------------------------------------------------------------------
   call initAccelerateCPA(mix_type, max_iter, cpa_mix, cpa_tol,       &
                          kmax_kkr_max, NumCPAMediums)
!  -------------------------------------------------------------------
!
   end subroutine initCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endCPAMedium()
!  ===================================================================
   use CrystalMatrixModule, only : endCrystalMatrix
   use AccelerateCPAModule, only : endAccelerateCPA
   use EmbeddedClusterModule, only : endEmbeddedCluster
!
   implicit none
!
   integer (kind=IntKind) :: n, ic
!
   deallocate(WORK0, WORK1, WORK2, MediumIndex)
!
   do n = 1, NumCPAMediums
      nullify(CPAMedium(n)%tau_c, CPAMedium(n)%Tcpa, CPAMedium(n)%TcpaInv)
      nullify(CPAMedium(n)%Tcpa_old, CPAMedium(n)%TcpaInv_old)
      do ic = 1, CPAMedium(n)%num_species
         nullify(CPAMedium(n)%CPAMatrix(ic)%tau_a, CPAMedium(n)%CPAMatrix(ic)%tmat_a)
      enddo
      deallocate(CPAMedium(n)%CPAMatrix)
   enddo
!
   if (NumCPAMediums > 0) then
      deallocate(Tcpa, TcpaInv, Tcpa_old, TcpaInv_old, CPAMedium)
!     deallocate(Tau, TauA)
      if (nSpinCant ==2) then
         deallocate(Tmat_global)
      endif
   endif
!
!  -------------------------------------------------------------------
   call endAccelerateCPA()
!  -------------------------------------------------------------------
   call endEmbeddedCluster()
!  -------------------------------------------------------------------
   call endCrystalMatrix()
!  -------------------------------------------------------------------
!
   NumCPAMediums = 0
   iteration = 0
!
   end subroutine endCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumCPAMediums() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = NumCPAMediums
!
   end function getNumCPAMediums
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeCPAMedium(e,cpa_medium)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
!
   use SSSolverModule, only : getScatteringMatrix
!
   use CrystalMatrixModule, only : calCrystalMatrix
!
   use AccelerateCPAModule, only : initializeAcceleration, accelerateCPA
!
   use EmbeddedClusterModule, only : getTau
!
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: cpa_medium
   integer (kind=IntKind) :: ia, id, dsize, nsize, n
   integer (kind=IntKind) :: site_config(LocalNumSites)
!
!  ===================================================================
!  CPA iteration acceleration parameters...
!  There parameters are taken from the mkkrcpa code
!  ===================================================================
   integer (kind=IntKind), parameter :: ipits = 4
   real (kind=RealKind), parameter ::  ctol=1.0d-08, cmix=0.15d0, cw0=5.0d-03
!  ===================================================================
!
   real (kind=RealKind) :: err
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: tau_a(:,:,:)
!  
   if (present(cpa_medium)) then
      n = cpa_medium
      if (n < 1 .or. n > NumCPAMediums) then
         call ErrorHandler('computeCPAMedium','The CPA medium index is out of bound',n)
      endif
   else
      n = 1
   endif
!
   if (nSpinCant == 1) then
      id = CPAMedium(n)%local_index
      do ia = 1, CPAMedium(n)%num_species
         CPAMedium(n)%CPAMatrix(ia)%tmat_a =>                      &
            getScatteringMatrix('T-Matrix',spin=1,site=id,atom=ia)
      enddo
   else
      call ErrorHandler('computeCPAMedium',                           &
               'Spin-canting calculation needs to be checked in CPA case')
   endif
!
!  -------------------------------------------------------------------
   call averageTMatrix(n)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Note that:
!       CPAMedium(n)%Tcpa(:,:) is an alias of Tcpa(:)
!       CPAMedium(n)%Tcpa_old(:,:) is an alias of Tcpa_old(:)
!       CPAMedium(n)%TcpaInv_old(:,:) is an alias of TcpaInv_old(:)
!  ===================================================================
   site_config = 0 ! Set each atomic site to be the CPA medium site
   iteration = 0
   LOOP_iter: do while (iteration < MaxIterations)
      Tcpa_old = Tcpa
      TcpaInv_old = TcpaInv
      iteration = iteration + 1
print *,'iteration = ',iteration
!     ----------------------------------------------------------------
      call initializeAcceleration(Tcpa,ndim_Tmat,iteration)
!     ----------------------------------------------------------------
!
!     ================================================================
!     calculate Tau00, Tauij, of the medium made of tmat_c's
!     ----------------------------------------------------------------
print *,'Before calCrystalMatrix'
      call calCrystalMatrix(e,getSingleSiteTmat,use_tmat=.true.,      &
                            configuration=site_config)
print *,'After calCrystalMatrix'
!     ----------------------------------------------------------------
!
!     ================================================================
!     Assume that the CPA mediums on each sub-lattice are not correlated,
!     i.e., we are taking the single site approximation.
!     ----------------------------------------------------------------
      call iterateCPAMedium(n)
!     ----------------------------------------------------------------
!
!     ----------------------------------------------------------------
      call accelerateCPA(Tcpa,ndim_Tmat,iteration)
!     ----------------------------------------------------------------
      dsize = CPAMedium(n)%dsize
      nsize = dsize*nSpinCant
      CPAMedium(n)%TcpaInv = CPAMedium(n)%Tcpa
!     ----------------------------------------------------------------
      call MtxInv_LU(CPAMedium(n)%TcpaInv,nsize)
!     ----------------------------------------------------------------
!
!     ----------------------------------------------------------------
      call checkCPAMedium(n,err)
!     ----------------------------------------------------------------
      if (print_instruction >= 0) then
         write(6,'(a,2i4,2x,d15.8)')'In computeCPAMedium: iter, medium, err = ', &
               iteration, n, err   
      endif
!
      if (err < CPA_tolerance) then
         exit LOOP_iter
      endif
   enddo LOOP_iter
!  
!  ===================================================================
!  Calculate Tau_a for each species
!  ===================================================================
   id = CPAMedium(n)%local_index
   do ia = 1, CPAMedium(n)%num_species
      tau_a => getTau(local_id=1) ! Associated tau_a space with the space
                                  ! allocated in EmbeddedCluster module
!     ================================================================
!     Substitute one CPA medium site by a real atom. The returning
!     xmat_a is the tau_a matrix.
!     ----------------------------------------------------------------
      call substituteTcByTa(id,ia,spin=1,mat_a=tau_a(:,:,1))
!     ----------------------------------------------------------------
      CPAMedium(n)%CPAMatrix(ia)%tau_a => tau_a
   enddo
!
   end subroutine computeCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine averageTMatrix(n)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
!
   use SpinRotationModule, only : rotateLtoG
!
   use SSSolverModule, only : getScatteringMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n ! n = CPA medium index
   integer (kind=IntKind) :: i, ic, ia, dsize, nsize
!
   complex (kind=CmplxKind), pointer :: tm1(:,:), tm2(:,:), tm0(:,:)
   complex (kind=CmplxKind), pointer :: tavr(:,:)
   complex (kind=CmplxKind) :: cfac
!
   tavr => CPAMedium(n)%Tcpa
   tavr = CZERO
   ia = CPAMedium(n)%local_index
   dsize = CPAMedium(n)%dsize
   nsize = dsize*nSpinCant
   do ic = 1, CPAMedium(n)%num_species
      tm0 => CPAMedium(n)%CPAMatrix(ic)%tmat_a
      if (nSpinCant == 2) then
         tm1 => getScatteringMatrix('T-Matrix',spin=1,site=ia,atom=ic)
         tm2 => getScatteringMatrix('T-Matrix',spin=2,site=ia,atom=ic)
         nsize = size(tm1,1)
!        -------------------------------------------------------------
         call rotateLtoG(i, nsize, nsize, tm1, tm2, tm0)
!        -------------------------------------------------------------
      endif
      cfac = CPAMedium(n)%CPAMatrix(ic)%content
!     ----------------------------------------------------------------
      call zaxpy(nsize*nsize,cfac,tm0,1,tavr,1)
!     ----------------------------------------------------------------
   enddo
   CPAMedium(n)%TcpaInv = CPAMedium(n)%Tcpa
!  -------------------------------------------------------------------
   call MtxInv_LU(CPAMedium(n)%TcpaInv,nsize)
!  -------------------------------------------------------------------
!
   end subroutine averageTMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumSpeciesInCPAMedium(site) result(num)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site
   integer (kind=IntKind) :: num, n
!
   if (site < 1 .or. site > LocalNumSites) then
      call ErrorHandler('getNumSpeciesInCPAMedium',                   &
                        'The local atom index is out of range',site)
   endif
!
   n = MediumIndex(site)
   if (n == 0) then
      num = 1
   else
      num = CPAMedium(n)%num_species
   endif
!
   end function getNumSpeciesInCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCPAMatrix(matrix_type,site,atom,dsize,err) result(mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site
   integer (kind=IntKind), intent(in), optional :: atom
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind), intent(out), optional :: err
   integer (kind=IntKind) :: n, ia
!
   character (len=*), intent(in) :: matrix_type
!
   complex (kind=CmplxKind), pointer :: mat(:,:)
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (site < 1 .or. site > LocalNumSites) then
      call ErrorHandler('getCPAMatrix','The local atom index is out of range',site)
   else if (len(matrix_type) /= 3 .or. len(matrix_type) /= 4) then
      call ErrorHandler('getCPAMatrix','The matrix type is invalid',matrix_type)
   endif
!
   n = MediumIndex(site)
   if (n == 0) then
      nullify(mat)
      if (present(dsize)) then
         dsize = 0
      endif
      if (present(err)) then
         err = 1
      endif
      return
   endif
!
   if (present(atom)) then
      if (atom < 0 .or. atom > CPAMedium(n)%num_species) then
         call ErrorHandler('getCPAMatrix','The species index is out of range',atom)
      else
         ia = atom
      endif
   else
      ia = 0
   endif
!
   if (present(dsize)) then
      dsize = CPAMedium(n)%dsize
   endif
   if (present(err)) then
      err = 0
   endif
   if (nocaseCompare(matrix_type,'Tcpa')) then
      mat => CPAMedium(n)%Tcpa
   else if (nocaseCompare(matrix_type,'Tau')) then
      if (ia == 0) then
         mat => CPAMedium(n)%tau_c(:,:,1)
      else
         mat => CPAMedium(n)%CPAMatrix(ia)%tau_a(:,:,1)
      endif
   else if (nocaseCompare(matrix_type,'Tmat_a') .and. ia > 0) then
      mat => CPAMedium(n)%CPAMatrix(ia)%tmat_a
   else 
      nullify(mat)
      call WarningHandler('getCPAMatrix','The requested matrix does not exist', &
                          matrix_type)
      if (present(dsize)) then
         dsize = 0
      endif
      if (present(err)) then
         err = 2
      endif
   endif
!
   end function getCPAMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine iterateCPAMedium(n)
!  ===================================================================
   use SSSolverModule, only : getScatteringMatrix
!
   use CrystalMatrixModule, only : getCPAMediumTau => getTau
!
   use EmbeddedClusterModule, only : embedScatterInHostMedium,        &
                                     setupHostMedium, beginEmbedding, &
                                     endEmbedding
   use EmbeddedClusterModule, only : calEmbeddedSiteMatrix
!
   use MatrixInverseModule, only : MtxInv_LU
!
   use MatrixModule, only : computeAprojB
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: id, ia, dsize, i, nsize
!
   complex (kind=CmplxKind) :: cfac
   complex (kind=CmplxKind), pointer :: xmat_a(:,:), xmat_c(:,:), tau_c(:,:,:)
   complex (kind=CmplxKind), pointer :: xmat_proj(:,:), tmat_a(:,:)
!
!  -------------------------------------------------------------------
   call setupHostMedium(t_host=CPAMedium(n)%TcpaInv,t_inverse=.true., &
                        getMediumTau=getCPAMediumTau)
!  -------------------------------------------------------------------
!
   id = CPAMedium(n)%local_index
   dsize = CPAMedium(n)%dsize
   nsize = dsize*nSpinCant
   xmat_a => aliasArray2_c(WORK1,nsize,nsize)
   xmat_c => aliasArray2_c(WORK2,nsize,nsize)
   xmat_c = CZERO
   do ia = 1, CPAMedium(n)%num_species
!     ================================================================
!     Substitute one CPA medium site by a real atom. The returning
!     xmat_a is the X_a matrix if compute_X=.true., or the tau_a matrix
!     if compute_X=.false.
!     ----------------------------------------------------------------
      call substituteTcByTa(id,ia,spin=1,mat_a=xmat_a,compute_X=.true.)
!     ----------------------------------------------------------------
!
!     ================================================================
!     add content*Xmat_a to Xmat_c
!     ================================================================
      cfac = CPAMedium(n)%CPAMatrix(ia)%content
!     ----------------------------------------------------------------
      call zaxpy(nsize*nsize,cfac,xmat_a,1,xmat_c,1)
!     ----------------------------------------------------------------
   enddo
   CPAMedium(n)%tau_c => getCPAMediumTau(id)
   xmat_proj => aliasArray2_c(WORK0,nsize,nsize)
!  ===================================================================
!  This only works for non-spin-canted case.
!  -------------------------------------------------------------------
   call computeAprojB('L',nsize,xmat_c,CPAMedium(n)%tau_c(:,:,1),xmat_proj)
!  -------------------------------------------------------------------
   CPAMedium(n)%TcpaInv = CPAMedium(n)%TcpaInv - xmat_proj
   CPAMedium(n)%Tcpa = CPAMedium(n)%TcpaInv
!  -------------------------------------------------------------------
   call MtxInv_LU(CPAMedium(n)%Tcpa,nsize)
!  -------------------------------------------------------------------
!
   end subroutine iterateCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine substituteTcByTa(id,ia,spin,mat_a,compute_X)
!  ===================================================================
   use SSSolverModule, only : getScatteringMatrix
!
   use EmbeddedClusterModule, only : embedScatterInHostMedium,        &
                                     setupHostMedium, beginEmbedding, &
                                     endEmbedding
   use EmbeddedClusterModule, only : calEmbeddedSiteMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, spin
!
   logical, intent(in), optional :: compute_X
   logical :: X
!
   complex (kind=CmplxKind), intent(out) :: mat_a(:,:)
   complex (kind=CmplxKind), pointer :: tmat_a(:,:)
!
   if (present(compute_X)) then
      X = compute_X
   else
      X = .false.
   endif
!
   tmat_a => getScatteringMatrix('T-Matrix',spin=spin,site=id,atom=ia)
!
!  ===================================================================
!  substitute tmat_c by tmat_a in the medium
!  -------------------------------------------------------------------
   call beginEmbedding()
   call embedScatterInHostMedium(id,tmat_a)
   call endEmbedding()
   call calEmbeddedSiteMatrix(id,mat_a,compute_X=X)
!  -------------------------------------------------------------------
!
   end subroutine substituteTcByTa
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSingleSiteTmat(sm_type,spin,site,atom,nsize) result(tmat)
!  ===================================================================
   use SpinRotationModule, only : rotateLtoG
!
   use SSSolverModule, only : getScatteringMatrix
!
   implicit none
!
   character (len=*), intent(in) :: sm_type
   integer (kind=IntKind), intent(in), optional :: spin, site, atom  ! Local atom index
   integer (kind=IntKind), intent(out), optional :: nsize
   integer (kind=IntKind) :: msize, ic, n
!
   complex (kind=CmplxKind), pointer :: tmat(:,:)  ! t-matrix
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (.not.nocaseCompare(sm_type,'T-Matrix')) then
      call ErrorHandler('getSingleSiteTmat',                          &
           'The scattering matrix type is not recognized in this case',sm_type)
   endif
!
   n = MediumIndex(site)
   if (n > 0) then
      if (atom < 0 .or. atom > CPAMedium(n)%num_species) then
         call ErrorHandler('getSingleSiteTmat','The species index is out of range',atom)
      else if (atom == 0) then ! For a CPA site, it returns Tcpa if atom = 0
         tmat => CPAMedium(n)%Tcpa
      else
         tmat => CPAMedium(n)%CPAMatrix(atom)%tmat_a
      endif
      if (present(nsize)) then
         nsize = CPAMedium(n)%dsize*nSpinCant
      endif
   else 
      if (atom == 0) then
         ic = 1   ! For a non-CPA site, index atom is irrelevant
      else
         ic = atom  ! For a CPA site, set ic to a real species index
      endif
      if (present(spin)) then
         tmat => getScatteringMatrix('T-Matrix',spin=spin,site=site,atom=ic,dsize=msize)
      else
         tmat => getScatteringMatrix('T-Matrix',spin=1,site=site,atom=ic,dsize=msize)
      endif
      if (present(nsize)) then
         nsize = msize
      endif
   endif
!
   nullify(tmat)
!
   end function getSingleSiteTmat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkCPAMedium(n,err)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i, j, dsize, nsize
!
   real (kind=RealKind), intent(out) :: err
   real (kind=RealKind) :: trace
!
   complex (kind=CmplxKind), pointer :: tc_new(:,:), tc_old(:,:)
!
   tc_new => CPAMedium(n)%Tcpa
   tc_old => CPAMedium(n)%Tcpa_old
   dsize = CPAMedium(n)%dsize
   nsize = dsize*nSpinCant
!
   err = ZERO; trace = ZERO
   do j = 1, nsize
      do i = 1, nsize
         err = err + abs(tc_new(i,j) - tc_old(i,j))
      enddo
      trace = trace + abs(tc_old(j,j))
   enddo
   if (trace < TEN2m6) then
!     ----------------------------------------------------------------
      call ErrorHandler('checkCPAMedium','trace of Tcpa is too small',trace)
!     ----------------------------------------------------------------
   endif
!
   err = err/trace
!
   end subroutine checkCPAMedium
!  ===================================================================
end module CPAMediumModule
