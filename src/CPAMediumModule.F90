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
      complex (kind=CmplxKind), pointer :: tmat(:,:)
      complex (kind=CmplxKind), pointer :: tau(:,:)
   end type CPAMatrixStruct
!
   type CPAMediumStruct
      integer (kind=IntKind) :: dsize
      integer (kind=IntKind) :: local_index
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: num_species
      complex (kind=CmplxKind), pointer :: tau(:,:)
      complex (kind=CmplxKind), pointer :: Tcpa(:,:)
      complex (kind=CmplxKind), pointer :: Tcpa_old(:,:)
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
   complex (kind=CmplxKind), allocatable, target :: Tcpa_old(:)
   complex (kind=CmplxKind), allocatable, target :: Tau(:)
   complex (kind=CmplxKind), allocatable, target :: TauA(:)
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
   subroutine initCPAMedium(cant, rel, mix_type, max_iter, cpa_mix, cpa_tol,   &
                            istop, iprint)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   use MediumMatrixModule, only : initMediumMatrix
   use MediumHostModule, only  : getNumSites, getLocalNumSites, getNumSpecies, &
                                 getLmaxKKR, getGlobalIndex, getSpeciesContent
!
   use AccelerateCPAModule, only : initAccelerateCPA
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: cant, rel, mix_type, max_iter
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: i, ig, lmaxi, kmaxi, ic, n, NumImpurities
   integer (kind=IntKind) :: aid, num, dsize
!
   real (kind=RealKind), intent(in) :: cpa_mix, cpa_tol
!
   character (len=14) :: sname = "initCPAMedium"
!
   GlobalNumSites = getNumSites()
   LocalNumSites = getLocalNumSites()
   nSpinCant = cant
!
   GroupID = getGroupID('Medium Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
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
      ig = getGlobalIndex(i)
      if (getNumSpecies(ig) > 1) then
         NumCPAMediums = NumCPAMediums + 1
         NumImpurities = NumImpurities + getNumSpecies(ig)
      endif
   enddo
!
   if (NumCPAMediums > 0) then
      allocate(CPAMedium(NumCPAMediums))
      allocate(Tcpa(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(Tcpa_old(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(Tau(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumCPAMediums))
      allocate(TauA(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant*NumImpurities))
      if (nSpinCant == 2) then
         allocate(tmat_global(kmax_kkr_max*nSpinCant*kmax_kkr_max*nSpinCant,NumImpurities))
      endif
      Tcpa = CZERO; Tau = CZERO
   endif
!
   n = 0
   NumImpurities = 0
   MediumIndex = 0
   ndim_Tmat = 0; aid = 0
   do i = 1, LocalNumSites
      ig = getGlobalIndex(i)
      num = getNumSpecies(ig)
      if (num > 1) then
         n = n + 1
         MediumIndex(i) = n
         CPAMedium(n)%num_species = num
         CPAMedium(n)%local_index = i
         CPAMedium(n)%global_index = ig
         lmaxi = getLmaxKKR(ig)
         dsize = (lmaxi+1)**2*nSpinCant
         CPAMedium(n)%dsize = dsize
         CPAMedium(n)%Tcpa => aliasArray2_c(Tcpa(ndim_Tmat+1:),dsize,dsize)
         CPAMedium(n)%Tcpa_old => aliasArray2_c(Tcpa_old(ndim_Tmat+1:),dsize,dsize)
         CPAMedium(n)%tau => aliasArray2_c(Tau(ndim_Tmat+1:),dsize,dsize)
         ndim_Tmat = ndim_Tmat + dsize*dsize
         allocate(CPAMedium(n)%CPAMatrix(num))
         do ic = 1, num
            NumImpurities = NumImpurities + 1
            CPAMedium(n)%CPAMatrix(ic)%content = getSpeciesContent(ic,ig)
            if (nSpinCant == 2) then ! Spin-canted case, use the locally allocated space
               CPAMedium(n)%CPAMatrix(ic)%tmat => aliasArray2_c(tmat_global(:,NumImpurities),dsize,dsize)
            else ! In (non-)spin-polarized case, use the Tmat space in SSSolverModule
               nullify(CPAMedium(n)%CPAMatrix(ic)%tmat)
            endif
            CPAMedium(n)%CPAMatrix(ic)%tau => aliasArray2_c(TauA(aid+1:),dsize,dsize)
            aid = aid + dsize*dsize
         enddo
      endif
   enddo
!
!  -------------------------------------------------------------------
   call initMediumMatrix(cant, rel, istop, iprint)
!  -------------------------------------------------------------------
   call initAccelerateCPA(mix_type, max_iter, cpa_mix, cpa_tol,       &
                          kmax_kkr_max, NumCPAMediums)
!  -------------------------------------------------------------------
!
   iteration = 0
   print_instruction = iprint
!
   end subroutine initCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endCPAMedium()
!  ===================================================================
   use MediumMatrixModule, only : endMediumMatrix
   use AccelerateCPAModule, only : endAccelerateCPA
!
   implicit none
!
   integer (kind=IntKind) :: n, ic
!
   deallocate(WORK0, WORK1, WORK2, MediumIndex)
!
   do n = 1, NumCPAMediums
      nullify(CPAMedium(n)%tau, CPAMedium(n)%Tcpa, CPAMedium(n)%Tcpa_old)
      do ic = 1, CPAMedium(n)%num_species
         nullify(CPAMedium(n)%CPAMatrix(ic)%tau, CPAMedium(n)%CPAMatrix(ic)%tmat)
      enddo
      deallocate(CPAMedium(n)%CPAMatrix)
   enddo
!
   if (NumCPAMediums > 0) then
      deallocate(Tcpa, Tcpa_old, Tau, TauA, CPAMedium)
      if (nSpinCant ==2) then
         deallocate(tmat_global)
      endif
   endif
!
!  -------------------------------------------------------------------
   call endMediumMatrix()
!  -------------------------------------------------------------------
   call endAccelerateCPA()
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
   subroutine computeCPAMedium(e)
!  ===================================================================
   use SSSolverModule, only : getTMatrix
!
   use MediumMatrixModule, only : calMediumMatrix
!
   use AccelerateCPAModule, only : initializeAcceleration, accelerateCPA
!
   implicit none
!
   integer (kind=IntKind) :: ia, id, n, dsize
!
!  ===================================================================
!  CPA iteration acceleration parameters...
!  There parameters are taken from the mkkrcpa code
!  ===================================================================
   integer (kind=IntKind), parameter :: ipits = 4
   real (kind=RealKind), parameter ::  ctol=1.0d-08, cmix=0.15d0, cw0=5.0d-03
!  ===================================================================
!
   real (kind=RealKind) :: err, total_err
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: cfac
   complex (kind=CmplxKind), pointer :: tainv(:,:), tcinv(:,:), tau_a(:,:), p_Tau(:,:)
!  
   if (nSpinCant == 1) then
      do n = 1, NumCPAMediums
         id = CPAMedium(n)%local_index
         do ia = 1, CPAMedium(n)%num_species
            CPAMedium(n)%CPAMatrix(ia)%tmat => getTMatrix(spin=1,site=id,atom=ia)
         enddo
      enddo
   endif
!
   do n = 1, NumCPAMediums
!     ----------------------------------------------------------------
      call averageTMatrix(n)
!     ----------------------------------------------------------------
   enddo
!
   iteration = 0
   LOOP_iter: do while (iteration < MaxIterations)
      Tcpa_old = Tcpa
      iteration = iteration + 1
!     ----------------------------------------------------------------
      call initializeAcceleration(Tcpa,ndim_Tmat,iteration)
!     ----------------------------------------------------------------
!
!     ================================================================
!     calculate Tau00, Tauij, of the medium made of tmat_c's
!     ----------------------------------------------------------------
      call calMediumMatrix(e,getSingleSiteTmat)
!     ----------------------------------------------------------------
!
!     ================================================================
!     Assume that the CPA mediums on each sub-lattice are not correlated,
!     i.e., we are taking the single site approximation.
!     ================================================================
      do n = 1, NumCPAMediums  ! For each Medium, calculate new tmat_c
!        -------------------------------------------------------------
         call iterateCPAMedium(n)
!        -------------------------------------------------------------
      enddo
!
!     ----------------------------------------------------------------
      call accelerateCPA(Tcpa,ndim_Tmat,iteration)
!     ----------------------------------------------------------------
!
      total_err = ZERO
      do n = 1, NumCPAMediums
!        -------------------------------------------------------------
         call checkCPAMedium(n,err)
!        -------------------------------------------------------------
         total_err = max(total_err,err)
         if (print_instruction >= 0) then
            write(6,'(a,2i4,2x,d15.8)')'In computeCPAMedium: iter, medium, err = ', &
                  iteration, n, err   
         endif
      enddo
!
      if (total_err < CPA_tolerance) then
         exit LOOP_iter
      endif
   enddo LOOP_iter
!
   end subroutine computeCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine averageTMatrix(n)
!  ===================================================================
   use SpinRotationModule, only : rotateLtoG
!
   use SSSolverModule, only : getTMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n ! n = CPA medium index
   integer (kind=IntKind) :: i, ic, ia, dsize
!
   complex (kind=CmplxKind), pointer :: tm1(:,:), tm2(:,:), tm0(:,:)
   complex (kind=CmplxKind), pointer :: tavr(:,:)
   complex (kind=CmplxKind) :: cfac
!
   tavr => CPAMedium(n)%Tcpa
   tavr = CZERO
   ia = CPAMedium(n)%local_index
   dsize = CPAMedium(n)%dsize
   do ic = 1, CPAMedium(n)%num_species
      tm0 => CPAMedium(n)%CPAMatrix(ic)%tmat
      if (nSpinCant == 2) then
         tm1 => getTMatrix(spin=1,site=ia,atom=ic)
         tm2 => getTMatrix(spin=2,site=ia,atom=ic)
         dsize = size(tm1,1)
!        -------------------------------------------------------------
         call rotateLtoG(i, dsize, dsize, tm1, tm2, tm0)
!        -------------------------------------------------------------
      endif
      cfac = CPAMedium(n)%CPAMatrix(ic)%content
!     ----------------------------------------------------------------
      call zaxpy(dsize*dsize,cfac,tm0,1,tavr,1)
!     ----------------------------------------------------------------
   enddo
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
      call ErrorHandler('getNumSpeciesInCPAMedium','The local atom index is out of range',site)
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
   function getCPAMatrix(site,atom,matrix_type,dsize,err) result(mat)
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
      dsize = 0
      err = 1
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
   dsize = CPAMedium(n)%dsize
   err = 0
   if (nocaseCompare(matrix_type,'Tcpa')) then
      mat => CPAMedium(n)%Tcpa
   else if (nocaseCompare(matrix_type,'Tau')) then
      if (ia == 0) then
         mat => CPAMedium(n)%tau
      else
         mat => CPAMedium(n)%CPAMatrix(ia)%tau
      endif
   else 
      nullify(mat)
      dsize = 0; err = 2
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
   use MediumMatrixModule, only : embedScatterToMedium, getTmatInvDiff, &
                                  projectMatrix, beginEmbedding, endEmbedding
!
   use MatrixInverseModule, only : MtxInv_LU
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: id, ia, dsize, i
!
   complex (kind=CmplxKind) :: cfac
   complex (kind=CmplxKind), pointer :: xmat_a(:,:), xmat_c(:,:)
   complex (kind=CmplxKind), pointer :: tmat_a(:,:), diff_tminv(:,:)
   complex (kind=CmplxKind), pointer :: xmat_proj(:,:), dmat(:,:), tc_new(:,:)
!
   id = CPAMedium(n)%local_index
   dsize = CPAMedium(n)%dsize
   xmat_a => aliasArray2_c(WORK1,dsize,dsize)
   xmat_c => aliasArray2_c(WORK2,dsize,dsize)
   xmat_c = CZERO
   do ia = 1, CPAMedium(n)%num_species
      tmat_a => getSingleSiteTmat(id,ia)
!     ================================================================
!     substitute tmat_c by tmat_a in the medium
!     ================================================================
      call beginEmbedding()
      call embedScatterToMedium(id,tmat_a,dsize)
      call endEmbedding()
!     ----------------------------------------------------------------
      diff_tminv => getTmatInvDiff(id)
!     ----------------------------------------------------------------
      call projectMatrix('L',id,diff_tminv,xmat_a,dsize)
!     ----------------------------------------------------------------
!
!     ================================================================
!     add content*Tau_a to Tau_c
!     ================================================================
      cfac = CPAMedium(n)%CPAMatrix(ia)%content
!     ----------------------------------------------------------------
      call zaxpy(dsize*dsize,cfac,xmat_a,1,xmat_c,1)
!     ----------------------------------------------------------------
   enddo
   xmat_proj => aliasArray2_c(WORK0,dsize,dsize)
!  -------------------------------------------------------------------
   call projectMatrix('L',id,xmat_c,xmat_proj,dsize)
!  -------------------------------------------------------------------
   dmat => aliasArray2_c(WORK1,dsize,dsize)
   dmat = CZERO
   do i = 1, dsize
      dmat(i,i) = CONE
   enddo
!  -------------------------------------------------------------------
   call zgemm('n','n',dsize,dsize,dsize,-CONE,xmat_proj,dsize,        &
               CPAMedium(n)%Tcpa,dsize,CONE,dmat,dsize)
!  -------------------------------------------------------------------
   call MtxInv_LU(dmat,dsize)
!  -------------------------------------------------------------------
   tc_new => aliasArray2_c(WORK0,dsize,dsize)
!  -------------------------------------------------------------------
   call zgemm('n','n',dsize,dsize,dsize,CONE,CPAMedium(n)%Tcpa,dsize, &
              dmat,dsize,CZERO,tc_new,dsize)
!  -------------------------------------------------------------------
   CPAMedium(n)%Tcpa = tc_new
!  -------------------------------------------------------------------
!
   end subroutine iterateCPAMedium
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSingleSiteTmat(site,atom,dsize) result(tmat)
!  ===================================================================
   use SpinRotationModule, only : rotateLtoG
!
   use SSSolverModule, only : getTMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom  ! Local atom index
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind) :: msize, ic, n
!
   complex (kind=CmplxKind), pointer :: tm1(:,:), tm2(:,:)
   complex (kind=CmplxKind), pointer :: tmat(:,:)  ! t-matrix
!
   n = MediumIndex(site)
   if (n > 0) then
      dsize = CPAMedium(n)%dsize
      if (atom < 0 .or. atom > CPAMedium(n)%num_species) then
         call ErrorHandler('getSingleSiteTmat','The species index is out of range',atom)
      else if (atom == 0) then ! For a CPA site, it returns Tcpa if atom = 0
         tmat => CPAMedium(n)%Tcpa
      else
         tmat => CPAMedium(n)%CPAMatrix(atom)%tmat
      endif
   else 
      if (atom == 0) then
         ic = 1   ! For a non-CPA site, index atom is irrelevant
      else
         ic = atom  ! For a CPA site, set ic to a real species index
      endif
      if (nSpinCant == 2) then
         tm1 => getTMatrix(spin=1,site=site,atom=ic,dsize=msize)
         tm2 => getTMatrix(spin=2,site=site,atom=ic)
         dsize = 2*msize
         tmat => aliasArray2_c(WORK0,dsize,dsize)
!        =============================================================
!        NOTE:
!        site is the local index for the medium host, therefore it needs to
!        find appropriate index for the following call. This will be looked
!        into further in the future (12/31/2016 - YW)
!        -------------------------------------------------------------
         call rotateLtoG(site, msize, msize, tm1, tm2, tmat)
!        -------------------------------------------------------------
      else
         tmat => getTMatrix(spin=1,site=site,atom=ic,dsize=dsize)
      endif
   endif
!
   nullify(tm1, tm2)
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
   integer (kind=IntKind) :: i, j, dsize
!
   real (kind=RealKind), intent(out) :: err
   real (kind=RealKind) :: trace
!
   complex (kind=CmplxKind), pointer :: tc_new(:,:), tc_old(:,:)
!
   tc_new => CPAMedium(n)%Tcpa
   tc_old => CPAMedium(n)%Tcpa_old
   dsize = CPAMedium(n)%dsize
!
   err = ZERO; trace = ZERO
   do j = 1, dsize
      do i = 1, dsize
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
