program testProcMapping
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE
   use MPPModule, only : initMPP, endMPP, syncAllPEs
   use GroupCommModule
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter, &
                                       isDataStorageExisting, &
                                       getDataStorage, RealMark
   use InputModule, only : initInput, endInput
!
   use SystemModule, only : initSystem, endSystem
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber
!
   use ScfDataModule, only : initScfData, endScfData
   use ScfDataModule, only : getPotentialTypeParam
   use ScfDataModule, only : isReadEmesh, getEmeshFileName
   use ScfDataModule, only : isReadKmesh, getKmeshFileName
   use ScfDataModule, only : NumEs, ContourType, eGridType
   use ScfDataModule, only : NumKMeshs, kGenScheme, Kdiv, Symmetrize
   use ScfDataModule, only : isScreenKKR_LSMS, isKKR, Temperature
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType, &
                                   isFullPotential
!
   use ContourModule, only : initContour, endContour, getNumEs
!
   use BZoneModule, only : initBZone, printBZone, endBZone, getNumKs
!
   use ProcMappingModule, only : initProcMapping, endProcMapping,            &
                                 createParallelization
!
   implicit none
!
   integer (kind=IntKind) :: def_id, info_id, MyPE, NumPEs, nk, ne
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: i, ia, ie, ik
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
!
   real (kind=RealKind), pointer :: bravais(:,:)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
!
!  -------------------------------------------------------------------
   call initMPP()
   call initGroupComm()
   call initDataServiceCenter()
   call initInput()
   call readInputs(def_id,info_id)
   call initScfData(def_id)
   call initPotentialType(getPotentialTypeParam())
   call initSystem(def_id)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Initialize the contour in energy complex plane to find the total
!  number of energy mesh needed
!  ===================================================================
   if (isReadEmesh()) then
!     ----------------------------------------------------------------
      call initContour(getEmeshFileName(), 'none', -1)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call initContour( ContourType, eGridType, NumEs, Temperature, 'none', -1)
!     ----------------------------------------------------------------
   endif
!  
!  -------------------------------------------------------------------
   ne = getNumEs()
!  -------------------------------------------------------------------
!  
!  ===================================================================
!  After the number of energy points is obtained, we need to call endContour
!  since initContour will be called again in ValenceStatesModule within
!  SCF loop
!  -------------------------------------------------------------------
   call endContour()
!  -------------------------------------------------------------------
!  
   if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
      bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('testProcMapping','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  
   NumAtoms = getNumAtoms()
!  
   allocate(AtomPosition(3,NumAtoms), AtomicNumber(NumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
      AtomicNumber(i)=getAtomicNumber(i)
   enddo
!  
!  ===================================================================
!  Initialize the Brillouin zone mesh for k-space integration
!  ===================================================================
   if (isKKR() .or. isScreenKKR_LSMS()) then
      if (isReadKmesh()) then
!        -------------------------------------------------------------
         call initBZone(getKmeshFileName(),'none',-1)
!        -------------------------------------------------------------
      else if (NumKMeshs > 0) then
!        -------------------------------------------------------------
         call initBZone(NumKMeshs,kGenScheme,Kdiv,Symmetrize,bravais, &
                        NumAtoms,AtomPosition,AtomicNumber,'none',-1)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call WarningHandler('main','No K mesh is initialized')
!        -------------------------------------------------------------
      endif
      call printBZone()
      nk = getNumKs()
   else
      nk = 0
   endif
!
!  ===================================================================
!  Initialize the processes mapping module that determines how the
!  parallization will be performed
!  -------------------------------------------------------------------
   call initProcMapping(NumAtoms, ne, nk, isFullPotential(), 'none', 0)
!  -------------------------------------------------------------------
   call createParallelization()
!  -------------------------------------------------------------------
!
!  ===================================================================
!
!
!  ===================================================================
!
!  -------------------------------------------------------------------
   call endProcMapping()
   call endSystem()
   call endPotentialType()
   call endScfData()
   call endInput()
   call endDataServiceCenter()
   call endGroupComm()
   call endMPP()
!  -------------------------------------------------------------------
   stop
!
end program testProcMapping
!  ===================================================================
