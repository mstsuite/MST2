module IBZRotationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO, HALF, ONE, CONE, PI, TEN2m10, SQRTm1
!
public :: initIBZRotation,          &
          endIBZRotation,           &
          computeRotationMatrix,    &
          isProperRotation,         &
          getNumIBZRotations,       &
          getNumProperIBZRotations, &
          getIBZRotationMatrix,     &
          getIBZRotationMatrix3D,   &
          printIBZRotationMatrix
!
private
!
      logical :: symmetrize = .true.
!
      character (len=12) :: lattice_name
!
      integer (kind=IntKind) :: l_only
      integer (kind=IntKind) :: lmax, kkrsz
      integer (kind=IntKind) :: NumRotations
      integer (kind=IntKind) :: NumProperRotations
      integer (kind=IntKind), parameter :: maxrot = 48
!
      real (kind=RealKind), target :: rot3d(3,3,maxrot)
      real (kind=RealKind) :: euler(3,maxrot/2)
      real (kind=RealKind), allocatable :: fact(:)
!
      complex (kind=CmplxKind), allocatable :: reflex(:,:)
      complex (kind=CmplxKind), allocatable, target :: dj(:,:,:)
      complex (kind=CmplxKind), allocatable, target :: djc(:,:,:)
!
      logical :: Initialized = .false.
!
contains
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initIBZRotation(relativistic,latsys,lmax_in,isym)
!     ================================================================
      use IntegerFactorsModule, only : isIntegerFactorsInitialized,   &
                                       initIntegerFactors
      implicit none
!
      logical, intent(in) :: relativistic
!
      character (len=*), intent(in) :: latsys
!
      integer (kind=IntKind), intent(in) :: lmax_in, isym
      integer (kind=IntKind) :: l
!
      if (relativistic) then
         l_only = 0
      else
         l_only = 1
      endif
!
      if (isym == 0) then
         symmetrize = .false.
      else
         symmetrize = .true.
      endif
!
      lattice_name = latsys
      lmax = lmax_in
      kkrsz = (lmax+1)**2
!
      if (.not.isIntegerFactorsInitialized()) then
         call initIntegerFactors(lmax)
      endif
!
      allocate( fact(0:2*lmax+1) ) ! The upper limit (2*lmax+1) may need to be increased.
      allocate( reflex(kkrsz,kkrsz) )
      allocate( dj(kkrsz,kkrsz,maxrot) )
      allocate( djc(kkrsz,kkrsz,maxrot) )
!
      rot3d = ZERO
      reflex = CZERO
      dj = CZERO
      djc = CZERO
!
!     Set factorials....................................................
      fact(0)=one
      do l=1,2*lmax+1   ! The upper limit (2*lmax+1) may need to be increased.
         fact(l)=fact(l-1)*l
      enddo
!
      Initialized = .true.
!
      end subroutine initIBZRotation
!     ================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine endIBZRotation()
!     ================================================================
      implicit none
!
      deallocate( fact, reflex, dj, djc )
!
      Initialized = .false.
!
      end subroutine endIBZRotation
!     ================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeRotationMatrix()
!     ================================================================
      implicit none
!
      integer (kind=IntKind) :: i, kl, klp
!
      if (symmetrize) then
!        -------------------------------------------------------------
         call getrot(lattice_name,l_only,NumRotations,NumProperRotations)
!        -------------------------------------------------------------
         do i = 1, NumRotations
            do kl = 1, kkrsz
               do klp = 1, kkrsz
                  djc(klp,kl,i) = conjg(dj(klp,kl,i))
               enddo
            enddo
         enddo
      else
         NumRotations = 1
         NumProperRotations = 1
         rot3d = ZERO
         rot3d(1,1,1) = ONE
         rot3d(2,2,1) = ONE
         rot3d(3,3,1) = ONE
         dj = CZERO
         do i = 1, kkrsz
            dj(i,i,1) = CONE
         enddo
         djc = dj
      endif
!
      end subroutine computeRotationMatrix
!     ================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function isProperRotation(i) result(y)
!     ================================================================
      use ErrorHandlerModule, only : ErrorHandler
      implicit none
!
      integer (kind=intKind), intent(in) :: i
!
      logical :: y
!
      if (i < 1 .or. i > NumRotations) then
         call ErrorHandler('isProperRotation','invalid rotation index',i)
      else if (i <= NumProperRotations) then
         y = .true.
      else
         y = .false.
      endif
!
      end function isProperRotation
!     ================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function getNumIBZRotations() result(n)
!     ================================================================
      implicit none
!
      integer (kind=intKind) :: n
!
      n = NumRotations
!
      end function getNumIBZRotations
!     ================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function getNumProperIBZRotations() result(n)
!     ================================================================
      implicit none
!
      integer (kind=intKind) :: n
!
      n = NumProperRotations
!
      end function getNumProperIBZRotations
!     ================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function getIBZRotationMatrix(c,i) result(r)
!     ================================================================
      use ErrorHandlerModule, only : ErrorHandler
      implicit none
!
      character (len=1), intent(in) :: c
!
      integer (kind=intKind), intent(in) :: i
!
      complex (kind=CmplxKind), pointer :: r(:,:)
!
      if (i < 1 .or. i > NumRotations) then
         call Errorhandler('getIBZRotationMatrix','Invalid rotation index',i)
      endif
!
      if (c == 'N' .or. c == 'n') then
         r => dj(:,:,i)
      else if (c == 'c' .or. c == 'C') then
         r => djc(:,:,i)
      else
         call Errorhandler('getIBZRotationMatrix','Invalid matrix kind',c)
      endif
!
      end function getIBZRotationMatrix
!     ================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function getIBZRotationMatrix3D(i) result(r)
!     ================================================================
      use ErrorHandlerModule, only : ErrorHandler
      implicit none
!
      integer (kind=intKind), intent(in) :: i
!
      real (kind=RealKind) :: r(3,3)
!
      if (i < 1 .or. i > NumRotations) then
         call Errorhandler('getIBZRotationMatrix3D','Invalid rotation index',i)
      endif
!
      r = rot3d(:,:,i)
!
      end function getIBZRotationMatrix3D
!     ================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine printIBZRotationMatrix(Rot3D_Only)
!     ================================================================
      use WriteMatrixModule, only : writeMatrix
!
      implicit none
!
      logical, intent(in), optional :: Rot3D_Only
      logical :: print_rot3d_only
!
      integer (kind=IntKind) :: irot, i, j
!
      print_rot3d_Only = .false.
      if (present(Rot3D_Only)) then
         if (Rot3D_Only) then
            print_rot3d_only = .true.
         endif
      endif
!
      write(6,'(/,80(''-''))')
      write(6,'(/,23x,a)')'**************************************'
      write(6,'( 23x,a )')'* Output from printIBZRotationMatrix *'
      write(6,'(23x,a,/)')'**************************************'
!
      write(6,'('' The lattice system is '',a)') lattice_name
      do irot=1,NumProperRotations
!        ===============================================================
         write(6,'(/,'' Proper rotation no.'',i4)') irot
         write(6,'('' Euler Angles '',3f10.2)')                         &
                   euler(1,irot)/PI,euler(2,irot)/PI,euler(3,irot)/PI
         write(6,'('' Geometrical Proper Rotation'')')
         do i=1,3
            write(6,'(3f15.6)') (rot3d(i,j,irot),j=1,3)
         enddo
         if (.not.print_rot3d_only) then
            write(6,'('' Angular Mom. Proper Rotation'')')
!           ------------------------------------------------------------
            call writeMatrix(' dj(m,m'') ',dj(1:kkrsz,1:kkrsz,irot),kkrsz,kkrsz,TEN2m10)
!           ------------------------------------------------------------
         endif
         write(6,'(''        Rotation no.'',i4)') irot+NumProperRotations
         write(6,'('' Geometrical Unproper Rotation'')')
         do i=1,3
            write(6,'(3f15.6)') (rot3d(i,j,irot+NumProperRotations),j=1,3)
         enddo
         if (.not.print_rot3d_only) then
            write(6,'('' Angular Mom. Unproper Rotation'')')
!           ------------------------------------------------------------
            call writeMatrix(' dj(m,m'') ',dj(1:kkrsz,1:kkrsz,irot+NumProperRotations), &
                             kkrsz,kkrsz,TEN2m10)
!           ------------------------------------------------------------
         endif
      enddo
!
      end subroutine printIBZRotationMatrix
!     ==================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getrot(latsys,intmom,nrot,nprop)
!     ==================================================================
!
!     This routine drives the calculation of the Complex Spherical
!     harmonics calculation of rotation matrices........................
!     It works for integer or odd half integer values of the Angular
!     momentum..........................................................
!     bg & eb .......................................March 1995.........
!
!     //////////////////////////////////////////////////////////////////
!     Calls:
!     spacerot: given latsys produces the set of the Euler Angles
!               corresponding to the rotation group and the
!               geometrical rotations of the group
!
!     djmatrx:  for a given j (integer or odd half integer) constucts
!               the j-th block of the proper rotation matrix
!               corresponding to 3 given Euler angles
!     //////////////////////////////////////////////////////////////////
!
!     INPUT:    latsys: Lattice Type (to select the symmetry ops.)
!               intmom:   if > 0 only integer j=l calculations
!
!     output:   rot(3,3,nrot) = geometrical rotations group
!               dj(kkrsz,kkrsz,nrot) = <j'm'|Rot|jm>
!
!     ******************************************************************
      use IntegerFactorsModule, only : m1m
!
      implicit none
!
      character (len=*), intent(in) :: latsys
!
      integer (kind=IntKind), intent(in) :: intmom
      integer (kind=IntKind), intent(out) :: nrot
      integer (kind=IntKind), intent(out) :: nprop
      integer (kind=IntKind) :: irot, jj, m, mm, mmp
      integer (kind=IntKind) :: jstart
      integer (kind=IntKind) :: jtimes
      integer (kind=IntKind) :: marray
      integer (kind=IntKind) :: kkrindex
      integer (kind=IntKind) :: kkr1, kkr2, kkrw1, kkrw2
      integer (kind=IntKind) :: l, ss, i, j, k
!
      complex (kind=CmplxKind) :: djm(kkrsz,kkrsz) ! Working space
!
!     ******************************************************************
!
!     Build the Reflection matrix to set detR=-1 rotations..............
!
      if(intmom.ne.0) then
         do l=0,lmax
            do m=l,-l,-1
               kkrindex=l*(l+1)+m+1
               reflex(kkrindex,kkrindex)=m1m(l)
            enddo
         enddo
      else
!        ---------------------------------------------------------------
         call seminv()
!        ---------------------------------------------------------------
      endif
!
!     get Euler Angles for the given lattice system.....................
!     ------------------------------------------------------------------
      call spacerot(latsys,nprop,nrot)
!     ------------------------------------------------------------------
!
!     make Angular Momentum rotation matrices (proper rotations)........
      do irot=1,nprop
!
         jstart=0
         do l=0,lmax
!
!           ************************************************************
!           Owing to the fact that dj depends only on j, not on l, the
!           j-th block evaluated for the couple (l,j=l+1/2) is written
!           also onto the block corresponding to the couple
!           (l+1,j=l-1/2) (if l < lmax).................................
!           ************************************************************
!
!           determine the 2*j value and the no. of possible values for
!           a given l...................................................
            if(intmom.ne.0) then
               jj=2*l
               jtimes=1
               marray=-1
            else
               jj=2*l+1
               marray=1
               if(l.eq.lmax) then
                  jtimes=1
               else
                  jtimes=2
               endif
            endif
!
!           evaluate the j-th block of the rotation matrix..............
!           ------------------------------------------------------------
            call djmatrx(euler(1,irot),euler(2,irot),euler(3,irot),jj,djm)
!           ------------------------------------------------------------
            do ss=1,jtimes
               do mm=jj,-jj,-2
                  kkr1=(jj-marray*mm)/2+1
                  kkrw1=(jj-mm)/2+1
                  do mmp=jj,-jj,-2
                     kkr2=(jj-marray*mmp)/2+1
                     kkrw2=(jj-mmp)/2+1
                     dj(kkr1+jstart,kkr2+jstart,irot)=djm(kkrw1,kkrw2)
                  enddo
               enddo
               jstart=jstart+jj+1
            enddo
!
         enddo
!
!        the unproper rotation matrices are evaluated as the product of
!        the proper rotation matrices and the reflection...............
!        ---------------------------------------------------------------
         do j=1,kkrsz
            do i=1,kkrsz
               dj(i,j,irot+nprop) = CZERO
            enddo
            do k=1,kkrsz
               do i=1,kkrsz
                  dj(i,j,irot+nprop)=dj(i,j,irot+nprop)+dj(i,k,irot)*reflex(k,j)
               enddo
            enddo
         enddo
!        ---------------------------------------------------------------
!
      enddo
!
      end subroutine getrot
!     ==================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine seminv()
!     ==================================================================
!
!     subroutine to build the reflection matrix for odd half integer
!     angular momenta...................................................
!
!     calls: vcoup is a function producing the general vector coupling
!            coefficients...............................................
!
!
      implicit none
!
      integer (kind=IntKind) :: nrkkrsz
      integer (kind=IntKind) :: knd
      integer (kind=IntKind) :: jj
      integer (kind=IntKind) :: mm
      integer (kind=IntKind) :: ss
      integer (kind=IntKind) :: l
      integer (kind=IntKind) :: ll
      integer (kind=IntKind) :: jstart
      integer (kind=IntKind) :: index
      integer (kind=IntKind) :: jndex
      integer (kind=IntKind) :: kndp
      integer (kind=IntKind) :: jjp
      integer (kind=IntKind) :: mmp
      integer (kind=IntKind) :: ssp
      integer (kind=IntKind) :: lp
      integer (kind=IntKind) :: llp
      integer (kind=IntKind) :: jpstart
      integer (kind=IntKind) :: indexp
      integer (kind=IntKind) :: jndexp
!
      complex (kind=CmplxKind) :: wkmat(kkrsz,kkrsz)
!
      wkmat = CZERO
!
!     set up the inversion matrix in lms representation.................
      nrkkrsz=kkrsz/2
      do l=0,lmax
         do mm=l,-l,-1
            index=l*(l+1)-mm+1
            wkmat(index,index)=sqrtm1**(2*l)
            wkmat(index+nrkkrsz,index+nrkkrsz)=sqrtm1**(2*l)
         enddo
      enddo
!
      jstart=0
      do knd=1,2*lmax+1
         l=knd/2
         ll=2*l
         if(mod(knd,2).ne.0) then
            jj=ll+1
         else
            jj=ll-1
         endif
         jpstart=0
         do kndp=1,2*lmax+1
            lp=kndp/2
            llp=2*lp
            if(mod(kndp,2).ne.0) then
               jjp=llp+1
            else
               jjp=llp-1
            endif
!
            do mm=jj,-jj,-2
               index=(jj-mm)/2+1+jstart
               do mmp=jjp,-jjp,-2
                  indexp=(jjp-mmp)/2+1+jpstart
                  do ss=1,-1,-2
                     jndex=l*(l+1)-(mm-ss)/2+1-(ss-1)*nrkkrsz/2
                     do ssp=1,-1,-2
                        jndexp=lp*(lp+1)-(mmp-ssp)/2+1-(ssp-1)*nrkkrsz/2
                        if(abs(mm-ss).le.ll .and. abs(mmp-ssp).le.llp) then
                           reflex(index,indexp)=reflex(index,indexp)+               &
                               vcoup(ll,1,jj,mm-ss,ss,mm)*wkmat(jndex,jndexp)* &
                               vcoup(llp,1,jjp,mmp-ssp,ssp,mmp)
                        endif
                     enddo
                  enddo
               enddo
            enddo
            jpstart=jpstart+jjp+1
         enddo
         jstart=jstart+jj+1
      enddo
!
      end subroutine seminv
!     ==================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spacerot(latsys,nprop,nrot)
!     ==================================================================
!
!     given the lattice system LATSYS this routine produces the Euler
!     Angles corresponding to the proper rotation matrices and the full
!     set of proper and improper geometrical rotation matrices..........
!
!     NOTE: convention on Euler Angles are taken following Edmonds
!     fig. 1.1 rather them Messiah App. C fig. 3, i.e.:
!
!     __________________________________________________________________
!     //////////////////////////////////////////////////////////////////
!
      use KindParamModule, only : IntKind, RealKind, CmplxKind
      use MathParamModule, only : ZERO, ONE, THIRD, HALF, TWO, FOUR, FIVE, PI
      use ErrorHandlerModule, only : ErrorHandler
!
      implicit none
!
      character (len=*), intent(in) :: latsys
!
      integer (kind=IntKind), intent(out) :: nprop
      integer (kind=IntKind), intent(out) :: nrot
      integer (kind=IntKind) :: irot
      integer (kind=IntKind) :: i,j
!
      real (kind=RealKind) :: alpha
      real (kind=RealKind) :: beta
      real (kind=RealKind) :: gamma
!
      interface
         function nocaseCompare(s1,s2) result(t)
            character (len=*), intent(in) :: s1
            character (len=*), intent(in) :: s2
            logical :: t
         end function nocaseCompare
      end interface
!
!     first zeroout all possible Euler Angles...........................
      euler = ZERO
!
!     define the non zero Euler Angles depending on the Lattice system.
!     Euler Angles are defined in units of pi...........................
!     ==================================================================
      if (nocaseCompare(latsys,'cubic')) then
!
!        Cubic System (Oh) 24 proper rotations...........................
         nprop=24
         euler(1, 2)=half
         euler(1, 3)=one
         euler(1, 4)= ONE+HALF ! threehalves
         euler(2, 5)=half
         euler(2, 6)=one
         euler(2, 7)= ONE+HALF ! threehalves
         euler(1, 8)=half
         euler(1, 9)=half
         euler(1,10)=half
         euler(1,11)=half
         euler(2, 8)=half
         euler(2, 9)=half
         euler(2,10)=half
         euler(2,11)=half
         euler(3, 9)=half
         euler(3,10)=one
         euler(3,11)= ONE+HALF ! threehalves
         euler(1,12)=one
         euler(1,13)=one
         euler(2,12)=half
         euler(2,13)=half
         euler(3,13)=half
         euler(1,14)= ONE+HALF ! threehalves
         euler(1,15)= ONE+HALF ! threehalves
         euler(1,16)= ONE+HALF ! threehalves
         euler(1,17)= ONE+HALF ! threehalves
         euler(2,14)=half
         euler(2,15)=half
         euler(2,16)=half
         euler(2,17)=half
         euler(3,15)=half
         euler(3,16)=one
         euler(3,17)= ONE+HALF ! threehalves
         euler(2,18)=half
         euler(2,19)=half
         euler(2,20)=half
         euler(3,18)=half
         euler(3,19)=one
         euler(3,20)= ONE+HALF ! threehalves
         euler(2,21)=one
         euler(2,22)=one
         euler(2,23)=one
         euler(3,21)=half
         euler(3,22)=one
         euler(3,23)= ONE+HALF ! threehalves
         euler(2,24)= ONE+HALF ! threehalves
         euler(3,24)=half
!
!     ==================================================================
      else if (nocaseCompare(latsys,'tetragonal')) then
!
!        Tetragonal System (D4h) 8 proper rotations.....................
         nprop=8
         euler(1,2)=half
         euler(1,3)=one
         euler(1,4)= ONE+HALF ! threehalves
         euler(2,5)=one
         euler(2,6)=one
         euler(2,7)=one
         euler(2,8)=one
         euler(3,6)=half
         euler(3,7)=one
         euler(3,8)= ONE+HALF ! threehalves
!
!     ==================================================================
      else if (nocaseCompare(latsys,'orthorhombic')) then
!
!        Orthorhombic System (D2h) 4 proper rotations...................
         nprop=4
         euler(1,2)=one
         euler(2,3)=one
         euler(2,4)=one
         euler(3,4)=one
!
!     ==================================================================
      else if (nocaseCompare(latsys,'monoclinic')) then
!
!        Monoclinic System (C2h) 2 proper rotations.....................
         nprop=2
         euler(1,2)=one
!
!     ==================================================================
      else if (nocaseCompare(latsys,'triclinic')) then
!
!        Triclinic System (C1) 1 proper rotations.......................
         nprop=1 
!
!     ==================================================================
      else if (nocaseCompare(latsys,'hexagonal')) then
!
!        Hexagonal System (D6h) 12 proper rotations......................
         nprop=12
         euler(1,2)=third
         euler(1,3)=TWO*THIRD ! twothird
         euler(1,4)=one
         euler(1,5)=FOUR*THIRD ! fourthird
         euler(1,6)=FIVE*THIRD ! fivethird
         euler(2,7)=one
         euler(2,8)=one
         euler(2,9)=one
         euler(2,10)=one
         euler(2,11)=one
         euler(2,12)=one
         euler(3,8)=third
         euler(3,9)=TWO*THIRD ! twothird
         euler(3,10)=one
         euler(3,11)=FOUR*THIRD ! fourthird
         euler(3,12)=FIVE*THIRD ! fivethird
!
!     ==================================================================
      else if (nocaseCompare(latsys,'trigonal')) then
!
!        Trigonal System (D3h) 6 proper rotations.......................
         nprop=6 
         euler(1,2)=TWO*THIRD ! twothird
         euler(1,3)=FOUR*THIRD ! fourthird
         euler(2,4)=one
         euler(2,5)=one
         euler(2,6)=one
         euler(3,5)=TWO*THIRD ! twothird
         euler(3,6)=FOUR*THIRD ! fourthird
!
!     ==================================================================
      else                                 
         call ErrorHandler('SPACEROT','Lattice system unknown',latsys)
      endif
!
!     all the above lattice systems do have the Inversion Symmetry......
      nrot=nprop*2
!
!     Change units of Euler Angles to radians...........................
      do irot=1,nprop
         do i=1,3
            euler(i,irot)=euler(i,irot)*PI
         enddo
      enddo
!
      do irot=1,nprop
!
!        make proper rotation matrices (Messiah, 2nd vol., App. C)......
         alpha=euler(1,irot)
         beta=euler(2,irot)
         gamma=euler(3,irot)
!
         rot3d(1,1,irot)=cos(gamma)*cos(beta)*cos(alpha)-sin(gamma)*sin(alpha)
         rot3d(1,2,irot)=-sin(gamma)*cos(beta)*cos(alpha)-cos(gamma)*sin(alpha)
         rot3d(1,3,irot)=sin(beta)*cos(alpha)
         rot3d(2,1,irot)=cos(gamma)*cos(beta)*sin(alpha)+sin(gamma)*cos(alpha)
         rot3d(2,2,irot)=-sin(gamma)*cos(beta)*sin(alpha)+cos(gamma)*cos(alpha)
         rot3d(2,3,irot)=sin(beta)*sin(alpha)
         rot3d(3,1,irot)=-cos(gamma)*sin(beta)
         rot3d(3,2,irot)=sin(gamma)*sin(beta)
         rot3d(3,3,irot)=cos(beta)
!
!        store the improper rotation matrices...........................
         do i=1,3
            do j=1,3
               rot3d(i,j,nprop+irot)=-rot3d(i,j,irot)
            enddo
         enddo
!
      enddo
!
      end subroutine spacerot
!     ==================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine djmatrx(alpha,beta,gamma,jj,djm)
!     ==================================================================
!
!     routine to set the Jth block of the rotation matrix corresponding
!     to the Euler angles alpha,beta,gamma..............................
!
!     bg... Feb 1995 ............ and Edmonds' book Chapter 4...........
!
!      INPUT:          alpha,beta gamma= Euler Angles in radiants
!                      jj=twice j
!                      j= half integer or integer angular momentum q.n.
!                      fact= factorials array
!
!      OUTPUT:         djm = jth block (size=2j+1) of the Rotation Matrix
!     ******************************************************************
      use IntegerFactorsModule, only : m1m
!
      implicit none
!
      integer (kind=IntKind), intent(in) :: jj
!
      integer (kind=IntKind) :: mm
      integer (kind=IntKind) :: mmp
      integer (kind=IntKind) :: mindex
      integer (kind=IntKind) :: mpindex
      integer (kind=IntKind) :: mmindex
      integer (kind=IntKind) :: mmpindex
      integer (kind=IntKind) :: isign
      integer (kind=IntKind) :: isign1
      integer (kind=IntKind) :: isign2
!
      real (kind=RealKind), intent(in) :: alpha
      real (kind=RealKind), intent(in) :: beta
      real (kind=RealKind), intent(in) :: gamma
!
      real (kind=RealKind) :: j
      real (kind=RealKind) :: m
      real (kind=RealKind) :: mp
      real (kind=RealKind) :: csi
      real (kind=RealKind) :: eta
      real (kind=RealKind) :: smalld
!
      complex (kind=CmplxKind), intent(out) :: djm(kkrsz,jj+1)
!
      j=jj*half
      djm = CZERO
!
      csi=cos(beta*half)
      eta=sin(beta*half)
!
!     according to Messiah's book convention, array is filled on
!     decreasing m from m=j to m=-j.....................................
!
!     Symmetry of djm and smalld is exploited, and the calculation of
!     smalld is limited only to m ge 0 and ge mp........................
      do mm=jj,0,-2
         m=mm*half
         mindex=(jj-mm)/2+1
         mmindex=(jj+mm)/2+1
         do mmp=mm,0,-2
            mp=mmp*half
            isign=m1m((mm-mmp)/2)
            isign1=m1m((mmp-jj)/2)
            isign2=m1m((mm+mmp)/2)
            mpindex=(jj-mmp)/2+1
            mmpindex=(jj+mmp)/2+1
!           ------------------------------------------------------------
            call edmonds(smalld,csi,eta,jj,mm,mmp)
!           ------------------------------------------------------------
!
            djm(mindex,mpindex)=exp(sqrtm1*m*alpha)*smalld*exp(sqrtm1*mp*gamma)
!
!           store the negative indeces elements.(Edmonds, eq. (4.2.7))..
            djm(mmindex,mmpindex)=conjg(djm(mindex,mpindex))*isign
!
            if(mm.ne.mp) then
!              store the transposes (Edmonds, eq. (4.2.6))..............
               djm(mpindex,mindex)=exp(sqrtm1*mp*alpha)*smalld*isign*exp(sqrtm1*m*gamma)
               djm(mmpindex,mmindex)=conjg(djm(mpindex,mindex))*isign
            endif
!
            if(mm.ne.0 .and. mmp.ne.0) then
!              change the sign of mp (Edmonds, eq. 4.2.3)...............
!              ---------------------------------------------------------
               call edmonds(smalld,-eta,csi,jj,mm,mmp)
!              ---------------------------------------------------------
!
               djm(mindex,mmpindex)=exp(sqrtm1*m*alpha)*smalld*isign1*exp(-sqrtm1*mp*gamma)
!
!              Edmonds, eq. (4.2.7).....................................
               djm(mmindex,mpindex)=conjg(djm(mindex,mmpindex))*isign2
               if(mm.ne.mmp) then
!                 Edmonds, eq. (4.2.6) and next.........................
                  djm(mmpindex,mindex)=exp(-sqrtm1*mp*alpha)*smalld*isign1*isign2*exp(sqrtm1*m*gamma)
                  djm(mpindex,mmindex)=conjg(djm(mmpindex,mindex))*isign2
               endif
            endif
!
         enddo
      enddo
!
      end subroutine djmatrx
!     ==================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine edmonds(smalld,csi,eta,j,m,mp)
!     ==================================================================
!
!     Edmonds's Formula  (4.1.15) for rotations.........................
!
      use IntegerFactorsModule, only : m1m
!
      implicit none
!
      integer (kind=IntKind), intent(in) :: j
      integer (kind=IntKind), intent(in) :: m
      integer (kind=IntKind), intent(in) :: mp
!
      integer (kind=IntKind) :: maxs
      integer (kind=IntKind) :: s
      integer (kind=IntKind) :: jpm
      integer (kind=IntKind) :: jmm
      integer (kind=IntKind) :: jpmp
      integer (kind=IntKind) :: jmmp
      integer (kind=IntKind) :: mpmp
!
      real (kind=RealKind), intent(out) :: smalld
      real (kind=RealKind), intent(in) :: csi
      real (kind=RealKind), intent(in) :: eta
!
      real (kind=RealKind) :: betafac
      real (kind=RealKind) :: sum
      real (kind=RealKind) :: factor
!
      jpm=(j+m)/2
      jmm=(j-m)/2
      jpmp=(j+mp)/2
      jmmp=(j-mp)/2
      mpmp=(m+mp)/2
!
      factor=sqrt( fact(jpmp)*fact(jmmp)/(fact(jpm)*fact(jmm)) )
!
      sum=zero
      maxs=min(jmmp,jmm)
      do s=0,maxs
         if(s+s+mpmp.eq.0) then
            betafac=one
         else
            betafac=csi**(s+s+mpmp)
         endif
         if(j-mpmp-s-s.ne.0) then
            betafac=betafac*eta**(j-mpmp-s-s)
         endif
         sum=sum+fact(jpm)*fact(jmm)/(fact(mpmp+s)*fact(jmm-s)*fact(jmmp-s)*fact(s) )*betafac*m1m(jmmp-s)
      enddo
      smalld=factor*sum
!
      end subroutine edmonds
!     ==================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function vcoup(j1,j2,j,m1,m2,m) result(vc)
!     ==================================================================
!
!     routine to evaluate Vector Coupling Coefficients
!     make use of Edmonds's equation (3.6.10) p. 44.....................
!
!     INPUT: j1,j2,j =  2* total angular momentum quantum numbers
!            m1,m2,m =  2* azimuthal angular momentum quantum umbers
!            fact= factorials array
!
!     OUTPUT: vcoup = <j1,m1,j2,m2|j1,j2,j,m>, Edmonds eq. 3.6.10
!             the same as <j1,m1,j2,m2|j,m>, Messiah's, App. C eq. (11)
!
      use KindParamModule, only : IntKind, RealKind
      use MathParamModule, only : ZERO
      use IntegerFactorsModule, only : m1m
!
      implicit   none
!
      integer (kind=IntKind), intent(in) :: j1
      integer (kind=IntKind), intent(in) :: j2
      integer (kind=IntKind), intent(in) :: j
      integer (kind=IntKind), intent(in) :: m1
      integer (kind=IntKind), intent(in) :: m2
      integer (kind=IntKind), intent(in) :: m
!
      real (kind=RealKind) :: vc
!
      integer (kind=IntKind) :: maxs
      integer (kind=IntKind) :: s
!
      real (kind=RealKind) :: sum
      real (kind=RealKind) :: factor
!
!     Consistency tests.................................................
      if(iabs(m1).gt.j1 .or. iabs(m2).gt.j2 .or. iabs(m).gt.j) then
         write(6,'('' VCOUP:: Error condition'',/,'' (j1,m1)='',2i5,    &
     &             '' (j2,m2)='',2i5,'' (j,m)='',2i5)') j1,m1,j2,m2,j,m
         write(6,'('' VCOUP:: 3-j symbol set zero'')')
         vc=zero
         return
      endif
!
!     "Triangular inequality"...........................................
      if(abs(j1-j2).gt.j .or. abs(j1+j2).lt.j) then
         vc=zero
      else
!        "Selection rule"...............................................
         if(m1+m2.ne.m) then
            vc=zero
         else
!           Calculation of 3-j symbol...................................
            factor= sqrt( real(j+1)*fact((j1+j2-j)/2)*fact((j1-m1)/2)*   &
                                    fact((j2-m2)/2)*fact((j+m)/2)*       &
                                    fact((j-m)/2)/(fact((j1+j2+j)/2+1)*  &
                                    fact((j1-j2+j)/2)*fact((j2-j1+j)/2)* &
                                    fact((j1+m1)/2)*fact((j2+m2)/2)) )
!
            maxs=max( (j1-m1)/2, (j2-m2)/2, (j1+m1)/2, (j2+m2)/2, (j-m)/2, &
                      (j+m)/2, (j1+j2-j)/2, (j2+j-j1)/2, (j+j1-j2)/2 ) + 1
            sum=zero
            do s=0,maxs
               if( ((j1+m1)/2+s).ge.0 .and. ((j2+j-m1)/2-s).ge.0         &
                                      .and. ((j1-m1)/2-s).ge.0 .and.     &
                   ((j-m)/2-s).ge.0 .and. ((j2-j+m1)/2+s).ge.0) then
                  sum=sum + m1m((j1-m1)/2+s)*fact((j1+m1)/2+s)*fact((j2+j-m1)/2-s)/ &
                      (fact(s)*fact((j1-m1)/2-s)*fact((j-m)/2-s)*fact((j2-j+m1)/2+s))
               endif
            enddo
            vc=factor*sum
         endif
      endif
!
      end function vcoup
!     ==================================================================
end module IBZRotationModule
