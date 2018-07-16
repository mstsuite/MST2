module ExchCorrFunctionalModule
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, TWO, THREE, THIRD, PI4, TEN2m6, &
                               TEN2m8, TEN2m10
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler,        &
                                  StopHandler
!
public :: initExchCorrFunctional,         &
          endExchCorrFunctional,          &
          getExchCorrPot,                 &
          getExchCorrEnDen,               &
          calExchangeCorrelation
!
   interface getExchCorrPot
      module procedure getExchCorrPot_s, getExchCorrPot_v
   end interface
!
   interface getExchCorrEnDen
      module procedure getExchCorrEnDen_s, getExchCorrEnDen_v
   end interface
!
   interface calExchangeCorrelation
      module procedure calExchangeCorrelation_s, calExchangeCorrelation_v
   end interface
!
private 
   logical :: Initialized = .false.
!
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: typeFunctional
!
   real (kind=RealKind) ::    Vexc_s
   real (kind=RealKind) ::    Eexc_s
   real (kind=RealKind), allocatable, target :: Vexc_v(:)
   real (kind=RealKind), allocatable, target :: Eexc_v(:)
!
   integer :: nV, nE
!
contains
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initExchCorrFunctional( pola, excorr )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: excorr
   integer (kind=IntKind), intent(in) :: pola
!
   n_spin_pola = pola
   typeFunctional = excorr
!
   if ( typeFunctional < 0 .or. typeFunctional >1 ) then
      call ErrorHandler( "initExchCorrFunctional",                    &
                         "Invalid LDA functional", typeFunctional )
   endif
!
   nV = 0
   nE = 0
   Initialized = .true.
!
   end subroutine initExchCorrFunctional
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExchCorrEnDen_s( rho_den, mag_den, is )       result(E)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is 
!
   real (kind=RealKind), intent(in) :: rho_den
   real (kind=RealKind), intent(in) :: mag_den
!
   real (kind=RealKind) :: E
!
   if (.not.Initialized) then
      call ErrorHandler('getExchCorrEnDen_s',                         &
                        'ExchCorrFunctional is not initialized')
   else if ( is<0 .or. is>2 ) then
      call ErrorHandler("getExchCorrEnDen_s","Wrong spin index",is)
   endif
!
!  ------------------------------------------------------------------
   call calExchCorr_s( rho_den, mag_den, is )
!  ------------------------------------------------------------------
!
   E = Eexc_s
!
   end function getExchCorrEnDen_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExchCorrEnDen_v(n_Rpts, rho_den, mag_den, is)  result(pE)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n_Rpts
   integer (kind=IntKind), intent(in) :: is 
!
   real (kind=RealKind), intent(in) :: rho_den(n_Rpts)
   real (kind=RealKind), intent(in) :: mag_den(n_Rpts)
!
   real (kind=RealKind), pointer :: pE(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getExchCorrEnDen_v',                         &
                        'ExchCorrFunctional is not initialized')
   else if ( is < 0 .or. is > 2 ) then
      call ErrorHandler("getExchCorrEnDen_v","Wrong spin index",is)
   else if ( n_Rpts <= 0 ) then
      call ErrorHandler("getExchCorrEnDen_v","n_Rpts out of range",n_Rpts)
   endif
!
   if ( allocated(Eexc_v) .and. nE < n_Rpts ) then
      deallocate( Eexc_v )
   endif
   if ( .not.allocated(Eexc_v) ) then
      allocate( Eexc_v(1:n_Rpts) )
      nE = n_Rpts
   endif
!
   if ( allocated(Vexc_v) .and. nV < n_Rpts ) then
      deallocate( Vexc_v )
   endif
   if ( .not.allocated(Vexc_v) ) then
      allocate( Vexc_v(1:n_Rpts) )
      nV = n_Rpts
   endif
!
!  ------------------------------------------------------------------
   call calExchCorr_v( n_Rpts, rho_den, mag_den, is )
!  ------------------------------------------------------------------
!
   pE => Eexc_v(1:n_Rpts)
!
   end function getExchCorrEnDen_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExchCorrPot_s( rho_den, mag_den, is )         result(V)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is 
!
   real (kind=RealKind), intent(in) :: rho_den
   real (kind=RealKind), intent(in) :: mag_den
!
   real (kind=RealKind) :: V
!
   if (.not.Initialized) then
      call ErrorHandler('getExchCorrPot_s',                         &
                        'ExchCorrFunctional is not initialized')
   else if ( is<0 .or. is>2 ) then
      call ErrorHandler("getExchCorrPot_s","Wrong spin index",is)
   endif
!
!  ------------------------------------------------------------------
   call calExchCorr_s( rho_den, mag_den, is )
!  ------------------------------------------------------------------
!
   V = Vexc_s
!
   end function getExchCorrPot_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calExchangeCorrelation_s( rho_den, mag_den, is, pot, en )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is 
!
   real (kind=RealKind), intent(in) :: rho_den
   real (kind=RealKind), intent(in) :: mag_den
!
   real (kind=RealKind), intent(out) :: pot
   real (kind=RealKind), intent(out) :: en
!
   if (.not.Initialized) then
      call ErrorHandler('calExchangeCorrelation_s',                         &
                        'ExchCorrFunctional is not initialized')
   else if ( is<0 .or. is>2 ) then
      call ErrorHandler("calExchangeCorrelation_s","Wrong spin index",is)
   endif
!
!  ------------------------------------------------------------------
   call calExchCorr_s( rho_den, mag_den, is )
!  ------------------------------------------------------------------
!
   pot = Vexc_s
   en  = Eexc_s
!
   end subroutine calExchangeCorrelation_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calExchangeCorrelation_v( n_Rpts, rho_den, mag_den, is, pot, en )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n_Rpts
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind), intent(in) :: rho_den(n_Rpts)
   real (kind=RealKind), intent(in) :: mag_den(n_Rpts)
!
   real (kind=RealKind), pointer :: pot(:)
   real (kind=RealKind), pointer :: en(:)
!
   if (.not.Initialized) then
      call ErrorHandler('calExchangeCorrelation_v',                 &
                        'ExchCorrFunctional is not initialized')
   else if ( is<0 .or. is>2 ) then
      call ErrorHandler("gcalExchangeCorrelation_v","Wrong spin index",is)
   endif
!
!  ------------------------------------------------------------------
   call calExchCorr_v( n_rpts, rho_den, mag_den, is )
!  ------------------------------------------------------------------
!
   pot => Vexc_v(1:n_Rpts)
   en  => Eexc_v(1:n_Rpts)
!
   end subroutine calExchangeCorrelation_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExchCorrPot_v(n_Rpts, rho_den, mag_den, is)  result(pV)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n_Rpts
   integer (kind=IntKind), intent(in) :: is 
!
   real (kind=RealKind), intent(in) :: rho_den(n_Rpts)
   real (kind=RealKind), intent(in) :: mag_den(n_Rpts)
!
   real (kind=RealKind), pointer :: pV(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getExchCorrPot_v',                         &
                        'ExchCorrFunctional is not initialized')
   else if ( is < 0 .or. is > 2 ) then
      call ErrorHandler("getExchCorrPot_v","Wrong spin index",is)
   else if ( n_Rpts <= 0 ) then
      call ErrorHandler("getExchCorrPot_v","n_Rpts out of range",n_Rpts)
   endif
!
   if ( allocated(Eexc_v) .and. nE < n_Rpts ) then
      deallocate( Eexc_v )
   endif
   if ( .not.allocated(Eexc_v) ) then
      allocate( Eexc_v(1:n_Rpts) )
      nE = n_Rpts
   endif
!
   if ( allocated(Vexc_v) .and. nV < n_Rpts ) then
      deallocate( Vexc_v )
   endif
   if ( .not.allocated(Vexc_v) ) then
      allocate( Vexc_v(1:n_Rpts) )
      nV = n_Rpts
   endif
!
!  ------------------------------------------------------------------
   call calExchCorr_v( n_Rpts, rho_den, mag_den, is )
!  ------------------------------------------------------------------
!
   pV => Vexc_v(1:n_Rpts)
!
   end function getExchCorrPot_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calExchCorr_s( rho_den, mag_den, is )
!  ===================================================================
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind), intent(in) :: rho_den
   real (kind=RealKind), intent(in) :: mag_den
!
   real (kind=RealKind) :: dz
   real (kind=RealKind) :: r_s
   real (kind=RealKind) :: sp
!
!  ===================================================================
!  exchange-correlation potential and density at a point
!  ===================================================================
!
   sp = THREE-TWO*is 
   if ( rho_den > ZERO ) then
      dz = mag_den/rho_den
      if (abs(dz) > ONE+TEN2m10) then
!        =============================================================
!        call ErrorHandler('calExchCorr','abs(dz) > 1',dz,.true.)
! modified here to help the vaccume case 05/29/17
         call WarningHandler('calExchCorr','abs(dz) > 1',dz,.true.)
         if ( dz > ONE ) then
            dz = ONE
         else
            dz = -ONE
         endif
!     else if ( dz > ONE ) then
!        dz = ONE
!     else if (dz < -ONE) then
!        dz = -ONE
!        =============================================================
      endif
      r_s = (THREE/(PI4*rho_den))**THIRD
!     ----------------------------------------------------------------
      call LDAfunctional( r_s, dz, sp, Eexc_s, Vexc_s )
!     ----------------------------------------------------------------
   else
      Vexc_s = ZERO
      Eexc_s = ZERO
   endif
!
   end subroutine calExchCorr_s 
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calExchCorr_v( n_Rpts, rho_den, mag_den, is )
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: n_Rpts
!
   real (kind=RealKind), intent(in) :: rho_den(n_Rpts)
   real (kind=RealKind), intent(in) :: mag_den(n_Rpts)
!
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind) :: dz
   real (kind=RealKind) :: r_s
   real (kind=RealKind) :: sp
!
!  ===================================================================
!  exchange-correlation potential and density for an array of points
!  ===================================================================
!
   sp = THREE-TWO*is 
!
   do ir = 1,n_Rpts
      if ( rho_den(ir) > ZERO ) then
         dz = mag_den(ir)/rho_den(ir)
         if (abs(dz) > ONE+TEN2m10) then
!           ==========================================================
!           call ErrorHandler('calExchCorr','abs(dz) > 1',dz,.true.)
! modified here to help the vaccume case 05/29/17
            call WarningHandler('calExchCorr','abs(dz) > 1',dz,.true.)
            if ( dz > ONE ) then
               dz = ONE
            else
               dz = -ONE
            endif
!        else if ( dz > ONE ) then
!           dz = ONE
!        else if (dz < -ONE) then
!           dz = -ONE
!           ==========================================================
         endif
         r_s = (THREE/(PI4*rho_den(ir)))**THIRD
!        -------------------------------------------------------------
         call LDAfunctional( r_s, dz, sp, Eexc_v(ir), Vexc_v(ir) )
!        -------------------------------------------------------------
      else
         Vexc_v(ir) = ZERO
         Eexc_v(ir) = ZERO
      endif
   enddo
!
   end subroutine calExchCorr_v 
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine LDAfunctional(r_s, dz, sp, exchg, vxchg)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind) ::  i
!
   real (kind=RealKind), intent(in) ::  r_s
   real (kind=RealKind), intent(in) ::  dz
   real (kind=RealKind), intent(in) ::  sp
   real (kind=RealKind), intent(out) ::  exchg
   real (kind=RealKind), intent(out) ::  vxchg
!
   real (kind=RealKind) ::  vxx(2)
   real (kind=RealKind) ::  g(3)
   real (kind=RealKind) ::  dg(3)
   real (kind=RealKind) ::  tbq(3)
   real (kind=RealKind) ::  tbxq(3)
   real (kind=RealKind) ::  bxx(3)
   real (kind=RealKind) ::  q(3)
   real (kind=RealKind) ::  bb(3)
   real (kind=RealKind) ::  cx(3)
   real (kind=RealKind) ::  fm
   real (kind=RealKind) ::  fdz
   real (kind=RealKind) ::  ex
   real (kind=RealKind) ::  exf
   real (kind=RealKind) ::  xp
   real (kind=RealKind) ::  xf
   real (kind=RealKind) ::  gp
   real (kind=RealKind) ::  gf
   real (kind=RealKind) ::  exc
   real (kind=RealKind) ::  excf
   real (kind=RealKind) ::  dedz
   real (kind=RealKind) ::  gpp
   real (kind=RealKind) ::  gfp
   real (kind=RealKind) ::  depd
   real (kind=RealKind) ::  defd
   real (kind=RealKind) ::  decd
   real (kind=RealKind) ::  bfc
   real (kind=RealKind) ::  zp1
   real (kind=RealKind) ::  zm1
   real (kind=RealKind) ::  xr
   real (kind=RealKind) ::  pex
   real (kind=RealKind) ::  xrsq
   real (kind=RealKind) ::  qi
   real (kind=RealKind) ::  txb
   real (kind=RealKind) ::  fx
   real (kind=RealKind) ::  arct
   real (kind=RealKind) ::  dxs
   real (kind=RealKind) ::  vcc
   real (kind=RealKind) ::  facc
   real (kind=RealKind) ::  ecp
   real (kind=RealKind) ::  zp3
   real (kind=RealKind) ::  zm3
   real (kind=RealKind) ::  zp3m3
   real (kind=RealKind) ::  fx1
   real (kind=RealKind) ::  z4
   real (kind=RealKind) ::  fz
   real (kind=RealKind) ::  beta
   real (kind=RealKind) ::  ec
   real (kind=RealKind) ::  f3ex
!
!  ==================================================================
!  data for von Barth-Hedin
!  ==================================================================
   real (kind=RealKind), parameter ::  ccp = 0.0450d+00
   real (kind=RealKind), parameter ::  rp = 21.0d+00
   real (kind=RealKind), parameter ::  ccf = 0.02250d+00
   real (kind=RealKind), parameter ::  rf = 52.9166820d+00
!
!  ==================================================================
!  data for Vosko-Wilks-Nusair
!  ==================================================================
!
   real (kind=RealKind), parameter ::  a(3) = (/-0.0337740d+00,      &
                                       0.06218140d+00,0.03109070d+00/)
   real (kind=RealKind), parameter ::  b(3) = (/1.131070d+00,        &
                                       3.727440d+00,7.060420d+00/)
   real (kind=RealKind), parameter ::  c(3) = (/13.00450d+00,        &
                                       12.93520d+00,18.05780d+00/)
   real (kind=RealKind), parameter ::  x0(3) = (/-0.00475840d+00,    &
                                       -0.104980d+00,-0.32500d+00/)
   real (kind=RealKind), parameter ::  cst = 1.923661050d+00
   real (kind=RealKind), parameter ::  aip = 0.916330590d+00
   real (kind=RealKind), parameter ::  fnot = 1.709920950d+00
   real (kind=RealKind), parameter ::  bip = 0.259921050d+00
   real (kind=RealKind), parameter ::  for3 = ONE + THIRD
   real (kind=RealKind), parameter ::  thrd = THIRD
!
   if ( typeFunctional == 0 ) then  
!     ================================================================
!     von Barth-Hedin  exch-corr potential
!     j. phys. c5,1629(1972)
!     ================================================================
      fm  = 2.0d+00**(4.0d+00/3.0d+00)-2.0d+00
      fdz = ((1.0d+00+dz)**(4.0d+00/3.0d+00)                         &
           +(1.0d+00-dz)**(4.0d+00/3.0d+00)-2.0d+00)/fm
      ex  = -0.916330d+00/r_s
      exf = ex*2.0d+00**0.333333330d+00
      xp  = r_s/rp
      xf  = r_s/rf
      gp = (1.0d+00+xp**3)*log(1.0d+00+1.0d+00/xp)                   &
          -xp*xp +xp/2.0d+00 - 0.333333330d+00
      gf = (1.0d+00+xf**3)*log(1.0d+00+1.0d+00/xf)                   &
          -xf*xf +xf/2.0d+00 - 0.333333330d+00
      exc  = ex-ccp*gp
      excf = exf-ccf*gf
      dedz = (4.0d+00/3.0d+00)*(excf-exc)                            &
            *((1.0d+00+dz)**(1.0d+00/3.0d+00)                        &
            -(1.0d+00-dz)**(1.0d+00/3.0d+00))/fm
      gpp = 3.0d+00*xp*xp*log(1.0d+00+1.0d+00/xp)-1.0d+00/xp         &
           +1.50d+00-3.0d+00*xp
      gfp = 3.0d+00*xf*xf*log(1.0d+00+1.0d+00/xf)-1.0d+00/xf         &
           +1.50d+00-3.0d+00*xf
      depd = -ex/r_s-ccp/rp*gpp
      defd = -exf/r_s-ccf/rf*gfp
      decd = depd+(defd-depd)*fdz
!     ================================================================
!     exchange-correlation energy
!     ================================================================
      exchg = exc + (excf-exc)*fdz
!     ================================================================
!     exchange-correlation potential
!     ================================================================
      vxchg = exc+(excf-exc)*fdz-r_s*decd/3.0d+00                     &
              +sp*(1.0d+00-sp*dz)*dedz
!
   else if ( typeFunctional == 1 ) then
!     ================================================================
!     Vosko-Wilks-Nusair exch-corr potential
!     From G.S. Painter : Phys. Rev. B24 4264,1981
!     ================================================================
!
!     ================================================================
!     generate constant coefficients for the parameterization (v-w-n)
!     ================================================================
      do i = 1,3
         cx(i)  = x0(i)**2 + b(i)*x0(i) + c(i)
         bfc    = 4.0d+00*c(i) - b(i)**2.0d+00
         q(i)   = sqrt(bfc)
         bxx(i) = b(i)*x0(i)/cx(i)
         tbq(i) = 2.0d+00*b(i)/q(i)
         tbxq(i)= tbq(i) + 4.0d+00*x0(i)/q(i)
         bb(i)  = 4.0d+00*b(i)*( 1 - x0(i)*(b(i) + 2.0d+00*x0(i))/cx(i) )
      enddo
!
      zp1  = 1.0d+00 + dz
      zm1  = 1.0d+00 - dz
      xr   = sqrt(r_s)
      pex  = -aip/r_s
      xrsq = r_s
!
!     ===============================================================
!     generate g(i)=alpha,epsilon fct.s
!     and their derivatives dg(i)
!     1=alpha(spin stiffness)  2=ecp  3=ecf
!     ===============================================================
      do i = 1,3
         qi   = q(i)
         txb  = 2.0d+00*xr + b(i)
         fx   = xrsq + xr*b(i) + c(i)
         arct = atan2(qi,txb)
         dxs  = (xr-x0(i))**2/fx
         g(i) = a(i)*( log(xrsq/fx) + tbq(i)*arct-bxx(i)*(log(dxs)  &
                +tbxq(i)*arct) )
         dg(i) = a(i)*( 2.0d+00/xr - txb/fx                         &
                -bxx(i)*(2.0d+00/(xr-x0(i))-txb/fx)                 &
                -bb(i)/(qi**2 + txb**2) )
      enddo
!
      ecp = g(2)
      zp3 = zp1**thrd
      zm3 = zm1**thrd
      zp3m3 = zp3-zm3
!     ===============================================================
!     part of last term in vx   eq(13)
!     ===============================================================
      fx1  = .50d+00*for3*pex*zp3m3
      z4   = dz**4
      fz   = cst*(zp1**for3 + zm1**for3 - 2.0d+00)
      beta = fnot*( g(3)-g(2) )/g(1) -1.0d+00
      ec   = ecp + fz*g(1)*( 1.0d+00 + z4*beta )/fnot
      ex   = pex*( 1.0d+00 + fz*bip )
      f3ex = for3*ex
!     ===============================================================
!     echange-correlation energy
!     ===============================================================
      exchg = ec + ex
!     ===============================================================
!     exchange potential
!     ===============================================================
      vxx(1) = f3ex + fx1*zm1
      vxx(2) = f3ex - fx1*zp1
!     ===============================================================
!     correlation potential
!     ===============================================================
      vcc = ec - xr*( (1.0d+00-z4*fz)*dg(2) + z4*fz*dg(3)           &
           +(1.0d+00 - z4)*fz*dg(1)/fnot )/6.0d+00
!
      facc = 4.0d+00*g(1)*( dz**3*fz*beta            &
            +( 1.0d+00 + beta*z4 )*zp3m3/(6.0d+00*bip) )/fnot
!
!     ===============================================================
!     exch-corr. potential for each spin as called in newpot
!     ===============================================================
!
      if ( sp >= 0 ) then 
         vxchg = vcc + zm1*facc + vxx(1)
      else
         vxchg = vcc - zp1*facc + vxx(2)
      endif
   endif
!
   end subroutine LDAfunctional
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endExchCorrFunctional()
!  ===================================================================
!
   implicit none
!
   if ( allocated(Vexc_v) ) then
      deallocate( Vexc_v )
   endif
!
   if ( allocated(Eexc_v) ) then
      deallocate( Eexc_v )
   endif
!
   nV = 0
   nE = 0
!
   Initialized = .false.
!
   end subroutine endExchCorrFunctional
!  ===================================================================
!
end module ExchCorrFunctionalModule
