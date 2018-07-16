!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getpotg(local_id, specie, num_atoms, num_species, nspin, n_spin_pola, &
                      efpot,evec,                                      &
                      nr,nrmax,jmt,jws,r_mesh,jmax_pot,                &
                      vr,vdif,rhotot,xvalws,                           &
                      ztotss,zcorss,numc,nc,lc,kc,ec,lst,p_numcmax,    &
                      pot_l, iform,                                    &
                      header,vunit,iprint,istop)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, TEN2m6, HALF, TWO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use MPPModule, only : sendMessage, recvMessage, waitMessage
   use MPPModule, only : setCommunicator, resetCommunicator, syncAllPEs
!
   use ParallelIOModule, only : isInputProc, getMyInputProc,           &
                                getNumInputClients, getInputClient,    &
                                getIOCommunicator, getMyPEinIOGroup, getNumPEsInIOGroup
!
   use SystemModule, only : getAtomicNumber
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use LdaCorrectionModule, only : checkLdaCorrection, insertDataPackFromInput
!
   use InterpolationModule, only : getInterpolation
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: p_numcmax
   integer (kind=IntKind), intent(in) :: n_spin_pola
!
   character (len=*), intent(in) :: istop
   character (len=*), intent(out) :: header
   character(len=5), intent(out) :: lst(p_numcmax,n_spin_pola)
!
   character (len=5) :: ctmp
   character(len=1) :: cmsgbuf((160+40*p_numcmax)*2)
!
   logical :: checking_atoms = .true.
   logical :: checking_spin = .true.
!
   integer (kind=IntKind), intent(in) :: local_id
   integer (kind=IntKind), intent(in) :: specie
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: num_species
   integer (kind=IntKind), intent(in) :: nspin
   integer (kind=IntKind), intent(in) :: jmt
   integer (kind=IntKind), intent(in) :: jws
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), intent(in) :: nrmax
   integer (kind=IntKind), intent(in) :: vunit
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind), intent(in) :: jmax_pot
   integer (kind=IntKind), intent(out) :: numc
   integer (kind=IntKind), intent(out) :: iform
   integer (kind=IntKind), intent(out) :: nc(p_numcmax,n_spin_pola)
   integer (kind=IntKind), intent(out) :: lc(p_numcmax,n_spin_pola)
   integer (kind=IntKind), intent(out) :: kc(p_numcmax,n_spin_pola)
!
   integer (kind=IntKind) :: jmtmax
   integer (kind=IntKind) :: jwsmax
   integer (kind=IntKind) :: numcmax
   integer (kind=IntKind) :: is, ns, nspin_in, n_spin_pola_in
   integer (kind=IntKind) :: present_atom
   integer (kind=IntKind) :: fp_pos
   integer (kind=IntKind) :: msg_bytes
   integer (kind=IntKind) :: itmp
   integer (kind=IntKind) :: i, ig, ia
   integer (kind=IntKind) :: jl, nr_old, ir, jmt0, nj
   integer (kind=IntKind) :: slen, pad_bytes
   integer (kind=IntKind) :: imsgbuf(20)
!
   integer (kind=IntKind) :: num_clients
   integer (kind=IntKind) :: proc_client
   integer (kind=IntKind) :: integer4_size
   integer (kind=IntKind) :: real8_size
   integer (kind=IntKind) :: dsize_ldapu, dsize_nspot, dsize_min
   integer (kind=IntKind) :: comm, MyPEinGroup, NumPEsInGroup, na_in
!
   integer (kind=IntKind), parameter :: fsize = 10000 ! > (9+2*nrmax)*n_spin_pola
!
   real (kind=RealKind), intent(in) :: r_mesh(nr)
   real (kind=RealKind), intent(in) :: ztotss
   real (kind=RealKind), intent(in) :: zcorss
   real (kind=RealKind), intent(out) :: evec(3)
   real (kind=RealKind), intent(out) :: efpot
   real (kind=RealKind), intent(out) :: vr(nr,n_spin_pola)
   real (kind=RealKind), intent(out) :: rhotot(nr,n_spin_pola)
   real (kind=RealKind), intent(out) :: xvalws(n_spin_pola)
   real (kind=RealKind), intent(out) :: vdif
   real (kind=RealKind), intent(out) :: ec(p_numcmax,n_spin_pola)
   real (kind=RealKind) :: fspace(fsize)
!  real (kind=RealKind), intent(out) :: fspace((9+2*nrmax)*n_spin_pola)
!
   real (kind=RealKind) :: alat, za
   real (kind=RealKind) :: rtmp, xst, xmt, hh, dvr
   real (kind=RealKind), allocatable :: data_ldapu(:)
   real (kind=RealKind), allocatable, target :: data_nspot(:)
   real (kind=RealKind), allocatable :: r_mesh_old(:)
!
   complex (kind=CmplxKind), target :: pot_l(nr,jmax_pot,n_spin_pola)
   complex (kind=CmplxKind), allocatable, target :: c_data_nspot(:)
   complex (kind=CmplxKind), pointer :: p_data_nspot(:)
!
!  *******************************************************************
!  cccccccccc read in the potentials for current sublattice  ccccccccc
!  *******************************************************************
!
   MyPEinGroup = getMyPEinIOGroup()
   NumPEsInGroup = getNumPEsInIOGroup()
   comm = getIOCommunicator()
   call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!  -------------------------------------------------------------------
   rhotot(1:nr,1:n_spin_pola) = ZERO
   imsgbuf(1:20) = 0
   call c_dtsize(integer4_size,real8_size)
   dsize_ldapu = 0
   dsize_nspot = 0
!  -------------------------------------------------------------------
!
!  ===================================================================
!  using gopen to perform the global openning......................
!  -------------------------------------------------------------------
   if( isInputProc() ) then
!     ================================================================
!     check the number of atoms in the potential file.................
!     ================================================================
      if (checking_atoms) then
         na_in = 0
         LOOP_do: do
            na_in = na_in + 1
            fp_pos=(na_in-1)*integer4_size*10+1
            call c_fseek(vunit,fp_pos,0)
            call c_read_integer(vunit,imsgbuf,10)
            if (imsgbuf(1) < 0) then
!              write(6,'(a)')'The potential contains modified data format.........!'
               imsgbuf(1) = -imsgbuf(1)
            endif
            if (imsgbuf(1) /= na_in) then
               na_in = na_in - 1
               exit LOOP_do
            endif
         enddo LOOP_do
         if (na_in /= num_atoms) then
            write(6,'(a,i6)')'The number atoms in the potential file = ',na_in
         endif
      else
         na_in = num_atoms
      endif
!     ================================================================
!     read in the potentials..........................................
!     ================================================================
      num_clients = getNumInputClients()
      do i = 1,num_clients
         proc_client = getInputClient(i)
         ig = getGlobalIndex(local_id,proc_client)
!        =============================================================
         if (na_in /= num_atoms) then
            present_atom = -1
            LOOP_ia: do ia = 1, na_in   ! Loop over the atoms in the potential file to find the matching atom
               fp_pos=(ia-1)*integer4_size*10+1
!              -------------------------------------------------------
               call c_fseek(vunit,fp_pos,0)
               call c_read_integer(vunit,imsgbuf,10)
               call c_string_padsize(vunit,imsgbuf(2),pad_bytes)
!              -------------------------------------------------------
               msg_bytes=imsgbuf(2)+pad_bytes+(imsgbuf(3)+3)*real8_size
               if (imsgbuf(1) > 0) then
                  fp_pos=na_in*integer4_size*10+(ia-1)*msg_bytes+1
               else
                  fp_pos=na_in*integer4_size*10+fp_pos
!                 ----------------------------------------------------
                  call c_fseek(vunit,fp_pos,0)
                  call c_read_integer(vunit,imsgbuf(11:20),10)
!                 ----------------------------------------------------
                  fp_pos=imsgbuf(12)*real8_size +                     &
                         na_in*integer4_size*20+(ia-1)*msg_bytes+1
               endif
               fp_pos = fp_pos + imsgbuf(2) + 3*real8_size
!              -------------------------------------------------------
               call c_fseek(vunit,fp_pos,0)
               call c_read_double(vunit,za,1)
!              -------------------------------------------------------
               if (int(za) /= getAtomicNumber(ig)) then
                  write(6,'(a,2i5)')'Za in file /= Za of my atom: ',int(za),getAtomicNumber(ig)
               else
                  present_atom = ia
                  exit LOOP_ia
               endif
            enddo LOOP_ia
         else
            present_atom = ig
         endif
         if (present_atom < 1) then
            call ErrorHandler('getpotg','My atom potential is not found in the potential data',ig)
         endif
!        =============================================================
         fp_pos=(present_atom-1)*integer4_size*10+1
!        -------------------------------------------------------------
         call c_fseek(vunit,fp_pos,0)
         call c_read_integer(vunit,imsgbuf,10)
!        -------------------------------------------------------------
         if (imsgbuf(3) > fsize) then
            call ErrorHandler('getpotg','imsgbuf(3) > size(fspace): ', &
                              present_atom,imsgbuf(3),fsize)
         endif
         if (imsgbuf(1) < 0) then
            fp_pos=na_in*integer4_size*10+fp_pos
            call c_fseek(vunit,fp_pos,0)
            call c_read_integer(vunit,imsgbuf(11:20),10)
            if (dsize_ldapu < imsgbuf(13)) then
               if ( allocated(data_ldapu) ) then
                  deallocate( data_ldapu )
               endif
               allocate( data_ldapu( 1:imsgbuf(13) ) )
               dsize_ldapu = imsgbuf(13)
            endif
            if (dsize_nspot < imsgbuf(14)) then
               if ( allocated(data_nspot) ) then
                  deallocate( data_nspot )
               endif
               allocate( data_nspot( 1:imsgbuf(14) ) )
               dsize_nspot = imsgbuf(14)
            endif
         endif
!        -------------------------------------------------------------
         jmtmax=imsgbuf(8)
         jwsmax=imsgbuf(9)
         numcmax=imsgbuf(10)
!        =============================================================
!        reading cmsgbuf, fspace and evec arrays from vfile........
!        =============================================================
         call c_string_padsize(vunit,imsgbuf(2),pad_bytes)
!             print *,'imsgbuf = ',imsgbuf(2),', pad_bytes = ',pad_bytes
         msg_bytes=imsgbuf(2)+pad_bytes+(imsgbuf(3)+3)*real8_size
         if (imsgbuf(1) > 0) then
            fp_pos=na_in*integer4_size*10+(present_atom-1)*msg_bytes+1
         else
            fp_pos=imsgbuf(12)*real8_size +                           &
                   na_in*integer4_size*20+(present_atom-1)*msg_bytes+1
         endif
!        -------------------------------------------------------------
         call c_fseek(vunit,fp_pos,0)
         call c_read_string(vunit,cmsgbuf,imsgbuf(2),slen)
         call c_read_double(vunit,fspace,imsgbuf(3))
         call c_read_double(vunit,evec,3)
         if (imsgbuf(13) > 0) then  
            call c_read_double(vunit,data_ldapu,imsgbuf(13))
         endif
         if (imsgbuf(14) > 0) then
            call c_read_double(vunit,data_nspot,imsgbuf(14))
         endif
!        -------------------------------------------------------------
         if (imsgbuf(4) < 3) then
            evec(1)=0.0d0
            evec(2)=0.0d0
            evec(3)=1.0d0
         endif
!        -------------------------------------------------------------
         call sendMessage(imsgbuf,20,23457,proc_client)
         call sendMessage(cmsgbuf,imsgbuf(2),23458,proc_client)
         call sendMessage(fspace,imsgbuf(3),23459,proc_client)
         call sendMessage(evec,3,23465,proc_client)
         if (imsgbuf(13) > 0) then  
!           ----------------------------------------------------------
            call sendMessage(data_ldapu,imsgbuf(13),23466,proc_client)
!           ----------------------------------------------------------
         endif
         if (imsgbuf(14) > 0) then  
!           ----------------------------------------------------------
            call sendMessage(data_nspot,imsgbuf(14),23467,proc_client)
!           ----------------------------------------------------------
         endif
      enddo
!
      present_atom = getGlobalIndex(local_id)
      fp_pos=(present_atom-1)*integer4_size*10+1
!     ----------------------------------------------------------------
      call c_fseek(vunit,fp_pos,0)
      call c_read_integer(vunit,imsgbuf,10)
      if (imsgbuf(3) > fsize) then
         call ErrorHandler('getpotg','imsgbuf(3) > size(fspace): ',   &
                           present_atom,imsgbuf(3),fsize)
      endif
      if (imsgbuf(1) < 0) then
         fp_pos=na_in*integer4_size*10+fp_pos
         call c_fseek(vunit,fp_pos,0)
         call c_read_integer(vunit,imsgbuf(11:20),10)
         if (dsize_ldapu < imsgbuf(13)) then
            if ( allocated(data_ldapu) ) then
               deallocate( data_ldapu )
            endif
            allocate( data_ldapu( 1:imsgbuf(13) ) )
            dsize_ldapu = imsgbuf(13)
         endif
         if (dsize_nspot < imsgbuf(14)) then
            if ( allocated(data_nspot) ) then
               deallocate( data_nspot )
            endif
            allocate( data_nspot( 1:imsgbuf(14) ) )
            dsize_nspot = imsgbuf(14)
         endif
      endif
!     ----------------------------------------------------------------
      jmtmax=imsgbuf(8)
      jwsmax=imsgbuf(9)
      numcmax=imsgbuf(10)
!
      if (numcmax > p_numcmax) then
         call ErrorHandler('getpotg','numcmax > p_numcmax',numcmax,p_numcmax)
      endif
!     ================================================================
!     reading cmsgbuf, fspace and evec arrays from vfile..............
!     ================================================================
      call c_string_padsize(vunit,imsgbuf(2),pad_bytes)
!          print *,'pad_bytes = ',pad_bytes
      msg_bytes=imsgbuf(2)+pad_bytes+(imsgbuf(3)+3)*real8_size
      if (imsgbuf(1) > 0) then
         fp_pos=na_in*integer4_size*10+(present_atom-1)*msg_bytes+1
      else
         fp_pos=imsgbuf(12)*real8_size +                              &
             na_in*integer4_size*20+(present_atom-1)*msg_bytes+1
      endif
!     ----------------------------------------------------------------
      call c_fseek(vunit,fp_pos,0)
      call c_read_string(vunit,cmsgbuf,imsgbuf(2),slen)
      call c_read_double(vunit,fspace,imsgbuf(3))
      call c_read_double(vunit,evec,3)
!     ----------------------------------------------------------------
      if(imsgbuf(4) < 3) then
         evec(1)=0.0d0
         evec(2)=0.0d0
         evec(3)=1.0d0
      endif
      if (imsgbuf(13) > 0) then  
         call c_read_double(vunit,data_ldapu,imsgbuf(13))
      endif
      if (imsgbuf(14) > 0) then
         call c_read_double(vunit,data_nspot,imsgbuf(14))
      endif
!
      if (checkLdaCorrection(local_id,1) .and. imsgbuf(13) > 1) then
         call insertDataPackFromInput(local_id,1,imsgbuf(13),data_ldapu)
      endif
   else
!     ----------------------------------------------------------------
      call recvMessage(imsgbuf,20,23457,getMyInputProc())
      call recvMessage(cmsgbuf,imsgbuf(2),23458,getMyInputProc())
      call recvMessage(fspace,imsgbuf(3),23459,getMyInputProc())
      call recvMessage(evec,3,23465,getMyInputProc())
!     ----------------------------------------------------------------
      if (imsgbuf(1) < 0) then
         if (dsize_ldapu < imsgbuf(13)) then
            if ( allocated(data_ldapu) ) then
               deallocate( data_ldapu )
            endif
            allocate( data_ldapu( 1:imsgbuf(13) ) )
            dsize_ldapu = imsgbuf(13)
         endif
         if (dsize_nspot < imsgbuf(14)) then
            if ( allocated(data_nspot) ) then
               deallocate( data_nspot )
            endif
            allocate( data_nspot( 1:imsgbuf(14) ) )
            dsize_nspot = imsgbuf(14)
         endif
      endif
      if (imsgbuf(13) > 0) then  
!        -------------------------------------------------------------
         call recvMessage(data_ldapu,imsgbuf(13),23466,getMyInputProc())
!        -------------------------------------------------------------
      endif
      if (imsgbuf(14) > 0) then  
!        -------------------------------------------------------------
         call recvMessage(data_nspot,imsgbuf(14),23467,getMyInputProc())
!        -------------------------------------------------------------
      endif
      if (checkLdaCorrection(local_id,1) .and. imsgbuf(13) > 1) then
         call insertDataPackFromInput(local_id,1,imsgbuf(13),data_ldapu)
      endif
!
      jmtmax=imsgbuf(8)
      jwsmax=imsgbuf(9)
      numcmax=imsgbuf(10)
   endif
!
   iform = imsgbuf(1)
!
   if (checking_spin) then
      nspin_in = imsgbuf(4)
      if (nspin_in /= nspin) then
         write(6,'(a,i6)')'The spin parameter in the potential file = ',nspin_in
         n_spin_pola_in = min(2,nspin_in)
      else
         n_spin_pola_in = n_spin_pola
      endif
   else
      n_spin_pola_in = n_spin_pola
   endif
!
!  ===================================================================
!  decoding cmsgbuf and fspace.....................................
!  ===================================================================
   do is=1,n_spin_pola_in
!     In case of n_spin_pola_in /= n_spin_pola in the potential data file
      if (n_spin_pola_in /= n_spin_pola) then
         ns = 1
      else
         ns = is
      endif
!     ----------------------------------------------------------------
      call potredg(jmtmax, jwsmax, nr, numcmax, imsgbuf(5:7),               &
                   cmsgbuf((160+40*numcmax)*(ns-1)+1:(160+40*numcmax)*ns),  &
                   fspace((9+jmtmax+jwsmax)*(ns-1)+1:(9+jmtmax+jwsmax)*ns), &
                   alat,efpot,                                              &
                   jmt,jws,r_mesh,vr(1:nr,ns),vdif,rhotot(1:nr,ns),         &
                   xvalws(ns),ztotss,zcorss,numc,                           &
                   nc(1:numcmax,ns),lc(1:numcmax,ns),                       &
                   kc(1:numcmax,ns),ec(1:numcmax,ns),                       &
                   lst(1:numcmax,ns),header,iprint,istop)
!     ----------------------------------------------------------------
   enddo
!
   if (n_spin_pola_in < n_spin_pola) then
      cmsgbuf((160+40*numcmax)+1:(160+40*numcmax)*2) = cmsgbuf(1:160+40*numcmax)
      fspace((9+jmtmax+jwsmax)+1:(9+jmtmax+jwsmax)*2)= fspace(1:9+jmtmax+jwsmax)
      vr(1:nr,2) = vr(1:nr,1)
      rhotot(1:nr,1) = HALF*rhotot(1:nr,1)
      rhotot(1:nr,2) = rhotot(1:nr,1)
      xvalws(1) = HALF*xvalws(1)
      xvalws(2) = xvalws(1)
      nc(1:numcmax,2) = nc(1:numcmax,1)
      lc(1:numcmax,2) = lc(1:numcmax,1)
      kc(1:numcmax,2) = kc(1:numcmax,1)
      ec(1:numcmax,2) = ec(1:numcmax,1)
   else if (n_spin_pola_in > n_spin_pola) then
      rhotot(1:nr,1) = TWO*rhotot(1:nr,1)
      xvalws(1) = TWO*xvalws(1)
   endif
!
!  ===================================================================
!  check if the moment is negative.................................
!  if necessary, make the moment positive, put evec to its opposite 
!  direction, and rearrange the potential and density..............
!  ===================================================================
   if(xvalws(1)-xvalws(n_spin_pola).lt.-0.001 .and. nspin.eq.3)then
      vdif=-vdif
      rtmp=xvalws(1)
      xvalws(1)=xvalws(n_spin_pola)
      xvalws(n_spin_pola)=rtmp
      evec(1)=-evec(1)
      evec(2)=-evec(2)
      evec(3)=-evec(3)
!     ----------------------------------------------------------------
!     call dcopy(jmt,vr(1,1),1,fspace,1)
      fspace(1:jmt) = vr(1:jmt,1)
!     ----------------------------------------------------------------
!     call dcopy(jmt,vr(1,n_spin_pola),1,vr(1,1),1)
      vr(1:jmt,1) = vr(1:jmt,n_spin_pola)
!     ----------------------------------------------------------------
!     call dcopy(jmt,fspace,1,vr(1,n_spin_pola),1)
      vr(1:jmt,n_spin_pola) = fspace(1:jmt)
!     ----------------------------------------------------------------
!     call dcopy(jws,rhotot(1,1),1,fspace,1)
      fspace(1:jws) = rhotot(1:jws,1)
!     ----------------------------------------------------------------
!     call dcopy(jws,rhotot(1,n_spin_pola),1,rhotot(1,1),1)
      rhotot(1:jws,1) = rhotot(1:jws,n_spin_pola)
!     ----------------------------------------------------------------
!     call dcopy(jws,fspace,1,rhotot(1,n_spin_pola),1)
      rhotot(1:jws,n_spin_pola) = fspace(1:jws)
!     ----------------------------------------------------------------
      do i=1,numc
         itmp=nc(i,1)
         nc(i,1)=nc(i,n_spin_pola)
         nc(i,n_spin_pola)=itmp
         itmp=lc(i,1)
         lc(i,1)=lc(i,n_spin_pola)
         lc(i,n_spin_pola)=itmp
         itmp=kc(i,1)
         kc(i,1)=kc(i,n_spin_pola)
         kc(i,n_spin_pola)=itmp
         rtmp=ec(i,1)
         ec(i,1)=ec(i,n_spin_pola)
         ec(i,n_spin_pola)=rtmp
         ctmp=lst(i,1)
         lst(i,1)=lst(i,n_spin_pola)
         lst(i,n_spin_pola)=ctmp
      enddo
   endif
!
   call resetCommunicator()
!
!  ===================================================================
!  Needs to take care the radial mesh interpolation...
!  ===================================================================
   if ( dsize_nspot > 0) then
!     ================================================================
!     perform interpolation: added by Yang @08/06/14
!     Note : imsgbuf(4) = nspin_old; imsgbuf(17) = jmax_pot_old
!     ================================================================
      nr_old = dsize_nspot/(2*imsgbuf(17)*n_spin_pola)
      allocate(r_mesh_old(nr_old), c_data_nspot(dsize_nspot/2))
      jmt0=imsgbuf(5)
      xst=fspace(8)
      xmt=fspace(9)
      hh=(xmt-xst)/real(jmt0-1,kind=RealKind) ! Assuming hin = hout
      do ir=1,nr_old
         r_mesh_old(ir)=exp(xst+(ir-1)*hh)
      enddo
      c_data_nspot = transfer(data_nspot,c_data_nspot)
      nj = nr_old*imsgbuf(17)
      do ns = 1, n_spin_pola
         if (nr /= nr_old .or. abs(r_mesh(jmt)-r_mesh_old(jmt0)) > TEN2m6) then
            do jl = 1, min(imsgbuf(17),jmax_pot)
               p_data_nspot => c_data_nspot(nj*(ns-1)+nr_old*(jl-1)+1:nj*(ns-1)+nr_old*jl)
               do ir = 1, nr
!                 ----------------------------------------------------
                  pot_l(ir,jl,ns) = getInterpolation(nr_old,r_mesh_old,p_data_nspot, &
                                                     r_mesh(ir),dvr)      
!                 ----------------------------------------------------
               enddo
            enddo
         else
            p_data_nspot => c_data_nspot(nj*(ns-1)+1:nj*ns)
!           ----------------------------------------------------------
            call zcopy(nr*min(imsgbuf(17),jmax_pot),p_data_nspot,1,pot_l(1,1,ns),1)
!           ----------------------------------------------------------
         endif
      enddo
      deallocate(r_mesh_old, c_data_nspot)
   endif
!
   if ( allocated(data_ldapu) ) then
      deallocate( data_ldapu )
   endif
   if ( allocated(data_nspot) ) then
      deallocate( data_nspot )
   endif
!
   end subroutine getpotg
