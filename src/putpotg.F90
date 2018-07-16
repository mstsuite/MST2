!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine putpotg(local_id,specie,num_atoms,num_species,nspin,n_spin_pola, &
                      efermi,evec,                                      &
                      jmt,xmt,jws,xstart,jmax_pot,                      &
                      vr,vdif,rhotot,xvalws,                            &
                      atname,ztotss,zcorss,zsemss,zvalss,               &
                      numc,nc,lc,kc,ec,lst,                             &
                      jmtmax,jwsmax,nr,numcmax,                         &
                      da_ldapu, da_nsp, pot_l,                          &
                      header,wunit,fmsgbuf,iprint,istop)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO
!
   use MPPModule, only : nbsendMessage, recvMessage, waitMessage
   use MPPModule, only : setCommunicator, resetCommunicator, syncAllPEs
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use ParallelIOModule, only : isOutputProc, getMyOutputProc,          &
                                getNumOutputClients, getOutputClient,   &
                                getIOCommunicator, getMyPEinIOGroup, getNumPEsInIOGroup
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use LdaCorrectionModule, only : checkLdaCorrection, getDataPack4Output, &
                                   getNumCorrOrbitals
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: numcmax
   integer (kind=IntKind), intent(in) :: n_spin_pola
!
   character (len=*), intent(in) :: header
   character (len=*), intent(in) :: istop
   character (len=*), intent(in) :: atname
   character (len=5), intent(in) :: lst(numcmax,n_spin_pola)
!
   character (len=80) :: jtitle
   character (len=40) :: temp_buf
   character (len=1) :: cmsgbuf((160+40*numcmax)*n_spin_pola)
!
   integer (kind=IntKind), intent(in) :: local_id
   integer (kind=IntKind), intent(in) :: specie
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: num_species
   integer (kind=IntKind), intent(in) :: nspin
   integer (kind=IntKind), intent(in) :: jmt
   integer (kind=IntKind), intent(in) :: jws
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), intent(in) :: numc
   integer (kind=IntKind), intent(in) :: jmtmax
   integer (kind=IntKind), intent(in) :: jwsmax
   integer (kind=IntKind), intent(in) :: jmax_pot
   integer (kind=IntKind), intent(in) :: nc(numcmax,n_spin_pola)
   integer (kind=IntKind), intent(in) :: lc(numcmax,n_spin_pola)
   integer (kind=IntKind), intent(in) :: kc(numcmax,n_spin_pola)
   integer (kind=IntKind), intent(in) :: da_ldapu(num_atoms)
   integer (kind=IntKind), intent(in) :: da_nsp(num_atoms)
   integer (kind=IntKind), intent(in) :: wunit
   integer (kind=IntKind), intent(in) :: iprint
!
   integer (kind=IntKind) :: integer4_size
   integer (kind=IntKind) :: real8_size
   integer (kind=IntKind) :: ns
   integer (kind=IntKind) :: i, j, n
   integer (kind=IntKind) :: present_atom
   integer (kind=IntKind) :: fp_pos
   integer (kind=IntKind) :: msg_bytes
   integer (kind=IntKind) :: imsgbuf(20)
   integer (kind=IntKind) :: slen, pad_bytes
   integer (kind=IntKind) :: num_clients, proc_client, msgid(6)
   integer (kind=IntKind) :: comm, MyPEinGroup, NumPEsInGroup
!
   real (kind=RealKind), intent(in) :: efermi
   real (kind=RealKind), intent(in) :: evec(3)
   real (kind=RealKind), intent(in) :: xmt
   real (kind=RealKind), intent(in) :: xstart
   real (kind=RealKind), intent(in) :: vr(nr,n_spin_pola)
   real (kind=RealKind), intent(in) :: rhotot(nr,n_spin_pola)
   real (kind=RealKind), intent(in) :: xvalws(n_spin_pola)
   real (kind=RealKind), intent(in) :: ztotss
   real (kind=RealKind), intent(in) :: zcorss
   real (kind=RealKind), intent(in) :: zsemss
   real (kind=RealKind), intent(in) :: zvalss
   real (kind=RealKind), intent(in) :: ec(numcmax,n_spin_pola)
   real (kind=RealKind), intent(in) :: vdif
!  real (kind=RealKind), target, intent(in) :: pot_l(1:2*nr*jmax_pot*n_spin_pola)
   real (kind=RealKind), target, intent(inout) :: pot_l(1:2*jwsmax*jmax_pot*n_spin_pola)
   real (kind=RealKind), intent(out) :: fmsgbuf((9+jmtmax+jwsmax)*n_spin_pola)
!
   real (kind=RealKind) :: evec_othern(3)
   real (kind=RealKind), pointer :: data_ldapu(:), data_nspot(:)
!
   interface
      subroutine copyString2CharArray(a,b,n)
         use KindParamModule, only : IntKind
         integer (kind=IntKind), intent(in) :: n
         character (len=*), intent(in) :: a
         character (len=1), intent(out) :: b(:)
      end subroutine copyString2CharArray
   end interface
!
!  *******************************************************************
!  cccccccccc write out the potentials  ccccccccccc                   
!  *******************************************************************
!
   if (getMyOutputProc() < 0) then
      return
   endif
!
   MyPEinGroup = getMyPEinIOGroup()
   NumPEsInGroup = getNumPEsInIOGroup()
   comm = getIOCommunicator()
   call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!
!  -------------------------------------------------------------------
   call c_dtsize(integer4_size,real8_size)
!  -------------------------------------------------------------------
!  write(6,'(a,2i5)')'integer4_size, real8_size = ',integer4_size,real8_size
!
   present_atom = getGlobalIndex(local_id)
!  write(6,'(a,i5)')'present_atom = ',present_atom
!
   imsgbuf(1:20) = 0
!
   imsgbuf(1)=present_atom
   imsgbuf(4)=nspin
   imsgbuf(5)=jmt
   imsgbuf(6)=jws
   imsgbuf(7)=numc
   imsgbuf(8)=jmtmax
   imsgbuf(9)=jwsmax
   imsgbuf(10)=numcmax
!
   if (checkLdaCorrection(local_id,1) .or. maxval(da_nsp) > 0 ) then
      imsgbuf(1)=-imsgbuf(1)   ! negative sign indicates a modified data format
      imsgbuf(11) = present_atom
      if (present_atom == 1) then
         imsgbuf(12) = 0
         imsgbuf(13) = da_ldapu(present_atom)
         imsgbuf(14) = da_nsp(present_atom)
      else
         imsgbuf(12) = da_ldapu(present_atom-1) + da_nsp(present_atom-1)
         imsgbuf(13) = da_ldapu(present_atom) - da_ldapu(present_atom-1)
         imsgbuf(14) = da_nsp(present_atom) - da_nsp(present_atom-1)
      endif
   endif
!
   if (checkLdaCorrection(local_id,1)) then
      imsgbuf(16) = getNumCorrOrbitals(local_id,1)
      data_ldapu => getDataPack4Output(local_id,1,n)
      if (n /= imsgbuf(13)) then
         call ErrorHandler('putpotg','Inconsistent data size for LDA+U', &
                           n,imsgbuf(13))
      endif
   endif
!
   if ( imsgbuf(14)>0 ) then
      if (2*nr*jmax_pot*n_spin_pola /= imsgbuf(14)) then
         call ErrorHandler('putpotg','Inconsistent data size for pot_l', &
                           nr*jmax_pot*n_spin_pola,imsgbuf(14))
      endif
!     data_nspot => pot_l(1:2*nr*jmax_pot*n_spin_pola)
      data_nspot => pot_l(1:imsgbuf(14))
   endif
!
   imsgbuf(17) = jmax_pot
!
   do ns=1,n_spin_pola
!     ================================================================
!     setup cmsgbuf................................................
!     ================================================================
      write(jtitle,'('' PUTPOTG:'',i5,3x,a2,''  zt='',f3.0,           &      
  &                  '' zc='',f3.0,'' zs='',f3.0,'' zv='',f3.0,       &
  &                  '' xv='',f10.5)') present_atom,atname,           &
                                       ztotss,zcorss,zsemss,zvalss,xvalws(ns)
!     ----------------------------------------------------------------
      call copyString2CharArray(header,cmsgbuf(imsgbuf(2)+1:),80)
      call copyString2CharArray(jtitle,cmsgbuf(imsgbuf(2)+81:),80)
!     ----------------------------------------------------------------
      imsgbuf(2)=imsgbuf(2)+160
      i=imsgbuf(2)
      do j=1,numc
         write(temp_buf,'(3i5,1e20.13,1a5)')                          &
               nc(j,ns),lc(j,ns),kc(j,ns),ec(j,ns),lst(j,ns)
!        -------------------------------------------------------------
         call copyString2CharArray(temp_buf,cmsgbuf(i+(j-1)*40+1:),40)
!        -------------------------------------------------------------
      enddo
      if (numcmax > numc) then
         cmsgbuf(i+numc*40+1:i+numcmax*40) = ' '
      endif
      imsgbuf(2)=imsgbuf(2)+numcmax*40
!
!     ================================================================
!     setup fmsgbuf................................................
!     ================================================================
      i=imsgbuf(3)
      fmsgbuf(i+1)=exp(xmt)      ! alat is replaced with rmt
      fmsgbuf(i+2)=efermi
      fmsgbuf(i+3)=vdif
      fmsgbuf(i+4)=ztotss
      fmsgbuf(i+5)=zcorss
      fmsgbuf(i+6)=xvalws(ns)
      fmsgbuf(i+7)=ZERO
      fmsgbuf(i+8)=xstart
      fmsgbuf(i+9)=xmt
      imsgbuf(3)=imsgbuf(3)+9
!     ----------------------------------------------------------------
      call dcopy(jmt,vr(1,ns),1,fmsgbuf(imsgbuf(3)+1),1)
!     ----------------------------------------------------------------
      imsgbuf(3)=imsgbuf(3)+jmtmax
!     ----------------------------------------------------------------
      call dcopy(jws,rhotot(1,ns),1,fmsgbuf(imsgbuf(3)+1),1)
!     ----------------------------------------------------------------
      imsgbuf(3)=imsgbuf(3)+jwsmax
   enddo
!
!  ===================================================================
!  using gopen to perform the global openning......................
!  ===================================================================
   if( isOutputProc() ) then
!     ================================================================
!     write out imsgbuf...............................................
!     ================================================================
      fp_pos=(present_atom-1)*integer4_size*10+1
!     ----------------------------------------------------------------
!     write(6,'(a,i10)')'#1: fp_pos = ',fp_pos
!     write(6,'(a,10i6)')'local IMSGBUF(1-10) = ',imsgbuf(1:10)
      call c_fseek(wunit,fp_pos,0)
      call c_write_integer(wunit,imsgbuf,10)
      if (imsgbuf(1) < 0) then
         fp_pos=num_atoms*integer4_size*10+fp_pos
         call c_fseek(wunit,fp_pos,0)
         call c_write_integer(wunit,imsgbuf(11:20),10)
      endif
!     ----------------------------------------------------------------
!
!     ================================================================
!     write out cmsgbuf and fmsgbuf.................................
!     ================================================================
      call c_string_padsize(wunit,imsgbuf(2),pad_bytes)
      msg_bytes=imsgbuf(2)+pad_bytes+(imsgbuf(3)+3)*real8_size
      if (imsgbuf(1) > 0) then
         fp_pos=num_atoms*integer4_size*10+(present_atom-1)*msg_bytes+1
      else
         fp_pos=imsgbuf(12)*real8_size +                              &
                num_atoms*integer4_size*20+(present_atom-1)*msg_bytes+1
      endif
!     ----------------------------------------------------------------
!     write(6,'(a,2i10)')'#2: fp_pos = ',fp_pos,msg_bytes
!     write(6,'(a,10d12.5)')'local FMSGBUF(1-10) = ',fmsgbuf(1:10)
      call c_fseek(wunit,fp_pos,0)
      call c_write_string(wunit,cmsgbuf,imsgbuf(2),slen)
!     write(6,'(a,2i8)')'imsgbuf, slen = ',imsgbuf(2),slen
      call c_write_double(wunit,fmsgbuf,imsgbuf(3))
      call c_write_double(wunit,evec,3)
      if (imsgbuf(13) > 0) then
         call c_write_double(wunit,data_ldapu,imsgbuf(13))
      endif
      if (imsgbuf(14) > 0) then
         call c_write_double(wunit,data_nspot,imsgbuf(14))
      endif
!     ----------------------------------------------------------------
!
      num_clients = getNumOutputClients()
      do i = 1, num_clients
         proc_client = getOutputClient(i)
         present_atom = getGlobalIndex(local_id,proc_client)
!        -------------------------------------------------------------
         call recvMessage(imsgbuf,20,23466,proc_client)
         call recvMessage(cmsgbuf,imsgbuf(2),23467,proc_client)
         call recvMessage(fmsgbuf,imsgbuf(3),23468,proc_client)
         call recvMessage(evec_othern,3,23469,proc_client)
         if (imsgbuf(13) > 1) then
            call recvMessage(data_ldapu,imsgbuf(13),23470,proc_client)
         endif
         if (imsgbuf(14) > 1) then
            data_nspot => pot_l(1:imsgbuf(14))
            call recvMessage(data_nspot,imsgbuf(14),23471,proc_client)
         endif
!        -------------------------------------------------------------
         fp_pos=(present_atom-1)*integer4_size*10+1
!        -------------------------------------------------------------
         call c_fseek(wunit,fp_pos,0)
         call c_write_integer(wunit,imsgbuf(1:10),10)
         if (imsgbuf(1) < 0) then
            fp_pos=num_atoms*integer4_size*10+fp_pos
            call c_fseek(wunit,fp_pos,0)
            call c_write_integer(wunit,imsgbuf(11:20),10)
         endif
!        -------------------------------------------------------------
!
!        =============================================================
!        write out cmsgbuf and fmsgbuf...............................
!        =============================================================
         call c_string_padsize(wunit,imsgbuf(2),pad_bytes)
         msg_bytes=imsgbuf(2)+pad_bytes+(imsgbuf(3)+3)*real8_size
         if (imsgbuf(1) > 0) then
            fp_pos=num_atoms*integer4_size*10+(present_atom-1)*msg_bytes+1
         else
            fp_pos=imsgbuf(12)*real8_size +                           &
                   num_atoms*integer4_size*20+(present_atom-1)*msg_bytes+1
         endif
!        -------------------------------------------------------------
         call c_fseek(wunit,fp_pos,0)
         call c_write_string(wunit,cmsgbuf,imsgbuf(2),slen)
         call c_write_double(wunit,fmsgbuf,imsgbuf(3))
         call c_write_double(wunit,evec_othern,3)
         if (imsgbuf(13) > 0) then
            call c_write_double(wunit,data_ldapu,imsgbuf(13))
         endif
         if (imsgbuf(14) > 0) then
            call c_write_double(wunit,data_nspot,imsgbuf(14))
         endif
!        -------------------------------------------------------------
      enddo
   else
!     ----------------------------------------------------------------
      msgid(1) = nbsendMessage(imsgbuf,20,23466,getMyOutputProc())
      msgid(2) = nbsendMessage(cmsgbuf,imsgbuf(2),23467,getMyOutputProc())
      msgid(3) = nbsendMessage(fmsgbuf,imsgbuf(3),23468,getMyOutputProc())
      msgid(4) = nbsendMessage(evec,3,23469,getMyOutputProc())
      if (imsgbuf(13) > 1) then
         msgid(5) = nbsendMessage(data_ldapu,imsgbuf(13),23470,getMyOutputProc())
      endif
      if (imsgbuf(14) > 1) then
         msgid(6) = nbsendMessage(data_nspot,imsgbuf(14),23471,getMyOutputProc())
      endif
!     ----------------------------------------------------------------
      call waitMessage(msgid(1))
      call waitMessage(msgid(2))
      call waitMessage(msgid(3))
      call waitMessage(msgid(4))
      if (imsgbuf(13) > 1) then
         call waitMessage(msgid(5))
      endif
      if (imsgbuf(14) > 1) then
         call waitMessage(msgid(6))
      endif
   endif
!
   call resetCommunicator()
!
   end subroutine putpotg
