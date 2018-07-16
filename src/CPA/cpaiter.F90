!****q* main/cpaiter
!
! NAME
!   cpaiter
!
! COPYRIGHT
! 
! DESCRIPTION
!
! REFERENCES
!
! INPUTS
!
! OUTPUT
!
! NOTES
!
! WARNING
!
! PARENTS
!
! CHILDREN
!   giterbr, giterat
! 
! SOURCE
!..............................................................................
subroutine cpaiter(matrix,ndim,iter,cdel,iscf_cpa,cdim,citer,cmix,cw0,ctol)
! use mpi
! use options
  implicit none 
!io............................................................................
  integer iter,ndim
  complex (kind(0.d0)) matrix(ndim)
!local.........................................................................
  integer :: cdel,iscf_cpa,cdim,citer
  real (kind(0.d0)) :: cmix,cw0,ctol
  integer i,ibr,m,k,n
  real(kind(0.d0)) pq,p,ppq,tp,u
!allocation....................................................................
  real(kind(0.d0)), allocatable :: t(:),x(:)
  real(kind(0.d0)), allocatable :: nm(:),nml(:),fm(:),fml(:),delta(:,:,:),&
       & bkni(:,:)
  save m,k,pq,p,ppq,tp,u,t,x,nm,nml,fm,fml,delta,bkni
!calculate length of the vector................................................
  n=2*ndim
!==============================================================================
  ibr=0
  if(abs(cdel).gt.0) ibr=mod(iter,cdel)
  if(iter.lt.0) then

     if(allocated(t))     deallocate(t)
     if(allocated(x))     deallocate(x)

     if(allocated(nm))    deallocate(nm)
     if(allocated(nml))   deallocate(nml)
     if(allocated(fm))    deallocate(fm)
     if(allocated(fml))   deallocate(fml)
     if(allocated(delta)) deallocate(delta)
     if(allocated(bkni))  deallocate(bkni)

     return
  elseif((abs(cdel).gt.0.and.ibr.eq.0).or.iter.eq.0) then
     
     if(allocated(x)) deallocate(x)
     allocate (x(n))
     x(1:ndim)  = real(matrix(1:ndim))
     x(1+ndim:n)= aimag(matrix(1:ndim))

     if(iscf_cpa.eq.1) then

        if(allocated(nm)) deallocate(nm)
        if(allocated(nml)) deallocate(nml)
        if(allocated(fm)) deallocate(fm)
        if(allocated(fml)) deallocate(fml)
        if(allocated(delta)) deallocate(delta)
        if(allocated(bkni)) deallocate(bkni)

        allocate (nm(n),fm(n),nml(n),fml(n))
     
        nm(1:n)=x(1:n)
        fm(1:n)=0.d0
        m=0
        if(cdim.gt.1) allocate (delta(n,cdim,2),bkni(cdim,cdim))
        delta=0.d0
        bkni=0.d0

        return
     elseif(iscf_cpa.eq.2) then
        i=2*n*(1+cdim)+(1+cdim)**2+2*cdim**2
        if(allocated(t)) deallocate(t)
        if(allocated(x)) deallocate(x)
        allocate (t(i),x(n))
        t=0.d0
        x(1:ndim)  = real(matrix(1:ndim))
        x(1+ndim:n)= aimag(matrix(1:ndim))
     
        t=0.d0
        pq=1.d0
        p=cmix
        ppq=cmix
        k=0
        m=0
        tp=0.d0
        u=0.d0
        do i=1,n
           t(i)=x(i)
        enddo
        ppq=max(ppq,cw0)

        return
     endif
  endif

  x(1:ndim)   = real(matrix(1:ndim))
  x(1+ndim:n) = aimag(matrix(1:ndim))

  if(iscf_cpa.eq.1) then
     call giterbr(cmix,cw0,cdim,x,n,m,nm,nml,fm,fml,delta,bkni,tp)
  elseif(iscf_cpa.eq.2) then
     call giterat(citer,cdim,cmix,x,n,t,u,m,k,p,ppq,pq,cw0,tp,ctol)
  endif

  matrix(1:ndim)=cmplx(x(1:ndim),x(ndim+1:n),kind(0.d0))
  return
end subroutine cpaiter

!***


