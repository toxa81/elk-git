module mod_expigqr
use mod_wannier
implicit none

! if wave-function is a 2-component spinor then e^{-i(G+q)x} must be 
! a 2x2 matrix in spin space; expigqr22 controls the valuse of this 2x2 matrix
! expigqr22=1 : diagonal matrix for charge response
integer expigqr22
data expigqr22/1/

! total number of matrix elements <nk|e^{-i(G+q)x}|n'k+q> in the Bloch basis
! for a given local k-point
integer, allocatable :: nmegqblh(:)

! bands (n,n') for matrix elements <nk|e^{-i(G+q)x}|n'k+q>  
!   1-st index :  1: n at k
!                 2: n' at k+q
!   2-nd index : global index of pair of bands (n,n')
!   3-rd index : k-point
integer, allocatable :: bmegqblh(:,:,:)

! matrix elements <nk|e^{-i(G+q)x}|n'k+q> in the Bloch basis
!   1-st index : local index of pair of bands (n,n')
!   2-nd index : G-vector
!   3-rd index : k-point
complex(8), allocatable :: megqblh(:,:,:)

! adjoint matrix elements <n,k-q|e^{-i(G+q)x}|n'k> in the Bloch basis
complex(8), allocatable :: amegqblh(:,:,:)

! number of adjoint matrix elements 
integer, allocatable :: namegqblh(:)

! band indices of adjoint matrix elements
integer, allocatable :: bamegqblh(:,:,:)

! interval of bands to take for matrix elements <nk|e^{-i(G+q)x}|n'k+q>
real(8) megq_include_bands(2)
data megq_include_bands/-100.1d0,100.1d0/

! minimum interband transition energy
real(8) lr_min_e12

! low level switch: compute matrix elements of e^{i(G+q)x} in the basis of
!   Wannier functions; depends on crpa and wannier_chi0_chi
logical wannier_megq
data wannier_megq/.false./

! minimum and maximum cutoff values for matrix elements in Wannier basis
real(8) megqwan_cutoff(2)
data megqwan_cutoff/-0.0d0,1000.d0/

real(8) megqwan_mindist
data megqwan_mindist/-0.0d0/
real(8) megqwan_maxdist
data megqwan_maxdist/0.1d0/

complex(8), allocatable :: megqwan(:,:)

integer nwann_include
data nwann_include/0/
integer, allocatable :: iwann_include(:)

integer nmegqblhwanmax
integer, allocatable :: nmegqblhwan(:)
integer, allocatable :: imegqblhwan(:,:)

complex(8), allocatable :: wann_c_jk(:,:,:)

integer ngntujumax
integer, allocatable :: ngntuju(:,:)
integer(2), allocatable :: igntuju(:,:,:,:)
complex(8), allocatable :: gntuju(:,:,:)

! array for k+q points
!  1-st index: index of k-point in BZ
!  2-nd index: 1: index of k'=k+q-K
!              2: index of K-vector which brings k+q to first BZ
!              3: index of k'=k-q point
integer, allocatable :: idxkq(:,:)

type(wannier_transitions) :: megqwantran

contains

!==============================================================================!
! the subroutine computes <psi_{n,k}|e^{-i(G+q)x}|psi_{n',k+q}>  and           !
!  <W_n|e^{-i(G+q)x}|W_{n'T}>                                                  !  
!==============================================================================!
subroutine genmegq(iq,tout,tg0q,allibt)
use modmain
use mod_nrkp
use mod_addons_q
use mod_wannier
implicit none
integer, intent(in) :: iq
logical, intent(in) :: tout
logical, intent(in) :: tg0q
logical, intent(in) :: allibt
! allocatable arrays
integer, allocatable :: igkignr_jk(:)
complex(8), allocatable :: wfsvmt_jk(:,:,:,:,:)
complex(8), allocatable :: compact_wfsvmt_jk(:,:,:)
complex(8), allocatable :: wfsvit_jk(:,:,:)
integer ngknr_jk
integer i,ikstep,sz,ig
integer nkstep,ik,ist1,ist2,ikloc
real(8) t1,t2,t3,t4,t5,dn1
integer lmaxexp,lmmaxexp
integer np
character*100 :: qnm,qdir,fout
integer, allocatable :: waninc(:)
call papi_timer_start(pt_megq)

! maximum l for exponent expansion
lmaxexp=lmaxvr
lmmaxexp=(lmaxexp+1)**2
call getqdir(iq,vqm(:,iq),qdir)
call getqname(vqm(:,iq),qnm)
qnm=trim(qdir)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k/)).and.tout) then
  wproc=.true.
  fout=trim(qnm)//"_ME.OUT"
  open(150,file=trim(fout),form="formatted",status="replace")
endif

if (wproc) then
  write(150,*)
  write(150,'("Calculation of matrix elements:")')
  write(150,'("  <n,k|e^{-i(G+q)x}|n'',k+q>")')
endif

if (wannier_megq) then
  call deletewantran(megqwantran)
  allocate(waninc(nwantot))
  if (nwann_include.eq.0) then
    waninc=1
  else
    waninc=0
    do i=1,nwann_include
      waninc(iwann_include(i))=1
    enddo
  endif
  call genwantran(megqwantran,megqwan_mindist,megqwan_maxdist,waninc=waninc)
  deallocate(waninc)
  !call printwantran(megqwantran)
  if (wproc) then
    write(150,*)
    write(150,'("Number of Wannier transitions : ",I6)')megqwantran%nwt
    write(150,'("Translation limits : ",6I6)')megqwantran%tlim(:,1), &
      &megqwantran%tlim(:,2),megqwantran%tlim(:,3)
    call flushifc(150)
  endif
endif

call timer_start(t_init_kgq,reset=.true.)
! initialize G+q vector arays
call init_gq(iq,lmaxexp,lmmaxexp,tg0q)
! initialize k+q array
call init_kq(iq)
! initialize interband transitions
call init_band_trans(allibt)
! initialize Gaunt-like coefficients 
call init_gntuju(iq,lmaxexp)
call timer_stop(t_init_kgq)
if (wproc) then
  write(150,*)
  write(150,'("maximum |G+q| [1/a.u.]                        : ",G18.10)')gqmax  
  write(150,'("number of G-vectors                           : ",I4)')ngq(iq)   
  write(150,'("number of G+q-shells                          : ",I4)')ngqsh(iq)   
  write(150,*)
  write(150,'("q-vector (lat.coord.)                         : ",&
    & 3G18.10)')vqlnr(:,iq)
  write(150,'("q-vector (Cart.coord.) [1/a.u.]               : ",&
    & 3G18.10)')vqcnr(:,iq)
  t1=sqrt(vqcnr(1,iq)**2+vqcnr(2,iq)**2+vqcnr(3,iq)**2)
  write(150,'("q-vector length [1/a.u.]                      : ",&
    & G18.10)')t1
  write(150,'("q-vector length [1/A]                         : ",&
    & G18.10)')t1/au2ang
  write(150,'("G-vector to reduce q to first BZ (lat.coord.) : ",&
    & 3I4)')ivg(:,ig0q(iq))
  write(150,'("global index of Gq-vector                     : ",&
    & I4)')ig0q(iq)
  write(150,'("relative index of Gq-vector                   : ",&
    & I4)')iig0q    
  write(150,'("reduced q-vector (lat.coord.)                 : ",&
    & 3G18.10)')vql(:,iq)
  write(150,'("reduced q-vector (Cart.coord.) [1/a.u.]       : ",&
    & 3G18.10)')vqc(:,iq)
  write(150,*)
  write(150,'("Bloch functions band interval (N1,N2 or E1,E2) : ",2F8.3)')&
    &megq_include_bands(1),megq_include_bands(2)
  write(150,*)
  write(150,'("Minimal energy transition (eV) : ",F12.6)')lr_min_e12*ha2ev    
  write(150,*)
  write(150,'("Approximate number of interband transitions : ",I5)')nmegqblh(1)
  if (wannier_megq) then
    write(150,*)
    write(150,'("Maximum number of interband transitions for megqwan : ",I5)')nmegqblhwanmax
  endif
  sz=int(16.d0*nstsv*nstsv*ngq(iq)*nkptnrloc/1048576.d0)
  write(150,*)
  write(150,'("Array size of matrix elements in Bloch basis (MB) : ",I6)')sz
  if (wannier_megq) then
    sz=int(16.d0*megqwantran%nwt*ngq(iq)/1048576.d0)
    write(150,*)
    write(150,'("Array size of matrix elements in Wannier basis (MB) : ",I6)')sz
  endif   
  sz=int(24.d0*ngntujumax*natmcls*ngq(iq)/1048576.d0)
  write(150,*)
  write(150,'("Maximum number of Gaunt-like coefficients : ",I8)')ngntujumax
  write(150,'("Array size of Gaunt-like coefficients (MB) : ",I6)')sz
  write(150,*)
  write(150,'("Init done in ",F8.2," seconds")')timer_get_value(t_init_kgq)
  call flushifc(150)
endif

if (allocated(megqblh)) deallocate(megqblh)
allocate(megqblh(nstsv*nstsv,ngq(iq),nkptnrloc))
megqblh(:,:,:)=zzero
allocate(wfsvmt_jk(lmmaxapw,nufrmax,natmtot,nspinor,nstsv))
allocate(compact_wfsvmt_jk(compact_wf_index%totsizemt,nspinor,nstsv))
allocate(wfsvit_jk(ngkmax,nspinor,nstsv))
allocate(igkignr_jk(ngkmax))
if (wannier_megq) then
  if (allocated(megqwan)) deallocate(megqwan)
  allocate(megqwan(megqwantran%nwt,ngq(iq)))
  megqwan(:,:)=zzero
  if (allocated(wann_c_jk)) deallocate(wann_c_jk)
  allocate(wann_c_jk(nwantot,nstsv,nkptnrloc))
endif

i=0
nkstep=mpi_grid_map(nkptnr,dim_k,x=i)
call timer_reset(t_getwfkq)
call timer_reset(t_genmegqblh)
call timer_reset(t_megqblh_mt)
call timer_reset(t_megqblh_it)
call timer_reset(t_megqblh_prod)
do ikstep=1,nkstep
! transmit wave-functions
  call timer_start(t_getwfkq)
  call getwfkq(ikstep,ngknr_jk,igkignr_jk,wfsvmt_jk,wfsvit_jk,compact_wfsvmt_jk)
  call timer_stop(t_getwfkq)
! compute matrix elements  
  call timer_start(t_genmegqblh)
  if (ikstep.le.nkptnrloc) then
    !call genmegqblh(iq,ikstep,ngknr(ikstep),ngknr_jk,igkignr(1,ikstep),&
    !               &igkignr_jk,wfsvmtnrloc(1,1,1,1,1,ikstep),wfsvmt_jk,&
    !               &wfsvitnrloc(1,1,1,ikstep),wfsvit_jk)
    call genmegqblh_v3(iq,ikstep,ngknr(ikstep),ngknr_jk,igkignr(1,ikstep),&
                      &igkignr_jk,compact_wfsvmtnrloc(1,1,1,ikstep),compact_wfsvmt_jk,&
                      &wfsvitnrloc(1,1,1,ikstep),wfsvit_jk)
  endif !ikstep.le.nkptnrloc
  call timer_stop(t_genmegqblh)
enddo !ikstep
if (wannier_megq) then
  call timer_start(t_genmegqwan,reset=.true.)
! compute matrix elements of e^{-i(G+q)x} in the basis of Wannier functions
  call genmegqwan(iq)
! sum over all k-points and interband transitions to get <n,T=0|e^{-i(G+q)x}|n',T'>
  call mpi_grid_reduce(megqwan(1,1),megqwantran%nwt*ngq(iq),dims=(/dim_k/),&
                       &all=.true.)
  megqwan=megqwan/nkptnr
  call timer_stop(t_genmegqwan)
  if (wproc) then
    write(150,*)
    write(150,'("Time for megqwan : ",F8.2)')timer_get_value(t_genmegqwan)
  endif
  !call printmegqwan(iq)
endif
!call printmegqblh(iq)
! for G=q=0: e^{iqx}=1+iqx
! from "Formalism of Bnad Theory" by E.I. Blount:
!   v=p/m
!   v=-i/\hbar [x,H]
! -i(xH-Hx)=p 
! xH-Hx=ip
! <nk|xH|n'k>-<nk|Hx|n'k>=E_{n'k}<nk|x|n'k>-E_{nk}<nk|x|n'k>
! <nk|x|n'k>*(E_{n'k}-E_{nk})=i<nk|p|n'k>
! <nk|x|n'k>=i<nk|p|n'k>/(E_{n'k}-E_{nk})
! <nk|e^{iqx}|n'k>=<nk|1+iqx|n'k>=\delta_{nn'}+iq<nk|x|n'k>
! <nk|e^{iqx}|n'k>=\delta_{nn'}-q*<nk|p|n'k>/(E_{n'k}-E_{nk})
if (vq_gamma(iq).and.allocated(pmatnrloc)) then
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do ig=1,ngq(iq)
      if (igqig(ig,iq).eq.1) then
        megqblh(:,ig,ikloc)=zzero
        do i=1,nmegqblh(ikloc)
          ist1=bmegqblh(1,i,ikloc)
          ist2=bmegqblh(2,i,ikloc)
          t1=evalsvnr(ist2,ik)-evalsvnr(ist1,ik)
          if (ist1.eq.ist2) megqblh(i,ig,ikloc)=zone
          if (abs(t1).gt.1d-8) then
            !if (t1.gt.0.d0) then
            !  t1=t1+swidth
            !else
            !  t1=t1-swidth
            !endif
            megqblh(i,ig,ikloc)=megqblh(i,ig,ikloc)-&
              &dot_product(vqc(:,iq),pmatnrloc(:,ist1,ist2,ikloc))/t1
          endif
        enddo
      endif
    enddo
  enddo
endif
!call printmegqblh(iq)
! time for wave-functions send/recieve
t1=timer_get_value(t_getwfkq)
call mpi_grid_reduce(t1,dims=(/dim_k/))
! total time for matrix elements calculation
t2=timer_get_value(t_genmegqblh)
call mpi_grid_reduce(t2,dims=(/dim_k/))
! time to precompute MT
t3=timer_get_value(t_megqblh_mt)
call mpi_grid_reduce(t3,dims=(/dim_k/))
! time to precompute IT
t4=timer_get_value(t_megqblh_it)
call mpi_grid_reduce(t4,dims=(/dim_k/))
! time to compute ME
t5=timer_get_value(t_megqblh_prod)
call mpi_grid_reduce(t5,dims=(/dim_k/))
! approximate number of matrix elements
dn1=1.d0*nmegqblh(1)*ngq(iq)*nkptnr
if (wannier_megq) dn1=dn1+1.d0*megqwantran%nwt*ngq(iq)
np=mpi_grid_dim_size(dim_k)
if (wproc) then
  write(150,*)
  write(150,'("Average time (seconds/proc)")')
  write(150,'("  send and receive wave-functions  : ",F8.2)')t1/np
  write(150,'("  compute matrix elements          : ",F8.2)')t2/np
  write(150,'("    precompute muffin-tin part     : ",F8.2)')t3/np
  write(150,'("    precompute interstitial part   : ",F8.2)')t4/np
  write(150,'("    multiply wave-functions        : ",F8.2)')t5/np
  write(150,'("Speed (me/sec/proc)                : ",F10.2)')dn1/t2
  call flushifc(150)
endif
deallocate(wfsvmt_jk)
deallocate(wfsvit_jk)
deallocate(compact_wfsvmt_jk)
deallocate(igkignr_jk)
deallocate(ngntuju)
deallocate(igntuju)
deallocate(gntuju)
call papi_timer_stop(pt_megq)
call mpi_grid_barrier((/dim_k/))
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
  close(150)
endif
return
end subroutine

subroutine get_adjoint_megqblh(iq)
use modmain
use mod_addons_q
implicit none
!
integer, intent(in) :: iq
!
integer ik,jk,ikstep,nkstep,jkloc,i,j,tag
logical need_to_recieve 
integer, allocatable :: jkmap(:,:)
!
if (allocated(amegqblh)) deallocate(amegqblh)
allocate(amegqblh(nstsv*nstsv,ngq(iq),nkptnrloc))
if (allocated(namegqblh)) deallocate(namegqblh)
allocate(namegqblh(nkptnrloc))
if (allocated(bamegqblh)) deallocate(bamegqblh)
allocate(bamegqblh(2,nstsv*nstsv,nkptnrloc))
!
nkstep=nkptnrloc
call mpi_grid_bcast(nkstep,dims=(/dim_k/))
allocate(jkmap(2,0:mpi_grid_dim_size(dim_k)-1))
do ikstep=1,nkstep
  jkmap=-1
  need_to_recieve=.false.
  ! if this processor has a k-point for this step
  if (ikstep.le.nkptnrloc) then
    ! k-point global index
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikstep)
    ! k-q point global index
    jk=idxkq(3,ik)
    ! find index of processor and a local jk index
    jkloc=mpi_grid_map(nkptnr,dim_k,x=j,glob=jk)
    ! save index of processor from which k-q point is recieved and local index of k-q point
    jkmap(1,mpi_grid_dim_pos(dim_k))=j
    jkmap(2,mpi_grid_dim_pos(dim_k))=jkloc
    ! make a local copy if jk is on the same processor
    if (j.eq.mpi_grid_dim_pos(dim_k)) then
      amegqblh(:,:,ikstep)=megqblh(:,:,jkloc)
      namegqblh(ikstep)=nmegqblh(jkloc)
      bamegqblh(:,:,ikstep)=bmegqblh(:,:,jkloc)
    else
      need_to_recieve=.true.
    endif
  endif
  call mpi_grid_reduce(jkmap(1,0),2*mpi_grid_dim_size(dim_k),dims=(/dim_k/),all=.true.,op=op_max)
  ! check who needs k-point which is stored on this processor
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (jkmap(1,i).eq.mpi_grid_dim_pos(dim_k).and.mpi_grid_dim_pos(dim_k).ne.i) then
      jkloc=jkmap(2,i)
      ! send to proc i
      tag=(ikstep*mpi_grid_dim_size(dim_k)+i)*10
      call mpi_grid_send(megqblh(1,1,jkloc),nstsv*nstsv*ngq(iq),(/dim_k/),(/i/),tag)
      call mpi_grid_send(nmegqblh(jkloc),1,(/dim_k/),(/i/),tag+1)
      call mpi_grid_send(bmegqblh(1,1,jkloc),2*nstsv*nstsv,(/dim_k/),(/i/),tag+2)
    endif
  enddo
  if (need_to_recieve) then
    j=jkmap(1,mpi_grid_dim_pos(dim_k))
    tag=(ikstep*mpi_grid_dim_size(dim_k)+mpi_grid_dim_pos(dim_k))*10
    call mpi_grid_recieve(amegqblh(1,1,ikstep),nstsv*nstsv*ngq(iq),(/dim_k/),(/j/),tag)
    call mpi_grid_recieve(namegqblh(ikstep),1,(/dim_k/),(/j/),tag+1)
    call mpi_grid_recieve(bamegqblh(1,1,ikstep),2*nstsv*nstsv,(/dim_k/),(/j/),tag+2)
  endif
enddo
deallocate(jkmap)
call mpi_grid_barrier((/dim_k/))
return
end subroutine

subroutine getwfkq(ikstep,ngknr_jk,igkignr_jk,wfsvmt_jk,wfsvit_jk,compact_wfsvmt_jk)
use modmain
use mod_nrkp
use mod_wannier
implicit none
integer, intent(in) :: ikstep
integer, intent(out) :: ngknr_jk
integer, intent(out) :: igkignr_jk(ngkmax)
complex(8), intent(out) :: wfsvmt_jk(lmmaxapw,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: wfsvit_jk(ngkmax,nspinor,nstsv)
complex(8), intent(out) :: compact_wfsvmt_jk(compact_wf_index%totsizemt,nspinor,nstsv)
!
integer i,ik,jk,nkptnrloc1,jkloc,j,tag

! each proc knows that it needs wave-functions at jk=idxkq(1,ik) (at k'=k+q)
!
! the distribution of k-points could look like this
!                p0          p1          p2
!          +-----------+-----------+-----------+
! ikstep=1 | ik=1 jk=3 | ik=4 jk=2 | ik=7 jk=5 |
! ikstep=2 | ik=2 jk=4 | ik=5 jk=7 | ik=8 jk=6 |
! ikstep=3 | ik=3 jk=1 | ik=6 jk=8 |  -        |
!          +-----------+-----------+-----------+

! two actions are required:
! 1) each processor scans trough other processors and determines, which
!    processors require its part of k-points; during this phase it
!    executes non-blocking 'send'
! 2) each processor must know the index of other processor, from which 
!    it gets jk-point; during this phase it executes blocking 'recieve'

do i=0,mpi_grid_dim_size(dim_k)-1
! number of k-points on the processor i
  nkptnrloc1=mpi_grid_map(nkptnr,dim_k,x=i)
  if (ikstep.le.nkptnrloc1) then
! for the step ikstep processor i computes matrix elements between k-point ik 
    ik=mpi_grid_map(nkptnr,dim_k,x=i,loc=ikstep)
! and k-point jk
    jk=idxkq(1,ik)
! find the processor j and local index of k-point jkloc for the k-point jk
    jkloc=mpi_grid_map(nkptnr,dim_k,glob=jk,x=j)
    if (mpi_grid_dim_pos(dim_k).eq.j.and.mpi_grid_dim_pos(dim_k).ne.i) then
! send to i
      tag=(ikstep*mpi_grid_dim_size(dim_k)+i)*10

      call mpi_grid_send(wfsvmtnrloc(1,1,1,1,1,jkloc),&
                        &lmmaxapw*nufrmax*natmtot*nspinor*nstsv,&
                        &(/dim_k/),(/i/),tag)
      
      call mpi_grid_send(wfsvitnrloc(1,1,1,jkloc),ngkmax*nspinor*nstsv,&
                        &(/dim_k/),(/i/),tag+1)
      
      call mpi_grid_send(ngknr(jkloc),1,&
                        &(/dim_k/),(/i/),tag+2)
      
      call mpi_grid_send(igkignr(1,jkloc),ngkmax,&
                        &(/dim_k/),(/i/),tag+3)
      
      if (wannier_megq) then
        call mpi_grid_send(wanncnrloc(1,1,jkloc),nwantot*nstsv,(/dim_k/),&
                          &(/i/),tag+4)
      endif
      
      call mpi_grid_send(compact_wfsvmtnrloc(1,1,1,jkloc),&
                        &compact_wf_index%totsizemt*nspinor*nstsv,&
                        &(/dim_k/),(/i/),tag+5)
    endif
    if (mpi_grid_dim_pos(dim_k).eq.i) then
      if (j.ne.i) then
! recieve from j
        tag=(ikstep*mpi_grid_dim_size(dim_k)+i)*10
        
        call mpi_grid_recieve(wfsvmt_jk(1,1,1,1,1),&
                             &lmmaxapw*nufrmax*natmtot*nspinor*nstsv,&
                             &(/dim_k/),(/j/),tag)
        
        call mpi_grid_recieve(wfsvit_jk(1,1,1),&
                             &ngkmax*nspinor*nstsv,&
                             &(/dim_k/),(/j/),tag+1)
        
        call mpi_grid_recieve(ngknr_jk,1,&
                             &(/dim_k/),(/j/),tag+2)
        
        call mpi_grid_recieve(igkignr_jk(1),ngkmax,&
                             &(/dim_k/),(/j/),tag+3)
        
        if (wannier_megq) then
          call mpi_grid_recieve(wann_c_jk(1,1,ikstep),nwantot*nstsv,&
                               &(/dim_k/),(/j/),tag+4)
        endif
        
        call mpi_grid_recieve(compact_wfsvmt_jk(1,1,1),&
                             &compact_wf_index%totsizemt*nspinor*nstsv,&
                             &(/dim_k/),(/j/),tag+5)
      else
! local copy
        wfsvmt_jk(:,:,:,:,:)=wfsvmtnrloc(:,:,:,:,:,jkloc)
        wfsvit_jk(:,:,:)=wfsvitnrloc(:,:,:,jkloc)
        ngknr_jk=ngknr(jkloc)
        igkignr_jk(:)=igkignr(:,jkloc)
        if (wannier_megq) wann_c_jk(:,:,ikstep)=wanncnrloc(:,:,jkloc)
        compact_wfsvmt_jk(:,:,:)=compact_wfsvmtnrloc(:,:,:,jkloc)
      endif
    endif
  endif   
enddo
call mpi_grid_barrier((/dim_k/))
return
end subroutine

subroutine genmegqblh_v3(iq,ikloc,ngknr1,ngknr2,igkignr1,igkignr2,wfsvmt1,wfsvmt2,&
                        &wfsvit1,wfsvit2)
use modmain
use mod_addons_q
use mod_nrkp
implicit none
!
integer, intent(in) :: iq
integer, intent(in) :: ikloc
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: wfsvmt1(compact_wf_index%totsizemt,nspinor,nstsv)
complex(8), intent(in) :: wfsvmt2(compact_wf_index%totsizemt,nspinor,nstsv)
complex(8), intent(in) :: wfsvit1(ngkmax,nspinor,nstsv)
complex(8), intent(in) :: wfsvit2(ngkmax,nspinor,nstsv)
!
integer ivg1(3)
integer i,j,ik,jk,igkq,n1,ispn1,ispn2,ist1,ist2,ic
integer ig,ig1,ig2,ias,ifg,ir,offs,sz
logical l1
complex(8), allocatable :: wf1gq(:,:)
complex(8), allocatable :: wftmp2(:,:)
complex(8), allocatable :: wfir1(:)
complex(8) b1(lmmaxapw*nufrmax),b2(lmmaxapw*nufrmax)
integer sizewf
integer num_ifg
integer, allocatable :: ifg_list(:)
integer, allocatable :: map_ifg(:)
logical found
complex(8), allocatable :: wf1pw(:,:,:)
!
! global k-point
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! jk=k+q-G_q
jk=idxkq(1,ik)
! G_q vector 
igkq=idxkq(2,ik)

sizewf=compact_wf_index%totsizemt+ngknr2
allocate(wf1gq(sizewf,ngq(iq)))
allocate(wftmp2(sizewf,nstsv))

allocate(ifg_list(ngrtot))
allocate(map_ifg(ngrtot))
num_ifg=0
do ig=1,ngq(iq)
  do ig2=1,ngknr2
    ivg1(:)=ivg(:,igkignr2(ig2))-ivg(:,igqig(ig,iq))-ivg(:,igkq)
    ifg=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
    found=.false.
    do i=1,num_ifg
      if (ifg_list(i).eq.ifg) then
        map_ifg(ifg)=i
        found=.true.
      endif
    enddo
    if (.not.found) then
      num_ifg=num_ifg+1
      ifg_list(num_ifg)=ifg
      map_ifg(ifg)=num_ifg
    endif
  enddo
enddo

call papi_timer_start(pt_megqblh_it)
call timer_start(t_megqblh_it)
allocate(wf1pw(num_ifg,nspinor,nstsv))
!$OMP PARALLEL DEFAULT(shared) PRIVATE(wfir1,ispn1,ig1,ifg,ir,ig,ig2,ivg1)
allocate(wfir1(ngrtot))
!$OMP DO
do ist1=1,nstsv
  do ispn1=1,nspinor
    wfir1=zzero
    do ig1=1,ngknr1
      ifg=igfft(igkignr1(ig1))
      wfir1(ifg)=wfsvit1(ig1,ispn1,ist1)
    enddo
    call zfftifc(3,ngrid,1,wfir1)
    do ir=1,ngrtot
      wfir1(ir)=wfir1(ir)*cfunir(ir)
    enddo
    call zfftifc(3,ngrid,-1,wfir1)
    do ig=1,ngq(iq)
      do ig2=1,ngknr2
  ! G1=G2-G-Gkq
        ivg1(:)=ivg(:,igkignr2(ig2))-ivg(:,igqig(ig,iq))-ivg(:,igkq)
        ifg=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
        wf1pw(map_ifg(ifg),ispn1,ist1)=dconjg(wfir1(ifg))
      enddo
    enddo
  enddo !ispnn
enddo !ist1
!$OMP END DO
deallocate(wfir1)
!$OMP END PARALLEL
call timer_stop(t_megqblh_it)      
call papi_timer_stop(pt_megqblh_it)

call papi_timer_start(pt_megqblh)
do ispn1=1,nspinor
  if (expigqr22.eq.1) ispn2=ispn1
! index of the interband transitions
  i=1
! go through the interband transitions    
  do while (i.le.nmegqblh(ikloc))
! left <bra| state 
    ist1=bmegqblh(1,i,ikloc)
    l1=.true.
    if (spinpol) then
      if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
    endif
    if (l1) then
      call timer_start(t_megqblh_mt)
      call papi_timer_start(pt_megqblh_mt)
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(ias,ic,b1,b2,j,sz,offs)
      do ig=1,ngq(iq)
! precompute muffint-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
        do ias=1,natmtot
          ic=ias2ic(ias)
          offs=compact_wf_index%offsetmt(ias)
          sz=compact_wf_index%sizemt(ias)
          b1(1:sz)=dconjg(wfsvmt1(offs+1:offs+sz,ispn1,ist1)*sfacgq(ig,ias))
          b2=zzero
          do j=1,ngntuju(ic,ig)
            b2(igntuju(2,j,ic,ig))=b2(igntuju(2,j,ic,ig))+&
                                  &b1(igntuju(1,j,ic,ig))*gntuju(j,ic,ig)
          enddo
          wf1gq(offs+1:offs+sz,ig)=b2(1:sz)
        enddo !ias
      enddo !ig  
!$OMP END PARALLEL DO
      call timer_stop(t_megqblh_mt)
      call papi_timer_stop(pt_megqblh_mt)
      do ig=1,ngq(iq)
        do ig2=1,ngknr2
! G1=G2-G-Gkq
          ivg1(:)=ivg(:,igkignr2(ig2))-ivg(:,igqig(ig,iq))-ivg(:,igkq)
          ifg=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
          wf1gq(compact_wf_index%totsizemt+ig2,ig)=wf1pw(map_ifg(ifg),ispn1,ist1)
        enddo
      enddo
    endif !l1
    call timer_start(t_megqblh_prod)
    n1=0
! collect right |ket> states into matrix wftmp2
    do while ((i+n1).le.nmegqblh(ikloc))
      if (bmegqblh(1,i+n1,ikloc).ne.bmegqblh(1,i,ikloc)) exit
      ist2=bmegqblh(2,i+n1,ikloc)
      n1=n1+1
      wftmp2(1:compact_wf_index%totsizemt,n1)=wfsvmt2(1:compact_wf_index%totsizemt,ispn2,ist2)
      do ig2=1,ngknr2
        wftmp2(compact_wf_index%totsizemt+ig2,n1)=wfsvit2(ig2,ispn2,ist2)
      enddo
    enddo !while
! update several matrix elements by doing matrix*matrix operation
!  me(ib,ig)=wftmp2(ig2,ib)^{T}*wftmp1(ig2,ig)
    call zgemm('T','N',n1,ngq(iq),sizewf,zone,wftmp2,sizewf,wf1gq,sizewf,&
              &zone,megqblh(i,1,ikloc),nstsv*nstsv)
    i=i+n1
    call timer_stop(t_megqblh_prod)
  enddo !while
enddo !ispn
deallocate(wf1gq)
deallocate(wftmp2)
deallocate(ifg_list)
deallocate(map_ifg)
deallocate(wf1pw)
call papi_timer_stop(pt_megqblh)

return
end
end module
