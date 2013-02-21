subroutine genu4(iq,ngq_,nwloc,gqcutoff,u4_)
use modmain
use mod_addons_q
use mod_wannier
use mod_expigqr
use mod_linresp
implicit none
!
integer, intent(in) :: iq
integer, intent(inout) :: ngq_
integer, intent(in) :: nwloc
real(8), intent(in) :: gqcutoff
complex(8), intent(out) :: u4_(megqwantran%nwt,megqwantran%nwt,megqwantran%ntr,nwloc)
!
integer iwloc,iw,n,n1,i,ig,vtl(3),j,it
real(8) v2(3),vtc(3),vqc1(3)
complex(8), allocatable :: vscrn(:,:)
complex(8), allocatable :: megqwan2(:,:)
complex(8), allocatable :: megqwan3(:,:)
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: zm2(:,:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: epsilon(:,:)
complex(8) zt1
integer igq_,ig0q_
complex(8), allocatable :: megqblh_(:,:,:)

ngq_=0
do ig=1,ngq(iq)
  if (gq(ig,iq).le.gqcutoff) ngq_=ngq_+1
enddo

allocate(megqblh_(nstsv*nstsv,ngq_,nkptnrloc))
igq_=0
do ig=1,ngq(iq)
  if (gq(ig,iq).le.gqcutoff) then
    igq_=igq_+1
    megqblh_(:,igq_,:)=megqblh(:,ig,:)
    if (igqig(ig,iq).eq.1) ig0q_=igq_
  endif
enddo

if (screenu4) call genchi0_v2(iq,ngq_,megqblh_)

if (vq_gamma(iq)) then
  vqc1=0.d0
else
  vqc1(:)=vqc(:,iq)
endif

call papi_timer_start(pt_uscrn)
allocate(vscrn(ngq_,ngq_))
allocate(krnl(ngq_,ngq_))
allocate(epsilon(ngq_,ngq_))
allocate(zm1(megqwantran%nwt,ngq_))
allocate(zm2(megqwantran%nwt,megqwantran%nwt))
krnl=zzero
igq_=0
do ig=1,ngq(iq)
  if (gq(ig,iq).le.gqcutoff) then
    igq_=igq_+1
    krnl(igq_,igq_)=vhgq(ig,iq)
  endif
enddo
allocate(megqwan2(ngq_,megqwantran%nwt))   
allocate(megqwan3(ngq_,megqwantran%nwt))   
! compute megqwan2=<W_n|e^{+i(G+q)x}|W_n'T'> and also rearrange megqwan
do i=1,megqwantran%nwt
  n=megqwantran%iwt(1,i)
  n1=megqwantran%iwt(2,i)
  vtl(:)=megqwantran%iwt(3:5,i)
  v2=dble(vtl)
  call r3mv(avec,v2,vtc)
  zt1=exp(-zi*dot_product(vqc1,vtc))
  j=megqwantran%iwtidx(n1,n,-vtl(1),-vtl(2),-vtl(3))
  if (j.le.0) then
    write(*,'("Error(genu4) wrong index of matrix element")')
    write(*,'(" n,n1,vtl : ",5I4)')n,n1,vtl
    call pstop
  endif
  igq_=0
  do ig=1,ngq(iq)
    if (gq(ig,iq).le.gqcutoff) then
      igq_=igq_+1
      megqwan2(igq_,i)=dconjg(megqwan(j,ig))*zt1
      megqwan3(igq_,i)=megqwan(i,ig)
    endif
  enddo
enddo
! compute 4-index U
! TODO: comments with formulas
do iwloc=1,nwloc
  iw=mpi_grid_map(lr_nw,dim_k,loc=iwloc)
! broadcast chi0
  if (screenu4) then
    call genvscrn_v2(iq,ngq_,ig0q_,chi0loc(1,1,iwloc),krnl,vscrn,epsilon)
  else
    vscrn=krnl
  endif
  call zgemm('T','N',megqwantran%nwt,ngq_,ngq_,zone,megqwan2,ngq_,&
            &vscrn,ngq_,zzero,zm1,megqwantran%nwt)
  call zgemm('N','N',megqwantran%nwt,megqwantran%nwt,ngq_,zone,zm1,&
            &megqwantran%nwt,megqwan3,ngq_,zzero,zm2,megqwantran%nwt)
  do it=1,megqwantran%ntr
    v2=dble(megqwantran%vtr(:,it))
    call r3mv(avec,v2,vtc)
    zt1=exp(-zi*dot_product(vqc1,vtc))/omega/nkptnr
    call zaxpy((megqwantran%nwt)**2,zt1,zm2(1,1),1,u4_(1,1,it,iwloc),1)
  enddo
enddo
deallocate(megqblh_)
deallocate(megqwan2)
deallocate(megqwan3)
deallocate(zm1)
deallocate(zm2)
deallocate(krnl,epsilon)
deallocate(vscrn)
call papi_timer_stop(pt_uscrn)
return
end

