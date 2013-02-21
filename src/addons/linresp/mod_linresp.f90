module mod_linresp
implicit none

! dimension for q-vectors
integer, parameter :: dim_q=2

! type of linear response calculation
!   0 : charge response
!   1 : magnetic response (not yet implemented)
integer lrtype
data lrtype/0/

! temporary array used in genchi0blh
complex(8), allocatable :: megqblh2(:,:)

! interval of bands to include in chi0
real(8) chi0_include_bands(2)
data chi0_include_bands/-100.1d0,100.1d0/

! interval of bands to exclude from chi0
real(8) chi0_exclude_bands(2)
data chi0_exclude_bands/100.1d0,-100.1d0/

complex(8), allocatable :: wann_cc(:,:,:)
complex(8), allocatable :: wann_cc2(:,:)

complex(8), allocatable :: chi0loc(:,:,:)
complex(8), allocatable :: chi0wanloc(:,:,:)

! number of energy-mesh points
integer lr_nw
data lr_nw/201/
! first energy point (Ha)
real(8) lr_w0
data lr_w0/0.d0/
! last energy point (Ha)
real(8) lr_w1
data lr_w1/1.d0/
! energy step
real(8) lr_dw
! energy mesh
complex(8), allocatable :: lr_w(:)
! broadening parameter (Ha)
real(8) lr_eta
data lr_eta/0.01d0/
! QP correction: 0 for G0W0, 1 for GW0 and 2 for GW
integer gw_mode
data gw_mode/0/
! energy step used in G0W0 calculation
real(8) del_e
data del_e/0.01d0/
! inverse temperature for the matsubara frequency in eV^-1
real(8) lr_beta
data lr_beta/30.d0/
! .true. if imaginary frequency mesh is required
logical timgw
data timgw/.false./
! first imaginary frequency
real(8) lr_iw0
data lr_iw0/0.d0/
! last imaginary frequency
real(8) lr_iw1
data lr_iw1/80.d0/

real(8) fxca0
data fxca0/0.d0/
real(8) fxca1
data fxca1/0.d0/
integer nfxca
data nfxca/1/
integer fxctype
data fxctype/0/

! high-level switch: compute chi0 and chi in Wannier functions basis
logical wannier_chi0_chi 
data wannier_chi0_chi/.false./
! high-level switch: .true. if chi0 should be multiplied by 2
logical wannier_chi0_afm
data wannier_chi0_afm/.false./

! indices of response functions in global array f_response(:,:,:)
integer, parameter :: f_chi0                 = 1
integer, parameter :: f_chi                  = 2
integer, parameter :: f_chi_scalar           = 3
integer, parameter :: f_chi_pseudo_scalar    = 4
integer, parameter :: f_epsilon_matrix_GqGq  = 5
integer, parameter :: f_epsilon_scalar_GqGq  = 6
integer, parameter :: f_inv_epsilon_inv_GqGq = 7
integer, parameter :: f_epsilon_eff          = 8
integer, parameter :: f_epsilon_eff_scalar   = 9
integer, parameter :: f_sigma                = 10
integer, parameter :: f_sigma_scalar         = 11
integer, parameter :: f_loss                 = 12
integer, parameter :: f_loss_scalar          = 13
integer, parameter :: f_chi0_wann            = 14
integer, parameter :: f_chi_wann             = 15
integer, parameter :: f_epsilon_eff_wann     = 16
integer, parameter :: f_sigma_wann           = 17
integer, parameter :: f_loss_wann            = 18
integer, parameter :: f_epsilon_inv_GqGq     = 19

integer, parameter :: nf_response            = 19
complex(8), allocatable :: f_response(:,:,:)

complex(8), allocatable :: u4(:,:,:,:,:)
logical screenu4
data screenu4/.true./

complex(8), allocatable :: gw_self_energy(:,:,:)
complex(8), allocatable :: self_energy_x(:,:)
contains

subroutine genchi0blh(ikloc,ngq,w,chi0w)
use modmain
use mod_nrkp
use mod_expigqr
implicit none
! arguments
integer, intent(in) :: ikloc
integer, intent(in) :: ngq
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0w(ngq,ngq)
! local variables
logical l1
integer i,ist1,ist2,ik,jk,ig
real(8) t1,t2
complex(8), allocatable :: wt(:)
logical, external :: bndint
! 
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
allocate(wt(nmegqblh(ikloc)))
wt(:)=zzero
do i=1,nmegqblh(ikloc)
  ist1=bmegqblh(1,i,ikloc)
  ist2=bmegqblh(2,i,ikloc)
! default : include all interband transitions         
  l1=.true.
! cRPA case : don't include bands in energy window [crpa_e1,crpa_e2]
  if (bndint(ist1,evalsvnr(ist1,ik),chi0_exclude_bands(1),&
      chi0_exclude_bands(2)).and.bndint(ist2,evalsvnr(ist2,jk),&
      chi0_exclude_bands(1),chi0_exclude_bands(2))) l1=.false.
  if (l1) then
    t1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
    if (abs(t1).gt.1d-6) then
      t2=sign(scissor,t1)
      wt(i)=t1/(evalsvnr(ist1,ik)-evalsvnr(ist2,jk)-t2+w)
    endif
  endif
enddo !i
call papi_timer_start(pt_megqblh2)
do ig=1,ngq
  megqblh2(1:nmegqblh(ikloc),ig)=dconjg(megqblh(1:nmegqblh(ikloc),ig,ikloc))*wt(1:nmegqblh(ikloc))
enddo
call papi_timer_stop(pt_megqblh2)
call papi_timer_start(pt_chi0_zgemm)
call zgemm('T','N',ngq,ngq,nmegqblh(ikloc),zone,megqblh(1,1,ikloc),nstsv*nstsv,&
  &megqblh2(1,1),nstsv*nstsv,zone,chi0w(1,1),ngq)
call papi_timer_stop(pt_chi0_zgemm)
deallocate(wt)
return
end subroutine

subroutine genchi0blh_v2(ikloc,ngq_,w,megqblh_,megqblh2_,chi0w)
use modmain
use mod_nrkp
use mod_expigqr
implicit none
! arguments
integer, intent(in) :: ikloc
integer, intent(in) :: ngq_
complex(8), intent(in) :: w
complex(8), intent(in) :: megqblh_(nstsv*nstsv,ngq_)
complex(8), intent(out) :: megqblh2_(nstsv*nstsv,ngq_)
complex(8), intent(out) :: chi0w(ngq_,ngq_)
! local variables
logical l1
integer i,ist1,ist2,ik,jk,ig
real(8) t1,t2
complex(8), allocatable :: wt(:)
logical, external :: bndint
! 
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
allocate(wt(nmegqblh(ikloc)))
wt(:)=zzero
do i=1,nmegqblh(ikloc)
  ist1=bmegqblh(1,i,ikloc)
  ist2=bmegqblh(2,i,ikloc)
! default : include all interband transitions         
  l1=.true.
! cRPA case : don't include bands in energy window [crpa_e1,crpa_e2]
  if (bndint(ist1,evalsvnr(ist1,ik),chi0_exclude_bands(1),&
      chi0_exclude_bands(2)).and.bndint(ist2,evalsvnr(ist2,jk),&
      chi0_exclude_bands(1),chi0_exclude_bands(2))) l1=.false.
  if (l1) then
    t1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
    if (abs(t1).gt.1d-6) then
      t2=sign(scissor,t1)
      wt(i)=t1/(evalsvnr(ist1,ik)-evalsvnr(ist2,jk)-t2+w)
    endif
  endif
enddo !i
call papi_timer_start(pt_megqblh2)
do ig=1,ngq_
  megqblh2_(1:nmegqblh(ikloc),ig)=dconjg(megqblh_(1:nmegqblh(ikloc),ig))*wt(1:nmegqblh(ikloc))
enddo
call papi_timer_stop(pt_megqblh2)
call papi_timer_start(pt_chi0_zgemm)
call zgemm('T','N',ngq_,ngq_,nmegqblh(ikloc),zone,megqblh_,nstsv*nstsv,&
          &megqblh2_,nstsv*nstsv,zone,chi0w,ngq_)
call papi_timer_stop(pt_chi0_zgemm)
deallocate(wt)
return
end subroutine

subroutine genchi0_v2(iq,ngq_,megqblh_)
use modmain
use mod_addons_q
use mod_expigqr
implicit none
! arguments
integer, intent(in) :: iq
integer, intent(in) :: ngq_
complex(8), intent(in) :: megqblh_(nstsv*nstsv,ngq_,nkptnrloc)
! local variables
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: chi0wan_k(:,:,:)
complex(8), allocatable :: chi0wan(:,:)
complex(8), allocatable :: mexp(:,:,:)
complex(8), allocatable :: megqwan_tmp(:,:)
integer i,iw,i1,i2,ikloc,n,j
integer ist1,ist2,nwloc,jwloc
integer it1(3),it2(3),it(3)
character*100 qnm,qdir,fout,fstat
integer ik,n1,n2
integer nwt
integer, allocatable :: iwt_tmp(:,:)
real(8) vtc(3)
complex(8), allocatable :: megqblh2_(:,:)
real(8) t1,t2,t3

call getqdir(iq,vqm(:,iq),qdir)
call getqname(vqm(:,iq),qnm)
qnm=trim(qdir)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k/))) then
  wproc=.true.
  fout=trim(qnm)//"_CHI0.OUT"
  open(150,file=trim(fout),form="FORMATTED",status="REPLACE")
  fstat=trim(qnm)//"_chi0_stat.txt"
endif

call papi_timer_start(pt_chi0)

if (wproc) then
  write(150,*)
  write(150,'("Calculation of chi0")')
  write(150,*)
  write(150,'("Energy mesh parameters:")')
  if (timgw) then
    write(150,'("  number of imaginary frequencies : ",I4)')lr_nw
    write(150,'("  frequency interval [eV] : ", 2F9.4)')lr_iw0,lr_iw1
  else
    write(150,'("  energy interval [eV] : ", 2F9.4)')lr_w0*ha2ev,lr_w1*ha2ev
    write(150,'("  energy step     [eV] : ", F9.4)')lr_dw*ha2ev
    write(150,'("  eta             [eV] : ", F9.4)')lr_eta*ha2ev
  endif
  write(150,*)  
  write(150,'("Included band interval (Ha)        : ",2F8.2)')&
    &chi0_include_bands(1),chi0_include_bands(2)
  write(150,'("Excluded band interval (Ha)        : ",2F8.2)')&
    &chi0_exclude_bands(1),chi0_exclude_bands(2) 
  call flushifc(150)
endif
  
allocate(chi0(ngq_,ngq_))
allocate(megqblh2_(nstsv*nstsv,ngq_))

! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)
if (allocated(chi0loc)) deallocate(chi0loc)
allocate(chi0loc(ngq_,ngq_,nwloc))

call timer_start(t_chi0_tot,reset=.true.)
! loop over energy points
do iw=1,lr_nw
  call timer_start(t_chi0_w,reset=.true.)
! sum over fraction of k-points
  chi0=zzero
  do ikloc=1,nkptnrloc
    if (nmegqblh(ikloc).gt.0) then
! for each k-point : sum over interband transitions
      call genchi0blh_v2(ikloc,ngq_,lr_w(iw),megqblh_(1,1,ikloc),megqblh2_,chi0)
    endif
  enddo
! find the processor j which will get the full chi0 and chi0wan matrices
  jwloc=mpi_grid_map(lr_nw,dim_k,glob=iw,x=j)
! sum over k-points and band transitions
  call mpi_grid_reduce(chi0(1,1),ngq_*ngq_,chi0loc(1,1,jwloc),dims=(/dim_k/),&
                      &root=(/j/))
  call timer_stop(t_chi0_w)
  if (wproc) then
    open(160,file=trim(fstat),status="REPLACE",form="FORMATTED")
    write(160,'(I8)')iw
    write(160,'(F8.2)')timer_get_value(t_chi0_w)   
    close(160)
  endif
enddo !iw
chi0loc=chi0loc/nkptnr/omega
call timer_stop(t_chi0_tot)
t1=timer_get_value(t_chi0_tot)
if (wproc) then
  write(150,*)
  write(150,'("Total time per frequency point   : ",F8.2)')t1/lr_nw
  call flushifc(150)
endif

call papi_timer_stop(pt_chi0)

call mpi_grid_barrier(dims=(/dim_k/))

deallocate(chi0)
deallocate(megqblh2_)
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
  close(150)
endif
return
end subroutine

subroutine genvscrn_v2(iq,ngq_,ig0q_,chi0,krnl,vscrn,epsilon)
use modmain
use mod_addons_q
implicit none
! arguments
integer, intent(in) :: iq
integer, intent(in) :: ngq_
integer, intent(in) :: ig0q_
complex(8), intent(in) :: chi0(ngq_,ngq_)
complex(8), intent(in) :: krnl(ngq_,ngq_)
complex(8), intent(out) :: vscrn(ngq_,ngq_)
complex(8), intent(out) :: epsilon(ngq_,ngq_)
! local variables
integer ig,ig1,ig2

call papi_timer_start(pt_vscrn)

epsilon=zzero
do ig=1,ngq_
  epsilon(ig,ig)=zone
enddo
call zgemm('N','N',ngq_,ngq_,ngq_,-zone,chi0,ngq_,krnl,ngq_,zone,epsilon,ngq_)
call invzge(epsilon,ngq_)
call zgemm('N','N',ngq_,ngq_,ngq_,zone,krnl,ngq_,epsilon,ngq_,zzero,vscrn,ngq_)
if (vq_gamma(iq)) then
  vscrn(ig0q_,ig0q_)=epsilon(ig0q_,ig0q_)*fourpi*q0wt*nkptnr*omega/(twopi**3)
  vscrn(:,:)=vscrn(:,:)/dble(nvq0)
endif
call papi_timer_stop(pt_vscrn)

return
end subroutine

end module
