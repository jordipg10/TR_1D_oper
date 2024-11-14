module m_surface_dl
!-------------------------------------------------------------------------
!
!   $Description: Represent one surface with diffuse layer model
!
!   $Use:
!
!   $Author: Sergio Andrés Bea Jofré
!
!   $License:
!
!-------------------------------------------------------------------------
use m_parentsurface
use m_phase 
use m_species
use flib_xpath
use flib_sax
use m_general_tools_cheproo
use m_constants_cheproo
!%------------------------------------------------------------
!%------------------------------------------------------------
private                             ::
!%------------------------------------------------------------
!%------------------------------------------------------------
public                              :: &
create_ &
,destroy_ &
,set_ &
,set_parent_ &
,compute_txoh_ &
,compute_dtxoh_ &
,compute_sigma0_ &
,compute_dsigma0_ &
,init_ &
,update_sk_ &
,get_num_sk_ &
,get_num_xoh_ &
,get_xoh_index_ &
,get_name_xoh_ &
,update_ &
,write_ &
,assignment(=)
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
integer, parameter          :: &
numsk=2
!%------------------------------------------------------------
!%------------------------------------------------------------
type, public::t_surface_dl
 
private                                :: 

type (t_parentsurface), pointer        :: pp              ! Pointer to parent surface
 
real*8, pointer, dimension (:,:)       :: txoh

real*8, pointer, dimension (:,:)       :: stq0
 
real*8                                 :: cte_8rtee0

real*8                                 :: cte_f_rt
 
logical                                :: locktxoh

logical                                :: lockstq0
 
end type t_surface_dl
 
!%----------------------------------------------------
!%----------------------------------------------------
interface create_
 
module procedure create_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface destroy_
 
module procedure destroy_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface set_
 
module procedure set_surface_dl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface set_parent_
 
module procedure set_parent_surface_dl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_txoh_
 
module procedure compute_txoh_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_dtxoh_
 
module procedure compute_dtxoh_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_sigma0_
 
module procedure compute_sigma0_surfdl1
module procedure compute_sigma0_surfdl2
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_dsigma0_
 
module procedure compute_dsigma0_surfdl1
module procedure compute_dsigma0_surfdl2
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface init_
 
module procedure init_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface update_sk_
 
module procedure update_sk_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_num_sk_
 
module procedure get_num_sk_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface write_
 
module procedure write_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_num_xoh_
 
module procedure get_num_xoh_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_xoh_index_
 
module procedure get_xoh_index_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_name_xoh_
 
module procedure get_name_xoh_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface update_
 
module procedure update_temp_param_surfdl
module procedure update_xoh_in_cd_surfdl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface assignment(=)
 
module procedure copy_surfdl
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
contains
!%------------------------------------------------------------------------
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_surfdl &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create the surface object with double layer model
!
!   $Arguments:
!
 
type(t_surface_dl), intent(inout) :: this 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, parameter                 :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
this%pp => null ()
!%------------------------------------------------------------
this%txoh => null ()
this%stq0 => null ()
!%------------------------------------------------------------
this%cte_8rtee0=r0
this%cte_f_rt=r0
this%locktxoh=.false.
this%lockstq0=.false.
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_surfdl &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy the surface object with double layer model
!
!   $Arguments:
!
 
type(t_surface_dl), intent(inout) :: this 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, parameter :: &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
this%pp => null ()
this%locktxoh=.false.
this%lockstq0=.false.
this%cte_8rtee0=r0
this%cte_f_rt=r0
!%------------------------------------------------------------
!% Deallocate pointers
!%------------------------------------------------------------
call check_pointer_ (this%txoh,1,1,.false.)
call check_pointer_ (this%stq0,1,1,.false.)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_surfdl &
   (this, &
    temp, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(inout)   :: this

real*8, intent(in)                  :: temp

logical, intent(out)                :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)   :: &
 nameprop
real*8               :: &
 epsiw, &
 tk
character(len=100)   :: &
 msg 
real*8, parameter    :: &
 r8=8.0d0, &
 zerokelvin=273.15d0, &
 r3=1.0d3  
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
tk=temp+zerokelvin
if (this%pp%lockpaqph) then
 nameprop='dielectriccte'
 call update_ (this%pp, temp,iserror)
 if (iserror) return
 call get_prop_ (this%pp%paqph,epsiw,nameprop,iserror)
 if (iserror) goto 10
 this%cte_8rtee0=r8*r3*tk*rgas*epsiz*epsiw 
else
 msg='Error, pointer to aqueous phase not locked'
 goto 10
end if
!%-----------------------------------------------------------
this%cte_f_rt=faraday/(rgas*tk)
!%-----------------------------------------------------------
return
10 print *,'********************************'
print *,'Surface:'
print *,'Name:',this%pp%name
print *,'Service: update_'
print *, msg
 print *,'********************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_xoh_in_cd_surfdl &
   (this, &
    cd, &
    sk)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)   :: this

real*8, intent(in)               :: sk(this%pp%numsite*numsk)

real*8, intent(inout)            :: cd(this%pp%numsp) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer             :: &
 i, &
 isk 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
isk=0
do i=1,this%pp%numsite
 cd(this%pp%idxoh(i))=sk(isk+1)
 isk=isk+numsk
end do
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_num_sk_surfdl &
   (this, &
    numskloc)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in) :: this

integer, intent(out)           :: numskloc 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
numskloc=numsk*this%pp%numsite
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_num_xoh_surfdl &
   (this, &
    numxoh)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)    :: this

integer, intent(out)              :: numxoh 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
numxoh=this%pp%numsite
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_name_xoh_surfdl &
   (this, &
    namexoh, &
    numxoh)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)             :: this

integer, intent(out)                       :: numxoh

character(len=100), pointer                :: namexoh(:)  
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 i, &
 ipri 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
call check_pointer_ (namexoh,this%pp%numsite,.true.)
!%------------------------------------------------------------
do i=1,this%pp%numsite
 ipri=this%pp%idxoh(i)
 namexoh(i)=this%pp%pspecies(ipri)%ptr%name
end do
numxoh=this%pp%numsite
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_surface_dl &
   (this, &
    parent, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(inout)         :: this

type(t_parentsurface), intent(in), target :: parent 

logical, intent(out)                      :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer             :: &
 i, &
 j, &
 isp
real*8              :: &
 z
character(len=100)  :: &
 nameprop, &
 msg
real*8, parameter   :: &
 r1=1.0d0 
!-------------------------------------------------------------------------
!
!   $code
!

!%---------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------
this%pp => parent
!%---------------------------------------------------------
if (this%pp%locksp) then
 call check_pointer_ (this%stq0,this%pp%numsite,this%pp%numsp,.true.)
 call check_pointer_ (this%txoh,this%pp%numsite,this%pp%numsp,.true.)
else
  msg='Error, species must be previously locked'
  goto 10
end if
!%----------------------------------------------------------
isp=0
do i=1,this%pp%numsite
 do j=1,this%pp%numspsite(i)
  isp=isp+1
  nameprop='charge'
  call get_prop_ (this%pp%pspecies(isp)%ptr,z,nameprop,msg,iserror)
  if (iserror) goto 10 

  this%stq0(i,isp)=z

  this%txoh(i,isp)=r1

 end do
end do
!%----------------------------------------------------------
this%lockstq0=.true.
this%locktxoh=.true.
!%----------------------------------------------------------
return
 
10 continue
print *,'********************'
print *,'Surface:'
print *,'Name:', this%pp%name
print *,'Service: set_'
print *, msg
print *,'********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_parent_surface_dl &
   (this, &
    parent)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
type(t_surface_dl), intent(inout)           :: this

type(t_parentsurface), intent(in), target   :: parent 

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%---------------------------------------------------------
this%pp => parent
!%----------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_xoh_index_surfdl &
   (this, &
    idxoh, &
    numxoh)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)     :: this

integer, intent(out)               :: numxoh

integer, pointer                   :: idxoh(:) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
call check_pointer_ (idxoh,this%pp%numsite,.true.)
idxoh=this%pp%idxoh
numxoh=this%pp%numsite
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine init_surfdl &
   (this, &
    sk, &
    cd, &
    txoh, &
    spsurfarea, &
    ionstr)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)   :: this

real*8, intent(in)               :: txoh(this%pp%numsite)

real*8, intent(in)               :: spsurfarea(this%pp%numsite)  

real*8, intent(in)               :: ionstr 

real*8, intent(inout)            :: cd(this%pp%numsp)

real*8, pointer                  :: sk(:)   
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer             :: &
 isk, &
 i 
real*8, pointer     :: &
 x(:) => null (), &
 sigma0(:) => null ()
real*8              :: &
 cte, &
 value 
real*8, parameter   :: &
 r0=0.0d0, &
 r1=1.0d0, &
 r2=2.0d0
character(len=100)  :: &
 msg  
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
msg=''
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
value=minval(spsurfarea)
if (value<=r0) then
 msg='Error, specific surface are must be > 0'
 goto 10 
end if
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
cte=dsqrt(this%cte_8rtee0*ionstr)
call check_pointer_ (sk,this%pp%numsite*numsk,.true.)
call check_pointer_ (sigma0,this%pp%numsite,.true.)
call check_pointer_ (x,this%pp%numsite,.true.)
call compute_sigma0_ (this,sigma0,cd,spsurfarea,this%pp%numsite,this%pp%numsp)
sk=r1
isk=0
do i=1,this%pp%numsite
 sk(isk+1)=cd(this%pp%idxoh(i))
 if (sk(isk+1)==r0) sk(isk+1)=txoh(i)
 x(i)=sigma0(i)/cte
 x(i)=-r2*dlog(x(i)+dsqrt(x(i)*x(i)+r1))
 sk(isk+2)=dexp(x(i))
 isk=isk+numsk
end do
!%------------------------------------------------------------
!% Update the sk in cd vector 
!%------------------------------------------------------------
call update_ (this,cd,sk)
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (x,1,.false.)
call check_pointer_ (sigma0,1,.false.)
!%------------------------------------------------------------
return
10 continue
print *,'********************'
print *,'Surface:'
print *,'Name:', this%pp%name
print *,'Service: set_'
print *, msg
print *,'********************'
stop 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_surfdl &
   (this, &
    ioutput, &
    msg, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_surface_dl), intent(in)    :: this

integer, intent(in)                :: ioutput    

logical, intent(out)               :: iserror

character(len=*), intent(out)      :: msg 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                            :: &
 i, &
 j 
real*8                             :: &
 zexch
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
write (ioutput,*) '----------------------------------------------'
write (ioutput,*) '           Surface information                '
write (ioutput,*) '           Model: Double Layer                '
!%-----------------------------------------------------------
!% Write the parent surface 
!%-----------------------------------------------------------
call write_ (this%pp,ioutput,msg,iserror)
!%-----------------------------------------------------------
!% Write the TXOH definition 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
write (ioutput,*) '----------------------------------------------'
write (ioutput,*) '      TXOH definition by sites       '
write(ioutput,1) (this%pp%pspecies(j)%ptr%name,j=1,this%pp%numsp)
do i=1,this%pp%numsite
  write(ioutput,2) (this%txoh(i,j),j=1,this%pp%numsp)
end do
!%-----------------------------------------------------------
if (iserror) goto 10
!%-----------------------------------------------------------
!% Write the STQ0 definition 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
write (ioutput,*) '----------------------------------------------'
write (ioutput,*) '      STQ0 definition by sites       '
write(ioutput,1) (this%pp%pspecies(j)%ptr%name,j=1,this%pp%numsp)
do i=1,this%pp%numsite
  write(ioutput,2) (this%stq0(i,j),j=1,this%pp%numsp)
end do
!%-----------------------------------------------------------
if (iserror) goto 10
!%-----------------------------------------------------------
write (ioutput,*) '----------------------------------------------'
!%-----------------------------------------------------------
return
10 continue  
print *,'***************************'
print *,'Surface:'
print *,'Name:',this%pp%name
print *,'Service: write_'
print *, msg
print*, '***************************'
iserror=.true.
return
1 format (7x,<this%pp%numsp>a10) 
2 format (<this%pp%numsp>f10.1) 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_sk_surfdl &
   (this, &
    sk, &
    cd, &
    dcd, &
    txoh, &
    spsurfarea, &
    str2, &
    isconvergence, &
    upmaxiter, &
    iter)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)  :: this

real*8, intent(out)             :: sk(this%pp%numsite*numsk)

real*8, intent(inout)           :: cd(this%pp%numsp)

real*8, intent(in)              :: dcd(this%pp%numsp,this%pp%numsite*numsk)

real*8, intent(in)              :: txoh(this%pp%numsite)

real*8, intent(in)              :: spsurfarea(this%pp%numsite)

real*8, intent(in)              :: str2

logical, intent(out)            :: isconvergence 

logical, intent(out)            :: upmaxiter

integer, intent(in)             :: iter 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer :: &
 jacobian(:,:) => null (), &
 residual(:) => null (), &
 expfi0(:) => null (), &
 sigma1(:) => null (), &
 dsigma1(:,:) => null (), &
 sigma2(:) => null (), &
 dsigma2(:,:) => null (), &
 txohc(:) => null (), &
 dtxohc(:,:) => null (), &
 error(:) => null ()
integer, pointer:: &
 indx(:) => null ()
integer         :: &
 numtotsk, &
 isk, &
 i
real*8          :: &
 dd, &
 value, &
 maxerror, &
 maxres 
character(len=100)  :: &
 msg
logical             :: &
 iserror
real*8, parameter   :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
if (iter>maxiterads) goto 10
isconvergence=.false.
upmaxiter=.false.
!%-----------------------------------------------------------
value=minval(spsurfarea)
if (value==r0) then
  msg='Error, specific surface area must be different than 0'
  goto 10
end if 
!%-----------------------------------------------------------
numtotsk=this%pp%numsite*numsk
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (jacobian,numtotsk,numtotsk,.true.)
call check_pointer_ (residual,numtotsk,.true.)
call check_pointer_ (indx,numtotsk,.true.)
call check_pointer_ (expfi0,this%pp%numsite,.true.)
call check_pointer_ (sigma1,this%pp%numsite,.true.)
call check_pointer_ (sigma2,this%pp%numsite,.true.)
call check_pointer_ (txohc,this%pp%numsite,.true.)
call check_pointer_ (dtxohc,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (dsigma1,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (dsigma2,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (error,numtotsk,.true.)
!%-----------------------------------------------------------
!% Initialice 
!%-----------------------------------------------------------
isk=0
do i=1,this%pp%numsite
 expfi0(i)=sk(isk+2)
 isk=isk+numsk
end do
!%-----------------------------------------------------------
!% Compute the jacobian and residual 
!%-----------------------------------------------------------
call compute_txoh_ (this,txohc,cd,this%pp%numsite)
call compute_dtxoh_ (this,dtxohc,sk,cd,dcd,this%pp%numsite,this%pp%numsp)
call compute_sigma0_(this,sigma1,cd,spsurfarea,this%pp%numsite,this%pp%numsp)
call compute_dsigma0_(this,dsigma1,cd,dcd,sk,spsurfarea,this%pp%numsite,this%pp%numsp)
call compute_sigma0_(this,sigma2,expfi0,str2,this%pp%numsite)
call compute_dsigma0_(this,dsigma2,expfi0,str2,this%pp%numsite)
isk=0
do i=1,this%pp%numsite
    residual(isk+1)=txohc(i)-txoh(i)
    jacobian(isk+1,:)=dtxohc(i,:)
    residual(isk+2)=sigma1(i)-sigma2(i)
    jacobian(isk+2,:)=dsigma1(i,:)-dsigma2(i,:)
    isk=isk+numsk
end do
residual=-residual
!%-----------------------------------------------------------
!% Compute the maximum error in residual
!%-----------------------------------------------------------
maxres=maxval(dabs(residual))
!%-----------------------------------------------------------
!% Solve the linear system 
!%-----------------------------------------------------------
call ludcmp (jacobian,numtotsk,numtotsk,indx,dd,msg,iserror)
if (iserror) goto 20
call lubksb (jacobian,numtotsk,numtotsk,indx,residual)
!%-----------------------------------------------------------
!% Compute error in unknowns
!%-----------------------------------------------------------
error=residual/sk
error=dabs(error)
maxerror=maxval(error)
if (maxerror>facmaxads) then
  residual = residual * facmaxads/maxerror
end if
sk = sk + residual
!%-----------------------------------------------------------
!% Check convergence 
!%-----------------------------------------------------------
isconvergence=(maxres<=tolresads.and.maxerror<=tolunkads)
upmaxiter=(iter>maxiterads)
!%-----------------------------------------------------------
!% Update the solution 
!%-----------------------------------------------------------
call update_ (this,cd,sk)
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (jacobian,1,1,.false.)
call check_pointer_ (residual,1,.false.)
call check_pointer_ (indx,1,.false.)
call check_pointer_ (expfi0,1,.false.)
call check_pointer_ (sigma1,1,.false.)
call check_pointer_ (sigma2,1,.false.)
call check_pointer_ (txohc,1,.false.)
call check_pointer_ (dtxohc,1,1,.false.)
call check_pointer_ (dsigma1,1,1,.false.)
call check_pointer_ (dsigma2,1,1,.false.)
call check_pointer_ (error,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
 
10 print *,'********************************'
print *,'Surface:'
print *,'Name:',this%pp%name
print *,'Service: update_'
print *, msg
print *,'********************************'
stop 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%*************************************************************
subroutine compute_dsigma0_surfdl2 &
   (this, &
    dsigma0, &
    expfi0, &
    str2, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)      :: this

integer, intent(in)                 :: numsite

real*8, intent(out)                 :: dsigma0(numsite,numsite*numsk)

real*8, intent(in)                  :: str2

real*8, intent(in)                  :: expfi0(numsite) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer         :: &
 ipri, &
 i
real*8          :: &
 fi0(numsite)
real*8          :: &
 cte, &
 cte2 
real*8, parameter :: &
 r0=0.0d0, &
 r2=2.0d0 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
cte=dsqrt(this%cte_8rtee0*str2)
cte2=this%cte_f_rt/r2
!%-----------------------------------------------------------
if (numsite/=this%pp%numsite) goto 10
fi0=-dlog(expfi0)/this%cte_f_rt
dsigma0=r0
ipri=0
do i=1,this%pp%numsite
 dsigma0(i,ipri+2)=-cte*cosh(cte2*fi0(i))/(r2*expfi0(i))
 ipri=ipri+numsk
end do
!%-----------------------------------------------------------
return
 
10 print *,'Surface:'
print *,'Service: compute_sigma0_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sigma0_surfdl2 &
   (this, &
    sigma0, &
    expfi0, &
    str2, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute tha charge in the layer 0 acccording
!%   sigma0 = (8*r*t*e*e0*I)1/2 * sinh( f *fi0)
!%                                     ---
!%                                     2rt
!%
!%   where
!%    I= Ionic strength of the solution
!%    r= Universal gas constant
!%    f= Faraday constant
!%    fi0= Electric potential in layer0
!%    e0= Permitivity of free space
!%    e= Dilectric constant of water
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)  :: this

integer, intent(in)             :: numsite

real*8, intent(out)             :: sigma0(numsite)

real*8, intent(in)              :: expfi0(numsite)

real*8, intent(in)              :: str2 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8          :: &
 cte, &
 cte2 
real*8, parameter :: &
 r0=0.0d0, &
 r2=2.0d0
real*8, pointer   :: &
 fi0(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
!% Initialice variables 
!%-----------------------------------------------------------
sigma0=r0
!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
cte=dsqrt(this%cte_8rtee0*str2)
cte2=this%cte_f_rt/r2
call check_pointer_ (fi0,this%pp%numsite,.true.)
fi0=-dlog(expfi0)/this%cte_f_rt
sigma0=cte*sinh(cte2*fi0)
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (fi0,1,.false.)
!%-----------------------------------------------------------
return
 
10 print *,'Surface:'
print *,'Service: compute_sigma0_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sigma0_surfdl1 &
   (this, &
    sigma0, &
    cd, &
    spsurfarea, &
    numsite, &
    numsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in) :: this

integer, intent(in)            :: numsite

integer, intent(in)            :: numsp

real*8, intent(out)            :: sigma0(numsite)

real*8, intent(in)             :: spsurfarea(numsite)

real*8, intent(in)             :: cd(numsp) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8          :: &
 spcarea(numsite)
integer         :: &
 i
character(len=100):: &
 nameprop 
real*8, parameter :: &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
!%-----------------------------------------------------------
sigma0=r0
sigma0=matmul(this%stq0,cd)
do i=1,this%pp%numsite
 sigma0(i)=(faraday/spsurfarea(i))*sigma0(i)
end do
!%-----------------------------------------------------------
return
 
10 print *,'Surface:'
print *,'Service: compute_sigma0_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dsigma0_surfdl1 &
   (this, &
    dsigma0, &
    cd, &
    dcd, &
    sk, &
    spsurfarea, &
    numsite, &
    numsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)  :: this

integer, intent(in)             :: numsite 

integer, intent(in)             :: numsp

real*8, intent(out)             :: dsigma0(numsite,numsite*numsk)

real*8, intent(in)              :: sk(numsite*numsk)

real*8, intent(in)              :: spsurfarea(numsite)

real*8, intent(in)              :: cd(numsp) 

real*8, intent(in)              :: dcd(numsp,numsite*numsk)
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer         :: &
 i, &
 j, &
 ipri, &
 ixoh
real*8             :: &
 spcarea 
real*8, parameter :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
dsigma0=r0
!%-----------------------------------------------------------
ipri=0
!%-----------------------------------------------------------
do i=1,numsite
 ixoh=this%pp%idxoh(i)
 do j=1,numsp
   !dsigma0(i,ipri+1)=dsigma0(i,ipri+1)+ this%stq0(i,j)*this%txoh(i,j)*cd(j)/sk(ipri+1)
   if (ixoh==j) then
    dsigma0(i,ipri+1)=this%stq0(i,j) 
   else 
    dsigma0(i,ipri+1)=this%stq0(i,j) * dcd(j,ipri+1)
   end if
   !dsigma0(i,ipri+2)=dsigma0(i,ipri+2)+ this%stq0(i,j)*this%stq0(i,j)*cd(j)/sk(ipri+2)
   dsigma0(i,ipri+2)=dsigma0(i,ipri+2)+ this%stq0(i,j) * dcd(j,ipri+2)
 end do
 
 dsigma0(i,:)=(faraday/spsurfarea(i))*dsigma0(i,:)
 
 ipri=ipri+numsk
 
end do
!%-----------------------------------------------------------
return
 
10 print *,'Surface:'
print *,'Service: compute_dsigma0_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_txoh_surfdl &
   (this, &
    txohc, &
    cd, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)  :: this

integer, intent(in)             :: numsite

real*8, intent(out)             :: txohc(numsite)

real*8, intent(in)              :: cd(this%pp%numsp) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, parameter :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
if (numsite/=this%pp%numsite) goto 10
txohc=r0
txohc=matmul(this%txoh,cd)
!%------------------------------------------------------------
return
 
10 print *,'Surface:'
print *,'Service: compute_txoh_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dtxoh_surfdl &
   (this, &
    dtxoh, &
    sk, &
    cd, &
    dcd, &
    numsite, &
    numsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the derivate of mass balnce equation with
!%    respect to xoh, exp((-f/rt)*sigam0), exp((-f/rt)*sigamb)
!%    exp((-f/rt)*sigamd)
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)    :: this

integer, intent(in)               :: numsite

integer, intent(in)               :: numsp

real*8, intent(out)               :: dtxoh(numsite,numsite*numsk)

real*8, intent(in)                :: cd(numsp)

real*8, intent(in)                :: dcd(numsp,numsite*numsk)

real*8, intent(in)                :: sk(numsite*numsk) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer            :: &
 j, &
 ipri, &
 i, &
 ixoh
real*8, parameter :: &
 r0=0.0d0, &
 r1=1.0d0
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
dtxoh=r0
ipri=0
do i=1,this%pp%numsite
 ixoh=this%pp%idxoh(i)
 do j=1,numsp
  !dtxoh(i,ipri+1)=dtxoh(i,ipri+1)+this%txoh(i,j)*this%txoh(i,j)*cd(j)/sk(ipri+1)
  !dtxoh(i,ipri+2)=dtxoh(i,ipri+2)+this%txoh(i,j)*this%stq0(i,j)*cd(j)/sk(ipri+2)
  if (ixoh==j) then
   dtxoh(i,ipri+1)=this%txoh(i,j)
  else
   dtxoh(i,ipri+1)=dtxoh(i,ipri+1)+this%txoh(i,j)*dcd(j,ipri+1)
  end if
  dtxoh(i,ipri+2)=dtxoh(i,ipri+2)+this%txoh(i,j)*dcd(j,ipri+2) 
 end do
 ipri=ipri+numsk
end do
!%------------------------------------------------------------
return
 
10 print *,'Surface:'
print *,'Service: compute_dtxoh_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_surfdl &
   (copied, &
    this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_dl), intent(in)  :: this

type(t_surface_dl), intent(out) :: copied 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 ndim1, &
 ndim2 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
if (this%locktxoh) then
 copied%locktxoh=this%locktxoh
 ndim1=size(this%txoh,1)
 ndim2=size(this%txoh,2)
 call check_pointer_ (copied%txoh,ndim1,ndim2,.true.)
 copied%txoh=this%txoh
end if
!%------------------------------------------------------------
if (this%lockstq0) then
 copied%lockstq0=this%lockstq0
 ndim1=size(this%stq0,1)
 ndim2=size(this%stq0,2)
 call check_pointer_ (copied%stq0,ndim1,ndim2,.true.)
 copied%stq0=this%stq0
end if
!%------------------------------------------------------------
copied%cte_8rtee0=this%cte_8rtee0
copied%cte_f_rt=this%cte_f_rt
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_surface_dl
