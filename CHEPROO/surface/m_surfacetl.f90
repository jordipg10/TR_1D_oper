module m_surface_tl
!-------------------------------------------------------------------------
!
!   $Description: Represent the surface with triple layer model 
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
private    ::
!%------------------------------------------------------------
!%------------------------------------------------------------
public:: &
create_ &
,destroy_ &
,set_ &
,set_parent_ &
,compute_txoh_ &
,compute_dtxoh_ &
,compute_sigma0_ &
,compute_sigmab_ &
,compute_sigmad_ &
,compute_dsigma0_ &
,compute_dsigmab_ &
,compute_dsigmad_ &
,init_ &
,update_sk_ &
,get_num_tot_sk_ &
,get_xoh_index_ &
,get_num_xoh_ &
,get_name_xoh_ &
,update_ &
,write_ &
,assignment(=)
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
integer, parameter          :: &
numsk=4
!%------------------------------------------------------------
!%------------------------------------------------------------
type, public::t_surface_tl
 
private                             ::
 
type (t_parentsurface), pointer     :: pp              ! Pointer to parent Interface
 
real*8, pointer, dimension(:,:)     :: txoh

real*8, pointer, dimension(:,:)     :: stq0

real*8, pointer, dimension(:,:)     :: stqb
 
real*8                              :: cte_8rtee0

real*8                              :: cte_f_rt
 
logical                             :: locktxoh

logical                             :: lockstq0

logical                             :: lockstqb
 
end type t_surface_tl
 
!%----------------------------------------------------
!%----------------------------------------------------
interface create_
 
module procedure create_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface destroy_
 
module procedure destroy_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface set_
 
module procedure set_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface set_parent_
 
module procedure set_parent_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_txoh_
 
module procedure compute_txoh_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_dtxoh_
 
module procedure compute_dtxoh_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_sigma0_
 
module procedure compute_sigma01_surftl
module procedure compute_sigma02_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_sigmab_
 
module procedure compute_sigmab1_surftl
module procedure compute_sigmab2_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_sigmad_
 
module procedure compute_sigmad1_surftl
module procedure compute_sigmad2_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_dsigma0_
 
module procedure compute_dsigma01_surftl
module procedure compute_dsigma02_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_dsigmab_
 
module procedure compute_dsigmab1_surftl
module procedure compute_dsigmab2_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_dsigmad_
 
module procedure compute_dsigmad1_surftl
module procedure compute_dsigmad2_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface init_
 
module procedure init_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface update_sk_
 
module procedure update_sk_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_num_tot_sk_
 
module procedure get_num_tot_sk_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface write_
 
module procedure write_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_xoh_index_
 
module procedure get_xoh_index_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_num_xoh_
 
module procedure get_num_xoh_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_name_xoh_
 
module procedure get_name_xoh_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface update_
 
module procedure update_tempdepparam_surftl
module procedure update_xoh_in_cd_surftl
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface assignment(=)
 
module procedure copy_surftl
 
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
subroutine create_surftl &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create triple layer object
!
!   $Arguments:
!
 
type(t_surface_tl), intent(inout) ::  this  ! Type triple layer model. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, parameter            :: &
 r0=0.0d0  
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
this%pp => null ()
this%txoh => null ()
this%stq0 => null ()
this%stqb => null ()
!%------------------------------------------------------------
this%cte_8rtee0=r0
this%cte_f_rt=r0
this%locktxoh=.false.
this%lockstq0=.false.
this%lockstqb=.false.
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_surftl &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy triple layer object
!
!   $Arguments:
!
 
type(t_surface_tl), intent(inout)  :: this 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, parameter            :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
this%pp => null ()
this%txoh => null ()
this%stq0 => null ()
this%stqb => null ()
!%------------------------------------------------------------
this%cte_8rtee0=r0
this%cte_f_rt=r0
this%locktxoh=.false.
this%lockstq0=.false.
this%lockstqb=.false.
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_parent_surftl &
   (this, &
    parent)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
type(t_surface_tl), intent(inout)          :: this

type(t_parentsurface), intent(in), target  :: parent 
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
subroutine compute_sigmad1_surftl &
   (this, &
    sigmad, &
    ionstr, &
    expfid, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute tha charge in the layer 0 acccording
!%   sigma0 = (8*r*t*e*e0*I)1/2 * sinh( f *fid)
!%                                     ---
!%                                     2rt
!%
!%   where
!%    I= Ionic strength of the solution
!%    r= Universal gas constant
!%    f= Faraday constant
!%    fid= Electric potential
!%    e0= Permitivity of free space
!%    e= Dilectric constant of water
!
!   $Arguments:
!
 
type(t_surface_tl),intent(in) :: this

integer, intent(in)           :: numsite

real*8, intent(out)           :: sigmad(numsite)

real*8, intent(in)            :: ionstr

real*8, intent(in)            :: expfid(numsite) 
 
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
 fid(numsite)
real*8          :: &
 cte, &
 cte2 
real*8, parameter :: &
 r2=2.0d0, &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
cte=dsqrt(this%cte_8rtee0*ionstr)
cte2=this%cte_f_rt/r2
!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
fid=r0
fid=-dlog(expfid)/this%cte_f_rt
sigmad=r0
sigmad=cte*sinh(cte2*fid)
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_sigmad_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sigmad2_surftl &
   (this, &
    sigmad, &
    expfib, &
    expfid, &
    capext, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_tl), intent(in)   :: this

integer, intent(in)              :: numsite

real*8, intent(out)              :: sigmad(numsite)
 
real*8, intent(in)               :: expfib(numsite)

real*8, intent(in)               :: expfid(numsite)

real*8, intent(in)               :: capext(numsite) 
 
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
 fib(numsite), &
 fid(numsite)
integer         :: &
 i
real*8          :: &
 c2(numsite)
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
fib=r0
fid=r0
fib=-dlog(expfib)/this%cte_f_rt
fid=-dlog(expfid)/this%cte_f_rt
!%-----------------------------------------------------------
sigmad=r0
sigmad=capext*(fid-fib)
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_sigmad_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dsigmad1_surftl &
   (this, &
    dsigmad, &
    ionstr, &
    expfid, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
type(t_surface_tl), intent(in)                 :: this

integer, intent(in)                            :: numsite

real*8, intent(out)                            :: dsigmad(numsite,numsite*numsk)

real*8, intent(in)                             :: expfid(numsite)

real*8, intent(in)                             :: ionstr
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
 fid(numsite)
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
if (numsite.ne.this%pp%numsite) goto 10
cte=dsqrt(this%cte_8rtee0*ionstr)
cte2=this%cte_f_rt/r2
fid=r0
fid=-dlog(expfid)/this%cte_f_rt
dsigmad=r0
ipri=0
do i=1,this%pp%numsite
 dsigmad(i,ipri+4)=-cte*cosh(cte2*fid(i))/(r2*expfid(i))
 ipri=ipri+numsk
end do
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_sigmad_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dsigmad2_surftl &
   (this, &
    dsigmad, &
    expfib, &
    expfid, &
    capext, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
type(t_surface_tl), intent(in)                 :: this

integer, intent(in)                            :: numsite

real*8, intent(out)                            :: dsigmad(numsite,numsite*numsk)

real*8, intent(in)                             :: expfid(numsite)

real*8, intent(in)                             :: expfib(numsite)

real*8, intent(in)                             :: capext(numsite)
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
 ipri
real*8          :: &
 c2, &
 cte
character(len=100):: &
 nameprop
real*8, parameter :: &
 r0=0.0d0, &
 r1=1.0d0
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
cte=r1/this%cte_f_rt
!%-----------------------------------------------------------
dsigmad=r0
ipri=0
do i=1,this%pp%numsite
 dsigmad(i,ipri+3)=capext(i)*cte/expfib(i)
 dsigmad(i,ipri+4)=-capext(i)*cte/expfid(i)
 ipri=ipri+numsk
end do
 
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_dsigmad_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sigma01_surftl &
   (this, &
    sigma0, &
    cd, &
    spsurfarea, &
    numsite, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_tl), intent(in)           :: this

integer, intent(in)                      :: numsite          ! Number od sites 

real*8, intent(out)                      :: sigma0(numsite)

real*8, intent(in)                       :: spsurfarea(numsite)

real*8, intent(in)                       :: cd(this%pp%numsp)

logical, intent(out)                     :: iserror ! if true, then there was an error 
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
integer         :: &
 i
character(len=100):: &
 msg  
real*8, parameter :: &
 r0=0.0d0
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
iserror=.false. 
msg='' 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
sigma0=r0
sigma0=matmul(this%stq0,cd)
sigma0=(faraday/spsurfarea)*sigma0
!%-----------------------------------------------------------
return
 
10 continue  
print *,'***************************'
print *,'Surface:'
print *,'Name:',this%pp%name
print *,'Service: compute_sigma0_'
print *, msg
print*, '***************************'
iserror=.true.
return 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sigma02_surftl &
   (this, &
    sigma0, &
    expfi0, &
    expfib, &
    capint, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: 
!%   Compute the charge in the layer 0 according triple layer
!%   model
!%
!%   sigma0 = c1 (fi0-fib)
!%
!%   where
!%   c1 = capacitance in [f/m2]
!%   fi0= Potential in layer 0
!%   fib= Potential in layer beta
!
!   $Arguments:
!
type(t_surface_tl), intent(in)                 :: this

integer, intent(in)                            :: numsite

real*8, intent(out)                            :: sigma0(numsite)

real*8, intent(in)                             :: expfi0(numsite)

real*8, intent(in)                             :: expfib(numsite)

real*8, intent(in)                             :: capint(numsite)
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
 fi0(numsite), &
 fib(numsite)
integer         :: &
 i
character(len=100) :: &
 nameprop 
real*8, parameter :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
fi0=r0
fib=r0
fi0=-dlog(expfi0)/this%cte_f_rt
fib=-dlog(expfib)/this%cte_f_rt
!%-----------------------------------------------------------
sigma0=r0
do i=1,this%pp%numsite
 sigma0(i)=capint(i)*(fi0(i)-fib(i))
end do
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_sigmab_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dsigma01_surftl &
   (this, &
    dsigma0, &
    cd, &
    sk, &
    dcd, &
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
 
type(t_surface_tl), intent(in)                 :: this

integer, intent(in)                            :: numsite

integer, intent(in)                            :: numsp

real*8, intent(out)                            :: dsigma0(numsite,numsite*numsk)

real*8, intent(in)                             :: spsurfarea(numsite)

real*8, intent(in)                             :: dcd(numsp,numsite*numsk) 

real*8, intent(in)                             :: cd(numsp)

real*8, intent(in)                             :: sk(numsite*numsk)  
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
   if (ixoh==j) then
    dsigma0(i,ipri+1)=this%stq0(i,j) 
   else 
    dsigma0(i,ipri+1)=this%stq0(i,j) * dcd(j,ipri+1)
   end if
   dsigma0(i,ipri+2)=dsigma0(i,ipri+2)+ this%stq0(i,j) * dcd(j,ipri+2)
   dsigma0(i,ipri+3)=dsigma0(i,ipri+3)+ this%stq0(i,j) * dcd(j,ipri+3)
 end do
 
 dsigma0(i,:)=(faraday/spsurfarea(i))*dsigma0(i,:)
 
 ipri=ipri+numsk
 
end do
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_dsigma0_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dsigma02_surftl &
   (this, &
    dsigma0, &
    expfi0, &
    expfib, &
    capint, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_tl), intent(in)          :: this

integer, intent(in)                     :: numsite

real*8, intent(out)                     :: dsigma0(numsite,numsite*numsk)

real*8, intent(in)                      :: expfi0(numsite)

real*8, intent(in)                      :: expfib(numsite)

real*8, intent(in)                      :: capint(numsite)  
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
real*8          :: &
 cte
integer         :: &
 i, &
 ipri
 real*8, parameter :: &
 r0=0.0d0, &
 r1=1.0d0 
!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
cte=r1/this%cte_f_rt
dsigma0=r0
ipri=0
do i=1,this%pp%numsite
 dsigma0(i,ipri+2)=-capint(i)*cte/expfi0(i)
 dsigma0(i,ipri+3)=capint(i)*cte/expfib(i)
 ipri=ipri+numsk
end do
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_dsigma0_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sigmab1_surftl &
   (this, &
    sigmab, &
    cd, &
    spsurfarea, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_tl), intent(in)           :: this

integer, intent(in)                      :: numsite

real*8, intent(out)                      :: sigmab(numsite)

real*8, intent(in)                       :: cd(this%pp%numsp)

real*8, intent(in)                       :: spsurfarea(numsite)
 
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
 i 
real*8, parameter :: &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
sigmab=r0
sigmab=matmul(this%stqb,cd)
sigmab=(faraday/spsurfarea)*sigmab
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_sigmab_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sigmab2_surftl &
   (this, &
    sigmab, &
    expfi0, &
    expfib, &
    expfid, &
    capint, &
    capext, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
type(t_surface_tl), intent(in)           :: this

integer, intent(in)                      :: numsite

real*8, intent(out)                      :: sigmab(numsite)

real*8, intent(in)                       :: expfi0(numsite)

real*8, intent(in)                       :: expfib(numsite)

real*8, intent(in)                       :: expfid(numsite)

real*8, intent(in)                       :: capint(numsite)

real*8, intent(in)                       :: capext(numsite) 
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
 fi0(numsite), &
 fib(numsite), &
 fid(numsite)
integer         :: &
 i
real*8, parameter :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
sigmab=r0
fi0=r0
fib=r0
fid=r0
fi0=-dlog(expfi0)/this%cte_f_rt
fib=-dlog(expfib)/this%cte_f_rt
fid=-dlog(expfid)/this%cte_f_rt
!%-----------------------------------------------------------
do i=1,this%pp%numsite
 sigmab(i)=capint(i)*(fib(i)-fi0(i))+capext(i)*(fib(i)-fid(i))
end do
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_sigmab_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dsigmab1_surftl &
   (this, &
    dsigmab, &
    cd, &
    sk, &
    dcd, &
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
 
type(t_surface_tl), intent(in)              :: this

integer, intent(in)                         :: numsite

integer, intent(in)                         :: numsp

real*8, intent(out)                         :: dsigmab(numsite,numsite*numsk)
 
real*8, intent(in)                          :: dcd(numsp,numsite*numsk)

real*8, intent(in)                          :: cd(numsp) 

real*8, intent(in)                          :: sk(numsite*numsk)

real*8, intent(in)                          :: spsurfarea(numsite)
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
real*8, parameter :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
dsigmab=r0
ipri=0
do i=1,numsite
 ixoh=this%pp%idxoh(i)
 do j=1,numsp
   if (ixoh==j) then
    dsigmab(i,ipri+1)=this%stqb(i,j) 
   else 
    dsigmab(i,ipri+1)=this%stqb(i,j) * dcd(j,ipri+1)
   end if
   dsigmab(i,ipri+2)=dsigmab(i,ipri+2)+ this%stqb(i,j) * dcd(j,ipri+2)
   dsigmab(i,ipri+3)=dsigmab(i,ipri+3)+ this%stqb(i,j) * dcd(j,ipri+3)
 end do
 
 dsigmab(i,:)=(faraday/spsurfarea(i))*dsigmab(i,:)
 
 ipri=ipri+numsk
 
end do
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_dsigmab_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dsigmab2_surftl &
   (this, &
    dsigmab, &
    expfi0, &
    expfib, &
    expfid, &
    capint, &
    capext, &
    numsite)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_tl), intent(in)  :: this

integer, intent(in)             :: numsite

real*8, intent(out)             :: dsigmab(numsite,numsite*numsk)

real*8, intent(in)              :: expfi0(numsite)

real*8, intent(in)              :: expfib(numsite)

real*8, intent(in)              :: expfid(numsite) 

real*8, intent(in)              :: capint(numsite)

real*8, intent(in)              :: capext(numsite) 
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
 cte
integer         :: &
 i, &
 ipri
real*8, parameter :: &
 r0=0.0d0, &
 r1=1.0d0  
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
dsigmab=r0
cte=r1/this%cte_f_rt
!%-----------------------------------------------------------
ipri=0
do i=1,this%pp%numsite
 dsigmab(i,ipri+2)=capint(i)*cte/expfi0(i)
 dsigmab(i,ipri+3)=-(capint(i)+capext(i))*cte/expfib(i)
 dsigmab(i,ipri+4)=capext(i)*cte/expfid(i)
 ipri=ipri+numsk
end do
!%-----------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_dsigmab_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine compute_txoh_surftl &
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
 
type(t_surface_tl), intent(in)  :: this

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
if (numsite.ne.this%pp%numsite) goto 10
txohc=r0
txohc=matmul(this%txoh,cd)
!%------------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_txoh_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%    Update all parameter in triple layer model that depend
!%    of the temperature
!%    In general here speak about constant
!%    cte1= 8RTEE0
!%    cte2= F/RT
!%
!%************************************************************
subroutine update_tempdepparam_surftl &
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
 
type(t_surface_tl) :: &
 this
real*8              :: &
 temp
logical             :: &
 iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100) :: &
 nameprop, &
 msg
real*8             :: &
 epsiw, &
 tk 
real*8, parameter  :: &
 zerokelvin=273.15d0, &
 r3=1.0d3, &
 r8=8.0d0 
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
 call update_ (this%pp%paqph,temp,iserror)
 if (iserror) return
 call get_prop_ (this%pp%paqph, epsiw, nameprop, iserror)
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
print *,'surface:'
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
subroutine copy_surftl &
   (copied, &
    this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_tl), intent(in)  :: this

type(t_surface_tl), intent(out) :: copied 
 
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
if (this%lockstqb) then
 copied%lockstqb=this%lockstqb
 ndim1=size(this%stqb,1)
 ndim2=size(this%stqb,2)
 call check_pointer_ (copied%stqb,ndim1,ndim2,.true.)
 copied%stqb=this%stqb
end if
!%------------------------------------------------------------
copied%cte_8rtee0=this%cte_8rtee0
copied%cte_f_rt=this%cte_f_rt
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine compute_dtxoh_surftl &
   (this, &
    dtxoh, &
    cd, &
    sk, &
    dcd, &
    numsite, &
    numsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_tl) :: &
 this
integer            :: &
 numsite, &
 numsp
real*8             :: &
 dtxoh(numsite,numsite*numsk), &
 dcd(numsp,numsite*numsk)
real*8             :: &
 cd(numsp), &
 sk(numsite*numsk) 
 
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
 ipri, &
 i, &
 j, &
 ixoh  
real*8, parameter  :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
if (numsite.ne.this%pp%numsite) goto 10
dtxoh=r0
!%------------------------------------------------------------
ipri=0
do i=1,this%pp%numsite
 ixoh=this%pp%idxoh(i)
 do j=1,numsp
  if (ixoh==j) then
   dtxoh(i,ipri+1)=this%txoh(i,j)
  else
   dtxoh(i,ipri+1)=dtxoh(i,ipri+1)+this%txoh(i,j)*dcd(j,ipri+1)
  end if
  dtxoh(i,ipri+2)=dtxoh(i,ipri+2)+this%txoh(i,j)*dcd(j,ipri+2) 
  dtxoh(i,ipri+3)=dtxoh(i,ipri+3)+this%txoh(i,j)*dcd(j,ipri+3) 
 end do
 ipri=ipri+numsk
end do

!%------------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: compute_dtxoh_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine get_num_tot_sk_surftl &
   (this, &
    numskloc)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
type(t_surface_tl) :: &
 this
integer       :: &
 numskloc
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
!%
!%************************************************************
subroutine get_xoh_index_surftl &
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
 
type(t_surface_tl) :: &
 this
integer             :: &
 numxoh
integer, pointer    :: &
 idxoh(:)
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
!%
!%************************************************************
subroutine update_xoh_in_cd_surftl &
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
 
type(t_surface_tl) :: &
 this
integer             :: &
 numsite
real*8              :: &
 sk(this%pp%numsite*numsk), &
 cd(this%pp%numsp)
 
integer             :: &
 i, &
 isk
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
!%
!%************************************************************
subroutine init_surftl &
   (this, &
    sk, &
    cd, &
    txoh)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
type(t_surface_tl) :: &
 this
integer             :: &
 numsite
real*8              :: &
 txoh(this%pp%numsite), &
 cd(this%pp%numsp)
real*8, pointer     :: &
 sk(:)
 
integer             :: &
 isk, &
 i
real*8, parameter   :: &
 r1=1.0d0 
!%------------------------------------------------------------
call check_pointer_ (sk,this%pp%numsite*numsk,.true.)
sk=r1
isk=1
do i=1,this%pp%numsite
 sk(isk)=txoh(i)
 isk=isk+numsk
end do
!%------------------------------------------------------------
call update_ (this, cd, sk)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine get_num_xoh_surftl &
   (this, &
    numxoh)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
type(t_surface_tl) :: &
 this
integer       :: &
 numxoh
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
!%
!%************************************************************
subroutine get_name_xoh_surftl &
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
 
type(t_surface_tl) :: &
 this
integer       :: &
 numxoh
character(len=100), pointer:: &
 namexoh(:)
 
integer                :: &
 i, &
 ipri
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
!%
!%************************************************************
subroutine set_surftl &
   (this, &
    pp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_tl), intent(inout)           :: this

type(t_parentsurface), intent(in), target   :: pp 

logical, intent(out)                        :: iserror
 
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
real*8             :: &
 z
character(len=100) :: &
 nameprop, &
 namesp, &
 msg
logical            :: &
 be 
real*8, parameter  :: &
 r1=1.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%---------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------
this%pp => pp
!%---------------------------------------------------------
if (this%pp%locksp) then
 call check_pointer_ (this%stq0,this%pp%numsite,this%pp%numsp,.true.)
 call check_pointer_ (this%stqb,this%pp%numsite,this%pp%numsp,.true.)
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
  namesp=this%pp%pspecies(isp)%ptr%name
  call get_prop_ (this%pp%pspecies(isp)%ptr,z,nameprop,msg,iserror)
  if (iserror) goto 10  
  call find(namesp,'h2+',be)
!%------------------
  select case (be)
  case (.true.)
   this%stq0(i,isp)=z
  case (.false.)
   call find(namesp,'o-',be)
   select case (be)
   case (.true.)
    this%stq0(i,isp)=z
   case (.false.)
   call find(namesp,'oh',be)
    select case (be)
      case (.false.)
       this%stq0(i,isp)=-z
       this%stqb(i,isp)=z
      end select
   end select
  end select
 
!%------------------
 this%txoh(i,isp)=r1
!%------------------
 end do
end do
!%----------------------------------------------------------
this%lockstq0=.true.
this%lockstqb=.true.
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
subroutine update_sk_surftl &
   (this, &
    sk, &
    cd, &
	dcd, &
    txoh, &
    capint, &
    capext, &
    spsurfarea, &
    ionstr, &
    isconvergence, &
	isupmxiter, &
    iter, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_tl), intent(in) :: this

real*8, intent(out)            :: sk(this%pp%numsite*numsk)

real*8, intent(in)             :: cd(this%pp%numsp)

real*8, intent(in)             :: dcd(this%pp%numsp,this%pp%numsite*numsk)

real*8, intent(in)             :: txoh(this%pp%numsite)

real*8, intent(in)             :: capint(this%pp%numsite)

real*8, intent(in)             :: capext(this%pp%numsite)

real*8, intent(in)             :: spsurfarea(this%pp%numsite)

real*8, intent(in)             :: ionstr

logical, intent(out)           :: isconvergence

logical, intent(out)           :: isupmxiter

integer, intent(in)            :: iter 

logical, intent(out)           :: iserror 
 
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
 expfib(:) => null (), &
 expfid(:) => null (), &
 sigma01(:) => null (), &
 dsigma01(:,:) => null (), &
 sigma02(:) => null (), &
 dsigma02(:,:) => null (), &
 sigmab1(:) => null (), &
 dsigmab1(:,:) => null (), &
 sigmab2(:) => null (), &
 dsigmab2(:,:) => null (), &
 sigmad1(:) => null (), &
 dsigmad1(:,:) => null (), &
 sigmad2(:) => null (), &
 dsigmad2(:,:) => null (), &
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
 maxerror, &
 maxres 
character(len=100)  :: &
 msg
real*8, parameter   :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false. 
msg='' 
isconvergence=.false.
isupmxiter=.false.
!%-----------------------------------------------------------
numtotsk=this%pp%numsite*numsk
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (jacobian,numtotsk,numtotsk,.true.)
call check_pointer_ (residual,numtotsk,.true.)
call check_pointer_ (indx,numtotsk,.true.)
call check_pointer_ (expfi0,this%pp%numsite,.true.)
call check_pointer_ (expfib,this%pp%numsite,.true.)
call check_pointer_ (expfid,this%pp%numsite,.true.)
call check_pointer_ (txohc,this%pp%numsite,.true.)
call check_pointer_ (dtxohc,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (sigma01,this%pp%numsite,.true.)
call check_pointer_ (dsigma01,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (sigma02,this%pp%numsite,.true.)
call check_pointer_ (dsigma02,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (sigmab1,this%pp%numsite,.true.)
call check_pointer_ (dsigmab1,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (sigmab2,this%pp%numsite,.true.)
call check_pointer_ (dsigmab2,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (sigmad1,this%pp%numsite,.true.)
call check_pointer_ (dsigmad1,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (sigmad2,this%pp%numsite,.true.)
call check_pointer_ (dsigmad2,this%pp%numsite,numtotsk,.true.)
call check_pointer_ (error,numtotsk,.true.)
!%-----------------------------------------------------------
isk=0
do i=1,this%pp%numsite
 expfi0(i)=sk(isk+2)
 expfib(i)=sk(isk+3)
 expfid(i)=sk(isk+4)
 isk=isk+numsk
end do
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
call compute_txoh_ (this,txohc,cd,this%pp%numsite)
call compute_dtxoh_(this,dtxohc,cd,sk,dcd,this%pp%numsite,this%pp%numsp)
call compute_sigma0_ (this,sigma01,cd,spsurfarea,this%pp%numsite,iserror)
if (iserror) goto 20 
call compute_dsigma0_ (this,dsigma01,cd,sk,dcd,spsurfarea,this%pp%numsite,this%pp%numsp)
call compute_sigma0_ (this,sigma02,expfi0,expfib,capint,this%pp%numsite)
call compute_dsigma0_(this,dsigma02,expfi0,expfib,capint,this%pp%numsite)
call compute_sigmab_(this,sigmab1,cd,spsurfarea,this%pp%numsite)
call compute_dsigmab_(this,dsigmab1,cd,sk,dcd,spsurfarea,this%pp%numsite,this%pp%numsp)
call compute_sigmab_(this,sigmab2,expfi0,expfib,expfid,capint,capext,this%pp%numsite)
call compute_dsigmab_(this,dsigmab2,expfi0,expfib,expfid,capint,capext,this%pp%numsite)
call compute_sigmad_(this,sigmad1,ionstr,expfid,this%pp%numsite)
call compute_dsigmad_(this,dsigmad1,ionstr,expfid,this%pp%numsite)
call compute_sigmad_(this,sigmad2,expfib,expfid,capext,this%pp%numsite)
call compute_dsigmad_(this,dsigmad2,expfib,expfid,capext,this%pp%numsite)
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
isk=0
do i=1,this%pp%numsite
  residual(isk+1)=txohc(i)-txoh(i)
  jacobian(isk+1,:)=dtxohc(i,:)
  residual(isk+2)=sigma01(i)-sigma02(i)
  jacobian(isk+2,:)=dsigma01(i,:)-dsigma02(i,:)
  residual(isk+3)=sigmab1(i)-sigmab2(i)
  jacobian(isk+3,:)=dsigmab1(i,:)-dsigmab2(i,:)
  residual(isk+4)=sigmad1(i)-sigmad2(i)
  jacobian(isk+4,:)=dsigmad1(i,:)-dsigmad2(i,:)
  isk=isk+numsk
end do
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
residual=-residual
!%-----------------------------------------------------------
!% Compute the maximum error in residual 
!%-----------------------------------------------------------
maxres=maxval(dabs(residual))
!%-----------------------------------------------------------
call ludcmp (jacobian,numtotsk,numtotsk,indx,dd,msg,iserror)
if (iserror) goto 20
call lubksb (jacobian,numtotsk,numtotsk,indx,residual)
error=r0
error=residual/sk
error=dabs(error)
maxerror=maxval(error)
if (maxerror.gt.facmaxads) then
  residual = residual * facmaxads/maxerror
end if

sk = sk + residual
 
isconvergence=(maxres.le.tolresads.and.maxerror.le.tolunkads)

isupmxiter=(iter>maxiterads)
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
call check_pointer_ (expfib,1,.false.)
call check_pointer_ (expfid,1,.false.)
call check_pointer_ (txohc,1,.false.)
call check_pointer_ (dtxohc,1,1,.false.)
call check_pointer_ (sigma01,1,.false.)
call check_pointer_ (dsigma01,1,1,.false.)
call check_pointer_ (sigma02,1,.false.)
call check_pointer_ (dsigma02,1,1,.false.)
call check_pointer_ (sigmab1,1,.false.)
call check_pointer_ (dsigmab1,1,1,.false.)
call check_pointer_ (sigmab2,1,.false.)
call check_pointer_ (dsigmab2,1,1,.false.)
call check_pointer_ (sigmad1,1,.false.)
call check_pointer_ (dsigmad1,1,1,.false.)
call check_pointer_ (sigmad2,1,.false.)
call check_pointer_ (dsigmad2,1,1,.false.)
call check_pointer_ (error,1,.false.)
if (iserror) goto 10 
!%-----------------------------------------------------------
return
 
10 continue  
print *,'***************************'
print *,'Surface:'
print *,'Name:',this%pp%name
print *,'Service: update_sk_'
print *, msg
print*, '***************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_surftl &
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
 
type (t_surface_tl), intent(in)    :: this

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
write (ioutput,*) '           Model: Triple Layer                '
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
!% Write the STQb definition 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
write (ioutput,*) '----------------------------------------------'
write (ioutput,*) '      STQB definition by sites       '
write(ioutput,1) (this%pp%pspecies(j)%ptr%name,j=1,this%pp%numsp)
do i=1,this%pp%numsite
  write(ioutput,2) (this%stqb(i,j),j=1,this%pp%numsp)
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
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
 end module m_surface_tl
