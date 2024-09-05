module m_surface
!-------------------------------------------------------------------------
!
!   $Description: Represent a surface
!
!   $Use:
!
!   $Author: Sergio Andrés Bea Jofré
!
!   $License:
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_parentsurface
use m_surface_cexch
use m_surface_dl
use m_surface_tl
use m_surface_constcap
use m_species
use m_phase
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%------------------------------------------------------
!%------------------------------------------------------
private  ::
!%------------------------------------------------------
!%------------------------------------------------------
public:: &
create_ &                   ! Create the surface object. 
,destroy_ &                 ! Destroy the surface object. 
,init_ &                    ! Init the surface object. 
,update_sk_ &
,get_name_ &
,get_numsites_ &
,get_num_tot_sk_ &
,get_name_model_ &
,get_namesp_ &
,get_idxoh_ &
,get_name_xoh_ &
,get_numsp_ &
,get_num_xoh_ &
,get_if_sp_is_present_ &
,get_if_active_ &
,set_ &
,update_ &                 ! Update parameters that depend of the temperature
,get_pspecies_ &
,write_ &
,assignment(=) &
,change_to_mol_
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
integer, parameter :: &
kd=1, &
langmuir=2, &
freundlich=3, &
cexch=4, &
constantcapacitance=5, &
diffuselayer=6, &
triplelayer=7
!%------------------------------------------------------
!%------------------------------------------------------
!% Type pointer to surafce object 
!%------------------------------------------------------
!%------------------------------------------------------
type, public::t_psurface

type(t_surface), pointer::ptr

end type
!%------------------------------------------------------
!%------------------------------------------------------
!%------------------------------------------------------
!%------------------------------------------------------
!%------------------------------------------------------
type, public::t_surface
 
private                                ::
 
type (t_parentsurface), pointer        :: pp                  ! Parent surface 
 
type (t_surface_cexch), pointer        :: t_psurf_cexch
 
type (t_surface_dl), pointer           :: t_psurf_dl
 
type (t_surface_constcap), pointer     :: t_psurf_constcap
 
type (t_surface_tl), pointer           :: t_psurf_tl
 
integer                                :: itype
 
 
end type t_surface
 
!%----------------------------------------------------
!%----------------------------------------------------
interface create_
 
module procedure create_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface init_
 
module procedure init_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface update_sk_
 
module procedure update_sk_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_name_
 
module procedure get_name_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_pspecies_
 
module procedure get_pspecies_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_numsites_
 
module procedure get_numsites_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_num_tot_sk_
 
module procedure get_num_tot_sk_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_name_model_
 
module procedure get_name_model_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_namesp_
 
module procedure get_namesp_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_idxoh_
 
module procedure get_idxoh_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_name_xoh_
 
module procedure get_name_xoh_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_numsp_
 
module procedure get_numsp_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_num_xoh_
 
module procedure get_num_xoh_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_if_sp_is_present_
 
module procedure get_if_sp_is_present_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_if_active_
 
module procedure get_if_active_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface set_
 
module procedure set_surf
module procedure set_paqph_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface destroy_
 
module procedure destroy_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface update_
 
module procedure update_tempdepparam_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface write_
 
module procedure write_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface assignment(=)
 
module procedure copy_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface change_to_mol_
 
module procedure change_to_mol_surf
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
contains
!%----------------------------------------------------
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine create_surf &
   (this)
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
 
type(t_surface) :: &
 this
!%------------------------------------------------------------
this%itype=0
this%pp => null ()
this%t_psurf_cexch => null ()
this%t_psurf_dl => null ()
this%t_psurf_tl => null ()
this%t_psurf_constcap => null ()
!%------------------------------------------------------------
allocate(this%pp)
call create_ (this%pp)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%    Destroy surface object
!%************************************************************
subroutine destroy_surf &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface), intent(inout)    :: this 
 
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
if (associated(this%pp)) then
 call destroy_ (this%pp)
 deallocate (this%pp)
 this%pp => null ()
end if
!%------------------------------------------------------------
select case (this%itype)
 
case (cexch)
 
 call destroy_ (this%t_psurf_cexch)
 deallocate (this%t_psurf_cexch)
 
case (diffuselayer)
 
 call destroy_ (this%t_psurf_dl)
 deallocate (this%t_psurf_dl)
 
case (constantcapacitance)
 
 call destroy_(this%t_psurf_constcap)
 deallocate (this%t_psurf_constcap)
 
case (triplelayer)
 
 call destroy_(this%t_psurf_tl)
 deallocate (this%t_psurf_tl)
 
end select
!%------------------------------------------------------------
this%itype=0
this%t_psurf_cexch => null ()
this%t_psurf_dl => null ()
this%t_psurf_constcap => null ()
this%t_psurf_tl => null ()
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_idxoh_surf &
   (this, &
    idxoh, &
    nsite, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set local index of sorption primary species.
!
!   $Arguments:
!
 
type(t_surface), intent(in)           :: this

integer, intent(out)                  :: nsite

integer, pointer, dimension(:)        :: idxoh 

logical, intent(out)                  :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
!
!   $code
!
call get_idxoh_ (this%pp,idxoh,nsite,iserror)
!%-----------------------------------------------------------
return
 

end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_num_tot_sk_surf &
   (this, &
    numsk)
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
 
type(t_surface), intent(in)     :: this

integer, intent(out)            :: numsk
!%------------------------------------------------------------
select case (this%itype)
 
case (cexch)
 
 call get_num_sk_ (this%t_psurf_cexch,numsk)
 
case (diffuselayer)
 
 call get_num_sk_ (this%t_psurf_dl,numsk)
 
case (constantcapacitance)
 
 call get_num_sk_ (this%t_psurf_constcap,numsk)
 
case (triplelayer)
 
 call get_num_tot_sk_ (this%t_psurf_tl,numsk)
 
end select
 
 
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
subroutine get_pspecies_surf &
   (this, &
    pspecies, &
    ithsp, &
    iserror)
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
 
type(t_surface), intent(in) :: &
 this
type(t_species), pointer:: &
 pspecies
integer, intent(in)        :: &
 ithsp
logical, intent(out)       :: &
 iserror
!%------------------------------------------------------------
call get_pspecies_ (this%pp, pspecies, ithsp,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_surf &
   (this, &
    name, &
    namemodel, &
    species, &
    numsp, &
    numsite, &
    numspsite, &
    namepropsite, &
    propsite, &
    numpropsite, &
    idxoh, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface), intent(inout)             :: this

integer, intent(in)                        :: numsp

integer, intent(in)                        :: numsite

integer, intent(in)                        :: numpropsite

character(len=*), intent(in)               :: namemodel

integer, intent(in)                        :: numspsite(numsite)

character(len=*), intent(in)               :: name

type(t_species), intent(in)                :: species(numsp)

real*8, intent(in)                         :: propsite(numsite,numpropsite)

character(len=*), intent(in)               :: namepropsite(numpropsite)

integer, intent(in)                        :: idxoh(numsite)

logical, intent(out)                       :: iserror 
 
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
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
select case (namemodel)
case ('cexch','ex','EX')
 this%itype = cexch
case ('dl','DL')
 this%itype = diffuselayer
case ('cc','CC')
 this%itype = constantcapacitance
case ('tl','TL')
 this%itype = triplelayer
case default
 msg='Error, not recognized the sorption model:'
 call add_ (msg,namemodel)
 goto 10
end select

!%-------------------------------------------------------------
!% Set in the parent surface object 
!%-------------------------------------------------------------
call set_ (this%pp,name)
call set_ &
   (this%pp, &
    species, &
    numsp, &
    numspsite, &
    numsite, &
    msg, &
    iserror)
if (iserror) goto 10
call set_ (this%pp,idxoh,numsite,msg,iserror)
if (iserror) goto 10
!%------------------------------------------------------------
call set_(this%pp,namepropsite,propsite,numsite,numpropsite,msg,iserror)
if (iserror) goto 10
!%-------------------------------------------------------------
select case (this%itype)
!%------------
case (cexch)
  allocate (this%t_psurf_cexch)
  call create_ (this%t_psurf_cexch)
  call set_ &
   (this%t_psurf_cexch, &
    this%pp, &
    iserror)
!%------------
case (diffuselayer)
  allocate (this%t_psurf_dl)
  call create_ (this%t_psurf_dl)
  call set_ (this%t_psurf_dl,this%pp,iserror)
!%------------
case (constantcapacitance)
  allocate (this%t_psurf_constcap)
  call create_ (this%t_psurf_constcap)
  call set_ (this%t_psurf_constcap,this%pp,iserror)
!%------------
case (triplelayer)
  allocate (this%t_psurf_tl)
  call create_ (this%t_psurf_tl)
  call set_(this%t_psurf_tl,this%pp,iserror)
end select
!%-------------------------------------------------------------
if (iserror) goto 10
!%-------------------------------------------------------------
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
!%
!%************************************************************
subroutine set_paqph_surf &
   (this, &
    paqphase, &
    iserror)
 
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
 
type(t_surface) :: &
 this
type(t_phase)      :: &
 paqphase
logical            :: &
 iserror
!%-------------------------------------------------------------
call set_ (this%pp,paqphase,iserror)
!%-------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine update_tempdepparam_surf &
   (this, &
    temp, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update parameters that depend of the temperature
!
!   $Arguments:
!
 
type(t_surface), intent(inout) :: this

real*8, intent(in)             :: temp

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
character(len=100) :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
select case (this%itype)
 
case (cexch)
  call update_ (this%t_psurf_cexch,temp,iserror)
case (constantcapacitance)
  call update_ (this%t_psurf_constcap,temp,iserror)
case (diffuselayer)
  call update_ (this%t_psurf_dl,temp,iserror)
case (triplelayer)
  call update_ (this%t_psurf_tl,temp,iserror)
case default
  msg='Error, public service not implemented'
  goto 10
end select
!%-----------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Surface:'
print *,'Name:',this%pp%name
print *,'Service: update_'
print *, msg
 print *,'************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_name_surf &
   (this, &
    name)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface), intent(in)   :: this

character(len=*), intent(out) :: name 
 
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
call get_name_ (this%pp,name)
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
subroutine get_numsites_surf &
   (this, &
    numsite)
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
 
type(t_surface) :: &
 this
integer            :: &
 numsite
!%------------------------------------------------------------
call get_numsites_ (this%pp,numsite)
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
subroutine get_name_model_surf &
   (this, &
    namemodel)
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
 
type(t_surface), intent(in)     :: &
 this
character(len=*), intent(out) :: &
 namemodel
!%------------------------------------------------------------
namemodel=' '
!%------------------------------------------------------------
select case (this%itype)
 
case (cexch)
  namemodel='cationexchange'
case (triplelayer)
  namemodel='triplelayer'
case (constantcapacitance)
  namemodel='constantcapacitance'
case (diffuselayer)
  namemodel='diffuselayer'
case (langmuir)
  namemodel='langmuir'
end select
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
subroutine get_namesp_surf &
   (this, &
    namesp)
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
 
type(t_surface) :: &
 this
 
character(len=100), pointer       :: &
 namesp(:)
integer                           :: &
 i
!%------------------------------------------------------------
call get_namesp_ &
   (this%pp, &
    namesp)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
subroutine get_name_xoh_surf &
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
 
type(t_surface) :: &
 this
character(len=100), pointer :: &
 namexoh(:)
integer                     :: &
 numxoh
!%-----------------------------------------------------------
select case (this%itype)
 
case (cexch)
 call get_name_xna_ (this%t_psurf_cexch,namexoh,numxoh)
case (constantcapacitance)
 call get_name_xoh_ (this%t_psurf_constcap,namexoh,numxoh)
case (diffuselayer)
 call get_name_xoh_ (this%t_psurf_dl,namexoh,numxoh)
case (triplelayer)
 call get_name_xoh_ (this%t_psurf_tl,namexoh,numxoh)
end select
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_numsp_surf &
   (this, &
    numsp)
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
 
type(t_surface) :: &
 this
integer            :: &
 numsp
 
call get_numsp_(this%pp,numsp)
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_num_xoh_surf &
   (this, &
    numxoh)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the number of primary sorption species
!
!   $Arguments:
!
 
type(t_surface), intent(in)  :: this

integer, intent(out)         :: numxoh 
 
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
select case (this%itype)
 
case (cexch)
  call get_num_xna_(this%t_psurf_cexch,numxoh)
case (constantcapacitance)
  call get_num_xoh_(this%t_psurf_constcap,numxoh)
case (diffuselayer)
  call get_num_xoh_(this%t_psurf_dl,numxoh)
case (triplelayer)
  call get_num_xoh_(this%t_psurf_tl,numxoh)
end select
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%    Change from equivalent fraction to molality (only
!%    for cation exchange model because now work in
!%    equivalente fraction
!%************************************************************
subroutine change_to_mol_surf &
   (this, &
    beta, &
    txoh, &
    iserror)
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
 
type (t_surface), intent(in)       :: &
 this
real*8                             :: &
 beta(this%pp%numsp)
real*8, intent(in)                 :: &
 txoh(this%pp%numsite)
logical                 :: &
 iserror
!%-----------------------------------------------------------
select case (this%itype)
 
case (cexch)
 
  call change_to_mol_ &
   (this%t_psurf_cexch, &
    beta, &
    txoh, &
    iserror)
 
end select
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_sk_surf &
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
 
type(t_surface), intent(in)              :: this

real*8, intent(out), dimension(:)        :: sk

real*8, intent(inout), dimension(:)      :: cd

real*8, intent(in), dimension(:,:)       :: dcd

real*8, intent(in), dimension(:)         :: txoh

real*8, intent(in), dimension(:)         :: capint

real*8, intent(in), dimension(:)         :: capext

real*8, intent(in), dimension(:)         :: spsurfarea

real*8, intent(in)                       :: ionstr          ! ionic strength 

integer, intent(in)                      :: iter            ! Number of iterations 

logical, intent(out)                     :: isconvergence

logical, intent(out)                     :: isupmxiter

logical, intent(out)                     :: iserror         ! if true, then there was an error
 
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
 

!%-------------------------------------------------------------
!% Select the specialization 
!%-------------------------------------------------------------
select case (this%itype)
 
case (cexch)
 
 call update_sk_(this%t_psurf_cexch,sk,cd,dcd,txoh,isconvergence,isupmxiter,iter,iserror)
 
case (constantcapacitance)
 
 call update_sk_(this%t_psurf_constcap,sk,cd,dcd,txoh,capint,spsurfarea,isconvergence,isupmxiter,iter)
 
case (diffuselayer)
 
 call update_sk_(this%t_psurf_dl,sk,cd,dcd,txoh,spsurfarea,ionstr,isconvergence,isupmxiter,iter)
 
case (triplelayer)
 
 call update_sk_(this%t_psurf_tl,sk,cd,dcd,txoh,capint,capext,spsurfarea,ionstr,isconvergence,isupmxiter,iter,iserror)
 
end select
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_if_sp_is_present_surf &
   (this, &
    namesp, &
    be)
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
 
type(t_surface) :: &
 this
character(len=*)   :: &
 namesp
logical            :: &
 be
!%------------------------------------------------------------
call get_if_sp_is_present_ (this%pp, namesp, be)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_surf &
  (copied, &
   this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface), intent(in)  :: this

type(t_surface), intent(out) :: copied 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer        :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 


!%------------------------------------------------------------
copied%itype=this%itype
if (associated(copied%pp)) then
 call destroy_ (copied%pp)
 deallocate (copied%pp)
end if
allocate (copied%pp)
call create_ (copied%pp)
copied%pp=this%pp
!%------------------------------------------------------------
select case (copied%itype)
 
case (cexch)
 
 allocate(copied%t_psurf_cexch)
 call create_ (copied%t_psurf_cexch)
 copied%t_psurf_cexch=this%t_psurf_cexch
 call set_parent_ (copied%t_psurf_cexch,copied%pp)
 
case (diffuselayer)
 
 allocate(copied%t_psurf_dl)
 call create_ (copied%t_psurf_dl)
 copied%t_psurf_dl=this%t_psurf_dl
 call set_parent_ (copied%t_psurf_dl,copied%pp)
 
case (constantcapacitance)
 
 allocate(copied%t_psurf_constcap)
 call create_ (copied%t_psurf_constcap)
 copied%t_psurf_constcap=this%t_psurf_constcap
 call set_parent_ (copied%t_psurf_constcap,copied%pp) 
 
case (triplelayer)
 
 allocate(copied%t_psurf_tl)
 call create_ (copied%t_psurf_tl)
 copied%t_psurf_tl=this%t_psurf_tl
 call set_parent_ (copied%t_psurf_tl,copied%pp) 
 
end select
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine init_surf &
  (this, &
   sk, &
   cd, &
   txoh, &
   capint, &
   spsurfarea, &
   ionstr)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface), intent(inout)     :: this

real*8, intent(in)                 :: txoh(this%pp%numsite)

real*8, intent(in)                 :: capint(this%pp%numsite)

real*8, intent(in)                 :: spsurfarea(this%pp%numsite)

real*8, intent(inout)              :: cd(this%pp%numsp)

real*8, intent(in)                 :: ionstr                        ! Ionic strength 

real*8, pointer                    :: sk(:)  
 
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
select case (this%itype)
 
case (cexch)
 
 call init_ (this%t_psurf_cexch, sk, cd, txoh)
 
case (triplelayer)
 
 call init_ (this%t_psurf_tl, sk, cd, txoh)
 
case (constantcapacitance)
 
 call init_ (this%t_psurf_constcap,sk,cd,txoh,capint,spsurfarea)
 
case (diffuselayer)
 
 call init_ (this%t_psurf_dl,sk,cd,txoh,spsurfarea,ionstr)
 
end select
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%   Check if all sites is active
!%************************************************************
subroutine get_if_active_surf &
  (this, &
   txoh, &
   beactive)
 
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
 
type(t_surface), intent(in)   :: &
 this
real*8, intent(in)            :: &
 txoh(this%pp%numsite)
logical, intent(out)          :: &
 beactive
!%------------------------------------------------------------
select case (this%itype)
 
case (cexch,triplelayer,langmuir,diffuselayer,constantcapacitance)
 
 call get_if_active_(this%pp,txoh,beactive)
 
end select
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_surf &
  (this, &
   ioutput, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface), intent(in)     :: this

integer, intent(in)             :: ioutput

logical, intent(out)            :: iserror 
 
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
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
select case (this%itype)
 
case (cexch)
 
 call write_ (this%t_psurf_cexch,ioutput,msg,iserror)
 if (iserror) goto 10
 
case (triplelayer)
 
 call write_ (this%t_psurf_tl,ioutput,msg,iserror)
 if (iserror) goto 10
 
case (constantcapacitance)
 
 call write_ (this%t_psurf_constcap,ioutput,msg,iserror)
 if (iserror) goto 10
 
case (diffuselayer)
 
 call write_ (this%t_psurf_dl,ioutput,msg,iserror)
 if (iserror) goto 10
 
end select
!%-----------------------------------------------------------
return
10 continue 
print *,'*******************'
print *,'Surface:           '
print *,'Name:',this%pp%name
print *,'Service: write_'
print *, msg
print *,'*******************'
iserror=.true.
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
end module m_surface
