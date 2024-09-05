module m_parentsurface
!-------------------------------------------------------------------------
!
!   $Description:
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
private                     ::
!%------------------------------------------------------
!%------------------------------------------------------
public                      :: &
create_ &                       ! Create a parent surface object. 
,set_ &                         ! Destroy a parent surface object. 
,destroy_ &
,get_name_ &
,get_numsites_ &
,get_namesp_ &
,get_idxoh_ &
,get_prop_ &
,get_numsp_ &
,get_if_sp_is_present_ &
,assignment(=) &
,update_ &
,get_pspecies_ &
,get_if_active_ &
,write_
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
real*8, public, parameter     :: &
facmaxads=0.5d0, &        ! Factor for to scale the Newton-Raphson solution 
tolresads=1.0d-10, &      ! Tolerance in residual
tolunkads=1.0d-5, &       ! Tolerance in unknowns 
maxiterads=100            ! Maximum number of iterations 
!%------------------------------------------------------
!%------------------------------------------------------
type, public::t_parentsurface
 
character(len=100)                       :: name          ! Name of interphase
 
type(t_pspecies), pointer, dimension(:)  :: pspecies      ! Pointer to chemical species
 
type(t_phase), pointer                   :: paqph
 
real*8, pointer, dimension(:,:)          :: propsite      ! Value of the properties of the sites
 
character(len=100), pointer, dimension(:):: namepropsite  ! Name of sites properties
 
integer                                  :: numsite       ! Number of site of the interphase

integer                                  :: numsp         ! Total number of species of the interp

integer                                  :: numpropsite   ! Number of properties of sites
 
integer, pointer, dimension(:)           :: numspsite     ! Number of species for site [nsite]

integer, pointer, dimension(:)           :: idxoh

logical                                  :: locksp        ! .true. if species vector was allocate

logical                                  :: lockpaqph     ! .true. if aqueous phase class was

logical                                  :: lockpropsite

logical                                  :: lockidxoh
 
end type t_parentsurface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
 
module procedure create_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
 
module procedure set_name_psurf
module procedure set_species_psurf
module procedure set_paqph_psurf
module procedure set_propsite_psurf
module procedure set_idxoh_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
 
module procedure destroy_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_name_
 
module procedure get_name_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_numsites_
 
module procedure get_numsites_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_pspecies_
 
module procedure get_pspecies_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_namesp_
 
module procedure get_namesp_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_prop_
 
module procedure get_prop_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_numsp_
 
module procedure get_numsp_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_if_sp_is_present_
 
module procedure get_if_sp_is_present_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_if_active_
 
module procedure get_if_active_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_idxoh_
 
module procedure get_idxoh_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_
 
module procedure update_temp_param_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
 
module procedure write_psurf
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment(=)
 
module procedure copy_psurf
 
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
subroutine create_psurf &
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
 
type(t_parentsurface) :: &
 this
!%------------------------------------------------------------
this%numsp = 0
this%name = ''
this%numsite=0
this%numsp=0
this%numpropsite=0
this%locksp=.false.
this%lockpaqph=.false.
this%lockpropsite=.false.
this%lockidxoh=.false.
this%pspecies => null ()
this%paqph => null ()
this%numspsite => null ()
this%namepropsite => null ()
this%propsite => null ()
this%idxoh => null ()
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_psurf &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy parent surface
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
 
type(t_parentsurface) :: &
 this
 
integer          :: &
 i
!%------------------------------------------------------------
this%name = ''
call check_pointer_ (this%propsite,1,1,.false.)
call check_pointer_ (this%namepropsite,1,.false.)
call check_pointer_ (this%numspsite,1,.false.)
call check_pointer_ (this%idxoh,1,.false.)
!%------------------------------------------------------------
if (this%numsp>0) then
 do i=1,this%numsp
    call destroy_ (this%pspecies(i)%ptr)
    deallocate (this%pspecies(i)%ptr)
    this%pspecies(i)%ptr => null ()
 end do
 deallocate (this%pspecies)
 this%pspecies => null ()
end if
!%------------------------------------------------------------
this%paqph => null ()
!%------------------------------------------------------------
this%numsp=0
this%locksp=.false.
this%lockpaqph=.false.
this%lockpropsite=.false.
this%lockidxoh=.false.
!%------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_psurf &
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
 
type(t_parentsurface), intent(inout) :: this

real*8, intent(in)                   :: temp

logical, intent(out)                 :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)  :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call update_ (this%paqph,temp,iserror)
if (iserror) then
 msg='Error when calling update_'
 goto 10
end if
!%---------------- --------------------------------------------
return
 
10 continue 
print *,'******************'
print *,'surface:'
print *,'Name:',this%name
print *,'Service: update_'
print *, msg
 print *,'******************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_name_psurf &
   (this, &
    name)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_parentsurface), intent(inout) :: this

character(len=*), intent(in)         :: name 
 
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
 

!%-----------------------------------------------------------
this%name=name
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_species_psurf &
   (this,      & ! i/o
    species,    & ! i
    numsp,     & ! i
    numspsite, & ! i
    numsite, &
    msg, &
    iserror)       ! i
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_parentsurface), intent(inout)     :: this

integer, intent(in)                      :: numsp

integer, intent(in)                      :: numsite

integer, intent(in), dimension(numsite)  :: numspsite

type(t_species), intent(in), dimension(:):: species

character(len=*), intent(out)            :: msg

logical, intent(out)                     :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                   :: &
 i
logical                   :: &
 isrepeated
character(len=100)        :: &
 name 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%----------------------------------------------------------
msg=''
iserror=.false.
!%----------------------------------------------------------
this%numsp=numsp
this%numsite=numsite
!%----------------------------------------------------------
!% Check if the total species in the sites is equal to total
!% number of species 
!%----------------------------------------------------------
if (numsp/=sum(numspsite)) then
 msg='Error, the total number of species is different in the sites'
 goto 10
end if
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
call check_pointer_ (this%numspsite,this%numsite,.true.)
this%numspsite=numspsite
!%----------------------------------------------------------
if (numsp==0) return
!%----------------------------------------------------------
!% Find repeated species
!%----------------------------------------------------------
call find_repeated_ &
    (species(1:numsp)%name, &
     isrepeated, &
     numsp, &
     name)
if (isrepeated) then
    msg='Error, species repeated:'
    call add_ (msg,name)
    goto 10
end if
!%----------------------------------------------------------
allocate (this%pspecies(numsp))
do i=1,this%numsp
 allocate (this%pspecies(i)%ptr)
 call create_ (this%pspecies(i)%ptr)
 this%pspecies(i)%ptr = species(i)
end do
this%locksp=.true.
!%----------------------------------------------------------
return
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_paqph_psurf &
   (this, &
    paqph, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set pointer to aqueous phase
!
!   $Arguments:
!
 
type(t_parentsurface)           :: &
 this
type(t_phase), target         :: &
 paqph
logical                       :: &
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
logical                       :: &
 be
character(len=100)            :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
call get_if_aqueous_ (paqph,be)
if (.not.be) then
 msg='Error, the associated phase is not aqueous'
 goto 10
end if
!%-----------------------------------------------------------
this%paqph => paqph
this%lockpaqph=.true.
!%-----------------------------------------------------------
return
 
10 continue 
print *,'********************'
print *,'surface:'
print *,'Name:',this%name
print *,'Service: set_paqph_'
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
subroutine set_propsite_psurf &
   (this, &
    namepropsite, &
    propsite, &
    numsite, &
    numpropsite, &
    msg, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set surface attributes. 
!
!   $Arguments:
!
 
type(t_parentsurface), intent(inout)  :: this

integer, intent(in)                   :: numsite

integer, intent(in)                   :: numpropsite

character(len=*), intent(in)          :: namepropsite(numpropsite)

character(len=*), intent(out)         :: msg

real*8, intent(in)                    :: propsite(numsite,numpropsite)

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
 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
if (this%lockpropsite) return
if (.not.this%locksp) then
 msg='Error, not associated species in the surface'
 goto 10
end if
if (numsite.ne.this%numsite) then
 msg='Error in number of sites'
 goto 10
end if
!%-----------------------------------------------------------
this%numpropsite=numpropsite
call check_pointer_ (this%namepropsite,numpropsite,.true.)
this%namepropsite=namepropsite
!%-----------------------------------------------------------
call check_pointer_ (this%propsite,this%numsite,this%numpropsite,.true.)
this%propsite=propsite
!%-----------------------------------------------------------
if (this%numsite>0.and.this%numpropsite>0) then
 this%lockpropsite=.true.
end if
!%-----------------------------------------------------------
return
 
10 iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_idxoh_psurf &
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
 
type(t_parentsurface), intent(in)     :: this

integer, intent(out)                  :: nsite

integer, pointer                      :: idxoh(:) 

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
character(len=100)     :: & 
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
nsite=this%numsite
!%-----------------------------------------------------------
call check_pointer_ (idxoh,nsite,.true.)
idxoh=this%idxoh
!%-----------------------------------------------------------
return
 
10 iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_idxoh_psurf &
   (this, &
    idxoh, &
    numsite, &
    msg, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set local index of sorption primary species.
!
!   $Arguments:
!
 
type(t_parentsurface), intent(inout)  :: this

integer, intent(in)                   :: numsite

integer, intent(in)                   :: idxoh(numsite)

character(len=*), intent(out)         :: msg

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
 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
if (this%lockpropsite) return
if (.not.this%locksp) then
 msg='Error, not associated species in the surface'
 goto 10
end if
if (numsite.ne.this%numsite) then
 msg='Error in number of sites'
 goto 10
end if
!%-----------------------------------------------------------
call check_pointer_ (this%idxoh,numsite,.true.)
this%idxoh=idxoh
!%-----------------------------------------------------------
this%lockidxoh=.true.
!%-----------------------------------------------------------
return
 
10 iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_name_psurf &
   (this, &
    name)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_parentsurface), intent(in) :: this

character(len=*), intent(out)     :: name 
 
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
 

 
name=this%name
 
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_pspecies_psurf &
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
 
type(t_parentsurface), intent(in)    :: this

integer, intent(in)                  :: ithsp

type(t_species), pointer             :: pspecies 

logical, intent(out)                 :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                  :: &
msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
if (ithsp<=0.or.ithsp>this%numsp) then
 msg='Error in species index, index='
 call add_ (msg,ithsp)
 goto 10
end if  
!%------------------------------------------------------------
pspecies => this%pspecies(ithsp)%ptr
!%------------------------------------------------------------
return
10 continue 	
print *,'**********************'
print *,'Surface:'
print *,'Name:', this%name
print *,'Service: get_pspecies_'
print *, msg
print *,'**********************'
iserror=.true. 
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine get_numsites_psurf &
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
 
type(t_parentsurface) :: &
 this
integer          :: &
 numsite
!%----------------------------------------------------------
numsite=this%numsite
!%----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%   Get the property of interphase according to name of
!%   property
!%
!%************************************************************
subroutine get_prop_psurf &
   (this, &
    value, &
    nameprop, &
    isite)
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
 
type(t_parentsurface) :: &
 this
integer          :: &
 isite
real*8           :: &
 value
character(len=*) :: &
 nameprop
 
logical          :: &
 be
integer          :: &
 i
!%------------------------------------------------------------
if (isite>this%numsite) goto 10
!%------------------------------------------------------------
be=.false.
value=0.0d0
do i=1,this%numpropsite
 if (nameprop.eq.this%namepropsite(i)) then
   value=this%propsite(isite,i)
   be=.true.
   exit
 end if
end do
if (.not.be) goto 10
!%------------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: get_prop_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%   Check if all sites is active
!%************************************************************
subroutine get_if_active_psurf &
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
 
type(t_parentsurface), intent(in)   :: &
 this
real*8, intent(in)                  :: &
 txoh(this%numsite)
logical, intent(out)                :: &
 beactive
 
integer                             :: &
 i
!%------------------------------------------------------------
beactive=.true.
!%------------------------------------------------------------
do i=1,this%numsite
  if (txoh(i)==0.0d0) then
   beactive=.false.
   exit
  end if
end do
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine get_namesp_psurf &
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
 
type(t_parentsurface) :: &
 this
 
character(len=100), pointer       :: &
 namesp(:)
integer                           :: &
 i
 
if (associated(namesp)) deallocate (namesp)
allocate(namesp(this%numsp))
 
do i=1,this%numsp
 namesp(i) = this%pspecies(i)%ptr%name
end do
 
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine get_numsp_psurf &
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
 
type(t_parentsurface) :: &
 this
integer          :: &
 numsp
 
numsp=this%numsp
 
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
Subroutine get_if_sp_is_present_psurf &
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
 
type(t_parentsurface) :: &
 this
character(len=*)    :: &
 namesp
logical             :: &
 be
 
integer             :: &
 i
character(len=100)  :: &
 namesploc
!%----------------------------------------------------------
if (.not.this%locksp) goto 10
be = .false.
do i=1,this%numsp
 namesploc=this%pspecies(i)%ptr%name
 if (namesploc.eq.namesp) then
 
   be = .true.
 
   exit
 
 end if
 
end do
!%---------------------------------------------------------
return
 
10 print *,'surface:'
print *,'Service: get_if_sp_is_present_'
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
subroutine write_psurf &
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
 
type(t_parentsurface)    :: &
 this
integer             :: &
 ioutput
logical             :: &
 iserror
 
integer             :: &
 i
character(len=*)    :: &
 msg
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
write(ioutput,*) '-------------------------------------'
write(ioutput,*) 'Surface information'
write(ioutput,*) 'Name:',this%name
write(ioutput,*) 'Numsp:',this%numsp
write(ioutput,*) '-------------------------------------'
do i=1,this%numsp
 call write_ (this%pspecies(i)%ptr,ioutput,iserror)
end do
write(ioutput,*) '-------------------------------------'
!%------------------------------------------------------------
return
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine copy_psurf &
  (copied, &
   this)
 
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
 
type(t_parentsurface), intent(in) :: &
 this
type(t_parentsurface), intent(out) :: &
 copied
integer        :: &
 i
!%------------------------------------------------------------
copied%name=this%name
!%------------------------------------------------------------
if (this%locksp) then
 copied%numsite=this%numsite
 copied%numsp=this%numsp
 copied%locksp=this%locksp
 allocate (copied%numspsite(copied%numsite))
 copied%numspsite=this%numspsite
 allocate (copied%pspecies(copied%numsp))
 do i=1,this%numsp
  allocate (copied%pspecies(i)%ptr)
  call create_ (copied%pspecies(i)%ptr)
  copied%pspecies(i)%ptr=this%pspecies(i)%ptr
 end do
end if
!%------------------------------------------------------------
if (this%lockpaqph) then
 copied%lockpaqph=this%lockpaqph
 copied%paqph => this%paqph
end if
!%------------------------------------------------------------
if (this%lockpropsite) then
 copied%lockpropsite=this%lockpropsite
 copied%numpropsite=this%numpropsite
 allocate(copied%namepropsite(copied%numpropsite))
 allocate(copied%propsite(copied%numpropsite,copied%numsite))
 copied%namepropsite=this%namepropsite
 copied%propsite=this%propsite
end if
!%-----------------------------------------------------------
if (this%lockidxoh) then
 copied%lockidxoh=this%lockidxoh
 allocate(copied%idxoh(copied%numsite))
 copied%idxoh=this%idxoh
end if
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_parentsurface
