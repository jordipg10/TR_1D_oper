module m_aqueousphaseideal
!-------------------------------------------------------------------------
!
!   $Description:Aqueous phase with ideal model (concentration=activity)
!
!   $Use: 
! m_parentaqueousphase.
! m_species
! m_general_tools_cheproo
! m_constants
!
!   $Author: Sergio Andr�s Bea Jofr� 
!
!   $License: UPC-CSIC (2008)
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_parentaqueousphase
use m_species
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%------------------------------------------------------------
!%------------------------------------------------------------
private::
!%------------------------------------------------------------
!%------------------------------------------------------------
public:: &
create_ &                  ! Create ideal object 
,destroy_ &                ! Destroy the object 
,set_ &                    ! Set the object  
,set_pparent_ &            ! Set the parent phase in the phase object 
,compute_act_coeff_ &      ! Compute activity coefficients
,compute_dact_coeff_ &     ! Compute derivatives of the activity coefficients
,update_ &                 ! Update parameters that are depending of the temperature 
,write_                    ! Write in ascii the attributes ancapsulated in the object 
!%------------------------------------------------------------
!%------------------------------------------------------------
!% Type definition 
!%------------------------------------------------------------
!%------------------------------------------------------------
Type, public::t_aqueousphaseideal
 
private                                ::
 
type (t_parentaqueousphase), pointer   :: pp   ! Parent aqueous phase pointer
 
 
end type t_aqueousphaseideal
!%------------------------------------------------------------
!%------------------------------------------------------------
!% Interfaces 
!%------------------------------------------------------------
!%------------------------------------------------------------
interface create_
 
module procedure create_aqphideal
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface destroy_
 
module procedure destroy_aqphideal
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface set_
 
module procedure set_aqphideal
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface set_pparent_
 
module procedure set_pparent_aqphideal
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_act_coeff_
 
module procedure compute_act_coeff_aqphideal
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_dact_coeff_
 
module procedure compute_dact_coeff_aqphideal
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface update_
 
module procedure update_temp_param_aqphideal
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface write_
 
module procedure write_aqphideal
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
Contains
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_aqphideal &
(this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create aqueous phase object with ideal thermodynamic model
!
!   $Arguments:
!
 
type(t_aqueousphaseideal), intent(inout):: this ! Ideal aqueous phase type
 
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
this%pp => null ()
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_aqphideal &
(this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy aqueous phase object with ideal thermodynamic model
!
!   $Arguments:
!
 
type(t_aqueousphaseideal), intent(inout):: this ! Ideal aqueous phase type
 
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
this%pp => null ()
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_aqphideal &
(this, &
pp, &
iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set attributes in the object 
!
!   $Arguments:
!
 
type(t_aqueousphaseideal), intent(inout):: this     ! Ideal aqueous phase type

type(t_parentaqueousphase), target      :: pp       ! Parent aqueous phase

logical, intent(out)                    :: iserror  ! iserror=true, then there was an error
 
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
iserror=.false. 
msg=''
!%-----------------------------------------------------------
this%pp => pp
!%------------------------------------------------------------
!% Update the parent phase according to reference temperature 
!%------------------------------------------------------------
call update_ (this,this%pp%pp%tempref,iserror)
if (iserror) goto 10
!%-----------------------------------------------------------
return
10 continue 
print *,'********************************'
print *,'Phase:'
print *,'Name:',this%pp%pp%name
print *,'Service: set_' 
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
subroutine set_pparent_aqphideal &
(this, &
pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_aqueousphaseideal), intent(inout)            :: this ! Ideal aqueous phase type

type(t_parentaqueousphase), intent(in), target      :: pp   ! Parent aqueous phase object 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)  ::  &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
!%-----------------------------------------------------------
this%pp => pp
!%-----------------------------------------------------------
return
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:',this%pp%pp%name
print *,'Service: specia_from_cpri_'
print *, msg
print *,'********************************'
stop 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_aqphideal &
(this, &
temp, &
iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update in the object parameters that depend of the 
! temperature 
!
!   $Arguments:
!
 
type(t_aqueousphaseideal), intent(inout) :: this      ! Ideal aqueous phase type

real*8, intent(in)                       :: temp      ! Temperature [C]

logical, intent(out)                     :: iserror   ! iserror=true, then there was an error
 
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
iserror=.false.
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_aqphideal &
(this, &
ioutput, &
iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii the attributes ancapsulated in the object 
!
!   $Arguments:
!
 
type (t_aqueousphaseideal), intent(in) :: this     ! Ideal aqueous phase type

integer, intent(in)                    :: ioutput  ! Output unit 

logical, intent(out)                   :: iserror  ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)      :: &
msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
write (ioutput,*) '************************************'
write (ioutput,*) '          Aqueous Phase             '
write (ioutput,*) '      Activity model: ideal         '
!%-----------------------------------------------------------
!% Write the aqueous parent phase
!%-----------------------------------------------------------
call write_ (this%pp,ioutput,iserror)
!%-----------------------------------------------------------
write (ioutput,*) '************************************'
!%-----------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Phase:'
print *,'Name:',this%pp%pp%name
print *,'Service: write_'
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
subroutine compute_act_coeff_aqphideal &
(this, &
g, &
c, &
iserror, &
ionstr)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute activity coefficients
!
!   $Arguments:
!
 
type(t_aqueousphaseideal), intent(in)             :: this     ! Ideal aqueous phase type 

real*8, pointer, dimension(:)                     :: g        ! Activity coefficients vector

real*8, intent(in), dimension(this%pp%pp%numsp)   :: c        ! Molality vector 

logical, intent(out)                              :: iserror  ! iserror=true, then there was an error

real*8, intent(out), optional                     :: ionstr   ! Ionic strength 
 
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
haveionstr 
character(len=100)            :: &
msg
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Check optional arguments 
!%-----------------------------------------------------------
haveionstr = present (ionstr)
!%-----------------------------------------------------------
call check_pointer_ (g,this%pp%pp%numsp,.true.)
g=1.0d0
!%-----------------------------------------------------------
if (haveionstr) then
   call compute_ionstr_ (this%pp,ionstr,c,iserror)
   if (iserror) goto 10 
end if
!%------------------------------------------------------------
return
10 continue 
print *,'***************************'
print *,'Phase:'
print *,'Name:',this%pp%pp%name
print *,'Service: compute_act_coeff_'
print *,msg
print *,'***************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dact_coeff_aqphideal &
(this, & 
dg,  & 
c,  &  
dc, &
dionstr)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute derivatives of the activity coefficients
!
!   $Arguments:
!
    
type(t_aqueousphaseideal), intent(in)             :: this    ! Ideal aqueous phase type 
 
real*8, pointer, dimension(:,:)                   :: dg      ! Derivatives of the activity coefficients

real*8, intent(in), dimension(this%pp%pp%numsp)   :: c       ! Molality vector 

real*8, intent(in), dimension(:,:)                :: dc      ! Derivatives of the molalities

real*8, pointer, dimension(:), optional           :: dionstr ! Derivative of the ionic strength 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer               :: &
dionstrloc(:) => null ()
logical                       :: &
havedionstr, &
iserror
integer                       :: &
ndimder 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
havedionstr=present(dionstr)
!%------------------------------------------------------------
ndimder=size(dc,2)
!%-----------------------------------------------------------
call check_pointer_ (dg,this%pp%pp%numsp,ndimder,.true.)
!%-----------------------------------------------------------
!% Compute derivatives of the ionic strength
!%-----------------------------------------------------------
call compute_dionstr_ (this%pp,dionstrloc,dc,iserror)
!%-----------------------------------------------------------
if (havedionstr) then
 call check_pointer_ (dionstr,ndimder,.true.)
 dionstr=dionstrloc
end if
!%-----------------------------------------------------------
call check_pointer_ (dionstrloc,1,.false.)
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
end module 