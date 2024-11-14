module m_gasidealphase
!-------------------------------------------------------------------------
!
!   $Description: This class represents a gas phase with ideal 
! thermodynamic behavior
!
!   $Use: use m_parentphase
! use m_general_tools_cheproo
! use m_constants
! use flib_xpath
! use flib_sax
!
!   $Author: Sergio Andrés Bea Jofré
!
!   $License: CSIC-UPC
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_parentphase
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
private            ::
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!% Public services 
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
public             :: &
create_ &                   ! Create an object 
,set_ &                     ! Set general attributes in a gas phase with ideal thermodynamic behavior.
,set_pparent_ &             ! Set the parent phase 
,destroy_ &                 ! Destroy an object.
,compute_act_coeff_ &       ! Compute activity coefficients.
,compute_dact_coeff_ &      ! Compute derivatives of the activity coefficients.
,update_ &                  ! Update parameters that depends of the temperature in the ideal gas phase object. 
,write_ &                   ! Write in ascii all attributes encapsulated in the gas phase object. 
,assignment(=)              ! Copy an object in other object. 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!% Type definition 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
type, public::t_gasidealphase
 
type (t_parentphase), pointer   :: pp    ! Type parent phase 
 
end type t_gasidealphase
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
 
module procedure create_gasidealph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
 
module procedure destroy_gasidealph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
 
module procedure set_gasidealph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_pparent_
 
module procedure set_pparent_gasidealph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_act_coeff_
 
module procedure compute_act_coeff_gasidealph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dact_coeff_
 
module procedure compute_dact_coeff_gasidealph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_
 
module procedure update_temp_param_gasidealph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
 
module procedure write_gasidealph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment(=)
 
module procedure copy_gasidealph
 
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
subroutine create_gasidealph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create an object. 
!
!   $Arguments:
!
 
type(t_gasidealphase), intent(inout) :: this   ! Type ideal gas phase variable.
 
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
subroutine destroy_gasidealph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy gas ideal object. 
!
!   $Arguments:
!
 
type(t_gasidealphase), intent(inout) :: this  ! Type ideal gas phase variable.
 
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
subroutine set_gasidealph &
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
 
type(t_gasidealphase), intent(inout)    :: this     ! Type ideal gas phase variable.

type(T_ParentPhase), intent(in), target :: pp       ! Type parent phase variable. 

logical, intent(out)                    :: iserror  ! iserror=true, then there was an error. 
 
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
 
iserror=.false. 
!%------------------------------------------------------------
this%pp => pp 
!%------------------------------------------------------------
!% Update according reference temperature
!%------------------------------------------------------------
call update_ (this,this%pp%tempref,iserror)
if (iserror) goto 10
!%------------------------------------------------------------
return
10 continue
iserror=.true. 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_pparent_gasidealph &
   (this, &
    pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_gasidealphase), intent(inout)     :: this   ! Type ideal gas phase variable.

type(T_ParentPhase), pointer             :: pp     ! Type parent phase variable 
 
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
this%pp => pp
!%------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_gasidealph &
   (this, &
    temp, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update parameters that depends of the temperature in the 
! ideal gas phase object. 
!
!   $Arguments:
!
 
type(t_gasidealphase), intent(inout)          :: this     ! Type ideal gas phase variable.

real*8, intent(in)                            :: temp     ! Temperature (in celcius). 

logical, intent(out)                          :: iserror  ! iserror=true, then there was an error. 
 
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
subroutine write_gasidealph &
   (this, &
    ioutput, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii all attributes encapsulated in the gas phase object. 
!
!   $Arguments:
!
 
type(t_gasidealphase), intent(in)     :: this     ! Type ideal gas phase variable.

integer, intent(in)                   :: ioutput  ! Output unit 

logical, intent(out)                  :: iserror  ! iserror=true, then there was an error. 
 
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
 
iserror=.false.
!%------------------------------------------------------------
write (ioutput,*) '------------------------------------'
write (ioutput,*) '-          Gas Phase               -'
write (ioutput,*) '-     Activity model: Ideal        -'
!%-----------------------------------------------------------
!% Write the parent phase 
!%-----------------------------------------------------------
call write_ (this%pp,ioutput,iserror)
!%-----------------------------------------------------------
write (ioutput,*) '------------------------------------'
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_act_coeff_gasidealph &
   (this, &
    g, &
    c, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute activity coefficients.
!
!   $Arguments:
!
 
type(t_gasidealphase), intent(in)             :: this      ! Type ideal gas phase variable.

real*8, pointer, dimension(:)                 :: g         ! Activity coefficients vector.

real*8, intent(in), dimension(this%pp%numsp)  :: c         ! Concentration vector. 

logical, intent(out)                          :: iserror   ! iserror=true, then there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                      :: &
msg 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false. 
msg=''
!%------------------------------------------------------------
!% Allocate pointer 
!%------------------------------------------------------------
call check_pointer_ (g,this%pp%numsp,.true.)
!%------------------------------------------------------------
g = 1.0d0
!%------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_act_coeff_'
print *, msg
print *,'*******************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
subroutine compute_dact_coeff_gasidealph &
   (this, &
    dg, &
    c, &
    dc, &
	iserror, &
    g,&
    dparam)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute derivatives of the activity coefficients.
!
!   $Arguments:
!
 
 
type(t_gasidealphase), intent(in)                        :: this      ! Type ideal gas phase variable.

real*8,pointer, dimension(:,:)                           :: dg        ! Derivatives of the activity coefficients.

real*8, intent(in), dimension(this%pp%numsp)             :: c         ! Concentration vector. 

real*8, intent(in), dimension(:,:)                       :: dc        ! Derivatives of the concentrations. 

real*8, pointer, optional, dimension(:)                  :: dparam    ! Derivatives of some parameters. 

logical, intent(out)                                     :: iserror   ! iserror=true, then there was an error. 

real*8, intent(in), optional, dimension(this%pp%numsp)   :: g         ! Activity coefficients vector.
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical                  :: &
 havedparam, &
 haveg
integer                  :: &
 ndim 
character(len=100)       :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg='' 
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
havedparam=present(dparam)
haveg=present(g)
!%------------------------------------------------------------
ndim=size(dc,2)
!%------------------------------------------------------------
!% Allocate pointer 
!%------------------------------------------------------------
call check_pointer_ (dg,this%pp%numsp,ndim,.true.)
!%------------------------------------------------------------
if (havedparam) then
  call check_pointer_ (dparam,ndim,.true.)
end if
!%------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_dct_coeff_'
print *, msg
print *,'*******************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_gasidealph &
  (tarjetobj, &
   sourceobj)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy an object in other object. 
!
!   $Arguments:
!
 
type(t_gasidealphase), intent(in) :: sourceobj   ! Type ideal gas phase variable.

type(t_gasidealphase), intent(out):: tarjetobj   ! Type ideal gas phase variable.
 
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
tarjetobj%pp => null ()
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
end module m_gasidealphase
