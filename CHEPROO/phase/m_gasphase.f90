module m_gasphase
!-------------------------------------------------------------------------
!
!   $Description:  This class represents a gas phase and its 
! thermodynamic behavior 
!
!   $Use: use m_parentphase
!use m_gasidealphase
!use m_general_tools_cheproo
!use m_constants
!use flib_xpath
!use flib_sax
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
use m_gasidealphase
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
public             :: &
create_ &                  ! Create a gas phase object 
,destroy_ &                ! Destroy a gas phase object 
,set_ &                    ! Set attributes in the gas phase. 
,read_xml_ &               ! Read a gas phase object from xml file 
,set_pparent_ &            ! Set the parent phase 
,get_pparent_ &            ! Return pointer to parent phase 
,compute_act_coeff_ &
,compute_dact_coeff_ &     ! Compute derivatives of the activity coefficients 
,update_ &
,write_ &                  ! Write in ascii the attributes encapsulated in the gas phase 
,assignment(=)
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
private             :: &
read_xml_loc &
,begin_element_handler 
!%----------------------------------------------------
!% Parameters constant
!%----------------------------------------------------
integer, parameter   :: &
gasidealph = 1
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
type, public::t_gasphase

private                               :: 

type (t_parentphase), pointer         :: pp           ! Parent phase. 

type (t_gasidealphase), pointer       :: pgasidealph  ! Gas phase with ideal behavior.

integer                               :: itype        ! Specialization indice.
 
end type t_gasphase
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
 
module procedure create_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
 
module procedure destroy_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
 
module procedure set_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_pparent_
 
module procedure set_pparent_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_pparent_
 
module procedure get_pparent_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_act_coeff_
 
module procedure compute_act_coeff_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dact_coeff_
 
module procedure compute_dact_coeff_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_
 
module procedure update_temp_param_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
 
module procedure write_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment(=)
 
module procedure copy_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface read_xml_
 
module procedure read_xml_gasph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface read_xml_loc
 
module procedure read_xml_loc
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface begin_element_handler
 
module procedure begin_element_handler
 
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
subroutine create_gasph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create a gas phase object 
!
!   $Arguments:
!
 
type(t_gasphase), intent(inout) :: this  ! Type gas phase variable 
 
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
this%pgasidealph => null ()
this%itype=0 
!%------------------------------------------------------------
allocate (this%pp)
call create_ (this%pp)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_gasph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy a gas phase object
!
!   $Arguments:
!
 
type(t_gasphase), intent(inout) :: this   ! Type gas phase variable 
 
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
!% Destroy and deallocate the specialization 
!%------------------------------------------------------------
select case (this%itype)
case (gasidealph)  
 call destroy_ (this%pgasidealph)
 deallocate (this%pgasidealph)
 this%pgasidealph => null ()
end select
!%-----------------------------------------------------------
this%itype=0
!%-----------------------------------------------------------
return 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_gasph &
   (this, &
    namefile, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read a gas phase object from xml file 
!
!   $Arguments:
!
 
type(t_gasphase), intent(inout)         :: this      ! Type gas phase variable 

character(len=*), intent(in)            :: namefile  ! Name of the xml file 

logical, intent(out)                    :: iserror
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer:: &
 iostat
type(xml_t):: &
 fxml
type(dictionary_t)          :: &
 attributes
character(len=100)  :: &
 nameactmodel, &
 msg, &
 name  
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%----------------------------------------------------------------
call open_xmlfile(namefile, fxml, iostat)
 
if (iostat.ne.0) then
 msg='Error opening file'
 call add_ (msg,namefile)
 goto 10
end if
!%----------------------------------------------------------------
call xml_parse(fxml, &
         begin_element_handler = begin_element_handler, &
         verbose = .false.)
!%--------------------------------------------------------
! End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%--------------------------------------------Read parent phase
call read_xml_ &
   (this%pp, &
    namefile, &
    iserror)
!%-------------------------------------------------------------
call read_xml_loc &
  (name, &
   attributes, &
   minph=this, &
   iserror=iserror, &
   nameactmodel=nameactmodel)
!%------------------------------------------------------------
if (iserror) goto 10
!%-------------------------------------------------------------
call set_gasph &
   (this, &
    this%pp, &
	nameactmodel, &
	iserror) 
!%------------------------------------------------------------
if (iserror) goto 10
!%-----------------------------------------------------------
! End and close xml file
!%-----------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Phase:'
print *,'Name:',this%pp%name
print *,'Service: read_xml_'
print *, msg
print *,'***********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_gasph &
   (this, &
    pp, &
	nameactmodel, &
	iserror)
implicit none
!-------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_gasphase), intent(inout)           :: this         ! Type gas phase variable 

type(T_ParentPhase), intent(in), target   :: pp

character(len=*), intent(in)              :: nameactmodel ! Name of the thermodynamic model 

logical, intent(out)                      :: iserror 
 
!-------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!--------------------------------------------------------------
character(len=100):: &
 msg
!--------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
iserror=.false.
msg='' 
!%------------------------------------------------------------
this%pp => pp
!%------------------------------------------------------------
select case(nameactmodel)
case('ideal','IDEAL')
 this%itype=gasidealph
 allocate (this%pgasidealph) 
 call create_ (this%pgasidealph)
 call set_ (this%pgasidealph,this%pp,iserror)
case default 
 msg='Error, not recognized thermodynamic model'
 call add_ (msg,nameactmodel)
 goto 10
end select 
!%------------------------------------------------------------ 
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: set_'
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
subroutine set_pparent_gasph &
   (this, &
    pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the parent phase 
!
!   $Arguments:
!
 
type(t_gasphase), intent(inout)         :: this ! Type gas phase variable 

type(T_ParentPhase), intent(in), target :: pp
 
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
select case (this%itype)
case (gasidealph)
 call set_pparent_ (this%pgasidealph,this%pp)
end select 
!%------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_pparent_gasph &
   (this, &
    pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return pointer to parent phase 
!
!   $Arguments:
!
 
type(t_gasphase), intent(in)     :: this   ! Type gas phase variable 

type(t_parentphase), pointer     :: pp     ! Type parent phase
 
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
pp => this%pp
!%------------------------------------------------------------ 
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_gasph &
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
 
type(t_gasphase), intent(inout)          :: this     ! Type gas phase variable 

real*8, intent(in)                       :: temp

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
 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
iserror=.false.
!%------------------------------------------------------------
if (this%pp%tempref==temp) return 
!%------------------------------------------------------------
this%pp%tempref=temp 
!%------------------------------------------------------------
select case (this%itype)
case (gasidealph)
 call update_ (this%pgasidealph,temp,iserror)
end select 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_gasph &
   (this, &
    ioutput, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii the attributes encapsulated in the gas phase 
!
!   $Arguments:
!
 
type(t_gasphase), intent(in)     :: this     ! Type gas phase variable 

integer, intent(in)              :: ioutput

logical, intent(out)             :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)         :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
iserror=.false.
msg='' 
!%------------------------------------------------------------
!%  Write the specialization 
!%------------------------------------------------------------
select case (this%itype)
case (gasidealph)
 call write_ (this%pgasidealph,ioutput,iserror)
case default
 msg='Error, service not implemented'
 goto 10 
end select
!%-----------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: write_'
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
subroutine compute_act_coeff_gasph &
   (this, &
    g, &
    c, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_gasphase), intent(in)                   :: this

real*8, pointer, dimension(:)                  :: g

real*8, intent(in), dimension(this%pp%numsp)   :: c

logical, intent(out)                           :: iserror 
 
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
select case (this%itype)
case (gasidealph)
 
 call compute_act_coeff_ (this%pgasidealph,g,c,iserror)

end select 
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
subroutine compute_dact_coeff_gasph &
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
!   $Description: Compute derivatives of the activity coefficients 
!
!   $Arguments:
!
 
 
type(t_gasphase), intent(in)                               :: this

real*8, intent(in), dimension(this%pp%numsp)               :: c

real*8, intent(in), dimension(:,:)                         :: dc

real*8,pointer, dimension(:,:)                             :: dg

real*8, intent(in), optional, dimension(this%pp%numsp)     :: g
 
real*8, pointer, optional, dimension(:)                    :: dparam

logical, intent(out)                                       :: iserror 
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
select case (this%itype)
case (gasidealph)
 call compute_dact_coeff_ &
   (this%pgasidealph, &
    dg, &
    c, &
    dc, &
	iserror, &
    g,&
    dparam)

end select 
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
subroutine copy_gasph &
  (copied, &
   this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy an gas phase object in other gas phase object.  
!
!   $Arguments:
!
 
type(t_gasphase), intent(in) :: this

type(t_gasphase), intent(out):: copied 
 
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
copied%itype=this%itype
!%------------------------------------------------------------
!% Copy the parent phase 
!%------------------------------------------------------------
copied%pp = this%pp
!%------------------------------------------------------------
select case (this%itype)
 
case (gasidealph)
 
 allocate(copied%pgasidealph)
 call create_ (copied%pgasidealph)
 copied%pgasidealph=this%pgasidealph
 call set_pparent_ (copied%pgasidealph,copied%pp)
 
end select 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
subroutine begin_element_handler(name,attributes)
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
integer        :: status
 
call read_xml_loc (name,attributes)
 
return
end subroutine begin_element_handler
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc &
  (name, &
   attributes, &
   minph, &
   iserror, &
   nameactmodel)
implicit none 
!-------------------------------------------------------------------------
!
!   $Description: Read phase from xml file
!
!   $Arguments:
!
 
type(dictionary_t), intent(in)  :: attributes

type(t_gasphase), optional      :: minph

logical, intent(out), optional  :: iserror 

character(len=*), optional      :: nameactmodel    ! Name of the thermodynamic model 

character(len=*)                :: name            ! Name of the gas phase      
 
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
 status
logical         :: &
 haveminph, &
 haveiserror
character(len=100)     :: &
 id
character(len=100),save:: &
 namemodel
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
haveminph=present(minph)
haveiserror=present(iserror)
!%-----------------------------------------------------------
if (haveiserror) iserror=.false.
!%-----------------------------------------------------------
select case (haveminph)
case (.true.)
  
 nameactmodel=namemodel
 
case default
 
  select case (name)
 
  case ('phase')
    id='' 
	call get_value (attributes,"model", id, status)
    namemodel=id
  end select  


end select
 
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
end module m_gasphase
