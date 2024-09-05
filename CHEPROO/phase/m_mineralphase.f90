module m_mineralphase
!-------------------------------------------------------------------------
!
!   $Description: This class represents a mineral phase and its 
! thermodynamic behavior 
!
!   $Use: use m_parentphase
!use m_mineralidealphase
!use m_binarysolutionphase
!use m_general_tools_cheproo
!use m_constants
!use flib_xpath
!use flib_sax
!
!   $Author: Sergio Andrés Bea Jofré
!
!   $License:
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_parentphase
use m_mineralidealphase
use m_binarysolutionphase
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
create_ &              ! Create phase object.
,set_ &                ! Set attributes in the phase object 
,read_xml_ &           ! Read the object from xml file 
,set_pparent_ &        ! Set pointer to parent phase 
,get_pparent_ &        ! Return pointer to parent phase 
,destroy_ &            ! Destroy phase object.
,compute_act_coeff_ &  ! Compute activity coefficients
,compute_dact_coeff_ & ! Compute the derivatives of the activity coefficients.  
,update_ &             ! Update paramaters that depends of the temperature in the phase
,write_ &              ! Write in ascii the attributes encapsulated in the object 
,assignment(=)         ! Copy a phase object in other phase object. 
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
private             :: &
read_xml_loc &
,begin_element_handler 
!%----------------------------------------------------
!%---------------------------------Parameters constant
!%----------------------------------------------------
integer, parameter   :: &
minidealph = 1, &
binsolph = 2
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!% Type definition 
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
type, public::t_mineralphase

private                               :: 

type (t_parentphase), pointer         :: pp            ! Parent phase 

type (t_mineralidealphase), pointer   :: pminidealph   ! Mineral phase with ideal thermodynamic behavior 

type (t_binarysolutionphase), pointer :: pbinsolph     ! Binary solution 

integer                               :: itype         ! Specialization indice 
 
end type t_mineralphase
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
 
module procedure create_minph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
 
module procedure destroy_minph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
 
module procedure set_minph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_pparent_
 
module procedure set_pparent_minlph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_pparent_
 
module procedure get_pparent_minlph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_act_coeff_
 
module procedure compute_act_coeff_minph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dact_coeff_
 
module procedure compute_dact_coeff_minph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_
 
module procedure update_temp_param_minph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
 
module procedure write_minph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment(=)
 
module procedure copy_minph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface read_xml_
 
module procedure read_xml_minph
 
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
subroutine create_minph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create phase object.
!
!   $Arguments:
!
 
type(t_mineralphase), intent(inout) :: this  ! Type mineral phase variable
 
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
this%pminidealph => null ()
this%pbinsolph => null ()
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
subroutine destroy_minph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy phase object.
!
!   $Arguments:
!
 
type(t_mineralphase), intent(inout) :: this  ! Type mineral phase variable
 
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
select case (this%itype)
case (minidealph)  
 call destroy_ (this%pminidealph)
 deallocate (this%pminidealph)
 this%pminidealph => null ()
case (binsolph)   ! Pitzer model
 call destroy_ (this%pbinsolph)
 deallocate (this%pbinsolph)
 this%pbinsolph => null ()
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
subroutine read_xml_minph &
   (this, &
    namefile, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read phase from xml file
!
!   $Arguments:
!
 
type(t_mineralphase), intent(inout)     :: this     ! Type mineral phase variable

character(len=*), intent(in)            :: namefile

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
call set_minph &
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
subroutine set_minph &
   (this, &
    pp, &
	nameactmodel, &
	iserror, &
	namedatabase)
implicit none
!-------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_mineralphase), intent(inout)       :: this             ! Type mineral phase variable

type(t_parentphase), intent(in), target   :: pp               ! Type parent phase 

character(len=*), intent(in)              :: nameactmodel     ! Name of the thermodynamic model

logical, intent(out)                      :: iserror          ! iserror=true, then there was an error

character(len=*), intent(in), optional    :: namedatabase     ! Name of some thermodynamic data base 
 
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
 this%itype=minidealph
 allocate (this%pminidealph) 
 call create_ (this%pminidealph)
 call set_ (this%pminidealph,this%pp,iserror)
case('nonideal','NONIDEAL')
 this%itype=binsolph
 allocate (this%pbinsolph) 
 call create_ (this%pbinsolph)
 call set_ (this%pbinsolph,this%pp,iserror)
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
subroutine set_pparent_minlph &
   (this, &
    pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_mineralphase), intent(inout)     :: this    ! Type mineral phase variable

type(t_parentphase), intent(in), target :: pp      ! Type parent phase 
 
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
case (minidealph)
 call set_pparent_ (this%pminidealph,this%pp)
case (binsolph)
 call set_pparent_ (this%pbinsolph,this%pp)
end select 
!%------------------------------------------------------------ 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_pparent_minlph &
   (this, &
    pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the parent phase pointer 
!
!   $Arguments:
!
 
type(t_mineralphase), intent(in)     :: this    ! Type mineral phase variable

type(t_parentphase), pointer         :: pp      ! Type parent phase 
 
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
subroutine update_temp_param_minph &
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
 
type(t_mineralphase), intent(inout)      :: this       ! Type mineral phase variable

real*8, intent(in)                       :: temp       ! Temperature in celcius

logical, intent(out)                     :: iserror    ! iserror=true, then there was an error
 
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
case (minidealph)
 call update_ (this%pminidealph,temp,iserror)
case (binsolph)
 call update_ (this%pbinsolph,temp,iserror)
end select 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_minph &
   (this, &
    ioutput, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii the attributes encapsulated in the object 
!
!   $Arguments:
!
 
type(t_mineralphase), intent(in) :: this    ! Type mineral phase variable

integer, intent(in)              :: ioutput ! Output unit 

logical, intent(out)             :: iserror ! iserror=true, then there was an error
 
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
select case (this%itype)
case (minidealph)
 call write_ (this%pminidealph,ioutput,iserror)
case (binsolph)
 call write_ (this%pbinsolph,ioutput,iserror)
end select
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_act_coeff_minph &
   (this, &
    g, &
    c, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute activity coefficients
!
!   $Arguments:
!
 
type(t_mineralphase), intent(in)              :: this      ! Type mineral phase variable 

real*8, pointer, dimension(:)                 :: g         ! Activity coefficient vector 

real*8, intent(in), dimension(this%pp%numsp)  :: c         ! Concentration vector 

logical, intent(out)                          :: iserror   ! iserror=true, then there was an error
 
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
case (minidealph)
 call compute_act_coeff_ (this%pminidealph,g,c,iserror)
case (binsolph)
 msg='Warning, not implemented service'
 goto 10
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
subroutine compute_dact_coeff_minph &
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
!   $Description: ! Compute the derivatives of the activity coefficients.  
!
!   $Arguments:
!
 
 
type(t_mineralphase), intent(in)                        :: this     ! Type mineral phase variable 

real*8, intent(in), dimension(this%pp%numsp)            :: c        ! Concentration vector

real*8, intent(in), dimension(:,:)                      :: dc       ! Derivatives of the concentrations  

real*8,pointer, dimension(:,:)                          :: dg       ! Derivatives of the activity  coefficients 

real*8, intent(in), optional, dimension(this%pp%numsp)  :: g        ! Activity coefficients
 
real*8, pointer, optional, dimension(:)                 :: dparam   ! Derivatives of some parameters 

logical, intent(out)                                    :: iserror  ! iserror=true, then there was an error
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
case (minidealph)
 call compute_dact_coeff_ (this%pminidealph,dg,c,dc,iserror,g,dparam)
case (binsolph)
 msg='Warning, not implemented service'
 goto 10
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
subroutine copy_minph &
  (copied, &
   this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy a object in other species object.  
!
!   $Arguments:
!
 
type(t_mineralphase), intent(in) :: this     ! Type mineral phase variable  

type(t_mineralphase), intent(out):: copied   ! Type mineral phase variable 
 
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
!% Copy the specialization 
!%------------------------------------------------------------
select case (this%itype)
 
case (minidealph)
 
 allocate(copied%pminidealph)
 call create_ (copied%pminidealph)
 copied%pminidealph=this%pminidealph
 call set_pparent_ (copied%pminidealph,copied%pp)
 
case (binsolph)
 
 allocate(copied%pbinsolph)
 call create_ (copied%pbinsolph)
 copied%pbinsolph=this%pbinsolph
 call set_pparent_ (copied%pbinsolph,copied%pp)

end select 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%
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
!   $Description: Read phase from xml file (private service)
!
!   $Arguments:
!
 
type(dictionary_t), intent(in)  :: attributes

type(t_mineralphase), optional  :: minph

logical, intent(out), optional  :: iserror 

character(len=*), optional      :: nameactmodel 

character(len=*)                :: name                  
 
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
end module m_mineralphase
