module m_parentphase
!-------------------------------------------------------------------------
!
!   $Description: Parent phase
!
!   $Use: use m_species
! use flib_xpath
! use flib_sax
! use m_general_tools_cheproo
! use m_constants
!
!   $Author: Sergio Andr�s Bea Jofr�. 
!
!   $License: 
!
!-------------------------------------------------------------------------
use m_species
use flib_xpath
use flib_sax
use m_general_tools_cheproo
use m_constants_cheproo
!%------------------------------------------------------
!%------------------------------------------------------
private::
!%------------------------------------------------------
!% Public services 
!%------------------------------------------------------
public:: &
create_ &                    ! Create the parent phase object. 
,destroy_ &                  ! Destroy the parent phase object. 
,read_xml_ &                 ! Read parent phase object from xml file. 
,set_ &                      ! Set attributes in the parent phase object 
,get_pspecies_ &             ! Return the pointer to species object. 
,get_namesp_ &               ! Return the name of the species. 
,get_name_ &                 ! Return the name of the phase 
,get_numsp_ &                ! Return the number of species 
,get_prop_ &                 ! Return some property according its name. 
,get_if_sp_is_present_ &     ! Return if the species is present in the phase
,compute_charge_balance_ &   ! Compute the charge balance. 
,compute_dcharge_balance_ &  ! Compute the derivatives of the charge balance. 
,compute_sum_c_ &            ! Copute SUM(c)
,assignment(=) &             ! Copy an object in other object. 
,write_                      ! Write in ascii the attributtes of parent phase object
!%------------------------------------------------------
!%------------------------------------------------------
private:: &
read_xml_loc_ &
,begin_element_handler
!%------------------------------------------------------
!%------------------------------------------------------
!% Type definition 
!%------------------------------------------------------
!%------------------------------------------------------
type, public::t_parentphase
 
character(len=100)                        :: name           ! Name of the phase

character(len=100)                        :: typephase      ! Type of phase 
 
type(t_pspecies), pointer, dimension(:)   :: pspecies       ! List of species pointers
 
integer                                   :: numsp          ! Number of species

integer                                   :: numprop        ! Number of properties
 
character(len=100), pointer, dimension(:) :: nameprop       ! Name of the properties [numprop]

character(len=100), pointer, dimension(:) :: modelprop      ! Name of the model property [numprop]
 
real*8, pointer, dimension(:)             :: prop           ! Value of the properties [numprop]

real*8                                    :: tempref        ! Temperature of reference 

real*8                                    :: pressref       ! Pressure of reference 
 
end type t_parentphase
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!% Public interfaces 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
 
module procedure create_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
 
module procedure destroy_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface read_xml_
 
module procedure read_xml_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
 
module procedure set_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_namesp_
 
module procedure get_namesp_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_prop_
 
module procedure get_prop_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_pspecies_
 
module procedure get_pspecies_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_numsp_
 
module procedure get_numsp_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_name_
 
module procedure get_name_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_if_sp_is_present_
 
module procedure get_if_sp_is_present_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_charge_balance_
 
 
module procedure compute_charge_balance_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_sum_c_
 
module procedure compute_sum_c_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dcharge_balance_
 
module procedure compute_dcharge_balance_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
 
module procedure write_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment(=)
 
module procedure copy_pph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!% Private interfaces 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface read_xml_loc_
 
module procedure read_xml_loc_pph
 
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
subroutine create_pph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create parent phase object 
!
!   $Arguments:
!
 
type(t_parentphase), intent(inout) :: this    ! Type parent phase variable. 
 
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
this%numsp = 0
this%numprop = 0
this%name = ''
this%typephase = ' '
this%tempref=25.0d0
this%pressref=0.0d0
!%-------------------------------------------------------------
this%pspecies => null ()
this%prop => null ()
this%nameprop => null ()
this%modelprop => null ()
!%-------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_pph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy parent phase object 
!
!   $Arguments:
!
 
type(t_parentphase), intent(inout) :: this   ! Type parent phase variable. 
 
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
 i 
!-------------------------------------------------------------------------
!
!   $code
!
this%numprop = 0
this%name = ''
this%typephase = ' '
this%tempref=0.0d0
this%pressref=0.0d0 
!%-------------------------------------------------------------
if (this%numsp>0) then
 do i=1,this%numsp
  call destroy_ (this%pspecies(i)%ptr)
  deallocate (this%pspecies(i)%ptr)
  this%pspecies(i)%ptr => null ()
 end do
 deallocate (this%pspecies)
 this%pspecies => null ()
end if
!%-------------------------------------------------------------
this%numsp=0
!%-------------------------------------------------------------
call check_pointer_ (this%prop,1,.false.)
call check_pointer_ (this%nameprop,1,.false.)
call check_pointer_ (this%modelprop,1,.false.)
!%-------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_pph &
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
 
type(t_parentphase), intent(inout)      :: this      ! Type parent phase variable. 

character(len=*), intent(in)            :: namefile  ! Name of the xml file 

logical, intent(out)                    :: iserror   ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: The namefile must contain the path
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
 name
character(len=100)  :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%----------------------------------------------------------------
!% Open xml file
!%----------------------------------------------------------------
call open_xmlfile(namefile,fxml, iostat)
if (iostat/=0) then
 msg='Error opening file'
 call add_ (msg,namefile)
 goto 10
end if
!%----------------------------------------------------------------
call xml_parse(fxml, &
         begin_element_handler = begin_element_handler, &
         verbose = .false.)
!%--------------------------------------------------------
!% End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%----------------------------------------------------------------
! Read xml file 
!%----------------------------------------------------------------
call read_xml_loc_ (name,attributes,pph=this,iserror=iserror)
if (iserror) goto 10 
!%-----------------------------------------------------------
!% Close xml file
!%-----------------------------------------------------------
call close_xmlfile(fxml)
!%------------------------------------------------------------
return
 
10 continue 
print *,'*******************'
print *,'Phase:'
print *,'Name:',this%name
print *,'Service: read_xml_'
print *, msg
print *,'*******************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_pph &
   (this, &
    species, &
    numsp, &
    numprop, &
    namephase, &
    typephase, &
    nameprop, &
	modelprop, &
    prop, &
    msg, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set attributes in the parent phase object 
!
!   $Arguments:
!
 
 
type(t_parentphase), intent(inout)                    :: this       ! Type parent phase variable. 

integer, intent(in)                                   :: numsp      ! Number of species 

integer, intent(in)                                   :: numprop    ! Number of properties 

type(t_species), intent(in), dimension(numsp), target :: species    ! List of the species objects. 

character(len=*), intent(in)                          :: namephase  ! Name of the phase

character(len=*), intent(in)                          :: typephase  ! Type of phase 

character(len=*), intent(in), dimension(numprop)      :: nameprop   ! Name of the properties [numprop]

character(len=*), intent(in), dimension(numprop)      :: modelprop  ! Name of the properties [numprop]

real*8, intent(in), dimension(numprop)                :: prop       ! Values of the properties [numprop]

character(len=*), intent(out)                         :: msg        ! Error message 

logical, intent(out)                                  :: iserror    ! iserror=true, then there was an error

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                             :: &
 i, &
 j, &
 ndim, &
 isp
character(len=100)                  :: &
 name
logical                             :: &
 isrepeated 
!-------------------------------------------------------------------------
!
!   $code
!
this%numsp=numsp
!%-----------------------------------------------------------
!% Check if some species is repeated in the list 
!%-----------------------------------------------------------
if (this%numsp>0) then
 allocate (this%pspecies(this%numsp))
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
end if
!%-----------------------------------------------------------
!% Allocate, create and copy the species 
!%-----------------------------------------------------------
do i=1,this%numsp
 allocate (this%pspecies(i)%ptr)
 call create_ (this%pspecies(i)%ptr)
 this%pspecies(i)%ptr = species (i)
end do
!%-----------------------------------------------------------
this%name = namephase
this%typephase = typephase
!%-----------------------------------------------------------
this%numprop=numprop
!%-----------------------------------------------------------
!% Set the properties 
!%-----------------------------------------------------------
if (this%numprop/=0) then
 call check_pointer_ (this%nameprop,ndim,.true.)
 call check_pointer_ (this%modelprop,ndim,.true.)
 this%nameprop = nameprop
 this%prop = prop
 this%modelprop = modelprop
end if
!%-----------------------------------------------------------
 
return
10 continue
print *,'********************************'
print *,'Phase:'
print *,'Name:', this%name
print *,'service: set_'
print *, msg
print *,'********************************'
iserror=.true.
return
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_namesp_pph &
   (this, &
    namesp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_parentphase), intent(in)           :: this     ! Type parent phase variable. 
 
character(len=*), pointer, dimension(:)   :: namesp   ! Name of the species 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                           :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%------------------------------------------------------------
call check_pointer_ (namesp,this%numsp,.true.)
!%------------------------------------------------------------
do i=1,this%numsp
 namesp(i) = this%pspecies(i)%ptr%name
end do
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_numsp_pph &
   (this, &
    numsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the number of species
!
!   $Arguments:
!
 
type(t_parentphase), intent(in)   :: this      ! Type parent phase variable. 

integer, intent(out)              :: numsp     ! Number of species 
 
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
numsp=this%numsp 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_pspecies_pph &
   (this, &
    pspecies, &
    ithsp, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return pointer to ith species object 
!
!   $Arguments:
!
 
type(t_parentphase), intent(in) :: this         ! Type parent phase variable. 

type(t_species), pointer        :: pspecies     ! Species pointer 

integer, intent(in)             :: ithsp        ! Local indice of the species 

logical, intent(out)            :: iserror      ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 character(len=100)              :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg='' 
!%------------------------------------------------------------
!% Check the index of the species 
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
print *,'********************************'
print *,'Phase:'
print *,'Name:', this%name
print *,'service: get_pspecies_'
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
subroutine get_prop_pph &
   (this, &
    value, &
    nameprop, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return some property according its name. 
!
!   $Arguments:
!
type(t_parentphase), intent(in)          :: this      ! Type parent phase variable. 

character(len=*), intent(in)             :: nameprop  ! Name of the property 
 
real*8, intent(out)                      :: value     ! Value of the property 

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

logical             :: &
 isbe
integer             :: &
 i
character(len=100)  :: &
 msg
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
value=0.0d0
!%------------------------------------------------------------
isbe=.false.
do i=1,this%numprop
 if (this%nameprop(i)==nameprop) then
  isbe=.true.
  value=this%prop(i)
  exit 
 end if
end do
!%------------------------------------------------------------
if (.not.isbe) goto 10
!%------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%name
print *,'Service: get_prop_'
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
subroutine get_if_sp_is_present_pph &
   (this, &
    namesp, &
    isbe)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
type(t_parentphase), intent(in)    :: this     ! Type parent phase variable. 

character(len=*), intent(in)       :: namesp   ! Name of the species 

logical, intent(out)               :: isbe     ! isbe=true, the species is present in the phase 
 
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
 i 
!-------------------------------------------------------------------------
!
!   $code
!

 

!%------------------------------------------------------------
isbe = .false.
!%------------------------------------------------------------
do i=1,this%numsp
 if (this%pspecies(i)%ptr%name==namesp) then
   isbe = .true.
   exit
 end if
end do
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_name_pph &
   (this, &
    name)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of the phase 
!
!   $Arguments:
!
 
type(t_parentphase), intent(in)    :: this    ! Type parent phase variable. 

character(len=*), intent(out)      :: name    ! Name of the phase.
 
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
name=this%name
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_charge_balance_pph &
   (this, &
    chgbal, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the charge balance.
!
!   $Arguments:
!
 
type(t_parentphase), intent(in)             :: this    ! Type parent phase variable. 

real*8, intent(out)                         :: chgbal  ! Charge balance 

real*8, intent(in), dimension(this%numsp)   :: c       ! Concentration vector 

logical, intent(out)                        :: iserror ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8        :: &
 z
integer       :: &
 i
character(len=100) :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

iserror=.false.
msg=''
!%------------------------------------------------------------
chgbal=0.0d0
z=0.0d0
!%------------------------------------------------------------
do i=1,this%numsp
 
   call get_prop_ (this%pspecies(i)%ptr,z,'charge',msg,iserror)
 
   if (iserror) goto 10
 
   chgbal = chgbal + z * c(i)
 
end do
!%------------------------------------------------------------
return
 
10 continue
print *,'********************************'
print *,'Phase:'
print *,'Name:', this%name
print *,'service: compute_charge_balance_'
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
subroutine compute_sum_c_pph &
   (this, &
    sumc, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copute SUM(c)
!
!   $Arguments:
!
 
type(t_parentphase), intent(in)             :: this     ! Type parent phase variable. 

real*8, intent(out)                         :: sumc     ! SUM(c)

real*8, intent(in), dimension(this%numsp)   :: c        ! Concentration vector 

logical, intent(out)                        :: iserror  ! iserror=true, then there was an error
 
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

iserror=.false.
msg=''
!%------------------------------------------------------------
sumc=sum(c)
!%------------------------------------------------------------
return
 
10 continue
print *,'********************************'
print *,'Phase:'
print *,'Name:', this%name
print *,'service: compute_charge_balance_'
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
subroutine compute_dcharge_balance_pph &
   (this, &
    dchgbal, &
    dc, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the derivatives of the charge balance (chgbal).
! 
! dchgbal=SUM(zi*dci)
!
!   $Arguments:
!
 
type(t_parentphase), intent(in)          :: this     ! Type parent phase variable. 

real*8, pointer, dimension(:)            :: dchgbal  ! Derivatives of the charge balance

real*8, intent(in), dimension(:,:)       :: dc       ! Derivatives of the concentrations 

logical, intent(out)                     :: iserror  ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!dchgbal
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!
 

real*8        :: &
 z
integer       :: &
 i, &
 ndimder
character(len=100) :: &
 msg
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
ndimder=size(dc,2)
call check_pointer_ (dchgbal,ndimder,.true.)
!%------------------------------------------------------------
z=0.0d0
do i=1,this%numsp
 
  call get_prop_ (this%pspecies(i)%ptr,z,'charge',msg,iserror)
 
  if (iserror) goto 10
 
  dchgbal = dchgbal + z * dc(i,:)
 
end do
!%------------------------------------------------------------
return
10 continue 
print *,'********************************'
print *,'Phase:'
print *,'Name:', this%name
print *,'service: compute_dcharge_balance_'
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
subroutine write_pph &
   (this, &
    ioutput, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii the attributtes of parent phase object
!
!   $Arguments:
!
 
 
type(t_parentphase), intent(in)  :: this    ! Type parent phase variable. 

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
integer             :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!

 

!%------------------------------------------------------------
iserror=.false.
!%------------------------------------------------------------
write(ioutput,*) '-------------------------------------------------------------------'
write(ioutput,*) 'Phase information'
write(ioutput,*) 'Name:',this%name
write(ioutput,1) 'Numsp:',this%numsp
write(ioutput,2) 'Reference temperature:',this%tempref,'C'
write(ioutput,2) 'Reference pressure:',this%pressref,'MPa'
write(ioutput,*) '-------------------------------------------------------------------'
!%------------------------------------------------------------
!% Write tthe phase properties
!%------------------------------------------------------------
if (this%numprop>0) then
  do i=1,this%numprop
   write(ioutput,*) 'Phase property:',this%nameprop(i)
   write(ioutput,*) 'Model property:',this%modelprop(i)
   write(ioutput,*) 'Value property:',this%prop(i)
   write(ioutput,*) '-------------------------------------------------------------------'
  end do
end if 
!%------------------------------------------------------------
!% Write species objects 
!%------------------------------------------------------------
do i=1,this%numsp
 call write_ (this%pspecies(i)%ptr,ioutput,iserror)
end do
write(ioutput,*) '-------------------------------------------------------------------'
!%------------------------------------------------------------
return
1 format(a6,i5)
2 format(a22,f6.2,a3)
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_pph &
  (copied, &
   this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy a object in other species object.  
!
!   $Arguments:
!
 
type(t_parentphase), intent(in)  :: this    ! Type parent phase variable. 

type(t_parentphase), intent(out) :: copied  ! Type parent phase variable. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: The target object must be previously created
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
 

integer        :: &
 i
!%-----------------------------------------------------------
call destroy_ (copied)
!%-----------------------------------------------------------
copied%name = this%name
copied%typephase = this%typephase
copied%numsp = this%numsp
copied%numprop = this%numprop
copied%tempref=this%tempref
copied%pressref=this%pressref 
!%-----------------------------------------------------------
!% Copy the species 
!%-----------------------------------------------------------
allocate (copied%pspecies(this%numsp))
do i=1,this%numsp
 allocate(copied%pspecies(i)%ptr)
 call create_ (copied%pspecies(i)%ptr)
 copied%pspecies(i)%ptr = this%pspecies(i)%ptr
end do
!%-----------------------------------------------------------
!% Copy the properties 
!%-----------------------------------------------------------
if (this%numprop>0) then
 call check_pointer_ (copied%prop,copied%numprop,.true.)
 call check_pointer_ (copied%nameprop,copied%numprop,.true.)
 call check_pointer_ (copied%modelprop,copied%numprop,.true.)
 copied%prop=this%prop
 copied%nameprop=this%nameprop
 copied%modelprop=this%modelprop
end if
!%-----------------------------------------------------------
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
 
call read_xml_loc_ (name,attributes)
 
return
end subroutine begin_element_handler 
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc_pph &
  (name, &
   attributes, &
   pph,&
   iserror)
implicit none 
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
type(t_parentphase), optional  :: pph 
logical, intent(out), optional :: iserror  
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
type(t_species), pointer       :: &
 species(:)
integer:: &
 numsp
character(len=100), pointer:: &
 namepropph(:)
real*8, pointer :: &
 propph(:)
 
integer         :: &
 n, &
 status, &
 i  
logical         :: &
 havepph
real*8          :: &
 reallocal(1)
character(len=100)     :: &
 id, &
 msg 
logical, save          :: &
 isbeph, &
 isbesp 
character(len=100),save:: &
 namephase, &
 typephase, &
 nameactmodel, &
 nameconvention, &
 convph
integer, save          :: &
 isp, &
 iprop, &
 ipropph
integer, pointer, save :: &
 numprop(:)
character(len=100),pointer, save:: &
 namesp(:), &
 nameprop(:,:), &
 unitprop(:,:), &
 namepropphase(:), &
 modelpropphase(:)
real*8, pointer, save           :: &
 valueprop(:,:), &
 valuepropphase(:)
integer, parameter     :: &
 mxdim=50 
!-------------------------------------------------------------------------
!
!   $code
! 
!%----------------------------------------------------------------
havepph=present(pph)
 
select case (havepph)
case (.true.)
 numsp=isp
 allocate(species(numsp))
 do i=1,numsp
  call create_ (species(i))
  call set_name_ (species(i),namesp(i))
  call set_prop_ (species(i),valueprop(i,1:numprop(i)),nameprop(i,1:numprop(i)),unitprop(i,1:numprop(i)),numprop(i))
 end do
 call set_ (pph,species,numsp,ipropph,namephase,typephase,namepropphase(1:ipropph), &
    modelpropphase(1:ipropph),valuepropphase(1:ipropph),msg,iserror)
 if (iserror) goto 10 

 do i=1,numsp
  call destroy_ (species(i))
 end do
 deallocate(species)
 deallocate(numprop)
 deallocate(namesp)
 deallocate(nameprop)
 deallocate(unitprop)
 deallocate(namepropphase)
 deallocate(modelpropphase)
 deallocate(valueprop)
 deallocate(valuepropphase)
 
!%---------------------------------------------------------
case default
 
  select case (name)
!%---------------------------------------------------------
!% phase
!%---------------------------------------------------------
  case ('phase')
   isbesp=.false.
   isbeph=.true.
   allocate(numprop(mxdim))
   allocate(namesp(mxdim))
   allocate(nameprop(mxdim,mxdim))
   allocate(unitprop(mxdim,mxdim))
   allocate(namepropphase(mxdim))
   allocate(modelpropphase(mxdim))
   allocate(valueprop(mxdim,mxdim))
   allocate(valuepropphase(mxdim))
   isp=0
   iprop=0
   ipropph=0
   numprop=0
   id=''
   call get_value (attributes,"name", id, status)
   namephase=id
   id=''
   call get_value (attributes,"type", id, status)
   typephase=id
   id=''
   call get_value (attributes,"model", id, status)
   nameactmodel=id
   id=''
   call get_value (attributes,"convention", id, status)
   nameconvention=id
!%---------------------------------------------------------
!% species
!%---------------------------------------------------------  
  case ('species')
   isbesp=.true.
   isbeph=.false.
   isp=isp+1
   iprop=0
   id=''
   call get_value (attributes,"name", id, status)
   namesp(isp)=id
!%---------------------------------------------------------
!% properties
!%---------------------------------------------------------
  case ('property')
   if (isbesp) then
      iprop=iprop+1
      numprop(isp)=iprop
      id=''
      call get_value (attributes,"name", id, status)
      nameprop(isp,iprop)=id
      id=''
      call get_value (attributes,"unit", id, status)
      unitprop(isp,iprop)=id
      id=''
      call get_value (attributes,"value", id, status)
      n=0
      call build_data_array (id,reallocal,n)
      valueprop(isp,iprop)=reallocal(1)
   end if 
   if (isbeph) then
      ipropph=ipropph+1
      call get_value (attributes,"name", id, status)
      namepropphase(ipropph)=id
      call get_value (attributes,"model", id, status)
      modelpropphase(ipropph)=id
	  call get_value (attributes,"value", id, status)
      n=0
      call build_data_array (id,reallocal,n)
      valuepropphase(ipropph)=reallocal(1)
   end if 
  case default
 
   goto 10
 
  end select
 
end select
!%------------------------------------------------------------
return
 
10 print *,'Phase:'
print *,'Service: read_'
print *,'Error in tag name:', name
stop
 
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_parentphase