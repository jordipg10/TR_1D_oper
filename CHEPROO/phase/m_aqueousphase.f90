module m_aqueousphase
!-------------------------------------------------------------------------
!
!   $Description: Represent an aqueous phase and its thermodynamic behavior
!
!   $Use: use m_parentaqueousphase
! use m_aqueousphasedavis
! use m_aqueousphasepitzer
! use m_aqueousphasebdot
! use m_aqueousphaseideal
! use m_parentphase
! use m_species
! use flib_xpath
! use flib_sax
! use m_general_tools_cheproo
! use m_constants
!
!   $Author: Sergio Andr�s Bea Jofr� 
!
!   $License: UPC-CSIC
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_parentaqueousphase
use m_aqueousphasedavis
use m_aqueousphasepitzer
use m_aqueousphasebdot
use m_aqueousphaseideal
use m_parentphase
use m_species
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
implicit none
!%-------------------------------------------------------------
!%-------------------------------------------------------------
private    ::
!%-------------------------------------------------------------
!%-------------------------------------------------------------
public:: &
create_ &                    ! Create aqueous phase object. 
,destroy_ &                  ! Destroy aqueous phase object. 
,read_xml_ &                 ! Read the aqueous phase from xml file          
,set_ &                      ! Set the aqueous phase object 
,set_pparent_ &              ! Set the parent phase 
,get_pparent_ &              ! Return the pointer to parent phase 
,compute_density_ &          ! Compute density (and derivatives if asked)
,compute_sum_c_ &            ! Compute SUM(caqti)
,compute_mass_ &             ! Compute the total mass of solutes in the aqueous phase [kgr]
,compute_dmass_dc_ &         ! Compute the derivative of mass of solutes
,compute_act_coeff_ &        ! Compute activity coefficients vector 
,compute_dact_coeff_ &       ! Compute derivatives of the activity coefficients
,update_ &                   ! Update parameters in the aqueous phase that depends of the temperature
,write_xml_ &                ! Write in xml the attributes in the aqueous phase (not implemented)
,assignment(=) &             ! Copy an aqueous phase object in other aqueous phase object 
,write_                      ! Write in ascii the attributes encapsulated in the aqueous phase 
!%----------------------------------------------------
!% Constant parameters for specialization 
!%----------------------------------------------------
integer, parameter   :: &
aqphdavis = 1, &
aqphpz    = 2, &
aqphbdot  = 3, &
aqphideal = 4
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!% Type definition 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%-------------------------------------------------------------
type, public::t_aqueousphase
 
private                              ::
 
type (t_parentaqueousphase), pointer :: pp            ! Parent aqueous phase. 
 
type(t_aqueousphaseideal), pointer   :: paqphideal    ! Aqueous phase with ideal model. 
 
type(t_aqueousphasedavis), pointer   :: paqphdavis    ! Aqueous phase with Debye-H�ckel (extended) model.
 
type(t_aqueousphasepitzer), pointer  :: paqphpitzer   ! Aqueous phase with Pitzer model.
 
type(t_aqueousphasebdot), pointer    :: paqphbdot     ! Aqueous phase with Bdot model.
 
integer                              :: itype         ! Specialization indice.
 
end type t_aqueousphase
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%-------------------------------------------------------------
interface create_
 
module procedure create_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface set_
 
module procedure set_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface set_pparent_
 
module procedure set_pparent_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_density_
 
module procedure compute_density_aqph
 
end interface

!%----------------------------------------------------
!%----------------------------------------------------
interface get_pparent_
 
module procedure get_pparent_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface destroy_
 
module procedure destroy_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface read_xml_
 
module procedure read_xml_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_sum_c_
 
module procedure compute_sum_c_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_mass_
 
module procedure compute_mass_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_dmass_dc_
 
module procedure compute_dmass_dc_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_act_coeff_
 
module procedure compute_act_coeff_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_dact_coeff_
 
module procedure compute_dact_coeff_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface update_
 
module procedure update_temp_param_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface write_xml_
 
module procedure write_xml_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface write_
 
module procedure write_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface assignment(=)
 
module procedure copy_aqph
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
!%----------------------------------------------------
!%----------------------------------------------------
!%----------------------------------------------------
!%----------------------------------------------------
!%----------------------------------------------------
!%----------------------------------------------------
!%----------------------------------------------------
!%----------------------------------------------------
interface read_xml_loc_
 
module procedure read_xml_loc_aqph
 
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
subroutine create_aqph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create aqueous phase object
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(inout) :: this    ! Type aqueous phase variable
 
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
this%paqphdavis => null ()
this%paqphpitzer => null ()
this%paqphideal => null ()
this%paqphbdot => null ()
!%------------------------------------------------------------
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
subroutine destroy_aqph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy aqueous phase object.
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(inout) :: this   ! Type aqueous phase variable
 
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
!% If the parent aqueous phase is associated, then destroy it
!%------------------------------------------------------------
if (associated(this%pp)) then
 call destroy_ (this%pp)
 deallocate(this%pp)
end if 
this%pp => null ()
!%------------------------------------------------------------
!% Destroy the specialization 
!%------------------------------------------------------------
select case (this%itype)
case (aqphdavis)   
 call destroy_ (this%paqphdavis)
 deallocate (this%paqphdavis)
 this%paqphdavis => null ()
case (aqphpz)   
 call destroy_ (this%paqphpitzer)
 deallocate (this%paqphpitzer)
 this%paqphpitzer => null ()
case (aqphbdot)   
 call destroy_ (this%paqphbdot)
 deallocate (this%paqphbdot)
 this%paqphbdot => null ()
case (aqphideal) 
 call destroy_ (this%paqphideal)
 deallocate (this%paqphideal)
 this%paqphideal => null ()
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
subroutine read_xml_aqph &
   (this, &
    namefile, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read the aqueous phase from xml file     
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(inout)     :: this       ! Type aqueous phase variable

character(len=*), intent(in)            :: namefile   ! Name of the xml file

logical, intent(out)                    :: iserror    ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: The namefile must include the path 
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
 nameconvention, &
 msg, &
 namedatabase, &
 name  
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%----------------------------------------------------------------
!% Open the xml file 
!%----------------------------------------------------------------
call open_xmlfile(namefile, fxml, iostat)
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
! End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%----------------------------------------------------------------
!% Read the parent phase
!%----------------------------------------------------------------
call read_xml_ &
   (this%pp%pp, &
    namefile, &
    iserror)
if (iserror) goto 10
!%-------------------------------------------------------------
call read_xml_loc_ &
  (name, &
   attributes, &
   aqph=this, &
   iserror=iserror, &
   nameactmodel=nameactmodel, &
   nameconvention=nameconvention,&
   namedatabase=namedatabase)
!%-------------------------------------------------------------
!% Set the aqueous phase 
!%-------------------------------------------------------------
call set_ &
   (this, &
    this%pp%pp, &
    nameactmodel, &
    nameconvention, &
    iserror, &
    namedatabase) 
if (iserror) goto 10
!%------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Phase:'
print *,'Name:',this%pp%pp%name
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
subroutine set_aqph &
   (this, &
    pp, &
    nameactmodel, &
    nameconvention, &
    iserror, &
    namedatabase)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the aqueous phase object 
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(inout)        :: this              ! Type aqueous phase variable

type(t_parentphase), intent(in), target    :: pp                ! Type parent phase variable 

character(len=*), intent(in)               :: nameactmodel      ! Name of the thermodynamic model

character(len=*), intent(in)               :: nameconvention    ! Name of convention scaled

logical, intent(out)                       :: iserror           ! iserror=true, then there was an error

character(len=*), intent(in), optional     :: namedatabase      ! Name of thermodynamic database 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical                                    :: &
 havenamedatabase
character(len=300)                         :: &
 namedatabaseloc
character(len=100)                         :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Check the optional arguments
!%-----------------------------------------------------------
havenamedatabase=present(namedatabase)
!%-----------------------------------------------------------
! Set the parent phase in the parent aqueous phase
!%-----------------------------------------------------------
call set_ (this%pp,pp,iserror)
if (iserror) goto 10 
!%-----------------------------------------------------------
!% Create and set the specialization 
!%-----------------------------------------------------------
select case (nameactmodel)
case ('DHEXT', 'dhext', 'D-H-EXt', 'd-h-ext', 'davis', 'DAVIS','dh','DH')
   
   this%itype = aqphdavis
   allocate (this%paqphdavis)
   call create_ (this%paqphdavis)
   call set_ (this%paqphdavis,this%pp,iserror) 
   if (iserror) goto 10
   call update_ (this%paqphdavis,this%pp%pp%tempref,iserror) 
   if (iserror) goto 10
   
case ('Pitzer', 'PITZER', 'pitzer', 'pz', 'PZ')

   this%itype = aqphpz
   if (havenamedatabase.and.namedatabase/='') then
    namedatabaseloc=namedatabase
   else 
	namedatabaseloc='pitzer.xml'
   end if 
   allocate (this%paqphpitzer)
   call create_ (this%paqphpitzer)
   call set_ (this%paqphpitzer,this%pp,namedatabaseloc,nameconvention,iserror)
   if (iserror) goto 10
   call update_ (this%paqphpitzer,this%pp%pp%tempref,iserror) 
   if (iserror) goto 10
     
case ('TJ', 'tj', 'Truesdell-Jones', 'truesdell-jones')
   
   this%itype = aqphbdot  
   allocate (this%paqphbdot)
   call create_ (this%paqphbdot)
   call set_ (this%paqphbdot,this%pp,iserror)
   if (iserror) goto 10
   call update_ (this%paqphbdot,this%pp%pp%tempref,iserror) 
   if (iserror) goto 10

case ('ideal', 'IDEAL')
   
   this%itype = aqphideal
   allocate (this%paqphideal)
   call create_ (this%paqphideal)
   call set_ (this%paqphideal,this%pp,iserror)
   if (iserror) goto 10
   call update_ (this%paqphideal,this%pp%pp%tempref,iserror) 
   if (iserror) goto 10

case default 
 msg='Error, not recognized thermodynamic model:'
 call add_ (msg,nameactmodel)
 goto 10
end select
!%------------------------------------------------------------
return
 
10 continue 
print *,'*********************'
print *,'Phase:'
print *,'Name:', this%pp%pp%name
print *,'Service: specia_'
print *, msg
print *,'*********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_pparent_aqph &
   (this, &
    pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the parent phase 
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(inout)    :: this   ! Type aqueous phase variable

type(t_parentphase), intent(in), target:: pp     ! Type parent phase variable 
 
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
!% Set the parent phase in the parent aqueous phase 
!%------------------------------------------------------------
call set_pparent_ (this%pp,pp)
!%------------------------------------------------------------
!% Set the parent aqueous phase in the specialization 
!%------------------------------------------------------------
select case (this%itype)
 
case (aqphdavis)
 
 call set_pparent_ (this%paqphdavis,this%pp)
 
case (aqphpz)
 
 call set_pparent_ (this%paqphpitzer,this%pp)
 
case (aqphbdot)
 
 call set_pparent_ (this%paqphbdot,this%pp)
 
case (aqphideal)
 
 call set_pparent_ (this%paqphideal,this%pp)
 
end select
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_pparent_aqph &
   (this, &
    pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the parent phase pointer 
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(in)         :: this    ! Type aqueous phase variable

type(t_parentphase), pointer             :: pp      ! Type parent phase variable 
 
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
pp => this%pp%pp 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_aqph &
   (this, &
    temp, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update parameters in the aqueous phase that depends of the temperature
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(inout)  :: this    ! Type aqueous phase variable

real*8, intent(in)                   :: temp    ! Temperature in celcius 

logical, intent(out)                 :: iserror ! iserror=true, then there was an error
 
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
if (this%pp%pp%tempref==temp) return 
!%------------------------------------------------------------
this%pp%pp%tempref=temp 
!%------------------------------------------------------------
!% Update the parent aqueous phase
!%------------------------------------------------------------
call update_ (this%pp,temp,iserror)
if (iserror) goto 20  
!%------------------------------------------------------------
!% Update the specialization 
!%------------------------------------------------------------
select case (this%itype)
 
case (aqphdavis)
 
 call update_ (this%paqphdavis,temp,iserror)
 
case (aqphpz)
 
 call update_ (this%paqphpitzer,temp,iserror)
 
case (aqphbdot)
 
 call update_ (this%paqphbdot,temp,iserror)
 
case (aqphideal)
 
 call update_ (this%paqphideal,temp,iserror)
 
end select
!%------------------------------------------------------------
20 continue 
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
subroutine compute_act_coeff_aqph &
   (this, &
    g, &
    c, &
    iserror, &
    ionstr, &
	faccap)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute activity coefficients vector 
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(in)                  :: this    ! Type aqueous phase variable             

real*8, intent(in), dimension(this%pp%pp%numsp)   :: c       ! Concentration vector 

real*8, pointer, dimension(:)                     :: g       ! Activity coefficients vector 

real*8, intent(out), optional                     :: ionstr  ! Ionic strength 

logical, intent(out)                              :: iserror ! iserror=true, then there was an error

real*8, intent(in), optional                      :: faccap  ! Capillary correction for water activity 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical :: havefaccap 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false. 
!%----------------------------------------------------------- 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
havefaccap=present(faccap)
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
select case (this%itype)
 
case (aqphdavis)
 
 call compute_act_coeff_ (this%paqphdavis,g,c,iserror,ionstr)
 if(iserror) return 
 if (havefaccap.and.this%pp%ithw>0) then
   g(this%pp%ithw)=g(this%pp%ithw)*faccap 
 end if 
 
case (aqphpz)
 
 call compute_act_coeff_ (this%paqphpitzer,g,c,iserror,ionstr)
 if(iserror) return 
 if (havefaccap.and.this%pp%ithw>0) then
   g(this%pp%ithw)=g(this%pp%ithw)*faccap 
 end if 
 
case (aqphbdot)
 
 call compute_act_coeff_ (this%paqphbdot,g,c,iserror,ionstr)
 if(iserror) return 
 if (havefaccap.and.this%pp%ithw>0) then
   g(this%pp%ithw)=g(this%pp%ithw)*faccap 
 end if 
 
case (aqphideal)
 
 call compute_act_coeff_ (this%paqphideal,g,c,iserror,ionstr)
 
end select
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_mass_aqph &
   (this, &
    mass, &
    mol, &
    nsp, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the total mass of solutes in the aqueous phase [kgr]
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(in)    :: this     ! Type aqueous phase variable    

integer, intent(in)                 :: nsp      ! Number of species 

real*8, intent(in), dimension(nsp)  :: mol      ! Molality vector

real*8, intent(out)                 :: mass     ! Mass of solutes [kgr]

logical, intent(out)                :: iserror  ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: The concentrations must be expressed in molalities 
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
!% Call the corresponding service implemented in the parent
!% aqueou phase 
!%------------------------------------------------------------
call compute_mass_ (this%pp,mass,mol,nsp,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dmass_dc_aqph &
   (this, &
    dm, &
    dc)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Compute derivative of Mass of salt
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(in)              :: this        !  aqueous phase type

real*8, intent(in), dimension(:,:)             :: dc          ! Derivatives of the molalities.     

real*8, pointer, dimension(:)                  :: dm          ! Derivative of Mass

 
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
 
!%------------------------------------------------------------
!% Call the corresponding service implemented in the parent
!% aqueou phase 
!%------------------------------------------------------------
call compute_dmass_dc_(this%pp,dm,dc)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sum_c_aqph &
   (this, &
    sum, &
    c)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute SUM(caqti)
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(in)                  :: this  ! Type aqueous phase variable            

real*8, intent(in), dimension(this%pp%pp%numsp)   :: c     ! Concentrations vector 

real*8, intent(out)                               :: sum   ! SUM(caqti)
 
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
call compute_sum_c_ (this%pp,sum,c)
!%------------------------------------------------------------
return
10 continue 
print *,'Phase:'
print *,'Service:compute_sum_mass_c_'
print *,'Error'
stop
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dact_coeff_aqph &
   (this, &
    dg, &
    c, &
    dc, &
	iserror, &
	g, & 
	dtemp, &
    dionstr, &
	faccap)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute derivatives of the activity coefficients
!
!   $Arguments:
!
 
 
type(t_aqueousphase), intent(in)                             :: this    ! Type aqueous phase variable  

real*8, intent(in), dimension(this%pp%pp%numsp)              :: c       ! Concentration vector

real*8, intent(in), dimension(:,:)                           :: dc      ! Derivatives of the concentrations 

real*8, pointer, dimension(:,:)                              :: dg      ! Derivatives of the activity coefficients 

logical, intent(out)                                         :: iserror ! iserror=true, then there was an error

real*8, pointer, optional, dimension(:)                      :: dionstr ! Derivatives of the ionic strength 

real*8, intent(in), optional                                 :: faccap  ! Capillary correction for water activity 

real*8, intent(in), optional, dimension(:)                   :: dtemp   ! Derivatives of the ionic strength 

real*8, intent(in), optional, dimension(this%pp%pp%numsp)    :: g       ! Activity coefficients vector 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical         :: havefaccap 
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false. 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%-------------------------------------------------------------
havefaccap=present(faccap)
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%-------------------------------------------------------------
select case (this%itype)
 
case (aqphdavis)
 
 call compute_dact_coeff_ &
   (this%paqphdavis, &
    dg, &
    c, &
    dc, &
	iserror, & 
	g=g, &
	dtemp=dtemp, &
    dionstr=dionstr)
 
 
case (aqphpz)
 
 call compute_dact_coeff_ &
   (this%paqphpitzer, &
    dg, &
    c, &
    dc, &
	iserror, &
	g=g, &
	dtemp=dtemp, &
    dionstr=dionstr)
 
case (aqphbdot)
 
 call compute_dact_coeff_ &
   (this%paqphbdot, &
    dg, &
    c, &
    dc, &
	iserror, &
    dionstr=dionstr)
 
case (aqphideal)
 
 call compute_dact_coeff_ &
   (this%paqphideal, &
    dg, &
    c, &
    dc, &
    dionstr=dionstr)
 
end select
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
if(iserror) return 
!%------------------------------------------------------------
!% Multiply by capillary correction on water activity 
!%------------------------------------------------------------
if (havefaccap.and.this%pp%ithw>0) then
   dg(this%pp%ithw,:)=dg(this%pp%ithw,:)*faccap 
end if 
!%------------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_xml_aqph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in xml the attributes in the aqueous phase (not implemented)
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(in) :: this   ! Type aqueous phase variable  
 
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
 

 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_aqph &
   (this, &
    ioutput, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii the attributes encapsulated in the aqueous phase 
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(in) :: this    ! Type aqueous phase variable  

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
character(len=100)   :: &
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
 
case (aqphdavis)
 
 call write_ (this%paqphdavis,ioutput,iserror)
 
case (aqphpz)
 
 call write_ (this%paqphpitzer,ioutput,iserror)
 
case (aqphbdot)
 
 call write_ (this%paqphbdot,ioutput,iserror) 
 
case (aqphideal)
 
 call write_ (this%paqphideal,ioutput,iserror)
 
end select
!%------------------------------------------------------------
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
subroutine copy_aqph &
  (copied, &
   this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy a object in other species object.  
!
!   $Arguments:
!
 
type(t_aqueousphase), intent(in) :: this    ! Type aqueous phase variable  

type(t_aqueousphase), intent(out):: copied  ! Type aqueous phase variable  
 
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
!% Copy the parent aqueous phase 
!%------------------------------------------------------------
copied%pp = this%pp
!%------------------------------------------------------------
!% Create and copy the specialization 
!%------------------------------------------------------------
select case (this%itype)
 
case (aqphdavis)
 
 allocate(copied%paqphdavis)
 call create_ (copied%paqphdavis)
 copied%paqphdavis=this%paqphdavis
 call set_pparent_ (copied%paqphdavis,copied%pp)
 
case (aqphpz)
 
 allocate(copied%paqphpitzer)
 call create_ (copied%paqphpitzer)
 copied%paqphpitzer=this%paqphpitzer
 call set_pparent_ (copied%paqphpitzer,copied%pp)
 
case (aqphbdot)
 
 allocate(copied%paqphbdot)
 call create_ (copied%paqphbdot)
 copied%paqphbdot=this%paqphbdot
 call set_pparent_ (copied%paqphbdot,copied%pp)
 
 
case (aqphideal)
 
 allocate(copied%paqphideal)
 call create_ (copied%paqphideal)
 call set_pparent_ (copied%paqphideal,copied%pp)
 
 
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
 
call read_xml_loc_ (name,attributes)
 
return
end subroutine begin_element_handler
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc_aqph &
  (name, &
   attributes, &
   aqph, &
   iserror, &
   nameactmodel, &
   nameconvention, &
   namedatabase)
implicit none 
!-------------------------------------------------------------------------
!
!   $Description: Read phase from xml file
!
!   $Arguments:
!
 
type(dictionary_t), intent(in)  :: attributes

type(t_aqueousphase), optional  :: aqph            ! Type aqueous phase variable  

logical, intent(out), optional  :: iserror         ! iserror=true, then there was an error

character(len=*), optional      :: nameconvention  ! Name of convention scaled 

character(len=*), optional      :: nameactmodel    ! Name of the thermodynamic model

character(len=*), optional      :: namedatabase    ! Name of the thermodynamic database 
 
character(len=*)                :: name            ! Name of the phase           
 
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
 haveaqph, &
 haveiserror 
character(len=100)     :: &
 id
character(len=100),save:: &
 namemodel, &
 nameconv, &
 namebase
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
haveaqph=present(aqph)
haveiserror=present(iserror)
!%-----------------------------------------------------------
if (haveiserror) iserror=.false. 
!%-----------------------------------------------------------
select case (haveaqph)
case (.true.)
  
 nameactmodel=namemodel
 nameconvention=nameconv 
 namedatabase='pitzer.xml'
 
case default
 
  select case (name)
 
  case ('phase')
    namebase=' ' 
    id='' 
	call get_value (attributes,"model", id, status)
    namemodel=id
    id='' 
	call get_value (attributes,"convention", id, status)
    nameconv=id
  end select  


end select
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_density_aqph &
    (this        ,&
     pliq        ,&
     temp        ,&
     c           ,&
     density     ,&
     ismolality  ,&
     iserror     ,&
     dc          ,&
     ddensitydc)
implicit none      
!-------------------------------------------------------------------------
!
!   $Description: Compute density from concentrations 
!
!   $Arguments:
!

type(t_aqueousphase), intent(in)               :: this         !Type aqueous phase 

real*8,intent(in)                              :: pliq         !liquid pressure

real*8,intent(in)                              :: temp         !temperature

real*8,dimension(this%pp%pp%numsp), intent(in) :: c            !concentration array

real*8,intent(out)                             :: density      !Density     

logical,intent(in)                             :: ismolality   !

logical,intent(out)                            :: IsError      !iserror=true, then there was an error

real*8,dimension(:,:),optional                 :: dc           !concentration derivatives (needed to calculate density derivatives)

real*8,pointer,dimension(:),optional           :: ddensitydc   !Density derivative

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                               :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
msg=' ' 
iserror=.false. 
!%------------------------------------------------------------
if ( present(ddensitydc)) then
   
    if (.not.present(dc)) then
        msg="dc are needed for calcutate density derivatives"
        isError=.true.
        goto 10
    endif
    
!%------------------------------------------------------------
!% Select the specialization 
!%------------------------------------------------------------
    select case (this%itype)
 
    case (aqphpz)
 
       call compute_density_(this%paqphpitzer,c,density,ismolality,iserror)
  
    case default
 
       call compute_density_(this%pp,pliq,temp,c,density,iserror,dc,ddensitydc) 
 
    end select
   
else
    
!%------------------------------------------------------------
!% Select the specialization 
!%------------------------------------------------------------
    select case (this%itype)
 
    case (aqphpz)
 
       call compute_density_(this%paqphpitzer,c,density,ismolality,iserror)
  
    case default
 
       call compute_density_(this%pp,pliq,temp,c,density,iserror) 
 
    end select
    
endif
!%----------------------------------------------------------
if (iserror) goto 10 
!%----------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'AqueousPhase:'
print *,'Service: compute_density_'
print *, msg
print *,'***********************'
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
end module m_aqueousphase