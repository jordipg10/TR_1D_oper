module m_phase
!-------------------------------------------------------------------------
!
!   $Description: The phase class describes a phase and its thermodynamic behavior. 
! It represents the 'non-ideality' of the chemical system. 
!Phases can be defined as homogeneous parts of the chemical system having distinct 
!boundaries with adjacent phases and being mechanically separable (an exception is 
! made for the surface phases).  
!The phase class contains the species class and the attributes that depend on the nature 
!of the phase (e.g.  
!the dielectric constant for the aqueous phase). 
!The main function of the phase class is to compute the thermodynamic activity coefficients 
!of involved species.
!The thermodynamic models also depend on the nature of the phase. 
!For example, both diluted water and brine are specializations of 'aqueous phase' but their 
!thermodynamic behavior 
!is different. 
!The thermodynamic behavior of dilute water can be approximated as ideal or according to 
!Debye-H�ckel-based formulas 
!(Davis and Bdot Classes, Debye and H�ckel, 1923). 
!Brines behave according to the ionic interaction model (Pitzer Class).
!A similar analysis can be followed for mineral and gas phases. 
!The thermodynamic behavior of solid and gas phases can be approximated as a pure phase 
!or as an ideal or non-ideal 
!solution.
!Different approaches to treat non-ideal solid solutions (i.e. regular, sub regular) and 
!a different equations of 
!state for non-ideal pure gases and gas mixtures can be incorporated as separate classes.
!
!   $Use: use m_parentphase
! use m_aqueousphase
! use m_mineralphase
! use m_gasphase
! use m_species
! use flib_xpath
! use flib_sax
! use m_general_tools_cheproo
! use m_constants
!
!   $Author: Sergio Andr�s Bea Jofr�
!
!   $License: CSIC-UPC
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_parentphase
use m_aqueousphase
use m_mineralphase
use m_gasphase
use m_species
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%-----------------------------------------------------------
!%-----------------------------------------------------------
private             ::
!%-----------------------------------------------------------
!% Public services 
!%-----------------------------------------------------------
public              :: &
create_ &                        ! Create phase object.
,read_xml_ &                     ! Read aqueous phase from xml file. 
,set_ &                          ! Set attributes in the phase object                
,destroy_ &                      ! Destroy phase object.
,compute_sum_c_ &                ! Copute SUM(c)
,compute_mass_ &                 ! Compute total mass in the phase. 
,compute_dmass_dc_ &             ! Compute the derivative of mass in the phase
,compute_act_coeff_ &            ! Compute activity coefficients
,compute_dact_coeff_ &           ! Compute the derivatives of the activity coefficients.  
,compute_charge_balance_ &       ! Compute the charge balance. 
,compute_dcharge_balance_ &      ! Compute derivatives of the charge balance. 
,compute_property_  &            ! Compute property of the phase  
,compute_density_ &              ! Compute density (and derivatives if asked)
,update_ &                       ! Update paramaters that depends of the temperature in the phase
,get_numsp_ &                    ! Return the number of species in the phase. 
,get_namesp_ &                   ! Return the name of the species 
,get_name_ &                     ! Return the name of the phase 
,get_prop_ &                     ! Return some property according its name
,get_if_aqueous_ &               ! Return if the phase is aqueous 
,get_if_gas_ &                   ! Return if the phase is a gas
,get_if_mineral_ &               ! Return if the phase is a mineral 
,get_if_sp_is_present_ &         ! Return if the species is present in the phase. 
,write_xml_ &                    ! Write in xml file the attributes encapsulated in the phase object (not implemented).  
,write_ &                        ! Write in ascii file the attributes encapsulated in the phase object. 
,verify_ &                       ! Verify attributes encapsulated in the object. 
,assignment(=) &                 ! Copy a phase object in other phase object. 
,get_pspecies_ &                 ! Return the pointer of ith species 
,test_cpu_time_act_coeff_        ! Compute the average cpu time consumed when the activity coefficients are computed. Also return the activity coefficients vector. 
!%-----------------------------------------------------------
!% Private services 
!%-----------------------------------------------------------
private              :: &
read_xml_loc &
,begin_element_handler
!%---------------------------------------------------
!% Constant parameters
!%----------------------------------------------------
integer, parameter   :: &
aqph  = 1, &     ! Aqueous phase specialization 
minph = 2, &     ! Mineral phase specialization 
gasph = 3        ! Gas phase specialization 
!%------------------------------------------------------------
!%------------------------------------------------------------
!% Type pointer to phase object 
!%------------------------------------------------------------
!%------------------------------------------------------------
type, public::t_pphase

type(t_phase), pointer::ptr

end type
!%------------------------------------------------------------
!%------------------------------------------------------------
!% Type definition 
!%------------------------------------------------------------
!%------------------------------------------------------------
type, public::t_phase
 
private                           ::
 
type (t_parentphase), pointer     :: pp       ! Type parent phase 
 
type(t_aqueousphase), pointer     :: paqph    ! Type aqueous phase 
 
type(t_mineralphase), pointer     :: pminph   ! Type mineral phase 

type(t_gasphase), pointer         :: pgasph   ! Type gas phase 
 
integer                           :: itype    ! Itype index 
 
 
end type t_phase
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
interface create_
 
module procedure create_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface destroy_
 
module procedure destroy_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface read_xml_
 
module procedure read_xml_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface set_
 
module procedure set_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_sum_c_
 
module procedure compute_sum_c_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_mass_
 
module procedure compute_mass_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_dmass_dc_
 
module procedure compute_dmass_dc_ph
 
end interface

!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_act_coeff_
 
module procedure compute_act_coeff_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_dact_coeff_
 
module procedure compute_dact_coeff_ph
module procedure compute_dact_coeff_numeric_ph 
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_charge_balance_
 
module procedure compute_charge_balance_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_property_
 
module procedure compute_property_ph 
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_dcharge_balance_
 
module procedure compute_dcharge_balance_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface update_
 
module procedure update_temp_param_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface get_numsp_
 
module procedure get_numsp_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface get_pspecies_
 
module procedure get_pspecies_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface get_namesp_
 
module procedure get_namesp_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface get_name_
 
module procedure get_name_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface get_prop_
 
module procedure get_prop_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface get_if_aqueous_
 
module procedure get_if_aqueous_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface get_if_gas_
 
module procedure get_if_gas_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface get_if_mineral_
 
module procedure get_if_mineral_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface get_if_sp_is_present_
 
module procedure get_if_sp_is_present_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface write_xml_
 
module procedure write_xml_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface write_
 
module procedure write_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface verify_
 
module procedure verify_ph
module procedure verify1_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface test_cpu_time_act_coeff_
 
module procedure test_cpu_time_act_coeff_ph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface compute_density_
 
module procedure compute_density_ph
 
end interface

!%------------------------------------------------------------
!%------------------------------------------------------------
interface assignment(=)
 
module procedure copy_ph

end interface


!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
interface read_xml_loc
 
module procedure read_xml_loc 
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
interface begin_element_handler
 
module procedure begin_element_handler
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
contains
!%------------------------------------------------------------
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_ph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create the phase object. 
!
!   $Arguments:
!
 
type(t_phase), intent(inout) :: this     ! Type phase variable 
 
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
this%itype=0
!%-----------------------------------------------------------
this%pp => null ()
this%paqph => null ()
this%pminph => null ()
this%pgasph => null ()
!%-----------------------------------------------------------
!% Allocate and create parent phase 
!%-----------------------------------------------------------
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
subroutine destroy_ph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy the phase object
!
!   $Arguments:
!
 
 
type(t_phase), intent(inout) :: this  ! Type phase variable 
 
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
!% If the parent phase was created, then destroy it and 
!% nullify pointer 
!%------------------------------------------------------------
if (associated(this%pp)) then
 call destroy_ (this%pp)
 deallocate (this%pp)
end if
this%pp => null ()
!%------------------------------------------------------------
!% Destroy and deallocate the specialization 
!%------------------------------------------------------------
select case (this%itype)
case (aqph)
 call destroy_ (this%paqph)
 deallocate(this%paqph)
 this%paqph => null ()
case (minph)
 call destroy_ (this%pminph)
 deallocate(this%pminph)
 this%pminph => null ()
case (gasph)
 call destroy_ (this%pgasph)
 deallocate(this%pgasph)
 this%pgasph => null ()
end select
!%------------------------------------------------------------
this%itype=0
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_ph &
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
 
type(t_phase), intent(inout):: this       ! Type phase variable 

character(len=*), intent(in):: namefile   ! Name of the xml file

logical, intent(out)        :: iserror    ! iserror=true, then there was an error
 
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
 name, &
 typephase, &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------------
! Open xml file
!%--------------------------------------------------------------------
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
! End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%--------------------------------------------------------------------
! Read phase from xml file 
!%--------------------------------------------------------------------
call read_xml_loc(name,attributes,phase=this,iserror=iserror,typephase=typephase)
if (iserror) goto 10 
!%--------------------------------------------------------------------
!% If the parent phase is associated then destroy and deallocate
!%--------------------------------------------------------------------
if (associated(this%pp)) then
 call destroy_ (this%pp)
 deallocate(this%pp)
end if 
this%pp => null ()
!%--------------------------------------------------------------------
!% Select phase type according physical properties and create the 
!% specialization  
!%--------------------------------------------------------------------
select case (typephase)
case ('aqueous','AQUEOUS')
 
 this%itype=aqph
 allocate (this%paqph)
 call create_ (this%paqph) 
 call read_xml_(this%paqph,namefile,iserror)
 call get_pparent_ (this%paqph,this%pp)

case ('mineral','MINERAL')
 
 this%itype=minph
 allocate (this%pminph)
 call create_ (this%pminph) 
 call read_xml_ (this%pminph,namefile,iserror)
 call get_pparent_ (this%pminph,this%pp)

case ('gas','GAS')
 
 this%itype=gasph
 allocate (this%pgasph)
 call create_ (this%pgasph) 
 call read_xml_ (this%pgasph,namefile,iserror)
 call get_pparent_ (this%pgasph,this%pp)

case default 
 msg='Error, not recognized phase type: '
 call add_ (msg,typephase) 
 goto 10  
end select 
!%------------------------------------------------------------
return
 
10 continue 
print *,'*******************'
print *,'Phase:'
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
subroutine set_ph &
   (this, &
    species, &
    numsp, &
    numprop, &
    nameph, &
    phasetype, &
    nameactmodel, &
    namepropph, &
    modelpropph, &
    propph, &
    nameconvention, &
    iserror, &
    namedatabase)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set attributes in the phase object 
!
!   $Arguments:
!
 
type(t_phase), intent(inout)                              :: this             ! Type phase variable 

integer, intent(in)                                       :: numsp            ! Number of species 

integer, intent(in)                                       :: numprop          ! Number of properties 

type(t_species), intent(in), dimension(numsp), target     :: species          ! List of the 

character(len=*), intent(in)                              :: nameph           ! Name of the phase

character(len=*), intent(in)                              :: nameactmodel     ! Name of the thermodynamic model

character(len=*), intent(in)                              :: phasetype        ! Type of phase 

character(len=*), intent(in)                              :: nameconvention   ! Name of the convention scaled (only for MacInnes scaled in the Pitzer specialization)

real*8, intent(in), dimension(numprop)                    :: propph           ! Values of the properties [numprop]

character(len=*), intent(in), dimension(numprop)          :: namepropph       ! Name of the properties [numprop]

character(len=*), intent(in), dimension(numprop)          :: modelpropph      ! Name of the model properties [numprop]

character(len=*), intent(in), optional                    :: namedatabase     ! Name of the thermodynamic database (e.g. virial coefficients database)

logical, intent(out)                                      :: iserror          ! iserror=true, then there was an error
 
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
!%--------------------------------------------------------------------
! Set the parent phase
!%--------------------------------------------------------------------
call set_ &
   (this%pp, &
    species, &
    numsp, &
    numprop, &
    nameph, &
    phasetype, &
    namepropph, &
    modelpropph, &
    propph, &
    msg, &
    iserror)
if (iserror) goto 10
!%--------------------------------------------------------------------
! Check and create phase type 
!%--------------------------------------------------------------------
select case (phasetype)
case ('aq','AQ','Aqueous','AQUEOUS','aqueous')
 
 this%itype=aqph
 allocate(this%paqph)
 call create_ (this%paqph)
 call set_ (this%paqph,this%pp,nameactmodel,nameconvention,iserror,namedatabase=namedatabase) 

case ('mineral','MINERAL','Mineral','MIN','min')

 this%itype=minph
 allocate(this%pminph)
 call create_ (this%pminph)
 call set_ (this%pminph,this%pp,nameactmodel,iserror) 

case ('gas','GAS','Gas')

 this%itype=gasph
 allocate(this%pgasph)
 call create_ (this%pgasph)
 call set_ (this%pgasph,this%pp,nameactmodel,iserror) 

case default

 msg='Error, not recognized phase type:'
 call add_ (msg,phasetype)
 goto 10 

end select 
!%---------------------------------------------------------------------
 
return
 
10 continue 
print *,'************************'
print *,'Phase:'
print *,'Name:',this%pp%name
print *,'Service: set_'
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
subroutine update_temp_param_ph &
   (this, &
    temp, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update paramaters that depends of the temperature in the phase. 
!
!   $Arguments:
!

type(t_phase), intent(inout)      :: this      ! Type phase variable 

real*8, intent(in)                :: temp      ! Temperature in celsius 

logical, intent(out)              :: iserror   ! iserror=true, then there was an error 
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
 
case (aqph)
 
 call update_ (this%paqph,temp,iserror)
 
case (minph)
 
 call update_ (this%pminph,temp,iserror)

case (gasph)
 
 call update_ (this%pgasph,temp,iserror)
 
end select
!%-----------------------------------------------------------
return
10 continue 
print *, '******************************'
print *, 'Phase:'
print *, 'Name:',this%pp%name
print *, 'Service: update_ (temperature)'
print *,  msg
print *, '******************************'
iserror=.true.
return 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_act_coeff_ph &
   (this, &
    g, &
    c, &
    iserror, & 
    param, &
	factor)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute activity coefficients.
!
!   $Arguments:
!
 
type(t_phase), intent(in)                      :: this             ! Type phase variable 

real*8, intent(in), dimension(this%pp%numsp)   :: c                ! Vector of concentrations 

real*8, pointer, dimension(:)                  :: g                ! Vector of activity coefficients

real*8, intent(out), optional                  :: param            ! May be ionic strength. 

real*8, intent(in), optional                   :: factor           ! Correction factor 

logical, intent(out)                           :: iserror          ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
msg=''
iserror=.false. 
!%-----------------------------------------------------------
select case (this%itype)
 
case (aqph)
 
 call compute_act_coeff_ (this%paqph,g,c,iserror,ionstr=param,faccap=factor)
 
case (minph)
 
 call compute_act_coeff_(this%pminph,g,c,iserror)

case (gasph)
 
 call compute_act_coeff_(this%pgasph,g,c,iserror)
 
end select
!%------------------------------------------------------------
return
10 continue 
print *, '******************************'
print *, 'Phase:'
print *, 'Name:',this%pp%name
print *, 'Service: compute_act_coeff_'
print *,  msg
print *, '******************************'
return 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sum_c_ph &
   (this, &
    sumc, &
    c)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copute SUM(c)
!
!   $Arguments:
!
 
type(t_phase), intent(in)                            :: this   ! Type phase variable 

real*8, intent(out)                                  :: sumc   ! Copute SUM(c)

real*8, intent(in), dimension(this%pp%numsp)         :: c      ! Concentration vector 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical                    :: &
iserror 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
select case (this%itype)
case (aqph)
 call compute_sum_c_ (this%paqph,sumc,c)
case default 
 call compute_sum_c_ (this%pp,sumc,c,iserror)
end select
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_mass_ph &
   (this, &
    mass, &
    mol, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the total mass in the phase from molality of the 
! species (actually only implmented for aqueous phase).
!
!   $Arguments:
!
 
type(t_phase), intent(in)            :: this     ! Type phase variable 

integer, intent(in)                  :: nsp      ! Number of species. 

real*8, intent(out)                  :: mass     ! Mass [kg]

real*8, intent(in), dimension(nsp)   :: mol      ! Molality vector

logical, intent(out)                 :: iserror  ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100):: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Call the specific service in the specialization 
!%------------------------------------------------------------
select case (this%itype)

case (aqph)

 call compute_mass_(this%paqph,mass,mol,nsp,iserror)

case default

 msg='Error, not implemented service for this specialization'

 goto 10

end select
!%-----------------------------------------------------------
return
10 continue 
print *, '*************************'
print *, 'Phase:'
print *, 'Name:',this%pp%name
print *, 'Service: compute_mass_'
print *,  msg
print *, '*************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
  subroutine compute_dmass_dc_ph &
   (this, &
    dm, &
    dc,iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Compute derivative of Mass of salt
!
!   $Arguments:
!
 
type(t_phase), intent(in)                      :: this       !  phase type

real*8, intent(in), dimension(:,:)             :: dc          ! Derivatives of the molalities.     

real*8, pointer, dimension(:)                  :: dm          ! Derivative of Mass

logical                                        :: iserror    ! iserror=true, then there was an error 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100):: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
isError=.false.
!%------------------------------------------------------------
!% Call the specific service in the specialization 
!%------------------------------------------------------------
select case (this%itype)

case (aqph)

 call compute_dmass_dc_(this%paqph,dm,dc)

case default

 msg='Error, not implemented service for this specialization'

 goto 10

end select
!%-----------------------------------------------------------
return
10 continue 
print *, '*************************'
print *, 'Phase:'
print *, 'Name:',this%pp%name
print *, 'Service: compute_dmass_dc_'
print *,  msg
print *, '*************************'
iserror=.true.
return
end subroutine



!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine test_cpu_time_act_coeff_ph &
 (this, &
  cputime, &
  g, &
  c, &
  ntime, &
  iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the average cpu time consumed when the activity 
! coefficients are computed. Also return the activity coefficients vector. 
! 
!
!   $Arguments:
!
 
type (t_phase), intent(in)                    :: this          ! Type phase variable 

integer, intent(in)                           :: ntime         ! Number of times 

real*8, intent(out)                           :: cputime       ! Average cpu time

real*8, intent(in), dimension(this%pp%numsp)  :: c             ! Concentration vector

real*8, pointer, dimension(:)                 :: g             ! Activity coefficients vector

logical, intent(out)                          :: iserror       ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                           :: &
 msg
integer                                      :: &
 i
real*8                                       :: &
 cputime1, &
 cputime2 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%------------------------------------------------------------
call cpu_time (cputime1)
!%------------------------------------------------------------
do i=1,ntime
 
 call compute_act_coeff_ (this,g,c,iserror)
 if (iserror) goto 10  

end do
!%------------------------------------------------------------
call cpu_time (cputime2)
!%------------------------------------------------------------
cputime=(cputime2-cputime1)/real(ntime)
!%------------------------------------------------------------
return
10 iserror=.true.
return 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dact_coeff_ph &
   (this, &
    dg, &
    c, &
    dc, &
    iserror, &
    g, &
	dtemp, &
    dparam, &
	factor)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the derivates of the activity coefficients.  
!
!   $Arguments:
!
 
type(t_phase), intent(in)                               :: this       ! Type phase variable 

real*8, pointer, dimension(:,:)                         :: dg         ! Derivatives of activity coefficients

real*8, intent(in), dimension(this%pp%numsp)            :: c          ! Molality vector. 

real*8, intent(in), dimension(:,:)                      :: dc         ! Derivatives of concentration species. 

logical, intent(out)                                    :: iserror    ! iserror=true, then there was an error

real*8, intent(in), optional, dimension(this%pp%numsp)  :: g          ! Activity coefficients vector

real*8, intent(in), optional, dimension(:)              :: dtemp      ! Derivatives of the temperature. 

real*8, intent(in), optional                            :: factor     ! Correction factor for activity coefficients 

real*8, pointer, optional, dimension(:)                 :: dparam     ! Derivatives a someone parameters (e.g. derivatives of the ionic strength)
 
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
iserror=.false.
!%---------------------------------------------------------
!% Call the specialization 
!%---------------------------------------------------------
select case (this%itype)
 
case (aqph)
 
 call compute_dact_coeff_ &
   (this%paqph, &
    dg, &
    c, &
    dc, &
	iserror, &
	g=g, &
	dtemp=dtemp, &
    dionstr=dparam, &
	faccap=factor)
 
case (minph)
 
 call compute_dact_coeff_ &
   (this%pminph, &
    dg, &
    c, &
    dc, &
	iserror, &
    dparam=dparam)

case (gasph)
 
 call compute_dact_coeff_ &
   (this%pgasph, &
    dg, &
    c, &
    dc, &
	iserror, &
    dparam=dparam)
 
end select
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dact_coeff_numeric_ph &
   (this, &
    dg, &
    c, &
    dc, &
    pert, &
    iserror, &
    dparam, &
	dtemp, &
	factor)
    
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute numerical analytic derivatives 
!
!   $Arguments:
!
 
type(t_phase), intent(in)                               :: this      ! Type phase variable 

real*8, intent(in), dimension(this%pp%numsp)            :: c         ! Concentration vector 

real*8, intent(in), dimension(:,:)                      :: dc        ! Derivatives of the concentrations 

real*8, pointer, dimension(:,:)                         :: dg        ! Derivatives of the activity coefficients 

real*8, intent(in)                                      :: pert      ! Perturbation factor

logical, intent(out)                                    :: iserror   ! iserror=true, then there was an error

real*8, pointer, optional, dimension(:)                 :: dparam    ! Derivatives of some parameters 

real*8, intent(in), optional                            :: factor    ! Correction factor 

real*8, intent(in), optional, dimension(:)              :: dtemp     ! Derivatives of the temperature 
 
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
 isps, &
 ndimder
real*8             :: &
 delta 
logical            :: &
 havedparam, &
 havedtemp  
real*8, pointer    :: &
 cpert(:) => null (), &
 dgloc(:,:) => null (), &
 g(:) => null (), &
 gpert(:) => null () 
character(len=100)   :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%---------------------------------------------------------
iserror=.false.
msg= ' ' 
!%---------------------------------------------------------
!% Check optional arguments 
!%---------------------------------------------------------
havedparam=present(dparam) 
havedtemp=present(dtemp) 
!%---------------------------------------------------------
!% Allocate pointers 
!%---------------------------------------------------------
ndimder=size(dc,2)
if (havedparam) then
 call check_pointer_ (dparam,ndimder,.true.)
end if
call check_pointer_ (dg,this%pp%numsp,ndimder,.true.)
!%---------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------
call check_pointer_ (dgloc,this%pp%numsp,this%pp%numsp,.true.)
call check_pointer_ (cpert,this%pp%numsp,.true.)
!%---------------------------------------------------------
!% Compute non perturbed activity coefficients 
!%---------------------------------------------------------
call compute_act_coeff_(this,g,c,iserror,factor=factor)
!%---------------------------------------------------------
if (iserror) goto 20 
!%---------------------------------------------------------
do isps=1,this%pp%numsp 
 
 cpert=c
 delta=cpert(isps)*pert
 cpert(isps)=cpert(isps)+delta
 
 call compute_act_coeff_(this,gpert,cpert,iserror,factor=factor)
 
 if (iserror) goto 20 
 
 dgloc(:,isps)=(gpert-g)/delta
 
end do
!%---------------------------------------------------------
!% Temperature ????
!%---------------------------------------------------------
if (havedtemp) then
 
end if
!%---------------------------------------------------------
!% 
!%---------------------------------------------------------
dg=matmul(dgloc,dc)
!%---------------------------------------------------------
20 continue 
!%---------------------------------------------------------
!% Deallocate local pointers 
!%---------------------------------------------------------
call check_pointer_ (dgloc,1,1,.false.)
call check_pointer_ (cpert,1,.false.)
call check_pointer_ (g,1,.false.)
call check_pointer_ (gpert,1,.false.)
!%---------------------------------------------------------
if (iserror) goto 10 
!%---------------------------------------------------------
return
10 continue 
print *,'****************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_dact_coeff_'
print *, msg
print *,'****************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_numsp_ph &
   (this, &
    numsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Return the number of species in the phase. 
!
!   $Arguments:
!
 
type(t_phase), intent(in)  :: this     ! Type phase variable 

integer, intent(out)       :: numsp    ! Number of species 
 
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
!% Call the corresponding service in the parent phase 
!%------------------------------------------------------------
call get_numsp_ (this%pp,numsp)
!%------------------------------------------------------------ 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_pspecies_ph &
   (this, &
    pspecies, &
    ithsp, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the pointer of ith species 
!
!   $Arguments:
!

type(t_phase), intent(in)  :: this       ! Type phase variable 

integer, intent(in)        :: ithsp      ! Local index of the species

type(t_species), pointer   :: pspecies   ! Pointer to species 

logical, intent(out)       :: iserror    ! iserror=true, then there was an error
 
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
call get_pspecies_ (this%pp,pspecies,ithsp,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_namesp_ph &
   (this, &
    namesp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Get the name of the phase
!
!   $Arguments:
!
 
type(t_phase), intent(in)                     :: this         ! Type phase variable 

character(len=*), pointer, dimension(:)       :: namesp       ! Name of the species 
 
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
call get_namesp_ (this%pp,namesp)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_if_aqueous_ph &
   (this, &
    isaqueous)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_phase), intent(in) :: this       ! Type phase variable 

logical, intent(out)      :: isaqueous  ! isaqueous=true, the phase is aqueous
 
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

isaqueous = .false.
!%----------------------------------------------------------- 
select case (this%itype)
case (aqph)
 isaqueous = .true.
end select
!%----------------------------------------------------------- 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_property_ph &
   (this, &
    prop, &
	nameprop, &
	c, &
	temp, &
	pressure, &
	nsp, &
	iserror)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute property of the phase 
!
!   $Arguments:
!
 
type(t_phase), intent(in)             :: this        ! Type phase variable 

real*8, intent(out)                   :: prop        ! Value of the property 

character(len=*), intent(in)          :: nameprop    ! Name of the property 

integer, intent(in)                   :: nsp         ! Number of species

real*8, intent(in), dimension(nsp)    :: c           ! Concentrations

real*8, intent(in)                    :: temp        ! Temperature

real*8, intent(in)                    :: pressure    ! Pressure

logical, intent(out)                  :: iserror     ! iserror=true, then there was an error
 
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
prop=prop*1.0d0
iserror=.false.

!%----------------------------------------------------------- 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_if_gas_ph &
   (this, &
    isgas)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_phase), intent(in) :: this    ! Type phase variable 

logical, intent(out)      :: isgas   ! isgas=true, the phase is a gas 
 
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
 

 
isgas = .false.
!%----------------------------------------------------------- 
select case (this%itype)
case (gasph)
 isgas = .true.
end select
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_if_mineral_ph &
   (this, &
    ismin)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_phase), intent(in) :: this    ! Type phase variable 

logical, intent(out)      :: ismin   ! ismin=true, the phase is a mineral
 
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
 

 
ismin = .false.
!%----------------------------------------------------------- 
select case (this%itype)
case (minph)
 ismin = .true.
end select
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_name_ph &
   (this, &
    name)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of the phase 
!
!   $Arguments:
!
 
type(t_phase), intent(in)         :: this   ! Type phase variable 

character(len=*), intent(out)     :: name   ! Name of the phase 
 
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
call get_name_ (this%pp,name) 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_prop_ph &
   (this, &
    value, &
    nameprop, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return some pproperty according its name
!
!   $Arguments:
!
 
type(t_phase), intent(in)        :: this       ! Type phase variable 

real*8, intent(out)              :: value      ! Value of the property 

character(len=*), intent(in)     :: nameprop   ! Name of the property 

logical, intent(out)             :: iserror    ! iserror=true, then there was an error
 
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
call get_prop_ (this%pp, value, nameprop, iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_if_sp_is_present_ph &
   (this, &
    namesp, &
    isbe)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return if the species is present in the phase. 
!
!   $Arguments:
!
 
type(t_phase), intent(in)           :: this      ! Type phase variable 

character(len=*), intent(in)        :: namesp    ! Name of the species 

logical, intent(out)                :: isbe      ! isbe=true, the species is present
 
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
call get_if_sp_is_present_ (this%pp, namesp, isbe)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_charge_balance_ph &
   (this, &
    chgbal, &
    c, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute charge balance in the phase. 
!
!   $Arguments:
!
 
type(t_phase), intent(in)                            :: this    ! Type phase variable. 

real*8, intent(out)                                  :: chgbal  ! Charge balance. 

real*8, intent(in), dimension(this%pp%numsp)         :: c       ! Concentration vector. 

logical, intent(out)                                 :: iserror ! iserror=true, there was an error. 
 
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
!% Call the corresponding service in the parent phase 
!%------------------------------------------------------------
call compute_charge_balance_(this%pp,chgbal,c,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dcharge_balance_ph &
   (this, &
    dchgbal, &
    dc, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_phase), intent(in)               :: this     ! Type phase variable. 

real*8, pointer, dimension(:)           :: dchgbal  ! Derivatives of the charge balance 

real*8, intent(in), dimension(:,:)      :: dc       ! Derivatives of the concentrations 

logical, intent(out)                    :: iserror  ! iserror=true, there was an error. 
 
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
!% Call the corresponding service in the parent phase 
!%------------------------------------------------------------
call compute_dcharge_balance_ (this%pp,dchgbal,dc,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_xml_ph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in xml file the attributes encapsulated 
! in the phase object (not implemented).
!
!   $Arguments:
!
 
type(t_phase), intent(in)  :: this  ! Type phase variable. 
 
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
subroutine write_ph &
   (this, &
    ioutput, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii file the attributes encapsulated in the 
! phase object. 
!
!   $Arguments:
!
 
type(t_phase), intent(in)  :: this     ! Type phase variable. 

integer, intent(in)        :: ioutput  ! Output unit 

logical, intent(out)       :: iserror  ! iserror=true, there was an error. 
  
 
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
 

character(len=100) :: &
 msg
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
select case (this%itype)
 
case (aqph)
 
 call write_ (this%paqph,ioutput,iserror)
 
case (minph)
 
 call write_ (this%pminph,ioutput,iserror)
 
case (gasph)
 
 call write_ (this%pgasph,ioutput,iserror)
 
end select
!%------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Phase:'
print *,'Name:',this%pp%name
print *,'Service: write_'
print *, msg
print *,'************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
subroutine begin_element_handler(name,attributes)
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
 
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
   phase, & 
   iserror, &
   typephase)
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
character(len=*), intent(in)                :: name             ! Name of the tag 

type(dictionary_t), intent(in)              :: attributes

type(t_phase), optional                     :: phase            ! Type phase variable. 

character(len=100), intent(out), optional   :: typephase        ! Type phase (e.g. aqueous, mineral, gas)

logical, intent(out), optional              :: iserror          ! iserror=true, there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical         :: &
 havephase, &
 haveiserror, &
 havetypephase
character(len=100)     :: &
 id
character(len=100),save:: &
 typephaseloc
integer        :: &
 status 
!-------------------------------------------------------------------------
!
!   $code
!


!%----------------------------------------------------------------
!% Check optional arguments
!%----------------------------------------------------------------
havephase=present(phase)
haveiserror=present(iserror)
havetypephase=present(typephase) 
 
select case (havephase)
case (.true.)
 if (haveiserror) iserror=.false. 
 if (havetypephase) typephase=typephaseloc
!%---------------------------------------------------------
case default
 
  select case (name)
 
  case ('phase')
    call get_value (attributes,"type", id, status)
    typephaseloc=id
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
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine verify_ph &
  (this, &
   name, &
   numsp, &
   species, &
   itype, &
   nerror, &
   errormsg)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_phase), intent(in)                          :: this         ! Type phase variable. 

integer, intent(out)                               :: nerror       ! Number of errors. 

integer, intent(in), optional                      :: numsp        ! Number of species 

integer, optional                                  :: itype        ! 

type(t_species), pointer, optional, dimension(:)   :: species 

character(len=*), optional                         :: name 

character(len=100), pointer, dimension(:)          :: errormsg
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                  :: &
 i, &
 ndim, &
 nerrorlocal
logical                  :: &
 havenumsp, &
 havespecies, &
 haveitype, &
 havename
character(len=100)              :: &
 errormsgloc(50) 
!-------------------------------------------------------------------------
!
!   $code
!

 
nerror=0
 
if (associated(errormsg)) deallocate (errormsg)
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------ 
havenumsp=present(numsp)
havespecies=present(species)
haveitype=present (itype)
havename=present(name)
 
 
!%------------------------------------------------------------
if (havenumsp) then
 
 if(numsp.ne.this%pp%numsp) then
  nerror=nerror+1
  errormsgloc(nerror)='error en numsp'
 end if
 
end if
!%------------------------------------------------------------
if (havename) then
 
 if(name.ne.this%pp%name) then
   nerror=nerror+1
   errormsgloc(nerror)='error en name'
 end if
 
end if
!%------------------------------------------------------------
if (havespecies) then
 
  if (associated(species)) then
 
    select case (associated(this%pp%pspecies))
 
    case (.true.)
      ndim=size(species)
      if(ndim.ne.this%pp%numsp) then
        nerror=nerror+1
        errormsgloc(nerror)='error en species'
      else
        do i=1,ndim
         call verify_ &
           (this%pp%pspecies(i)%ptr, &
            species(i), &
            nerror=nerrorlocal, &
            errormsg=errormsg)
            if (nerrorlocal.gt.0) then
                   nerror=nerror+1
             errormsgloc(nerror)='error en species'
                  end if
        end do
        end if
 
    case default
     nerror=nerror+1
     errormsgloc(nerror)='error en species'
    end select
 
 
 
  end if
 
end if
!%------------------------------------------------------------
if (haveitype) then
 if(this%itype.ne.itype) nerror=nerror+1
 if (itype.eq.0) then
  if (associated(this%paqph)) then
    nerror=nerror+1
    errormsgloc(nerror)='error en aqphdh'
   end if
  if (associated(this%pminph)) then
    nerror=nerror+1
    errormsgloc(nerror)='error en ideal'
   end if
  if (associated(this%pminph)) then
    nerror=nerror+1
    errormsgloc(nerror)='error en mineral phase'
   end if
 
 else
 if (.not.associated(this%pp)) then
    nerror=nerror+1
    errormsgloc(nerror)='error en puntero parent'
   end if
 select case (itype)
 case(aqph)
   if (.not.associated(this%paqph)) then
    nerror=nerror+1
    errormsgloc(nerror)='error en aqphdh'
   end if
 case(minph)
   if (.not.associated(this%pminph)) then
    nerror=nerror+1
    errormsgloc(nerror)='error en ideal'
   end if
 
 
 case(gasph)
   if (.not.associated(this%pminph)) then
    nerror=nerror+1
    errormsgloc(nerror)='error en binarysolution'
   end if
 end select
 
 endif
 
end if
!%------------------------------------------------------------
 
if (nerror.ne.0) then
 allocate (errormsg(nerror))
 errormsg=errormsgloc(1:nerror)
end if
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine verify1_ph &
  (this, &
   phase, &
   nerror, &
   errormsg)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_phase), intent(in)                     :: this          ! Type phase variable. 

type(t_phase), intent(in)                     :: phase

integer, intent(out)                          :: nerror

character(len=100), pointer, dimension(:)     :: errormsg
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
type(t_species), pointer              :: &
 pspecies(:)
integer                              :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
 
 
if (phase%pp%numsp.gt.0) then
 allocate (pspecies(phase%pp%numsp))
end if
 
do i=1,phase%pp%numsp
 pspecies (i)= phase%pp%pspecies(i)%ptr
end do
 
call verify_ &
  (this, &
   name=phase%pp%name, &
   numsp=phase%pp%numsp, &
   species=pspecies, &
   itype=phase%itype, &
   nerror=nerror, &
   errormsg=errormsg)
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_ph &
  (copied, &
   this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy a object in other species object.  
!
!   $Arguments:
!
 
type(t_phase), intent(in)     :: this    ! Type phase variable. 

type(t_phase), intent(out)    :: copied  ! Type phase variable. 
 
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
 
 

!%-------------------------------------------------------------
copied%itype=this%itype
!%-------------------------------------------------------------
!% Copy the parent phase 
!%-------------------------------------------------------------
copied%pp=this%pp
!%-------------------------------------------------------------
!% Copy the specialization 
!%-------------------------------------------------------------
select case (this%itype)
 
case (aqph)
 
 allocate(copied%paqph)
 call create_ (copied%paqph)
 copied%paqph = this%paqph
 call set_pparent_ (copied%paqph,copied%pp)
 
case (minph)
 
 allocate (copied%pminph)
 call create_ (copied%pminph)
 copied%pminph = this%pminph
 call set_pparent_ (copied%pminph,copied%pp)
 
case (gasph)
 
 allocate (copied%pgasph)
 call create_ (copied%pgasph)
 copied%pgasph = this%pgasph
 call set_pparent_ (copied%pgasph,copied%pp) 
 
 
end select
!%------------------------------------------------------------
 
return
end subroutine




!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_density_ph &
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
!   $Description: Compute density 
!
!   $Arguments:
!
type(t_phase), intent(in)                        :: this         !Type phase 

real*8,intent(in)                                :: pliq         !liquid pressure

real*8,intent(in)                                :: temp         !temperature

real*8,dimension(this%pp%numsp), intent(in)      :: c            !concentration array

real*8,intent(out)                               :: density      !Density     

logical,intent(in)                               :: ismolality   !iserror=true, there was an error 

logical,intent(out)                              :: IsError      !iserror=true, there was an error 

real*8,dimension(:,:),optional                   :: dc           !concentration derivatives (needed to calculate density derivatives)

real*8,pointer,dimension(:),optional             :: ddensitydc   !Density derivative

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


if (this%itype.eq.aqph) then

    if  (present(ddensitydc) ) then
        if (.not.present(dc)) then
            msg="dc are needed for calcutate density derivatives"
            isError=.true.
            goto 10
        endif
       
        call compute_density_(this%paqph,pliq,temp,c,density,ismolality,iserror,dc,ddensitydc)
        if (isError) goto 10
    else
        call compute_density_(this%paqph,pliq,temp,c,density,ismolality,iserror)
        if (isError) goto 10
    endif

else

    msg='Density function only supported for liquid phase'
    isError=.true.
    goto 10
    
endif

!%----------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Phase:'
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
end module m_phase