module m_reactionratelaw
!-------------------------------------------------------------------------
!
!   $Description: Reaction rate law class represents the mathematical expression of one kinetic law of one chemical reaction (see figure \ref{fig:chemsys}b). 
! In the same way that the previous classes, reaction rate law class is also associated to species class and in addition it contains all parameters that define a law (e.g. rate constants, half saturation constants, inhibition constants, activation energy). 
! Thus, an reaction rate law object can have none or a lot of species objects.
! These species participate in the kinetic law as catalytic and inhibiting and not necessarily these species should be participating in
! the stoichiometry of the reaction. 
! Internally reaction rate law class is organized by terms being these monod, inhibitors, or affinity terms.
! This last is called also like the far-from-equilibrium term, decreasing the reaction rate in a non-linear way, as the solution approaches to equilibrium. 
!
!   $Use: m_species
! m_parentreactionratelaw
! m_reactionratelaw_lasaga
! m_reactionratelaw_monod
! flib_xpath
! flib_sax
! m_general_tools_cheproo
! m_constants
!
!   $Author: Sergio Andrés Bea Jofré 
!
!   $License: UPC-CSIC
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_species
use m_parentreactionratelaw
use m_reactionratelaw_lasaga
use m_reactionratelaw_monod
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%------------------------------------------------------------------
!%------------------------------------------------------------------
private                     ::
!%------------------------------------------------------------------
!%------------------------------------------------------------------
public                      :: &
create_ &                ! Create the reaction rate law object 
,read_xml_ &             ! Read the reaction rate law object from xml file. 
,destroy_ &              ! Destroy the reaction rate law object. 
,set_ &                  ! Set the different attributes in the reaction rate law object. 
,set_pspecies_ &         ! Set the 
,compute_rk_ &
,compute_drk_ &
,update_ &
,get_namesp_ &           ! Return the name of species in the reaction rate law
,get_if_activities_ &    ! 
,write_ &                ! Write in ascii format the reaction rate law object. 
,rewrite_ &              ! Rewrite the stoichiometric matrix of catalytic species
,assignment(=)           ! Copy the reaction rate law object in other reaction rate law object. 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
private :: &
read_xml_loc_ &
,begin_element_handler
!%-----------------------------------------------------------------
!% Constant parameters 
!%-----------------------------------------------------------------
integer, parameter                       :: &
lasaga=1, &
monod=2
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Type pointer to reaction rate law object 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
type, public::t_prrlaw

type(t_reactionratelaw), pointer::ptr

end type
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
type, public:: t_reactionratelaw
 
private                                  ::
 
type (t_parentreactionratelaw), pointer  :: pp              ! Parent reaction rate law pointer 
 
type (t_reactionratelaw_lasaga), pointer :: prrlaw_lasaga   ! Specialization type lasaga
 
type (t_reactionratelaw_monod), pointer  :: prrlaw_monod    ! Specialization type monod
 
integer                                  :: itype           ! Index of specialization 
 
end type t_reactionratelaw
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface create_
 
module procedure create_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface read_xml_
 
module procedure read_xml_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface destroy_
 
module procedure destroy_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_
 
module procedure set_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_pspecies_
 
module procedure set_pspecies_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_rk_
 
module procedure compute_rk_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_drk_
 
module procedure compute_drk_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface update_
 
module procedure update_temp_param_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_namesp_
 
module procedure get_namesp_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_if_activities_
 
module procedure get_if_activities_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface write_
 
module procedure write_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface rewrite_
 
module procedure rewrite_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface assignment(=)
 
module procedure copy_rrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface read_xml_loc_
 
module procedure read_xml_loc_rrlawaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface begin_element_handler
 
module procedure begin_element_handler
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
contains
!%-----------------------------------------------------------------
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_rrlaw &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create reaction rate law object
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(inout) :: this 
 
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
this%prrlaw_lasaga => null ()
this%prrlaw_monod => null ()
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
subroutine destroy_rrlaw &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy the reaction rate law object
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(inout)   :: this 
 
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
end if
!%------------------------------------------------------------
this%pp => null ()
!%------------------------------------------------------------
select case (this%itype)
case (lasaga)
 call destroy_ (this%prrlaw_lasaga)
 deallocate (this%prrlaw_lasaga)
 this%prrlaw_lasaga => null ()
case (monod)
 call destroy_ (this%prrlaw_monod)
 deallocate (this%prrlaw_monod)
 this%prrlaw_monod => null ()
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
subroutine read_xml_rrlaw &
  (this, &
   namefile, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read reaction rate law from xml file
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(inout)    :: this

character(len=*), intent(in), optional    :: namefile

logical, intent(out)                      :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                               :: &
 iostat
type(xml_t):: &
 fxml
type(dictionary_t)                    :: &
 attributes
character(len=100)                    :: &
 name, &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%--------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------Open xml file
call open_xmlfile(namefile,fxml, iostat)
if (iostat /= 0) goto 10
!%--------------------------------------------------------
call xml_parse(fxml,begin_element_handler = begin_element_handler)
!%--------------------------------------------------------
call read_xml_loc_ (name,attributes,this,msg,iserror)
if (iserror) goto 10
!%--------------------------------------------------------
! End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%-------------------------------------------------------
select case (this%itype)
case (lasaga)
 allocate (this%prrlaw_lasaga)
 call create_ (this%prrlaw_lasaga)
 call read_xml_(this%prrlaw_lasaga,this%pp,namefile,iserror)
case (monod)
 allocate (this%prrlaw_monod)
 call create_ (this%prrlaw_monod)
 call read_xml_(this%prrlaw_monod,this%pp,namefile,iserror)
end select
!%--------------------------------------------------------
return
10 continue 
print *,'**********************'
print *,'Reaction Rate Law:'
print *,'Name:',this%pp%name
print *,'Service: read_xml_'
print *, msg
print *,'**********************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_namesp_rrlaw &
  (this, &
   namesp, &
   nsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of species in the reaction rate law
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(in)    :: this

character(len=*), dimension(:), pointer:: namesp

integer, intent(out)                   :: nsp 
 
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
call get_namesp_ (this%pp,namesp,nsp)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_if_activities_rrlaw &
  (this, &
   isbeact)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(in)    :: this

logical, intent(out)                   :: isbeact 
 
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
case (lasaga)
 isbeact=.true.  
case (monod)
 isbeact=.false.  
end select
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_rrlaw &
  (this, &
   itype, &
   name, &
   nspterm, &
   attrterm, &
   attrspterm, &
   namespterm, &
   typeterm, &
   nterm, &
   nattrterm, &
   mxsp, &
   ea, &
   attrrrlaw1, &     
   attrrrlaw2, &
   attrrrlaw3, &
   attrrrlaw4, &     
   isareadep, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the reaction rate law object
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(inout)  :: this                             ! Reaction rate law type

integer, intent(in)                     :: itype                            ! Specialization index 

character(len=*), intent(in)            :: name                             ! Name of the reaction rate law object 

integer, intent(in)                     :: nterm                            ! Number of terms 

integer, intent(in)                     :: nattrterm                        ! Number of attributes per term 

integer, intent(in)                     :: mxsp                             ! Maximum number of species 

real*8, intent(in)                      :: attrrrlaw1                       ! Maximum rate (monod)
                                                                            ! omegthr (lasaga)

real*8, intent(in)                      :: attrrrlaw2                       ! Minimum rate (monod)

real*8, intent(in)                      :: attrrrlaw3                       ! theta (monod)

real*8, intent(in)                      :: attrrrlaw4                       ! eta (monod)

character(len=*), intent(in)            :: namespterm(mxsp,nterm)

character(len=*), intent(in)            :: typeterm(nterm)

real*8, intent(in)                      :: ea                               ! Activation energy 

real*8, intent(in)                      :: attrspterm(mxsp,nterm)

real*8, intent(in)                      :: attrterm(nattrterm,nterm)        

integer, intent(in)                     :: nspterm(nterm)                   ! Number of species per term 

logical, intent(in)                     :: isareadep                        ! isareadep=true, then the reaction rate law depends of the reactive surface

logical, intent(out)                    :: iserror                          ! iserror=true, then there was an error
 
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
!%----------------------------------------------------------
msg=''
iserror=.false.
!%----------------------------------------------------------
!% Set the parent reaction rate law 
!%----------------------------------------------------------
 call set_ &
  (this%pp, &
   name, &
   nspterm, &
   attrterm, &
   attrspterm, &
   namespterm, &
   typeterm, &
   nterm, &
   nattrterm, &
   mxsp, &
   ea, &
   isareadep, &
   iserror)
!%------------------------------------------------------------
!% Store the specialization 
!%------------------------------------------------------------
this%itype=itype
!%------------------------------------------------------------
!% Create and set the specialization 
!%------------------------------------------------------------
select case (this%itype)
case (lasaga)
 allocate (this%prrlaw_lasaga)
 call create_(this%prrlaw_lasaga)
 call set_ (this%prrlaw_lasaga,this%pp,attrrrlaw1,iserror)
case (monod)
 allocate (this%prrlaw_monod)
 call create_(this%prrlaw_monod)
 call set_ (this%prrlaw_monod,this%pp,attrrrlaw1,attrrrlaw2,attrrrlaw3,attrrrlaw4,iserror)
case default
 msg='Error, not defined reaction rate law specialization'
 goto 10
end select
!%------------------------------------------------------------
return
10 continue 
print *,'**********************'
print *,'Reaction Rate Law:'
print *,'Name:',this%pp%name
print *,'Service: set_'
print *, msg
print *,'**********************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_pspecies_rrlaw &
  (this, &
   species)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(inout)      :: this

type(t_species), intent(in)                 :: species 
 
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
call set_pspecies_ (this%pp,species)
!%------------------------------------------------------------
return
 
end subroutine
!%********************
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_rrlaw &
  (this, &
   temp, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:Update termperature parameters in the reaction rate law objects
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(inout) :: this      ! Type reaction rate law 

real*8, intent(in)                     :: temp      ! Temperature in celcius 

logical, intent(out)                   :: iserror   ! iserror=true, then there was an error 
 
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
if (this%pp%tempref==temp) return 
this%pp%tempref=temp 
!%------------------------------------------------------------
call update_ (this%pp,this%pp%tempref,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_rrlaw &
  (this, &
   ioutput, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: write in txt format attributtes in the reaction rate
!%    law
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(in) :: this

integer, intent(in)                 :: ioutput 

logical, intent(out)                :: iserror 
 
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
case (lasaga)
 call write_ (this%prrlaw_lasaga,ioutput,iserror)
case (monod)
 call write_ (this%prrlaw_monod,ioutput,iserror)
end select
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_rk_rrlaw &
  (this, &
   rk, &
   omega, &
   c, &
   g, &
   nsp, &
   cold, &
   area, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute reaction rate. 
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(in)            :: this     ! Reaction rate law object 

real*8, intent(out)                            :: rk       ! reaction rate [mol m-3 s-1]

integer, intent(in)                            :: nsp      ! Number of species

real*8, intent(in)                             :: omega    ! IAP/Ke

real*8, intent(in)                             :: area

real*8, intent(in)                             :: cold     ! Previous concentrations

logical, intent(out)                           :: iserror  ! iserror=true, then there was an error

real*8, intent(in), optional, dimension(nsp)   :: c        ! Concentration vector 

real*8, intent(in), optional, dimension(nsp)   :: g        ! Activity coefficient vector 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                 :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%--------------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------------
!% Select the specialization 
!%--------------------------------------------------------------------
select case (this%itype)
case (lasaga)
 call compute_rk_ (this%prrlaw_lasaga,rk,omega,c,g,nsp,cold,area,iserror)
case (monod)
 call compute_rk_ (this%prrlaw_monod,rk,omega,c,nsp,iserror)
end select
!%-------------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Reaction Rate Law:'
print *,'Name:',this%pp%name
print *,'Service: compute_rk_'
print *, msg
print *,'**********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_drk_rrlaw &
  (this, &
   drk, &
   omega, &
   domega, &
   c, &
   dc, &
   g, &
   dg, &
   nsp, &
   ndimder, &
   area, &
   isbe, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute reaction ratederivatives
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(in):: this

integer, intent(in)                :: nsp              ! number of species 

integer, intent(in)                :: ndimder          ! number of derivatives 

real*8, pointer                    :: drk(:)           

real*8, intent(in)                 :: domega(ndimder)  ! d (IAP/Ke) / dc1 

real*8, intent(in)                 :: omega            ! IAP/Ke

real*8, intent(in)                 :: area             ! Area  

logical, intent(in)                :: isbe             ! isbe=true, mineral is present 

logical, intent(out)               :: iserror       

real*8, intent(in)                 :: c(nsp)           ! molality vector

real*8, intent(in)                 :: g(nsp)           ! activity coefficients vector 

real*8, intent(in)                 :: dc(nsp,ndimder)  ! derivatives of molality vector 

real*8, intent(in)                 :: dg(nsp,ndimder)  ! derivatives of activity coefficients 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                 :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%--------------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------------
!% Select the specialization 
!%--------------------------------------------------------------------
select case (this%itype)
case (lasaga)
 call compute_drk_(this%prrlaw_lasaga,drk,omega,domega,c,dc,g,dg,nsp,ndimder,area,isbe,iserror)
case (monod)
 call compute_drk_ (this%prrlaw_monod,drk,omega,domega,c,dc,nsp,ndimder,iserror)
end select
!%-------------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Reaction Rate Law:'
print *,'Name:',this%pp%name
print *,'Service: compute_drk_'
print *, msg
print *,'**********************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_rrlaw &
  (copied, &
   this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:Copy reaction rate law object
!
!   $Arguments:
!
 
 
type(t_reactionratelaw), intent(in) :: this

type(t_reactionratelaw), intent(out):: copied 
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
copied%itype=this%itype
!%---------------------------------------------------------
select case (copied%itype)
case (lasaga)
 if (associated(copied%prrlaw_lasaga)) then
  call destroy_ (copied%prrlaw_lasaga)
  deallocate (copied%prrlaw_lasaga)
  copied%prrlaw_lasaga => null ()
 end if
case (monod)
 if (associated(copied%prrlaw_monod)) then
  call destroy_ (copied%prrlaw_monod)
  deallocate (copied%prrlaw_monod)
  copied%prrlaw_monod => null ()
 end if
end select
!%---------------------------------------------------------
if (associated(copied%pp)) then
 call destroy_ (copied%pp)
else
 allocate (copied%pp)
 call create_ (copied%pp)
end if
!%---------------------------------------------------------
copied%pp=this%pp
!%---------------------------------------------------------
select case (this%itype)
case (lasaga)
 allocate (copied%prrlaw_lasaga)
 call create_ (copied%prrlaw_lasaga)
 call set_ (copied%prrlaw_lasaga,copied%pp)
 copied%prrlaw_lasaga=this%prrlaw_lasaga
case (monod)
 allocate (copied%prrlaw_monod)
 call create_ (copied%prrlaw_monod)
 call set_ (copied%prrlaw_monod,copied%pp)
 copied%prrlaw_monod=this%prrlaw_monod
end select
!%---------------------------------------------------------
return
end subroutine
!%***********************************************************
!%***************Public subroutine***************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
subroutine rewrite_rrlaw &
  (this, &
   namesp, &
   nsp, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_reactionratelaw), intent(inout)        :: this

integer, intent(in)                           :: nsp

character(len=*), intent(in), dimension(nsp)  :: namesp

logical, intent(out)                          :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------


!%---------------------------------------------------------
select case (this%itype)
case (lasaga)
  call rewrite_ (this%prrlaw_lasaga,namesp,nsp,iserror) 
case (monod)
  call rewrite_ (this%prrlaw_monod,namesp,nsp,iserror) 
end select
!%--------------------------------------------------------------- 
return
 
end subroutine
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine begin_element_handler (name, attributes)
 
 
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
 
 
 call read_xml_loc_ (name,attributes)
 
return
end subroutine
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc_rrlawaw &
 (name, &
  attributes, &
  this, &
  msg, &
  iserror)
 
 
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
type (t_reactionratelaw ),        optional:: &
 this
logical, optional              :: iserror
character(len=*), optional     :: msg
 
character(len=100)             :: &
 id
character(len=100), save       :: &
 typerrlaw
logical                        :: &
 havethis
integer                        :: &
 status
!%--------------------------------------------------------------
havethis=present(this)
!%--------------------------------------------------------------
select case (havethis)
 
case (.true.)
 
 msg=''
 
 select case (typerrlaw)
 case ('LASAGA','lasaga','Lasaga')
  this%itype=lasaga
 case ('MONOD','Monod','monod')
  this%itype=monod
 case default
  msg='Error, not recognized reaction rate law type'
  return
 end select
 
case default
 
 select case (name)
 
 case ('reactionratelaw')
   id=''
   call get_value (attributes,"type", id, status)
   typerrlaw=id
 end select
 
end select
!%-------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_reactionratelaw
