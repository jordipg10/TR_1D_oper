module m_reaction
!-------------------------------------------------------------------------
!
!   $Description: The reaction class represents a chemical reaction in a geochemical system. 
! It establishes the relationship between species in the same phase (homogeneous
! reaction), or species in different phases (heterogeneous reaction). The reaction
! class is associated to species and reaction rate law classes. In addition, it encapsulates other data like the stoichiometric coefficients,
! the equilibrium constant and the power function coefficients for the temperature
! dependence of the equilibrium constant. The relationship between species
! (Sps) of jth reaction can be formalized as
!
! Spsj <==> SUM(Ns) ( stqji Spi)   Kj(T)
!
! where Ns is the number of reacting species, stqij the stoichiometric coefficient of
! the ith reacting species in the jth reaction, and Kj(T) the equilibrium constant
! of the jth reaction, which is a function of temperature (T). At equilibrium,
! the relation (3) can be expressed by the mass action equation
!
!                
!                
! xj = Kj^(-1) gj^(-1) MULT(Ns)(ci gi)^(stqji) MULT(Nk) (sk)^(stqadsjk)
!
!where g are activity coefficients, and c and x are concentrations of the participating
!species. It is always possible to write the reactions in terms of a
!subset of c called primary species (dimension [Ns-Nr], where Nr is the number
!of reactions), so that the secondary species can be obtained explicitly
!from equation (4) and do not participate in any other reaction. The main
!method implemented in reaction class solve the equation (4) (see Figure 1c,
!’compute xj ’) and calculates its derivatives with respect to the participating
!species. Note that the meaning of xj in equation (4) depends of the type
!of reaction. Thus, it could be the concentration of an aqueous or a surface
!complex, a ratio between ionic activity product and equilibrium constant, or
!a partial pressure of a gas species. If a reaction object is defined as kinetic,
!the class can additionally calculate the reaction rate (and its derivatives with
!respect to all species). Internally the reaction rate calculation is carried out
!by the associated reaction rate law object.
!
!   $Use: use m_species
! use m_reactionratelaw
! use flib_xpath
! use flib_sax
! use m_general_tools_cheproo
! use m_constants
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
use m_reactionratelaw
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%------------------------------------------------------------------
!%------------------------------------------------------------------
private                   ::
!%------------------------------------------------------------------
!% Public services 
!%------------------------------------------------------------------
public                    :: &
create_ &                 ! Create and initialice reaction object. 
,destroy_ &               ! Destroy reaction object 
,compute_x_ &             ! Solve the explicitly exprexssion of the mass action equation.
,compute_dx_ &            ! Compute derivatives of the xj  
,compute_rk_ &            ! Compute reaction rate (only if the reaction was defined in kinetic) 
,compute_drk_ &           ! Compute derivatives of the kinetic rate 
,get_namesp_ &            ! Return the name of the species involved in the reaction. 
,get_stq_ &               ! Return the stoichiometric coefficients of the species involved in the reaction. 
,get_lnk_ &               ! Return the logarithm of the equilibrium constant. 
,get_logk_ &              ! Return the logarithm of the equilibrium constant. 
,get_numsp_ &             ! Return the number of species involved in the reaction. 
,get_if_sp_is_present_ &  ! Return if a species participates in the reaction. 
,get_coeff_stq_ &         ! Return the stoichiometric coefficient corresponding to a species. 
,get_name_ &              ! Return the name of the reaction (in general is associated with the product species). 
,get_name_sec_ &          ! Return the name of the secondary species associated to reaction. 
,get_if_kinetic_ &        ! Return if the reaction is kinetic. 
,get_polinom_logk_ &      ! Return the polinomial coefficients corresponding to the dependence of the equilibrium constant with the temperature. 
,get_namesp_rrlaw_ &      ! Return the species involved in reaction rate law. 
,get_prrlaw_ &            ! Return the pointer to reaction rate law object (only if the reaction is kinetic). 
,read_xml_ &              ! Read and set the reaction object from xml file. 
,rewrite_ &               ! Rewrite the stoichiometric matrix according to species list 
,rewrite_rrlaw_ &         ! Rewrite the stoichiometric matrix of catalytic species in the reaction rate law object according to species list (only if the reaction is kinetic)
,set_ &                   ! Set general attributes in the reaction object 
,set_pspecies_ &          ! Set the pointer to species object
,set_prrlaw_ &            ! Set the pointer to reaction rate law object
,set_stqads_ &            ! Set the stoichimetric matrix corresponding to sorption primary species according to sorption model. 
,set_coeff_logk_ &        ! Set the polinomial coefficients of the Log K = F(T). 
,update_ &                ! Update all attributtes in the reaction object that depend of the temperature 
,verify_ &                ! Verify different attributes encapsulated in the reaction object. 
,write_ &                 ! Write in ioutput unit all attributtes encapsulated in the reaction object
,assignment(=)            ! Copy the reaction object in other reaction object. 
!%-----------------------------------------------------------------
!% Private services
!%-----------------------------------------------------------------
private                   :: &
LEAST_SQUARES_FIT &
,begin_element_reaction &
,read_xml_loc_ &
,compute_dx_dc1_numeric_reaction
!%-----------------------------------------------------------------
!% Parameters 
!%-----------------------------------------------------------------
real*8, parameter     :: &
tempref=25.0d0               ! Reference temperature 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Type pointer to reaction object 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
type, public::t_preaction

type(t_reaction), pointer::ptr

end type
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Type definition 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
type, public::   t_reaction
 
private                                  ::
 
character(len=20)                        :: name         ! Name of reaction (in general is the nam product chemical species)
  
type (t_pspecies), pointer, dimension(:) :: pspecies     ! Pointer to species objects
 
type(t_reactionratelaw), pointer         :: prrlaw       ! Pointer to reaction rate law
 
real*8, pointer, dimension(:)            :: coefflogk    ! Polinomial coefficients for log(k)=f(T)
 
real*8                                   :: lnk          ! Log(e) of the equilibrium constant to different temperatures

real*8, pointer, dimension(:)            :: logktemp     ! Equilibrium constant to different temperatures [ntemp]

real*8, pointer, dimension(:)            :: templogk     ! Temperatures corresponding to logktemp vector [numtemp]
 
real*8, pointer, dimension(:)            :: stq          ! Stoichiometric coefficients for the species in the reaction [numsp]

real*8, pointer, dimension(:)            :: stqaqpri     ! Stoichiometric vector corresponding to primary species

real*8, pointer, dimension(:)            :: stqadspri    ! Stoichiometric vector corresponding to adicional adsorption primary species
 
integer                                  :: numsp        ! Number of species of the reaction

integer                                  :: numcoefflogk ! Number of coefficient for log (K) = f(T)

integer                                  :: numtemp      ! Number of temperatures corresponding to log k vector 

integer                                  :: ithsecsp     ! Local index of the secondary species associated to the reaction 

real*8                                   :: tempref      ! Last temperature that was updated the reaction 

real*8                                   :: pressref     ! Last pressure that was updated the reaction 
 
logical                                  :: iskin        ! If .true. the reaction object is defined kinetic

logical, pointer, dimension(:)           :: islocksps    ! If .true. the species object was allocated for reaction object [numsp]

logical                                  :: islockrrlaw  ! If .true. the reaction rate law object was allocated for reaction object 
 
end type t_reaction
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface create_
 
module procedure create_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_x_
 
module procedure compute_x_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_dx_
 
module procedure compute_dx_reaction
module procedure compute_dx_dsk_reaction
module procedure compute_dx_numeric_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_rk_
 
module procedure compute_rk_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_drk_
 
module procedure compute_drk_dx_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface destroy_
 
module procedure destroy_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_namesp_
 
module procedure get_namesp_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_stq_
 
module procedure get_stq_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_lnk_
 
module procedure get_lnk_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_logk_
 
module procedure get_logk_temp_reaction

end interface 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_numsp_
 
module procedure get_numsp_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_prrlaw_
 
module procedure get_prrlaw_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_if_sp_is_present_
 
module procedure get_if_sp_is_present_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_coeff_stq_
 
module procedure get_coeff_stq_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_name_
 
module procedure get_name_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_name_sec_
 
module procedure get_name_sec_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_if_kinetic_
 
module procedure get_if_kinetic_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_polinom_logk_
 
module procedure get_polinom_logk_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_namesp_rrlaw_
 
module procedure get_namesp_rrlaw_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface read_xml_
 
module procedure read_xml_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_stqads_
 
module procedure set_stqads_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface rewrite_
 
module procedure rewrite_stq_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface rewrite_rrlaw_
 
module procedure rewrite_rrlaw_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_
 
module procedure set_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_logk_
 
module procedure set_logk_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_prrlaw_
 
module procedure set_prrlaw_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_coeff_logk_
 
module procedure set_coeff_logk_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_pspecies_
 
module procedure set_pspecies_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface update_
 
module procedure update_temp_param_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface verify_
 
module procedure verify_reaction
module procedure verify2_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface write_
 
module procedure write_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface assignment(=)
 
module procedure copy_reaction
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------Private services--------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface LEAST_SQUARES_FIT
 
module procedure LEAST_SQUARES_FIT
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface begin_element_handler
 
module procedure begin_element_handler
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface read_xml_loc_
 
module procedure read_xml_loc_reaction 
 
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
subroutine create_reaction &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create reaction object
!
!   $Arguments:
!
 
type(t_reaction), intent(inout) :: this ! Type reaction variable. 
 
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
this%coefflogk => null ()
this%pspecies => null ()
this%stq => null ()
this%stqaqpri => null ()
this%stqadspri => null ()
this%prrlaw => null ()
this%islocksps => null ()
this%logktemp => null () 
this%templogk => null ()
!%------------------------------------------------------------
this%lnk=0.0d0 
this%numtemp=0 
this%numcoefflogk=0
this%name=''
this%ithsecsp=0
this%numsp=0
this%iskin=.false.
this%tempref=-1.0d3
this%pressref=0.0d0 
this%islockrrlaw=.false.
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_reaction &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy reaction object.
!
!   $Arguments:
!
 
type(t_reaction), intent(inout) :: this   ! Type reaction variable.
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer          :: i 
!-------------------------------------------------------------------------
!
!   $code
!
this%numcoefflogk = 0
this%numtemp=0 
this%name=' '
this%ithsecsp=0
this%lnk=0.0d0 
this%tempref=-1.0d3
this%pressref=0.0d0
!%--------------------------------------------------------------
!% Destroy reaction rate law object (only if was locked)
!%--------------------------------------------------------------
if (this%islockrrlaw) then
 call destroy_ (this%prrlaw)
 deallocate (this%prrlaw)
end if  
this%prrlaw => null ()
this%islockrrlaw=.false.
this%iskin=.false.
!%--------------------------------------------------------------
!% Destroy the species objects (only if was locked)
!%--------------------------------------------------------------
if (this%numsp>0) then
 do i=1,this%numsp
  if(this%islocksps(i)) then
    call destroy_ (this%pspecies(i)%ptr)
    deallocate (this%pspecies(i)%ptr) 
  end if 
  this%pspecies(i)%ptr => null ()
 end do
 deallocate(this%pspecies) 
 this%pspecies => null ()
 this%numsp=0
end if
!%------------------------------------------------------------
call check_pointer_ (this%stq,1,.false.)
call check_pointer_ (this%templogk,1,.false.)
call check_pointer_ (this%logktemp,1,.false.)
call check_pointer_ (this%stqaqpri,1,.false.)
call check_pointer_ (this%stqadspri,1,.false.)
call check_pointer_ (this%coefflogk,1,.false.)
call check_pointer_ (this%islocksps,1,.false.)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_reaction &
  (this, &
   namefile, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read the reaction object from xml file 
!
!   $Arguments:
!
 
type(t_reaction), intent(inout)   :: this     ! Type reaction variable. 

character(len=*), intent(in)      :: namefile ! Name (including path) of xml file.

logical, intent(out)              :: iserror  ! iserror=true, then there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer::iostat
type(xml_t)::fxml
character(len=100)             :: &
 name, &
 msg
type(dictionary_t)             :: &
 attributes 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
iserror=.false.
msg='' 
!%------------------------------------------------------------
!% Open the xml file 
!%------------------------------------------------------------
call open_xmlfile(namefile, fxml, iostat)
if (iostat/=0) then
 msg='Error when open file: '
 call add_ (msg,namefile) 
 goto 10
end if 
!%------------------------------------------------------------
call xml_parse(fxml, &
         begin_element_handler = begin_element_handler, &
         verbose = .false.)
!%--------------------------------------------------------
! End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%------------------------------------------------------------
call read_xml_loc_ (name,attributes,this)
!%------------------------------------------------------------
return
10 continue 
print *,'**********************'
print *,'Reaction:'
print *,'Name:',this%name
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
subroutine get_if_kinetic_reaction &
  (this, &
   iskin)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return if the reaction is kinetic
!
!   $Arguments:
!
 
type(t_reaction), intent(in) :: this    ! Type reaction variable.

logical, intent(out)         :: iskin   ! iskin=true, the reaction is kinetis. 
 
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
iskin=this%iskin
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_prrlaw_reaction &
  (this, &
   rrlaw)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the pointer to reaction rate law object. If the reaction
! rate object was previously allocated for the reaction object then destroy 
! and deallocate the reaction rate object
!
!   $Arguments:
!
 
type(t_reaction), intent(inout)                  :: this   ! Type reaction object 

type(t_reactionratelaw), intent(in), target      :: rrlaw  ! Type reaction rate law object  
 
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
this%iskin=.true.
!%------------------------------------------------------------
if (this%islockrrlaw) then
 call destroy_ (this%prrlaw)
 deallocate (this%prrlaw)
 this%islockrrlaw=.false.  
end if  
!%------------------------------------------------------------
this%prrlaw => rrlaw
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_coeff_logk_reaction &
  (this, &
   coefflogk, &
   ncoeff, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the polinomial coefficients of the Log K = F(T)
!
!   $Arguments:
!
 
type(t_reaction), intent(inout)              :: this        ! Type reaction variable.

integer, intent(in)                          :: ncoeff      ! Number of polinomila coefficients.

real*8, intent(in), dimension(ncoeff)        :: coefflogk   ! Polinomial coefficients [ncoeff]

logical, intent(out)                         :: iserror     ! iserror=true, then there was an error
 
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
!%------------------------------------------------------------
msg='' 
iserror=.false. 
!%------------------------------------------------------------
this%numcoefflogk=ncoeff
!%------------------------------------------------------------
call check_pointer_ (this%coefflogk,this%numcoefflogk,.true.)
!%------------------------------------------------------------
this%coefflogk=coefflogk
!%------------------------------------------------------------
!% Compute the equilibrium constant to reference temperature
!% (generally 25oC)
!%------------------------------------------------------------
call update_ (this,tempref,iserror)
if (iserror) then
 msg='Error when calling update_'
 goto 10
end if 
!%------------------------------------------------------------
return
10 continue 
print *,'************************'
print *,'Reaction:'
print *,'Name:',this%name
print *,'Service: set_coeff_logk_'
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
subroutine set_logk_reaction &
  (this, &
   logk, &
   temp, &
   ntemp, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the list of equilibrium constants. 
!
!   $Arguments:
!
 
type(t_reaction), intent(inout)        :: this      ! Type reaction variable. 

integer, intent(in)                    :: ntemp     ! Number of temperatures (the same dimension of log10 equilibrium constant vector). 
 
real*8, intent(in), dimension(ntemp)   :: logk      ! log10 of equilibrium constants [ntemp] 

real*8, intent(in), dimension(ntemp)   :: temp      ! Temperatures [ntemp]

logical, intent(out)                   :: iserror   ! iserror=true, there was an error. 
 
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
 tk
real*8, pointer        :: &
 atbsfn(:,:) => null (), &
 anrmt(:,:) => null ()
integer       :: &
 i, &
 j, &
 k, &
 ifail, &
 iprint, &
 ndnrmt, &
 ioutput
integer, pointer           :: &
 indx(:) => null()
integer               :: &
 id 
character(len=100)    :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
!%----------------------------------------------------------------------
!% Definition of temp. dependence-basis functions
!%----------------------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------------------
!% Storage the information in the reaction object 
!% (storage in base 10 logaritm)
!%----------------------------------------------------------------------
if (ntemp>0) then 
 this%numtemp=ntemp
 call check_pointer_ (this%logktemp,this%numtemp,.true.)
 call check_pointer_ (this%templogk,this%numtemp,.true.)
 this%logktemp=logk
 this%templogk=temp 
end if 
!%----------------------------------------------------------------------
!% If ntemp>1, then compute polinomial coefficients for temperature 
!% dependence of the equilibrium constant 
!%----------------------------------------------------------------------
if (ntemp>1) then
!%----------------------------------------------------------------------
!% Allocate local pointers 
!%----------------------------------------------------------------------
  call check_pointer_ (atbsfn,this%numcoefflogk,this%numtemp,.true.)
  call check_pointer_ (anrmt,this%numcoefflogk,this%numcoefflogk,.true.)
  call check_pointer_ (indx,this%numcoefflogk,.true.)
  tk=273.15
  do i=1,ntemp
        atbsfn(1,i) = dlog(temp(i) + tk)
        atbsfn(2,i) = 1.0d0
        atbsfn(3,i) = temp(i) + tk
        atbsfn(4,i) = 1.0d0/(temp(i) + tk)
        atbsfn(5,i) = 1.0d0/((temp(i) + tk)*(temp(i) + tk))
  end do
  ndnrmt = 0
!%----------------------------------------------------------------------
!% Compute polynomial coefficient using least square method. 
!%----------------------------------------------------------------------  
  call least_squares_fit &
 (this%numcoefflogk, &
  this%numcoefflogk, &
  ntemp, &
  ndnrmt, &
  ifail, &
  iprint, &
  ioutput, &
  logk, &
  atbsfn, &
  anrmt, &
  this%coefflogk, &
  indx, &
  msg, &
  iserror)
!%----------------------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------------------
  call check_pointer_ (atbsfn,1,1,.false.)
  call check_pointer_ (anrmt,1,1,.false.)
  call check_pointer_ (indx,1,.false.)
  if (iserror) goto 10 
!%----------------------------------------------------------------------
!% Update the equilibrium constant to temperature reference
!% (generally 25 oC)
!%----------------------------------------------------------------------
  call update_ (this,tempref,iserror)
  if (iserror) then
     msg='Error when calling update_'
     goto 10
  end if
else 
 this%tempref=this%templogk(1)
 this%lnk = logk(1)*dlog(10.0d0)
end if
!%----------------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Reaction:'
print *,'Name:',this%name
print *,'Service: set_logk_'
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
subroutine set_reaction &
  (this, &
   name, &
   iskin, &
   stq, &
   namesp, &
   logk, &
   temp, &
   ntemp, &
   ncoefflogk, &
   nsp, &
   iserror, &
   species, &
   rrlaw)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set general attributes in the reaction object. 
!
!   $Arguments:
!
 
type(t_reaction), intent(inout)                             :: this       ! Type reaction variable. 

integer, intent(in)                                         :: ntemp      ! Number of temperatures 

integer, intent(in)                                         :: ncoefflogk ! Number of polinomial coefficients for log K=f(T)

integer, intent(in)                                         :: nsp        ! Number of species associated to reaction object. 

character(len=*), intent(in)                                :: name       ! Name of the reaction.

real*8, intent(in), dimension(nsp)                          :: stq        ! Stoichiometric coefficients 

logical, intent(in)                                         :: iskin      ! iskin=true, the reaction is considered in kinetic

real*8, intent(in), dimension(ntemp)                        :: logk       ! Logaritm of equilibrum constants [ntemp]

real*8, intent(in), dimension(ntemp)                        :: temp       ! Temperature corresponding to vector of equilibrium constants [ntemp]

character(len=*), intent(in), dimension(nsp)                :: namesp     ! Name of the species. 

type(t_species), intent(in), dimension(:), optional, target :: species    ! Species objects (optional)

type(t_reactionratelaw), optional                           :: rrlaw      ! Reaction rate law object

logical, intent(out)                                        :: iserror    ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer       :: &
 i, &
 j, &
 sum, &
 nspext
integer               :: &
 id
logical                     :: &
 havespecies, &
 haverrlaw, &
 isrepeated
character(len=100)          :: &
 msg, &
 namerepeated 
!-------------------------------------------------------------------------
!
!   $code
!

iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check optional arguments
!%-------------------------------------------------------------
havespecies = present(species)
haverrlaw = present (rrlaw)
!%-------------------------------------------------------------
if(name==' ') then
 msg='Error, not defined name of the reaction object'
 goto 10
end if 
!%-------------------------------------------------------------
!% Check if the species are repeated 
!%-------------------------------------------------------------
if (this%numsp>0) then
 call find_repeated_(namesp,isrepeated,this%numsp,namerepeated)
 if (isrepeated) then
  msg='Error, species repeated:'
  call add_ (msg,namerepeated)
  goto 10
 end if
end if
!%-------------------------------------------------------------
if (ntemp>1) then
 this%numcoefflogk=ncoefflogk
 call check_pointer_ (this%coefflogk,this%numcoefflogk,.true.)
end if
!%-------------------------------------------------------------
this%name=name
this%iskin=iskin
this%numsp=nsp
allocate(this%pspecies(this%numsp))
call check_pointer_ (this%islocksps,this%numsp,.true.)
call check_pointer_ (this%stq,this%numsp,.true.)
call check_pointer_ (this%stqaqpri,this%numsp-1,.true.)
this%stq=stq
!%-------------------------------------------------------------
! If it has kinetic rate law
!%-------------------------------------------------------------
if (iskin) then
    if(.not.haverrlaw) goto 10
    if (this%islockrrlaw) then
	  call destroy_ (this%prrlaw)
      deallocate (this%prrlaw)
	  this%islockrrlaw=.false. 
    end if
    this%islockrrlaw=.true. 
	this%prrlaw => null ()  
    allocate (this%prrlaw)
    call create_ (this%prrlaw)
    this%prrlaw = rrlaw
	if (ntemp==1) call update_ (this%prrlaw,temp(1),iserror)
	if (iserror) goto 10  
end if
!%--------------------------------------------------------------
select case (havespecies)
case (.false.) 
 
 this%islocksps=.true.
 do i=1,this%numsp
   allocate (this%pspecies(i)%ptr)
   call create_(this%pspecies(i)%ptr)
   call set_name_ (this%pspecies(i)%ptr,namesp(i))
   if(namesp(i)==this%name) this%ithsecsp = i
 end do
 
case default
 this%islocksps=.false. 
 nspext = size(species)
 sum=0
 do i=1,this%numsp
  do j=1,nspext
     if(species(j)%name==namesp(i)) then
        sum=sum+1
        this%pspecies(i)%ptr => species(j)
     end if
  end do
  if(this%pspecies(i)%ptr%name==this%name) this%ithsecsp = i
 end do
 if (sum/=this%numsp) then
  msg='Error, not defined all species in the species vector'
  goto 10
 end if
end select
!%------------------------------------------------------------
if (this%ithsecsp==0) goto 10
!%------------------------------------------------------------
sum=0
do i=1,this%numsp
 
  if(this%pspecies(i)%ptr%name/=this%name) then
   sum=sum+1
   if (sum>this%numsp-1) goto 10
   this%stqaqpri(sum)=this%stq(i)
  end if
 
end do
!%-------------------------------------------------------------
!% Compute log(k)=f(T,a1,a2,a3,...ai) and storage the 
!% information 
!%-------------------------------------------------------------
call set_logk_(this,logk,temp,ntemp,iserror)
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Reaction:'
print *,'Name:',this%name
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
subroutine update_temp_param_reaction &
  (this, &
   temp, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update all attributtes in the reaction object that 
! depend of the temperature 
!
!   $Arguments:
!
 
type(t_reaction), intent(inout) :: this      ! Type reaction variable. 

real*8, intent(in)              :: temp      ! Temperature [celcius] 

logical, intent(out)            :: iserror   ! iserror=true, there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8             :: &
 tk
integer            :: &
 i 
character(len=100) :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%------------------------------------------------------------
!% If the tempref=temp, then return 
!%------------------------------------------------------------
if (this%tempref==temp) return
!%------------------------------------------------------------
!% Copy the temperature 
!%------------------------------------------------------------
this%tempref=temp 
!%------------------------------------------------------------
!% Compute the equilibrium constant according polynomial 
!% expression 
!%------------------------------------------------------------
if (this%numcoefflogk>0) then
  tk=temp+273.15d0
  this%lnk=this%coefflogk(1)*dlog(tk)+this%coefflogk(2)+ &
  this%coefflogk(3)*tk+this%coefflogk(4)/tk+ &
  this%coefflogk(5)/tk/tk
  !%----------------------------------------------------------
  !% Change log10 K => ln K
  !%----------------------------------------------------------
  this%lnk = 10.0d0**this%lnk 
  this%lnk = dlog(this%lnk)
end if
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Update parameters in the reaction rate law object 
!%------------------------------------------------------------
if (this%iskin) then
  call update_ (this%prrlaw,temp,iserror)
  if (iserror) goto 10  
end if
!%------------------------------------------------------------
return
 
10 continue
print *,'*******************************'
print *,'Reaction:'
print *,'Name:', this%name
print *,'Service: update_'
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
subroutine compute_x_reaction &
  (this, &
   x, &
   ci, &
   gi, &
   g, &
   iserror, &
   sk)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Solve the explicitly exprexssion of the mass action equation.
! A reaction at equilibrium can be expressed by the mass action equation
!                
! xj = Kj^(-1) gj^(-1) MULT(Ns)(ci gi)^(stqji) MULT(Nk) (sk)^(stqadsjk)   (1)
!
!where g are activity coefficients, and c and x are concentrations of the participating
!species. It is always possible to write the reactions in terms of a
!subset of c called primary species (dimension [Ns-Nr], where Nr is the number
!of reactions), so that the secondary species can be obtained explicitly
!from equation (1)
!
!   $Arguments:
!
 
type(t_reaction), intent(in)                  :: this      ! Type reaction variable.

real*8,intent(out)                            :: x         ! Secondary calculation.

real*8, intent(in)                            :: g         ! Activity coefficient of the jth species.

real*8, intent(in), dimension(:)              :: gi        ! Vector of the activity coefficients of the reacting species [npri]

real*8, intent(in), dimension(:)              :: ci        ! Vector of concentrations of the reacting species [npri]

logical, intent(out)                          :: iserror   ! iserror=true, then there was an error.

real*8, intent(in), optional, dimension(:)    :: sk        ! Vector of concentrations of the additional primary species (for adsorption)
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
 ndim1,i, &
 ndim2
logical                             :: &
 havesk
real*8                              :: &
 adsterm, &
 value 
character(len=100)                  :: &
 msg 
real*8, parameter                   :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%---------------------------------------------------
!% Check the number of species 
!%---------------------------------------------------
if(this%numsp==0) then
 msg='Error, not defined species'
 goto 10
end if
!%---------------------------------------------------
!% Initialice variables 
!%---------------------------------------------------
adsterm=r0
!%---------------------------------------------------
!% Check optional arguments 
!%---------------------------------------------------
havesk=present(sk)
!%---------------------------------------------------
!% Check negative concentrations in c1 vector
!%---------------------------------------------------
value=minval(ci)
if (value<r0) then
  msg='Error, negative concentration in reactive species'
  goto 10 
end if
!%---------------------------------------------------
!% Check negative activity coefficient in g1 vector
!%---------------------------------------------------
value=minval(gi)
if (value<r0) then
  msg='Error, negative activity coefficient in reactive species'
  goto 10 
end if
!%---------------------------------------------------
!% Check dimension of stqaqpri vector 
!%---------------------------------------------------
ndim1=size(this%stqaqpri)
ndim2=size(ci)
if (ndim1/=ndim2) then
 msg='Error in number of reacting species'
 goto 10
end if
!%---------------------------------------------------
!% 
!%---------------------------------------------------
if (associated(this%stqadspri)) then
 ndim1=size(this%stqadspri)
else
 ndim1=0
end if
!%---------------------------------------------------
!% Compute adsorption term
!%---------------------------------------------------
if (havesk) then
  ndim2=size(sk)
  if (ndim1/=ndim2) then
    msg='Error in number of additional primary species'
    goto 10
  end if
  adsterm=dot_product(this%stqadspri,dlog(sk))
else
  if (ndim1/=0) then
    msg='Error, not defined aditional primary species (sk)'
    goto 10
  end if
end if
!%----------------------------------------------------
!% Compute xj
!%----------------------------------------------------
x = dot_product (this%stqaqpri,dlog(ci*gi)) - dlog(g) - this%lnk + adsterm
x = dexp(x)
!%---------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Reaction:'
print *,'Name:',this%name
print *,'Service: compute_x_'
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
subroutine compute_dx_dsk_reaction &
  (this, &
   dxj, &
   xj, &
   sk, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_reaction), intent(in)   :: this

real*8, pointer, dimension(:)  :: dxj

real*8, intent(in)             :: xj

real*8, intent(in),dimension(:):: sk

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
integer                        :: &
 ndim1, &
 ndim2, &
 i
real*8, pointer                :: &
 prod(:,:) => null ()
character(len=100)             :: &
 msg 
real*8, parameter              :: &
 r1=1.0d0
!-------------------------------------------------------------------------
!
!   $code
!
 

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
ndim1=size(sk)
ndim2=size(this%stqadspri)
if (ndim1==0) return
!%----------------------------------------------------------
!% Check dimmensions 
!%----------------------------------------------------------
if (ndim1/=ndim2) then
 msg='Error in number of additional primary species'
 goto 10
end if
!%----------------------------------------------------------
!% Allocate pointers 
!%----------------------------------------------------------
call check_pointer_ (dxj,ndim2,.true.)
call check_pointer_ (prod,ndim2,ndim2,.true.)
!%----------------------------------------------------------
do i=1,ndim2
 prod(i,i)=r1/sk(i)
end do
!%----------------------------------------------------------
dxj=matmul(transpose(prod),this%stqadspri)
dxj=xj*dxj
!%----------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------
call check_pointer_ (prod,1,1,.false.)
!%----------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Reaction:'
print *,'Name:',this%name
print *,'Service: compute_dx_'
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
subroutine compute_rk_reaction &
  (this, &
   rk, &
   c1, &
   g1, &
   c, &
   g, &
   npri, &
   nsp, &
   namesp, &
   cold, &
   alpha, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute reaction rate using reaction rate law
!
!   $Arguments:
!
 
type(t_reaction), intent(in)                 :: this    ! Type reaction variable 

real*8, intent(out)                          :: rk      ! Reaction rate [mol s-1 kgw-1]

integer, intent(in)                          :: npri    ! Number of primary species

integer, intent(in)                          :: nsp     ! Number of species

real*8, intent(in), dimension(npri)          :: c1      ! Concentration of primary species

real*8, intent(in), dimension(npri)          :: g1      ! Activity coefficients of primary species

real*8, intent(in), dimension(nsp)           :: c       ! Concentration of species
   
real*8, intent(in), dimension(nsp)           :: g       ! 

real*8, intent(in)                           :: alpha   ! Reactive surface  

character(len=*), intent(in), dimension(nsp) :: namesp  ! Name of the species 

real*8, intent(in)                           :: cold     

logical, intent(out)                         :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                        :: &
 omega
character(len=100)            :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%--------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------
!% Zeroing reaction rate 
!%--------------------------------------------------------
rk=0.0d0
!%--------------------------------------------------------
!% If the reaction is in equilibrium, then return 
!%--------------------------------------------------------
if (.not.this%iskin) then
 print *,'Warning: the reaction was not defined in kinetic'
 return
end if 
!%--------------------------------------------------------
!% Check if the reaction rate law pointer is associated
!%--------------------------------------------------------
if (.not.associated(this%prrlaw)) then
 msg='Error, not associated reaction rate law in the reaction object'
 goto 10
end if
!%---------------------------------------------------------------
!% Compute omega
!%---------------------------------------------------------------
call compute_x_(this,omega,c1,g1,1.0d0,iserror)
if (iserror) then 
    msg='Error when calling compute_x_'
    goto 10
end if 
!%---------------------------------------------------------------
!% Compute rk 
!%---------------------------------------------------------------
call compute_rk_(this%prrlaw,rk,omega,c,g,nsp,cold,alpha,iserror)
!%---------------------------------------------------------------
if (iserror) then
   msg='Error when calling compute_rk_'
   goto 10 
end if
!%---------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Reaction:'
print *,'Name:',this%name
print *,'Service: compute_rk_'
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
subroutine compute_drk_dx_reaction &
  (this, &
   drk, &
   c1, &
   g1, &
   dc1, &
   dg1, &
   c, &
   dc, &
   g, &
   dg, &
   namesp, &
   npri, &
   nsp, &
   ndimder, &
   alpha, &
   isbe, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute reaction rate derivatives. 
!
!   $Arguments:
!
 
type(t_reaction), intent(in)                     :: this
 
real*8, pointer, dimension(:)                    :: drk 

integer, intent(in)                              :: npri

integer, intent(in)                              :: nsp

integer, intent(in)                              :: ndimder

real*8, intent(in), dimension(npri)              :: c1

real*8, intent(in), dimension(npri)              :: g1

real*8, intent(in), dimension(nsp)               :: c 

real*8, intent(in), dimension(nsp)               :: g

real*8, intent(in), dimension(nsp,ndimder)       :: dc

real*8, intent(in), dimension(nsp,ndimder)       :: dg

real*8, intent(in), dimension(npri,ndimder)      :: dc1

real*8, intent(in), dimension(npri,ndimder)      :: dg1

character(len=*), intent(in), dimension(nsp)     :: namesp

logical, intent(out)                             :: iserror

real*8, intent(in)                               :: alpha 

logical, intent(in)                              :: isbe 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                  :: &
 omega
real*8, pointer         :: &
 domega(:) => null (), &
 dgomega(:) => null ()
character(len=100)      :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call check_pointer_ (drk,ndimder,.true.)
!%------------------------------------------------------------
if (.not.this%iskin) return
!%------------------------------------------------------------
!% Check if the reaction rate law pointer is associated 
!%------------------------------------------------------------
if (.not.associated(this%prrlaw)) then
 msg='Error, not associated reaction rate law'
 goto 20
end if
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
call check_pointer_ (dgomega,ndimder,.true.)
!%------------------------------------------------------------
!% Compute omega
!%------------------------------------------------------------
call compute_x_(this,omega,c1,g1,1.0d0,iserror)
if (iserror) goto 20
!%------------------------------------------------------------
!% Compute derivatives of omega
!%------------------------------------------------------------
call compute_dx_(this,domega,c1,g1,1.0d0,dc1,dg1,dgomega,iserror,xj=omega)
if (iserror) then
   msg='Error when calling compute_dx_'
   goto 20
end if
!%------------------------------------------------------------
call compute_drk_(this%prrlaw,drk,omega,domega,c,dc,g,dg, &
                  nsp,ndimder,alpha,isbe,iserror)
if (iserror) then
   msg='Error when calling compute_drk_'
   goto 20
end if
!%------------------------------------------------------------
20 continue
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (domega,1,.false.)
call check_pointer_ (dgomega,1,.false.)
if (iserror) goto 10 
!%------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Reaction:'
print *,'Name:',this%name
print *,'Service: compute_drk_'
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
subroutine compute_dx_reaction &
  (this, &
   dxj, &
   c1, &
   g1, &
   gj, &
   dc1, &
   dg1, &
   dgj, &
   iserror, &
   xj, &
   dtemp, &
   sk, &
   dsk)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute derivatives of the xj with respect to concentration of primary species. 
!
!   $Arguments:
!
 
type(t_reaction), intent(in)                        :: this        ! Type reaction variable 

real*8, pointer, dimension(:)                       :: dxj         ! Derivatives 

real*8, intent(in), dimension(:)                    :: g1          ! Activity coefficients vector 

real*8, intent(in), dimension(:)                    :: c1          ! Concentration vector 

real*8, intent(in), dimension(:)                    :: dgj

real*8, intent(in), dimension(:,:)                  :: dg1

real*8, intent(in), dimension(:,:)                  :: dc1

real*8, intent(in)                                  :: gj

logical, intent(out)                                :: iserror     ! iserror=true, there was an error. 

real*8, intent(in), optional                        :: xj 

real*8, intent(in), dimension(:), optional          :: dtemp 

real*8, intent(in), dimension(:), optional          :: sk

real*8, intent(in), dimension(:,:), optional        :: dsk
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer :: &
 prod(:,:) => null ()
integer         :: &
 i, &
 ndim1, &
 ndim2
real*8          :: &
 xjloc
logical         :: &
 havexj, &
 havedtemp, &
 havesk, &
 havedsk  
character(len=100) :: &
 msg 
real*8, parameter  :: &
 r1=1.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------------
!% Check optional arguments 
!%-----------------------------------------------------------------
havexj=present(xj)
havedtemp=present(dtemp)
havesk=present(sk)
havedsk=present(dsk)
!%-----------------------------------------------------------------
ndim1=size(this%stqaqpri)
ndim2=size(c1)
if (ndim1/=ndim2) then
 msg='Error in number of reacting species' 
 goto 10
end if
!%-----------------------------------------------------------------
ndim1=size(dc1,1)
ndim2=size(dc1,2)
call check_pointer_ (dxj,ndim2,.true.)
!%-----------------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (prod,ndim1,ndim2,.true.)
!%----------------------------------------------------------------
if (havexj) then
 xjloc=xj
else
 call compute_x_(this,xjloc,c1,g1,gj,iserror)
 if (iserror) goto 20
end if
!%-----------------------------------------------------------------
do i=1,ndim2
 prod(i,:) = dc1(i,:)/c1(i) + dg1(i,:)/g1(i)
end do
!%-----------------------------------------------------------------
dxj = (matmul(transpose(prod),this%stqaqpri) - dgj/gj)
!%-----------------------------------------------------------------
!% Only if sk and dsk are present 
!%-----------------------------------------------------------------
if (havesk.and.havedsk) then
 ndim1=size(this%stqaqpri)
 ndim2=size(dsk,2)
 if (ndim1/=ndim2) then
   msg='Error in number of reacting species in dsk/dc1' 
   goto 10
 end if
 ndim1=size(dsk,1)
 ndim2=size(dsk,2)
 call check_pointer_ (prod,ndim1,ndim2,.true.)
 do i=1,ndim1
   prod(i,:) = dsk(i,:)/sk(i)
 end do
 dxj = dxj + matmul(transpose(prod),this%stqadspri)
end if 
!%-----------------------------------------------------------------
!% 
!%-----------------------------------------------------------------
dxj = xjloc * dxj
!%-----------------------------------------------------------------
20 continue 
!%-----------------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (prod,1,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Reaction:'
print *,'Name:',this%name
print *,'Service: compute_dx_'
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
subroutine compute_dx_numeric_reaction &
  (this, &
   dx, &
   ci, &
   gi, &
   g, &
   dci, &
   factor, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute numeric derivative 
!
!   $Arguments:
!
 
type(t_reaction), intent(in)        :: this ! Type reaction variable. 

real*8, pointer, dimension(:)       :: dx

real*8, intent(in), dimension(:)    :: gi

real*8, intent(in), dimension(:)    :: ci

real*8, intent(in), dimension(:,:)  :: dci

real*8, intent(in)                  :: g       ! Activity coefficient  

real*8, intent(in)                  :: factor  ! Perturbation factor. 

logical, intent(out)                :: iserror ! iserror=true, there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer :: &
 ciperturbed(:) => null (), &
 dxloc(:) => null ()
integer         :: &
 isps, &
 ndim1, &
 ndim2
real*8          :: &
 xperturbed, &
 x, &
 delta 
character(len=100) :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------------
ndim1=size(this%stqaqpri)
ndim2=size(ci)
!%-----------------------------------------------------------------
!% Check the number of product species 
!%-----------------------------------------------------------------
if (ndim1/=ndim2) then
 msg='Error in the number of reacting species' 
 goto 10
end if
!%-----------------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (ciperturbed,ndim1,.true.)
call check_pointer_ (dxloc,ndim1,.true.)
!%-----------------------------------------------------------------
!% Allocate dx 
!%-----------------------------------------------------------------
ndim1=size(dci,1)
ndim2=size(dci,2)
call check_pointer_ (dx,ndim2,.true.)
!%----------------------------------------------------------------
!% Compute xj (not perturbed)
!%----------------------------------------------------------------
call compute_x_(this,x,ci,gi,g,iserror)
!%----------------------------------------------------------------
do isps=1,ndim1
 delta=ci(isps)*factor
 ciperturbed(isps)=ci(isps)+delta
 call compute_x_(this,xperturbed,ciperturbed,gi,g,iserror)
 if (iserror) goto 20
 dxloc(isps)=(xperturbed-x)/delta
end do
!%-----------------------------------------------------------------
!% 
!%-----------------------------------------------------------------
dx=matmul(dxloc,dci)
!%-----------------------------------------------------------------
20 continue 
!%-----------------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (ciperturbed,1,.false.)
call check_pointer_ (dxloc,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Reaction:'
print *,'Name:',this%name
print *,'Service: compute_dx_'
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
subroutine get_namesp_reaction &
  (this, &
   namesp, &
   nsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of species. 
!
!   $Arguments:
!
 
type(t_reaction), intent(in)               :: this

character(len=*), pointer, dimension(:)    :: namesp

integer, intent(out)                       :: nsp 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                      :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%-----------------------------------------------------------
if (this%numsp>0) then
  nsp=this%numsp
!%-----------------------------------------------------------
  call check_pointer_ (namesp,this%numsp,.true.)
!%-----------------------------------------------------------
  do i=1,this%numsp
    namesp (i) = this%pspecies(i)%ptr%name
  end do
else
  nsp=0
  call check_pointer_ (namesp,1,.false.)
end if
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_numsp_reaction &
  (this, &
   numsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the number of species in the reaction object 
!
!   $Arguments:
!
 
type(t_reaction), intent(in)        :: this

integer, intent(out)                :: numsp 
 
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
numsp = this%numsp
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_stq_reaction &
  (this, &
   stq)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the stoichiometric matric corresponding to primary aqueous species. 
!
!   $Arguments:
!
type(t_reaction), intent(in)  :: this

real*8, pointer, dimension(:) :: stq
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer          :: &
  ndim 
!-------------------------------------------------------------------------
!
!   $code
!'
!%-------------------------------------------------------------
if (this%numsp>0) then
!%-------------------------------------------------------------
 ndim=size(this%stqaqpri)
!%-------------------------------------------------------------
 if (ndim>0) then
  call check_pointer_ (stq,ndim,.true.)
  stq=this%stqaqpri
 end if
end if
!%------------------------------------------------------------
return 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_lnk_reaction &
  (this, &
   lnk)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return logarithm of the equilibrium constant.
!
!   $Arguments:
!
 
type(t_reaction), intent(in)        :: this

real*8, intent(out)                 :: lnk 
 
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
lnk = this%lnk
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_logk_temp_reaction &
  (this, &
   logktemp, &
   templogk, &
   ntemp, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the logaritm of equilibrium constant to 
! different temperatures 
!
!   $Arguments:
!
 
type(t_reaction), intent(in)        :: this

real*8, pointer, dimension(:)       :: logktemp 

real*8, pointer, dimension(:)       :: templogk

integer, intent(out)                :: ntemp  ! Number of temperatures

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
iserror=.false. 
!%------------------------------------------------------------
if (this%numtemp>0) then
 ntemp=this%numtemp
 call check_pointer_ (logktemp,ntemp,.true.)
 call check_pointer_ (templogk,ntemp,.true.)
 logktemp=this%logktemp
 templogk=this%templogk
end if  
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_if_sp_is_present_reaction &
  (this, &
   namesp, &
   isbe)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return if the species is invlved in the reaction. 
!
!   $Arguments:
!
 
type(t_reaction), intent(in), target:: this

character(len=*), intent(in)        :: namesp

logical, intent(out)                :: isbe 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                 :: &
 i 
character(len=100), pointer     :: &
 name => null ()
 
!-------------------------------------------------------------------------
!
!   $code
!
isbe=.false.
!%------------------------------------------------------------
do i=1,this%numsp
  name => this%pspecies(i)%ptr%name
  if(name==namesp) then
    isbe=.true.
	name => null () 
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
subroutine copy_reaction &
  (targetobj, &
   sourceobj)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy the reaction object in other reaction object. 
!
!   $Arguments:
!
 
type (t_reaction), intent(in) :: sourceobj  ! Type reaction variable.

type (t_reaction), intent(out):: targetobj  ! Type reaction variable.
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: The target object must be previously created
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                       :: &
 ndim, &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------
targetobj%name=sourceobj%name
!%------------------------------------------------
targetobj%ithsecsp=sourceobj%ithsecsp
!%------------------------------------------------
targetobj%lnk = sourceobj%lnk
!%------------------------------------------------
targetobj%numsp=sourceobj%numsp
!%------------------------------------------------
targetobj%numtemp=sourceobj%numtemp
!%------------------------------------------------
targetobj%numcoefflogk=sourceobj%numcoefflogk
!%------------------------------------------------
targetobj%iskin=sourceobj%iskin
!%------------------------------------------------
targetobj%tempref=sourceobj%tempref
!%------------------------------------------------
targetobj%pressref=sourceobj%pressref
!%------------------------------------------------
targetobj%islockrrlaw=sourceobj%islockrrlaw
!%------------------------------------------------
if (associated(sourceobj%islocksps)) then
 ndim=size(sourceobj%islocksps)
 call check_pointer_ (targetobj%islocksps,ndim,.true.)
 targetobj%islocksps = sourceobj%islocksps
end if
!%------------------------------------------------
nullify(targetobj%pspecies)
allocate(targetobj%pspecies(sourceobj%numsp))
do i=1,sourceobj%numsp
  
  if (sourceobj%islocksps(i)) then
   allocate (targetobj%pspecies(i)%ptr)
   call create_ (targetobj%pspecies(i)%ptr)
   targetobj%pspecies(i)%ptr = sourceobj%pspecies(i)%ptr
  else
   targetobj%pspecies(i)%ptr => sourceobj%pspecies(i)%ptr
  end if 

end do
!%------------------------------------------------
if (sourceobj%numcoefflogk>0) then
 call check_pointer_ (targetobj%coefflogk,targetobj%numcoefflogk,.true.)
 targetobj%coefflogk = sourceobj%coefflogk
end if
!%------------------------------------------------
if (sourceobj%numtemp>0) then
 call check_pointer_ (targetobj%templogk,targetobj%numtemp,.true.)
 call check_pointer_ (targetobj%logktemp,targetobj%numtemp,.true.)
 targetobj%templogk = sourceobj%templogk
 targetobj%logktemp = sourceobj%logktemp
end if
!%-------------------------------------------------
if (associated(sourceobj%stqaqpri)) then
 ndim=size(sourceobj%stqaqpri)
 call check_pointer_ (targetobj%stqaqpri,ndim,.true.)
 targetobj%stqaqpri = sourceobj%stqaqpri
end if
!%-------------------------------------------------
if (associated(sourceobj%stqadspri)) then
 ndim=size(sourceobj%stqadspri)
 call check_pointer_ (targetobj%stqadspri,ndim,.true.)
 targetobj%stqadspri = sourceobj%stqadspri
end if
!%-------------------------------------------------
if (associated(sourceobj%stq)) then
 ndim=size(sourceobj%stq)
 call check_pointer_ (targetobj%stq,ndim,.true.)
 targetobj%stq = sourceobj%stq
end if
!%-------------------------------------------------
targetobj%prrlaw => null ()
!%-------------------------------------------------
if (sourceobj%iskin) then
 
 if (sourceobj%islockrrlaw) then 
  allocate (targetobj%prrlaw)
  call create_ (targetobj%prrlaw)
  targetobj%prrlaw=sourceobj%prrlaw
 else
  targetobj%prrlaw => sourceobj%prrlaw 
 end if 

end if
!%------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_name_reaction &
  (this, &
   name)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of the reaction object. 
!
!   $Arguments:
!
 
type(t_reaction), intent(in)     :: this    ! Type reaction variable.

character(len=100), intent(out)  :: name    ! Name of the reaction  
 
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
name = this%name
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_name_sec_reaction &
  (this, &
   name)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of the secondary species 
!
!   $Arguments:
!
 
type(t_reaction), intent(in)     :: this    ! Type reaction variable.

character(len=100), intent(out)  :: name 
 
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
if (this%ithsecsp>0.and.this%ithsecsp<this%numsp) then 
 name = this%pspecies(this%ithsecsp)%ptr%name
else
 name=''
end if 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_pspecies_reaction &
  (this, &
   species)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the pointer to species object. If the species object 
! was previously allocated for reaction object then destroy and deallocate the 
! species object.
!
!   $Arguments:
!
 
type(t_reaction), intent(inout)         :: this    ! Type reaction variable.

type(t_species), intent(in), target     :: species ! Type species object. 
 
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
 isp
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
do isp=1,this%numsp
 if (this%pspecies(isp)%ptr%name==species%name) then
  if (this%islocksps(isp)) then 
    call destroy_ (this%pspecies(isp)%ptr)
    deallocate (this%pspecies(isp)%ptr)
	this%pspecies(isp)%ptr => null () 
	this%islocksps(isp)=.false.  
  end if 
  this%pspecies(isp)%ptr => species     
  exit
 end if
end do
!%------------------------------------------------------------
!% Set the species in the reaction rate law object 
!% (only if the reaction is kinetic)
!%------------------------------------------------------------
if (this%iskin) then
 call set_pspecies_ (this%prrlaw,species)
end if
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_coeff_stq_reaction &
  (this, &
   namesp, &
   coeff)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Get stoichiometic coefficient corresponding of any species
!%  according the name of species
!%  If the species not participing in the reaction then coeff= 0
!
!   $Arguments:
!
 
type(t_reaction), intent(in), target  :: this

character(len=*), intent(in)          :: namesp

real*8, intent(out)                   :: coeff 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                 :: &
  i 
character(len=100), pointer :: &
 name => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%----------------------------------------------------------
coeff=0.0d0
!%----------------------------------------------------------
do i=1,this%numsp
 name => this%pspecies(i)%ptr%name
 if(name==namesp) then
   coeff = this%stq(i)
   name => null ()
   exit  
 end if
end do
!%----------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_polinom_logk_reaction &
  (this, &
   coefflogk, &
   ncoeff)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Retrn the polinomial term of the log k =f(T)
!
!   $Arguments:
!
 
type(t_reaction), intent(in)         :: this

real*8, pointer, dimension(:)        :: coefflogk

integer, intent(out)                 :: ncoeff 
 
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
ncoeff=0
if (this%numcoefflogk>0) then
 ncoeff=this%numcoefflogk
 call check_pointer_ (coefflogk,ncoeff,.true.)
 coefflogk=this%coefflogk
end if
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_namesp_rrlaw_reaction &
  (this, &
   namesp, &
   nsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the species in the reaction rate law
!
!   $Arguments:
!
 
type(t_reaction), intent(in)            :: this

character(len=*), pointer, dimension(:) :: namesp

integer, intent(out)                    :: nsp 
 
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
if (this%iskin) then
 call get_namesp_ (this%prrlaw,namesp,nsp)
else
 nsp=0
 call check_pointer_ (namesp,1,.false.) 
end if
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_prrlaw_reaction &
  (this, &
   prrlaw)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the pointer to reaction rate law object
!
!   $Arguments:
!
 
type(t_reaction), intent(in)     :: this     ! Type reaction object 

type(t_reactionratelaw), pointer :: prrlaw   ! Type reaction rate law object
 
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
prrlaw => this%prrlaw
!%-----------------------------------------------------------
return
end subroutine
!%***********************************************************
!%***************Public subroutine***************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
subroutine rewrite_stq_reaction &
  (this, &
   namepri, &
   npri, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Rewrite the stoichiometric matrix according to species list.
!
!   $Arguments:
!
 
type(t_reaction), intent(inout), target       :: this     ! Type reaction variable. 

integer, intent(in)                           :: npri     ! Number primary species 

character(len=*), intent(in), dimension(npri) :: namepri  ! Name of primary species (reactin species)

logical, intent(out)                          :: iserror  ! iserror=true, there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer               :: &
 i, &
 j, &
 isp
character(len=100)    :: &
 msg 
character(len=100), pointer :: &
 name => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%--------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------
call check_pointer_ (this%stqaqpri,npri,.true.)
!%---------------------------------------------------------------
isp=0
do i=1,npri
 
  do j=1,this%numsp
    name => this%pspecies(j)%ptr%name
    if(name==namepri(i)) then
       isp=isp+1
       this%stqaqpri (i) = this%stq(j)
	   exit  
    end if
 
  end do
 
  if(namepri(i)==this%name) then
   msg='Error, the name of reaction is equal to one primary species'
   goto 10
  end if
 
end do
!%---------------------------------------------------------------
name => null () 
!%--------------------------------------------------------------- 
return
 
10 continue 
print *,'*******************'
print *,'Reaction:'
print *,'Name:', this%name
print *,'Service: rewrite_'
print *, msg
print *,'*******************'
iserror=.true.
return
 
end subroutine
!%***********************************************************
!%***************Public subroutine***************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
subroutine rewrite_rrlaw_reaction &
  (this, &
   namesp, &
   nsp, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Rewrite the stoichiometric matrix of catalytic species 
! in the reaction rate law object according to species list (only if the reaction is kinetic)
!
!   $Arguments:
!
 
type(t_reaction), intent(inout), target       :: this    ! Type reaction variable. 

integer, intent(in)                           :: nsp     ! Number of species.

character(len=*), intent(in), dimension(nsp)  :: namesp  ! Name of the species.

logical, intent(out)                          :: iserror ! iserror=true, then there was an error. 
 
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
!%--------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Only if the reaction was defined in kinetic. 
!%--------------------------------------------------------------
if (this%iskin) then 
 call rewrite_ (this%prrlaw,namesp,nsp,iserror)
end if 
!%--------------------------------------------------------------- 
return
 
10 continue 
print *,'***********************'
print *,'Reaction:'
print *,'Name:', this%name
print *,'Service: rewrite_rrlaw_'
print *, msg
print *,'***********************'
iserror=.true.
return
 
end subroutine
!%***********************************************************
!%***************Public subroutine***************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
subroutine set_stqads_reaction &
  (this, &
   namexoh, &
   numsite, &
   nameadsmodel, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the stoichimetric matrix corresponding to sorption 
! primary species according to sorption model. 
!
!   $Arguments:
!
 
type(t_reaction), intent(inout), target           :: this          ! Type reaction variable 

integer, intent(in)                               :: numsite       ! Number of sites associated to any surface

character(len=*), intent(in), dimension(numsite)  :: namexoh       ! Name of sorption primary species [numsite]

character(len=*), intent(in)                      :: nameadsmodel  ! Name of adosorption model. 

logical, intent(out)                              :: iserror       ! iserror=true, there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer               :: &
 i, &
 j, &
 isite, &
 ipri
real*8                :: &
 z
integer               :: &
 numsk
logical               :: &
 isbexoh
character(len=100), pointer :: &
 name => null ()
character(len=100)    :: &
 msg
character(len=100), parameter :: &
 nameprop='charge'
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------
isbexoh=.false.
!%-------------------------------------------------------
!% Conmpute dimension of stqadspri [numsk] according type of
!% sorption model  
!%-------------------------------------------------------
select case (nameadsmodel)
!%-----------------------
case ('cationexchange','langmuir')
 numsk=numsite
!%-----------------------
case ('triplelayer')
 numsk=numsite*4
!%-----------------------
case ('diffuselayer','constantcapacitance')
 numsk=numsite*2
!%-----------------------
case ('freundlich','kd','equilibrium')
 return
case default
 msg='Error, not recognized sorption model:'
 call add_ (msg,nameadsmodel)
 goto 10
end select
!%-------------------------------------------------------
!% Allocate stqadspri vector
!%-------------------------------------------------------
call check_pointer_ (this%stqadspri,numsk,.true.)
!%-------------------------------------------------------
!% Rewrite xoh primary species
!%-------------------------------------------------------
ipri=1
do j=1,numsite
 do i=1,this%numsp
  name => this%pspecies(i)%ptr%name
  if (name==namexoh(j)) then
    this%stqadspri(ipri)=this%stq(i)
    isbexoh=.true.
    isite=j
    exit
  end if
 end do
!%-----------------
 select case (nameadsmodel)
 case ('cationexchange','langmuir')
  ipri=ipri+1
 case ('triplelayer')
  ipri=ipri+4
 case ('diffuselayer','constant capacitance')
  ipri=ipri+2
 end select
!%-----------------
end do
if (.not.isbexoh) then
  msg='Error, sorption primary not defined in the reaction'
  iserror=.true. 
  goto 20
end if
!%-----------------------------------------------------
!%-----------------------------------------------------
!%-----------------------------------------------------
select case (nameadsmodel)
!%-------------------------
case ('triplelayer')
 ipri=(isite-1)*4
 
 do i=1,this%numsp
    name => this%pspecies(i)%ptr%name
    if(name=='h+'.or.name=='H+') then
      this%stqadspri(ipri+2)=this%stq(i)
    else if((name/='h+'.or.name/='H+').and.name/=this%name.and.name/=namexoh(isite)) then
      call get_prop_ (this%pspecies(i)%ptr,z,'charge',msg,iserror)
      !call get_prop_ (this%pspecies(this%ithsecsp)%ptr,z,'charge',msg,iserror)
      if (iserror) goto 20
      this%stqadspri(ipri+3)=z
    end if
 end do
!%-------------------------
case ('diffuselayer','constantcapacitance')
!%--------------------------------------------------
!% Add the charge in the stoichiometric matrix 
!% of adsorption reactions
!%--------------------------------------------------
 ipri=(isite-1)*2
 call get_prop_ (this%pspecies(this%ithsecsp)%ptr, z,'charge',msg,iserror)
 if (iserror) goto 20 
 this%stqadspri(ipri+2)=z
end select
!%-----------------------------------------------------
20 continue 
!%-----------------------------------------------------
!% Nullify local pointer 
!%-----------------------------------------------------
name => null ()
!%-----------------------------------------------------
if (iserror) goto 10 
!%-----------------------------------------------------
return
 
10 continue 
print *,'******************'
print *,'Reaction:'
print *,'Name:', this%name
print *,'Service: rewrite_'
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
subroutine write_reaction &
  (this, &
   ioutput, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ioutput unit all attributtes encapsulated in the reaction object.
!
!   $Arguments:
!
 
type(t_reaction), intent(in)        :: this    ! Type reaction variable

integer, intent(in)                 :: ioutput ! Output unit

logical, intent(out)                :: iserror ! ierror=true, there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                 :: &
 i, &
 sum, &
 ipos, &
 sum1, &
 sum2, &
 isp1, &
 isp2
character(len=100)      :: &
 typereaction
character(len=200)      :: &
 msg, &
 string, &
 string1, &
 string2, &
 name 
real*8                  :: &
 value, &
 temp 
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
msg=''
iserror=.false.
!%-------------------------------------------------------------
string1=''
string2=''
string=''
sum=0
sum1=0
sum2=0
isp1=0
isp2=0
!%-------------------------------------------------------------
if(this%iskin) then
 typereaction='kinetic'
else
 typereaction='equilibrium'
end if
!%-------------------------------------------------------------
write (ioutput,*) "-------------------------------------------"
write (ioutput,*) "               Reaction object             "
write (ioutput,*) "-------------------------------------------"
write (ioutput,1) "type=",typereaction
write (ioutput,*)
!%------------------------------------------------------------ 
do i=1,this%numsp
 
   value = this%stq(i)
 
   name = this%pspecies(i)%ptr%name
 
   call lastletter_ (ipos,name)
 
   if (value>0.0d0) then
    isp1=isp1+1
    value=dabs(value)
    if (isp1.ne.1) then
        string1(sum1+1:sum1+3)=' + '
        sum1 = sum1 + 3
    end if
    if (value.ne.1.0d0) then
        write(unit=string1(sum1+1:sum1+4),fmt=2),value
        sum1 = sum1 + 5
    end if
    string1(sum1+1:sum1+ipos)=name(1:ipos)
    sum1 = sum1+ ipos
   else
    isp2=isp2+1
    value=dabs(value)
        if (isp2.ne.1) then
          string2(sum2+1:sum2+3)=' + '
          sum2 = sum2 + 3
        end if
      if (value.ne.1.0d0) then
        write(unit=string2(sum2+1:sum2+4),fmt=2),value
        sum2 = sum2 + 5
      end if
    string2(sum2+1:sum2+ipos)=name(1:ipos)
    sum2 = sum2 + ipos
   end if
 
end do
!%-----------------------------------------------------------------
sum = sum + sum2
string(1:sum)=string2(1:sum2)
sum = sum + 2
string(sum+1:sum+8)='<======>'
sum = sum + 10
string(sum+1:sum+sum1)=string1(1:sum1)
sum = sum + sum1
write (ioutput,*) string(1:sum)
!%-----------------------------------------------------------------
!% Write equilibrium constant
!%-----------------------------------------------------------------
 write (ioutput,*) "-------------------------------------------"
 write (ioutput,*) '            Equilibrium constant           '
 write (ioutput,*) "-------------------------------------------"
 write (ioutput,3) 'Temp=',this%tempref 
 write (ioutput,3) 'log K=',dlog10(exp(this%lnk))
 write (ioutput,*) "-------------------------------------------"
 write (ioutput,*) 'Equilibrium constants stored' 
 write (ioutput,*) "-------------------------------------------"
 write (ioutput,*) '      Temp         log K                   ' 
 write (ioutput,*) "-------------------------------------------"
 do i=1,this%numtemp  
   write (ioutput,5) this%templogk(i),this%logktemp(i)
 end do 
!%-----------------------------------------------------------------
!% Write polinomial coefficents of the equilibrium constant 
!%-----------------------------------------------------------------
if (this%numcoefflogk>0) then 
 write (ioutput,*) "-------------------------------------------"
 write (ioutput,*) '        Polinomial coefficients:           '
 write (ioutput,*) "-------------------------------------------"
 write (ioutput,*) 'log K(T) = at1 log T + at2 + at3 T + at4/T + at5/T2' 
 write (ioutput,*) 'where' 
 do i=1,this%numcoefflogk
  value=this%coefflogk(i)
  write (ioutput,4) 'at',i,'=',value
 end do 
end if 
!%-----------------------------------------------------------------
!% Write the kinetic reaction rate law 
!%-----------------------------------------------------------------
if (this%iskin) then
 call write_ (this%prrlaw,ioutput,iserror)
 if (iserror) goto 10  
end if
!%-----------------------------------------------------------------
write (ioutput,*) "-------------------------------------------"
!%------------------------------------------------------------- 
return
 
1 format(a5,a15)
2 format(f4.1)
3 format(a7,f7.3)
4 format(a2,i1,a1,f15.5)
5 format(2f10.4)
10 continue 
print *,'*************'
print *,'Reaction:'
print *,'Name:',this%name
print *,'Service: write_'
print *, msg
print *,'*************'
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
integer                        :: status
 
 
 
call read_xml_loc_ (name, attributes)
 
 
return
end subroutine
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc_reaction (name, attributes,this)
implicit none 
!-------------------------------------------------------------------------
!
!   $Description: Write in ioutput unit all information about the reaction
!
!   $Arguments:
!
 
type (t_reaction), optional, intent(inout)    :: this

character(len=*), intent(in)                  :: name

type(dictionary_t), intent(in)                :: attributes
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer        :: status,n,ncoeff
real*8         :: reallocal(1)
logical        :: havethis
character(len=100)      :: &
 id
integer, save           :: &
 numlogk, &
 ilogk, &
 isp, &
 iterm, &
 numsp, &
 itype, &
 typerrlaw, &
 nattrtermrrl 
real*8, save          :: &
 ea, &
 attrrrlaw1, &
 attrrrlaw2, &
 attrrrlaw3, &
 attrrrlaw4
integer, pointer, save:: &
 numspterm(:)
real*8, pointer, save :: &
 temp(:), &
 logkvalue(:), &
 coeff(:), &
 attrtermrrlaw(:,:), &
 attrsprrlaw(:,:)
character(len=100), save:: &
 namereact, &
 namerrlaw, &
 typereact
character(len=100), pointer, save:: &
 namesptermrrlaw(:,:), &
 namesp(:), &
 typetermrrlaw(:)
logical, save           :: &
 iskin, &
 isberrlaw, &
 isareadep
integer, parameter      :: &
 mxdim=20
logical                 :: &
 iserror
type(t_reactionratelaw), pointer :: &
 rrlaw 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------------
havethis=present(this)
!%------------------------------------------------------------------
if(havethis) then
 
     if(iskin.and.isberrlaw) then
         allocate (rrlaw)
         call create_ (rrlaw)
         call set_ &
        (rrlaw, &
         typerrlaw, &
         namerrlaw, &
         numspterm(1:iterm), &
         attrtermrrlaw(1:nattrtermrrl,1:iterm), &
         attrsprrlaw(1:10,1:iterm), &
         namesptermrrlaw(1:10,1:iterm), &
         typetermrrlaw(1:iterm), &
         iterm, &
         nattrtermrrl, &
         10, &
         ea, &
         attrrrlaw1, &
         attrrrlaw2, &
         attrrrlaw3, &
		 attrrrlaw4, &
         isareadep, &
         iserror)
     end if
     
	 call set_ &
    (this, &
     namereact, &
     iskin, &
     coeff(1:numsp), &
     namesp(1:numsp), &
     logkvalue(1:numlogk), &
     temp(1:numlogk), &
     numlogk, &
     5, &
     numsp, &
     iserror, &
	 rrlaw=rrlaw)

   if(iskin.and.isberrlaw) then
    call destroy_ (rrlaw)
    deallocate (rrlaw) 
   end if 
   deallocate(numspterm)
   deallocate(temp)
   deallocate(logkvalue)
   deallocate(coeff)
   deallocate(attrtermrrlaw)
   deallocate(attrsprrlaw)
   deallocate(namesptermrrlaw)
   deallocate(typetermrrlaw) 
   deallocate(namesp)
   if(iserror) goto 10
!%------------------------------------------------------------
else
 
 select case(name)
 
 
 case ('reaction')
   
   allocate(typetermrrlaw(mxdim))  
   allocate(numspterm(mxdim))
   allocate(temp(mxdim))
   allocate(logkvalue(mxdim))
   allocate(coeff(mxdim))
   allocate(namesp(mxdim))
   allocate(attrtermrrlaw(mxdim,mxdim))
   allocate(attrsprrlaw(mxdim,mxdim))
   allocate(namesptermrrlaw(mxdim,mxdim))
   isberrlaw=.false.
   ilogk=0
   isp=1
   id='' 
   call get_value (attributes,"name", id, status)
   if (status<0) goto 10
   namereact=id
   namesp(isp)=namereact
   coeff(isp)=-1.0d0
   id='' 
   call get_value (attributes,"type", id, status)
   typereact=id
   select case (id)
   case ('kinetic','KINETIC','Kinetic','kin','KIN')
    iskin=.true.
   case ('equilibrium','Equilibrium','EQUILIBRIUM','eq','EQ','')
    iskin=.false.
   case default
    goto 10  
   end select
 
 case ('logk')
   ilogk=ilogk+1
   numlogk=ilogk
   call get_value (attributes,"temp", id, status)
   n=0
   call build_data_array (id,reallocal,n)
   temp(ilogk)=reallocal(1)
   call get_value (attributes,"value", id, status)
   n=0
   call build_data_array (id,reallocal,n)
   logkvalue(ilogk)=reallocal(1)
 
 case ('species')
   isp=isp+1
   if(.not.isberrlaw) then
     numsp=isp
     id=''
     call get_value (attributes,"name", id, status)
     namesp(isp)=id
     id=''
     call get_value (attributes,"coeff", id, status)
     n=0
     call build_data_array (id,reallocal,n)
     coeff(isp)=reallocal(1)
!%---------------
   else if (isberrlaw.and.iskin) then 
     numspterm(iterm)=isp
     call get_value (attributes,"name", id, status)
     call lowercase (id)
     namesptermrrlaw(isp,iterm) = id
     call get_value (attributes,"attr1", id, status)
     n=0
     reallocal(1)=0.0d0
     call build_data_array (id,reallocal,n)
     attrsprrlaw(isp,iterm) = reallocal(1)
   end if
!%---------------
 case ('reactionratelaw')
  isberrlaw=.true. 
  isareadep=.true. 
  if (iskin) then
    isp=0
    iterm=0
    id=''
    call get_value (attributes,"name", id, status)
    call lowercase (id)
    namerrlaw=id
    id=''
    call get_value (attributes,"type", id, status)
    call lowercase (id)
!%-------------
    select case (id)
    case ('lasaga')
     typerrlaw=1
     nattrtermrrl=3
    case ('monod')
     typerrlaw=2
     nattrtermrrl=0
    case default
     typerrlaw=0
     nattrtermrrl=0
    end select
!%-------------
    call get_value (attributes,"ea", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    ea= reallocal(1)
!%--------
    id=''
    call get_value (attributes,"attr1", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    attrrrlaw1 = reallocal(1)
!%--------
    id=''
    call get_value (attributes,"attr2", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    attrrrlaw2 = reallocal(1)
!%--------
    id=''
    call get_value (attributes,"attr3", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    attrrrlaw3 = reallocal(1)
	!%--------
    id=''
    call get_value (attributes,"attr4", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    attrrrlaw4 = reallocal(1)
!%--------
    id=''
    call get_value (attributes,"areadep", id, status)
    select case (id)
    case('not','NOT','N','n')
     isareadep=.false.
    end select 
  end if
!%------------------------------------
!%------------------------------------
 case ('term')
  if (iskin) then   
	   iterm=iterm+1
       isp=0
     if (typerrlaw==1) then    ! For reaction rate law type
       typetermrrlaw(iterm)='catalyst'
       id=' '
       call get_value (attributes,"attr1", id, status)  ! rate constan
       if (status/=0) goto 10
       n=0
       reallocal(1)=0.0d0
       call build_data_array (id,reallocal,n)
       attrtermrrlaw(1,iterm) = reallocal(1)
       call get_value (attributes,"attr2", id, status)  ! theta
       if (status/=0) goto 10
       n=0
       reallocal(1)=0.0d0
       call build_data_array (id,reallocal,n)
       attrtermrrlaw(2,iterm) = reallocal(1)
       call get_value (attributes,"attr3", id, status)   ! eta
       if (status/=0) goto 10
       n=0
       reallocal(1)=0.0d0
       call build_data_array (id,reallocal,n)
       attrtermrrlaw(3,iterm) = reallocal(1)

     else if (typerrlaw==2) then  ! For reaction rate law typ
       id=' '
       call get_value (attributes,"type", id, status)
       if (status/=0) goto 10
       typetermrrlaw(iterm)=id
     end if
  end if 
 end select
 
end if
return
 
10 continue  
print *,'Reaction:'
print *,'Service: read_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine verify_reaction &
  (this, &
   name, &
   numsp, &
   species, &
   logk, &
   stq, &
   numcoefflogk, &
   idspstq, &
   ithsecsp, &
   nerror, &
   errormsg)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Verify different attributes encapsulated in the reaction object. 
!
!   $Arguments:
!
 
type(t_reaction)                            :: this

integer                                     :: nerror

integer, optional                           :: numcoefflogk
 
integer, optional                           :: ithsecsp
 
integer, optional                           :: numsp

integer, pointer, optional                  :: idspstq(:)

type(t_species), pointer, optional          :: species(:)

real*8, pointer, optional                   :: stq(:)

character(len=*), optional                  :: name

real*8, optional                            :: logk

character(len=100), pointer                 :: errormsg(:) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical                     :: &
 havename, &
 havespecies, &
 havelogk, &
 havestq, &
 havenumcoefflogk, &
 haveidspstq, &
 haveithsecsp, &
 havenumsp
integer                  :: &
 ndim, &
 i
real*8                   :: &
 error
character(len=100)       :: &
 errormsgloc(50) 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

 
if(associated(errormsg)) deallocate(errormsg)
 
nerror=0
errormsgloc=' '
 havename=present(name)
 havespecies=present(species)
 havelogk=present(logk)
 havestq=present(stq)
 havenumcoefflogk=present(numcoefflogk)
 haveidspstq=present(idspstq)
 haveithsecsp=present(ithsecsp)
 havenumsp=present(numsp)
!%---------------------------------------------
if (havename) then
 
 if (name.ne.this%name) then
  nerror=nerror+1
  errormsgloc(nerror)='error en name'
 end if
 
end if
!%---------------------------------------------
if (havenumsp) then
 
 if (numsp.ne.this%numsp) then
  nerror=nerror+1
  errormsgloc(nerror)='error en numsp'
 end if
 
end if
!%---------------------------------------------
if (havelogk) then
 
 if(this%lnk.ne.0.0) then
  error=(this%lnk-logk)/this%lnk
  if (abs(error).gt.1.0d-2) then
   nerror=nerror+1
   errormsgloc(nerror)='error en logk'
  end if
 else
  if (logk.ne.this%lnk) then
   nerror=nerror+1
   errormsgloc(nerror)='error en logk'
  end if
 end if
 
end if
!%---------------------------------------------
if (havenumcoefflogk) then
 
 if (numcoefflogk.ne.this%numcoefflogk) then
   nerror=nerror+1
   errormsgloc(nerror)='error en logk'
 end if
 
end if
!%---------------------------------------------
if (haveithsecsp) then
 
 if (ithsecsp.ne.this%ithsecsp) then
   nerror=nerror+1
   errormsgloc(nerror)='error en ithsecsp'
 end if
 
end if
!%---------------------------------------------
if (havespecies) then
 select case (associated(species))
 case (.true.)
  if (.not.associated(this%pspecies)) nerror=nerror+1
 
  ndim=size(species)
 
  if (ndim.ne.size(this%pspecies)) nerror=nerror+1
  if (ndim.ne.this%numsp) nerror=nerror+1
 
  do i=1,ndim
 
   if(species(i)%name.ne.this%pspecies(i)%ptr%name) nerror=nerror+1
 
  end do
 
 case default
 
  if(associated(this%pspecies)) nerror=nerror+1
 
 end select
end if
!%--------------------------------------------
if (havestq) then
 select case (associated(stq))
 case (.true.)
  ndim=size(stq)
  if (ndim.ne.size(this%stqaqpri)) then
   nerror=nerror+1
   errormsgloc(nerror)='error en stq'
  end if
  do i=1,ndim
   if(stq(i).ne.this%stqaqpri(i)) then
    nerror=nerror+1
    errormsgloc(nerror)='error en stq'
   end if
  end do
 
 case default
  if(associated(this%stqaqpri)) then
   nerror=nerror+1
   errormsgloc(nerror)='error en stq'
  end if
 end select
end if
 
!%---------------------------------------------
!%--------------------------------------------------
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
subroutine verify2_reaction &
  (this, &
   reaction, &
   nerror, &
   errormsg)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Verify different attributes encapsulated in the reaction object. 
!
!   $Arguments:
!
 
type(t_reaction)                  :: this

type(t_reaction), intent(in)      :: reaction

integer                           :: nerror

character(len=300), pointer       :: errormsg(:) 
 
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

if (reaction%numsp>0) then
 allocate (pspecies(reaction%numsp))
end if
 
do i=1,reaction%numsp
 pspecies (i)= reaction%pspecies(i)%ptr
end do
 
 
call verify_reaction &
  (this, &
   name=reaction%name, &
   numsp=reaction%numsp, &
   species=pspecies, &
   logk=reaction%lnk, &
   stq=reaction%stqaqpri, &
   numcoefflogk=reaction%numcoefflogk, &
   ithsecsp=reaction%ithsecsp, &
   nerror=nerror, &
   errormsg=errormsg)
 
if(associated(pspecies)) nullify(pspecies)
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
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
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
SUBROUTINE LEAST_SQUARES_FIT &
 (MXBASIS,NBASIS ,NTEMP  ,NDNRMT ,IFAIL  ,IPRINT ,IOUTPUT &
 ,FNCT   ,ATBSFN ,ANRMT  ,CFBSFN ,INDX,MSG,ISERROR)
 
!%************************************************************************
!%
!%PURPOSE
!%   Least squares fit
!%
!%DESCRIPTION
!%   This routine computes coefficients of asn interpolation function give
!%   matrix ATBSFN of basis functions evaluated at a set of NTEMP data poi
!%   It is designed for frequent calls with the same matrix ATBSFN, so tha
!%   normal matrix is saved isbetween calls. The routine may decide to recom
!%   the normal matrix if data are missing (see definition of FNCT)
!%
!%   It was originally designed for fitting the polynomial coefs. of logK
!%   dependence on temeprature.
!%   In this case, ATBSFN contains the values of these polynomials
!%   at the temp. data points; FNCT contains the values of log-K at the te
!%   points; and CFBSFN returns the coefficients of the functions of log K
!%   dependence on temperature.
!%
!%THEORY
!%   The interpolation function is given by:
!%           f(T)=sum_k(A_k*fk(T))
!%
!%   where
!%           fk  is the k-th basis function (see READ_TEMP_CHEM)
!%           A_k is the k-th coefficient
!%
!%   The A_k's(CFBSFN) are computed as the solution of the "normal" equati
!%
!%             AN(j,k)*A_k=B_j                      (Eq. 1)
!%
!%   where - B_j is the right hand side
!%
!%             B_j=sum_i(f(T_i)*AT(j,i))            (Eq. 2)
!%
!%         - AT(j,i) is the j-th basis function evaluated at data point T_
!%         - f(T_i) interpolated function evaluated at data points (FNCT)
!%         - AN(j,k) is the "normal" matrix                        (ANRMT)
!%
!%             AN(j,k)=sum_i(AT(j,i)*AT(k,i))       (Eq. 3)
!%                                                                 (ATBSFN
!%
!%ARGUMENTS : SCALARS
!%   MXBASIS       Max numisber of basis functions as specified in calling r
!%   NTEMP         Actual numisber of data points
!%   NBASIS        Actual numisber of "basis" functions
!%   NDNRMT        On input, 0 if normal matrix needs to isbe calculated.
!%                 On output, 0 if normal matrix will have to isbe recalcula
!%                 next time, and 1 otherwise. It should isbe initialized to
!%                 isbefore the first call to this routine and not touched
!%                 afterwards
!%   IFAIL:        On output, numisber of missing data points in FNCT
!%   IPRINT        Controls printout on messagges
!%   IOUTPUT       Output file unit
!%
!%ARGUMENTS : ARRAYS
!%   FNCT         Values of interpolated function at the NTEMP data points
!%                a component is greater than 499.0, them the value is unk
!%                (ie, missing) and the normal matrix will isbe recomputed.
!%                this case, NDNRMT is set to 1, so that the normal matrix
!%                isbe recomputed in the next call to this routine, and IFAI
!%                set to the numisber of "missing" data points
!%   ATBSFN       Matrix of NBASIS basis functions evaluated at the data p
!%   ANRMT        LU decomposition of normal matrix
!%   CFBSFN       Coefficients of interpolated function
!%   INDX         Workspace for matrix operations
!%
!%INTERNAL VARIABLES: SCALARS
!%   FI           Temporary storage of FNCT(I)
!%
!%INTERNAL VARIABLES: SCALARS
!%   INDICE(I)    zero if the Ith data point is missing
!%
!%HISTORY
!%   Created by J.Carrera and C.Ayora (Jan,1998)
!%************************************************************************
 
IMPLICIT REAL*8 (A-H,O-Z)
 
PARAMETER (MXTMP=20)
 
DIMENSION FNCT(NTEMP), ATBSFN(MXBASIS,NTEMP), &
          ANRMT(MXBASIS,MXBASIS), CFBSFN(MXBASIS), INDX(MXBASIS), &
          INDICE(MXTMP)
CHARACTER(LEN=*) MSG
LOGICAL ISERROR

MSG=''
ISERROR=.FALSE. 
 
!%---------------------------------------- Initial checks and debugging op
 
IF (NBASIS.GT.NTEMP) THEN
   WRITE (IOUTPUT,8001) NBASIS,NTEMP
8001    FORMAT(' ERROR from subroutine LEAST_SQUARES_FIT:', &
   ' cannot fit a',i3,' terms polynomial(NBASIS) to',i3, &
   'data (NTEMP)'/' either increase or reduce NBASIS')
   STOP
END IF
 
IF (NBASIS.GT.MXBASIS) THEN
   WRITE (IOUTPUT,8002) NBASIS,MXBASIS
8002    FORMAT(' Dimensions ERROR from subroutine LEAST_SQUARES_FIT:', &
   '  NBASIS=',I3,'  MXBASIS=',I3)
   STOP
END IF
 
!%-------------------------------- Computes right hand side of normal equa
 
 DO  J=1, NBASIS
     CFBSFN(J)= 0D0
 END DO
 IFAIL=0
 DO I=1, NTEMP
    IF(FNCT(I).GT.499.0) THEN  ! Data missing for I-th temp. data po
      IFAIL=IFAIL+1           ! Normal matrix needs to isbe recalculat
      INDICE(I)=0
    ELSE
      INDICE(I)=1
      FI=FNCT(I)
      DO J=1, NBASIS
         CFBSFN(J)=CFBSFN(J)+FI*ATBSFN(J,I)                ! (Eq. 2)
      END DO
    END IF
 END DO
 
!%-------------- If needed, computes normal equations matrix(LU decomposit
 
 IF (IFAIL.GE.1 .OR. NDNRMT.EQ.0) THEN
    DO J=1, NBASIS
       DO K=1, J
          ANRMT(K,J)=0D0
          DO I=1, NTEMP
             IF (INDICE(I).EQ.1) THEN
                ANRMT(K,J)=ANRMT(K,J)+ATBSFN(K,I)*ATBSFN(J,I)   ! (E
             END IF
          END DO
          ANRMT(J,K)=ANRMT(K,J)
       END DO
    END DO
 
!%-------------------------------------------- LU decomposition of normal
 
    CALL LUDCMP(ANRMT,NBASIS,MXBASIS,INDX,DD,MSG,ISERROR)
	IF (ISERROR) RETURN 
 END IF
 
!%----------------------- Solve (Eq. 1) to compute the CFBSFN (A_k) coeffi
 
 CALL LUBKSB(ANRMT,NBASIS,MXBASIS,INDX,CFBSFN)
 IF (IFAIL.GE.1) THEN
    NDNRMT=0         ! Next time will need to recompute normal matri
 ELSE
    NDNRMT=1         ! Next time will not need to recompute normal m
 END IF
 
 RETURN
end subroutine
!%***********************************************************************
!%***********************************************************************
!%***********************************************************************
!%***********************************************************************
!%***********************************************************************
!%***********************************************************************
end module m_reaction
