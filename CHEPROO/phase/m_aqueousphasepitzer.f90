module m_aqueousphasepitzer
!-------------------------------------------------------------------------
!
!   $Description: This class represents an aqueous phase with the 
! Pitzer model. This class compute activity coefficient and its derivatives 
! according to 
!
! S. A. Bea, J. Carrera and C. Ayora. Modeling of concentrated aqueous solutions: 
! Efficient implementation of Pitzer equations in geochemical and reactive
! transport models (under review). 
!
!
!   $Use: 
! use m_parentaqueousphase
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
!% Modules corresponding to CHEPROO package
!%-------------------------------------------------------------------------
use m_species
use m_parentaqueousphase
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%-------------------------------------------------------------
!%-------------------------------------------------------------
private                                ::
!%-------------------------------------------------------------
!% Public services
!%-------------------------------------------------------------
public:: &
create_ &                  ! Create the object
,destroy_ &                ! Destroy the object 
,set_ &                    ! Set the object from virial coefficients data base
,set_pparent_ &            ! Set pointr to parent aqueous phase
,compute_act_coeff_ &      ! Compute activity coefficients according to Pitzer equations
,compute_dact_coeff_ &     ! Compute derivatives of the activity coefficients
,compute_density_ &        ! Compute density accordig to Monnin (1994)
,update_ &                 ! Update the parameters that depends of the temperature in the object 
,write_ &                  ! Write in ascii the attributes encapsulated in the object
,assignment(=)             ! Copy an object in other object
!%-------------------------------------------------------------
!% Private services
!%-------------------------------------------------------------
private   :: &
set_virial_  &
,read_pitzer_data_base_ &
,compute_act_coeff_loc_ &
,compute_qmatrices_ &
,compute_volex_ &
,compute_dh_ &
,compute_qc_ &
,compute_q_ &
,compute_ql_ &
,compute_t_ &
,compute_fbeta_ &
,compute_fbetavol_ &
,compute_fj_ &
,assign_fbeta_ &
,assign_fbetavol_ &
,assign_fj_ &
,compute_ddh_ &
,compute_dq_ &
,compute_dqc_ &
,compute_dql_ &
,compute_dt_ &
,check_consistence_ &
,gclm_ &
,dgclm_
!%-------------------------------------------------------------
!% Constant parameters
!%-------------------------------------------------------------
real*8, parameter               :: &
!%-------------------------------------------------------------
!% Virial coefficients for KCl system (for MacInnes convention
!% scaled 
!%-------------------------------------------------------------
!% Virial coefficients for Cl- and K+ (at 25oC) used for 
!% MacInnes scaled   
!%-------------------------------------------------------------
mtb0kcl=0.04835d0,         & ! B0 virial coefficient K+ Cl-
mtb1kcl=0.2122d0,          & ! B1 virial coefficient K+ Cl-
mtc0kcl=-0.00084d0,        & ! C0 virial coefficient K+ Cl-
!%-------------------------------------------------------------
!% Coefficients of the temperature dependence of Afi 
!% (see Clegg and Whitfield (1991))
!%-------------------------------------------------------------
!% Constant Debye-H�ckel parameter
!%-------------------------------------------------------------
bdh=1.2d0,                  & 
!%-------------------------------------------------------------
!% Limiting Debye-H�ckel slope to 25�   0.39153  0.392
!%-------------------------------------------------------------
aphi25=0.392d0,             &  
!%-------------------------------------------------------------
!% Limiting Debye-H�ckel slope to 25� for density calculations 
!% dAphi/dP??
!% cm3 kg1/2/mol3/2
!% taken from Monnin (1994)
!%-------------------------------------------------------------
aphi25vol=1.8743d0,         &                               
!%-------------------------------------------------------------
!cprovi aphi25=0.39153d0,   & ! Limiting Debye-H�ckel slope to 25�   0.39153  0.392
!%-------------------------------------------------------------
c0aphi=0.13422d0,          & ! Temperature depending coefficients 
c1aphi=0.0368329d0,        & ! Temperature depending coefficients 
c2aphi=14.62718d0,         & ! Temperature depending coefficients 
c3aphi=1530.1474d0,        & ! Temperature depending coefficients 
c4aphi=80.40631d0,         & ! Temperature depending coefficients 
c5aphi=4.1725332d0,        & ! Temperature depending coefficients 
c6aphi=0.1481291d0,        & ! Temperature depending coefficients 
c7aphi=1.5188505d-5,       & ! Temperature depending coefficients 
c8aphi=1.8016317d-8,       & ! Temperature depending coefficients 
c9aphi=9.3816144d-10         ! Temperature depending coefficients 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Type definition 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
type, public::t_aqueousphasepitzer    
 
private                                ::
 
type (t_parentaqueousphase), pointer   :: pp            ! Parent aqueous phase 

real*8, pointer, dimension(:)          :: beta0         ! Beta_0 virial coefficients

real*8, pointer, dimension(:)          :: beta1         ! Beta_1 virial coefficients

real*8, pointer, dimension(:)          :: beta2         ! Beta_2 virial coefficients

real*8, pointer, dimension(:)          :: theta         ! Theta virial coefficients

real*8, pointer, dimension(:)          :: lambda        ! Lambda virial coefficients

real*8, pointer, dimension(:)          :: psi           ! Psi virial coefficients

real*8, pointer, dimension(:)          :: cpz           ! C pitzer matrix (computed according equation A8)

real*8, pointer, dimension(:)          :: zprod         ! charge products [nfunj]

real*8, pointer, dimension(:,:)        :: alpha         ! [2,nfunb]

real*8, pointer, dimension(:,:)        :: alphavol      ! [2,nfunbvol]

real*8                                 :: aphi          ! Debye-H�ckel limiting slope

integer, pointer, dimension(:,:)       :: indnzbeta     ! Non zero indices for beta virial coefficients

integer, pointer, dimension(:,:)       :: indnztheta    ! Non zero indices for theta virial coefficients

integer, pointer, dimension(:,:)       :: indnzlambda   ! Non zero indices for lambda virial coefficients

integer, pointer, dimension(:,:)       :: indnzpsi      ! Non zero indices for psi virial coefficients

integer, pointer, dimension(:,:)       :: indnzcpz      ! Non zero indices for c virial coefficients

integer                                :: nfunb         ! Number of beta functions

integer                                :: nfunbvol      ! Number of beta functions

integer                                :: nfunj         ! Number of j functions

integer                                :: nnzbeta       ! Number of non cero Beta matrix terms

integer                                :: nnztheta      ! Number of non cero theta matrix terms

integer                                :: nnzcpz        ! Number of non cero C_ca matrix terms

integer                                :: nnzlambda     ! Number of non cero lambda matrix terms

integer                                :: nnzpsipz      ! Number of non cero Psi matrix terms

integer                                :: nnzq          ! Number of non cero q matrix terms (q,q',q'',q-fi and q-fi')

integer                                :: ithcl         ! Local indice of Cl species  (usefull for macinnes convention)

integer                                :: ithk          ! Local indice of K species   (usefull for macinnes convention)

logical                                :: ismacinnes    ! ismacinnes=true, then activity coefficients will be scaled according macinnes convention
 
real*8, pointer, dimension(:)          :: beta0vol      ! Derivative of beta0 virial coefficient with respect to pressure

real*8, pointer, dimension(:)          :: beta1vol      ! Derivative of beta1 virial coefficient with respect to pressure

real*8, pointer, dimension(:)          :: beta2vol      ! Derivative of beta2 virial coefficient with respect to pressure

real*8, pointer, dimension(:)          :: cpzvol        ! Derivative of cphi virial coefficient with respect to pressure

real*8, pointer, dimension(:)          :: vol0          ! Standard partial volume of salt

integer, pointer, dimension(:,:)       :: indnzbetavol  ! Non zero indices for c virial coefficients

integer, pointer, dimension(:,:)       :: indnzcpzvol   ! Non zero indices for c virial coefficients

integer, pointer, dimension(:,:)       :: indnzvol0     ! Non zero indices for vol0 

integer                                :: nnzbetavol    ! Number of non cero derivaties of Beta virial coefficients with respect to pressure 

integer                                :: nnzcpzvol     ! Number of non cero derivaties of Cij virial coefficients with respect to pressure

integer                                :: nnzvol0       ! Number of non cero standard partial volume of salt 

end type t_aqueousphasepitzer
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!% Interfaces 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
 
module procedure create_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
 
module procedure destroy_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
 
module procedure set_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_pparent_
 
module procedure set_pparent_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_act_coeff_
 
module procedure compute_act_coeff_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dact_coeff_
 
module procedure compute_dact_coeff_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_density_
 
module procedure compute_density_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_
 
module procedure update_temp_param_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
 
module procedure write_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment (=)
 
module procedure copy_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%-----------------------------Private services--------------------------- 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_virial_
 
module procedure set_virial_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface read_pitzer_data_base_
 
module procedure read_pitzer_data_base_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_act_coeff_loc_
 
module procedure compute_act_coeff_loc_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_volex_
 
module procedure compute_volex_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_qmatrices_
 
module procedure compute_qmatrices_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dh_
 
module procedure compute_dh_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_qc_
 
module procedure compute_qc_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_q_
 
module procedure compute_q_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_ql_
 
module procedure compute_ql_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_t_
 
module procedure compute_t_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_fbeta_
 
module procedure compute_fbeta_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_fbetavol_
 
module procedure compute_fbetavol_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_fj_
 
module procedure compute_fj_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assign_fbeta_
 
module procedure assign_fbeta_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assign_fbetavol_
 
module procedure assign_fbetavol_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assign_fj_
 
module procedure assign_fj_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_ddh_
 
module procedure compute_ddh_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dq_
 
module procedure compute_dq_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dqc_
 
module procedure compute_dqc_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dql_
 
module procedure compute_dql_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dt_
 
module procedure compute_dt_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface check_consistence_
 
module procedure check_consistence_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface gclm_
 
module procedure gclm_aqphpitzer
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface dgclm_
 
module procedure dgclm_aqphpitzer
 
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
subroutine create_aqphpitzer &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create the object
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(inout) :: this  ! Type aqueous phase Pitzer 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical                :: &
 iserror 
real*8, parameter      :: &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
!% Nullify and zeroing the attributes
!%------------------------------------------------------------
this%pp => null ()
!%------------------------------------------------------------
this%beta0 => null ()
this%beta1 => null ()
this%beta2 => null ()
this%theta => null ()
this%lambda => null ()
this%psi => null ()
this%cpz => null ()
!%------------------------------------------------------------
!% Nullify derivatives of the virial coefficients with 
!% respect to the temperature 
!%------------------------------------------------------------
this%beta0vol => null ()
this%beta1vol => null ()
this%beta2vol => null ()
this%cpzvol => null ()
this%vol0 => null ()
this%indnzbetavol => null ()
this%indnzcpzvol => null ()
this%indnzvol0 => null ()
!%------------------------------------------------------------
!%
!%------------------------------------------------------------
this%zprod => null ()
this%alpha => null ()
this%alphavol => null ()
this%indnzbeta => null ()
this%indnztheta => null ()
this%indnzlambda => null ()
this%indnzpsi => null ()
this%indnzcpz => null ()
!%------------------------------------------------------------
!% Initialice variables 
!%------------------------------------------------------------
this%nfunb=0
this%nfunbvol=0
this%nfunj=0
this%nnzbeta=0
this%nnzbetavol=0
this%nnztheta=0
this%nnzlambda=0
this%nnzq=0
this%nnzcpz=0
this%nnzcpzvol=0
this%nnzvol0=0
this%nnzpsipz=0
this%ithcl=0
this%ithk=0
this%ismacinnes=.false.
this%aphi=r0
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_aqphpitzer &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy the object 
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(inout) :: this   ! Type aqueous phase Pitzer 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, parameter     :: &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
!% Nullify pointer to parent aqueous phase object 
!%------------------------------------------------------------
this%pp => null ()
!%------------------------------------------------------------
!% Deallocate pointers 
!%------------------------------------------------------------
call check_pointer_ (this%beta0,1,.false.)
call check_pointer_ (this%beta1,1,.false.)
call check_pointer_ (this%beta2,1,.false.)
call check_pointer_ (this%theta,1,.false.)
call check_pointer_ (this%lambda,1,.false.)
call check_pointer_ (this%psi,1,.false.)
call check_pointer_ (this%cpz,1,.false.)
!%-----------------------------------------------------------
!%
!%-----------------------------------------------------------
call check_pointer_ (this%beta0vol,1,.false.)
call check_pointer_ (this%beta1vol,1,.false.)
call check_pointer_ (this%beta2vol,1,.false.)
call check_pointer_ (this%cpzvol,1,.false.)
call check_pointer_ (this%vol0,1,.false.)
call check_pointer_ (this%alphavol,1,1,.false.)
call check_pointer_ (this%indnzbetavol,1,1,.false.)
call check_pointer_ (this%indnzcpzvol,1,1,.false.)
call check_pointer_ (this%indnzvol0,1,1,.false.)
!%-----------------------------------------------------------
!%
!%-----------------------------------------------------------
call check_pointer_ (this%zprod,1,.false.)
call check_pointer_ (this%alpha,1,1,.false.)
call check_pointer_ (this%indnzbeta,1,1,.false.)
call check_pointer_ (this%indnztheta,1,1,.false.)
call check_pointer_ (this%indnzlambda,1,1,.false.)
call check_pointer_ (this%indnzpsi,1,1,.false.)
call check_pointer_ (this%indnzcpz,1,1,.false.)
!%-----------------------------------------------------------
!% Initialice variables 
!%-----------------------------------------------------------
this%nfunb=0
this%nfunbvol=0
this%nfunj=0
this%nnzbeta=0
this%nnzbetavol=0
this%nnztheta=0
this%nnzlambda=0
this%nnzq=0
this%nnzcpz=0
this%nnzcpzvol=0
this%nnzvol0=0
this%nnzpsipz=0
this%ithcl=0
this%ithk=0
this%ismacinnes=.false.
this%aphi=r0
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_aqphpitzer &
   (this, &
    pp, &
    namedatabase, &
    nameconvention, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the object from virial coefficients data base.
! 
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(inout)       :: this            ! Type aqueous phase Pitzer 
 
type(t_parentaqueousphase), intent(in), target  :: pp              ! Type parent aqueous phase 

character(len=*), intent(in)                    :: namedatabase    ! Name of the virial coefficients data base

character(len=*), intent(in)                    :: nameconvention  ! Name of the convention scaled  

logical, intent(out)                            :: iserror         ! If iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: The namedatabase must contain the directory path
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                  :: &
 namevirial
integer        :: &
 i, &
 j, &
 nsp
integer::  &
 iostat
type(xml_t)::  &
 fxml
character(len=100) :: &
 name, &
 msg, &
 empty, &
 nameprop 
real*8             :: &
 value1, &
 value2
type(dictionary_t) :: &
 attributes 
real*8, pointer    :: &
 alpha(:,:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
!%-------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------
empty=' '
!%-------------------------------------------------------
this%pp => pp
!%-------------------------------------------------------
!% Check the convention for scaled 
!%-------------------------------------------------------
select case (nameconvention)
 
case ('macinnes','Macinnes','MACINNES')
 this%ismacinnes=.true.
case default
  if (nameconvention/=empty) then
   msg='Warning, not recognized convention:'
   call add_ (msg,nameconvention)
   print *,'****************************'
   print *,'Phase:'
   print *,'Name:', this%pp%pp%name
   print *,'service: set_'
   print *, msg
   print *,'****************************'
  end if
  this%ismacinnes=.false.
end select
!%------------------------------------------------------
!% If the MacInnes convention is implemented, then find
!% the local indice of Cl- and K+ 
!%------------------------------------------------------
if (this%ismacinnes) then
  do i=1,this%pp%pp%numsp
   select case (this%pp%pp%pspecies(i)%ptr%name)
   case ('cl-','CL-','Cl-','cl-1','Cl-1','CL-1')
    this%ithcl=i
   case ('k+','K+','k+1','K+1')
    this%ithk=i
   end select
  end do
  if (this%ithcl==0.or.this%ithk==0) then
    msg='Warning: MacInnes convention not implemented (Cl and K was not defined as components)'
    print *,'****************************'
    print *,'Phase:'
    print *,'Name:', this%pp%pp%name
    print *,'service: set_'
    print *, msg
    print *,'****************************'
    this%ismacinnes=.false.
  end if
end if
!%---------------------------------------------------
!% Open the xml file
!%---------------------------------------------------
call open_xmlfile (namedatabase,fxml,iostat)
if (iostat/=0) then
 msg='Not found Pitzer data base: '
 call add_ (msg,namedatabase)
 goto 10
end if
!%---------------------------------------------------
!% Reading Pitzer data base 
!%---------------------------------------------------
print *,'=======> Reading Pitzer virial coefficients data base (starting)'
!%---------------------------------------------------
! Reading Aphi 
!%---------------------------------------------------
namevirial='//aphi'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading B0 virial terms
!%---------------------------------------------------
namevirial='//b0'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror,alpha=alpha)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading B1 virial terms
!%---------------------------------------------------
namevirial='//b1'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading B2 virial terms
!%---------------------------------------------------
namevirial='//b2'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading Cfi virial terms
!%---------------------------------------------------
namevirial='//cfi'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading Theta virial terms
!%---------------------------------------------------
namevirial='//theta'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading Psi virial terms
!%---------------------------------------------------
namevirial='//psi'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading Lambda virial terms
!%---------------------------------------------------
namevirial='//lambda'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
print *,'=======> Reading Pitzer virial coefficients data base (finishing)'
!%---------------------------------------------------
!% Now read the virial coefficients for the prediction
!% of fluid properties based on Pitzer equations 
!% Read Monnin (1994)
!%---------------------------------------------------
print *,'=======> Reading Pitzer virial coefficients for fluid properties prediction (starting)'
!%---------------------------------------------------
! Reading dB0/dP virial terms
!%---------------------------------------------------
namevirial='//b0vol'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading dB1/dP virial terms
!%---------------------------------------------------
namevirial='//b1vol'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading dB2/dP virial terms
!%---------------------------------------------------
namevirial='//b2vol'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%---------------------------------------------------
! Reading dCfi/dP virial terms
!%---------------------------------------------------
namevirial='//cfivol'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10
!%---------------------------------------------------
! Reading vol0 terms
!%---------------------------------------------------
namevirial='//vol0'
call read_pitzer_data_base_ (this,fxml,namevirial,iserror)
if (iserror) goto 10 
!%------------------------------------------------------------
print *,'=======> Reading Pitzer virial coefficients for fluid properties prediction (finishing)'
!%------------------------------------------------------------
!% Check the consistence between geochemical system and 
!% virial coefficients data base. 
!%------------------------------------------------------------
call check_consistence_ (this,6,iserror)
if (iserror) goto 20 
!%------------------------------------------------------------
!%------------------------------------------------------------
!% Assign beta functions 
!%------------------------------------------------------------
call assign_fbeta_ (this,alpha(1:2,1:this%nnzbeta),iserror)
if (iserror) goto 20
!%------------------------------------------------------------
!% Assign betavol functions 
!% for density calculations 
!%------------------------------------------------------------
call assign_fbetavol_ (this,iserror)
if (iserror) goto 20
!%------------------------------------------------------------
!% Assign J functions 
!%------------------------------------------------------------
call assign_fj_ (this)
!%------------------------------------------------------------
this%nnzq = this%nnzbeta + this%nnztheta + this%nnzlambda
!%------------------------------------------------------------
!% Update the parent phase according to reference temperature 
!%------------------------------------------------------------
call update_ (this,this%pp%pp%tempref,iserror)
!%-----------------------------------------------------------
20 continue
!%-----------------------------------------------------------
! End and close xml file
!%-----------------------------------------------------------
call endfile_xmlfile (fxml)
call close_xmlfile(fxml)
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (alpha,1,1,.false.)
!%------------------------------------------------------------
if (iserror) goto 10
!%------------------------------------------------------------
!% Read molar volumes from data base file 
!%------------------------------------------------------------
open(unit=1,file='aqueousphase.dat',status='old',err=30)
!%------------------------------------------------------------
nameprop='molvol'
value1=r0
value2=r0
do i=1,this%pp%pp%numsp  
   nameprop='molvol'
   call add_ (this%pp%pp%pspecies(i)%ptr,value1,nameprop,'m3/mol',iserror)
   if (iserror) goto 30
   nameprop='molweight'
   call add_ (this%pp%pp%pspecies(i)%ptr,value1,nameprop,'m3/mol',iserror)
   if (iserror) goto 30   
end do  
read(1,*) nsp
do i=1,nsp
  
  read(1,*) name,value1,value2
    
  do j=1,this%pp%pp%numsp 
    if (this%pp%pp%pspecies(j)%ptr%name==name) then
       nameprop='molweight'
       call add_ (this%pp%pp%pspecies(j)%ptr,value1,nameprop,'gr/mol',iserror)
       if (iserror) goto 30
       nameprop='molvol'
       call add_ (this%pp%pp%pspecies(j)%ptr,value2,nameprop,'m3/mol',iserror)
       if (iserror) goto 30
       exit 
    end if 
  end do  
   
end do
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
40 continue 
if (iserror) goto 10 
!%------------------------------------------------------------
return
 
10 continue 
print *,'****************************'
print *,'Phase:'
print *,'Name:', this%pp%pp%name
print *,'Service: set_'
print *, msg
print *,'****************************'
iserror=.true.
return

30 goto 40 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_pparent_aqphpitzer &
   (this, &
    pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set pointr to parent aqueous phase
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(inout)             :: this    ! Type aqueous phase Pitzer

type(t_parentaqueousphase), intent(in), target        :: pp      ! Type parent aqueous phase
 
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
subroutine compute_act_coeff_aqphpitzer &
 (this, &
  g, &
  c, &
  iserror, &
  ionstr)
 
implicit none 
!-------------------------------------------------------------------------
!
!   $Description: Compute activity coefficients according to Pitzer equations.
! The implemented algorithm was based on Bea et al. (2008) (eq. 1, 2 and 3)
!  
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(in)          :: this      ! Type aqueous phase Pitzer

real*8, pointer, dimension(:)                    :: g         ! Activity coefficients vector [numsp]

real*8, intent(in), dimension(this%pp%pp%numsp)  :: c         ! Molality of species [numsp] 

logical, intent(out)                             :: iserror   ! iserror=true, then there was an error

real*8, intent(out), optional                    :: ionstr    ! Ionic strength [M]
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: 
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer                           :: &
 qphi(:) => null (), &
 q(:) => null (), &
 qpri(:) => null (), &
 dqpri(:) => null (), &
 dqphi(:)  => null ()
integer, pointer                          :: &
 indnzq(:,:) => null ()
real*8                                    :: &
 ionstrloc, &
 ionstrz, &
 m, &
 gclm, &
 gcl, &
 charge
logical                                   :: &
 haveionstr
character(len=100)                        :: &
 msg 
integer                                   :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
! 

 

!%--------------------------------------------------------------
msg=''
iserror=.false.
!%--------------------------------------------------------------
!% Check optional arguments
!%--------------------------------------------------------------
haveionstr=present(ionstr)
!%--------------------------------------------------------------
!% Activity coefficients are computed for this private 
!% subroutine
!%--------------------------------------------------------------
call compute_act_coeff_loc_ &
 (this, &
  g, &
  c, &
  ionstrloc, &
  ionstrz, &
  m, &
  qphi, &
  q, &
  qpri, &
  dqpri, &
  dqphi, &
  indnzq, &
  gclm, &
  iserror)
if (iserror) goto 20  
!%-------------------------------------------------------------
!% If ismacinnes=true, then activity coeffients are scaled 
!% according to MacInnes convention (MacInnes, 1919) 
!%-------------------------------------------------------------
if (this%ismacinnes) then
 gcl = g (this%ithcl)
 do i=1,this%pp%pp%numsp
  call get_prop_ (this%pp%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
  if (iserror) goto 10
  if (i/=this%pp%ithw) then
    g(i) = g(i) * (gcl / gclm)**charge
  end if
 end do
end if
!%--------------------------------------------------------------
!% If the ionstr is present as argument, then assign the value 
!%--------------------------------------------------------------
if (haveionstr) ionstr=ionstrloc
!%--------------------------------------------------------------
20 continue 
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (q,1,.false.)
call check_pointer_ (qphi,1,.false.)
call check_pointer_ (qpri,1,.false.)
call check_pointer_ (dqpri,1,.false.)
call check_pointer_ (dqphi,1,.false.)
call check_pointer_ (indnzq,1,1,.false.)
if (iserror) goto 10 
!%--------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%pp%name
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
!%************************************************************
subroutine compute_dact_coeff_aqphpitzer &
   (this, &
    dg, &
    c, &
    dc, &
	iserror, &
    dionstr, &
	dtemp, &
    g)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute derivatives of the activity coefficients.
! The implemented algorithm was based on Bea et al., (2008) (eq. 19 and 20)
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(in)                    :: this      ! Type aqueous phase Pitzer

real*8, pointer, dimension(:,:)                            :: dg        ! Derivatives of the activity coefficients  

real*8, intent(in), dimension(:,:)                         :: dc        ! Derivatives of the concentration of species

real*8, intent(in), dimension(this%pp%pp%numsp)            :: c         ! Molality vector

real*8, intent(in), optional, dimension(this%pp%pp%numsp)  :: g         ! Activity coefficients vector 

logical, intent(out)                                       :: iserror   ! iserror=true, then there was an error

real*8, pointer, optional, dimension(:)                    :: dionstr   ! Derivatives of the ionic strength

real*8, intent(in), optional, dimension(:)                 :: dtemp     ! Derivatives of the temperature
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
 ndimder, &
 i
real*8, pointer                     :: &
 dosco(:) => null(), &
 dionstrloc(:) => null(), &
 dm(:) => null(), &
 gloc(:) => null()
real*8, pointer                     :: &
 qphi(:) => null(), &
 q(:) => null(), &
 qpri(:) => null(), &
 dqpri(:) => null(), &
 dqphi(:) => null(), &
 dgcl(:) => null(), &
 dgclm(:) => null()
integer, pointer :: &
 indnzq(:,:) => null()
real*8           :: &
 ionstr, &
 ionstrz, &
 m, &
 gclm, &
 cte, &
 gcl, &
 charge, &
 gm, &
 osco
logical          :: &
 havedionstr, &
 haveg, &
 havedtemp 
character(len=100)   :: &
 msg
type(t_species), pointer :: &
 pspecies => null () 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Check optional arguments
!%--------------------------------------------------------------
havedionstr=present(dionstr)
haveg=present(g)
havedtemp=present(dtemp) 
!%--------------------------------------------------------------
!% Allocate dg
!%--------------------------------------------------------------
ndimder=size(dc,2)
call check_pointer_ (dg,this%pp%pp%numsp,ndimder,.true.)
!%--------------------------------------------------------------
!% Allocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (dosco,ndimder,.true.)
!%--------------------------------------------------------------
!% Computes the ionic strength (eq. 4)
!%--------------------------------------------------------------
call compute_ionstr_ (this%pp,ionstr,c,iserror)
if (iserror) goto 20  
!%--------------------------------------------------------------
!% Computes derivatives of the ionic strength (eq. 21)
!%--------------------------------------------------------------
call compute_dionstr_ (this%pp,dionstrloc,dc,iserror)
if (iserror) goto 20  
!%--------------------------------------------------------------
!% Computes dM (eq. 23)
!%--------------------------------------------------------------
call compute_sum_dc_ (this%pp,dm,dc)
!%--------------------------------------------------------------
if (haveg) then
   call check_pointer_ (gloc,this%pp%pp%numsp,.true.)
   if (this%ismacinnes) then 
     call compute_dh_ (this,osco,gloc,ionstr,gclm,iserror)
     if (iserror) goto 20
   end if 
   gloc=g
   call compute_qmatrices_ &
    (this, &
     qphi, &
     q, &
     qpri, &
     dqpri, &
     dqphi, &
     indnzq, &
     ionstr, &
     msg, &
     iserror)
   if(iserror) goto 20   
else
   call compute_act_coeff_loc_ &
    (this, &
     gloc, &
     c, &
     ionstr, &
     ionstrz, &
     m, &
     qphi, &
     q, &
     qpri, &
     dqpri, &
     dqphi, &
     indnzq, &
     gclm, &
     iserror)
   if (iserror) goto 20  
end if
!%--------------------------------------------------------------
!% Computes and adds eq. B1 and B2
!%--------------------------------------------------------------
call compute_ddh_ &
 (this, &
  dg, &
  dosco, &
  ionstr, &
  dionstrloc, &
  ndimder, &
  dgclm, &
  iserror)
if (iserror) goto 20
!%--------------------------------------------------------------
!% Computes and adds eq. 25 and 26
!%--------------------------------------------------------------
call compute_dq_ &
  (this, &
   dg, &
   dosco, &
   q, &
   qpri, &
   qphi, &
   dqpri, &
   dqphi, &
   indnzq, &
   c, &
   ionstr, &
   dionstrloc, &
   dc, &
   ndimder, &
   msg, &
   iserror)
!%--------------------------------------------------------------
!% Computes and adds eq. 27
!%--------------------------------------------------------------
call compute_dqc_ &
  (this, &
   dg, &
   dosco, &
   c, &
   ionstrz, &
   dc, &
   ndimder, &
   msg, &
   iserror)
 if (iserror) goto 20
!%--------------------------------------------------------------
!% Computes and adds eq. 28 
!%--------------------------------------------------------------
call compute_dql_ &
  (this, &
   dosco, &
   c, &
   dc, &
   ndimder, &
   msg, &
   iserror)
 if (iserror) goto 20
!%--------------------------------------------------------------
!% Computes and adds eq. 29
!%--------------------------------------------------------------
call compute_dt_ &
  (this, &
   dg, &
   dosco, &
   c, &
   gloc, &
   dc, &
   dm, &
   ndimder)
!%--------------------------------------------------------------
!% Allocate optional pointer 
!%--------------------------------------------------------------
if (havedionstr) then
 call check_pointer_ (dionstr,ndimder,.true.)
 dionstr=dionstrloc
end if
!%--------------------------------------------------------------
!% Scales the derivatives according to MacInnes 
!% convention (MacInnes, 1919)
!%--------------------------------------------------------------
if (this%ismacinnes) then
  call check_pointer_ (dgcl,ndimder,.true.)
  gcl = gloc (this%ithcl)
  dgcl = dg(this%ithcl,:)
  cte=gcl/gclm
  do i=1,this%pp%pp%numsp
    pspecies => this%pp%pp%pspecies(i)%ptr
    call get_prop_ (pspecies,charge,'charge',msg,iserror)
    if (iserror) goto 20
    if (i/=this%pp%ithw) then
        gm=gloc(i)*cte**charge
        dg(i,:)=dg(i,:)/gloc(i)
        dg(i,:)=dg(i,:)+charge*dgcl/gcl
        dg(i,:)=dg(i,:)+charge*dgclm/gclm
        dg(i,:)=dg(i,:)*gm
   end if
  end do
end if
!%--------------------------------------------------------------
20 continue 
!%--------------------------------------------------------------
!% Deallocate and nullify local pointers 
!%--------------------------------------------------------------
call check_pointer_ (gloc,1,.false.)
call check_pointer_ (q,1,.false.)
call check_pointer_ (qphi,1,.false.)
call check_pointer_ (dqpri,1,.false.)
call check_pointer_ (qpri,1,.false.)
call check_pointer_ (dqphi,1,.false.)
call check_pointer_ (dosco,1,.false.)
call check_pointer_ (dionstrloc,1,.false.)
call check_pointer_ (dm,1,.false.)
call check_pointer_ (indnzq,1,1,.false.)
if (this%ismacinnes) then
  call check_pointer_ (dgcl,1,.false.)
  call check_pointer_ (dgclm,1,.false.)
end if
pspecies => null ()
if (iserror) goto 10
!%--------------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%pp%name
print *,'Service: compute_dact_coeff_'
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
subroutine write_aqphpitzer &
   (this, &
    ioutput, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii the attributes encapsulated in the object.
! Also write the virial coefficients data base.   
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(in), target   :: this     ! Type aqueous phase Pitzer

integer, intent(in)                               :: ioutput  ! Output unit 

logical, intent(out)                              :: iserror  ! iserror=true, then there was an error
 
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
integer                 :: &
 i, &
 nvirial
integer, pointer        :: &
 isp1 => null (), &
 isp2 => null (), &
 isp3 => null (), &
 ifunb => null ()
real*8, pointer         :: &
 coeff1 => null (), &
 coeff2 => null (), &
 coeff3 => null (), &
 alpha1 => null (), &
 alpha2 => null ()
character(len=100), pointer :: &
 namesp1 => null (), &
 namesp2 => null (), &
 namesp3 => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
nvirial=0 
!%-----------------------------------------------------------
!% Write head 
!%-----------------------------------------------------------
write(ioutput,*)  '-------------------------------------------------------------------'
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) '-------------------------------------------------------------------'
write(ioutput,*)  'Aqueous Phase Object                                               '
write(ioutput,*)  'HMW model (Harvie et al., 1984) implemented according to           '
write(ioutput,*)  'Bea et al. (2009)                                                  '
write(ioutput,*)  '-------------------------------------------------------------------'
!%-----------------------------------------------------------
if (this%ismacinnes) then
 write (ioutput,*) 'MacInnes convention to scale activity coefficients is used'
else
 write (ioutput,*) 'MacInnes convention to scale activity coefficients is not used'
end if
!%-----------------------------------------------------------
!% Write the parent aqueous phase 
!%-----------------------------------------------------------
call write_ (this%pp,ioutput,iserror)
if (iserror) goto 20 
!%-----------------------------------------------------------
!% Write parameters concerning to Debye-H�ckel term 
!%-----------------------------------------------------------
write(ioutput,*)  '-------------------------------------------------------------------'
write (ioutput,9) 'Debye-Huckel limiting slope (aphi):',this%aphi
write(ioutput,*)  '-------------------------------------------------------------------'
!%-----------------------------------------------------------
!% Write B0, B1 and B2 virial coefficients
!%-----------------------------------------------------------
if (this%nnzbeta>0) then
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) 'Beta virial coefficients                                           '
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,1) 'sps1','sps2','B0','B1','B2','Alpha1','Alpha2'
write (ioutput,*) '-------------------------------------------------------------------'
do i=1,this%nnzbeta
 isp1 => this%indnzbeta(1,i)
 isp2 => this%indnzbeta(2,i)
 ifunb => this%indnzbeta(3,i)
 namesp1 => this%pp%pp%pspecies(isp1)%ptr%name
 namesp2 => this%pp%pp%pspecies(isp2)%ptr%name
 coeff1 => this%beta0(i)
 coeff2 => this%beta1(i)
 coeff3 => this%beta2(i)
 alpha1 => this%alpha(1,ifunb)
 alpha2 => this%alpha(2,ifunb)
 write (ioutput,2) namesp1,namesp2,coeff1,coeff2,coeff3,alpha1,alpha2
 if (coeff1/=0.0d0) nvirial=nvirial+1
 if (coeff2/=0.0d0) nvirial=nvirial+1
 if (coeff3/=0.0d0) nvirial=nvirial+1
end do
write (ioutput,*) '-------------------------------------------------------------------'
end if
!%-----------------------------------------------------------
!% Write C terms
!%-----------------------------------------------------------
if (this%nnzcpz>0) then
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) 'C terms                                                            '
write (ioutput,*) 'Cij =       Cij0  <======    these coefficients are read from      '     
write (ioutput,*) '      -------------------    data base                             '
write (ioutput,*) '      2 aqrt( (abs(zizj) )                                         '
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,7) 'sps1','sps2','C'
write (ioutput,*) '-------------------------------------------------------------------'
do i=1,this%nnzcpz
 isp1 => this%indnzcpz(1,i)
 isp2 => this%indnzcpz(2,i)
 namesp1 => this%pp%pp%pspecies(isp1)%ptr%name
 namesp2 => this%pp%pp%pspecies(isp2)%ptr%name
 coeff1 => this%cpz(i)
 write (ioutput,4) namesp1,namesp2,coeff1
end do
write (ioutput,*) '-------------------------------------------------------------------'
nvirial=nvirial+this%nnzcpz 
end if
!%-----------------------------------------------------------
!% Write THETA virial coefficients
!%-----------------------------------------------------------
if (this%nnztheta>0) then
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) 'Theta virial coefficients                                          '
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,8) 'sps1','sps2','theta'
write (ioutput,*) '-------------------------------------------------------------------'
do i=1,this%nnztheta
 isp1 => this%indnztheta(1,i)
 isp2 => this%indnztheta(2,i)
 namesp1 => this%pp%pp%pspecies(isp1)%ptr%name
 namesp2 => this%pp%pp%pspecies(isp2)%ptr%name
 coeff1 => this%theta(i)
 write (ioutput,4) namesp1,namesp2,coeff1
end do
write (ioutput,*) '-------------------------------------------------------------------'
nvirial=nvirial+this%nnztheta
end if
!%-----------------------------------------------------------
!% Write LAMBDA virial coefficients
!%-----------------------------------------------------------
if (this%nnzlambda>0) then
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) 'Lambda virial coefficients                                         '
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,3) 'sps1','sps2','lambda'
write (ioutput,*) '-------------------------------------------------------------------'
do i=1,this%nnzlambda
 isp1 => this%indnzlambda(1,i)
 isp2 => this%indnzlambda(2,i)
 namesp1 => this%pp%pp%pspecies(isp1)%ptr%name
 namesp2 => this%pp%pp%pspecies(isp2)%ptr%name
 coeff1 => this%lambda(i)
 write (ioutput,4) namesp1,namesp2,coeff1
end do
write (ioutput,*) '-------------------------------------------------------------------'
nvirial=nvirial+this%nnzlambda
end if
!%-----------------------------------------------------------
!% Write PSI virial coefficients
!%-----------------------------------------------------------
if (this%nnzpsipz>0) then
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) 'Psi virial coefficients                                            '
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,5) 'sps1','sps2','sps3','psi'
write (ioutput,*) '-------------------------------------------------------------------'
do i=1,this%nnzpsipz
 isp1 => this%indnzpsi(1,i)
 isp2 => this%indnzpsi(2,i)
 isp3 => this%indnzpsi(3,i)
 namesp1 => this%pp%pp%pspecies(isp1)%ptr%name
 namesp2 => this%pp%pp%pspecies(isp2)%ptr%name
 namesp3 => this%pp%pp%pspecies(isp3)%ptr%name
 coeff1 => this%psi(i)
 write (ioutput,6) namesp1,namesp2,namesp3,coeff1
end do
write (ioutput,*) '-------------------------------------------------------------------'
nvirial=nvirial+this%nnzpsipz
end if
!%-----------------------------------------------------------
write (ioutput,11) 'Total number of virial coefficients:',nvirial 
!%-----------------------------------------------------------
!% Write the standard partial volume of salt  
!%-----------------------------------------------------------
if (this%nnzvol0>0) then
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) ' Standard partial volume of salt [cm3 / mol]                       '
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,3) 'sps1','sps2','vol0'
write (ioutput,*) '-------------------------------------------------------------------'
do i=1,this%nnzvol0
 isp1 => this%indnzvol0(1,i)
 isp2 => this%indnzvol0(2,i)
 namesp1 => this%pp%pp%pspecies(isp1)%ptr%name
 namesp2 => this%pp%pp%pspecies(isp2)%ptr%name
 coeff1 => this%vol0(i)
 write (ioutput,4) namesp1,namesp2,coeff1
end do
write (ioutput,*) '-------------------------------------------------------------------'
end if
!%-----------------------------------------------------------
!% Write dB/dP coefficients 
!%-----------------------------------------------------------
if (this%nnzbetavol>0) then
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) ' dBeta/dP virial coefficients                                      '
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,12) 'sps1','sps2','B0vol','B1vol','B2vol'
write (ioutput,*) '-------------------------------------------------------------------'
do i=1,this%nnzbetavol
 isp1 => this%indnzbetavol(1,i)
 isp2 => this%indnzbetavol(2,i)
 namesp1 => this%pp%pp%pspecies(isp1)%ptr%name
 namesp2 => this%pp%pp%pspecies(isp2)%ptr%name
 coeff1 => this%beta0vol(i)
 coeff2 => this%beta1vol(i)
 coeff3 => this%beta2vol(i)
 write (ioutput,13) namesp1,namesp2,coeff1,coeff2,coeff3
end do
write (ioutput,*) '-------------------------------------------------------------------'
end if
!%-----------------------------------------------------------
!% Write dcfi/dP coefficients 
!%-----------------------------------------------------------
if (this%nnzcpzvol>0) then
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) ' dCfi/dP virial coefficients                                      '
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,7) 'sps1','sps2','c'
write (ioutput,*) '-------------------------------------------------------------------'
do i=1,this%nnzcpzvol
 isp1 => this%indnzcpzvol(1,i)
 isp2 => this%indnzcpzvol(2,i)
 namesp1 => this%pp%pp%pspecies(isp1)%ptr%name
 namesp2 => this%pp%pp%pspecies(isp2)%ptr%name
 coeff1 => this%cpzvol(i)
 write (ioutput,14) namesp1,namesp2,coeff1
end do
write (ioutput,*) '-------------------------------------------------------------------'
end if
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) '-------------------------------------------------------------------'
write (ioutput,*) '-------------------------------------------------------------------'
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Nullify local pointers 
!%-----------------------------------------------------------
isp1 => null ()
isp2 => null ()
isp3 => null ()
ifunb => null ()
namesp1 => null ()
namesp2 => null ()
namesp3 => null ()
coeff1 => null ()
coeff2 => null ()
coeff3 => null ()
alpha1 => null ()
alpha2 => null ()
!%-----------------------------------------------------------
if (iserror) goto 10  
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
 
1 format (x,a3,6x,a3,9x,a2,8x,a2,8x,a2,8x,a6,4x,a6)
2 format (2a9,5f10.5)
3 format (x,a3,6x,a3,9x,a6)
4 format (2a9,f9.5)
5 format (x,a3,6x,a3,6x,a3,9x,a5)
6 format (3a9,f9.5)
7 format (x,a3,6x,a3,9x,a1)
8 format (x,a3,6x,a3,8x,a5)
9 format (a36,f7.4)
11 format (a37,i5)
12 format (x,a3,6x,a3,9x,a5,4x,a5,5x,a5)
13 format (2a9,3e10.3)
14 format (2a9,e10.3)

end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_aqphpitzer &
 (this, &
  temp, &
  iserror)
 
implicit none  
!-------------------------------------------------------------------------
!
!   $Description: Update the parameters that depends of the temperature in 
! the object. At the present version, Debye-H�ckel limiting slope is updated 
! according to function defined by Clegg and Whitfield(1991) and 
! Zhang et al. (2005) (page 12)
!
!   $Arguments:
!

type (t_aqueousphasepitzer), intent(inout)       :: this     ! Type aqueous phase Pitzer

real*8, intent(in)                               :: temp     ! Temperature [C] 

logical, intent(out)                             :: iserror  ! iserror=true, then there was an error
!-------------------------------------------------------------------------
!
!   $Pre-cond: The temperature must be in celsius 
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)            :: &
 msg
real*8                        :: &
 tempk
real*8, parameter             :: &
 r0=0.0d0, &
 r100=100.0d0, &
 r60=60.0d0, &
 r25=25.0d0, &
 r0kelvin=273.15d0
!-------------------------------------------------------------------------
!
!   $code
!

 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
if (this%aphi==r0) then
 this%aphi=aphi25
end if 
return
!%------------------------------------------------------------
!% Update the temperature dependence of Aphi according to 
!% expression reported by Clegg and Whitfield (1991)
!%------------------------------------------------------------
if (temp==r25) then
  this%aphi=aphi25
else
 tempk = r0kelvin + temp
 if (temp>=-r60.and.temp<r0) then
  this%aphi=c0aphi*(c1aphi*tempk-c2aphi*dlog(tempk)-(c3aphi/tempk)+c4aphi)
 else if (temp>=r0.and.temp<=r100) then
  this%aphi=c0aphi*(c5aphi-c6aphi*dsqrt(tempk)+   &
  c7aphi*tempk*tempk-c8aphi*tempk*tempk*tempk+c9aphi*tempk**3.5d0)
 end if
end if
!%------------------------------------------------------------
return
10 continue 
print *,'**********************'
print *,'Phase:'
print *,'Name:',this%pp%pp%name
print *,'Service: update_ '
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
subroutine compute_density_aqphpitzer &
    (this        ,&
     c           ,&
     density     ,&
     ismolality, &
     iserror)
implicit none      
!-------------------------------------------------------------------------
!
!   $Description: Compute density according to Monnin (1994)
!
! Monnin, C., 1994. Density calculation and concentration scale conversions 
! for natural waters. Computers & Geosciences 20, 1435-1445.
!
!   $Arguments:
!

type (t_aqueousphasepitzer), intent(in)          :: this         ! Type aqueous phase Pitzer

real*8, dimension(this%pp%pp%numsp), intent(in)  :: c            ! Concentration vector 

real*8,intent(out)                               :: density      ! Liquid density [kgr/m3]

logical,intent(in)                               :: ismolality   !  

logical,intent(out)                              :: iserror      ! if iserror=true, then there was an error
!-------------------------------------------------------------------------
!
!   $Pre-cond: The temperature must be in celsius 
!              Concentration vector must be in molalities 
!              This service is not finished 
!              
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)            :: &
 msg
integer                       :: &
 i, &
 iter  
real*8                        :: &
 mass,    &   ! mass of solutes  
 mass0,   &   ! mass of solutes  
 vol,    &    ! total volume of the solution 
 volex,   &   ! total excess volume of the solution  
 volid, &     ! total volume of fresh water 
 ionstr,&     ! ionic strength 
 densityold 
real*8, pointer               :: &
 mol(:) => null ()
logical                       :: &
 isconvergence 
real*8, parameter             :: &
 r1=1.0d0, &
 r6=1.0d6, &
 r3=1.0d3, &
 tol=1.0d-4, &
 r0=0.0d0
integer, parameter            :: &
 mxiter=50
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=' '
!%------------------------------------------------------------
!% Initialice variables 
!%------------------------------------------------------------
density=r3 
volex=r0
iter=0 
isconvergence=.false. 
call check_pointer_ (mol,this%pp%pp%numsp,.true.)
mol=c
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
do while(.not.isconvergence)
densityold=density 
iter=iter+1
!%------------------------------------------------------------
!% Compute ideal volume of the solution  
!%------------------------------------------------------------
call compute_volid_  (this%pp,volid,mol,this%pp%pp%numsp,iserror)
if (iserror) then
 msg="Error when call compute_volid_ service"
 goto 20
end if 
!%--------------------------------------------------------------
!% Computes the ionic strength
!%--------------------------------------------------------------
call compute_ionstr_ (this%pp,ionstr,mol,iserror)
if (iserror) then
 msg="Error when call compute_ionstr_ service"
 goto 20
end if 
!%------------------------------------------------------------
!% Compute the total mass of solutes in thequeous phase 
!%------------------------------------------------------------
call compute_mass_ (this%pp,mass,mol,this%pp%pp%numsp,iserror)
if (iserror) then
    msg="Error when call compute_mass_ service"
    goto 20
end if 
if (iter==1) mass0=mass 
!%------------------------------------------------------------
!% Compute total excess volume 
!%------------------------------------------------------------
call compute_volex_ (this,volex,mol,ionstr,iserror)
if (iserror) then
 msg="Error when call compute_volex_ service"
 goto 20
end if 
!%------------------------------------------------------------
!% Compute the total volume of the solution 
!% (wich contains 1000 gr of water)
!%------------------------------------------------------------
vol=volid+volex 
!%------------------------------------------------------------
!% Compute the density [kgr / m3] 
!% vol => cm3
!% mass => kgr
!%------------------------------------------------------------
density=(r1+mass)*r6/vol
!%------------------------------------------------------------
!% Check the convergence for molalities 
!%------------------------------------------------------------
if (ismolality) then 
  isconvergence=(dabs(density-densityold)<=tol) 
else
  isconvergence=.true. 
end if 
!%------------------------------------------------------------
!% Check the maximum number of iterations 
!%------------------------------------------------------------
if (iter>mxiter) then
 msg='Error, not convergence in molalities'
 iserror=.true.
 goto 20 
end if 

mol=r3*c/(density-mass0*r3) 

end do 
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
call check_pointer_ (mol,1,.false.)
if (iserror) goto 10 
!%------------------------------------------------------------
return
10 continue 
print *,'************************'
print *,'Phase:'
print *,'Name:',this%pp%pp%name
print *,'Service: compute_density_'
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
subroutine compute_volex_aqphpitzer &
    (this        ,&
     volex, &
     c           ,&
     ionstr      ,&
     iserror)
implicit none      
!-------------------------------------------------------------------------
!
!   $Description: Compute density according to Monnin (1994). 
!
!   $Arguments:
!

type (t_aqueousphasepitzer), intent(in)          :: this         ! Type aqueous phase Pitzer

real*8,intent(out)                               :: volex        ! Total excess volume of a multicomponent 

real*8, dimension(this%pp%pp%numsp), intent(in)  :: c            ! Concentration vector 

real*8,intent(in)                                :: ionstr       ! Ionic strength 

logical,intent(out)                              :: iserror      ! if iserror=true, then there was an error
!-------------------------------------------------------------------------
!
!   $Pre-cond: The temperature must be in celsius 
!              
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)            :: &
 msg
real*8                        :: &
 dhvol, &
 ionstri, &
 bij, &
 b0ij, &
 b1ij, &
 b2ij, &
 cij, &
 mi, &
 mj, &
 z, &
 g1, &
 g2  
integer                       :: &
 i, &
 j, &
 k, &
 nz
character(len=50)             :: &
 nameprop 
real*8, pointer               :: &
 g(:,:) => null ()
real*8, parameter             :: &
 r10=1.0d1, &
 r1=1.0d0, &
 r2=2.0d0, &
 tempref=273.15d0
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=' '
!%------------------------------------------------------------
!% Initialice variables 
!%------------------------------------------------------------
volex=0.0d0 
!%------------------------------------------------------------
!% Compute zc
!%------------------------------------------------------------
call compute_sum_absz_c_ (this%pp,z,c,iserror)
z=z/r2
if (iserror) goto 20 
!%------------------------------------------------------------
!% Compute Debye Huckel term 
!%------------------------------------------------------------
dhvol = aphi25vol * ionstr / bdh
dhvol = dhvol * dlog(r1 + bdh * dsqrt(ionstr))
!%------------------------------------------------------------
!% Compute g functions 
!%------------------------------------------------------------
ionstri=dsqrt(ionstr)
call compute_fbetavol_ (this,g,ionstri)
!%------------------------------------------------------------
do nz=1,this%nnzbetavol
  i = this%indnzbetavol(1,nz)
  j = this%indnzbetavol(2,nz)
  k = this%indnzbetavol(3,nz)
  mi = c(i) 
  mj = c(j)
  b0ij = this%beta0vol(nz)
  b1ij = this%beta1vol(nz)
  b2ij = this%beta2vol(nz)
  g1 = g(k,1)
  g2 = g(k,2)
  bij = b0ij + b1ij * g1 + b2ij * g2 
  volex = volex + bij * mi * mj 
end do  
!%------------------------------------------------------------
!%------------------------------------------------------------
do nz=1,this%nnzcpzvol
  i = this%indnzcpzvol(1,nz)
  j = this%indnzcpzvol(2,nz)
  mi = c(i) 
  mj = c(j)
  cij = this%cpzvol(nz)
  volex = volex + z * cij * mi * mj 
end do  
volex = volex * r2 * rgas * r10 * (tempref + this%pp%pp%tempref)
volex = volex + dhvol 
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
if (iserror) goto 10 
!%------------------------------------------------------------
return
10 continue 
print *,'************************'
print *,'Phase:'
print *,'Name:',this%pp%pp%name
print *,'Service: compute_volex_'
print *, msg
print *,'************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_aqphpitzer &
  (targetobj, &
   sourceobj)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy an object in other object.
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(in) :: sourceobj     ! Type aqueous phase Pitzer

type(t_aqueousphasepitzer), intent(out):: targetobj     ! Type aqueous phase Pitzer
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: The target object must be previously created
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer        :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!%-----------------------------------------------------------
targetobj%nfunb=sourceobj%nfunb
targetobj%nfunj=sourceobj%nfunj
targetobj%nnzbeta=sourceobj%nnzbeta
targetobj%nnztheta=sourceobj%nnztheta
targetobj%nnzlambda=sourceobj%nnzlambda
targetobj%nnzq=sourceobj%nnzq
targetobj%nnzcpz=sourceobj%nnzcpz
targetobj%nnzpsipz=sourceobj%nnzpsipz
targetobj%ithcl=sourceobj%ithcl
targetobj%ithk=sourceobj%ithk
targetobj%ismacinnes=sourceobj%ismacinnes
targetobj%aphi=sourceobj%aphi
!%----------------------------------------------------------
!% Copy attibutes regarding derivatives of virial 
!% coefficients with respect to pressure 
!%----------------------------------------------------------
targetobj%nnzbetavol=sourceobj%nnzbetavol
targetobj%nnzcpzvol=sourceobj%nnzcpzvol
targetobj%nnzvol0=sourceobj%nnzvol0
targetobj%nfunbvol=sourceobj%nfunbvol
!%----------------------------------------------------------
!% Copy beta virial coefficients 
!%----------------------------------------------------------
if (targetobj%nnzbeta>0) then
 call check_pointer_ (targetobj%beta0,targetobj%nnzbeta,.true.)
 call check_pointer_ (targetobj%beta1,targetobj%nnzbeta,.true.)
 call check_pointer_ (targetobj%beta2,targetobj%nnzbeta,.true.)
 call check_pointer_ (targetobj%indnzbeta,3,targetobj%nnzbeta,.true.)
 targetobj%beta0=sourceobj%beta0
 targetobj%beta1=sourceobj%beta1
 targetobj%beta2=sourceobj%beta2
 targetobj%indnzbeta=sourceobj%indnzbeta
end if
!%----------------------------------------------------------
!% Copy betavol virial coefficients 
!%----------------------------------------------------------
if (targetobj%nnzbetavol>0) then
 call check_pointer_ (targetobj%beta0vol,targetobj%nnzbetavol,.true.)
 call check_pointer_ (targetobj%beta1vol,targetobj%nnzbetavol,.true.)
 call check_pointer_ (targetobj%beta2vol,targetobj%nnzbetavol,.true.)
 targetobj%beta0vol=sourceobj%beta0vol
 targetobj%beta1vol=sourceobj%beta1vol
 targetobj%beta2vol=sourceobj%beta2vol
end if
!%----------------------------------------------------------
!% Copy betavol virial coefficients 
!%----------------------------------------------------------
if (targetobj%nnzbetavol>0) then
 call check_pointer_ (targetobj%beta0vol,targetobj%nnzbetavol,.true.)
 call check_pointer_ (targetobj%beta1vol,targetobj%nnzbetavol,.true.)
 call check_pointer_ (targetobj%beta2vol,targetobj%nnzbetavol,.true.)
 call check_pointer_ (targetobj%indnzbetavol,3,targetobj%nnzbetavol,.true.)
 targetobj%beta0vol=sourceobj%beta0vol
 targetobj%beta1vol=sourceobj%beta1vol
 targetobj%beta2vol=sourceobj%beta2vol
 targetobj%indnzbetavol=sourceobj%indnzbetavol
end if
!%----------------------------------------------------------
!% Copy theta virial coefficients 
!%----------------------------------------------------------
if (targetobj%nnztheta>0) then
 call check_pointer_ (targetobj%theta,targetobj%nnztheta,.true.)
 call check_pointer_ (targetobj%indnztheta,5,targetobj%nnztheta,.true.)
  targetobj%theta=sourceobj%theta
 targetobj%indnztheta=sourceobj%indnztheta
end if
!%----------------------------------------------------------
!% Copy lambda virial coefficients 
!%----------------------------------------------------------
if (targetobj%nnzlambda>0) then
 call check_pointer_ (targetobj%lambda,targetobj%nnzlambda,.true.)
 call check_pointer_ (targetobj%indnzlambda,2,targetobj%nnzlambda,.true.)
 targetobj%lambda=sourceobj%lambda
 targetobj%indnzlambda=sourceobj%indnzlambda
end if
!%----------------------------------------------------------
!% Copy psi virial coefficients 
!%----------------------------------------------------------
if (targetobj%nnzpsipz>0) then
 call check_pointer_ (targetobj%psi,targetobj%nnzpsipz,.true.)
 call check_pointer_ (targetobj%indnzpsi,3,targetobj%nnzpsipz,.true.)
 targetobj%psi=sourceobj%psi
 targetobj%indnzpsi=sourceobj%indnzpsi
end if
!%----------------------------------------------------------
!% Copy Cij terms 
!%----------------------------------------------------------
if (targetobj%nnzcpz>0) then
 call check_pointer_ (targetobj%cpz,targetobj%nnzcpz,.true.)
 call check_pointer_ (targetobj%indnzcpz,2,targetobj%nnzcpz,.true.)
 targetobj%cpz=sourceobj%cpz
 targetobj%indnzcpz=sourceobj%indnzcpz
end if
!%----------------------------------------------------------
!% Copy Cvolij terms 
!%----------------------------------------------------------
if (targetobj%nnzcpzvol>0) then
 call check_pointer_ (targetobj%cpzvol,targetobj%nnzcpzvol,.true.)
 call check_pointer_ (targetobj%indnzcpzvol,2,targetobj%nnzcpzvol,.true.)
 targetobj%cpzvol=sourceobj%cpzvol
 targetobj%indnzcpzvol=sourceobj%indnzcpzvol
end if
!%----------------------------------------------------------
!% Copy vol0 terms 
!%----------------------------------------------------------
if (targetobj%nnzvol0>0) then
 call check_pointer_ (targetobj%vol0,targetobj%nnzvol0,.true.)
 call check_pointer_ (targetobj%indnzvol0,2,targetobj%nnzvol0,.true.)
 targetobj%vol0=sourceobj%vol0
 targetobj%indnzvol0=sourceobj%indnzvol0
end if
!%----------------------------------------------------------
if (targetobj%nfunj>0) then
 call check_pointer_ (targetobj%zprod,targetobj%nfunj,.true.)
 targetobj%zprod=sourceobj%zprod
end if
!%----------------------------------------------------------
if (targetobj%nfunb>0) then
 call check_pointer_ (targetobj%alpha,2,targetobj%nfunb,.true.)
 targetobj%alpha=sourceobj%alpha
end if
!%----------------------------------------------------------
if (targetobj%nfunbvol>0) then
 call check_pointer_ (targetobj%alphavol,2,targetobj%nfunbvol,.true.)
 targetobj%alphavol=sourceobj%alphavol
end if
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%*******************Private Services*************************
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
subroutine read_pitzer_data_base_aqphpitzer  &
 (this, &
  fxml, &
  typevirial, &
  iserror, &
  alpha)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read virial coefficients from xml file 
! Also read the virial coefficients for the prediction of
! fluid properties like the density 
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(inout)       :: this          ! Type aqueous phase Pitzer

type(xml_t), intent(in)                          :: fxml          ! Type xml_t 

character(len=*), intent(in)                     :: typevirial    ! Type of virial coefficients  

logical, intent(out)                             :: iserror       ! iserror=true, then there was an error

real*8, pointer, dimension(:,:), optional        :: alpha         ! Alpha coefficients for Beta virial coefficients
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
type(dictionary_t)                :: attributes 
type(xml_t)                       :: context
type(xml_t)                       :: ff
integer        :: status,n
real*8         :: reallocal(1)
character(len=100)      :: &
 id1, &
 id2, &
 id3, &
 id4, &
 namesp
integer                 :: &
 i, &
 j, &
 isp1, &
 isp2, &
 isp3, &
 ivirial
logical                 :: &
 besp1, &
 besp2, &
 besp3, &
 havealpha
integer, pointer        :: &
 idvirial(:,:) => null ()
real*8, pointer         :: &
 virial(:) => null (), &
 alpha1(:,:) => null ()
integer, parameter      :: &
 mxdim1=3, &
 mxdim2=500, &
 mxdim3=2 
real*8, parameter       :: &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false. 
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
havealpha=present(alpha)
!%------------------------------------------------------------
!% Allocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (idvirial,mxdim1,mxdim2,.true.)
call check_pointer_ (virial,mxdim2,.true.)
call check_pointer_ (alpha1,mxdim3,mxdim2,.true.)
!%------------------------------------------------------------
!% Initialice variables 
!%------------------------------------------------------------
ivirial=0
!%------------------------------------------------------------
context=fxml
!%------------------------------------------------------------
call sync_xmlfile (context,status)
!%------------------------------------------------------------
call mark_node(context,path=typevirial,status=status)
!%------------------------------------------------------------
if (status<0) goto 10 
!%------------------------------------------------------------
ff=context
call sync_xmlfile (ff,status)
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
select case (typevirial)
!%------------------------------------------------------------
case ('//b0','//b1','//b2','//cfi','//theta','//lamda','//lambda','//b0vol','//b1vol','//b2vol','//cfivol','//vol0')
 
do
    besp1=.false.
    besp2=.false.
    call get_node(ff,path="coeff",attributes=attributes,status=status)
    if (status<0) exit
 
    id1=''
	call get_value (attributes,"sps1", id1, status)
    id2=''
	call get_value (attributes,"sps2", id2, status)
	id3=''
	call get_value (attributes,"alpha1", id3, status)
    id4=''
	call get_value (attributes,"alpha2", id4, status)
 
    do i=1,this%pp%pp%numsp
      namesp=this%pp%pp%pspecies(i)%ptr%name
      if (namesp==id1) then
          isp1=i
          besp1=.true.
        exit
      end if
    end do
 
    if (besp1) then
     do i=1,this%pp%pp%numsp
      namesp=this%pp%pp%pspecies(i)%ptr%name
      if (namesp==id2) then
        isp2=i
        besp2=.true.
        exit
      end if
     end do
    end if
 
    if (besp1.and.besp2) then
     ivirial=ivirial+1
     id1=''
	 call get_value (attributes,"value", id1, status)
     n=0
     call build_data_array (id1,reallocal,n)
     virial(ivirial)=reallocal(1)
     idvirial(1,ivirial)=isp1
     idvirial(2,ivirial)=isp2
	 n=0
	 reallocal=r0
     call build_data_array (id3,reallocal,n)
	 alpha1(1,ivirial)=reallocal(1)
	 n=0
	 reallocal=r0
     call build_data_array (id4,reallocal,n)
	 alpha1(2,ivirial)=reallocal(1)
     besp1=.false.
     besp2=.false.
    end if
 
end do
!%------------------------------------------------------------
case ('//psi')
 
   do
    besp1=.false.
    besp2=.false.
    besp3=.false.
    call get_node(ff,path="coeff",attributes=attributes,status=status)
 
    if (status<0) exit
 
    id1=''
	call get_value (attributes,"sps1", id1, status)
	id2=''
    call get_value (attributes,"sps2", id2, status)
	id3=''
    call get_value (attributes,"sps3", id3, status)
 
    do i=1,this%pp%pp%numsp
      namesp=this%pp%pp%pspecies(i)%ptr%name
      if (namesp==id1) then
          isp1=i
          besp1=.true.
          exit
      end if
    end do
 
    if (besp1) then
      do i=1,this%pp%pp%numsp
        namesp=this%pp%pp%pspecies(i)%ptr%name
        if (namesp==id2) then
         isp2=i
         besp2=.true.
         exit
        end if
       end do
    end if
 
    if (besp1.and.besp2) then
      do i=1,this%pp%pp%numsp
        namesp=this%pp%pp%pspecies(i)%ptr%name
        if (namesp==id3) then
          isp3=i
          besp3=.true.
          exit
        end if
      end do
    end if
 
 
    if (besp1.and.besp2.and.besp3) then
      ivirial=ivirial+1
      id1=''
	  call get_value (attributes,"value", id1, status)
      n=0
      reallocal=r0
	  call build_data_array (id1,reallocal,n)
      virial(ivirial)=reallocal(1)
      idvirial(1,ivirial)=isp1
      idvirial(2,ivirial)=isp2
      idvirial(3,ivirial)=isp3
      besp1=.false.
      besp2=.false.
      besp3=.false.
    end if
 
end do


case ('//aphi')
    
    call get_node(ff,path="coeff",attributes=attributes,status=status)
    if (status<0) goto 10
    id1=''
	call get_value (attributes,"value", id1, status)
    n=0
    reallocal=r0
	call build_data_array (id1,reallocal,n)
    this%aphi=reallocal(1)
 
end select
!%----------------------------------------------------------------
!% Set the virial coefficients 
!%----------------------------------------------------------------
call set_virial_(this,virial,idvirial,ivirial,typevirial,iserror)
!%------------------------------------------------------------
!% Allocate alpha array if it is present as argument 
!%------------------------------------------------------------
if (havealpha) then 
 call check_pointer_ (alpha,mxdim3,mxdim2,.true.)
 alpha=alpha1
end if
10 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (idvirial,1,1,.false.)
call check_pointer_ (virial,1,.false.)
call check_pointer_ (alpha1,1,1,.false.)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_pitzer_data_base_phreeqc_aqphpitzer  &
 (this, &
  iunit, &
  typevirial, &
  iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read virial coefficients from xml file
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(inout)       :: this          ! Type aqueous phase Pitzer

integer, intent    (in)                          :: iunit         ! Phreeqc data base

character(len=*), intent(in)                     :: typevirial    ! Type of virial coefficients  

logical, intent(out)                             :: iserror       ! iserror=true, then there was an error
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
 id1, &
 id2, &
 id3, &
 namesp
integer                 :: &
 i, &
 j, &
 isp1, &
 isp2, &
 isp3, &
 ivirial
real*8                  :: &
 value 
logical                 :: &
 besp1, &
 besp2, &
 besp3
integer, pointer        :: &
 idvirial(:,:) => null ()
real*8, pointer         :: &
 virial(:) => null ()
integer, parameter      :: &
 mxdim1=3, &
 mxdim2=200 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false. 
!%------------------------------------------------------------
!% Allocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (idvirial,mxdim1,mxdim2,.true.)
call check_pointer_ (virial,mxdim2,.true.)
ivirial=0
!%------------------------------------------------------------
!% Select the type of mineral 
!%------------------------------------------------------------
select case (typevirial)
!%------------------------------------------------------------
case ('-B0','-B1','-B2','-C0','-THETA','-LAMBDA')
 
do
    besp1=.false.
    besp2=.false.
    
	read (iunit,*) id1,id2,value
	
    do i=1,this%pp%pp%numsp
      namesp=this%pp%pp%pspecies(i)%ptr%name
      if (namesp==id1) then
          isp1=i
          besp1=.true.
        exit
      end if
    end do
 
    if (besp1) then
     do i=1,this%pp%pp%numsp
      namesp=this%pp%pp%pspecies(i)%ptr%name
      if (namesp==id2) then
        isp2=i
        besp2=.true.
        exit
      end if
     end do
    end if
 
    if (besp1.and.besp2) then
     ivirial=ivirial+1
     virial(ivirial)=value
     idvirial(1,ivirial)=isp1
     idvirial(2,ivirial)=isp2
     besp1=.false.
     besp2=.false.
    end if
 
end do
!%------------------------------------------------------------
case ('-PSI')
 
   do
    besp1=.false.
    besp2=.false.
    besp3=.false.
    
    read (iunit,*) id1,id2,id3,value
 
    do i=1,this%pp%pp%numsp
      namesp=this%pp%pp%pspecies(i)%ptr%name
      if (namesp==id1) then
          isp1=i
          besp1=.true.
          exit
      end if
    end do
 
    if (besp1) then
      do i=1,this%pp%pp%numsp
        namesp=this%pp%pp%pspecies(i)%ptr%name
        if (namesp==id2) then
         isp2=i
         besp2=.true.
         exit
        end if
       end do
    end if
 
    if (besp1.and.besp2) then
      do i=1,this%pp%pp%numsp
        namesp=this%pp%pp%pspecies(i)%ptr%name
        if (namesp==id3) then
          isp3=i
          besp3=.true.
          exit
        end if
      end do
    end if
 
 
    if (besp1.and.besp2.and.besp3) then
      ivirial=ivirial+1
      virial(ivirial)=value
      idvirial(1,ivirial)=isp1
      idvirial(2,ivirial)=isp2
      idvirial(3,ivirial)=isp3
      besp1=.false.
      besp2=.false.
      besp3=.false.
    end if
 
end do
 
 
end select
!%----------------------------------------------------------------
!% Set the virial coefficients 
!%----------------------------------------------------------------
call set_virial_(this,virial,idvirial,ivirial,typevirial,iserror)
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (idvirial,1,1,.false.)
call check_pointer_ (virial,1,.false.)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_virial_aqphpitzer &
  (this, &
   virial, &
   idvirial, &
   ivirial, &
   typevirial, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the virial coefficients in the object 
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(inout) :: this            ! Type aqueous phase Pitzer

integer, intent(in)                       :: ivirial         ! Number of virial coefficients 

real*8, intent(in), dimension(ivirial)    :: virial          ! Value of the virial coefficients 

integer, intent(in), dimension(:,:)       :: idvirial        ! Indice of the virial coefficients 

character(len=*), intent(in)              :: typevirial      ! Type of virial coefficients (e.g. B0, B1, B2)

logical, intent(out)                      :: iserror         ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 i, &
 j, &
 k, &
 inew, &
 jnew, &
 iold, &
 jold
logical               :: &
 be
real*8                :: &
 z1, &
 z2
real*8, pointer       :: &
 arrayloc(:) => null ()
integer, pointer      :: &
 idarrayloc(:,:) => null ()
character(len=100)    :: &
 msg 
type(t_species), pointer  :: &
 ipspecies => null (), &
 jpspecies => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
if (ivirial==0) return
!%------------------------------------------------------------
select case (typevirial)
case ('//b0','b0','B0','//B0')
 be=.false.
 
 if(this%nnzbeta==0) then
 
  this%nnzbeta = ivirial
  call check_pointer_ (this%beta0,this%nnzbeta,.true.)
  call check_pointer_ (this%beta1,this%nnzbeta,.true.)
  call check_pointer_ (this%beta2,this%nnzbeta,.true.)
  call check_pointer_ (this%indnzbeta,3,this%nnzbeta,.true.)
  this%beta0 = virial (1:ivirial)
  this%indnzbeta = idvirial(:,1:ivirial)
 
 else
 
  call check_pointer_ (this%beta0,this%nnzbeta,.true.)
  do i=1,ivirial
    inew = idvirial(1,i)
    jnew = idvirial(2,i)
   do j=1,this%nnzbeta
     iold = this%indnzbeta (1,j)
     jold = this%indnzbeta (2,j)
     if (iold==inew.and.jold==jnew.or.iold==jnew.and.jold==inew) then
         this%beta0 (j) = virial (i)
         be=.true.
     end if
   end do
    if (.not.be) then
!%--------------
      call check_pointer_ (arrayloc,this%nnzbeta,.true.)
      call check_pointer_ (idarrayloc,3,this%nnzbeta,.true.)
      arrayloc = this%beta0
      idarrayloc = this%indnzbeta
      this%nnzbeta = this%nnzbeta + 1
      call check_pointer_ (this%beta0,this%nnzbeta,.true.)
      call check_pointer_ (this%indnzbeta,3,this%nnzbeta,.true.)
      this%beta0 (1:this%nnzbeta-1) = arrayloc
      this%indnzbeta (:,1:this%nnzbeta-1) = idarrayloc
      this%beta0 (this%nnzbeta) = virial (i)
      this%indnzbeta (1,this%nnzbeta) = idvirial (1,i)
      this%indnzbeta (2,this%nnzbeta) = idvirial (2,i)
!%--------------
      arrayloc = this%beta1
      call check_pointer_ (this%beta1,this%nnzbeta,.true.)
      this%beta1 (1:this%nnzbeta-1) = arrayloc
!%--------------
      arrayloc = this%beta2
      call check_pointer_ (this%beta2,this%nnzbeta,.true.)
      this%beta2 (1:this%nnzbeta-1) = arrayloc
!%--------------
    end if
 
  end do
 
 end if
!%------------------------------------------------------------
case ('//b1','b1','B1','//B1')
 
 be=.false.
 
 if(this%nnzbeta==0) then
 
  this%nnzbeta = ivirial
  call check_pointer_ (this%beta0,this%nnzbeta,.true.)
  call check_pointer_ (this%beta1,this%nnzbeta,.true.)
  call check_pointer_ (this%beta2,this%nnzbeta,.true.)
  call check_pointer_ (this%indnzbeta,3,this%nnzbeta,.true.)
  this%beta1 = virial (1:ivirial)
  this%indnzbeta = idvirial(:,1:ivirial)
 
 else
 
  call check_pointer_ (this%beta1,this%nnzbeta,.true.)
  do i=1,ivirial
    inew = idvirial(1,i)
    jnew = idvirial(2,i)
   do j=1,this%nnzbeta
     iold = this%indnzbeta (1,j)
     jold = this%indnzbeta (2,j)
     if (iold==inew.and.jold==jnew.or.iold==jnew.and.jold==inew) then
         this%beta1 (j) = virial (i)
         be=.true.
     end if
   end do
    if (.not.be) then
!%--------------
      call check_pointer_ (arrayloc,this%nnzbeta,.true.)
      call check_pointer_ (idarrayloc,3,this%nnzbeta,.true.)
      arrayloc = this%beta1
      idarrayloc = this%indnzbeta
      this%nnzbeta = this%nnzbeta + 1
      call check_pointer_ (this%beta1,this%nnzbeta,.true.)
      call check_pointer_ (this%indnzbeta,3,this%nnzbeta,.true.)
      this%beta1 (1:this%nnzbeta-1) = arrayloc
      this%indnzbeta (:,1:this%nnzbeta-1) = idarrayloc
      this%beta1 (this%nnzbeta) = virial (i)
      this%indnzbeta (1,this%nnzbeta) = idvirial (1,i)
      this%indnzbeta (2,this%nnzbeta) = idvirial (2,i)
!%--------------
      arrayloc = this%beta0
      call check_pointer_ (this%beta0,this%nnzbeta,.true.)
      this%beta0 (1:this%nnzbeta-1) = arrayloc
!%--------------
      arrayloc = this%beta2
      call check_pointer_ (this%beta2,this%nnzbeta,.true.)
      this%beta2 (1:this%nnzbeta-1) = arrayloc
!%--------------
    end if
 
  end do
 
 end if
!%------------------------------------------------------------
case ('//b2','b2','B2','//B2')
 
 be=.false.
 
 if(this%nnzbeta==0) then
 
  this%nnzbeta = ivirial
  call check_pointer_ (this%beta0,this%nnzbeta,.true.)
  call check_pointer_ (this%beta1,this%nnzbeta,.true.)
  call check_pointer_ (this%beta2,this%nnzbeta,.true.)
  call check_pointer_ (this%indnzbeta,3,this%nnzbeta,.true.)
  this%beta2 = virial (1:ivirial)
  this%indnzbeta = idvirial(:,1:ivirial)
 
 else
 
  call check_pointer_(this%beta2,this%nnzbeta,.true.)
  do i=1,ivirial
    inew = idvirial(1,i)
    jnew = idvirial(2,i)
   do j=1,this%nnzbeta
     iold = this%indnzbeta (1,j)
     jold = this%indnzbeta (2,j)
       if (iold==inew.and.jold==jnew.or.iold==jnew.and.jold==inew) then
           this%beta2 (j) = virial (i)
           be=.true.
       end if
   end do
    if (.not.be) then
!%--------------
      call check_pointer_ (arrayloc,this%nnzbeta,.true.)
      call check_pointer_ (idarrayloc,2,this%nnzbeta,.true.)
      arrayloc = this%beta2
      idarrayloc = this%indnzbeta
      this%nnzbeta = this%nnzbeta + 1
      call check_pointer_ (this%beta2,this%nnzbeta,.true.)
      call check_pointer_ (this%indnzbeta,3,this%nnzbeta,.true.)
      this%beta2 (1:this%nnzbeta-1) = arrayloc
      this%indnzbeta (:,1:this%nnzbeta-1) = idarrayloc
      this%beta2 (this%nnzbeta) = virial (i)
      this%indnzbeta (1,this%nnzbeta) = idvirial (1,i)
      this%indnzbeta (2,this%nnzbeta) = idvirial (2,i)
!%--------------
      arrayloc = this%beta0
      call check_pointer_ (this%beta0,this%nnzbeta,.true.)
      this%beta0 (1:this%nnzbeta-1) = arrayloc
!%--------------
      arrayloc = this%beta1
      call check_pointer_ (this%beta1,this%nnzbeta,.true.)
      this%beta1 (1:this%nnzbeta-1) = arrayloc
!%--------------
    end if
 
  end do
 
 end if
!%------------------------------------------------------------
case ('//theta','theta','THETA','//THETA')
 if (ivirial==0) return
 this%nnztheta = ivirial
 call check_pointer_ (this%theta,this%nnztheta,.true.)
 call check_pointer_ (this%indnztheta,5,this%nnztheta,.true.)
 this%indnztheta (1:2,:)= idvirial(:,1:ivirial)
 this%theta = virial (1:ivirial)
 do i=1,this%nnztheta
   inew=idvirial(1,i)
   jnew=idvirial(2,i)
   ipspecies => this%pp%pp%pspecies(inew)%ptr
   jpspecies => this%pp%pp%pspecies(jnew)%ptr
   call get_prop_ (ipspecies,z1,'charge',msg,iserror)
   if (iserror) goto 20
   call get_prop_ (jpspecies,z2,'charge',msg,iserror)
   if (iserror) goto 20
   if (dabs(z1)==0.0d0) then
    msg='Error, the electric charge must be different to 0 in species:'
	call add_ (msg,ipspecies%name)
	iserror=.true.
	goto 20
   end if  
   if (dabs(z2)==0.0d0) then
    msg='Error, the electric charge must be different to 0 in species:'
	call add_ (msg,jpspecies%name)
	iserror=.true.
	goto 20
   end if  
   this%indnztheta (3,i) = int(z1 * z2)
   this%indnztheta (4,i) = int(z1 * z1)
   this%indnztheta (5,i) = int(z2 * z2)
 end do
 
 be=.false.
!%---------------------------------------------------------------------
!% Add other anions-anions and cations-cations pair
!%---------------------------------------------------------------------
 do i=1,this%pp%pp%numsp
  do j=1,this%pp%pp%numsp
    if (i/=j) then
	   ipspecies => this%pp%pp%pspecies(i)%ptr
       jpspecies => this%pp%pp%pspecies(j)%ptr
       call get_prop_ (ipspecies,z1,'charge',msg,iserror)
       if (iserror) goto 20
       call get_prop_ (jpspecies,z2,'charge',msg,iserror)
       if (iserror) goto 20
       if ((z1>0.0d0.and.z2>0.0d0).or.(z1<0.0d0.and.z2<0.0d0)) then
          do k=1,this%nnztheta
            inew=this%indnztheta(1,k)
            jnew=this%indnztheta(2,k)
            if ((inew==i.or.jnew==i).and.(inew==j.or.jnew==j)) then
                be=.true.
                exit
            end if
           end do
 
              if (be) then
               be=.false.
              else
               call check_pointer_ (arrayloc,this%nnztheta,.true.)
               call check_pointer_ (idarrayloc,5,this%nnztheta, .true.)
               arrayloc=this%theta
               idarrayloc=this%indnztheta
               this%nnztheta =this%nnztheta +1
               call check_pointer_ (this%theta,this%nnztheta,.true.)
               call check_pointer_ (this%indnztheta,5,this%nnztheta,.true.)
               this%theta(1:this%nnztheta-1)=arrayloc
               this%indnztheta(:,1:this%nnztheta-1)=idarrayloc
               this%indnztheta(1,this%nnztheta)=i
               this%indnztheta(2,this%nnztheta)=j
               this%indnztheta(3,this%nnztheta)=int(z1 * z2)
               this%indnztheta(4,this%nnztheta)=int(z1 * z1)
               this%indnztheta(5,this%nnztheta)=int(z2 * z2)
              end if
 
       end if
    end if
  end do
 end do
 
case ('//cfi','cfi','CFI','//CFi')
 
 if (ivirial==0) return
 this%nnzcpz = ivirial
 call check_pointer_ (this%cpz,this%nnzcpz,.true.)
 call check_pointer_ (this%indnzcpz,2,this%nnzcpz,.true.)
 this%indnzcpz = idvirial(:,1:ivirial)
 do i=1,this%nnzcpz
   inew=idvirial(1,i)
   jnew=idvirial(2,i)
   ipspecies => this%pp%pp%pspecies(inew)%ptr
   jpspecies => this%pp%pp%pspecies(jnew)%ptr
   call get_prop_ (ipspecies,z1,'charge',msg,iserror)
   if (iserror) goto 20
   call get_prop_ (jpspecies,z2,'charge',msg,iserror)
   if (iserror) goto 20
   this%cpz(i) = virial (i)/(2.0d0*dsqrt(dabs(z1*z2)))
 end do
 
case ('//psi','psi','PSI','//PSI')
 
 this%nnzpsipz = ivirial
 call check_pointer_ (this%psi,this%nnzpsipz,.true.)
 call check_pointer_ (this%indnzpsi,3,this%nnzpsipz,.true.)
 this%indnzpsi = idvirial(:,1:ivirial)
 this%psi = virial (1:ivirial)
 
case ('//lambda','LAMBDA','Lambda','//LAMBDA')
 
 this%nnzlambda = ivirial
 call check_pointer_ (this%lambda,this%nnzlambda,.true.)
 call check_pointer_ (this%indnzlambda,2,this%nnzlambda,.true.)
 this%indnzlambda = idvirial(:,1:ivirial)
 this%lambda = virial (1:ivirial)
 
case ('//b0vol','b0vol','B0VOL','//B0Vol')
 
 be=.false.
 
 if(this%nnzbetavol==0) then
 
  this%nnzbetavol = ivirial
  call check_pointer_ (this%beta0vol,this%nnzbetavol,.true.)
  call check_pointer_ (this%beta1vol,this%nnzbetavol,.true.)
  call check_pointer_ (this%beta2vol,this%nnzbetavol,.true.)
  call check_pointer_ (this%indnzbetavol,3,this%nnzbetavol,.true.)
  this%beta0vol = virial (1:ivirial)
  this%indnzbetavol = idvirial(:,1:ivirial)
 
 else
 
  call check_pointer_ (this%beta0vol,this%nnzbetavol,.true.)
  do i=1,ivirial
    inew = idvirial(1,i)
    jnew = idvirial(2,i)
   do j=1,this%nnzbetavol
     iold = this%indnzbetavol (1,j)
     jold = this%indnzbetavol (2,j)
     if (iold==inew.and.jold==jnew.or.iold==jnew.and.jold==inew) then
         this%beta0vol (j) = virial (i)
         be=.true.
     end if
   end do
    if (.not.be) then
!%--------------
      call check_pointer_ (arrayloc,this%nnzbeta,.true.)
      call check_pointer_ (idarrayloc,3,this%nnzbeta,.true.)
      arrayloc = this%beta0vol
      idarrayloc = this%indnzbetavol
      this%nnzbetavol = this%nnzbetavol + 1
      call check_pointer_ (this%beta0vol,this%nnzbetavol,.true.)
      call check_pointer_ (this%indnzbetavol,3,this%nnzbetavol,.true.)
      this%beta0vol (1:this%nnzbetavol-1) = arrayloc
      this%indnzbetavol (:,1:this%nnzbetavol-1) = idarrayloc
      this%beta0vol (this%nnzbetavol) = virial (i)
      this%indnzbetavol (1,this%nnzbetavol) = idvirial (1,i)
      this%indnzbetavol (2,this%nnzbetavol) = idvirial (2,i)
!%--------------
      arrayloc = this%beta1vol
      call check_pointer_ (this%beta1vol,this%nnzbetavol,.true.)
      this%beta1vol (1:this%nnzbetavol-1) = arrayloc
!%--------------
      arrayloc = this%beta2vol
      call check_pointer_ (this%beta2vol,this%nnzbetavol,.true.)
      this%beta2vol (1:this%nnzbetavol-1) = arrayloc
!%--------------
    end if
 
  end do
 
 end if
!%------------------------------------------------------------
case ('//b1vol','b1vol','B1VOL','//B1Vol')
 
 be=.false.
 
 if(this%nnzbetavol==0) then
 
  this%nnzbetavol = ivirial
  call check_pointer_ (this%beta0vol,this%nnzbetavol,.true.)
  call check_pointer_ (this%beta1vol,this%nnzbetavol,.true.)
  call check_pointer_ (this%beta2vol,this%nnzbetavol,.true.)
  call check_pointer_ (this%indnzbetavol,3,this%nnzbeta,.true.)
  this%beta1vol = virial (1:ivirial)
  this%indnzbetavol = idvirial(:,1:ivirial)
 
 else
 
  call check_pointer_ (this%beta1vol,this%nnzbetavol,.true.)
  do i=1,ivirial
    inew = idvirial(1,i)
    jnew = idvirial(2,i)
   do j=1,this%nnzbetavol
     iold = this%indnzbetavol (1,j)
     jold = this%indnzbetavol (2,j)
     if (iold==inew.and.jold==jnew.or.iold==jnew.and.jold==inew) then
         this%beta1vol (j) = virial (i)
         be=.true.
     end if
   end do
    if (.not.be) then
!%--------------
      call check_pointer_ (arrayloc,this%nnzbetavol,.true.)
      call check_pointer_ (idarrayloc,3,this%nnzbetavol,.true.)
      arrayloc = this%beta1vol
      idarrayloc = this%indnzbetavol
      this%nnzbetavol = this%nnzbetavol + 1
      call check_pointer_ (this%beta1vol,this%nnzbetavol,.true.)
      call check_pointer_ (this%indnzbetavol,3,this%nnzbetavol,.true.)
      this%beta1vol (1:this%nnzbetavol-1) = arrayloc
      this%indnzbetavol (:,1:this%nnzbetavol-1) = idarrayloc
      this%beta1vol (this%nnzbetavol) = virial (i)
      this%indnzbetavol (1,this%nnzbetavol) = idvirial (1,i)
      this%indnzbetavol (2,this%nnzbetavol) = idvirial (2,i)
!%--------------
      arrayloc = this%beta0vol
      call check_pointer_ (this%beta0vol,this%nnzbetavol,.true.)
      this%beta0vol (1:this%nnzbetavol-1) = arrayloc
!%--------------
      arrayloc = this%beta2vol
      call check_pointer_ (this%beta2vol,this%nnzbetavol,.true.)
      this%beta2vol (1:this%nnzbetavol-1) = arrayloc
!%--------------
    end if
 
  end do
 
 end if
!%------------------------------------------------------------
case ('//b2vol','b2vol','B2VOL','//B2Vol')
 
 be=.false.
 
 if(this%nnzbetavol==0) then
 
  this%nnzbetavol = ivirial
  call check_pointer_ (this%beta0vol,this%nnzbetavol,.true.)
  call check_pointer_ (this%beta1vol,this%nnzbetavol,.true.)
  call check_pointer_ (this%beta2vol,this%nnzbetavol,.true.)
  call check_pointer_ (this%indnzbetavol,3,this%nnzbetavol,.true.)
  this%beta2vol = virial (1:ivirial)
  this%indnzbetavol = idvirial(:,1:ivirial)
 
 else
 
  call check_pointer_(this%beta2vol,this%nnzbetavol,.true.)
  do i=1,ivirial
    inew = idvirial(1,i)
    jnew = idvirial(2,i)
   do j=1,this%nnzbetavol
     iold = this%indnzbetavol (1,j)
     jold = this%indnzbetavol (2,j)
       if (iold==inew.and.jold==jnew.or.iold==jnew.and.jold==inew) then
           this%beta2vol (j) = virial (i)
           be=.true.
       end if
   end do
    if (.not.be) then
!%--------------
      call check_pointer_ (arrayloc,this%nnzbetavol,.true.)
      call check_pointer_ (idarrayloc,2,this%nnzbetavol,.true.)
      arrayloc = this%beta2vol
      idarrayloc = this%indnzbetavol
      this%nnzbetavol = this%nnzbetavol + 1
      call check_pointer_ (this%beta2vol,this%nnzbetavol,.true.)
      call check_pointer_ (this%indnzbetavol,3,this%nnzbetavol,.true.)
      this%beta2vol (1:this%nnzbetavol-1) = arrayloc
      this%indnzbetavol (:,1:this%nnzbetavol-1) = idarrayloc
      this%beta2vol (this%nnzbetavol) = virial (i)
      this%indnzbetavol (1,this%nnzbetavol) = idvirial (1,i)
      this%indnzbetavol (2,this%nnzbetavol) = idvirial (2,i)
!%--------------
      arrayloc = this%beta0vol
      call check_pointer_ (this%beta0vol,this%nnzbetavol,.true.)
      this%beta0vol (1:this%nnzbetavol-1) = arrayloc
!%--------------
      arrayloc = this%beta1vol
      call check_pointer_ (this%beta1vol,this%nnzbetavol,.true.)
      this%beta1vol (1:this%nnzbetavol-1) = arrayloc
!%--------------
    end if
 
  end do
 
 end if
 
 case ('//cfivol','cfivol','CFIVOL','//CFivol')
 
 if (ivirial==0) return
 this%nnzcpzvol = ivirial
 call check_pointer_ (this%cpzvol,this%nnzcpzvol,.true.)
 call check_pointer_ (this%indnzcpzvol,2,this%nnzcpzvol,.true.)
 this%indnzcpzvol = idvirial(:,1:ivirial)
 do i=1,this%nnzcpzvol
   inew=idvirial(1,i)
   jnew=idvirial(2,i)
   ipspecies => this%pp%pp%pspecies(inew)%ptr
   jpspecies => this%pp%pp%pspecies(jnew)%ptr
   call get_prop_ (ipspecies,z1,'charge',msg,iserror)
   if (iserror) goto 20
   call get_prop_ (jpspecies,z2,'charge',msg,iserror)
   if (iserror) goto 20
   this%cpzvol(i) = virial (i)
 end do
!%------------------------------------------------------------ 
case ('//vol0','vol0','VOL0','//Vol0')
 
 if (ivirial==0) return
 this%nnzvol0 = ivirial
 call check_pointer_ (this%vol0,this%nnzvol0,.true.)
 call check_pointer_ (this%indnzvol0,2,this%nnzvol0,.true.)
 this%indnzvol0 = idvirial(:,1:ivirial)
 do i=1,this%nnzvol0
   inew=idvirial(1,i)
   jnew=idvirial(2,i)
   ipspecies => this%pp%pp%pspecies(inew)%ptr
   jpspecies => this%pp%pp%pspecies(jnew)%ptr
   call get_prop_ (ipspecies,z1,'charge',msg,iserror)
   if (iserror) goto 20
   call get_prop_ (jpspecies,z2,'charge',msg,iserror)
   if (iserror) goto 20
   this%vol0(i) = virial (i)
 end do 
end select
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate and nullify local pointers 
!%------------------------------------------------------------
call check_pointer_ (arrayloc,1,.false.)
call check_pointer_ (idarrayloc,1,1,.false.)
ipspecies => null ()
jpspecies => null ()
if (iserror) goto 10
!%------------------------------------------------------------
return
10 continue 
print *,'************************'
print *,'Phase:'
print *,'Name:', this%pp%pp%name
print *,'Service: set_'
print *, msg
print *,'************************'
iserror=.true.
return 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_act_coeff_loc_aqphpitzer &
 (this, &
  g, &
  c, &
  ionstr, &
  ionstrz, &
  mt, &
  qphi, &
  q, &
  qpri, &
  dqpri, &
  dqphi, &
  indnzq, &
  gclm, &
  iserror)
 
implicit none 
!-------------------------------------------------------------------------
!
!   $Description: Compute activity coeffients vector 
!
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(in)         :: this       ! Type aqueous phase Pitzer

real*8, intent(out)                             :: gclm       ! MacInnes activity coefficients for Cl- 

real*8, intent(in), dimension(this%pp%pp%numsp) :: c          ! Molality vector 

real*8, pointer, dimension(:)                   :: g          ! Activity coefficients vector

real*8, pointer, dimension(:)                   :: qphi       ! qphi matrix (mat. 13)

real*8, pointer, dimension(:)                   :: q          ! Q matrix (mat. 11)

real*8, pointer, dimension(:)                   :: qpri       ! Q' matrix (mat. 12)

real*8, pointer, dimension(:)                   :: dqpri      ! dQ' \ dI matrix (mat. 30)

real*8, pointer, dimension(:)                   :: dqphi      ! dQphi \ dI matrix (mat. 31) 

integer, pointer, dimension(:,:)                :: indnzq     ! Local indices of non-cero components in Q matrices 

real*8, intent(out)                             :: ionstr     ! Ionic strength (eq. 4)

real*8, intent(out)                             :: ionstrz    ! Z (eq. 5)

logical, intent(out)                            :: iserror    ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: 
!
!   $Post-cond:  The matrices Q, Q' and qphi are stored in sparse scheme 
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                        :: &
 mt, &
 osco 
character(len=100)            :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
! 
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call check_pointer_ (g,this%pp%pp%numsp,.true.)
g=1.0d0
!%------------------------------------------------------------
!% Computes the ionic strength (eq. 4)
!%------------------------------------------------------------
call compute_ionstr_(this%pp,ionstr,c,iserror)
if (iserror) goto 10
!%------------------------------------------------------------------------
!% Computes Z (eq. 5) 
!%------------------------------------------------------------------------
call compute_sum_absz_c_(this%pp,ionstrz,c,iserror)
if (iserror) goto 10 
!%------------------------------------------------------------------------
!% Computes M (eq. 6)
!%------------------------------------------------------------------------
call compute_sum_c_(this%pp,mt,c)
!%------------------------------------------------------------------------
!% Computes matrices
!% Q (mat. 12)
!% Q' (mat. 13)
!% Qphi (mat. 14)
!% dQ'/dI (mat. 30)
!% dQphi/dI (mat. 31)
!%------------------------------------------------------------------------
call compute_qmatrices_ (this,qphi,q,qpri,dqpri,dqphi,indnzq,ionstr,msg,iserror)
if (iserror) goto 10 
!%------------------------------------------------------------------------
!% Computes and add
!% ln gamma_dh (eq. A1)
!% ln gamma'_dh (eq. A2)
!%------------------------------------------------------------------------
call compute_dh_ (this,osco,g,ionstr,gclm,iserror)
if (iserror) goto 10
!%------------------------------------------------------------------------
!% Computes and adds eq. 7 and 8 
!%------------------------------------------------------------------------
call compute_q_ (this,g,c,osco,q,qpri,qphi,indnzq,iserror)
if (iserror) goto 10
!%------------------------------------------------------------------------
!% Computes and adds eq. 10
!%------------------------------------------------------------------------
call compute_ql_ (this,osco,c,iserror)
if (iserror) goto 10
!%------------------------------------------------------------------------
!% Computes and adds eq. 9 
!%------------------------------------------------------------------------
call compute_qc_ (this,g,c,ionstrz,osco,iserror)
if (iserror) goto 10
!%------------------------------------------------------------------------
!% Computes and adds eq. 11
!%------------------------------------------------------------------------
call compute_t_(this,g,c,osco,mt,iserror)
if (iserror) goto 10
!%------------------------------------------------------------------------
 
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%pp%name
print *,'Service: compute_act_coeff_loc_'
print *, msg
print *,'*******************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_qmatrices_aqphpitzer &
 (this, &
  qphi, &
  q, &
  qpri, &
  dqpri, &
  dqphi, &
  indnzq, &
  str, &
  msg, &
  iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Compute Q, Q', Qphi matrices
!    
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(in) :: this              ! Type aqueous phase Pitzer

real*8, intent(in)                     :: str               ! Ionic strength 

real*8, pointer, dimension(:)          :: qphi              ! qphi matrix 
  
real*8, pointer, dimension(:)          :: q                 ! Q matrix

real*8, pointer, dimension(:)          :: qpri              ! Q' matrix 

real*8, pointer, dimension(:)          :: dqpri             ! dQ'/ dI matrix 

real*8, pointer, dimension(:)          :: dqphi             ! dQphi / dI matrix 

integer, pointer, dimension(:,:)       :: indnzq            ! Non-zero indices 

character(len=*), intent(out)          :: msg               ! Error message 

logical, intent(out)                   :: iserror           ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                             :: &
 stri, &
 str2, &
 str3, &
 zizj, &
 afi_stri, &
 zizi, &
 zjzj, &
 xii, &
 xij, &
 xjj, &
 dxii, &
 dxij, &
 dxjj, &
 eth, &
 eth_i, &
 eth_i2, &
 ethpri, &
 ethpri_i, &
 eth2pri, &
 A1_2stri, &
 A2_2stri, &
 t1, &
 t2, &
 B0ij, &
 B1ij, &
 B2ij, &
 expo, &
 expo1, &
 expo2, &
 jij, &
 jii, &
 jjj, &
 jpij, &
 jpii, &
 jpjj, &
 j2pij, &
 j2pii, &
 j2pjj, &
 thij, &
 G1, &
 G2, &
 GP1, &
 GP2, &
 G2P1, &
 G2P2, &
 ttij, &
 ttii, &
 ttjj
integer                              :: &
 naqt, &
 i, &
 j, &
 k, &
 nz, &
 nz_b, &
 nz_c, &
 z1, &
 z2, &
 nz_th, &
 nz_lam, &
 funii, &
 funjj, &
 funij
real*8, pointer                      :: &
 g_(:,:) => null(), &
 g_pri_(:,:) => null(), &
 g_2pri_(:,:) => null(), &
 f_(:,:) => null(), &
 j_(:) => null(), &
 j_pri_(:) => null(), &
 j_2pri_(:) => null() 
!-------------------------------------------------------------------------
!
!   $code
!
!%--------------------------------------------------------------------
iserror=.false. 
msg='' 
!%--------------------------------------------------------------------
!% Allocate Q matrices 
!%--------------------------------------------------------------------
call check_pointer_ (qphi,this%nnzq,.true.)
call check_pointer_ (q,this%nnzq,.true.)
call check_pointer_ (qpri,this%nnzq,.true.)
call check_pointer_ (dqpri,this%nnzq,.true.)
call check_pointer_ (dqphi,this%nnzq,.true.)
call check_pointer_ (indnzq,2,this%nnzq,.true.)
!%--------------------------------------------------------------------
stri=dsqrt(str)
str2=str*str
str3=str2*str
afi_stri=2.352d0 * stri
!%--------------------------------------------------------------------
!% Computes g functions (eq. A6, A7 and B5)
!%--------------------------------------------------------------------
call compute_fbeta_ (this, stri, f_, g_, g_pri_, g_2pri_)
!%--------------------------------------------------------------------
!% Computes j functions (eq. A14, A15, B9, B10, B11 and B12)
!%--------------------------------------------------------------------
call compute_fj_ (this, j_, j_pri_, j_2pri_, afi_stri)
!%--------------------------------------------------------------------
nz=0
!%--------------------------------------------------------------------
!% Computes Q, Q', Qphi, dQphi and dQ' matrices
!%--------------------------------------------------------------------
do nz_b=1, this%nnzbeta
 
      nz=nz+1
 
      i=this%indnzbeta(1,nz_b)
      j=this%indnzbeta(2,nz_b)
      k=this%indnzbeta(3,nz_b)
!%-------------
      B0ij=this%beta0(nz_b)
      B1ij=this%beta1(nz_b)
      B2ij=this%beta2(nz_b)
!%-------------
      A1_2stri=this%alpha(1,k)/(2.0d0*stri)
      A2_2stri=this%alpha(2,k)/(2.0d0*stri)
!%-------------
      expo1=f_(k,1)
      expo2=f_(k,2)
!%-------------
      G1=g_(k,1)
      G2=g_(k,2)
!%------------
      Gp1=g_pri_(k,1)
      Gp2=g_pri_(k,2)
!%-------------
      G2p1=g_2pri_(k,1)
      G2p2=g_2pri_(k,2)
!%-------------------------------
!% Eq. A5
!%-------------------------------
      qphi(nz)=B0ij+B1ij*expo1+B2ij*expo2
!%-------------------------------
!%  Eq. B4
!%-------------------------------
      dqphi(nz)=-B1ij*A1_2stri*expo1-B2ij*A2_2stri*expo2
!%-------------------------------
!%  Eq. A3
!%-------------------------------
      q(nz)=B0ij+B1ij*G1+B2ij*G2
!%-------------------------------
!%  Eq. A4
!%-------------------------------
      qpri(nz)=B1ij*(Gp1/str)+B2ij*(Gp2/str)
 
      t1=(str*G2p1*A1_2stri-Gp1)/str2
      t2=(str*G2p2*A2_2stri-Gp2)/str2
!%-------------------------------
!% Eq. B3
!%-------------------------------
      dqpri(nz)=B1ij*t1+B2ij*t2
      indnzq(1,nz)=i
      indnzq(2,nz)=j
 
 
end do
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!% Anions-Anions, Cations-Cations terms
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
do nz_th=1, this%nnztheta
 !%--------------------------------------------------------------------
      nz=nz+1
 !%--------------------------------------------------------------------
      i=this%indnztheta(1,nz_th)
      j=this%indnztheta(2,nz_th)
      funij=this%indnztheta(3,nz_th)
      funii=this%indnztheta(4,nz_th)
      funjj=this%indnztheta(5,nz_th)
      thij=this%theta(nz_th)
 !%--------------------------------------------------------------------
      jij=j_(funij)
      jii=j_(funii)
      jjj=j_(funjj)
      jpij=j_pri_(funij)
      jpii=j_pri_(funii)
      jpjj=j_pri_(funjj)
      j2pij=j_2pri_(funij)
      j2pii=j_2pri_(funii)
      j2pjj=j_2pri_(funjj)
 !%--------------------------------------------------------------------
      zizj=this%zprod(funij)
      zizi=this%zprod(funii)
      zjzj=this%zprod(funjj)
      xij=afi_stri*zizj
      xii=afi_stri*zizi
      xjj=afi_stri*zjzj
!%-----------------------------
!%  Eq. A12
!%-----------------------------
      eth=(zizj/(4.0d0*str))*(jij-0.5d0*jii-0.5d0*jjj)
      q(nz)=thij+eth
      eth_i=eth/str
      eth_i2=eth/str2
!%----------------------------
!% Eq. A13
!%----------------------------
      ethpri=-eth_i+(zizj/(8.0d0*str2))*(xij*jpij-0.5d0*xii*jpii-0.5d0*xjj*jpjj)
      qpri(nz)=ethpri
      qphi(nz)=thij+eth+str*ethpri
      ethpri_i=ethpri/str
!%----------------------------
!%  Eq. B8
!%----------------------------
      ttij=(jpij+xij*j2pij)*xij
      ttii=(jpii+xii*j2pii)*xii
      ttjj=(jpjj+xjj*j2pjj)*xjj
      eth2pri=-3.d0*ethpri_i-eth_i2+(zizj/(16.d0*str3))*(ttij-0.5d0*ttii-0.5d0*ttjj)
      dqpri(nz)=eth2pri
      dqphi(nz)=2.d0*ethpri+str*eth2pri
      indnzq(1,nz)=i
      indnzq(2,nz)=j
 
end do
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!% Neutral ions-cations and Neutral ions-Anions terms
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
do nz_lam=1,this%nnzlambda
 i=this%indnzlambda(1,nz_lam)
 j=this%indnzlambda(2,nz_lam)
 nz=nz+1
!%-------------------------------
!% Add in Q matrix (mat. 12) 
!% the Lij terms 
!%-------------------------------
 q(nz)=this%lambda(nz_lam)
 indnzq(1,nz)=i
 indnzq(2,nz)=j
end do
!%--------------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------------
call check_pointer_ (g_,1,1,.false.)
call check_pointer_ (g_pri_,1,1,.false.)
call check_pointer_ (g_2pri_,1,1,.false.)
call check_pointer_ (f_,1,1,.false.)
call check_pointer_ (j_,1,.false.)
call check_pointer_ (j_pri_,1,.false.)
call check_pointer_ (j_2pri_,1,.false.)
!%--------------------------------------------------------------------
return
10 continue
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dh_aqphpitzer &
  (this, &
   osco, &
   g, &
   ionstr, &
   gclm, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute DH contributions in the osmotic and activity 
! coefficients 
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(in)          :: this      ! Type aqueous phase Pitzer 

real*8, intent(out)                              :: osco      ! Osmotic coefficient

real*8, intent(out)                              :: gclm      ! Activity coefficient of Cl- in the KCl system 

real*8, intent(in)                               :: ionstr    ! Ionic strength 

real*8, intent(out), dimension(this%pp%pp%numsp) :: g         ! Activity coefficients vector 

logical, intent(out)                             :: iserror   ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 real*8                              :: &
 den, &
 dh, &
 sqri, &
 charge
integer                             :: &
 i
character(len=100)                  :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Zeroing 
!%-----------------------------------------------------------
gclm=0.0d0
!%-----------------------------------------------------------
!% Eq. A1 and A2
!%-----------------------------------------------------------
sqri = dsqrt(ionstr)
den = 1.0d0+bdh*sqri
dh = -this%aphi*sqri/den
if (this%pp%ithw/=0) osco = ionstr * dh  
dh = dh - (2.0d0/bdh)*this%aphi*dlog(den)    
!%-----------------------------------------------------------
do i=1,this%pp%pp%numsp
   call get_prop_ (this%pp%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
   if (iserror) goto 10
   g(i) = charge * charge * dh
end do
!%-------------------------------------------------------------
!% Compute activity coefficient for Cl in a KCl system  
!%-------------------------------------------------------------
if (this%ismacinnes) gclm = gclm_ (ionstr,dh)
!%-------------------------------------------------------------
return
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_qc_aqphpitzer &
 (this, &
  g, &
  c, &
  ionstrz, &
  osco, &
  iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes and adds terms that are involving C matrix 
! (eq. 1 and 3)
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(in)            :: this       ! Type aqueous phase Pitzer

real*8, intent(inout)                              :: osco       ! Osmotic coefficient 

real*8, intent(in)                                 :: ionstrz    ! Z (eq. 5)

real*8, intent(inout), dimension(this%pp%pp%numsp) :: g          ! Activity coefficients vector 

real*8, intent(in), dimension(this%pp%pp%numsp)    :: c          ! Molality vector 

logical, intent(out)                               :: iserror    ! iserror=true, then there was an error
 
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
 miCij, &
 mi, &
 mj, &
 Cij, &
 qc, &
 mjCij, &
 mimjCij, &
 charge
integer                 :: &
 i, &
 j, &
 nz
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
!% Zeroing  
!%-----------------------------------------------------------
qc=0.0d0
!%-----------------------------------------------------------
do  nz=1,this%nnzcpz
      i=this%indnzcpz (1,nz)
      j=this%indnzcpz (2,nz)
      mi=c(i)
      mj=c(j)
      Cij=this%cpz(nz)
      mjCij=mj*Cij
      miCij=mi*Cij
      mimjCij=mi*mj*Cij
      qc=qc+mimjCij
      g(i)=g(i)+mjCij*ionstrz
      g(j)=g(j)+miCij*ionstrz
 end do
!%-----------------------------------------------------------
if (this%pp%ithw/=0) then
 osco=osco+ionstrz*qc
end if
!%-----------------------------------------------------------
do i=1,this%pp%pp%numsp
   call get_prop_ (this%pp%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
   if (iserror) goto 10
   g(i)=g(i)+dabs(charge)*qc
end do
!%------------------------------------------------------------
return
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_t_aqphpitzer &
   (this, &
    g, &
    c, &
    osco, &
    m, &
	iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes and adds the 3rd order terms (eq. 11).
! Also converts the logarithm activity coefficients vector. 
!
!   $Arguments:
!
 
  
type (t_aqueousphasepitzer), intent(in)             :: this     ! Type aqueous phase Pitzer

real*8, intent(inout)                               :: osco     ! Osmotic coefficient

real*8, intent(in)                                  :: m        ! eq. 6

real*8, intent(inout), dimension(this%pp%pp%numsp)  :: g        ! Activity coefficients vector

real*8, intent(in), dimension(this%pp%pp%numsp)     :: c        ! Molality vector 

logical, intent(out)                                :: iserror  ! iserror=true, then there was an error          
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                              :: &
 miCij, &
 mi, &
 mj, &
 mk, &
 Psiijk, &
 tril, &
 ctw1, &
 ctw2
integer                             :: &
 i, &
 j, &
 k, &
 nz
real*8, pointer                     :: &
 vector(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false. 

!%------------------------------------------------------
ctw1=m/cwater
ctw2=2.0d0/m
!%------------------------------------------------------
call check_pointer_ (vector,this%pp%pp%numsp,.true.)
!%------------------------------------------------------
tril=0.0d0
!%------------------------------------------------------
!% Compute t and vect
!%------------------------------------------------------
do  nz=1,this%nnzpsipz
 
   Psiijk=this%psi(nz)
   i=this%indnzpsi(1,nz)       
   j=this%indnzpsi(2,nz)      
   k=this%indnzpsi(3,nz)      
   mi=c(i)
   mj=c(j)
   mk=c(k)
   tril=tril+mi*mj*mk*Psiijk
!%----------------------------------------------------------
!% este se utiliza para el c�lculo de los coeficientes de 
!% actividad de las especies acuosas que no sean H2O
!%----------------------------------------------------------
   vector(i)=vector(i)+Psiijk*mj*mk
   vector(j)=vector(j)+Psiijk*mi*mk
   vector(k)=vector(k)+Psiijk*mj*mi
end do
!%----------------------------------------------------
!% 
!%----------------------------------------------------
do i=1,this%pp%pp%numsp
   if (i==this%pp%ithw) then
    osco=osco+tril
    osco=(ctw2*osco)+1.0d0
    g(this%pp%ithw)=-osco*ctw1
    g(this%pp%ithw)=dexp(g(i))
   else
    g(i)=dexp(g(i)+vector(i))
   end if
end do
!%----------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------
call check_pointer_ (vector,1,.false.)
!%----------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_q_aqphpitzer &
(this, &
 g, &
 c, &
 osco, &
 q, &
 qpri, &
 qphi, &
 indnzq, &
 iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes and adds the terms containing Q matrices in osmotic and 
! activity coefficients 
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(in)             :: this    ! Type aqueous phase Pitzer

real*8, intent(inout)                               :: osco    ! Osmotic coefficient 

real*8, intent(inout), dimension(this%pp%pp%numsp)  :: g       ! Acitivity coefficients vector

real*8, intent(in), dimension(this%pp%pp%numsp)     :: c       ! Molality vector

real*8, intent(in), dimension(this%nnzq)            :: q       ! Q matrix

real*8, intent(in), dimension(this%nnzq)            :: qpri    ! Q' matrix 

real*8, intent(in), dimension(this%nnzq)            :: qphi    ! Qphi matrix

integer, intent(in), dimension(2,this%nnzq)         :: indnzq  ! Non-zero indices for Q matrices 

logical, intent(out)                                :: iserror ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                                :: &
 q2phi, &
 q2prim, &
 qij, &
 qpriij, &
 qphiij, &
 charge, &
 mi, &
 mj
integer                              :: &
 i, &
 j, &
 nz 
character(len=100)                  :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
!%-----------------------------------------------------------
q2phi=0.0d0
q2prim=0.0d0
!%-----------------------------------------------------------
do   nz=1,this%nnzq
 
      qij=q(nz)
      qpriij=qpri(nz)
      qphiij=qphi(nz)
      i=indnzq(1,nz)
      j=indnzq(2,nz)
      mi=c(i)
      mj=c(j)
!%------------------------------------------------------------
!% q2_PHi = sUMij(mi*mj*qPHiij)
!%------------------------------------------------------------
      q2phi=q2phi+mi*mj*qphiij
!%------------------------------------------------------------
!%q'2 = sUMij(mi*mj*q'ij)
!%------------------------------------------------------------
      q2prim=q2prim+mi*mj*qpriij
!%------------------------------------------------------------
!%GAMi = GAMi + sUMj(qij*mj)
!%------------------------------------------------------------
      g(i)=g(i)+2.0d0*qij*mj

      g(j)=g(j)+2.0d0*qij*mi
 
end do
!%------------------------------------------------------------
!% Add the above terms to osco and osco = osco + q2_PHi
!%------------------------------------------------------------
if (this%pp%ithw/=0) osco=osco+q2phi
!%------------------------------------------------------------
 do  i=1,this%pp%pp%numsp
     call get_prop_ (this%pp%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
     if (iserror) goto 10
!%------------------------------------------------------------
!% GAMi = GAMi + (chargei^2)*q'2
!%------------------------------------------------------------
     g(i)=g(i)+ charge * charge * q2prim
 
 end do
!%-----------------------------------------------------------
return
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_ql_aqphpitzer &
(this, &
 osco, &
 c, &
 iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes and adds eq. 10
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(in)            :: this    ! Type aqueous phase Pitzer type

real*8, intent(inout)                              :: osco    ! Osmotic coefficient 

real*8, intent(in), dimension(this%pp%pp%numsp)    :: c       ! Molality vector

logical, intent(out)                               :: iserror ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                                :: &
 Lij, &
 mi, &
 mj
integer                              :: &
 i, &
 j, &
 nz 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
!%-----------------------------------------------------------
!% Add contribution of neutral species (only if the water 
!% species was defined in the aqueous phase)
!%-----------------------------------------------------------
if (this%pp%ithw>0.and.this%nnzlambda>0) then 
 
  do   nz=1,this%nnzlambda
      
      i=this%indnzlambda(1,nz)
      j=this%indnzlambda(2,nz)
      Lij=this%lambda(nz)
      mi=c(i)
      mj=c(j)
      osco=osco+mi*mj*Lij
 
  end do
  
  
end if
!%-----------------------------------------------------------
return
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_fbeta_aqphpitzer &
  (this, &
   stri, &
   f_, &
   g_, &
   g_pri_, &
   g_2pri_)
 
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes eq. A6, A7 and B5
!
!   $Arguments:
!
type (t_aqueousphasepitzer), intent(in) :: this     ! Type aqueous phase Pitzer

real*8, pointer, dimension(:,:)         :: f_       ! 

real*8, pointer, dimension(:,:)         :: g_       ! g function

real*8, pointer, dimension(:,:)         :: g_pri_   ! g' function 

real*8, pointer, dimension(:,:)         :: g_2pri_  ! g'' function 

real*8, intent(in)                      :: stri     ! (Ionic strength)^(1/2)
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                              :: &
 alfa1, &
 alfa2, &
 x1, &
 x2, &
 x1q, &
 x2q, &
 x1c, &
 x2c, &
 expo1, &
 expo2
integer                        :: &
 j 
!-------------------------------------------------------------------------
!
!   $code
!
!%---------------------------------------------------------------------
!% Allocate pointers 
!%---------------------------------------------------------------------
call check_pointer_ (f_,this%nfunb,2,.true.)
call check_pointer_ (g_,this%nfunb,2,.true.)
call check_pointer_ (g_pri_,this%nfunb,2,.true.)
call check_pointer_ (g_2pri_,this%nfunb,2,.true.)
!%---------------------------------------------------------------------
do j=1,this%nfunb
 
   alfa1=this%alpha(1,j)
   alfa2=this%alpha(2,j)
 
   x1=alfa1*stri
   x1q=x1*x1
   x1c=x1q*x1
   
   expo1=dexp(-x1)
 
   f_(j,1)=expo1
!%---------------------------------------------------------------------
!% eq. A6
!%--------------------------------------------------------------------- 
   g_(j,1)=2.0d0*(1.0d0-(1.0d0+x1)*expo1)/x1q
!%---------------------------------------------------------------------
!% eq. A7
!%---------------------------------------------------------------------
 
   g_pri_(j,1)=-2.0d0*(1.0d0-(1.0d0+x1+(x1q/2.0d0))*expo1)/x1q

!%---------------------------------------------------------------------
!% eq. B5
!%---------------------------------------------------------------------
 
   g_2pri_(j,1)=4.0d0*(1.0d0-(1.0d0+x1+(x1q/2.0d0)+(x1C/4.0d0))*expo1)/x1c
 
 
!%---------------------------------------------------------------------
!% when the ions aren't monovalents
!%---------------------------------------------------------------------

   if(alfa2/=0.0d0) then
     x2=alfa2*stri
     x2q=x2*x2
     x2c=x2q*x2
 
     expo2=dexp(-x2)
     f_(j,2)=expo2
!%---------------------------------------------------------------------
!% eq. A6
!%---------------------------------------------------------------------
     g_(j,2)=2.0d0*(1.0d0-(1.0d0+x2)*expo2)/x2q
!%---------------------------------------------------------------------
!% eq. A7
!%---------------------------------------------------------------------
     g_pri_(j,2)=-2.0d0*(1.0d0-(1.0d0+x2+(x2q/2.0d0))*expo2)/x2q
!%---------------------------------------------------------------------
!% eq. B5
!%---------------------------------------------------------------------
     g_2pri_(j,2)=4.0d0*(1.0d0-(1.0d0+x2+(x2q/2.0d0)+(x2C/4.0d0))*expo2)/x2c
   end if
!%------------------------------------------------------------
end do
 
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_fbetavol_aqphpitzer &
  (this, &
   g_, &
   stri)
 
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes g functions 
!
!   $Arguments:
!
type (t_aqueousphasepitzer), intent(in) :: this     ! Type aqueous phase Pitzer

real*8, pointer, dimension(:,:)         :: g_       ! g function

real*8, intent(in)                      :: stri     ! (Ionic strength)^(1/2)
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                              :: &
 alfa1, &
 alfa2, &
 x1, &
 x2, &
 x1q, &
 x2q, &
 expo1, &
 expo2
integer                        :: &
 j 
!-------------------------------------------------------------------------
!
!   $code
!
!%---------------------------------------------------------------------
!% Allocate pointers 
!%---------------------------------------------------------------------
call check_pointer_ (g_,this%nfunbvol,2,.true.)
!%---------------------------------------------------------------------
do j=1,this%nfunbvol
 
   alfa1=this%alphavol(1,j)
   alfa2=this%alphavol(2,j)
 
   x1=alfa1*stri
   x1q=x1*x1
   expo1=dexp(-x1)

!%---------------------------------------------------------------------
!% eq. A6
!%--------------------------------------------------------------------- 
   g_(j,1)=2.0d0*(1.0d0-(1.0d0+x1)*expo1)/x1q
!%---------------------------------------------------------------------
!% when the ions aren't monovalents
!%---------------------------------------------------------------------
   if(alfa2/=0.0d0) then
     x2=alfa2*stri
     x2q=x2*x2
     expo2=dexp(-x2)
!%---------------------------------------------------------------------
!% eq. A6
!%---------------------------------------------------------------------
     g_(j,2)=2.0d0*(1.0d0-(1.0d0+x2)*expo2)/x2q
   end if
!%------------------------------------------------------------
end do
 
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_fj_aqphpitzer &
  (this, &
   j_, &
   j_pri_, &
   j_2pri_, &
   afi_stri)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes eq. A14, A15, B9, B10, B11 and B12
!
!   $Arguments:
!
 
type (t_aqueousphasepitzer), intent(in) :: this       ! Type aqueous phase Pitzer

real*8, pointer, dimension(:)           :: j_         ! j fuction 

real*8, pointer, dimension(:)           :: j_pri_     ! j' function 

real*8, pointer, dimension(:)           :: j_2pri_    ! j'' function 

real*8, intent(in)                      :: afi_stri   ! Afi*(Ionic Strength)^(1/2)
 
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
 x2, &
 x3, &
 x4, &
 x, &
 expo, &
 j2p, &
 termxc4, &
 jp_j, &
 jp2_j, &
 s1, &
 s2, &
 s3, &
 s1q, &
 s3q, &
 s1c, &
 zizj, &
 logar, &
 xc4, &
 xc2, &
 td, &
 td1
integer                      :: &
 i, &
 j, &
 k
real*8                       :: &
 d(6)
real*8, parameter            :: &
  e1=4.581d0, &
  e2=0.7237d0, &
  e3=0.0120d0, &
  e4=0.528d0, &
  e12=7.8963d0, &
  e13=0.029025d0, &
  e14=0.1957d0 
!-------------------------------------------------------------------------
!
!   $code
d(1)=4.118d0
d(2)=7.247d0
d(3)=-4.408d0
d(4)=1.837d0
d(5)=-0.251d0
d(6)=0.0164d0
!%------------------------------------------------------------------
!% Allocate pointers 
!%------------------------------------------------------------------
call check_pointer_ (j_,this%nfunj,.true.)
call check_pointer_ (j_pri_,this%nfunj,.true.)
call check_pointer_ (j_2pri_,this%nfunj,.true.)
!%------------------------------------------------------------------
do i=1,this%nfunj
 
zizj = this%zprod(i)
 
x=zizj*afi_stri  
x2=x*x
x3=x2*x
x4=x3*x
 
!cprovi if (x<=0.03d0) then    !For x<= 0.03
if (x<=80.0d0) then    !For x<= 0.03

!%------------------------------------------------------------------
!% Compute s1=sum1 (C(k)/x^k)
!% Compute s2=suma((d(k)*k*(k+1.d0)/(x**(k+2.d0))))
!% Compute s3=suma((d(k)*k/(x**(k+1.d0))))
!%------------------------------------------------------------------
 
     s1=d(6)/x
     s3=6.0d0*s1
     s2=7.0d0*s3
 
     do k=5,1,-1
       s1=(s1+d(k))/x
       s2=(s2+k*(k+1)*d(k))/x
       s3=(s3+k*d(k))/x
     end do
 
 
      s2=s2/x2
      s3=s3/x
      s1q=s1*s1
      s1c=s1q*s1
      s3q=s3*s3
 
      logar=dlog(x)
      expo=dexp(-10.0d0*x2)
!%-------------------------------------
!% A14
!%-------------------------------------
      j_(i)=-(1.0d0/6.0d0)*x2*logar*expo+(1/s1)
!%-------------------------------------
!% eq.B9
!%-------------------------------------
      j_pri_(i)=((10.d0*x2-1.d0)*logar-0.5d0)*(x/3.0d0)*expo+(s3/s1q)
!%-------------------------------------
!% eq.B10
!%-------------------------------------
      j2p=-(s2/s1q)+(2.0d0*s3q/s1c)
      j2p=j2p+((50.d0*x2-200.d0*x4-1.0d0)*logar+20.d0*x2-1.5d0)*expo/3.0d0
      j_2pri_(i)=j2p
!%----------------------------------
 
else       ! For x>= 0.03
 
      xc4=x**e4
      xc2=x**e2
      expo=dexp(-e3*xc4)
      td1=(e1/xc2)*expo
      td=4.0d0+td1
!%-------------------------------------
!% eq. A15
!%-------------------------------------
      j_(i)=x/td
!%-------------------------------------
!% eq. B11
!%-------------------------------------
      j_pri_(i)=(j_(i)/x2)*(x+td1*(e2+e3*e4*xc4)*j_(i))
!%-------------------------------------
!% eq. B12
!%-------------------------------------
      termxc4=e3*e4*xc4
      jp_j=j_pri_(i)/j_(i)
      jp2_j=jp_j*j_pri_(i)
      j2p=e4-e2-termxc4
      j2p=j2p*termxc4/(e2+termxc4)+(jp_j*x-e2)
      j2p=j2p*(jp_j*x-1.0d0)+1.0d0
      j_2pri_(i)=j2p*j_(i)/x2-2.0d0*j_pri_(i)/x+jp2_j
!%-------------------------------------
 
end if
 
 
 
end do
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine assign_fbeta_aqphpitzer &
  (this, &
  alpha, &
  iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Assign beta functions 
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(inout)      :: this     ! Type aqueous phase Pitzer

real*8, intent(in), dimension(2,this%nnzbeta)  :: alpha    ! Alpha coefficients read from Pitzer data base 

logical, intent(out)                           :: iserror  ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 nz, &
 l1, &
 l2, &
 l3, &
 i, &
 j, &
 n1, &
 n2, &
 n3
real*8                 :: &
 z1, &
 z2
real*8, pointer         :: &
 alphaloc(:,:) => null ()
character(len=100)     :: &
 msg 
integer, parameter     :: &
 ndim=1000
real*8, parameter      :: &
 r0=0.0d0, &
 r1=1.0d0, &
 r2=2.0d0, &
 r12=12.0d0, &
 r1_4=1.4d0, &
 r50=50.0d0
!-------------------------------------------------------------------------
!
!   $code
!


!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (alphaloc,2,ndim,.true.)
!%-----------------------------------------------------------
!% Initialice variables
!%-----------------------------------------------------------
this%nfunb=0
l1=0
l2=0
l3=0
!%-----------------------------------------------------------
! Asiggn Beta's functions
!%-----------------------------------------------------------
do nz=1,this%nnzbeta
	      
   i=this%indnzbeta(1,nz)
   j=this%indnzbeta(2,nz)
   call get_prop_ (this%pp%pp%pspecies(i)%ptr,z1,'charge',msg,iserror)
   if (iserror) goto 20 
   call get_prop_ (this%pp%pp%pspecies(j)%ptr,z2,'charge',msg,iserror)
   if (iserror) goto 20

   if (alpha(1,nz)/=r0.or.alpha(2,nz)/=r0) then
           
		   this%nfunb = this%nfunb + 1
		   this%indnzbeta(3,nz) = this%nfunb
           
		   alphaloc (1,this%nfunb) = alpha(1,nz)
		   alphaloc (2,this%nfunb) = alpha(2,nz)

   else
       
	       if (abs(z1)==r1.or.abs(z2)==r1) then
                   	if (l1==0) then 
 						   
						   this%nfunb = this%nfunb + 1
						   this%indnzbeta(3,nz) = this%nfunb
	                       
						   n1 = this%nfunb
                           
						   alphaloc (1,this%nfunb) = r2
						   alphaloc (2,this%nfunb) = r0
						   
						   l1=1
				     else
						   
						   this%indnzbeta(3,nz) = n1
						  
				     end if 	                   
 					   
		    else 
					   
					 if (abs(z1)/=abs(z2)) then
	                         
							 if (l2==0) then 
	                           this%nfunb = this%nfunb + 1
						       this%indnzbeta(3,nz) = this%nfunb
                              
							   n2 = this%nfunb
							   
							   alphaloc(1,this%nfunb)=r2
						       alphaloc(2,this%nfunb)=r50
						       
	                           l2=1
	                         else
							 
							   this%indnzbeta(3,nz)=n2
						  
						     end if 

					  else 
	                       
						     if (l3==0) then
						       this%nfunb = this%nfunb + 1
						       this%indnzbeta(3,nz) = this%nfunb
                               n3=this%nfunb
							   							   
							   alphaloc (1,this%nfunb) = r1_4
						       alphaloc (2,this%nfunb) = r12
						       
	                           l3=1
                             else
					       
                               this%indnzbeta(3,nz)=n3
 						  
						     end if 
					       
					  end if

			end if  
        
	  
   end if	        
	  
	  
end do
!%------------------------------------------------------------
!% Allocate pointer
!%------------------------------------------------------------
call check_pointer_ (this%alpha,2,this%nfunb,.true.)
!%------------------------------------------------------------
this%alpha = alphaloc (1:2,1:this%nfunb)
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (alphaloc,1,1,.false.)
!%------------------------------------------------------------
if (iserror) goto 10 
!%------------------------------------------------------------
return
10 iserror=.true.
return 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine assign_fbetavol_aqphpitzer &
  (this, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Assign betavol functions 
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(inout)         :: this     ! Type aqueous phase Pitzer

logical, intent(out)                              :: iserror  ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 nz, &
 l1, &
 l2, &
 l3, &
 i, &
 j, &
 n1, &
 n2, &
 n3
real*8                 :: &
 z1, &
 z2
real*8, pointer         :: &
 alphaloc(:,:) => null ()
character(len=100)     :: &
 msg 
integer, parameter     :: &
 ndim=1000
!-------------------------------------------------------------------------
!
!   $code
!


!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (alphaloc,2,ndim,.true.)
!%-----------------------------------------------------------
!% Initialice variables
!%-----------------------------------------------------------
this%nfunbvol=0
l1=0
l2=0
l3=0
!%-----------------------------------------------------------
! Asiggn Beta's functions
!%-----------------------------------------------------------
do nz=1,this%nnzbetavol
	      
   i=this%indnzbetavol(1,nz)
   j=this%indnzbetavol(2,nz)
   call get_prop_ (this%pp%pp%pspecies(i)%ptr,z1,'charge',msg,iserror)
   if (iserror) goto 20 
   call get_prop_ (this%pp%pp%pspecies(j)%ptr,z2,'charge',msg,iserror)
   if (iserror) goto 20

	       if (abs(z1)==1.0d0.or.abs(z2)==1.0d0) then
                   	if (l1==0) then 
 						   
						   this%nfunbvol = this%nfunbvol + 1
						   this%indnzbetavol(3,nz) = this%nfunbvol
	                       
						   n1 = this%nfunbvol
                           
						   alphaloc (1,this%nfunbvol) = 2.0d0
						   alphaloc (2,this%nfunbvol) = 0.0d0
						   
						   l1=1
				     else
						   
						   this%indnzbetavol(3,nz) = n1
						  
				     end if 	                   
 					   
		    else 
					   
					 if (abs(z1)/=abs(z2)) then
	                         
							 if (l2==0) then 
	                           this%nfunbvol = this%nfunbvol + 1
						       this%indnzbetavol(3,nz) = this%nfunbvol
                              
							   n2 = this%nfunbvol
							   
							   alphaloc(1,this%nfunbvol)=2.0d0
						       alphaloc(2,this%nfunbvol)=50.0d0
						       
	                           l2=1
	                         else
							 
							   this%indnzbetavol(3,nz)=n2
						  
						     end if 

					  else 
	                       
						     if (l3==0) then
						       this%nfunbvol = this%nfunbvol + 1
						       this%indnzbetavol(3,nz) = this%nfunbvol
                               n3=this%nfunbvol
							   							   
							   alphaloc (1,this%nfunbvol) = 1.4d0
						       alphaloc (2,this%nfunbvol) = 12.0d0
						       
	                           l3=1
                             else
					       
                               this%indnzbetavol(3,nz)=n3
 						  
						     end if 
					       
					  end if

			end if         
	  
	  
end do
!%------------------------------------------------------------
!% Allocate pointer
!%------------------------------------------------------------
call check_pointer_ (this%alphavol,2,this%nfunbvol,.true.)
!%------------------------------------------------------------
this%alphavol = alphaloc (1:2,1:this%nfunbvol)
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (alphaloc,1,1,.false.)
!%------------------------------------------------------------
if (iserror) goto 10 
!%------------------------------------------------------------
return
10 iserror=.true.
return 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine assign_fj_aqphpitzer &
  (this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(inout) :: this    ! Aqueous phase Pitzer type
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: Assign j functions 
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 nz, &
 k, &
 zz, &
 nn, &
 j
real*8, pointer        :: &
 zprodloc (:) => null () 
integer, parameter     :: &
 ndim=20
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
!% Compute number of functions j and save
!%-----------------------------------------------------------
call check_pointer_ (zprodloc,ndim,.true.)
!%------------------------------------------------------------
this%nfunj=1
 
do nz=1,this%nnztheta
 do k=3,5
    zz=this%indnztheta(k,nz)
    do j=1,this%nfunj
      if (zprodloc(j)==real(zz)) then
            nn=j
            goto 10
      end if
     end do
     zprodloc(this%nfunj) = real(zz) 
     nn = this%nfunj
     this%nfunj = this%nfunj + 1
  10 continue     
  this%indnztheta(k,nz)=nn
 end do
end do
!%------------------------------------------------------------
this%nfunj = this%nfunj - 1
!%------------------------------------------------------------
!% Allocate pointer
!%------------------------------------------------------------
call check_pointer_ (this%zprod,this%nfunj,.true.)
this%zprod = zprodloc (1:this%nfunj)
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (zprodloc,1,.false.)
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_ddh_aqphpitzer &
 (this, &
  dg, &
  dosco, &
  ionstr, &
  dionstr, &
  ndimder, &
  dgclm, &
  iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute debye-Huckel contribution to osmotic 
! and activity coefficients
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(in)                      :: this     ! Type aqueous phase Pitzer

integer, intent(in)                                         :: ndimder  ! Number of derivatives 

real*8 , intent(in), dimension(this%pp%pp%numsp)            :: dionstr  ! Derivatives of the ionic strength

real*8 , intent(inout), dimension(this%pp%pp%numsp,ndimder) :: dg       ! Derivatives of the activity coefficients 

real*8 , intent(inout), dimension(ndimder)                  :: dosco    ! Derivatives of the osmotic coefficients 

real*8 , intent(in)                                         :: ionstr   ! Ionic strength

real*8 , pointer, dimension(:)                              :: dgclm    ! Derivatives of the Cl- activity coefficient for KCl system 

logical, intent(out)                                        :: iserror  ! iserror=true, then there was an error

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 i, &
 l
real*8                 :: &
 sqri, &
 den, &
 cte, &
 ddh, &
 charge, &
 daux
character(len=100)                  :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
iserror=.false.
!%----------------------------------------------------------
! Compute eq. 22 and 23
!%----------------------------------------------------------
!%----------------------------------------------------------
!% Compute Debye-H�ckel contribution to osmotic coefficient 
!% derivative
!%----------------------------------------------------------
dosco=0.0d0
sqri = sqrt(ionstr)
den = 1.0d0+bdh*sqri
cte=-this%aphi/(2.0d0*den)
ddh=cte*((3.0d0/sqri)-(bdh/den))
!%----------------------------------------------------------
!% Add dosco = ddH1_dHi(ddh1/di)*(di/dc1)
!%----------------------------------------------------------
if (this%pp%ithw/=0) then
 dosco = ddh*ionstr*dionstr
end if
!%----------------------------------------------------------
!% Compute derivative of modifiel dH contrib. to activities
!%----------------------------------------------------------
do i=1,this%pp%pp%numsp
    call get_prop_ (this%pp%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
    if (iserror) goto 10
    dg(i,1:ndimder)= ddh*charge*charge*dionstr(1:ndimder)
end do
!%------------------------------------------------------------
!% Compute derivatives of the Cl- activity coefficient for 
!% the KCl system 
!%------------------------------------------------------------
if (this%ismacinnes) then
 call check_pointer_ (dgclm,ndimder,.true.)
 dgclm=dgclm_(ionstr,ddh)*dionstr
end if
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
subroutine compute_dq_aqphpitzer &
 (this, &
  dg, &
  dosco, &
  q, &
  qpri, &
  qphi, &
  dqpri, &
  dqphi, &
  indnzq, &
  c, &
  ionstr, &
  dionstr, &
  dc, &
  ndimder, &
  msg, &
  iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes and adds the terms containing Q and dQ
! matrices in the osmotic coefficients and activity coefficients
!
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(in)                      :: this       ! Type aqueous phase Pitzer

integer, intent(in)                                         :: ndimder    ! Number of derivatives 

real*8, intent(inout), dimension(this%pp%pp%numsp,ndimder)  :: dg         ! Derivatives of the activity coefficients 

real*8, intent(inout), dimension(ndimder)                   :: dosco      ! Derivatives of the osmotic coefficient 

real*8, intent(in)                                          :: ionstr     ! Ionic strength 

real*8, intent(in), dimension(ndimder)                      :: dionstr    ! Derivatives of the ionic strength 

real*8, intent(in), dimension(this%nnzq)                    :: q          ! Q matrix 

real*8, intent(in), dimension(this%nnzq)                    :: qpri       ! Q' matrix 

real*8, intent(in), dimension(this%nnzq)                    :: qphi       ! dQphi /dI matrix 

real*8, intent(in), dimension(this%nnzq)                    :: dqpri      ! dQpri /dI matrix 

real*8, intent(in), dimension(this%nnzq)                    :: dqphi      ! dQphi /dI matrix 

real*8, intent(in), dimension(this%pp%pp%numsp,ndimder)     :: dc         ! Derivatives of the molalities 

real*8, intent(in), dimension(this%pp%pp%numsp)             :: c          ! Concentration vector 

integer, intent(in), dimension(2,this%nnzq)                 :: indnzq     ! Non-zero indices in the Q matrices

character(len=*), intent(out)                               :: msg        ! Error message 

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
integer                            :: &
 i, &
 j, &
 l, &
 nz
real*8                             :: &
 mimjdqphi_di, &
 mi, &
 mj, &
 qij, &
 qpij, &
 q_p_i, &
 charge
real*8, pointer                    :: &
 vq(:) => null (), &
 vqphi(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
!% Allocate local pointers 
!%----------------------------------------------------------
call check_pointer_ (vq,this%pp%pp%numsp,.true.)
call check_pointer_ (vqphi,this%pp%pp%numsp,.true.)
q_p_i=0.0d0
mimjdqphi_di=0.0d0
!%---------------------------------------------------------
!% First, compute de q_fhi terms 
!%---------------------------------------------------------
do  nz=1,this%nnzq
    i=indnzq(1,nz)
    j=indnzq(2,nz)
    mi=c(i)
    mj=c(j)
    qij=q(nz)
    qpij=qpri(nz)
!%---------------------------------------------------------
!% Add terms dq_PHi/dC1 to dosco/dC1
!%---------------------------------------------------------
    dg(i,:)=dg(i,:)+2.0d0*(qij*dc(j,:)+qpij*mj*dionstr)
    dg(j,:)=dg(j,:)+2.0d0*(qij*dc(i,:)+qpij*mi*dionstr)
!%---------------------------------------------------------
!% Now compute q'-i=mi(dqij/di)mj2)
!% q_P_i es de F
!%---------------------------------------------------------
    q_p_i=q_p_i+mi*dqpri(nz)*mj
    vq(i)=vq(i)+qpij*mj
    vq(j)=vq(j)+qpij*mi
    mimjdqphi_di= mimjdqphi_di+mi*dqphi(nz)*mj
    vqphi(i)=vqphi(i)+qphi(nz)*mj
    vqphi(j)=vqphi(j)+qphi(nz)*mi
end do
!%------------------------------------------------------
if (this%pp%ithw/=0) then
     do l=1,ndimder
        do j=1,this%pp%pp%numsp
          dosco(l)=dosco(l)+vqphi(j)*dc(j,l)
        end do
        dosco(l)=dosco(l)+mimjdqphi_di*dionstr(l)
     end do
end if
!%---------------------------------------------------------
 do i=1,this%pp%pp%numsp
     call get_prop_ (this%pp%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
     if (iserror) goto 10
	 do l=1,ndimder
        dg(i,l)=dg(i,l)+charge*charge*q_p_i*dionstr(l)
        do j=1,this%pp%pp%numsp
          dg(i,l)=dg(i,l)+charge*charge*vq(j)*dc(j,l)
        end do
     end do
 end do
!%---------------------------------------------------------
!% Deallocate local pointers 
!%---------------------------------------------------------
call check_pointer_ (vq,1,.false.)
call check_pointer_ (vqphi,1,.false.)
!%---------------------------------------------------------
return
10 iserror=.true.
stop
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dqc_aqphpitzer  &
 (this, &
  dg, &
  dosco, &
  c, &
  strz, &
  dc, &
  ndimder, &
  msg, &
  iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes and adds terms containing C matrix in osmotic and 
! activity coefficients 
!
!
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(in)                      :: this     ! Type aqueous phase Pitzer

integer, intent(in)                                         :: ndimder  ! Number of derivatives 

real*8, intent(inout), dimension(this%pp%pp%numsp,ndimder)  :: dg       ! Derivatives of the activity coefficients 

real*8, intent(inout), dimension(ndimder)                   :: dosco    ! Derivatives of the osmotic coefficient 

real*8, intent(in), dimension(this%pp%pp%numsp)             :: c        ! Concentration vector 

real*8, intent(in), dimension(this%pp%pp%numsp,ndimder)     :: dc       ! Derivatives of the concentrations 

real*8, intent(in)                                          :: strz     ! eq. 5 

character(len=*), intent(out)                               :: msg      ! Error message 

logical, intent(out)                                        :: iserror  ! iserror=true, then there was an error 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 i, &
 j, &
 l, &
 nz
real*8                 :: &
 mimjCij, &
 mi, &
 mj, &
 ZCij, &
 Cij, &
 charge, &
 dqC, &
 ZdMk
real*8, pointer         :: &
 v_c(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
!% Allocate local pointers 
!%----------------------------------------------------------
call check_pointer_ (v_c,this%pp%pp%numsp,.true.)
!%----------------------------------------------------------
mimjCij=0.0d0
!%----------------------------------------------------------
do   nz=1,this%nnzcpz
        i=this%indnzcpz(1,nz)
        j=this%indnzcpz(2,nz)
        Cij=this%cpz(nz)
        mi=c(i)
        mj=c(j)
 
                                                        
        v_c(i)=v_c(i)+Cij*mj
        v_c(j)=v_c(j)+Cij*mi
 
        ZCij=strz*Cij
        mimjCij=mimjCij+mi*mj*Cij
        do l=1,ndimder                                  
          dg(i,l)=dg(i,l)+ZCij*dc(j,l)
          dg(j,l)=dg(j,l)+ZCij*dc(i,l)
        end do
 
end do
!%----------------------------------------------------------
do l=1,ndimder
 
 dqc=0.0d0
 zdmk=0.0d0
 
    do i=1,this%pp%pp%numsp     
       call get_prop_ (this%pp%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
       if (iserror) goto 10
       dqc=dqc+v_c(i)*dc(i,l)
       zdmk=zdmk+dabs(charge)*dc(i,l)
    end do
 
    do i=1,this%pp%pp%numsp
       call get_prop_ (this%pp%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
       if (iserror) goto 10
       dg(i,L)=dg(i,L)+v_c(i)*ZdMk+dabs(charge)*dqC
    end do
 
    if (this%pp%ithw/=0) then
	   dosco(L)=dosco(L)+(strz*dqC+Zdmk*mimjCij)
    end if

end do
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (v_c,1,.false.)
!%------------------------------------------------------------
return
10 iserror=.true.
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dql_aqphpitzer  &
 (this, &
  dosco, &
  c, &
  dc, &
  ndimder, &
  msg, &
  iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Computes and adds terms containing Ql matrix in the 
! osmotic coefficients 
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(in)                   :: this       ! Type aqueous phase Pitzer

integer, intent(in)                                      :: ndimder    ! Number of derivatives 

real*8, intent(inout), dimension(ndimder)                :: dosco      ! Derivatives of the osmotic coefficients 

real*8, intent(in), dimension(this%pp%pp%numsp)          :: c          ! Molality vector 

real*8, intent(in), dimension(this%pp%pp%numsp,ndimder)  :: dc         ! Derivatives of the molalities

logical, intent(out)                                     :: iserror    ! iserror=true, then there was an error

character(len=*), intent(out)                            :: msg        ! Error message 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 i, &
 j, &
 l, &
 nz
real*8                 :: &
 mi, &
 mj, &
 Lij
real*8, pointer         :: &
 dmi(:) => null (), &
 dmj(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
if (this%pp%ithw>0.and.this%nnzlambda>0) then
!%----------------------------------------------------------
!% Allocate local pointers 
!%----------------------------------------------------------
 call check_pointer_ (dmi,ndimder,.true.)
 call check_pointer_ (dmj,ndimder,.true.)
!%----------------------------------------------------------
 do  nz=1,this%nnzlambda
    i=this%indnzlambda(1,nz)
    j=this%indnzlambda(2,nz)
    mi=c(i)
    mj=c(j)
    dmi=dc(i,:)
    dmj=dc(j,:)
    Lij=this%lambda(nz)
    dosco=dosco+(dmi*mj+dmj*mi)*Lij
 end do
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
 call check_pointer_ (dmi,1,.false.)
 call check_pointer_ (dmj,1,.false.)
end if
!%------------------------------------------------------------
return
10 iserror=.true.
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dt_aqphpitzer  & 
 (this, &
  dg, &
  dosco, &
  c, &
  g, &
  dc, &
  dm, &
  ndimder)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute and adds 3th order terms in the osmotic and 
! activity coefficients. Also it changes the logarithm. 
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(in)                     :: this     ! Type aqueous phase Pitzer

integer, intent(in)                                        :: ndimder  ! Number of derivatives 

real*8, intent(inout), dimension(this%pp%pp%numsp,ndimder) :: dg       ! Derivatives of the activity coefficients 

real*8, intent(inout), dimension(ndimder)                  :: dosco    ! Derivatives of the osmotic coefficient

real*8, intent(in), dimension(this%pp%pp%numsp)            :: c        ! Concentration vector 

real*8, intent(in), dimension(this%pp%pp%numsp)            :: g        ! Activity coefficients 

real*8, intent(in), dimension(this%pp%pp%numsp,ndimder)    :: dc       ! Derivatives of the concentrations 

real*8, intent(in), dimension(ndimder)                     :: dm       ! dM (eq. 23)

 
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
 l, &
 k, &
 nz
real*8                              :: &
 tijk, &
 mi, &
 mj, &
 mk, &
 dil, &
 djl, &
 dkl 
!-------------------------------------------------------------------------
!
!   $code
!


!%---------------------------------------------------------
do  nz=1,this%nnzpsipz
 
     i=this%indnzpsi(1,nz)          
     j=this%indnzpsi(2,nz)          
     k=this%indnzpsi(3,nz)          
     mi=c(i)                 
     mj=c(j)                  
     mk=c(k)                  
     tijk=this%psi(nz)
 
 
       do l=1,ndimder
 
         dil=dc(i,l)
         djl=dc(j,l)
         dkl=dc(k,l)
 
 
         dg(i,l)=dg(i,l)+tijk*(mj*dkl+mk*djl)
 
         dg(j,l)=dg(j,l)+tijk*(Mi*dkl+Mk*dil)
 
         dg(k,l)=dg(k,l)+tijk*(Mj*dil+Mi*djl)
 
 
 
        if (this%pp%ithw/=0) then    
         dosco(l)=dosco(l)+(Mi*Mj*dkl+Mi*djl*Mk+dil*Mj*Mk)*tijk
        end if
 
 
 
       end do
 
 
 
 end do
!%-----------------------------------------------------------
if (this%pp%ithw/=0) then
   dg(this%pp%ithw,1:ndimder)=-(2.0d0/cwater)*(dosco+0.5d0 * dm)
end if
!%-----------------------------------------------------------
do i=1,this%pp%pp%numsp
   dg(i,1:ndimder) = g(i) * dg(i,1:ndimder)
end do
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_consistence_aqphpitzer  & 
 (this, &
  ioutput, &
  iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Check the consistence between chemical system and 
! virial coefficients data base.  
!
!   $Arguments:
!
 
type(t_aqueousphasepitzer), intent(in)                     :: this     ! Type aqueous phase Pitzer

integer, intent(in)                                        :: ioutput  ! Output unit 

logical, intent(out)                                       :: iserror  ! if true, then there was an error
 
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
 isps1, &
 isps2, &
 ilast1, &
 ilast2  
real*8                              :: &
 coeff
character(len=100)                  :: &
 msg, &
 label, &
 namesps1, &
 namesps2 
real*8, parameter                   :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Check BETA0 virial coefficients 
!%-----------------------------------------------------------
do i=1,this%nnzbeta
  isps1=this%indnzbeta(1,i)
  isps2=this%indnzbeta(2,i)
  namesps1=this%pp%pp%pspecies(isps1)%ptr%name
  namesps2=this%pp%pp%pspecies(isps2)%ptr%name
  coeff=this%beta0(i)
  label='WARNING, BETA0 virial coefficient between '
  if (coeff==r0) then
    call lastletter_ (ilast1,label)
    call lastletter_ (ilast2,namesps1)
    label(ilast1+2:ilast1+ilast2+2)=namesps1
    call lastletter_ (ilast1,label)
    label(ilast1+2:ilast1+5)='and'
    call lastletter_ (ilast1,label)
    call lastletter_ (ilast2,namesps2)
    label(ilast1+2:ilast1+ilast2+2)=namesps2
    call lastletter_(ilast1,label)
    label(ilast1+2:ilast1+29)='not found in the data base.'
    call lastletter_(ilast1,label)
    write (ioutput,'(<ilast1>a)') label(1:ilast1)
  end if
end do
!%-----------------------------------------------------------
!% Check BETA1 virial coefficients 
!%-----------------------------------------------------------
do i=1,this%nnzbeta
  isps1=this%indnzbeta(1,i)
  isps2=this%indnzbeta(2,i)
  namesps1=this%pp%pp%pspecies(isps1)%ptr%name
  namesps2=this%pp%pp%pspecies(isps2)%ptr%name
  coeff=this%beta1(i)
  label='WARNING, BETA1 virial coefficient between '
  if (coeff==r0) then
    call lastletter_ (ilast1,label)
    call lastletter_ (ilast2,namesps1)
    label(ilast1+2:ilast1+ilast2+2)=namesps1
    call lastletter_ (ilast1,label)
    label(ilast1+2:ilast1+5)='and'
    call lastletter_ (ilast1,label)
    call lastletter_ (ilast2,namesps2)
    label(ilast1+2:ilast1+ilast2+2)=namesps2
    call lastletter_(ilast1,label)
    label(ilast1+2:ilast1+29)='not found in the data base.'
    call lastletter_(ilast1,label)
    write (ioutput,'(<ilast1>a)') label(1:ilast1)
  end if
end do
!%-----------------------------------------------------------
!% Check BETA2 virial coefficients 
!%-----------------------------------------------------------  
do i=1,this%nnzbeta
  isps1=this%indnzbeta(1,i)
  isps2=this%indnzbeta(2,i)
  namesps1=this%pp%pp%pspecies(isps1)%ptr%name
  namesps2=this%pp%pp%pspecies(isps2)%ptr%name
  coeff=this%beta2(i)
  label='WARNING, BETA2 virial coefficient between '
  if (coeff==r0) then
    call lastletter_ (ilast1,label)
    call lastletter_ (ilast2,namesps1)
    label(ilast1+2:ilast1+ilast2+2)=namesps1
    call lastletter_ (ilast1,label)
    label(ilast1+2:ilast1+5)='and'
    call lastletter_ (ilast1,label)
    call lastletter_ (ilast2,namesps2)
    label(ilast1+2:ilast1+ilast2+2)=namesps2
    call lastletter_(ilast1,label)
    label(ilast1+2:ilast1+29)='not found in the data base.'
    call lastletter_(ilast1,label)
    write (ioutput,'(<ilast1>a)') label(1:ilast1)
  end if
end do
!%-----------------------------------------------------------
!% Check THETA virial coefficients 
!%-----------------------------------------------------------
do i=1,this%nnztheta  
  isps1=this%indnztheta(1,i)
  isps2=this%indnztheta(2,i)
  namesps1=this%pp%pp%pspecies(isps1)%ptr%name
  namesps2=this%pp%pp%pspecies(isps2)%ptr%name
  coeff=this%theta(i)
  label='WARNING, THETA virial coefficient between'
  if (coeff==r0) then
    call lastletter_ (ilast1,label)
    call lastletter_ (ilast2,namesps1)
    label(ilast1+2:ilast1+ilast2+2)=namesps1
    call lastletter_ (ilast1,label)
    label(ilast1+2:ilast1+5)='and'
    call lastletter_ (ilast1,label)
    call lastletter_ (ilast2,namesps2)
    label(ilast1+2:ilast1+ilast2+2)=namesps2
    call lastletter_(ilast1,label)
    label(ilast1+2:ilast1+29)='not found in the data base.'
    call lastletter_(ilast1,label)
    write (ioutput,'(<ilast1>a)') label(1:ilast1)
  end if
end do 
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private function*****************************
!%************************************************************
!%************************************************************
!%************************************************************
double precision function gclm_aqphpitzer (ionstr, dhterm)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Activity coefficient for KCl for MAcInnes scaled
!
!
!   $Arguments:
!
 
real*8, intent(in)                 :: ionstr  ! Ionic strength

real*8, intent(in)                 :: dhterm  ! Derivative of the Debye-H�ckel term with respect to ionic strength
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                 :: &
 x, &
 xxx, &
 yyy 
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
x=2.0d0*dsqrt(ionstr)
xxx=-2.0d0*(1.0d0-(1.0d0+x+0.5d0*x*x)*dexp(-x))/(x*x)
xxx=mtb1kcl*xxx/ionstr
yyy=2.0d0*(1.0d0-(1.0d0+x)*dexp(-x))/(x*x)
yyy=mtb0kcl+mtb1kcl*yyy
gclm_aqphpitzer = dhterm+ionstr*ionstr*xxx+ionstr*(2.0d0*yyy+ionstr*mtc0kcl)+ionstr*ionstr*mtc0kcl/2.0d0
gclm_aqphpitzer = dexp(gclm_aqphpitzer)
!%------------------------------------------------------------
return
end function
!%************************************************************
!%***************Private function*****************************
!%************************************************************
!%************************************************************
!%************************************************************
double precision function dgclm_aqphpitzer (ionstr,ddh)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the derivatives of the gamma +- clk for KCl 
! solution with respect to ionic strength (for MaInnes scaled)
!
!   $Arguments:
!
 
real*8, intent(in)     :: ionstr ! Ionic strength 

real*8, intent(in)     :: ddh    ! Derivative of the modified Debye-H�ckel 
                                 ! term with respect to ionic strength
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                 :: &
 x, &
 xx, &
 xxx, &
 yyy, &
 dxxx, &
 expx, &
 t1, &
 t2, &
 sqri 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
sqri=dsqrt(ionstr)
x=2.0d0*sqri
expx=dexp(-x)
xx=x*x
xxx=x*x*x
t1=2.0d0*x-(xx/2.0d0)
t2=1.0d0+x-(xx/2.0d0)
dxxx=(x*t1*expx-2.0d0*(1.0d0-t2*expx))/(xxx*sqri)
t1=1.0d0+2.0d0*x-xx/2.0d0
yyy=(1.0d0-t1*expx)/xx
t1=2.0d0*ionstr*mtb1kcl*dxxx
t2=2.0d0*(mtb0kcl+mtb0kcl*yyy)+3.0d0*ionstr*mtc0kcl
dgclm_aqphpitzer=ddh+t1+t2
!%------------------------------------------------------------
return
end function
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_aqueousphasepitzer