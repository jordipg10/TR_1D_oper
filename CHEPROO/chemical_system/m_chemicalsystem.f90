module m_chemicalsystem
!-------------------------------------------------------------------------
!
!   $Description: The entire description of the geochemical system is contained 
! in the chemical system class (except the composition data). 
! Chemical system class encapsulates the species, phases and reactions that define the chemical problem. 
! Within the chemical system class the different information about geochemical system is encapsulated 
! in other classes like species, phase, reaction and reaction rate law classes (see figure \ref{fig:chemsys}).  
! In addition, other data types as the chemical base and components definitions are also encapsulated in the chemical system class. 
! Within it, the species class constitutes the basic entity and it is shared by the other classes. 
! Thus, at the same time one species object could be shared by phase, reaction and reaction rate law objects.
! Chemical system class offers varied functions such as to give information about chemical base as of chemical 
! speciation or derivate calculations.
! State variables are not stored within the chemical system object and for this reason to carry out certain operations 
! implies input them like dummy arguments. 
! Several speciation algorithms have been implemented in the chemical system class and the necessary calculations are 
! carried out by the different classes performed in the chemical system class (e.g. phase class computes thermodynamic activity coefficients). 
!
!   $Use: m_reaction.
! m_reactionratelaw.
! m_surface.
! m_phase.
! m_species.
! m_parentchemicalsystem.
! m_chemicalsystemclassic.
! m_chemicalsystemretraso.
! m_chemicalsystemmolins.
! flib_xpath.
! flib_sax.
! m_general_tools_cheproo.
! m_constants. 
!
!   $Author: Sergio Andrés Bea Jofré 
!
!   $License:
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_reaction
use m_reactionratelaw
use m_surface
use m_phase
use m_species
use m_parentchemicalsystem
use m_chemicalsystemclassic
use m_chemicalsystemretraso
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%--------------------------------------------------------
!%--------------------------------------------------------
private                           ::
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
public                            :: &
create_ &                        ! Create the chemical system object. 
,read_xml_ &                     ! Read and initialice the chemical system object from a xml file. 
,read_txt_ &                     ! Read and initilice the chemical system object from an ascii file 
                                 ! (at the moment read the "che.inp", and thermodynamic databases "master25.dat", 
								 ! "pitzer.xml" and "kinetic.dat"
,destroy_ &                      ! Destroy the chemical system object. 
,set_ &                          ! Set attributes in the chemical system object. 
,specia_ &                       ! Specia an aqueous solution according different constraint. 
,select_cmob_ &                  ! Select the mobile species in a concentration vector. 
,select_caq_ &                   ! Select the aqueous species in a concentration vector. 
,select_cgas_ &                  ! Select the gas species in a concentration vector.  
,select_cads_ &                  ! Select the sorbed species in a concentration vector.   
,select_cmin_ &                  ! Select the mineral species in a concentration vector.  
,compute_usktrk_ &               ! Compute ukin=U*Skt*rk. 
,compute_dusktrk_ &              ! Compute dukin/dc1 = U*Skt* drk/dc1. 
,compute_umob_ &                 ! Compute umob=Umob*cmob.
,compute_uads_ &                 ! Compute uads=Ud*cads.
,compute_dcmob_ &                ! Compute dcmob/dc1.
,compute_dcads_ &                ! Compute dcads/dc1.
,compute_dumob_ &                ! Compute dumob/dc1=Umob*dcmob/dc1. 
,compute_duads_ &                ! Compute duads/dc1=Uads*dcads/dc1. 
,compute_iumob_ &                ! Compute umobi=Umobi*cmob
,compute_ith_total_ &            ! Compute ti=Ui*c
,compute_min_sat_index_ &        ! Compute vector of IAPmin/Kmin for minerals. 
,compute_mass_salt_ &            ! Compute the mass of salt [kgr] in the aqueous phase from molality vector. 
,compute_dmassSalt_dc_ &         ! Compute the derivative mass of salt [kgr] in the aqueous phase. 
,compute_omgwcryst_ &            ! Compute the mass of water fixed in the hydrated minerals. 
,compute_from_setre_ &           ! Compute x as: Transpose(Se)*x=b
,compute_density_ &              ! Compute density (and derivatives if asked) for any phase
,get_hashcompz_ &                 ! Return the hash index to identify the components definition. 
,swith_base_ &                   ! Switch the set of primary species in the chemical system object.          
,compute_alkalinity_ &           ! Compute alkalinity from concentration vector. 
,make_lin_trf_ &                 ! 1) Compute the matricial product vector2=U*vector1 
                                 ! 2) Compute the matricial product array2=U*array1.
,get_sp_index_ &                 ! Return the global indice of some species according its name. 
,get_new_chemsys_ &              ! Return a new chemical system from the name of phases and surfaces.  
,get_new_chemical_base_ &        ! Return a new set of chemical base determined from the concentration vector.  
,get_aq_sps_pos_ &               ! Return the first and last indice of the species corresponding to the aqueous phase. 
,get_gas_sp_index_ &             ! Return the global indices of gas species. 
,get_pha_sps_pos_ &              ! Return the first and last indice of the species corresponding to any phase. 
,get_surf_index_ &               ! Return the global indice of some surface (including the number of sites) according its name. 
,get_chem_info_ &                ! Return general chemical information about the chemical system.
,get_surf_model_ &               ! Return the thermodynamic model of the ith surface object. 
,change_chem_unit_ &             ! Change chemical units. 
,update_ &                       ! 1) Update the chemical system obejct to some temperature. 
                                 ! 2) Update the porosity from the mineral concentrations. 
,write_ &                        ! 1) Write in ascii all attributes ancapsulated in the chemical system object. 
                                 ! 2) Write in ascii the speciation information.
,assignment(=) &                 ! Copy a chemical system object in other chemical system object. 
,check_inv_points_ &             ! Check and determine the set of reactions that constraint the water activity. 
,set_iswcompbal_ &               ! Set if the water component mass balance must be evaluated in the speciation algorithms.
,set_iswriteinfo_ &              ! Set if the chemical system object must write the Newton-Raphson information during the speciation processes. 
,update_min_area_ &              ! Update the reactive surface of all minerals from its  concentration. 
,update_txoh_                    ! Update TXOH
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
integer, parameter               :: &
chemsysclassic =1, &
chemsysretraso =2
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!% Type chemical system pointer
!%--------------------------------------------------------
type, public:: t_pchemsys

type(t_chemicalsystem), pointer:: ptr 

end type
!%--------------------------------------------------------
!%--------------------------------------------------------
!% Type definition 
!%--------------------------------------------------------
!%--------------------------------------------------------
type, public::t_chemicalsystem
 
private                                    ::
 
type(t_parentchemicalsystem), pointer      :: pp               ! Parent chemical system 
 
type(t_chemicalsystemclassic), pointer     :: pchemsysclassic  ! Classic chemical system 
 
type(t_chemicalsystemretraso), pointer     :: pchemsysretraso  ! Retraso chemical system 
 
integer                                    :: itype            ! Specialization indice 
 
 
end type t_chemicalsystem
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_iswcompbal_
 
module procedure set_iswcompbal_chemsys 
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_iswriteinfo_
 
module procedure set_iswriteinfo_chemsys 
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface create_
 
module procedure create_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_xml_
 
module procedure read_xml_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_txt_
 
module procedure read_txt_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface specia_
 
module procedure equilibrate_chemsys
module procedure specia_from_solution_type_chemsys
module procedure add_and_specia_chemsys
module procedure specia_from_u_chemsys
module procedure specia_from_cpri_chemsys
module procedure specia_eqmin_from_setre_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_usktrk_
 
module procedure compute_usktrk1_chemsys
module procedure compute_usktrk2_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dusktrk_
 
module procedure compute_dusktrk1_chemsys
module procedure compute_dusktrk2_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_ith_total_
 
module procedure compute_ith_total_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_umob_
 
module procedure compute_umob_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_uads_
 
module procedure compute_uads_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dumob_
 
module procedure compute_dumob1_chemsys
module procedure compute_dumob2_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_duads_
 
module procedure compute_duads_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_iumob_
 
module procedure compute_iumob_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cmob_
 
module procedure select_cmob_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_caq_
 
module procedure select_caq_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cgas_
 
module procedure select_cgas_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cads_
 
module procedure select_cads_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cmin_
 
module procedure select_cmin_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dcmob_
 
module procedure compute_dcmob1_chemsys
module procedure compute_dcmob2_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dcads_
 
module procedure compute_dcads_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface switch_base_
 
module procedure switch_base_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface make_lin_trf_
 
module procedure make_lin_trf_vector_chemsys
module procedure make_lin_trf_array_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface destroy_
 
module procedure destroy_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_
 
module procedure set_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_chem_info_
 
module procedure get_chem_info_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_hashcompz_
 
module procedure get_hashcompz_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_sp_index_
 
module procedure get_sp_index_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_new_chemsys_
 
module procedure get_new_chemsys_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_new_chemical_base_
 
module procedure get_new_chemical_base_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_aq_sps_pos_
 
module procedure get_aq_sps_pos_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_min_sat_index_
 
module procedure compute_min_sat_index_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_mass_salt_
 
module procedure compute_mass_salt_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dmassSalt_dc_
 
module procedure compute_dmassSalt_dc_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_omgwcryst_
 
module procedure compute_omgwcryst_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_from_setre_
 
module procedure compute_from_setre_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_alkalinity_
 
module procedure compute_alkalinity_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_surf_index_
 
module procedure get_surf_index_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface change_chem_unit_
 
module procedure change_chem_unit_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface update_
 
module procedure update_temp_param_chemsys
module procedure update_porosity_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface write_
 
module procedure write_chemsys
module procedure write_speciation_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface assignment (=)
 
module procedure copy_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface check_inv_points_
 
module procedure check_inv_points_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface update_min_area_
 
module procedure update_min_area_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cmineq_
 
module procedure select_cmineq_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface update_txoh_
 
module procedure update_txoh_chemsys
 
end interface


!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_density_
 
module procedure compute_density_chemsys
 
end interface




!%--------------------------------------------------------
!%--------------------------------------------------------
!%********************************************************
!%********************************************************
!%********************************************************
!%********************************************************
!%******************Private subroutines*******************
!%********************************************************
!%********************************************************
!%********************************************************
!%********************************************************
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_xml_loc_
 
module procedure read_xml_loc_chemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface begin_element_handler
 
module procedure begin_element_handler
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
contains
!%--------------------------------------------------------
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_chemsys &
   (this)
 
 implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create the chemical system object and initialice attributes. 
! Also create the parent chemical system class.
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout) :: this 
 
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
this%pchemsysclassic => null ()
this%pchemsysretraso => null ()
!%-------------------------------------------------------------
this%itype = 0
allocate(this%pp)
!%-------------------------------------------------------------
call create_ (this%pp)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_chemsys &
   (this, &
    namefile, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read the chemical system according xml file. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout):: this ! Chemical system object 

character(len=*), intent(in)          :: namefile ! Name of.xml file

logical, intent(out)                  :: iserror ! If .true. there was error 
 
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
 iostat, &
 i, &
 j, &
 isp, &
 k, &
 isp2
type(xml_t):: &
 fxml
character(len=100)  :: &
 name, &
 msg
type(dictionary_t) :: &
 attributes 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------
call open_xmlfile(namefile, fxml, iostat)
if (iostat/=0) then
 msg='Error when open xml file:'
 call add_ (msg,namefile)
 goto 10
end if
!%---------------------------------------------------------
print *,'=======> Reading Chemical System object'
!%----------------------------------------------------------
!%-----------------------------------------------------Read
!%----------------------------------------------------------
call xml_parse &
 (fxml, &
  begin_element_handler=begin_element_handler)
!%----------------------------------------------------------
!% Read type of chemical system
!%----------------------------------------------------------
call read_xml_loc_ (name,attributes,this,iserror)
if (iserror) goto 10
!%----------------------------------------------------------
! End and close xml file
!%----------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%----------------------------------------------------------
!% Read parent chemical system
!%----------------------------------------------------------
call read_xml_ (this%pp,namefile, iserror)
if (iserror) goto 10
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
  allocate(this%pchemsysclassic)
  call create_ (this%pchemsysclassic)
  print *,'====> Set Chemical System (Classic)'
  call set_ (this%pchemsysclassic,this%pp,iserror)
 
case (chemsysretraso)
 
  allocate(this%pchemsysretraso)
  call create_ (this%pchemsysretraso)
  print *,'====> Set Chemical System (Retraso)'
  call set_ (this%pchemsysretraso,this%pp,iserror)
 
case default
 
  msg='Error, not implemented chemical system type'
  goto 10
 
end select
 
if (iserror) then
 msg='Error when calling set_'
 goto 10
end if
!%----------------------------------------------------------
print *,'=======> Reading Chemical System object finished' 
!%----------------------------------------------------------
!% Write the chemical system
!%----------------------------------------------------------
if (this%pp%iswriteinfo) then
 open(unit=999,file='chemical_system_info.dat',status='unknown')
 call write_ (this,999,iserror)
 if (iserror) goto 10  
 close(999)
end if
!%----------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Chemical System:'
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
subroutine read_txt_pchemsys &
   (this, &
    namefile, &
	namethdb, &
	namekindb, &
	filebase, &
	itype, &
	ioptxhl, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read the chemical system from txt files. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout) :: this

character(len=*), intent(in)           :: namefile 

character(len=*), intent(in)           :: namethdb

character(len=*), intent(in)           :: namekindb

character(len=*), intent(in)           :: filebase

integer, intent(in)                    :: itype 

integer, intent(in)                    :: ioptxhl 

logical, intent(out)                   :: iserror 
 
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
iserror=.false.
msg=''
!%--------------------------------------------------------
print *,'=======> Reading Chemical System object'
!%----------------------------------------------------------
!% Read parent chemical system 
!%----------------------------------------------------------
call read_txt_ (this%pp,namefile,namethdb,namekindb,filebase, &
                ioptxhl,iserror)
if (iserror) goto 10
!%----------------------------------------------------------
!% Set Chemical System (child)
!%----------------------------------------------------------
this%itype=itype 
!%----------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
  allocate(this%pchemsysclassic)
  call create_ (this%pchemsysclassic)
  call set_ (this%pchemsysclassic,this%pp,iserror)
 
case (chemsysretraso)
 
  allocate(this%pchemsysretraso)
  call create_ (this%pchemsysretraso)
  call set_ (this%pchemsysretraso,this%pp,iserror)

case default  
 
  msg='Error, not defined specialization, itype='
  call add_ (msg,itype)
  goto 10 
 
end select
 
if (iserror) then
 msg='Error when calling set_'
 goto 10
end if
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
print *,'=======> Reading Chemical System object finished'
!%--------------------------------------------------------
return
 
10 continue 
print *,'******************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: read_txt_'
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
subroutine destroy_chemsys &
   (this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy the chemical system object
!
!   $Arguments:
!
 
type (t_chemicalsystem) :: this 
 
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
if (associated(this%pp)) then
 call destroy_ (this%pp)
 deallocate (this%pp)
end if
this%pp => null ()
!%-------------------------------------------------------------
if (associated(this%pchemsysclassic)) then
 call destroy_ (this%pchemsysclassic)
 deallocate (this%pchemsysclassic)
end if
this%pchemsysclassic => null ()
!%-------------------------------------------------------------
if (associated(this%pchemsysretraso)) then
 call destroy_ (this%pchemsysretraso)
 deallocate (this%pchemsysretraso)
end if
this%pchemsysretraso => null ()
!%------------------------------------------------------------
this%itype=0
!%------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine switch_base_chemsys &
   (this, &
    naqprisp, &
	nadsprisp, &
    nameaqprisp, &
	nameadsprisp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Switch the set of primary species in the chemical system 
! object. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout), target         :: this           ! Type parent chemical system variable

integer, intent(in)                                    :: naqprisp       ! Number of aqueous primary species 

integer, intent(in)                                    :: nadsprisp      ! Number of adsorption primary species 

character(len=*), intent(in), dimension(naqprisp)      :: nameaqprisp    ! New set of aqueous primary species  

character(len=*), intent(in), dimension(nadsprisp)     :: nameadsprisp   ! New set of adsorption primary species 

logical, intent(out)                                   :: iserror        ! iserror=true, then there was an error 
 
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
 
!%---------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Call the corresponding service in the specialization 
!%------------------------------------------------------------
select case (this%itype)
 
case (chemsysretraso)
 
  call switch_base_ &
   (this%pchemsysretraso, &
    naqprisp, &
	nadsprisp, &
    nameaqprisp, &
	nameadsprisp, &
    iserror)
 
 
case default
 
  call switch_base_ &
   (this%pp, &
    naqprisp, &
	nadsprisp, &
    nameaqprisp, &
	nameadsprisp, &
    iserror)
 
end select
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: set_chemical_base_'
print *,'********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_iswriteinfo_chemsys &
   (this, &
    iswriteinfo, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set if the chemical system object must write the 
! Newton-Raphson information during the speciation processes. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout)       :: this

logical, intent(in)                          :: iswriteinfo

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
!-------------------------------------------------------------------------
!
!   $code
call set_iswriteinfo_ (this%pp,iswriteinfo,iserror) 
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_omgwcryst_chemsys &
   (this, &
    omgwcryst, &
    omgwfree, &
    c, &
    nsp, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the mass of water fixed in the hydrated minerals. 
!omgwcryst= PMh2o * 10^-3 * sum(coeff_j_h2o * cm * omgwfree)
!%
!%   where
!%
!%   PMh2o = Molecular weight of water
!%   coeff_j_h2o = stoichiometric coefficient of water in the
!%                 jth mineral
!%
!%   ---------------
!%   nsp    = number of species
!%   c(nsp) = vector of concentrations
!%   omgwfree = mass of water free per unit of volume
!%
!%   ----------------
!%   omgwcryst = mass of water per unit of volume contained in
!%            hydrated minerals
!%   iserror = .true. then there was one error
!%   msg=    message error
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent (in)      :: this

integer, intent(in)                       :: nsp

real*8, intent(out)                       :: omgwcryst

real*8, intent(in)                        :: omgwfree

real*8, intent(in), dimension(nsp)        :: c 

logical, intent(out)                      :: iserror

character(len=*), intent(out)             :: msg 
 
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
!%----------------------------------------------------------
call compute_omgwcryst_ &
   (this%pp, &
    omgwcryst, &
    omgwfree, &
    c, &
    nsp, &
    msg, &
    iserror)
!%-------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine specia_eqmin_from_setre_chemsys &
   (this, &
    temp, &
    ck1, &
    g, &
    ck, &
    iscompzchanged, &
    setre, &
    nsp, &
    dtime, &
    hashcompz, &
    omgwfreek1, &
    omgwfreek, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Specia according Setre=b 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout)  :: this

real*8, intent(in)                      :: temp 

real*8, intent(in)                      :: dtime

integer, intent(in)                     :: nsp

integer, intent(inout)                  :: hashcompz

real*8, intent(inout), dimension(nsp)   :: ck1                 ! Concentration vector in k+1

real*8, intent(inout), dimension(nsp)   :: setre

real*8, intent(in), dimension(nsp)      :: g                   ! Activity coefficients vector 

real*8, intent(in), dimension(nsp)      :: ck                  ! Concentration vector in k 

logical, intent(out)                    :: iserror

real*8, intent(in)                      :: omgwfreek1          ! Mass of free water in k+1

real*8, intent(in)                      :: omgwfreek           ! Mass of free water in k
 
logical, intent(out)                    :: iscompzchanged 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)          :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
select case (this%itype)
 
case (chemsysretraso)
 
  call specia_eqmin_from_setre_ &
   (this%pchemsysretraso, &
    ck1, &
    g, &
    ck, &
    iscompzchanged, &
    setre, &
    nsp, &
    dtime, &
    hashcompz, &
    omgwfreek1, &
    omgwfreek, &
    iserror)
 
case default
 
 msg='Error, public service not implemented'
 goto 10
 
end select
!%-------------------------------------------------------------
return
 
10 continue 
print *,'*********************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
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
subroutine update_temp_param_chemsys &
   (this, &
    temp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update the chemical system obejct to some temperature. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout) :: this    ! Type chemical system variable 

real*8, intent(in)                     :: temp    ! Temperature in celsius

logical, intent(out)                   :: iserror ! iserror=true, then there was an error 
 
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
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% If the reference temperature is equal to temp, then return 
!%-----------------------------------------------------------
if (this%pp%tempref==temp) return
!%-----------------------------------------------------------
!% Assign the new temperature 
!%-----------------------------------------------------------
this%pp%tempref=temp 
!%-----------------------------------------------------------
!% Update the depending temperature parameters in the parent 
!% chemical system 
!%-----------------------------------------------------------
call update_ (this%pp,this%pp%tempref,iserror)
if (iserror) goto 10
!%-----------------------------------------------------------
return
10 continue 
print *,'******************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: update_ (temperature)'
print *, msg
 print *,'******************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_porosity_chemsys &
   (this, &
    porosity, &
    ck1, &
    ck, &
    nsp, &
    omgwfreek1, &
    omgwfreek, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update the porosity from the mineral concentrations.
! por = porosity (in/out)
! 
!%    c = vector of concentrations in k+1
!%    cold = vector of concentrations in k
!%    nsp = number of species
!%    omgwfreek = mass of water free per unit of volume in k
!%    omgwfreek1 = mass of water free per unit of volume in k + 1
!%    iserror = .true. there was error
!%
!%   delta(m3min)= sum( (cminj(k+1)*omgwfree(k+1) - cminj(k+1)*omgwfree(k)
!%   (1-por(k+1)) - (1-por(k)) = delta(m3min)
!%    por(k+1) = por(k) - delta(m3min)
!%
!%   where
!%   cminj = concentration of jth mineral
!%   VMj = Molar volume of jth mineral
!%   por = porosity 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in):: this       ! Chemical system object 

integer, intent(in)                :: nsp        ! Number of species

real*8, intent(inout)              :: porosity   ! Porosity 

real*8, intent(in)                 :: omgwfreek  ! Mass of water in k 

real*8, intent(in)                 :: omgwfreek1 ! Mass of water in k+1

real*8, intent(in), dimension(nsp) :: ck1        ! Concentrations in k+1

real*8, intent(in), dimension(nsp) :: ck         ! Concentraions in k

logical, intent(out)               :: iserror    ! iserror=true, then there was an error. 
 
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

msg=''
iserror=.false.
!%------------------------------------------------------------
call update_(this%pp,porosity,ck1,ck,nsp,omgwfreek1,omgwfreek,iserror)
if (iserror) goto 10
!%------------------------------------------------------------
return
10 continue 
print *,'******************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: update_ (porosity)'
print *, msg
 print *,'******************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_chem_info_chemsys &
   (this, &
    iserror, &
    nummobph, &
    numbase, &
    numsp, &
	numreact, &
	numeqreact, &
	numkinreact, &
    numsites, &
    wcindex, &
    ecindex, &
    numph, &
    wspindex, &
    numsurf, &
    namesp, &
	nameminsp, &
    nameph, &
    namesurf, &
    idbase, &
    hashcompz, &
    ideqminsp, &
    neqminsp, &
    namebase, &
    nminsp, &
    idminsp, &
    ngassp, &
    idgassp,&
    nummobsp)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return general chemical information about the 
! chemical system.
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)         :: this

logical, intent(out)                        :: iserror

integer, intent(out), optional              :: numsp        ! Number of species

integer, intent(out), optional             :: nummobph      ! Number of mobile phases. 

integer, intent(out), optional             :: nummobsp      ! Number of mobile species

integer, intent(out), optional              :: numbase

integer, intent(out), optional              :: numsites

integer, intent(out), optional              :: numreact       ! Number of reactions

integer, intent(out), optional              :: numeqreact     ! Number of equilibrium reactions 

integer, intent(out), optional              :: numkinreact    ! Number of kinetic reactions 

integer, intent(out), optional              :: wcindex       

integer, intent(out), optional              :: ecindex

integer, intent(out), optional              :: numph

integer, intent(out), optional              :: numsurf

integer, intent(out), optional              :: wspindex

integer, intent(out), optional              :: neqminsp

integer, intent(out), optional              :: nminsp

integer, intent(out), optional              :: ngassp

integer, intent(in), optional               :: hashcompz       ! Hash index of the components zone

integer, pointer, optional                  :: idbase(:)

integer, pointer, optional                  :: ideqminsp(:)   !Index of mineral species in equilibrium

integer, pointer, optional                  :: idminsp(:)

integer, pointer, optional                  :: idgassp(:)

character(len=100), pointer, optional       :: namesp(:)

character(len=100), pointer, optional       :: nameminsp(:)

character(len=100), pointer, optional       :: nameph(:)

character(len=100), pointer, optional       :: namesurf(:)

character(len=100), pointer, optional       :: namebase(:)
  
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                              :: &
 msg
logical                                         :: &
 havehashcompz 
!-------------------------------------------------------------------------
!
!   $code
!

 

!%-------------------------------------------------------------
iserror=.false.
msg=' '
!%-------------------------------------------------------------
!% Check optional arguments 
!%-------------------------------------------------------------
havehashcompz=present(hashcompz)
!%-------------------------------------------------------------
call get_chem_info_ &
   (this%pp, &
    msg, &
    iserror, &
    nummobph=nummobph, &
    numsp=numsp, &
	ngassp=ngassp, &
    numsites=numsites, &
    ecindex=ecindex, &
    numph=numph, &
    wspindex=wspindex, &
    numsurf=numsurf, &
	numreact=numreact, &
	numeqreact=numeqreact, &
	numkinreact=numkinreact, &
    namesp=namesp, &
	nameminsp=nameminsp, &
    nameph=nameph, &
    namesurf=namesurf, &
    ideqminsp=ideqminsp, &
    neqminsp=neqminsp, &
    nminsp=nminsp, &
    idgassp=idgassp, &
    idminsp=idminsp,&
    nummobsp=nummobsp)
!%-------------------------------------------------------------
!% Call the corrresponding service in the specialization 
!%-------------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
  call get_chem_info_ &
   (this%pp, &
    msg, &
    iserror, &
    numbase=numbase, &
    wcindex=wcindex, &
    idbase=idbase, &
    namebase=namebase)
 
case (chemsysretraso)
 
 if (havehashcompz) then
 
  call get_chem_info_ &
      (this%pchemsysretraso, &
       hashcompz, &
       msg, &
       iserror, &
       numbase=numbase, &
       idbase=idbase, &
       namebase=namebase, &
       wcindex=wcindex)
 
 else
 
  call get_chem_info_ &
   (this%pp, &
    msg, &
    iserror, &
    numbase=numbase, &
    wcindex=wcindex, &
    idbase=idbase, &
    namebase=namebase)
 
 end if
 
end select
!%-------------------------------------------------------------
if (iserror) goto 10
!%-------------------------------------------------------------
return
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: get_chem_info_'
print *,msg
print *,'**************************'
iserror=.true.
return
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_iumob_chemsys &
   (this, &
    iumob, &
    naqcol, &
    ngascol, &
    nnonaqcol, &
    c, &
    nsp, &
    ithcomp, &
    hashcompz, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the icomp aqueous components in the mobile phases 
! (aqueous, gas, and non-aqueous phases
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)  :: this

real*8, pointer, dimension(:)        :: iumob 

integer, intent(in)                  :: nsp

integer, intent(in)                  :: ithcomp

integer, intent(in)                  :: hashcompz

integer, intent(out)                 :: naqcol

integer, intent(out)                 :: ngascol

integer, intent(out)                 :: nnonaqcol

real*8, intent(in), dimension(nsp)   :: c 

logical, intent(out)                 :: iserror 
 
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
 
case (chemsysclassic)
 
   call compute_iumob_ &
   (this%pchemsysclassic, &
    iumob, &
    naqcol, &
    ngascol, &
    nnonaqcol, &
    c, &
    nsp, &
    ithcomp, &
    iserror)
 
 
 
case (chemsysretraso)
 
  call compute_iumob_ &
   (this%pchemsysretraso, &
    iumob, &
    naqcol, &
    ngascol, &
    nnonaqcol, &
    c, &
    nsp, &
    ithcomp, &
	hashcompz, &
    iserror)
 
 
end select
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine select_cmob_chemsys &
   (this, &
    cmob, &
    nummobph, &
    numsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the icomp aqueous components in the mobile phases (aqueous, gas, and non-aqueous phases
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)   :: this

integer, intent(out)                  :: nummobph

integer, intent(in)                   :: numsp

logical, intent(out)                  :: iserror

real*8, intent(in), dimension(numsp)  :: c 

real*8, pointer, dimension(:,:)       :: cmob 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call select_cmob_(this%pp,cmob,nummobph,numsp,c,iserror)
!%------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: select_cmob_'
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
subroutine select_caq_chemsys &
   (this, &
    caq, &
    nsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)        :: this

integer, intent(in)                        :: nsp

logical, intent(out)                       :: iserror

real*8, intent(in), dimension(nsp)         :: c 

real*8, pointer, dimension(:)              :: caq
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)             :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!


!%----------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call select_caq_(this%pp,caq,nsp,c,iserror)
!%------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: select_caq_'
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
subroutine select_cgas_chemsys &
   (this, &
    cgas, &
    nsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)        :: this

integer, intent(in)                        :: nsp

logical, intent(out)                       :: iserror

real*8, intent(in), dimension(nsp)         :: c 

real*8, pointer, dimension(:)              :: cgas
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)             :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!


!%----------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call select_cgas_(this%pp,cgas,nsp,c,iserror)
!%------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: select_cgas_'
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
subroutine select_cads_chemsys &
   (this, &
    cads, &
    nsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)        :: this

integer, intent(in)                        :: nsp

logical, intent(out)                       :: iserror

real*8, intent(in), dimension(nsp)         :: c 

real*8, pointer, dimension(:)              :: cads 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=300)             :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call select_cads_(this%pp,cads,nsp,c,iserror)
!%------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: select_cads_'
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
subroutine select_cmin_chemsys &
   (this, &
    cmin, &
    nsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)        :: this

integer, intent(in)                        :: nsp

logical, intent(out)                       :: iserror

real*8, intent(in), dimension(nsp)         :: c 

real*8, pointer, dimension(:)              :: cmin
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)             :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!


!%----------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call select_cmin_(this%pp,cmin,nsp,c,iserror)
!%------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: select_cmin_'
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
subroutine compute_dcmob1_chemsys &
   (this, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    nsp, &
    hashcompz, &
    iserror, &
    dg, &
    g)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)               :: this

integer, intent(in)                               :: nsp

integer, intent(in)                               :: hashcompz

real*8, intent(in)                                :: c(nsp)

real*8, pointer, dimension(:,:)                   :: dcmob

integer, intent(out)                              :: nrow

integer, intent(out)                              :: ncol

integer, intent(out)                              :: nmobph

logical, intent(out)                              :: iserror

real*8, pointer, optional                         :: dg(:,:)

real*8, pointer, optional                         :: g(:) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 character(len=100)           :: &
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
 
case (chemsysclassic)
 
   call compute_dcmob_ &
   (this%pchemsysclassic, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    nsp, &
    iserror, &
    dg, &
    g)
 
 
 
case (chemsysretraso)
 
  call compute_dcmob_ &
   (this%pchemsysretraso, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    nsp, &
    hashcompz, &
    iserror, &
    dg, &
    g)
 
case default
 
   msg='Error, not implemented public service'
   goto 10
 
end select
!%------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dcmob_'
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
subroutine compute_dcmob2_chemsys &
   (this, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    nsp, &
    naqpri, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem)               :: &
 this
integer                               :: &
 nsp, &
 naqpri
real*8                                :: &
 dc(nsp,naqpri)
real*8, pointer                       :: &
 dcmob(:,:)
integer, intent(out)                  :: &
 nrow, &
 ncol, &
 nmobph
logical                               :: &
 iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                    :: &
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
 
case (chemsysclassic)
 
 call compute_dcmob_ &
   (this%pchemsysclassic, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    nsp, &
    naqpri, &
    iserror)
 
 
 
case (chemsysretraso)
 
 call compute_dcmob_ &
   (this%pchemsysretraso, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    nsp, &
    naqpri, &
    iserror)

case default
 
   msg='Error, not implemented public service'
   goto 10
 
end select
!%------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dcmob_'
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
subroutine compute_dcads_chemsys &
   (this, &
    dcads, &
    nrow, &
    ncol, &
    dc, &
    nsp, &
    naqpri, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute dcads/dc1
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)        :: this

integer, intent(in)                        :: nsp

integer, intent(in)                        :: naqpri

real*8, intent(in)                         :: dc(nsp,naqpri)

real*8, pointer                            :: dcads(:,:)

integer, intent(out)                       :: nrow

integer, intent(out)                       :: ncol

logical, intent(out)                       :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
 call compute_dcads_ &
   (this%pchemsysclassic, &
    dcads, &
    nrow, &
    ncol, &
    dc, &
    nsp, &
    naqpri, &
    iserror)
 
 
 
case (chemsysretraso)
 
 call compute_dcads_ &
   (this%pchemsysretraso, &
    dcads, &
    nrow, &
    ncol, &
    dc, &
    nsp, &
    naqpri, &
    iserror)
 
case default
 
   msg='Error, not implemented public service'
   goto 10
 
end select
!%------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dcads_'
print *, msg
print *,'*************************'
iserror=.true.
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_aq_sps_pos_chemsys &
   (this, &
    isps1, &
    isps2, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the first and last indice of the species 
! corresponding to the aqueous phase.
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in) :: this

integer, intent(out)                :: isps1

integer, intent(out)                :: isps2

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
!% Call the corresponding service in the parent chemical system
!%------------------------------------------------------------
call get_aq_sps_pos_ (this%pp,isps1,isps2,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public Subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_gas_sp_index_chemsys &
   (this, &
    iserror, &
    idgassp, &
    ngassp)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return global index corresponding to gas species and the total number of gas species.
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)         :: this

integer, pointer, optional, dimension(:)          :: idgassp

integer, intent(out), optional                    :: ngassp

logical, intent(out)                              :: iserror 
 
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
!% Call the corresponding service in the parent chemical system
!%------------------------------------------------------------
call get_gas_sp_index_(this%pp,iserror,idgassp,ngassp)
 !%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_new_chemical_base_chemsys &
   (this, &
    nameaqprisp, &
	nameadsprisp, &
    c, &
	nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return a new set of chemical base determined from the 
! concentration vector. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)                    :: this           ! Type chemical system variable

integer, intent(in)                                    :: nsp            ! Number of species. 

real*8, intent(in), dimension(nsp)                     :: c              ! Concentration vector.  

character(len=*), pointer, dimension(:)                :: nameaqprisp    ! Name of new set of primary aqueous species.

character(len=*), pointer, dimension(:)                :: nameadsprisp   ! Name of new set of primary sorption species.

logical, intent(out)                                   :: iserror        ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                         :: &
 i, &
 isps1, &
 isps2
real*8                          :: &
 cmax
character(len=100)              :: &
 name, &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Call the corresponding service in the parent chemical 
!% system 
!%------------------------------------------------------------
call get_new_chemical_base_ &
   (this%pp, &
    nameaqprisp, &
	nameadsprisp, &
    c, &
	nsp, &
    iserror)
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: get_new_chemical_base_'
print *,'********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_chemsys &
   (this, &
    name, &
    tempisoterm, &
    phase, &
    surface, &
    reaction, &
    numph, &
    numsurf, &
    numreact, &
    itype, &
    components, &
    tolunknr, &
    tolresnr, &
    deltasatmin, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the chemical system object
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout)      :: this

integer, intent(in)                         :: numph

integer, intent(in)                         :: numsurf

integer, intent(in)                         :: numreact

type (t_phase), intent(in)                  :: phase(numph)

type (t_surface), intent(in)                :: surface(numsurf)

type (t_reaction)                           :: reaction(numreact)

character(len=*), intent(in)                :: name

character(len=*), intent(in)                :: components

integer, intent(in)                         :: itype

logical, intent(out)                        :: iserror

real*8, intent(in)                          :: tempisoterm

real*8, intent(in)                          :: tolunknr

real*8, intent(in)                          :: tolresnr

real*8, intent(in)                          :: deltasatmin
 
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
iserror=.false.
!%--------------------------------------------------------
call set_ &
   (this%pp, &
    name, &
    tempisoterm, &
    phase, &
    surface, &
    reaction, &
    numph, &
    numsurf, &
    numreact, &
    components, &
    tolunknr, &
    tolresnr, &
    deltasatmin, &
    iserror)
 
if (iserror) return
!%--------------------------------------------------------
! Dispatching
!%--------------------------------------------------------
this%itype=itype
!%--------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
   allocate(this%pchemsysclassic)
   call create_ (this%pchemsysclassic)
   call set_ (this%pchemsysclassic,this%pp,iserror)
 
case (chemsysretraso)
 
   allocate(this%pchemsysretraso)
   call create_ (this%pchemsysretraso)
   call set_ (this%pchemsysretraso,this%pp,iserror)
 
end select
 
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_iswcompbal_chemsys &
   (this, &
    iswcompbal, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set if the water component mass balance must be evaluated 
! in the speciation algorithms. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout)      :: this

logical, intent(in)                         :: iswcompbal 

logical, intent(out)                        :: iserror
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
!%--------------------------------------------------------
call set_iswcompbal_ (this%pp,iswcompbal,iserror)
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine specia_from_cpri_chemsys &
   (this, &
    temp, &
    c, &
    g, &
	cold, &
    alpha, &
    nsp, &
    cpri, &
    npri, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    ntxoh, &
    nsurf, &
    dtime, &
	volgas, &
    ionstr, &
    hashcompz, &
	faccap, &
    isanomalous, &
    isupmxitergam, &
    isreset, &
    iserror, &
    dc, &
    dg, &
    sktrk, &
    dsktrk)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Make the speciation from concentration of primary aqueous
! species and compute derivates and kinetic terms (optional)
!%
!%
!%
!%    c[nsp] = Concentration vector
!%    g[nsp] = Activity coefficients vector
!%    area[nsp] = factor for multiply species
!%    area0[nsp] = initial factor
!%    nsp  = number of species
!%    cpri = concentration of base
!%    npri = number of base species
!%    txoh[ntxoh,nsurf] = total sites for adsorption
!%    mxtxoh = maximun number of sites
!%    nsurf  = number of sufaces
!%    temp  = temperature in celsius
!%    dt  = time increment
!%    ionstr = ionic strength
!%    hashcompz  = hash index of components (now only used for chemical sys
!%    isanomalous = .true. if there is isanomalous concentrations
!%    reset = if .true. then init c = 0
!%    iserror = if .true. there was error
!%
 
!%    ------------------
!%    dc      = dc/dc1
!%    dg      = dg/dc1
!%    sktrk   = sktrk
!%    dsktrk  = d sktrk / dc1
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)                :: this

real*8, intent(in)                                 :: temp               ! Temperaure in celcius 

integer, intent(in)                                :: hashcompz          ! Components zone index 

real*8, intent(in)                                 :: faccap             ! Capillary correction for water activity 

integer, intent(in)                                :: nsp                ! Number of species. 

integer, intent(in)                                :: ntxoh    

integer, intent(in)                                :: nsurf

integer, intent(in)                                :: npri               ! Number of primary species.

real*8, intent(inout), dimension(nsp)              :: c                  ! Concentration vector 

real*8, intent(inout), dimension(nsp)              :: g                  ! Activity coefficients vector 

real*8, intent(in), dimension(nsp)                 :: cold

real*8, intent(in), dimension(nsp)                 :: alpha

real*8, intent(in), dimension(npri)                :: cpri

real*8, intent(in), dimension(ntxoh,nsurf)         :: txoh

real*8, intent(in), dimension(ntxoh,nsurf)         :: cap1

real*8, intent(in), dimension(ntxoh,nsurf)         :: cap2

real*8, intent(in), dimension(ntxoh,nsurf)         :: spsurfarea 

logical, intent(out)                               :: iserror

logical, intent(out)                               :: isanomalous

logical, intent(out)                               :: isupmxitergam

logical, intent(in)                                :: isreset

real*8, intent(in)                                 :: dtime

real*8, intent(in)                                 :: volgas 

real*8, intent(out)                                :: ionstr

real*8, pointer, optional, dimension(:,:)          :: dc
 
real*8, pointer, optional, dimension(:,:)          :: dg

real*8, pointer, optional, dimension(:)            :: sktrk

real*8, pointer, optional, dimension(:,:)          :: dsktrk
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)          :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%-----------------------------------------------------------
select case (this%itype)

 
case (chemsysretraso)
 
 call specia_from_cpri_ &
   (this%pchemsysretraso, &
    temp, &
    c, &
    g, &
	cold, &
    alpha, &
    nsp, &
    cpri, &
    npri, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &   
    ntxoh, &
    nsurf, &
    dtime, &
	volgas, &
    ionstr, &
    hashcompz, &
	faccap, &
    isanomalous, &
    isupmxitergam, &
    isreset, &
    iserror, &
    dc=dc, &
    dg=dg, &
    sktrk=sktrk, &
    dsktrk=dsktrk)
 
case default
 
 call specia_from_cpri_ &
   (this%pp, &
    temp, &
    c, &
    g, &
	cold, &
    alpha, &
    nsp, &
    cpri, &
    npri, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    ntxoh, &
    nsurf, &
    dtime, &
	volgas, &
    ionstr, &
	faccap, &
    isanomalous, &
    isupmxitergam, &
    isreset, &
    iserror, &
    dc=dc, &
    dg=dg, &
    sktrk=sktrk, &
    dsktrk=dsktrk)
 
end select
 
 
return
 
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: specia_from_cpri_'
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
subroutine specia_from_u_chemsys &
   (this, &
    temp, &
    c, &
    g, &
	iseqmin, &
	iseqgas, &
    ionstr, &
    alpha, &
    t, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    nsp, &
    npri, &
    ntxoh, &
    nsurf, &
    dtime, &
    hashcompz, &
    isconvergence, &
    factoromgw, &
	volgas, &
	pgas, &
	msg, &
    iserror, &
    cguess, &
    dc, &
    sktrk, &
    dsktrk, &
    dg, &
    simin, &
    nchemiter, &
    thetakin)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout)        :: this

real*8, intent(in)                            :: temp   ! Temperature [C]

real*8, intent(in)                            :: dtime  ! Time increment 

real*8, intent(out)                           :: ionstr ! Ionic strength 

real*8, intent(out)                           :: factoromgw 

real*8, intent(in)                            :: volgas

real*8, intent(in)                            :: pgas

logical, intent(in)                           :: iseqmin

logical, intent(in)                           :: iseqgas

integer, intent(in)                           :: nsp

integer, intent(inout)                        :: hashcompz

integer, intent(in)                           :: npri

integer, intent(in)                           :: ntxoh

integer, intent(in)                           :: nsurf

real*8, intent(inout), dimension(nsp)         :: c

real*8, intent(inout), dimension(nsp)         :: g

real*8, intent(in), dimension(nsp)            :: alpha

real*8, intent(in), dimension(npri)           :: t

real*8, intent(in), dimension(ntxoh,nsurf)    :: txoh

real*8, intent(in), dimension(ntxoh,nsurf)    :: c1

real*8, intent(in), dimension(ntxoh,nsurf)    :: c2

real*8, intent(in), dimension(ntxoh,nsurf)    :: spsurfarea

logical, intent(out)                          :: iserror

character(len=*), intent(out)                 :: msg

logical, intent(out)                          :: isconvergence 

real*8, pointer, dimension(:,:), optional     :: dc

real*8, pointer, dimension(:), optional       :: sktrk

real*8, pointer, dimension(:,:), optional     :: dsktrk

real*8, pointer, dimension(:,:), optional     :: dg

real*8, pointer, dimension(:), optional       :: simin

integer, intent(out), optional                :: nchemiter 

real*8, intent(in), dimension(npri), optional :: cguess

real*8, intent(in), optional                  :: thetakin 
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
msg=''
iserror=.false.
!%-----------------------------------------------------------
select case (this%itype)
 
 
case (chemsysretraso)
 
  call specia_from_u_ &
   (this%pchemsysretraso, &
    temp, &
    c, &
    g, &
	iseqmin, &
	iseqgas, &
    ionstr, &
    alpha, &
    t, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    nsp, &
    npri, &
    ntxoh, &
    nsurf, &
    dtime, &
    hashcompz, &
    isconvergence, &
    factoromgw, &
	volgas, &
	pgas, &
	msg, &
    iserror, &
    cguess=cguess, &
    dc=dc, &
    sktrk=sktrk, &
    dsktrk=dsktrk, &
    dg=dg, &
    simin=simin, &
    nchemiter=nchemiter, &
    thetakin=thetakin)
 
case default
 
 call specia_from_u_ &
   (this%pp, &
    temp, &
    c, &
    g, &
	iseqmin, & 
	iseqgas, & 
    ionstr, &
    alpha, &
    t, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    nsp, &
    npri, &
    ntxoh, &
    nsurf, &
    dtime, &
    isconvergence, &
    factoromgw, &
	volgas, &
	pgas, &
	msg, &
    iserror, &
    cguess=cguess, &
    dc=dc, &
    sktrk=sktrk, &
    dsktrk=dsktrk, &
    dg=dg, &
    simin=simin, &
    nchemiter=nchemiter, &
    thetakin=thetakin)

end select
 
return
!%-----------------------------------------------------------
10  continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: specia_from_u_'
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
subroutine specia_from_solution_type_chemsys &
   (this, &
    temp, &
    c, &
    g, &
    alpha, &
    numsp, &
    nummin, &
    simin, &
    namemin, &
    numcomp, &
    namecomp, &
    icon, &
    ctot, &
    cguess, &
    constraint, &
    ionstr, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    numsites, &
    numsurf, &
    isconvergence, &
    hashcompz, &
    iserror, &
    dc, &
	dg, &
    sktrk, &
    dsktrk, &
	cputime, &
    nchemiter)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout)     :: this

real*8, intent(in)                         :: temp                      ! Temperature [C]

integer, intent(in)                        :: numsites

integer, intent(in)                        :: numcomp
 
integer, intent(in)                        :: numsp

integer, intent(in)                        :: numsurf

integer, intent(out)                       :: nummin

integer, intent(inout)                     :: hashcompz

integer, intent(in)                        :: icon(numcomp)

real*8, intent(out)                        :: c(numsp)

real*8, intent(out)                        :: g(numsp)

real*8, intent(in)                         :: alpha(numsp)

real*8, intent(in)                         :: txoh(numsites,numsurf)

real*8, intent(in)                         :: c1(numsites,numsurf)

real*8, intent(in)                         :: c2(numsites,numsurf)

real*8, intent(in)                         :: spsurfarea(numsites,numsurf)

real*8, intent(in)                         :: cguess(numcomp)

real*8, intent(in)                         :: ctot(numcomp)

real*8, pointer                            :: simin(:)

character(len=*), intent(in)               :: namecomp(numcomp)

character(len=*), intent(in)               :: constraint(numcomp)

character(len=100), pointer                :: namemin(:)

real*8, intent(out)                        :: ionstr

real*8, pointer, optional                  :: dc(:,:)

real*8, pointer, optional                  :: dg(:,:)

real*8, pointer, optional                  :: sktrk(:)

real*8, pointer, optional                  :: dsktrk(:,:)

logical, intent(out)                       :: isconvergence 

logical, intent(out)                       :: iserror 

integer, intent(out), optional             :: nchemiter

real*8, intent(out), optional              :: cputime      ! CPU time consumed during speciation 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)          :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------
msg=''
iserror=.false.
!%-----------------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
 call specia_from_solution_type_ &
   (this%pp, &
    c, &
    g, &
    alpha, &
    numsp, &
    nummin, &
    simin, &
    namemin, &
    numcomp, &
    namecomp, &
    icon, &
    ctot, &
    cguess, &
    constraint, &
    ionstr, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    numsites, &
    numsurf, &
    isconvergence, &
    iserror, &
    dc=dc, &
	dg=dg, & 
    sktrk=sktrk, &
    dsktrk=dsktrk, &
	cputime=cputime, &
    nchemiter=nchemiter)
 
 
case (chemsysretraso)
 
 call specia_from_solution_type_ &
   (this%pchemsysretraso, &
    c, &
    g, &
    alpha, &
    numsp, &
    nummin, &
    simin, &
    namemin, &
    numcomp, &
    namecomp, &
    icon, &
    ctot, &
    cguess, &
    constraint, &
    ionstr, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    numsites, &
    numsurf, &
    isconvergence, &
    hashcompz, &
    iserror, &
    dc=dc, &
	dg=dg, &
    sktrk=sktrk, &
    dsktrk=dsktrk, &
	cputime=cputime, &
    nchemiter=nchemiter)
 
 
case default
 msg='Error, not implemented public service'
 goto 10
end select
!%-----------------------------------------------------------
return
 
10 continue 
print *,'*********************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
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
subroutine add_and_specia_chemsys &
   (this, &
    temp, &
    c, &
    g, &
    ionstr, &
    simin, &
    namespadd, &
    moladd, &
    mh2o, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    nsp, &
    nsurf, &
    ntxoh, &
    omgwfree, &
    isconvergence, &
    iserror, &
    uadd, &
    dc, &
    sktrk, &
    dsktrk, &
    dg, &
    nchemiter)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Add mol of components (uad) or mass of water (mh2o) and make the speciation. 
! The equilibrium with other phases is not imposed.
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)      :: this

real*8, intent(in)                       :: temp  ! Temperature [C]

real*8, intent (in)                      :: mh2o  ! mass of water per unit of volume evaporated/adde
 
real*8, intent (in)                      :: moladd

real*8, intent(inout)                    :: omgwfree

real*8, intent (out)                     :: ionstr             ! ionic strength

integer, intent(in)                      :: nsp

integer, intent(in)                      :: nsurf

integer, intent(in)                      :: ntxoh

character(len=*)                         :: namespadd      ! Name of species added

real*8, intent(inout)                    :: c(nsp) ! concentration vector

real*8, intent(inout)                    :: g(nsp)         ! activity coefficients vector

real*8, intent(out)                      :: simin(nsp)

real*8, intent(in)                       :: txoh(ntxoh,nsurf)

real*8, intent(in)                       :: cap1(ntxoh,nsurf)

real*8, intent(in)                       :: cap2(ntxoh,nsurf)

real*8, intent(in)                       :: spsurfarea(ntxoh,nsurf)

logical, intent(out)                     :: iserror

logical, intent(out)                     :: isconvergence

real*8, intent(in), optional             :: uadd(this%pp%numaqprisp)

real*8, pointer, optional                :: dc(:,:)

real*8, pointer, optional                :: sktrk(:)

real*8, pointer, optional                :: dsktrk(:,:)

real*8, pointer, optional                :: dg(:,:) 

integer, intent(out), optional           :: nchemiter 
 
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
!%-----------------------------------------------------------------
call add_and_specia_ &
   (this%pp, &
    c, &
    g, &
    ionstr, &
    simin, &
    namespadd, &
    moladd, &
    mh2o, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    nsp, &
    nsurf, &
    ntxoh, &
    omgwfree, &
    isconvergence, &
    iserror, &
    uadd=uadd, &
    dc=dc, &
    sktrk=sktrk, &
    dsktrk=dsktrk, &
    dg=dg, &
    nchemiter=nchemiter)
!%------------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: add_and_specia_'
print *, msg
print *,'*************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine equilibrate_chemsys &
   (this, &
    temp, &
    c, &
    g, &
    si, &
    ionstr, &
    nameph, &
    siph, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    nsp, &
    nsurf, &
    ntxoh, &
    omgwfree, &
    isconvergence, &
    iserror, &
    dc, &
    sktrk, &
    dsktrk, &
    dg, &
    nchemiter)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Make the chemical speciation equilibrating wiht nameph.
!
!   $Arguments:
!
 
 
type (t_chemicalsystem), intent(in)       :: this

real*8, intent(in)                        :: temp     ! Temperature [C]

real*8, intent(inout)                     :: omgwfree ! Mass of free water 

real*8, intent (out)                      :: ionstr ! Ionic strength 

real*8, intent (in)                       :: siph

integer, intent(in)                       :: nsp ! Number of species 

integer, intent(in)                       :: nsurf ! Number of surfaces

integer, intent(in)                       :: ntxoh

character(len=*), intent(in)              :: nameph         ! Name of species added

real*8, intent(inout), dimension(nsp)     :: c ! concentration vector

real*8, intent(inout), dimension(nsp)     :: g  ! activity coefficients vector

real*8, intent(out), dimension(nsp)       :: si   ! saturation indices of mineral species
 
real*8, intent(in), dimension(ntxoh,nsurf):: txoh

real*8, intent(in), dimension(ntxoh,nsurf):: cap1

real*8, intent(in), dimension(ntxoh,nsurf):: cap2

real*8, intent(in), dimension(ntxoh,nsurf):: spsurfarea

logical, intent(out)                      :: iserror

logical, intent(out)                      :: isconvergence

real*8, pointer, optional, dimension(:,:) :: dc
 
real*8, pointer, optional, dimension(:)   :: sktrk

real*8, pointer, optional, dimension(:,:) :: dsktrk

real*8, pointer, optional, dimension(:,:) :: dg 

integer, intent(out), optional            :: nchemiter 
 
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
!%-----------------------------------------------------------------
msg='' 
!%-----------------------------------------------------------------
call equilibrate_ &
   (this%pp, &
    c, &
    g, &
    si, &
    ionstr, &
    nameph, &
    siph, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    nsp, &
    nsurf, &
    ntxoh, &
    omgwfree, &
    isconvergence, &
    iserror, &
    dc=dc, &
    sktrk=sktrk, &
    dsktrk=dsktrk, &
    dg=dg, &
    nchemiter=nchemiter)

return
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: equilibrate_'
print *, msg
print *,'*************************'
iserror=.true.
return
!%-----------------------------------------------------------------
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_hashcompz_chemsys &
   (this, &
    hashcompz, &
    c, &
    g, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the hash index to identify the components definition. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(inout)           :: this

integer, intent(out)                             :: hashcompz

integer, intent(in)                              :: nsp        ! Number of species

real*8, intent(in), dimension(nsp)               :: c          ! Concentration vector [nsp]

real*8, intent(in), dimension(nsp)               :: g

logical, intent(out)                             :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                   :: &
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

case (chemsysretraso)
 
 call get_hashcompz_ &
  (this%pchemsysretraso, &
   hashcompz, &
   c, &
   g, &
   nsp, &
   iserror)
 
case default
 
  hashcompz=0
 
end select
!%------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: get_hashcompz_'
print *, msg
print *,'*************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_chemsys &
   (this, &
    ioutput, &
	iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii all attributes ancapsulated in the 
! chemical system object. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)  :: this

integer, intent(in)                  :: ioutput

logical, intent(out)                 :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                   :: &
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
 
case (chemsysclassic)
 
 call write_ (this%pchemsysclassic,ioutput,iserror)
 
case (chemsysretraso)
 
 call write_ (this%pchemsysretraso,ioutput,iserror)

 
case default
 
 msg='Error, public service not implemented'
 goto 10
 
end select
!%------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: write_'
print *, msg
print *,'*************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_chemsys &
   (copied, &
    this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy a chemical system object in other chemical system
! object. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)  ::this

type (t_chemicalsystem), intent(out) ::copied 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical                            :: &
 iserror 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%-----------------------------------------------------------------------
copied%itype=this%itype
copied%pp=this%pp
!%----------------------------------------------------------------------- 
select case (this%itype)
 
case (chemsysclassic)
 
   allocate(copied%pchemsysclassic)
   call create_ (copied%pchemsysclassic)
   call set_parent_ (copied%pchemsysclassic,copied%pp,iserror)
 
case (chemsysretraso)
 
   allocate(copied%pchemsysretraso)
   call create_ (copied%pchemsysretraso)
   copied%pchemsysretraso = this%pchemsysretraso
   call set_parent_ (copied%pchemsysretraso,copied%pp,iserror)
 
end select
!%-----------------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_mass_salt_chemsys &
   (this, &
    mass, &
    mol, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the mass of salt [kgr] in the aqueous phase 
! from molality vector. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)       :: this

integer, intent(in)                       :: nsp

real*8, intent(out)                       :: mass

real*8, intent(in), dimension(nsp)        :: mol 

logical, intent(out)                      :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: 
!
!   $Post-cond:
!
!   $License: CSIC-UPC
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
iserror=.false.
!%------------------------------------------------------------
!% Call the corresponding service in the parent chemical 
!% system object 
!%------------------------------------------------------------
call compute_mass_salt_(this%pp,mass,mol,nsp,iserror)
!%------------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
 subroutine compute_dmassSalt_dc_chemsys &
   (this, &
    dc, &
    dmassdc, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute derivative of the mass of salt in the aqueous phase
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in) :: this

real*8,dimension(:,:),pointer             :: dc

real*8,dimension(:),pointer               :: dmassdc

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
 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
iserror=.false.
!%------------------------------------------------------------
!% Call the corresponding service in the parent chemical 
!% system object 
!%------------------------------------------------------------
call compute_dmassSalt_dc_(this%pp,dc,dmassdc,iserror)
!%------------------------------------------------------------
return
 
end subroutine

!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_new_chemsys_chemsys &
   (this, &
    pnewchemsys, &
    nameph, &
    namesurf, &
    numph, &
    numsurf, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return a new chemical system from the name of phases and 
! surfaces.  
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)   :: this

type (t_chemicalsystem), pointer      :: pnewchemsys

integer, intent(in)                   :: numph

integer, intent(in)                   :: numsurf

character(len=*), intent(in)          :: nameph(numph)

character(len=*), intent(in)          :: namesurf(numsurf)

logical, intent(out)                  :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------
if (associated(pnewchemsys)) nullify (pnewchemsys)
allocate (pnewchemsys)
call create_ (pnewchemsys)
pnewchemsys%itype=this%itype
!%--------------------------------------------------------
call get_new_chemsys_ &
   (this%pp, &
    pnewchemsys%pp, &
    nameph, &
    namesurf, &
    numph, &
    numsurf, &
    iserror)
if (iserror) then
 msg='Error when call get_new_chemsys_'
 goto 10
end if
!%--------------------------------------------------------
select case (pnewchemsys%itype)
 
case (chemsysclassic)
 
   allocate(pnewchemsys%pchemsysclassic)
   call create_ (pnewchemsys%pchemsysclassic)
   call set_(pnewchemsys%pchemsysclassic,pnewchemsys%pp,iserror)
 
case (chemsysretraso)
 
   allocate(pnewchemsys%pchemsysretraso)
   call create_ (pnewchemsys%pchemsysretraso)
   call set_(pnewchemsys%pchemsysretraso,this%pp,iserror)
 
case default
 
   msg='Error, not recognized chemical system type'
   call add_ (msg,pnewchemsys%itype)
   goto 10
 
end select
!%-------------------------------------------------------- 
if (iserror) then
 msg='Error when calling set_'
 goto 10
end if
!%------------------------------------------------------------
return
 
10 continue 
print *, '*************************'
print *, 'Chemical System:'
print *, 'Name:',this%pp%name
print *, 'Service: get_new_chemsys_'
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
subroutine compute_usktrk1_chemsys &
   (this, &
    usktrk, &
    ncomp, &
    c, &
    g, &
    alpha, &
    nsp, &
    hashcompz, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute U*Skt*rk 
!
!   $Arguments:
!
type (t_chemicalsystem), intent(in) :: this ! Type chemical system variable. 

real*8, pointer, dimension(:)       :: usktrk

integer, intent(in)                 :: nsp ! Number of species. 

real*8, intent(in), dimension(nsp)  :: c  ! Concentration vector. 

real*8, intent(in), dimension(nsp)  :: g  ! Activity coefficients vector. 

real*8, intent(in), dimension(nsp)  :: alpha

integer, intent(in)                 :: hashcompz ! Index of components zone. 

integer, intent(out)                :: ncomp ! Number of components. 

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
 
!-------------------------------------------------------------------------
!
!   $code
!
 

 
character(len=100)           :: &
 msg
!%----------------------------------------------
msg=''
iserror=.false.
!%----------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
   call compute_usktrk_ &
   (this%pchemsysclassic, &
    usktrk, &
    ncomp, &
    c, &
    g, &
    alpha, &
    nsp, &
    iserror)
 
case (chemsysretraso)
 
  call compute_usktrk_ &
   (this%pchemsysretraso, &
    usktrk, &
    ncomp, &
    c, &
    g, &
    alpha, &
    nsp, &
    hashcompz, &
    iserror)

case default
 
 msg='Error, public service not implemented'
 goto 10
 
end select
!%------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_usktrk_'
print *, msg
print *,'*************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_usktrk2_chemsys &
   (this, &
    usktrk, &
    ncomp, &
    sktrk, &
    nsp, &
    hashcompz, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)  :: this

real*8, pointer, dimension(:)        :: usktrk 

integer, intent(in)                  :: nsp

real*8, intent(in), dimension(nsp)   :: sktrk 

integer, intent(in)                  :: hashcompz

integer, intent(out)                 :: ncomp

logical, intent(out)                 :: iserror
 
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
 
 

!%----------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
   call compute_usktrk_ &
   (this%pchemsysclassic, &
    usktrk, &
    ncomp, &
    sktrk, &
    nsp, &
    iserror)
 
case (chemsysretraso)
 
  call compute_usktrk_ &
   (this%pchemsysretraso, &
    usktrk, &
    ncomp, &
    sktrk, &
    nsp, &
    hashcompz, &
    iserror)
 
 
end select
!%------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dusktrk1_chemsys &
   (this, &
    dusktrk, &
    naqpri, &
    c, &
    alpha, &
    nsp, &
    hashcompz, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem)                :: this

real*8, pointer, dimension(:,:)        :: dusktrk

integer                                :: hashcompz 
 
integer                                :: naqpri ! Number of aqueous primary species.

integer                                :: nsp ! Number of species.

real*8, dimension(nsp)                 :: c

real*8, dimension(nsp)                 :: alpha

logical                                :: iserror  ! iserror=true, there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)           :: &
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
 
case (chemsysclassic)
 
  call compute_dusktrk_ &
   (this%pchemsysclassic, &
    dusktrk, &
    naqpri, &
    c, &
    alpha, &
    nsp, &
    iserror)
 
case (chemsysretraso)
 
 
    call compute_dusktrk_ &
   (this%pchemsysretraso, &
    dusktrk, &
    naqpri, &
    c, &
    alpha, &
    nsp, &
    hashcompz, &
    iserror)
 
 
end select
!%--------------------------------------------------------------
return
 
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dusktrk2_chemsys &
   (this, &
    dusktrk, &
    naqpri, &
    dsktrk, &
    numsp, &
    hashcompz, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem) :: &
 this
real*8, pointer              :: &
 dusktrk(:,:)
integer                      :: &
 hashcompz, &
 naqpri, &
 numsp
real*8                       :: &
 dsktrk(numsp,naqpri)
logical                      :: &
 iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)           :: &
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
 
case (chemsysclassic)
 
 call compute_dusktrk_ &
   (this%pchemsysclassic, &
    dusktrk, &
    naqpri, &
    dsktrk, &
    numsp, &
    iserror)
 
 
case (chemsysretraso)
 
 call compute_dusktrk_ &
   (this%pchemsysretraso, &
    dusktrk, &
    naqpri, &
    dsktrk, &
    numsp, &
    hashcompz, &
    iserror)

end select
!%--------------------------------------------------------------
return
 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_umob_chemsys &
   (this, &
    umob, &
    naqpri, &
    nmobph, &
    numsp, &
    c, &
    hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem)    :: &
 this
integer                               :: &
 nmobph, &
 naqpri, &
 numsp, &
 hashcompz
real*8                                :: &
 c(numsp)
real*8, pointer                       :: &
 umob(:,:)
logical                               :: &
 iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
 call compute_umob_ &
   (this%pchemsysclassic, &
    umob, &
    naqpri, &
    nmobph, &
    numsp, &
    c, &
    iserror)
 
 
case (chemsysretraso)
 
 call compute_umob_ &
   (this%pchemsysretraso, &
    umob, &
    naqpri, &
    nmobph, &
    numsp, &
    c, &
    hashcompz, &
    iserror)
 
case default
 msg='Error, not implemented public service'
 goto 10
end select
!%------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_umob_'
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
subroutine compute_ith_total_chemsys &
   (this, &
    ithcomp, &
    nametotal, &
    c, &
    numsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute total of ith aqueous components.
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)   ::  this

integer, intent(in)                   :: numsp

character(len=*), intent(in)          :: nametotal

real*8, intent(in)                    :: c(numsp)

real*8, intent(out)                   :: ithcomp

logical, intent(out)                  :: iserror 
 
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
call compute_ith_total_ &
   (this%pp, &
    ithcomp, &
    nametotal, &
    c, &
    numsp, &
    iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_alkalinity_chemsys &
   (this, &
    alk, &
    c, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute alkalinity from concentration vector. 
!
!   In general the alcalinity is defined as
!   alkalinity = 2*[co3-2] + [hco3-]
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)               :: this

integer, intent(in)                               :: nsp

real*8, intent(in), dimension(nsp)                :: c

real*8, intent(out)                               :: alk

logical, intent(out)                              :: iserror 
 
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
 
 
 

!%------------------------------------------------------
!% Call the corresponding service in the parent chemical
!% system object. 
!%------------------------------------------------------
call compute_alkalinity_ (this%pp,alk,c,nsp,iserror)
!%------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_inv_points_chemsys &
   (this, &
    awconstr, &
    isconstr, &
    idreact, &
    nreact, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Check and determine the set of reactions that constraint 
! the water activity. The implemented method is based on Risacher anr Clement (2001)
! If there are lineally dependent reactions that constraint the water 
! activity, then compute the water activity value. 
!
!   $Arguments:
!!
 
type (t_chemicalsystem), intent(in)       :: this

integer, intent(in)                       :: nreact   ! Number of reactions 

integer, intent(in), dimension(nreact)    :: idreact  ! Global index of reactions 

logical, intent(out), dimension(nreact)   :: isconstr ! Global index of contrained reactions 

real*8, intent(out)                       :: awconstr ! Constrained water activity 

logical, intent(out)                      :: iserror  ! If true there was error 

character(len=*), intent(out)             :: msg      ! Error message  
 
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
!% Call the corresponding service in the parent chemical 
!% system object. 
!%------------------------------------------------------------
call check_inv_points_ &
   (this%pp, &
    awconstr, & 
    isconstr, &
    idreact, &
    nreact, &
    msg, &
    iserror)
!%------------------------------------------------------------
!%------------------------------------------------------------
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
subroutine update_min_area_chemsys &
   (this, &
    area, &
    area0, &
    c, &
    c0, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update the reactive surface of all minerals from its  
! concentration. 
!
!   $Arguments:
!
 
 
type (t_chemicalsystem), intent(in)             :: this

real*8, intent(inout), dimension(this%pp%numsp) :: area     ! Area [m2]. 

real*8, intent(in), dimension(this%pp%numsp)    :: area0    ! Initial area [m2].

real*8, intent(in), dimension(this%pp%numsp)    :: c        ! Concentration [mol].

real*8, intent(in), dimension(this%pp%numsp)    :: c0       ! Initial concentration [mol].

character(len=*), intent(out)                   :: msg      ! Massage error. 

logical, intent(out)                            :: iserror  ! iserror=true, then there was an error. 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond: Only are updated the reactive surface of the minerals 
!               defnined in kinetc. 
!
!   $License:
!
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
!% Call the corresponding service in the parent chemical 
!% system object. 
!%------------------------------------------------------------
call update_min_area_ &
   (this%pp, &
    area, &
    area0, &
    c, &
    c0, &
    msg, &
    iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_txoh_chemsys &
   (this, &
    txoh, &
    txoh0, &
	idtxohmin, &
    c, &
    c0, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update the reactive surface of kinetic minerals
!
!   $Arguments:
!
 
 
type (t_chemicalsystem), intent(in)                                :: this      ! Type chemical system variable. 

real*8, intent(inout), dimension(this%pp%numsites,this%pp%numsurf) :: txoh      ! Total sites 

real*8, intent(in), dimension(this%pp%numsites,this%pp%numsurf)    :: txoh0     ! Initial total sites

integer, intent(in), dimension(this%pp%numsurf)                    :: idtxohmin ! Initial total sites

real*8, intent(in), dimension(this%pp%numsp)                       :: c         ! Concentration 

real*8, intent(in), dimension(this%pp%numsp)                       :: c0        ! Initial concentration 
 
character(len=*), intent(out)                                      :: msg       ! Error message. 

logical, intent(out)                                               :: iserror   ! iserror=true, then there was an error. 
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
!% 
!%------------------------------------------------------------
call update_txoh_ (this%pp,txoh,txoh0,idtxohmin,c,c0,msg,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_uads_chemsys &
   (this, &
    uads, &
    npri, &
    c, &
    nsp, &
    hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)     :: this

integer, intent(in)                     :: nsp

integer, intent(out)                    :: npri

integer, intent(in)                     :: hashcompz

real*8, intent(in)                      :: c(nsp)

real*8, pointer                         :: uads(:)

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
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%---------------------------------------------------------
msg=''
iserror=.false.
!%---------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
 call compute_uads_ &
   (this%pchemsysclassic, &
    uads, &
    npri, &
    c, &
    nsp, &
    iserror)
 
 
case (chemsysretraso)
 
 call compute_uads_ &
   (this%pchemsysretraso, &
    uads, &
    npri, &
    c, &
    nsp, &
    hashcompz, &
    iserror)

 
case default
 
 msg='Error, not defined public service'
 goto 10
 
end select
 
 
!%------------------------------------------------------------
 
 
 
return
 
10 continue 
print *,'************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_uads_'
print *, msg
iserror=.true.
print *,'************************'
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dumob1_chemsys &
   (this, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    hashcompz, &
    iserror)
 
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
type (t_chemicalsystem)               :: &
 this
integer                               :: &
 numsp, &
 nmobph, &
 nrow, &
 ncol, &
 hashcompz
real*8                                :: &
 c(numsp)
real*8, pointer                       :: &
 dumob(:,:)
logical                               :: &
 iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%---------------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
 call compute_dumob_ &
   (this%pchemsysclassic, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    iserror)
 
 
case (chemsysretraso)
 
 call compute_dumob_ &
   (this%pchemsysretraso, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    hashcompz, &
    iserror)
 
case default
 
 msg='Error, not defined public service'
 goto 10
 
end select
!%---------------------------------------------------------------
 
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dumob_'
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
subroutine compute_dumob2_chemsys &
   (this, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute U*dcmob/dc1
!
!   $Arguments:
!
 
type (t_chemicalsystem)               :: &
 this
integer                               :: &
 numsp, &
 naqpri, &
 hashcompz
real*8                                :: &
 dc(numsp,naqpri)
real*8, pointer                       :: &
 dumob(:,:)
integer, intent(out)                  :: &
 nrow, &
 ncol, &
 nmobph
logical                               :: &
 iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                    :: &
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
 
case (chemsysclassic)
 
  call compute_dumob_ &
   (this%pchemsysclassic, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    iserror)
 
 
case (chemsysretraso)
 
  call compute_dumob_ &
   (this%pchemsysretraso, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    hashcompz, &
    iserror)
 
case default
 
 msg='Error, not defined public service'
 goto 10
 
end select
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dumob_'
print *, msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_duads_chemsys &
   (this, &
    duads, &
    nrow, &
    ncol, &
    numsp, &
    c, &
    hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute U*dcads/dc11ith
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)   :: this

integer, intent(in)                   :: numsp

real*8, intent(in)                    :: c(numsp)

real*8, pointer                       :: duads(:,:)

integer, intent(out)                  :: nrow

integer, intent(out)                  :: ncol

integer, intent(in)                   :: hashcompz

logical, intent(out)                  :: iserror  
 
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
 

 

!%---------------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
 call compute_duads_ &
   (this%pchemsysclassic, &
    duads, &
    nrow, &
    ncol, &
    numsp, &
    c, &
    iserror)
 
case (chemsysretraso)
 
 call compute_duads_ &
   (this%pchemsysretraso, &
    duads, &
    nrow, &
    ncol, &
    numsp, &
    c, &
    hashcompz, &
    iserror)
 
 
case default
 
 msg='Error, not defined public service'
 goto 10
 
end select
 
!%---------------------------------------------------------------
 
 
return
10 continue       
print *,'***************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_duads_'
print *,msg
print *,'***************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%****************Public subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine make_lin_trf_vector_chemsys &
   (this, &
    vnew, &
    vold, &
    hashcompz, &
    iserror, &
	isueq)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the matricial product vector2=U*vector1.
!
!   $Arguments:
!

type (t_chemicalsystem), intent(in)  :: this

real*8, pointer, dimension(:)        :: vnew

real*8, intent(in), dimension(:)     :: vold

integer, intent(in)                  :: hashcompz

logical, intent(out)                 :: iserror 

logical, intent(in), optional        :: isueq

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
!% Call the corresponding service in the 
!% specialization. 
!%-----------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
 call make_lin_trf_(this%pchemsysclassic,vnew,vold,iserror,isueq=isueq)
 
 
case (chemsysretraso)
 
 call make_lin_trf_(this%pchemsysretraso,vnew,vold,hashcompz,iserror,isueq=isueq)

 
case default
 
 msg='Error, not implemented public service'
 goto 10
 
end select
 
!%-----------------------------------------------------------
return
10 continue 
print *,'**********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: make_lin_trf_'
print *,msg
print *,'**********************'
return
end subroutine
!%************************************************************
!%****************Public subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine make_lin_trf_array_chemsys &
   (this, &
    anew, &
    aold, &
    hashcompz, &
    iserror, &
	isueq)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the matricial product array2=U*array1
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)   :: this

real*8, pointer, dimension(:,:)       :: anew

real*8, intent(in), dimension(:,:)    :: aold 

integer, intent(in)                   :: hashcompz

logical, intent(out)                  :: iserror  

logical, intent(in), optional         :: isueq
 
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
!%-----------------------------------------------------------
!% Call the corresponding service in the 
!% specialization. 
!%-----------------------------------------------------------
select case (this%itype)
 
case (chemsysclassic)
 
 call make_lin_trf_ &
   (this%pchemsysclassic, &
    anew, &
    aold, &
    iserror, &
	isueq=isueq)
 
 
case (chemsysretraso)
 
 call make_lin_trf_ &
   (this%pchemsysretraso, &
    anew, &
    aold, &
    hashcompz, &
    iserror, &
	isueq=isueq)
 
case default
 
 msg='Error, not implemented public service'
 goto 10
 
end select
 
 
!%------------------------------------------------------------
return
10 continue 
print *,'*****************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: make_lin_trf_'
print *, msg
iserror=.true.
print *,'*****************************'
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_min_sat_index_chemsys &
   (this, &
    simin, &
    idminsp, &
    nameminsp, &
    nminsp, &
    c, &
    g, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute mineral saturation indices
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)    :: this

real*8, pointer                        :: simin(:)

character(len=*), pointer              :: nameminsp(:)

integer, intent(out)                   :: nminsp

integer, intent(in)                    :: nsp

integer, pointer                       :: idminsp(:)

real*8, intent(in)                     :: c(nsp)

real*8, intent(in)                     :: g(nsp)

logical, intent(out)                   :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=300)      :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 
!%------------------------------------------------------------
call compute_min_sat_index_ &
   (this%pp, &
    simin, &
    idminsp, &
    nameminsp, &
    nminsp, &
    c, &
    g, &
    nsp, &
    iserror)
!%------------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_min_sat_index_'
print *, msg
iserror=.true.
print *,'*******************************'
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_sp_index_chemsys &
   (this, &
    namesp, &
    isps)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the global indice of some species according its name. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)    :: this

character(len=*), intent(in)           :: namesp

integer, intent(out)                   :: isps 
 
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
!% Call the corresponding service in the parent chemical 
!% system 
!%------------------------------------------------------------
call get_sp_index_(this%pp,namesp,isps)
!%------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_surf_index_chemsys &
   (this, &
    ithsurf, &
    nsites, &
    namesurf)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the global indice of some surface (including the number of sites) 
! according its name. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)    :: this

character(len=*), intent(in)           :: namesurf

integer, intent(out)                   :: ithsurf

integer, intent(out)                   :: nsites 
 
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
call get_surf_index_ &
   (this%pp, &
    ithsurf, &
    nsites, &
    namesurf)
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_surf_model_chemsys &
   (this, &
    namemodel, &
    ithsurf)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the thermodynamic model of the ith surface object. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)       :: this

character(len=*), intent(out)             :: namemodel

integer, intent(in)                       :: ithsurf
 
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
call get_surf_model_ &
   (this%pp, &
    namemodel, &
    ithsurf)
!%------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_speciation_chemsys &
   (this, &
    c, &
    g, &
    temp, &
    ionstr, &
    nsp, &
    ioutput, &
    omgwfree, &
    msg, &
    iserror, &
    simin, &
    sktrk)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii the speciation information. 
!
!   $Arguments:
!
 
 
type (t_chemicalsystem), intent(in)           :: this

integer, intent(in)                           :: nsp

integer, intent(in)                           :: ioutput

real*8, intent(in)                            :: ionstr

real*8, intent(in)                            :: temp

real*8, intent(in)                            :: omgwfree ! Mass of free water 

real*8, intent(in), dimension(nsp)            :: c

real*8, intent(in), dimension(nsp)            :: g

logical, intent(out)                          :: iserror

character(len=*), intent(out)                 :: msg

real*8, optional, intent(in), dimension(nsp)  :: simin

real*8, optional, intent(in), dimension(nsp)  :: sktrk
 
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
call write_ &
   (this%pp, &
    c, &
    g, &
    temp, &
    ionstr, &
    nsp, &
    ioutput, &
    omgwfree, &
    msg, &
    iserror, &
    simin=simin, &
    sktrk=sktrk)
!%-----------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine change_chem_unit_chemsys &
   (this, &
    c, &
    namesp, &
    unitin, &
    unitout, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Change chemical units. 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)          :: this

real*8, intent(inout)                        :: c

character(len=*), intent(in)                 :: namesp

character(len=*), intent(in)                 :: unitin

character(len=*), intent(in)                 :: unitout

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
 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

!%------------------------------------------------------------
call change_chem_unit_ &
   (this%pp, &
    c, &
    namesp, &
    unitin, &
    unitout, &
    iserror)
!%------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%************************************************************
subroutine begin_element_handler (name,attributes)
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
 
 
 
call read_xml_loc_ (name,attributes)
 
 
 
return
end subroutine
!%************************************************************
!%***************Private subroutine-**************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc_chemsys (name,attributes,this,iserror)
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
type(t_chemicalsystem), optional::this
logical, optional               ::iserror
integer        :: status, n
real*8         :: reallocal(1)
integer        :: integerlocal(1)
 
character(len=100), save:: &
 typechemsys
character(len=100):: &
 id
logical           :: &
 havethis
character(len=100)        :: &
 msg
!%-----------------------------------------------------------------
havethis = present(this)
call lowercase (name)
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
if (havethis) then
  msg=''
  iserror=.false.
  call lowercase (typechemsys)
!%-----------------
  select case (typechemsys)
  case ('classic')
   this%itype=chemsysclassic
  case ('retraso')
   this%itype=chemsysretraso
  case default
   msg='Error, not recognized chemical system type:'
   call add_ (msg,typechemsys)
  end select
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
else
 
 select case (name)
!%------------------------------------
  case ('chemicalsystem')
     id=' '
     call get_value (attributes,"type", id, status)
     call lowercase (id)
     if (id=='') then
        typechemsys='classic'
     else
        typechemsys=id
     end if
 end select
 
end if
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: read_xml_'
print *,msg
iserror=.true.
print *,'********************'
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine select_cmineq_chemsys &
   (this, &
    cmineq, &
    nsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the concentration of the non kinetic minerals
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in) :: this      ! Type parent chemical system variable

real*8, pointer, dimension(:)             :: cmineq    ! Concentration of mineral species in equilibrium

integer, intent(in)                       :: nsp       ! Number of species 

real*8, intent(in), dimension(nsp)        :: c         ! Concentration of species 

logical, intent(out)                      :: iserror   ! iserror=true, then there was an error 
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
call select_cmineq_(this%pp,cmineq,nsp,c,iserror)
if (iserror) goto 10
 
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: select_cmineq_'
print *,'*************************'
iserror=.true.
return
 


end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_density_chemsys &
    (this       ,&
    phasecode   ,&
    pliq        ,&
    temp        ,&
    c           ,&
    density     ,&
    iserror     ,&
    dc          ,&
    ddensitydc)


    

type(t_chemicalsystem), intent(in)      :: this         !Type chemical system. 

integer,intent(in)                      :: phasecode    !Phase for which the density will be calculated

real*8,intent(in)                       :: pliq         !liquid pressure

real*8,intent(in)                       :: temp         !temperature

real*8,pointer,dimension(:)             :: c         !concentration array

real*8,intent(out)                      :: density      !Density     

logical,intent(out)                     :: IsError

real*8,dimension(:,:),optional          :: dc  !concentration derivatives (needed to calculate density derivatives)

real*8,pointer,dimension(:),optional    :: ddensitydc !Density derivative

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

if  (present(ddensitydc) ) then
    if (.not.present(dc)) then
        msg="dc are needed for calcutate density derivatives"
        isError=.true.
        goto 10
    endif
    call compute_density_(this%pp,phasecode,pliq,temp,c,density,iserror,dc,ddensitydc)
    if (isError) goto 10
else
    call compute_density_(this%pp,phasecode,pliq,temp,c,density,iserror)
    if (isError) goto 10
endif


!%----------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical system:'
print *,'Service: compute_density_'
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
subroutine compute_from_setre_chemsys &
   (this, &
    ck, &
    setre, &
	nreact, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Specia according Setre=b 
!
!   $Arguments:
!
 
type (t_chemicalsystem), intent(in)     :: this

integer, intent(in)                     :: nsp

real*8, intent(inout), dimension(nsp)   :: setre

real*8, intent(in), dimension(nsp)      :: ck

logical, intent(out)                    :: iserror

Integer, intent(out)                    :: nreact

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)          :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------

call compute_from_setre_ &
   (this%pp, &
    ck, &
    setre, &
	nreact, &
    nsp, &
    iserror)

!%-------------------------------------------------------------
return
 
10 continue 
print *,'*********************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: compute_from_setre_chemsys'
print *, msg
print *,'*********************'
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
end module m_chemicalsystem
