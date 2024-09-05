module m_parentchemicalsystem
!-------------------------------------------------------------------------
!
!   $Description: Parent chemical system
!
!   $Use: m_phase, m_species, m_reaction, m_reactionratelaw, m_surface, flib_xpath, flib_sax, m_general_tools_cheproo, m_constants
!
!   $Author: Sergio Andrés Bea Jofré 
!
!   $License: UPC-CSIC
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_phase
use m_species
use m_reaction
use m_reactionratelaw
use m_surface
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%------------------------------------------------------
!%------------------------------------------------------
private                      ::
!%--------------------------------------------------------
!% Constant parameters per default 
!%--------------------------------------------------------
real*8, public, parameter    :: &
tolresnr=1.0d-10,  &    ! Tolerance in residual (for newton raphson implemented in any speciation algorithm) 
tolunknr=1.0d-5,   &    ! Tolerance in unknown (for newton raphson implemented in any speciation algorithm) (
facmax=0.5d0,        &  ! Factor of correction of unknowns
zero=1.0d-300,       &  ! Value for below is consider 0.0d0
maxanomc=1.0d20,     &  ! Maximum value for anomalous concentrations
minanomc=1.0d-300       ! Minimum value for anomalous concentrations
integer, public, parameter   :: &
mxiter=300,          &  ! Maximum interations for to solve non_linea
mxtry=20,            &  ! Maximum try for to find stable ste of mine
mxdiverg=30,         &  ! Maximum divergence considered 
mxitergam=100,       &  ! Maximum number of iteration of gammas
iouwiinfo=90909         ! Unit for to write chemical system info
!%------------------------------------------------------
!%------------------------------------------------------
public                            :: &
create_ &                      ! Create the parent chemical system. 
,read_xml_ &                   ! Read the parent chemical system from xml file. 
,read_txt_ &                   ! Read the parent chemical system from ascii input data (at the moment read the inputa file performed by retraso)
,destroy_ &                    ! Destroy the parent chemical system.
,set_ &                        ! Set attributes in the parent chemical system object according other types (e.g. reactions, phases, etc.)
,write_ &                      ! Write all atributes encapsulated in the parent chemical system in ascii format 
,get_aq_sps_pos_ &             ! Return the first and last indices of the aqueous spcecies. 
,get_idreaction_ &             ! Return the global indices of reactions in or between two phases.  
,get_min_sp_index_ &           ! Return the global indices of mineral species. 
,get_gas_sp_index_ &           ! Return the global indices of gas species. 
,get_chem_info_ &              ! Return general chemical information.
,get_sp_index_ &               ! Return the global indice of ith species according its name. 
,get_surf_index_ &             ! Return the global indice of surface according its name. 
,get_iposspsph_ &              ! Return the first and last global indices corresponding to species in the ith phase.  
,get_iposspsurf_ &             ! Return the first and last global indices corresponding to species in the ith surface.
,get_new_chemsys_ &            ! Return a new chemical system.
,get_new_chemical_base_ &      ! Get the name of the species of the new chemical base. 
,get_surf_model_ &             
,change_chem_unit_ &           ! Change chemical units. 
,switch_base_ &                ! Switch the set of primary species in the chemical system object. 
,check_mineral_ &              ! Check the mineral set 
,check_gas_ &                  ! Check the if there are negative gas species 
,add_and_specia_ &
,equilibrate_ &                ! Equilibrate solution with mineral phase.
,specia_from_solution_type_ &
,specia_from_cpri_ &           ! Compute the speciation from concentration of primary species
,specia_from_u_ &              ! Make the chemical speciation from total of components in all phases. 
,compute_min_sat_index_ &      ! Compute saturation indices of minerals 
,compute_molsalt_ &             
,compute_mass_salt_ &          ! Compute Mass of solute in a phase
,compute_dmassSalt_dc_ &       ! Compute derivative of Mass of solute in a phase
,compute_omgwcryst_ &          ! Compute mass of water fixed in hydrated minerals.
,compute_density_ &            ! Compute density (and derivatives if asked)
,compute_from_setre_ &         ! Compute x as: Transpose(Se)*x=b 
,select_cmob_ &
,select_caq_ &
,select_cads_ &
,select_cmin_ &
,select_cgas_ &
,select_cmineq_&               ! Select the concentration of the non kinetic minerals
,compute_dcmob_ &              ! Compute dcmob/dc1
,compute_dcads_ &              ! Compute dcads/dc1
,compute_dumob_ &              ! Compute U*dcmob/dc1
,compute_kinetic_ &            ! Compute the kinetic calculations. 
,update_min_area_ &            ! Update the reactive surface of minerals 
,update_txoh_ &                ! Update TXOH 
,compute_usktrk_ &             ! Compute U*Skt*rk
,compute_dusktrk_ &            ! Compute derivative of changes in components due kinetic reactions with respect to primary aqueous species. 
,compute_umob_ &
,compute_iumob_ &
,compute_ith_total_ &
,compute_uads_ &
,compute_r_from_stqtr_ &       ! Compute r solving the following system equations.
,compute_alkalinity_ &
,compute_secondaries_ &
,update_ &                     ! Update parameters depending of the temperature in the chemical system
,assignment(=) &               ! Copy a parent chemical system object in other parent chemical system object
,begin_xml_pchemsys &
,check_inv_points_ &           ! Check and determine the reaction set that constraint the water activity. 
,set_iswcompbal_ &
,set_iswriteinfo_
!%--------------------------------------------------------
!%--------------------------------------------------------
private                           :: &
build_chemical_base_ &
,build_stq_ &
,build_components_matrix_sing_values_ &
,read_prim_master25_ &
,read_scnd_master25_ &
,read_ads_sys_ &
,read_genr_master25_ &
,read_temp_master25_ &
,read_aqsp_sys_ &
,read_miga_sys_ &
,get_ph_index_ &
,check_convergence_ &
,get_eq_min_sp_index_ &
,get_num_eq_min_sp_ &
,read_xml_loc_ &
,begin_element_handler &
,check_ &
,set_from_master25_ &
,set_from_phreeqc_ &
,read_from_datakin_ 
!%------------------------------------------------------
!%------------------------------------------------------
!%------------------------------------------------------
!% Type definition 
!%------------------------------------------------------
!%------------------------------------------------------
!%------------------------------------------------------
type, public::t_parentchemicalsystem
 
character(len=100)                      :: name             ! Name of chemical system
 
type(t_pspecies), pointer, dimension(:) :: pspecies         ! List of pointers to specie type objects [numsp]
 
type(t_pphase), pointer, dimension(:)   :: pphase           ! List of pointers to phase type objects [numph]
 
type(t_psurface), pointer, dimension(:) :: psurface         ! List of pointers to surface type objects [numsurface]
 
type(t_preaction), pointer, dimension(:):: preaction        ! List of pointers to reaction type objects [numreact]
 
type(t_prrlaw), pointer, dimension(:)   :: prrlaw           ! List of pointers to reaction rate law objects
 
real*8, pointer, dimension(:,:)         :: ueq              ! Components matrix (are included only equilibrium reactions) [numcomp,numsp]

real*8, pointer, dimension(:,:)         :: u                ! Components matrix (are included all reactions) [numcomp,numsp]

real*8, pointer, dimension(:,:)         :: stq              ! Stoichiometric matrix (are included all reactions) [numreac,numsp]
 
integer, pointer, dimension(:)          :: idminph          ! Global index of mineral phases [numminph]

integer, pointer, dimension(:)          :: idgasph          ! Global index of gas phases [numgasph]

integer, pointer, dimension(:)          :: idreactsp        ! Indices of secondary species [numreact]

integer, pointer, dimension(:)          :: idaqprisp        ! Global index of primary aqueous species [numaqprisp]

integer                                 :: numaqprisp       ! Number of aqueous primary species

integer                                 :: numph            ! Number of phases 

integer                                 :: numsurf          ! Number of surfaces 

integer                                 :: numreact         ! Number of reactions objects
 
integer                                 :: numrrlaw         ! Number of reaction rate laws 

integer                                 :: numsp            ! Number of species

integer                                 :: numsites         ! Total number of sites for adsorption or cation exchange 

integer                                 :: numminph         ! Number of mineral phases 

integer                                 :: numgasph         ! Number of gas phases

integer                                 :: aqphindex        ! Global index of the aqueous phase

integer                                 :: wcindex          ! Global index of the water component

integer                                 :: ecindex          ! Global index of the electron component

integer                                 :: hcindex          ! Global index of the proton component
 
logical, pointer, dimension(:)          :: iskinreact       ! [numreact] 
                                                            ! iskinreact=true, the reaction is kinetic
                                                            ! iskinreact=false, the reaction is in equilibrium

logical, pointer, dimension(:,:)        :: ishomreact       ! .true. if the reaction is homogeneous [numph+numsurf,numreact]

logical, pointer, dimension(:,:)        :: ishetreact       ! .true. if the reaction is heterogenous [numph+numsurf,numreact]

logical                                 :: isgaussjordan    ! .true. if the components matrix is com
                                                            ! using gauss-jordan elimination.
                                                            ! .false. if the components matrix is co
                                                            ! using singular values descomposition
 
real*8                                  :: tempref          ! Reference temperature

real*8                                  :: pressref         ! Pressure of reference 
 
logical                                 :: iswriteinfo      ! If .true. then write Newton Raphson information 
                                                            ! when make the speciation
 
real*8                                  :: tolunknr         ! tolerance in unknown for newton raphson 

real*8                                  :: tolresnr         ! teolerance in residual for newton raphson

real*8                                  :: deltasatmin      ! Thershold of the saturation of minerals 

logical                                 :: iswcompbal       ! If true then evaluate the water component equation during speciation 
 
end type t_parentchemicalsystem
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_iswcompbal_  
 
module procedure set_iswcompbal_pchemsys 
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_iswriteinfo_
 
module procedure set_iswriteinfo_pchemsys 
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface create_
 
module procedure create_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_xml_
 
module procedure read_xml_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_txt_
 
module procedure read_txt_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface destroy_
 
module procedure destroy_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_
 
module procedure set_from_types_pchemsys
module procedure set_from_data_base_info_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface write_
 
module procedure write_pchemsys
module procedure write_speciation_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_chem_info_
 
module procedure get_chem_info_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_aq_sps_pos_
 
module procedure get_aq_sps_pos_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_idreaction_
 
module procedure get_idreaction_pchemsys
module procedure get_idreaction_from_sps_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_sp_index_
 
module procedure get_sp_index_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_surf_index_
 
module procedure get_surf_index_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_min_sp_index_
 
module procedure get_min_sp_index_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_gas_sp_index_
 
module procedure get_gas_sp_index_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_iposspsph_
 
module procedure get_iposspsph_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_iposspsurf_
 
module procedure get_iposspsurf_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_new_chemsys_
 
module procedure get_new_chemsys_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_surf_model_
 
module procedure get_surf_model_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_new_chemical_base_
 
module procedure get_new_chemical_base_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_iumob_
 
module procedure compute_iumob_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_ith_total_
 
module procedure compute_ith_total_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_kinetic_
 
module procedure compute_kinetic_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface update_min_area_
 
module procedure update_min_area_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface update_txoh_
 
module procedure update_txoh_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface change_chem_unit_
 
module procedure change_chem_unit_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface switch_base_
 
module procedure switch_base1_pchemsys
module procedure switch_base2_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cmob_
 
module procedure select_cmob_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_caq_
 
module procedure select_caq_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cads_
 
module procedure select_cads_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cmin_
 
module procedure select_cmin_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cgas_
 
module procedure select_cgas_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dcmob_
 
module procedure compute_dcmob1_pchemsys
module procedure compute_dcmob2_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dcads_
 
module procedure compute_dcads_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dumob_
 
module procedure compute_dumob1_pchemsys
module procedure compute_dumob2_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface check_mineral_
 
module procedure check_mineral_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface check_gas_
 
module procedure check_gas_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface specia_from_solution_type_
 
module procedure specia_from_solution_type_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface equilibrate_
 
module procedure equilibrate_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface add_and_specia_
 
module procedure add_and_specia_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface specia_from_cpri_
 
module procedure specia_from_cpri_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface specia_from_u_
 
module procedure specia_from_u_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_min_sat_index_
 
module procedure compute_min_sat_index_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_molsalt_
 
module procedure compute_molsalt_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_mass_salt_
 
module procedure compute_mass_salt_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dmassSalt_dc_
 
module procedure compute_dmassSalt_dc_pchemsys
 
end interface

!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_omgwcryst_
 
module procedure compute_omgwcryst_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface update_
 
module procedure update_temp_param_pchemsys
module procedure update_porosity_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_usktrk_
 
module procedure compute_usktrk1_pchemsys
module procedure compute_usktrk2_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dusktrk_
 
module procedure compute_dusktrk1_pchemsys
module procedure compute_dusktrk2_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_umob_
 
module procedure compute_umob_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_uads_
 
module procedure compute_uads_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_r_from_stqtr_
 
module procedure compute_r_from_stqtr_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_alkalinity_
 
module procedure compute_alkalinity_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_secondaries_
 
module procedure compute_secondaries_pchemsys
module procedure compute_adsorption_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface assignment (=)
 
module procedure copy_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface check_inv_points_
 
module procedure check_inv_points_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface select_cmineq_
 
module procedure select_cmineq_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_density_
 
module procedure compute_density_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_from_setre_
 
module procedure compute_from_setre_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------

!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%-------------------Private services---------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_chemical_base_
 
module procedure build_chemical_base_pchemsys

end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_stq_
 
module procedure build_stq_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_components_matrix_sing_values_
 
module procedure build_components_matrix_sing_values_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_components_matrix_gauss_jordan_
 
module procedure build_components_matrix_gauss_jordan_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_prim_master25_
 
module procedure read_prim_master25_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_scnd_master25_
 
module procedure read_scnd_master25_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_ads_sys_
 
module procedure read_ads_sys_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_miga_sys_
 
module procedure read_miga_sys_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_genr_master25_
 
module procedure read_genr_master25_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_temp_master25_
 
module procedure read_temp_master25_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_aqsp_sys_
 
module procedure read_aqsp_sys_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_ph_index_
 
module procedure get_ph_index_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface check_convergence_
 
module procedure check_convergence_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_eq_min_sp_index_
 
module procedure get_eq_min_sp_index_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_num_eq_min_sp_
 
module procedure get_num_eq_min_sp_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_xml_loc_
 
module procedure read_xml_loc_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface begin_element_handler
 
module procedure begin_element_handler
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface check_
 
module procedure check_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_from_master25_ 
 
module procedure set_from_master25_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_from_datakin_ 
 
module procedure read_from_datakin_pchemsys
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_from_phreeqc_ 
 
module procedure set_from_phreeqc_pchemsys
 
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
subroutine create_pchemsys &
   (this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create parent chemical system
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout) :: this   ! Type parent chemical system variable 
 
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
this%pphase => null ()
this%psurface => null ()
this%preaction=> null ()
this%pspecies => null ()
this%idreactsp => null ()
this%iskinreact => null ()
this%idminph => null ()
this%idgasph => null ()
this%ueq => null ()
this%u => null ()
this%stq => null ()
this%idaqprisp => null ()
this%ishomreact => null ()
this%ishetreact => null ()
this%pspecies => null ()
this%prrlaw => null ()
!%-----------------------------------------------------------
this%name=' '
this%numph = 0
this%numsurf = 0
this%numsp = 0
this%numreact = 0
this%numrrlaw = 0
this%aqphindex = 0
this%isgaussjordan = .false.
this%tempref = 25.0d0
this%pressref = 0.0d0
this%numminph = 0
this%numgasph = 0
this%numaqprisp = 0
this%numsites = 0
this%wcindex=0
this%ecindex=0
this%hcindex=0
this%iswriteinfo=.false.
this%iswcompbal =.false. 
!%------------------------------------------------------------
this%tolunknr=tolunknr
this%tolresnr=tolresnr
this%deltasatmin=0.0d0 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_pchemsys &
   (this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy parent chemical system object
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout) :: this 
 
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
if (this%numsp>0) then
 do i=1,this%numsp
  this%pspecies(i)%ptr => null ()
 end do
 deallocate (this%pspecies)
 this%pspecies => null ()
end if
!%----------------------------------------------------
if (this%numph>0) then
 do i=1,this%numph
   call destroy_ (this%pphase(i)%ptr)
   deallocate (this%pphase(i)%ptr)
   this%pphase(i)%ptr => null () 
 end do
 deallocate (this%pphase)
 this%pphase => null ()
end if
!%----------------------------------------------------
if (this%numsurf>0) then
 do i=1,this%numsurf
   call destroy_ (this%psurface(i)%ptr)
   deallocate (this%psurface(i)%ptr)
   this%psurface(i)%ptr => null ()
 end do
 deallocate (this%psurface)
 this%psurface => null ()
end if
!%----------------------------------------------------
if (this%numrrlaw>0) then
 do i=1,this%numrrlaw
  this%prrlaw(i)%ptr => null ()  
 end do
 deallocate (this%prrlaw)
 this%prrlaw => null ()
end if
!%----------------------------------------------------
if (this%numreact>0) then
 do i=1,this%numreact
  call destroy_ (this%preaction(i)%ptr)
  this%preaction(i)%ptr => null ()
 end do
 deallocate (this%preaction)
 this%preaction=> null ()
end if
!%----------------------------------------------------
!% Deallote pointers 
!%----------------------------------------------------
call check_pointer_ (this%ueq,1,1,.false.)
call check_pointer_ (this%u,1,1,.false.)
call check_pointer_ (this%stq,1,1,.false.)
call check_pointer_ (this%idreactsp,1,.false.)
call check_pointer_ (this%idaqprisp,1,.false.)
call check_pointer_ (this%iskinreact,1,.false.)
call check_pointer_ (this%idgasph,1,.false.)
call check_pointer_ (this%idminph,1,.false.)
call check_pointer_ (this%ishomreact,1,1,.false.)
call check_pointer_ (this%ishetreact,1,1,.false.)
!%----------------------------------------------------
!% Zeroing attributtes 
!%----------------------------------------------------
this%numaqprisp = 0
this%numph = 0
this%numsurf = 0
this%numreact = 0
this%numrrlaw = 0
this%numsp = 0
this%aqphindex = 0
this%numminph = 0
this%numgasph = 0
this%numsites = 0
this%wcindex = 0
this%ecindex = 0
this%hcindex = 0
this%isgaussjordan = .false.
this%iswriteinfo = .false.
this%tolresnr=0.0d0
this%tolunknr=0.0d0
this%pressref = 0.0d0
this%tempref = 0.0d0
this%deltasatmin=0.0d0 
this%iswcompbal =.false.
!%---------------------------------------------------
!% Close unit of write chemical system information
!%---------------------------------------------------
if (this%iswriteinfo) then
 close (unit=iouwiinfo)
end if
this%iswriteinfo=.false.
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_iswcompbal_pchemsys &
   (this, &
    iswcompbal, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: active the water component mass balance 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout) :: this

logical, intent(in)                          :: iswcompbal 

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
iserror=.false.
!%--------------------------------------------------------
this%iswcompbal=iswcompbal
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
subroutine set_iswriteinfo_pchemsys &
   (this, &
    iswriteinfo, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set if the chemical system object must write information 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout) :: this

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
iserror=.false.
!%--------------------------------------------------------
this%iswriteinfo=iswriteinfo 
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
subroutine read_xml_pchemsys &
   (this, &
    namefile, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read the chemical system from xml file 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout) :: this

character(len=*), intent(in)                 :: namefile

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
! Open xml file
!%--------------------------------------------------------
call open_xmlfile (namefile,fxml,iostat)
if (iostat/=0) then
 msg='Error when open xml file'
 call add_ (msg,namefile)
 goto 10
end if
!%--------------------------------------------------------
! Read parent chemical system
!%--------------------------------------------------------
call xml_parse &
 (fxml, &
  begin_element_handler=begin_element_handler)
!%--------------------------------------------------------
call read_xml_loc_ (name,attributes,this,iserror)
if (iserror) goto 10
!%--------------------------------------------------------
! End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%--------------------------------------------------------
if (this%iswriteinfo) then
  open(unit=iouwiinfo,file='chemsys_info.dat',status='unknown')
end if
!%--------------------------------------------------------
return
 
10 continue 
print *,'******************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: read_xml_'
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
subroutine read_txt_pchemsys &
   (this, &
    namefile, &
	namethdb, &
    namekindb, &
	filebase, &
	ioptxhl, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read the chemical system from txt file. 
! Read the chemical input file performed by Retraso 
!
!   $Arguments:
!
 
 
type (t_parentchemicalsystem), intent(inout):: this

character(len=*), intent(in)                :: namefile     ! Name of chemical system definition file

character(len=*), intent(in)                :: namethdb     ! Name of thermodynamic database

character(len=*), intent(in)                :: namekindb    ! Name of kinetic database

character(len=*), intent(in)                :: filebase

integer, intent(in)                         :: ioptxhl      ! Option to compute slat mass fraction 

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
character(len=100)  :: &
 msg
character(len=100)   :: &
 namtyp
character(len=100), pointer:: &
 napri(:) => null(), & 
 naaqx(:) => null(), &
 naaqt(:) => null(), &
 nagas(:) => null(), &
 naads(:) => null(), &
 namin(:) => null(), &
 label(:) => null(), &
 labelm(:) => null(), &
 naadsmod(:) => null()
integer, pointer          :: &
 idmeq(:) => null(), &
 nadsmod(:) => null()
real*8                    :: &
 tempisoterm 
integer                   :: &
 iact, &
 iconv, &
 ngas, &
 naqx, &
 naqt, &
 nads, &
 npri, &
 naqpri, &
 nadspri, &
 nmin, &
 itemp, &
 nmineq, &
 nminkin, &
 nmod, &
 lbase
integer, parameter        :: &
 mxsp=50, &
 mxlabel=5
character(len=100):: &
 name
!-------------------------------------------------------------------------
!
!   $code
!
!%--------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------
!%Open the unit
!%---------------------------------------------------------
name=namefile
call lastletter_ (lbase,filebase)
name=filebase(1:lbase)//name
open(unit=1,file=name,status='old',err=30)
!%---------------------------------------------------------
!%Allocate local pointers 
!%---------------------------------------------------------
call check_pointer_ (napri,mxsp,.true.)
call check_pointer_ (naaqx,mxsp,.true.)
call check_pointer_ (naaqt,mxsp,.true.)
call check_pointer_ (nagas,mxsp,.true.)
call check_pointer_ (namin,mxsp,.true.)
call check_pointer_ (label,mxlabel,.true.)
call check_pointer_ (naadsmod,mxsp,.true.)
call check_pointer_ (labelm,mxlabel,.true.)
call check_pointer_ (naads,mxsp,.true.)
call check_pointer_ (idmeq,mxsp,.true.)
call check_pointer_ (nadsmod,mxsp,.true.)
!%---------------------------------------------------------
!%Read the aqueous species defined in the chemical system
!%---------------------------------------------------------
call read_aqsp_sys_ (itemp,iact,iconv,naqt,npri,tempisoterm, &
                     label,naaqt,mxsp,mxlabel,iserror)
if (iserror) goto 20          
!%---------------------------------------------------------
!%Read mineral and gases defined in the chemical system
!%---------------------------------------------------------
call read_miga_sys_ &
  (ngas,nmin,nmineq,nminkin,naqt,naaqt,idmeq,labelm,nagas,namin, &
   mxlabel,mxsp,iserror)
if (iserror) goto 20  
!%---------------------------------------------------------
!% Compute the number of primary aqueous species 
!%---------------------------------------------------------
naqx=naqt-npri
napri(1:npri)=naaqt(1:npri)
naaqx(1:naqx)=naaqt(npri+1:npri+naqx)
!%---------------------------------------------------------
!%Read the surface complexes defined in the chemical system
!%---------------------------------------------------------
call read_ads_sys_ (nads,nmod,npri,napri,naads,naadsmod, &
                    nadsmod,mxsp,mxlabel,iserror)
if (iserror) goto 20   
!%---------------------------------------------------------
!% Set the chemical system object according thermodynamic
!% data base 
!%---------------------------------------------------------
nadspri=nmod 
naqpri=npri-nadspri
select case (namethdb)
case ('master25.dat','MASTER25.DAT')
  call set_from_master25_ &
   (this, &
    tempisoterm, &
    iact, &
    iconv, &
    ioptxhl, &
    naqpri, &
    nadspri, & 
    naqx, &
    nmin, &
    ngas, &
    nads, &
    nmod, &
    nadsmod, &
    napri, &
    naaqx, &
    namin, &
    nagas, &
    naads, &
    naadsmod, &
    idmeq, &
    namekindb, &
	filebase, &
    iserror)
case ('phreeqc.dat','PHREEQC.DAT')
  call set_from_phreeqc_ &
   (this, &
    tempisoterm, &
    iact, &
    iconv, &
    ioptxhl, &
    naqpri, &
    nadspri, & 
    naqx, &
    nmin, &
    ngas, &
    nads, &
    nmod, &
    nadsmod, &
    napri, &
    naaqx, &
    namin, &
    nagas, &
    naads, &
    naadsmod, &
    idmeq, &
    namekindb, &
	filebase, &
    iserror)    
 case default
  msg='Error, not defined thermodynamic data base:'
  call add_ (msg,namethdb)
  iserror=.true.
  goto 20 
 end select 
!%---------------------------------------------------------
!% Close unit
!%---------------------------------------------------------
close(unit=1)
!%---------------------------------------------------------
20 continue
!%---------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------- 
call check_pointer_ (napri,1,.false.)
call check_pointer_ (naaqx,1,.false.)
call check_pointer_ (naaqt,1,.false.)
call check_pointer_ (nagas,1,.false.)
call check_pointer_ (namin,1,.false.)
call check_pointer_ (label,1,.false.)
call check_pointer_ (naads,1,.false.)
call check_pointer_ (labelm,1,.false.)
call check_pointer_ (naadsmod,1,.false.)
call check_pointer_ (idmeq,1,.false.)
if (iserror) goto 10 
return
 
10 continue 
print *,'******************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: read_txt_'
print *, msg
print *,'******************'
iserror=.true.
return
 
30 continue
msg='Error when open file:'
call add_ (msg,namefile)
goto 10  


end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_from_master25_pchemsys &
   (this, &
    tempisoterm, &
    iact, &
    iconvact, &
    ioptxhl, &
    naqpri, &
    nadspri, &
    naqx, &
    nmin, &
    ngas, &
    nads, &
    nmod, &
    nadsmod, &
    napri, &
    naaqx, &
    namin, &
    nagas, &
    naads, &
    naadsmod, &
    idmeq, &
    namekindb, &
	filebase, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the chemical system object according master25.dat 
! (or mastertemp.dat) thermodynamic database. 
! The kinetic information is read of kinetic.dat file. 
!
!   $Arguments:
!
 
 
type (t_parentchemicalsystem), intent(inout)             :: this        ! Parent chemical system type 

real*8, intent(in)                                       :: tempisoterm ! Temperature 

integer, intent(in)                                      :: naqpri      ! Number of aqueous primary species

integer, intent(in)                                      :: nadspri     ! Number of adsorption species 

integer, intent(in)                                      :: naqx        ! Number of aqueous complexes species  

integer, intent(in)                                      :: nmin        ! Number of minerals species 

integer, intent(in)                                      :: ngas        ! Number of gas species 

integer, intent(in)                                      :: nads        ! Number of sorbed species 

integer, intent(in)                                      :: nmod        ! Number of surface complexes models. 

integer, intent(in), dimension(nmod)                     :: nadsmod     ! Number of species for surface complex model 

character(len=*), intent(in), dimension(naqpri+nadspri)  :: napri       ! Name of primary species

character(len=*), intent(in), dimension(naqx)            :: naaqx       ! Name of aqueous complexes

character(len=*), intent(in), dimension(nmin)            :: namin       ! Name of mineral species

character(len=*), intent(in), dimension(ngas)            :: nagas       ! Name of gas species

character(len=*), intent(in), dimension(nads)            :: naads       ! Name of sorbed species

character(len=*), intent(in), dimension(nmod)            :: naadsmod    ! Name of models where forming surface complexes (or cation exchange species)

integer, intent(in), dimension(nmin)                     :: idmeq       ! Index of equilibrium/kinetic minerals 

character(len=*), intent(in)                             :: namekindb   ! Name of kinetic database

character(len=*), intent(in)                             :: filebase    ! Path of input files

integer, intent(in)                                      :: iact        ! Name of aqueous activity model 

integer, intent(in)                                      :: iconvact    ! Index of convention 

integer, intent(in)                                      :: ioptxhl     ! Name of aqueous activity model  

logical, intent(out)                                     :: iserror     ! If true there was error 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)  :: &
 msg
character(len=100)   :: &
 namtyp
character(len=100), pointer:: &
 namsec(:) => null (), &
 label(:) => null (), &
 labelm(:) => null (), &
 namevalue(:) => null (), &
 nacat(:,:,:) => null ()
real*8, pointer           :: &
 zt(:) => null (), &
 a0t(:) => null (), &
 wmolt(:) => null (), &
 bdottot(:) => null (), &
 alogk(:,:) => null (), &
 stcrc(:,:) => null (), &
 stcpr(:,:) => null (), &
 stqt(:,:) => null (), &
 vmin(:,:) => null (), &
 propgas(:,:) => null (), &
 propads(:,:) => null (), &
 value(:) => null (), &
 s2lu(:,:) => null (), &
 temp(:)  => null (), &
 ea(:) => null (), &
 thresh(:) => null (), &
 rk(:,:) => null (), &
 pcat(:,:,:) => null (), &
 fora(:,:) => null (), &
 dins (:,:) => null ()
integer, pointer          :: &
 indxpri(:)  => null (), &
 indxsec(:) => null (), &
 indx(:) => null (), &
 indxgas(:) => null (), &
 indxmin(:) => null (), &
 indxads(:) => null (), &
 nkin(:) => null (), &
 ncat(:,:) => null ()
integer                   :: &
 i, &
 j, &
 isp, &
 ireact, &
 iph, &
 ntemp, &
 nbasis, &
 nreac, &
 nprop, &
 iopgastr, &
 ioutput, &
 itemp, &
 lbase, &
 nsp, &
 nkinmin 
integer, parameter        :: &
 npropgas=2, &
 mxlabel=5, &
 mxkinterm=5, &
 ncoeff=5, &
 iuthdb=5, &
 iprint=0, &
 mxdim=200
real*8, parameter         :: &
 tempref=25.0d0 
character(len=100):: &
 name
!-------------------------------------------------------------------------
!
!   $code
!
 
!%--------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------
!% Compute the total number of species 
!%--------------------------------------------------------
nsp=naqpri+nadspri+naqx+nmin+ngas+nads
!%--------------------------------------------------------
!% Check if there are kinetic minerals
!%--------------------------------------------------------
nkinmin=0 
do i=1,nmin
 if (idmeq(i)/=0) then
  nkinmin=nkinmin+1
 end if
end do
!%--------------------------------------------------------
!% Allocate local pointers 
!%--------------------------------------------------------
nprop=1
call check_pointer_ (label,mxlabel,.true.)
call check_pointer_ (temp,mxdim,.true.)
call check_pointer_ (bdottot,nsp,.true.)
call check_pointer_ (zt,nsp,.true.)
call check_pointer_ (a0t,nsp,.true.)
call check_pointer_ (wmolt,nsp,.true.)
call check_pointer_ (alogk,nsp,mxdim,.true.)
call check_pointer_ (stcpr,nsp,naqpri+nadspri,.true.)
call check_pointer_ (stcrc,nsp,nsp,.true.)
call check_pointer_ (label,mxlabel,.true.)
call check_pointer_ (labelm,mxlabel,.true.)
call check_pointer_ (vmin,nprop,nmin,.true.)
call check_pointer_ (propads,nprop,nads,.true.)
call check_pointer_ (propgas,npropgas,ngas,.true.)
call check_pointer_ (indxads,nads,.true.)
call check_pointer_ (indx,mxdim,.true.)
call check_pointer_ (indxmin,nmin,.true.)
call check_pointer_ (indxsec,naqx+nmin+ngas+nads,.true.)
call check_pointer_ (indxpri,naqpri+nadspri,.true.)
call check_pointer_ (indxgas,ngas,.true.)
call check_pointer_ (namsec,nsp,.true.)
!%--------------------------------------------------------
!% Allocate local pointers corresponding to kinetic 
!% information 
!%--------------------------------------------------------
call check_pointer_ (ea,nmin,.true.)
call check_pointer_ (nkin,nmin,.true.)
call check_pointer_ (thresh,nmin,.true.)
call check_pointer_ (dins,nmin,mxkinterm,.true.)
call check_pointer_ (fora,nmin,mxkinterm,.true.)
call check_pointer_ (rk,nmin,mxkinterm,.true.)
call check_pointer_ (pcat,nmin,mxkinterm,naqpri+naqx+nmin,.true.)
call check_pointer_ (ncat,nmin,mxkinterm,.true.)
call check_pointer_ (nacat,nmin,mxkinterm,naqpri+naqx+nmin,.true.) 
!%---------------------------------------------------------
!% Open file containing the chemical database (master25)
!%---------------------------------------------------------
if (tempisoterm==tempref) then 
 write(6,*) '=======> Reading database:','master25.dat'
 call lastletter_ (lbase,filebase)
!cprovi name=filebase(1:lbase)//'master25.dat'
 name='master25.dat'
else
 write(6,*) '=======> Reading database:','mastertemp.dat'
 call lastletter_ (lbase,filebase)
 name=filebase(1:lbase)//'mastertemp.dat'
end if 
!%---------------------------------------------------------
!% Open the thermodynamic data base
!%---------------------------------------------------------
open(unit=iuthdb,file=name,status='old',err=20)
!%---------------------------------------------------------
!%Read temperature points and build ATBSFN
!%---------------------------------------------------------
call read_temp_master25_(ntemp,temp,nbasis,mxdim,iuthdb,ioutput)
!%---------------------------------------------------------
!%Read primary species in the thermodynamic data base
!%---------------------------------------------------------
call read_prim_master25_ &
 (this,naqpri+nadspri,naqx,nads,nsp,iprint,iuthdb,ioutput,napri, &
  naaqx,naads,propads,zt,a0t,wmolt,bdottot,ioptxhl,iact,iserror)
if (iserror) then
  msg='Error when call internal service read_prim_master25_'
  goto 30
end if
!%---------------------------------------------------------
!%Read secundary species of the chemical system
!%---------------------------------------------------------
call read_scnd_master25_ &
 (naqpri+nadspri,naqx,ngas,naqx+nmin,ntemp,mxdim,nbasis, &
  nsp,nreac,iprint,6,iuthdb,alogk(1:nsp,1:ntemp), &
  napri,indxpri,naaqx,indxsec,nagas,indxgas, &
  indx,namsec,stcpr,stcrc,zt,a0t,wmolt,bdottot, &
  ioptxhl,iact)
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%Read thermodynamic data base (for minerals)
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
namtyp='minerals'
call read_genr_master25_ &
 (namtyp,nmin,nprop,naqpri+nadspri,naqx,ngas,naqx,ntemp, &
  mxdim,nbasis,nsp,nreac,iprint,6,iuthdb, &
  namin,indxmin,vmin,alogk(1:nsp,1:ntemp),napri,indxpri,naaqx, &
  indxsec,nagas,indxgas,namsec,stcpr, &
  stcrc,iopgastr,zt)
!%---------------------------------------------------------
!%Read thermodynamic data base (for gases)
!%---------------------------------------------------------
namtyp='gasses'
call read_genr_master25_ &
 (namtyp,ngas,npropgas,naqpri+nadspri,naqx,ngas,naqx+nmin,ntemp, &
  mxdim,nbasis,nsp,nreac,iprint,6,iuthdb, &
  nagas,indxgas,propgas,alogk(1:nsp,1:ntemp),napri,indxpri,naaqx, &
  indxsec,nagas,indxgas,namsec,stcpr, &
  stcrc,iopgastr,zt)
!%---------------------------------------------------------
!%Read thermodynamic data base (adsorbed species)
!%---------------------------------------------------------
namtyp='adsorbed species'
nprop=1
call read_genr_master25_ &
 (namtyp,nads,nprop,naqpri+nadspri,naqx,nads,naqx+nmin+ngas,ntemp, &
  mxdim,nbasis,nsp,nreac,iprint,6,iuthdb, &
  naads,indxads,propads,alogk(1:nsp,1:ntemp),napri,indxpri,naaqx, &
  indxsec,naads,indxads,namsec,stcpr, &
  stcrc,iopgastr,zt)    
!%---------------------------------------------------------
!% Change basis from master25 to chemical system 
!% (Compute STQT)
!%---------------------------------------------------------
if (nreac>0) then 
  call check_pointer_ (stqt,nreac,nsp,.true.)
  call check_pointer_ (s2lu,nreac,nreac,.true.)
  call check_pointer_ (indx,nreac,.true.)
  call switch_base_ &
  (naqpri+nadspri,nreac,ntemp,nsp, & 
   stcrc(1:nreac,1:nreac),stcpr(1:nreac,1:naqpri+nadspri), &
   alogk(1:nreac,1:ntemp),s2lu,stqt,indx)
end if
!%---------------------------------------------------------
!% Read kinetic data base  
!%---------------------------------------------------------
if (nkinmin>0) then 
 call read_from_datakin_ (naqpri,naqx,nmin,nkinmin,mxkinterm,napri(1:naqpri), &
  naaqx,idmeq,namin,ea,nkin,thresh,dins,fora,rk,pcat,ncat, &
  nacat,filebase,namekindb,iserror)
 if (iserror) then
    msg='Error reading kinetic data base'
	goto 30
 end if
end if 
!%---------------------------------------------------------
!% Set in the chemical system object  
!%---------------------------------------------------------
call set_ (this,tempisoterm,iact,iconvact,naqpri,nadspri,naqx, &
           nmin,nads,ngas,ntemp,nreac,ncoeff,nmod,mxkinterm,temp(1:ntemp), &
		   napri,naaqx,namin,naads, &
		   nagas,naadsmod,nadsmod,idmeq,a0t(1:naqpri+nadspri+naqx), &
		   zt(1:naqpri+nadspri+naqx),wmolt(1:naqpri+nadspri+naqx), &
		   bdottot(1:naqpri+nadspri+naqx), &
		   vmin(1,1:nmin),propgas(1,1:ngas),propads(1,1:nads), &
		   alogk(1:nreac,1:ntemp),stqt,ea,nkin,thresh,dins,fora,rk, &
		   pcat,ncat,nacat,iserror)		   
!%---------------------------------------------------------
30 continue 
!%---------------------------------------------------------
!% Deallocate local pointers
!%---------------------------------------------------------
call check_pointer_ (bdottot,1,.false.)
call check_pointer_ (zt,1,.false.)
call check_pointer_ (a0t,1,.false.)
call check_pointer_ (wmolt,1,.false.)
call check_pointer_ (indxsec,1,.false.)
call check_pointer_ (indxpri,1,.false.)
call check_pointer_ (indxgas,1,.false.)
call check_pointer_ (alogk,1,1,.false.)
call check_pointer_ (stcpr,1,1,.false.)
call check_pointer_ (stcrc,1,1,.false.)
call check_pointer_ (namsec,1,.false.)
call check_pointer_ (indxads,1,.false.)
call check_pointer_ (indxmin,1,.false.)
call check_pointer_ (indx,1,.false.)
call check_pointer_ (stqt,1,1,.false.)
call check_pointer_ (vmin,1,1,.false.)
call check_pointer_ (propgas,1,1,.false.)
call check_pointer_ (propads,1,1,.false.)
call check_pointer_ (label,1,.false.)
call check_pointer_ (label,1,.false.)
call check_pointer_ (labelm,1,.false.)
call check_pointer_ (value,1,.false.) 
call check_pointer_ (namevalue,1,.false.)
call check_pointer_ (s2lu,1,1,.false.)
call check_pointer_ (temp,1,.false.)
!%--------------------------------------------------------
!% Deallocate local pointers corresponding to kinetic 
!% information 
!%--------------------------------------------------------
call check_pointer_ (ea,1,.false.)
call check_pointer_ (nkin,1,.false.)
call check_pointer_ (thresh,1,.false.)
call check_pointer_ (dins,1,1,.false.)
call check_pointer_ (fora,1,1,.false.)
call check_pointer_ (rk,1,1,.false.)
call check_pointer_ (pcat,1,1,1,.false.)
call check_pointer_ (ncat,1,1,.false.)
call check_pointer_ (nacat,1,1,1,.false.)
!%---------------------------------------------------------
!% Close the thermodynamic data base
!%---------------------------------------------------------
close(unit=iuthdb)
!%--------------------------------------------------------
if (iserror) goto 10 
!%--------------------------------------------------------
return
 
10 continue 
print *,'***************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: set_from_master25_'
print *, msg
print *,'***************************'
iserror=.true.
return
20 msg='Error when trying to open thermodynamic data base'
goto 10
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_from_phreeqc_pchemsys &
   (this, &
    tempisoterm, &
    iact, &
    iconvact, &
    ioptxhl, &
    naqpri, &
    nadspri, &
    naqx, &
    nmin, &
    ngas, &
    nads, &
    nmod, &
    nadsmod, &
    napri, &
    naaqx, &
    namin, &
    nagas, &
    naads, &
    naadsmod, &
    idmeq, &
    namekindb, &
	filebase, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the chemical system object according master25.dat 
! (or mastertemp.dat) thermodynamic database. 
! The kinetic information is read of kinetic.dat file. 
!
!   $Arguments:
!
 
 
type (t_parentchemicalsystem), intent(inout):: this

real*8, intent(in)                          :: tempisoterm

integer, intent(in)                         :: naqpri ! Number of aqueous primary species

integer, intent(in)                         :: nadspri ! Number of adsorption species 

integer, intent(in)                         :: naqx ! Number of aqueous complexes species  

integer, intent(in)                         :: nmin ! Number of minerals species 

integer, intent(in)                         :: ngas ! Number of gas species 

integer, intent(in)                         :: nads ! Number of adsorbed species 

integer, intent(in)                         :: nmod ! Number of surface complexes models. 

integer, intent(in)                         :: nadsmod(nmod) ! Number of species for surface complex model 

character(len=*), intent(in)                :: napri(naqpri+nadspri) 

character(len=*), intent(in)                :: naaqx(naqx)

character(len=*), intent(in)                :: namin(nmin)

character(len=*), intent(in)                :: nagas(ngas)

character(len=*), intent(in)                :: naads(nads)

character(len=*), intent(in)                :: naadsmod(nmod) ! Name of models where forming surface complexes (or cation exchange species)

integer, intent(in)                         :: idmeq(nmin)

character(len=*), intent(in)                :: namekindb     ! Name of kinetic database

character(len=*), intent(in)                :: filebase

integer, intent(in)                         :: iact ! Name of aqueous activity model 

integer, intent(in)                         :: iconvact ! Index of convention 

integer, intent(in)                         :: ioptxhl ! Name of aqueous activity model  

logical, intent(out)                        :: iserror ! If true there was error 
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
!%------------------------------------------------------------------------
return
 
10 continue 
print *,'***************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: set_from_phreeac_'
print *, msg
print *,'***************************'
iserror=.true.
return
20 msg='Error when trying to open thermodynamic data base'
goto 10
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_from_types_pchemsys &
   (this, &
    name, &
    tempref, &
    phase, &
    surface, &
    reaction, &
    numph, &
    numsurf, &
    numreact, &
    typecomponent, &
    tolunknr, &
    tolresnr, &
    deltasatmin, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set attributes in the parent chemical system object according other types (e.g. reactions, phases, etc.)
!
!   $Arguments:
!
 
 
type (t_parentchemicalsystem), intent(inout)        :: this               ! Type parent chemical system variable

integer, intent(in)                                 :: numph              ! number of phases objects

integer, intent(in)                                 :: numsurf            ! number of surfaces objects

integer, intent(in)                                 :: numreact           ! number of reactions objects

type (t_phase), intent(in), dimension(numph)        :: phase              ! vector of phase type objects

type (t_surface), intent(in), dimension(numsurf)    :: surface            ! vector of surface type objects
 
type (t_reaction), intent(in), dimension(numreact)  :: reaction           ! vector of reaction type objects

character(len=100), intent(in)                      :: name               ! Name of parent chemical system 

logical, intent(out)                                :: iserror            ! If .true. there was error

character(len=*), intent(in)                        :: typecomponent      ! Method to compute matrix of components 

real*8, intent(in)                                  :: tempref            ! Temperature (only in case of isoterm chemical system)

real*8, intent(in)                                  :: tolunknr           ! Tolerance in unknowns for Newton-Raphson method

real*8, intent(in)                                  :: tolresnr           ! Tolerance in residual for Newton-Raphson method

real*8, intent(in)                                  :: deltasatmin        ! Tolerance when considering saturation of minerals
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                       :: &
 i, &
 sum, &
 isum
logical                       :: &
 iskin
character(len=100)            :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------------
this%name=name
this%tempref=tempref
!%------------------------------------------------------------------
!% Set the tolerance in unknowns 
!%------------------------------------------------------------------
if (tolunknr/=0.0d0) this%tolunknr=dabs(tolunknr)
!%----------------------------------------------------------
!% Set the tolerance in residual 
!%----------------------------------------------------------
if (tolresnr/=0.0d0) this%tolresnr=dabs(tolresnr)
!%----------------------------------------------------------
!% Set delta saturation indices considered 
!%----------------------------------------------------------
this%deltasatmin=dabs(deltasatmin)
!%----------------------------------------------------------
!% Set reactions
!%----------------------------------------------------------
sum=0 
this%numreact=numreact
if (this%numreact>0) then
 allocate(this%preaction(this%numreact))
 call check_pointer_ (this%iskinreact,this%numreact,.true.)
 call check_pointer_ (this%idreactsp,this%numreact,.true.)
 do i=1,this%numreact
   allocate (this%preaction(i)%ptr)
   call create_ (this%preaction(i)%ptr)
   call get_if_kinetic_ (reaction(i),iskin)
   this%preaction(i)%ptr=reaction(i)
   this%iskinreact(i)=iskin
   if (iskin) sum=sum+1
 end do
end if
!%-------------------------------------------------------------
!% Set number of reaction rate laws
!%-------------------------------------------------------------
this%numrrlaw=sum
!%-------------------------------------------------------------
if (this%numrrlaw>0) then
 sum=0
 allocate (this%prrlaw(this%numrrlaw))
 do i=1,this%numreact
  if (this%iskinreact(i)) then
    sum=sum+1
    allocate(this%prrlaw(sum)%ptr)
    call get_prrlaw_ (this%preaction(i)%ptr,this%prrlaw(sum)%ptr)
  end if
 end do
end if
!%--------------------------------------------------------------
!% Set phases
!%--------------------------------------------------------------
this%numph=numph
allocate(this%pphase(this%numph))
do i=1,this%numph
 allocate (this%pphase(i)%ptr)
 call create_ (this%pphase(i)%ptr)
 this%pphase(i)%ptr=phase(i)
end do
!%--------------------------------------------------------------
!% Set surfaces
!%--------------------------------------------------------------
sum=0
this%numsurf=numsurf
if (this%numsurf>0) then
 allocate(this%psurface(this%numsurf))
 do i=1,this%numsurf
   allocate (this%psurface(i)%ptr)
   call create_ (this%psurface(i)%ptr)
   this%psurface(i)%ptr=surface(i)
   call get_num_xoh_ (this%psurface(i)%ptr,isum)
   sum=sum+isum
 end do
 this%numsites = sum
end if
!%--------------------------------------------------------------
!% Check chemical system
!%--------------------------------------------------------------
call check_(this,msg,iserror)
if (iserror) goto 10
!%--------------------------------------------------------------
! Build the chemical base
!%--------------------------------------------------------------
call build_chemical_base_(this,msg,iserror)
if (iserror) goto 10
!%--------------------------------------------------------------
!% Build the stoichiometric matrix
!%--------------------------------------------------------------
call build_stq_ (this,msg,iserror)
if (iserror) goto 10
!%--------------------------------------------------------------
!% Compute components matrix
!%--------------------------------------------------------------
select case (typecomponent)
case ('gj','GJ','') ! Gauss-Jordan
 this%isgaussjordan=.true.
case ('sv','SV') ! Singular values descomposition
 this%isgaussjordan=.false.
case default
 msg='Error, not recognized the method for compute the components matrix:'
 call add_ (msg,typecomponent)
 goto 10
end select
!%----------------------------------------------------------------
!% Build the components matrix
!%----------------------------------------------------------------
if (this%isgaussjordan) then
  call build_components_matrix_gauss_jordan_(this,msg,iserror)
else
  call build_components_matrix_sing_values_(this,msg,iserror)
end if
if (iserror) goto 10
!%----------------------------------------------------------------
!% if the problem is isoterm update parameters that depend of the 
!% temperature
!%----------------------------------------------------------------
call update_ (this,this%tempref,iserror)
if (iserror) then
   msg='Error when calling update_'
   goto 10
end if
!%----------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Chemical System:'
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
subroutine set_from_data_base_info_pchemsys &
   (this, &
    tempisoterm, &
    iact, &
    iconv, &
    naqpri, &
    nadspri, &
    naqx, &
    nmin, &
    nads, &
    ngas, &
    ntemp, &
    nreac, &
    ncoeff, &
    nmod, &
    mxkinterm, &
    temp, &
    napri, &
    naaqx, &
    namin, &
    naads, &
    nagas, &
    naadsmod, &
    nadsmod, &
    idmeq, &
    a0t, &
    zt, &
    wmolt, &
    bdot, &
    vmin, &
    propgas, &
    propads, &
    logk, &
    stqt, &
    ea, &
    nkin, &
    thresh, &
    dins, &
    fora, &
    rk, &
	pcat, &
	ncat, &
	nacat, &
	iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the parent chemical system
!
!   $Arguments:
!
 
 
type (t_parentchemicalsystem), intent(inout)                            :: this

real*8, intent(in)                                                      :: tempisoterm

integer, intent(in)                                                     :: naqpri

integer, intent(in)                                                     :: nadspri

integer, intent(in)                                                     :: naqx

integer, intent(in)                                                     :: nmin

integer, intent(in)                                                     :: nads

integer, intent(in)                                                     :: ngas

integer, intent(in)                                                     :: ntemp

integer, intent(in)                                                     :: nreac

integer, intent(in)                                                     :: ncoeff

integer, intent(in)                                                     :: nmod 

integer, intent(in)                                                     :: mxkinterm 

real*8, intent(in), dimension(ntemp)                                    :: temp

character(len=*), intent(in), dimension(naqpri+nadspri)                 :: napri

character(len=*), intent(in), dimension(naqx)                           :: naaqx

character(len=*), intent(in), dimension(nmin)                           :: namin

character(len=*), intent(in), dimension(nmod)                           :: naadsmod
    
integer, intent(in), dimension(nmod)                                    :: nadsmod ! Number of adsorbed species by model 

character(len=*), intent(in), dimension(nads)                           :: naads

real*8, intent(in), dimension(ngas)                                     :: propgas

real*8, intent(in), dimension(nads)                                     :: propads

character(len=*), intent(in), dimension(ngas)                           :: nagas

integer, intent(in), dimension(nmin)                                    :: idmeq

real*8, intent(in), dimension(naqpri+nadspri+naqx)                      :: a0t

real*8, intent(in), dimension(naqpri+nadspri+naqx)                      :: zt

real*8, intent(in), dimension(naqpri+nadspri+naqx)                      :: wmolt

real*8, intent(in), dimension(naqpri+nadspri+naqx)                      :: bdot
    
real*8, intent(in), dimension(nmin)                                     :: vmin

real*8, intent(in), dimension(nreac,ntemp)                              :: logk
    
real*8, intent(in), dimension(nreac,naqpri+nadspri+naqx+nmin+ngas+nads) :: stqt

real*8, intent(in), dimension(nmin)                                     :: ea

integer, intent(in), dimension(nmin)                                    :: nkin

real*8, intent(in), dimension(nmin)                                     :: thresh

real*8, intent(in), dimension(nmin,mxkinterm)                           :: dins

real*8, intent(in), dimension(nmin,mxkinterm)                           :: fora

real*8, intent(in), dimension(nmin,mxkinterm)                           :: rk

real*8, intent(in), dimension(nmin,mxkinterm,naqpri+naqx+nmin)          :: pcat

character(len=*), intent(in), dimension(nmin,mxkinterm,naqpri+naqx+nmin):: nacat

integer, intent(in), dimension(nmin,mxkinterm)                          :: ncat

integer, intent(in)                                                     :: iact

integer, intent(in)                                                     :: iconv
 
logical, intent(out)                                                    :: iserror

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)       :: &
 msg, &
 namecomp, &
 namechemsys, &
 nameph, &
 namephtype, &
 nameactmodel, &
 nameconvention, &
 namekin 
logical                  :: &
 iskin, &
 isareadep  
integer                  :: &
 nsp, &
 nph, &
 nsurf, &
 ireact, &
 isps, &
 isps1, &
 i, &
 j, &
 k, &
 iph, &
 nprop, &
 ipri, &
 nterm, &
 ikin, &
 itypekin, &
 nattr
integer, pointer          :: &
 idxoh(:) => null ()   
type(t_species), pointer  :: &
 species(:) => null ()
type(t_phase), pointer  :: &
 phase(:) => null ()
type(t_surface), pointer  :: &
 surface(:) => null ()
type(t_reaction), pointer  :: &
 reaction(:) => null ()
type(t_reactionratelaw), pointer  :: &
 prrlaw
real*8, pointer          :: &
 stqloc(:) => null (), &
 prop(:) => null (), &
 propsite(:,:) => null (), &
 attrterm(:,:) => null (), &
 attrspterm(:,:) => null ()
character(len=100), pointer :: &
 namesp(:) => null (), &
 nameprop(:) => null (), &
 unitprop(:) => null (), &
 modelprop(:) => null (), &
 typeterm(:) => null (), &
 namespterm(:,:) => null ()
real*8, parameter        :: &
 deltasatmin=1.d-3, &       ! Saturation index should be grater that 1+deltasatmin to consider the mineral
 zero=1.0d-10
integer, parameter        :: &
 idealmodel=1, &
 davismodel=2, &
 bdotmodel=3, &
 pitzermodel=4, &
 mxdim=100 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg='' 
!%-------------------------------------------------------------------------
!% Check number of reactions
!%-------------------------------------------------------------------------
if (nreac/=naqx+nmin+ngas+nads) then
 msg='Error, number of reactions different to number of secondary species:'
 call add_ (msg,nreac)
 goto 10
end if  
!%-------------------------------------------------------------------------
!% Initialice local variables 
!%-------------------------------------------------------------------------
nsurf=0 
namecomp='gj'
namechemsys='Chemical System'
!%-------------------------------------------------------------------------
!% Allocate local pointers 
!%-------------------------------------------------------------------------
call check_pointer_ (namesp,mxdim,.true.)
call check_pointer_ (stqloc,mxdim,.true.)
!%-------------------------------------------------------------------------
!% Allocate and create reaction objects 
!%-------------------------------------------------------------------------
allocate (reaction(nreac))
do i=1,nreac
 call create_ (reaction(i))
end do 
ireact=0 
ikin=0 
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Homogeneous aqueous reactions 
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
do i=1,naqx
 ireact=ireact+1
 isps=1
 stqloc(isps)=-1.0d0
 namesp(isps)=naaqx(ireact) 
 
 do j=1,naqpri
  if (dabs(stqt(ireact,j))>zero) then
   isps=isps+1
   namesp(isps)=napri(j)
   stqloc(isps)=stqt(ireact,j)
  end if
 end do
!%-------------------------------------------------------------------------
!% Set in the reaction object 
!%------------------------------------------------------------------------- 
 call set_ &
  (reaction(ireact), &
   naaqx(ireact), &
   .false., &
   stqloc(1:isps), &
   namesp(1:isps), &
   logk(ireact,1:ntemp), &
   temp(1:ntemp), &
   ntemp, &
   ncoeff, &
   isps, &
   iserror)
end do
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Heterogeneous reaction objects (minerals)
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
nsp=naqpri+naqx+nmin
do i=1,nmin
 ireact=ireact+1
 isps=1
 stqloc(isps)=-1.0d0
 namesp(isps)=namin(i) 
 if (idmeq(i)/=0) then
  if (idmeq(i)==2.or.idmeq(i)==3) then
     isareadep=.false. 
  else
     isareadep=.true. 
  end if 
  ikin=ikin+1
  nterm=nkin(i)
  call check_pointer_ (attrterm,3,nterm,.true.)
  call check_pointer_ (typeterm,nterm,.true.)
  call check_pointer_ (namespterm,nsp,nterm,.true.)
  call check_pointer_ (attrspterm,nsp,nterm,.true.)
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
  if (idmeq(i)==1.or.idmeq(i)==2) then  ! lasaga type 
    itypekin=1
    nattr=3 
	do j=1,nterm
         isps1=0 
	     typeterm(j)='catalitic'  
         attrterm(1,j)=rk(i,j)
         attrterm(2,j)=dins(i,j)
         attrterm(3,j)=fora(i,j)
         do k=1,nsp
           if (nacat(i,j,k)/='') then
             isps1=isps1+1 
             namespterm(isps1,j)=nacat(i,j,k)
             attrspterm(isps1,j)=pcat(i,j,k)
           end if 
         end do  
    end do 
  
  else if (idmeq(i)==3) then   ! monod type
     
	 itypekin=2
	 nattr=1 
   
     do j=1,nterm
        
		  isps1=0 
		  
		  if (rk(i,j)>0.0d0) then
		   typeterm(j)='monod'  
		  else if (rk(i,j)<0.0d0) then
           typeterm(j)='inhibitor'  
          else
		   typeterm(j)='nonmonod'  
		  end if  			

          do k=1,nsp
            if (nacat(i,j,k)/='') then
               isps1=isps1+1 
               namespterm(isps1,j)=nacat(i,j,k)
               attrspterm(isps1,j)=pcat(i,j,k)
            end if 
          end do  

	 end do 
  
  end if 
!%-------------------------------------------------------------------------
!% Create and set the reaction rate law object 
!%-------------------------------------------------------------------------
  namekin='reaction rate law'
  call add_ (namekin,ikin)
  allocate (prrlaw)
  call create_ (prrlaw)
  call set_ &
  (prrlaw,&
   itypekin, &
   namekin, &
   ncat(i,1:nterm), &
   attrterm, &
   attrspterm, &
   namespterm, &
   typeterm, &
   nterm, &
   nattr, &
   nsp, &
   ea(i), &
   thresh(i), &
   0.0d0, &
   dins(i,1), &
   fora(i,1), &
   isareadep, &
   iserror) 
  iskin=.true.
  if (iserror) goto 20 
 else
  iskin=.false.
 end if 
 
 do j=1,naqpri
  if (dabs(stqt(ireact,j))>zero) then
   isps=isps+1
   namesp(isps)=napri(j)
   stqloc(isps)=stqt(ireact,j)
  end if
 end do
!%-------------------------------------------------------------------------
!% Set in the reaction object 
!%-------------------------------------------------------------------------
 call set_ &
  (reaction(ireact), &
   namin(i), &
   iskin, &
   stqloc(1:isps), &
   namesp(1:isps), &
   logk(ireact,1:ntemp), &
   temp(1:ntemp), &
   ntemp, &
   ncoeff, &
   isps, &
   iserror, &
   rrlaw=prrlaw)
 if (iskin) then
  call destroy_ (prrlaw)
  deallocate (prrlaw)
 end if 
 if (iserror) goto 20 
end do
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Heterogeneous reaction objects (gas)
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
do i=1,ngas
 ireact=ireact+1
 isps=1
 stqloc(isps)=-1.0d0
 namesp(isps)=nagas(i) 
 
 do j=1,naqpri
  if (dabs(stqt(ireact,j))>zero) then
   isps=isps+1
   namesp(isps)=napri(j)
   stqloc(isps)=stqt(ireact,j)
  end if
 end do
!%-------------------------------------------------------------------------
!% Set in the reaction object 
!%------------------------------------------------------------------------- 
 call set_ &
  (reaction(ireact), &
   nagas(i), &
   .false., &
   stqloc(1:isps), &
   namesp(1:isps), &
   logk(ireact,1:ntemp), &
   temp(1:ntemp), &
   ntemp, &
   ncoeff, &
   isps, &
   iserror)
end do
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Heterogeneous reaction objects (adsorption and cation exchange)
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
isps1=0
do i=1,nmod
 do j=1,nadsmod(i)
  ireact=ireact+1
  isps=1
  isps1=isps1+1
  stqloc(isps)=-1.0d0
  namesp(isps)=naads(isps1) 
  do k=1,naqpri+nadspri 
   if (dabs(stqt(ireact,k))>zero) then
    isps=isps+1
    namesp(isps)=napri(k)
    stqloc(isps)=stqt(ireact,k)
   end if
  end do
!%-------------------------------------------------------------------------
!% Set in the reaction object 
!%-------------------------------------------------------------------------  
  call set_ &
   (reaction(ireact), &
    naads(isps1), &
    .false., &
    stqloc(1:isps), &
    namesp(1:isps), &
    logk(ireact,1:ntemp), &
    temp(1:ntemp), &
    ntemp, &
    ncoeff, &
    isps, &
    iserror)
  end do 
end do
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Set species objects 
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Compute the total number of species objects 
!%-------------------------------------------------------------------------
nsp=naqpri+nadspri+naqx+nmin+ngas+nads
!%-------------------------------------------------------------------------
!% Allocate and create the species objects 
!%-------------------------------------------------------------------------
allocate (species(nsp))
do i=1,nsp
 call create_ (species(i))
end do
!%-------------------------------------------------------------------------
!% Aqueous species 
!%-------------------------------------------------------------------------
isps=0 
nprop=4
call check_pointer_ (nameprop,nprop,.true.)
call check_pointer_ (unitprop,nprop,.true.)
call check_pointer_ (prop,nprop,.true.)
nameprop(1)='charge'
unitprop(1)=' ' 
nameprop(2)='ionsize'
unitprop(2)='mum'
nameprop(3)='molweight'
unitprop(3)='gr/mol'
nameprop(4)='bdot'
unitprop(4)=' '
do i=1,naqpri
 isps=isps+1
 prop(1)=zt(isps)
 prop(2)=a0t(isps)
!%-------------------------------------------------------------------------
!% Change unit Kgr/mol => gr/mol 
!%-------------------------------------------------------------------------
 prop(3)=wmolt(isps)*1.0d3 
!%-------------------------------------------------------------------------
!%
!%-------------------------------------------------------------------------
 prop(4)=bdot(isps)
 call set_name_ (species(isps),napri(i))
 call set_prop_ (species(isps),prop,nameprop,unitprop,nprop) 
end do 
do i=1,naqx
 isps=isps+1
 prop(1)=zt(nadspri+isps)
 prop(2)=a0t(nadspri+isps)
!%-------------------------------------------------------------------------
!% Change unit Kgr/mol => gr/mol 
!%-------------------------------------------------------------------------
 prop(3)=wmolt(nadspri+isps)*1.0d3
!%-------------------------------------------------------------------------
!% 
!%-------------------------------------------------------------------------
 prop(4)=bdot(nadspri+isps)
 call set_name_ (species(isps),naaqx(i))
 call set_prop_ (species(isps),prop,nameprop,unitprop,nprop) 
end do 
!%-------------------------------------------------------------------------
!% Mineral species objects 
!%-------------------------------------------------------------------------
call check_pointer_ (nameprop,nprop,.true.)
call check_pointer_ (unitprop,nprop,.true.)
call check_pointer_ (prop,nprop,.true.)
nameprop(1)='molvol'
unitprop(1)='cm3/mol'
nprop=1
do i=1,nmin
 isps=isps+1
 prop(1)=vmin(i)
 call set_name_ (species(isps),namin(i))
 call set_prop_ (species(isps),prop(1:nprop),nameprop(1:nprop),unitprop(1:nprop),nprop) 
end do 
!%-------------------------------------------------------------------------
!% Gas species objects 
!%-------------------------------------------------------------------------
nprop=1
call check_pointer_ (nameprop,nprop,.true.)
call check_pointer_ (unitprop,nprop,.true.)
call check_pointer_ (prop,nprop,.true.)
nameprop(1)='molvol'
unitprop(1)='??????'
do i=1,ngas
 isps=isps+1
 prop(1)=propgas(i)
 call set_name_ (species(isps),nagas(i))
 call set_prop_ (species(isps),prop(1:nprop),nameprop(1:nprop),unitprop(1:nprop),nprop) 
end do 
!%-------------------------------------------------------------------------
!% Adsorbed species objects 
!%-------------------------------------------------------------------------
isps1=0 
do i=1,nmod
 select case (naadsmod(i))
 case ('ex','EX')
  nameprop(1)='chargecexch'
  unitprop(1)=' '
 case ('dl','cc','tl','DL','CC','TL')
  nameprop(1)='charge'
  unitprop(1)=' '
 end select 
 isps=isps+1
 prop(1)=zt(naqpri+i) 
 call create_ (species(isps))
 call set_name_ (species(isps),napri(naqpri+i))
 call set_prop_ (species(isps),prop(1:nprop),nameprop(1:nprop),unitprop(1:nprop),nprop) 
 do j=1,nadsmod(i)
  isps=isps+1
  isps1=isps1+1
  prop(1)=propads(isps1)
  call set_name_ (species(isps),naads(isps1))
  call set_prop_ (species(isps),prop(1:nprop),nameprop(1:nprop),unitprop(1:nprop),nprop)  
 end do
end do
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Build the Phases objects 
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Compute the number of phases 
!%-------------------------------------------------------------------------
nph=1+nmin+ngas
!%-------------------------------------------------------------------------
!% Allocate and create the phase objects 
!%-------------------------------------------------------------------------
allocate (phase(nph))
do i=1,nph
 call create_ (phase(i))
end do
!%-------------------------------------------------------------------------
!% Create and set the aqueous phase object
!%-------------------------------------------------------------------------
nameph='aqueous phase'
namephtype='aqueous'
nameactmodel=''
nameconvention='' 
!%-------------------------------------------------------------------------
!% Select the thermodynamic model (Davis, Pitzer, etc)
!%-------------------------------------------------------------------------
select case (iact)
case (idealmodel)
 nameactmodel='ideal'
case (davismodel)
 nameactmodel='davis'
case (bdotmodel)
 nameactmodel='bdot'
case (pitzermodel)
 nameactmodel='pitzer'
 if (iconv/=0) nameconvention='macinnes'
end select 

iph=1
call set_ &
   (phase(iph), &
    species(1:naqpri+naqx), &
    naqpri+naqx, &
    0, &
    nameph, &
    namephtype, &
    nameactmodel, &
    nameprop, &
    modelprop, &
    prop, &
    nameconvention, &
    iserror) 
if (iserror) goto 20  
!%-------------------------------------------------------------------------
!% Mineral phase objects 
!%-------------------------------------------------------------------------
namephtype='mineral'
nameactmodel='ideal'
nameconvention='' 
isps=naqpri+naqx
do i=1,nmin 
 iph=iph+1
 call set_ &
   (phase(iph), &
    species(isps+i:isps+i), &
    1, &
    0, &
    namin(i), &
    namephtype, &
    nameactmodel, &
    nameprop, &
    modelprop, &
    prop, &
    nameconvention, &
    iserror) 
  if (iserror) goto 20 
 end do 
!%-------------------------------------------------------------------------
!% Gas phase objects 
!%-------------------------------------------------------------------------
namephtype='gas'
nameactmodel='ideal'
nameconvention='' 
isps=naqpri+naqx+nmin 
do i=1,ngas
 iph=iph+1
 call set_ &
   (phase(iph), &
    species(isps+i:isps+i), &
    1, &
    0, &
    nagas(i), &
    namephtype, &
    nameactmodel, &
    nameprop, &
    modelprop, &
    prop, &
    nameconvention, &
    iserror) 
  if (iserror) goto 20 
 end do
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Surface objects 
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
nsurf=nmod
if (nsurf>0) then 
 allocate (surface(nsurf))
 do i=1,nsurf
  call create_ (surface(i))
 end do
 call check_pointer_ (idxoh,nmod,.true.)
 call check_pointer_ (propsite,nmod,nmod,.true.)
 idxoh=1.0d0
 !%-------------------------------
 !% Compute pointer to species 
 !%-------------------------------
 isps1=naqpri+naqx+nmin+ngas
 do i=1,nsurf
  nameph='surface'
  call add_ (nameph,i)
  isps=isps1+1
  isps1=isps1+1+nadsmod(i)
  call set_ &
   (surface(i), &
    nameph, &
    naadsmod(i), &
    species(isps:isps1), &
    nadsmod(i)+1, &
    1, &
    nadsmod(i:i)+1, &
    nameprop, &
    propsite, &
    0, &
    idxoh(i:i), &
    iserror) 
   if (iserror) goto 20 
  end do  
end if 
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Set the chemical system object according types 
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
call set_ &
   (this, &
    namechemsys, &
    tempisoterm, &
    phase, &
    surface, &
    reaction, &
    nph, &
    nsurf, &
    nreac, &
    namecomp, &
    tolunknr, &
    tolresnr, &
    deltasatmin, &
    iserror)
!%-------------------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------------------
!% Destroy species objects
!%-------------------------------------------------------------------------
do i=1,nsp
 call destroy_ (species(i))
end do
deallocate(species)
species => null ()
!%-------------------------------------------------------------------------
!% Destroy phase objects
!%-------------------------------------------------------------------------
do i=1,nph
 call destroy_ (phase(i))
end do
deallocate(phase)
phase => null ()
!%-------------------------------------------------------------------------
!% Destroy surface objects
!%-------------------------------------------------------------------------
if (nmod>0) then 
 do i=1,nmod
  call destroy_ (surface(i))
 end do
 deallocate(surface)
 surface => null ()
end if 
!%-------------------------------------------------------------------------
!% Destroy reaction objects
!%-------------------------------------------------------------------------
do i=1,nreac
 call destroy_ (reaction(i))
end do
deallocate(reaction)
reaction => null ()
!%-------------------------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------------------------
call check_pointer_ (nameprop,1,.false.)
call check_pointer_ (unitprop,1,.false.)
call check_pointer_ (modelprop,1,.false.)
call check_pointer_ (prop,1,.false.)
call check_pointer_ (idxoh,1,.false.)
call check_pointer_ (propsite,1,1,.false.)
call check_pointer_ (attrterm,1,1,.false.)
call check_pointer_ (typeterm,1,.false.)
call check_pointer_ (namespterm,1,1,.false.)
call check_pointer_ (attrspterm,1,1,.false.)
if (iserror) goto 10 
!%-------------------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Chemical System:'
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
subroutine select_cmob_pchemsys &
   (this, &
    cmob, &
    nmobph, &
    nsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Select mobile concentrations from c vector and storage in cmob
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

integer, intent(in)                       :: nsp ! Number of species 

integer, intent(out)                      :: nmobph ! Number of mobile phases 

logical, intent(out)                      :: iserror

real*8, intent(in), dimension(nsp)        :: c ! Concentrations vector 

real*8, pointer, dimension(:,:)           :: cmob ! Concentrations vector of mobile species 
 
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
 isps1, &
 isps2, &
 i
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%----------------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%----------------------------------------------------------------
!% Compute the number of mobile phases 
!%----------------------------------------------------------------
nmobph=1+this%numgasph
call check_pointer_ (cmob,nsp,nmobph,.true.)
!%----------------------------------------------------------------
!% For aqueous phase
!%----------------------------------------------------------------
call get_iposspsph_ (this, this%aqphindex,isps1,isps2)
cmob(isps1:isps2,1)=c(isps1:isps2)
!%----------------------------------------------------------------
!% For gas phases
!%----------------------------------------------------------------
do i=1,this%numgasph
 call get_iposspsph_ (this, this%idgasph(i),isps1,isps2)
 cmob(isps1:isps2,1+i)=c(isps1:isps2)
end do
!%----------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%name
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
subroutine select_caq_pchemsys &
   (this, &
    caq, &
    nsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Select the concentrations of aqueous species 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

integer, intent(in)                       :: nsp ! Number of species 

logical, intent(out)                      :: iserror

real*8, intent(in), dimension(nsp)        :: c ! Concentrations vector 

real*8, pointer, dimension(:)             :: caq ! Concentrations vector of aqueous species
 
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
 isps1, &
 isps2, &
 i
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!


!%----------------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%----------------------------------------------------------------
call check_pointer_ (caq,nsp,.true.)
!%----------------------------------------------------------------
!% For aqueous phase
!%----------------------------------------------------------------
call get_iposspsph_ (this,this%aqphindex,isps1,isps2)
caq(isps1:isps2)=c(isps1:isps2)
!%----------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: select_caq_'
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
subroutine select_cads_pchemsys &
   (this, &
    cads, &
    nsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Select adsorbed concentrations from c vector and storage in cads
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in):: this

integer, intent(in)                      :: nsp

logical, intent(out)                     :: iserror

real*8, intent(in), dimension(nsp)       :: c

real*8, pointer, dimension(:)            :: cads 
 
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
 isps1, &
 isps2, &
 isurf
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%----------------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%----------------------------------------------------------------
!% Allocate pointer 
!%----------------------------------------------------------------
call check_pointer_ (cads,this%numsp,.true.)
!%----------------------------------------------------------------
!% For surfaces
!%----------------------------------------------------------------
do isurf=1,this%numsurf
 call get_iposspsurf_ (this,isurf,isps1,isps2)
 cads(isps1:isps2)=c(isps1:isps2)
end do
!%----------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: select_cads_'
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
subroutine select_cmin_pchemsys &
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
 
type (t_parentchemicalsystem), intent(in) :: this      ! Type parent chemical system variable

real*8, pointer, dimension(:)             :: cmin      ! Concentration of mineral species 

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
integer                      :: &
 isps1, &
 isps2, &
 iph 
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%----------------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------------
!% Check the number of species 
!%----------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%----------------------------------------------------------------
call check_pointer_ (cmin,this%numsp,.true.)
!%----------------------------------------------------------------
!% For mineral
!%----------------------------------------------------------------
do iph=1,this%numminph
 call get_iposspsph_ (this,this%idminph(iph),isps1,isps2)
 cmin(isps1:isps2)=c(isps1:isps2)
end do
!%----------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: select_cmin_'
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
subroutine select_cgas_pchemsys &
   (this, &
    cgas, &
    nsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Select the concentrations of gas species
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

real*8, pointer, dimension(:)             :: cgas

integer, intent(in)                       :: nsp

real*8, intent(in), dimension(nsp)        :: c

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
integer                      :: &
 isps1, &
 isps2, &
 iph
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%----------------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------------
!% Check the number of species 
!%----------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%----------------------------------------------------------------
!% Allocate pointer 
!%----------------------------------------------------------------
call check_pointer_ (cgas,this%numsp,.true.)
!%----------------------------------------------------------------
!% For gases
!%----------------------------------------------------------------
do iph=1,this%numgasph
 call get_iposspsph_ (this,this%idgasph(iph),isps1,isps2)
 cgas(isps1:isps2)=c(isps1:isps2)
end do
!%----------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: select_cgas_'
print *, msg
print *,'*************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_pchemsys &
   (this, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Check chemical system consistence 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout)  :: this       ! Type parent chemical system variable. 

character(len=*), intent(out)                 :: msg        ! Message error. 

logical, intent(out)                          :: iserror    ! iserror=true, then there was an error. 
 
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
 sum, &
 i, &
 j, &
 iph, &
 isurf, &
 ireact, &
 isp1, &
 isp2, &
 numsplocal, &
 numminph, &
 numgasph, &
 numsurfph, &
 nreact
integer, pointer               :: &
 idminphloc(:) => null (), &
 idgasphloc(:) => null (), &
 idreact(:) => null ()
character(len=100), pointer    :: &
 namesp(:) => null ()
character(len=100)             :: &
 name1, &
 name2
logical                        :: &
 be, &
 besp 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
call check_pointer_ (idminphloc,this%numph,.true.)
call check_pointer_ (idgasphloc,this%numph,.true.)
sum = 0
iserror = .false.
numminph=0
numgasph=0
numsurfph=0
!%----------------------------------------------------------
!% Check if the phase is repeated and compute number of local 
!% chemical species
!%----------------------------------------------------------
do i=1,this%numph
 iph=0
 call get_name_ (this%pphase(i)%ptr,name1)
 do j=1,this%numph
   call get_name_ (this%pphase(j)%ptr,name2)
   if (name1==name2) then
    iph=iph+1
   end if
   if (iph>1) then
    msg='Error, the phase is repeated:'
    call add_ (msg,name1)
    goto 10
   end if
 end do
 call get_numsp_ (this%pphase(i)%ptr,numsplocal)
 call get_if_aqueous_ (this%pphase(i)%ptr,be)
 if (be) this%aqphindex=i
 call get_if_mineral_ (this%pphase(i)%ptr,be)
 if (be) then
  numminph=numminph+1
  idminphloc(numminph)=i
 end if
 call get_if_gas_ (this%pphase(i)%ptr,be)
 if (be) then
  numgasph=numgasph+1
  idgasphloc(numgasph)=i
 end if
 
 sum = sum + numsplocal
 
end do
if (this%aqphindex==0) then
 msg='Error, not defined aqueous phase in the chemical system'
 goto 20
end if
!%---------------------------------------------------------
do i=1,this%numsurf
!%---------------
 isurf=0
 call get_name_ (this%psurface(i)%ptr,name1)
 do j=1,this%numsurf
   call get_name_ (this%psurface(j)%ptr,name2)
   if (name1==name2) then
    isurf=isurf+1
   end if
   if (isurf>1) then
    msg='Error, the surface is repeated:'
    call add_ (msg,name1)
    goto 10
   end if
 end do
!%---------------
 call get_numsp_ (this%psurface(i)%ptr,numsplocal)
 sum = sum + numsplocal
end do
!%----------------------------------------------------------
!% Assign local number of chemical species
!%----------------------------------------------------------
this%numsp = sum
this%numminph = numminph
this%numgasph = numgasph
!%----------------------------------------------------------
!% Allocate idminph
!%----------------------------------------------------------
if (numminph>0) then
 call check_pointer_ (this%idminph,numminph,.true.)
 this%idminph=idminphloc(1:numminph)
end if
!%----------------------------------------------------------
!% Allocate idgasph
!%----------------------------------------------------------
if (numgasph>0) then
 call check_pointer_ (this%idgasph,numgasph,.true.)
 this%idgasph=idgasphloc(1:numgasph)
end if
!%----------------------------------------------------------
!% Allocate pointer to species objects 
!%----------------------------------------------------------
allocate(this%pspecies(this%numsp))
sum=0
!%----------------------------------------------------------
!% Pointer species in phase objects 
!%----------------------------------------------------------
do i=1,this%numph
 
 call get_numsp_ (this%pphase(i)%ptr,numsplocal)
 
 do isp1=1,numsplocal
 
   call get_pspecies_ (this%pphase(i)%ptr,this%pspecies(sum+isp1)%ptr,isp1,iserror)
   
   if (iserror) goto 20  
    
 end do

 sum=sum+numsplocal
 
end do
!%----------------------------------------------------------
!%----------------------------------------------------------
!% Pointer species in surface objects 
!%----------------------------------------------------------
!%----------------------------------------------------------
do i=1,this%numsurf
 
 call get_numsp_ (this%psurface(i)%ptr,numsplocal)
 
 do isp1=1,numsplocal
 
   call get_pspecies_  (this%psurface(i)%ptr,this%pspecies(sum+isp1)%ptr,isp1,iserror)

   if (iserror) goto 20 

 end do
 sum=sum+numsplocal
!%----------------------------------------------------------
!% Set pointer to aqueous phase in the surfaces 
!%----------------------------------------------------------
 call set_ (this%psurface(i)%ptr,this%pphase(this%aqphindex)%ptr,iserror)
 if (iserror) then
    msg='Error when calling set_'
    call add_ (msg,i)
    goto 20
 end if
end do
!%----------------------------------------------------------
!% Check species in reactions
!%----------------------------------------------------------
sum=0
do i=1,this%numreact
 call get_namesp_ (this%preaction(i)%ptr,namesp,numsplocal)
 do isp1=1,numsplocal
   besp=.false.
   do isp2=1,this%numsp
    if(this%pspecies(isp2)%ptr%name==namesp(isp1)) then
     besp=.true.
     call set_pspecies_ (this%preaction(i)%ptr,this%pspecies(isp2)%ptr)
     exit
    end if
   end do
   if(.not.besp) then
    msg='Error, not defined in the chemical system the species:'
    call add_ (msg,namesp(isp1))
      iserror=.true.
      goto 20
   end if
 end do
end do
!%----------------------------------------------------------
!% Check species in reaction rate laws
!%----------------------------------------------------------
sum=0
do i=1,this%numrrlaw
 call get_namesp_ (this%prrlaw(i)%ptr,namesp,numsplocal)
 do isp1=1,numsplocal
   besp=.false.
   do isp2=1,this%numsp
    if(this%pspecies(isp2)%ptr%name==namesp(isp1)) then
     besp=.true.
     call set_pspecies_ (this%prrlaw(i)%ptr,this%pspecies(isp2)%ptr)
     exit
    end if
   end do
   if(.not.besp) then
      msg='Error, not defined in the chemical system the species:'
      call add_ (msg,namesp(isp1))
      iserror=.true.
      goto 20
   end if
 end do
end do
!%----------------------------------------------------------
!% Build reaction networks
!%----------------------------------------------------------
call check_pointer_ (this%ishomreact,this%numreact,this%numph+this%numsurf,.true.)
call check_pointer_ (this%ishetreact,this%numreact,this%numph+this%numsurf,.true.)
 
do i=1,this%numph
  call get_idreactionb_pchemsys &
   (this, &
    idreact, &
    nreact, &
    i, &
    0, &
    .false.)
 do ireact=1,nreact
   this%ishomreact(idreact(ireact),i)=.true.
 end do
 
 do j=1,this%numph
  if (i.ne.j) then
   call get_idreactionb_pchemsys &
   (this, &
    idreact, &
    nreact, &
    i, &
    j, &
    .false.)
    do ireact=1,nreact
     this%ishetreact(idreact(ireact),i)=.true.
     this%ishetreact(idreact(ireact),j)=.true.
    end do
  end if
 end do
 
end do
!%----------------------------------------------------------
!% For surfaces
!%----------------------------------------------------------
do i=1,this%numsurf
   call get_idreactionb_pchemsys &
   (this, &
    idreact, &
    nreact, &
    this%aqphindex, &
    i, &
    .true.)
 
    do ireact=1,nreact
     this%ishetreact(idreact(ireact),this%aqphindex)=.true.
     this%ishetreact(idreact(ireact),this%numph+i)=.true.
    end do
 
 end do
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (idminphloc,1,.false.)
call check_pointer_ (idgasphloc,1,.false.)
call check_pointer_ (idreact,1,.false.)
call check_pointer_ (namesp,1,.false.)
if (iserror) goto 10
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
subroutine build_chemical_base_pchemsys &
 (this, &
 msg, &
 iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build the chemical base for aqueous species and rewrite
!%   the stoichiomeric of reactions.
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout) :: this     ! Type parent chemical system variable. 

character(len=*), intent(out)                :: msg      ! Error message. 

logical, intent(out)                         :: iserror  ! iserror=true, then there was an error. 
 
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
 namereact, &
 nameadsmodel
logical        :: &
 be
integer        :: &
 ipri, &
 i, &
 j, &
 k, &
 n, &
 sumsp, &
 numsp, &
 isp, &
 neq, &
 nkin, &
 naqreact, &
 neqaqreact, &
 naqsp, &
 numreact, &
 ipos1, &
 ipos2, &
 iph, &
 ireact, &
 numxoh, &
 loc(1)
integer       :: &
 idsp(this%numsp)
integer, pointer             :: &
 idaqreact(:) => null(), &
 idreact(:) => null()
character(len=100), pointer  :: &
 nameaqsp(:) => null(), &
 nameprisp(:) => null()
character(len=100):: &
 namespglobal(this%numsp) 
!-------------------------------------------------------------------------
!
!   $code
!
 

 


!%---------------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------------
ipri=0
idsp=0
isp=0
namespglobal=' '
nkin=0
neq=0
neqaqreact=0
!%--------------------------------------------------------------
!% Get aqueous homogeneous reactions
!%--------------------------------------------------------------
call get_idreaction_(this,idaqreact,naqreact,this%aqphindex,0,.false.)
!%---------------------------------------------------------------
 do j=1,naqreact
  select case (this%iskinreact(idaqreact(j)))
  case (.false.)
   neqaqreact = neqaqreact + 1
  end select
 end do
 call get_numsp_ (this%pphase(this%aqphindex)%ptr,naqsp)
!%--------------------------------------------------------------
! Determine dimension of the chemical base
!%--------------------------------------------------------------
this%numaqprisp = naqsp - neqaqreact 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
 do i=1,this%numreact
  select case (this%iskinreact(i))
  case(.true.)
   nkin=nkin+1
  case default
   neq=neq+1
  end select
 end do
 call check_pointer_ (this%idaqprisp,this%numaqprisp,.true.)
 call check_pointer_ (this%idreactsp,this%numreact,.true.)
 call check_pointer_ (this%stq,this%numreact,this%numsp,.true.)
!%--------------------------------------------------------------
!% Build idreactsp vector 
!%--------------------------------------------------------------
isp=0
do i=1,this%numsp
namespglobal(i) = this%pspecies(i)%ptr%name
end do
 
do i=1,this%numreact
 
 call get_name_ (this%preaction(i)%ptr,namereact)
 
 do j=1,this%numsp
 
    where(namespglobal==namereact)
     idsp=1
    elsewhere
     idsp=0
    end where
    n=sum(idsp)
 
    if(n>0) then
      isp=isp+1
      loc=maxloc(idsp,mask=idsp>0)
      this%idreactsp(i)=loc(1)
      exit
    else
      msg='Error, not defined in the chemical system the species'
      call add_ (msg,namereact)
      goto 20
    end if
 
 end do
 
end do
!%--------------------------------------------------------------
!% Build idaqprisp vector
!%--------------------------------------------------------------
call get_namesp_ (this%pphase(this%aqphindex)%ptr,nameaqsp)
call get_iposspsph_ (this, this%aqphindex, ipos1, ipos2)
ipri=0
do i=1,naqsp
  be=.false.
  do j=1,this%numreact
    call get_name_ (this%preaction(j)%ptr,namereact)
    if(nameaqsp(i)==namereact) be=.true.
  end do
  if (.not.be) then
    ipri=ipri+1
    this%idaqprisp(ipri)= ipos1 - 1 + i
  end if
end do
!%-------------------------------------------------------------- 
!% Rewrite reactions according primary species
!%--------------------------------------------------------------
call check_pointer_ (nameprisp,this%numaqprisp,.true.)
 
do i=1,this%numaqprisp
 nameprisp(i) = this%pspecies(this%idaqprisp(i))%ptr%name
 if (nameprisp(i)=='h2o') this%wcindex=i
 if (nameprisp(i)=='e-') this%ecindex=i
 if (nameprisp(i)=='h+') this%hcindex=i
end do
!%--------------------------------------------------------------
do i=1,this%numreact
 
 call rewrite_ (this%preaction(i)%ptr,nameprisp,this%numaqprisp,iserror)
 if (iserror) then
  msg='Error when calling rewrite_'
  call add_ (msg,i)
  goto 20
 end if
 call rewrite_rrlaw_ (this%preaction(i)%ptr,nameaqsp,naqsp,iserror)
 if (iserror) then
  msg='Error when calling rewrite_rrlaw_'
  call add_ (msg,i)
  goto 20
 end if
 
end do
!%-------------------------------------------------------------- 
!% Rewrite stqds of the sorption reactions
!%--------------------------------------------------------------
do i=1,this%numsurf
 
 call get_name_xoh_ (this%psurface(i)%ptr, nameprisp, numxoh)
 call get_name_model_ (this%psurface(i)%ptr, nameadsmodel)
 call get_idreaction_(this,idreact,numreact,this%aqphindex,i,.true.)
 if (numreact==0) then
   msg='Error, not defined reactions for surface'
   call add_ (msg,i)
   goto 20
 end if
 
 do j=1,numreact
   ireact=idreact(j)
   call set_stqads_(this%preaction(ireact)%ptr,nameprisp,numxoh,nameadsmodel,iserror)
   if (iserror) then
      msg='Error when calling set_stqads_'
      call add_ (msg,i)
      goto 20
   end if
 end do
 
 
end do
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (nameprisp,1,.false.)
call check_pointer_ (nameaqsp,1,.false.)
call check_pointer_ (idreact,1,.false.)
call check_pointer_ (idaqreact,1,.false.)
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
subroutine get_sp_index_pchemsys &
   (this, &
    namesp, &
    index)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the global index of one species
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

integer, intent(out)                      :: index ! Global index of species

character(len=*), intent(in)              :: namesp ! Name of the species 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=len(namesp))                :: name 

integer                                   :: isp

character(len=100)                        :: msg

logical                                   :: iserror 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

!%----------------------------------------------------------
index=0
!%----------------------------------------------------------
do isp=1,this%numsp
 
 name=this%pspecies(isp)%ptr%name
 if (name==namesp) then
   index=isp
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
subroutine get_surf_model_pchemsys &
   (this, &
    namemodel, &
    ithsurf)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of the model of ith surface
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

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
namemodel=' '
!%------------------------------------------------------------
if (ithsurf==0.or.ithsurf>this%numsurf) return
!%------------------------------------------------------------
call get_name_ (this%psurface(ithsurf)%ptr,namemodel)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_surf_index_pchemsys &
   (this, &
    ithsurf, &
    nsites, &
    namesurf)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the surface index and number of sites from namesurf
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)           :: this

integer, intent(out)                                :: ithsurf

integer, intent(out)                                :: nsites

character(len=*), intent(in)                        :: namesurf 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=len(namesurf))                        :: &
 name
integer                                             :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
ithsurf=0
nsites=0
!%------------------------------------------------------------
do i=1,this%numsurf
 
 call get_name_ (this%psurface(i)%ptr,name)
 if (name==namesurf) then
  ithsurf=i
  call get_numsites_ (this%psurface(i)%ptr,nsites)
  exit  
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
subroutine get_eq_min_sp_index_pchemsys &
   (this, &
    ideqminsp, &
    neqminsp)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the global pointers to equilibrium mineral species
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)   :: this

integer, intent(out)                        :: neqminsp

integer, pointer, dimension(:)              :: ideqminsp
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                 :: &
 i, &
 j, &
 iph, &
 nreact
integer, pointer                        :: &
 idreact(:) => null (), &
 ideqminsp1(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!
neqminsp=0
!%-------------------------------------------------------------
if (this%numminph==0) return
!%-------------------------------------------------------------
call check_pointer_ (ideqminsp1,this%numreact,.true.)
!%-------------------------------------------------------------
do i=1,this%numminph
 iph=this%idminph(i)
 call get_idreaction_(this,idreact,nreact,this%aqphindex,iph,.false.)
 do j=1,nreact
  if(.not.this%iskinreact(idreact(j))) then
   neqminsp=neqminsp+1
   ideqminsp1(neqminsp)=this%idreactsp(idreact(j))
  end if
 end do
end do
!%--------------------------------------------------------------
call check_pointer_ (ideqminsp,neqminsp,.true.)
ideqminsp(1:neqminsp)=ideqminsp1
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (idreact,1,.false.)
call check_pointer_ (ideqminsp1,1,.false.)
!%--------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_num_eq_min_sp_pchemsys &
   (this, &
    neqminsp)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the total number of equilibrium mineral species
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)    :: this

integer, intent(out)                         :: neqminsp
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                 :: &
 i, &
 j, &
 iph, &
 nreact
integer, pointer                        :: &
 idreact(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!
neqminsp=0
!%-------------------------------------------------------------
if (this%numminph==0) return
!%-------------------------------------------------------------
do i=1,this%numminph
 iph=this%idminph(i)
 call get_idreaction_(this,idreact,nreact,this%aqphindex,iph,.false.)
 do j=1,nreact
  if(.not.this%iskinreact(idreact(j))) neqminsp=neqminsp+1
 end do
end do
!%--------------------------------------------------------------
call check_pointer_ (idreact,1,.false.)
!%--------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_pchemsys &
   (this, &
    ioutput, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write parent chemical system attributes
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)    :: this

integer, intent(in)                          :: ioutput 

logical, intent(out)                         :: iserror ! If .true. there was error 
 
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
 j
character (len=100), pointer    :: &
 namesp(:) => null ()
character (len=100)             :: &
 namepri, &
 namesec, &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%-------------------------------------------------------------
write(ioutput,*)'-----------------------------------------------'
write(ioutput,*)'        Chemical System Information           '
write(ioutput,*)'-----------------------------------------------'
!%-------------------------------------------------------------
if (this%numph>0) then
 write(ioutput,*)'------------------------------------------'// &
                '----------------'
 write(ioutput,*)'========> Phases defined in the chemical system'
 write(ioutput,*)'------------------------------------------'// &
                '----------------'
!%-------------------------------------------------------------
 do i=1,this%numph
  call write_ (this%pphase(i)%ptr,ioutput,iserror)
  if (iserror) then
   msg='Error when calling write_'
   call add_ (msg,i)
  end if
 end do
end if
!%-------------------------------------------------------------
if (this%numsurf>0) then
 write(ioutput,*)'------------------------------------------'// &
                '----------------'
 write(ioutput,*)'========> Surfaces defined in the chemical system'
 write(ioutput,*)'------------------------------------------'// &
                '----------------'
!%-------------------------------------------------------------
  do i=1,this%numsurf
   call write_ (this%psurface(i)%ptr,ioutput,iserror)
   if (iserror) then
    msg='Error when calling write_'
    call add_ (msg,i)
   end if
  end do
end if
!%-------------------------------------------------------------
call get_chem_info_ (this,msg,iserror,namesp=namesp)
!%-------------------------------------------------------------
write(ioutput,*)'------------------------------------------'// &
                '----------------'
write(ioutput,*)'========> Reaction network defined in the chemical system'
write(ioutput,*)'------------------------------------------'// &
                '----------------'
!%-------------------------------------------------------------
!% Write the stoichiometric matrix 
!%-------------------------------------------------------------
write(ioutput,*)'-----------------------------------------------'
write(ioutput,*)'          Stoichiometric matrix                '
write(ioutput,*)'-----------------------------------------------'
!%------------------------------------------------------------------
write(ioutput,2) 'species','=',(namesp(j),j=1,this%numsp)
do i=1,this%numreact
 call get_name_ (this%preaction(i)%ptr,namesec)
 write(ioutput,1) namesec,'=',(this%stq(i,j),j=1,this%numsp)
end do
!%-------------------------------------------------------------
!% Write the reactions 
!%-------------------------------------------------------------
do i=1,this%numreact
 call write_ (this%preaction(i)%ptr, ioutput, iserror)
 if (iserror) then
  msg='Error when calling write_'
  call add_ (msg,i)
 end if
end do
!%------------------------------------------------------------------
write(ioutput,*)'-----------------------------------------------'
write(ioutput,*)'          Components matrix  (equilibrium)     '
if (this%isgaussjordan) then
 write(ioutput,*)'  (Computed using Gauss-Jordan Elimination) '
else
 write(ioutput,*)'(Computed using Singular Value Descomposition)'
end if
write(ioutput,*)'-----------------------------------------------'
!%------------------------------------------------------------------
 
 call get_chem_info_ (this,msg,iserror,namesp=namesp)
 
 write(ioutput,2) 'species','=',(namesp(j),j=1,this%numsp)
 
 
do i=1,this%numaqprisp
 call get_name_ (this%pspecies(this%idaqprisp(i))%ptr,namepri)
 write(ioutput,1) namepri,'=',(this%ueq(i,j),j=1,this%numsp)
end do
!%------------------------------------------------------------------
!%------------------------------------------------------------------
!%------------------------------------------------------------------ 
write(ioutput,*)'-----------------------------------------------'
write(ioutput,*)'          Components matrix  (general)         '
if (this%isgaussjordan) then
 write(ioutput,*)'  (Computed using Gauss-Jordan Elimination) '
else
 write(ioutput,*)'(Computed using Singular Value Descomposition)'
end if
write(ioutput,*)'-----------------------------------------------'
!%------------------------------------------------------------------
 
 call get_chem_info_ (this,msg,iserror,namesp=namesp)
 
 write(ioutput,2) 'species','=',(namesp(j),j=1,this%numsp)
 
 
do i=1,this%numaqprisp
 call get_name_ (this%pspecies(this%idaqprisp(i))%ptr,namepri)
 write(ioutput,1) namepri,'=',(this%u(i,j),j=1,this%numsp)
end do
!%------------------------------------------------------------
call check_pointer_ (namesp,1,.false.)
!%------------------------------------------------------------
return
 
10 continue       
print *,'*************************'
print *,'Chemical System:'
print *,'Name:', this%name
print *,'Service: write_'
print *, msg
print *,'*************************'
iserror=.true.
return
 
1 format (a5,a3,<this%numsp>f10.1)
2 format (a10,a3,<this%numsp>a10)
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_omgwcryst_pchemsys &
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
!   $Description: Compute mass of water fixed in hydrated minerals.
! 
! omgwcryst= PMh2o * 10^-3 * sum(coeff_j_h2o * cm * omgwfree)
!
!   where
!
!   PMh2o = Molecular weight of water
!   coeff_j_h2o = stoichiometric coefficient of water in the
!                 jth mineral
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this       ! Type parent chemical system variable

integer, intent(in)                       :: nsp        ! Number of species

real*8, intent(out)                       :: omgwcryst  ! Mass of fixed water [kgw]

real*8, intent(in)                        :: omgwfree   ! Mass of free water [kgw]

real*8, intent(in), dimension(nsp)        :: c          ! Concentration of the species

logical, intent(out)                      :: iserror    ! If .true. there was error 

character(len=*), intent(out)             :: msg        ! Message error 
 
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
 isps, &
 iph, &
 isps1, &
 isps2 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Initialice variables
!%-------------------------------------------------------------
omgwcryst=0.0d0
!%-------------------------------------------------------------
!% Check the number of the species
!%-------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species,nsp:'
 call add_ (msg,nsp)
 goto 10 
end if  
!%-------------------------------------------------------------
!% if the water not was considered in the chemical system, then
!% return. 
!%-------------------------------------------------------------
if (this%wcindex==0) return

do iph=1,this%numminph
   
    call get_iposspsph_ (this,this%idminph(iph),isps1,isps2)
    
	do isps=isps1,isps2
     omgwcryst=omgwcryst+this%ueq(this%wcindex,isps)*c(isps)
    end do

end do
omgwcryst=omgwcryst*omgwfree*kgwmol
!%-------------------------------------------------------------
return
10 continue       
print *,'***************************'
print *,'Chemical System:'
print *,'Name:', this%name
print *,'Service: compute_omgwcryst_'
print *, msg
print *,'***************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_chem_info_pchemsys &
   (this, &
    msg, &
    iserror, &
    nummobph, &
    numbase, &
    numsp, &
    numsites, &
	numreact, &
	numeqreact, &
	numkinreact, &
    wcindex, &
	hcindex, &
    ecindex, &
    numph, &
    wspindex, &
    numsurf, &
    namesp, &
	nameminsp, &
    nameph, &
    namesurf, &
    idbase, &
    ideqminsp, &
    neqminsp, &
    namebase, &
    nminsp, &
    idminsp, &
    idgassp, &
    ngassp,&
    nummobsp)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return general chemical information 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)   :: this

logical, intent(out)                        :: iserror

character(len=*), intent(inout)             :: msg

integer, intent(out), optional              :: numsp        ! Number of species

integer, intent(out), optional              :: nummobph      ! Number of mobile phases. 

integer, intent(out), optional             :: nummobsp      ! Number of mobile species

integer, intent(out), optional              :: numbase

integer, intent(out), optional              :: numsites

integer, intent(out), optional              :: numreact     ! Number of reactions

integer, intent(out), optional              :: numeqreact   ! Number of equilibrium reactions

integer, intent(out), optional              :: numkinreact  ! Number of kinetic reactions

integer, intent(out), optional              :: wcindex      ! Index of water component

integer, intent(out), optional              :: hcindex      ! Index of proton component

integer, intent(out), optional              :: ecindex

integer, intent(out), optional              :: numph

integer, intent(out), optional              :: numsurf

integer, intent(out), optional              :: wspindex

integer, intent(out), optional              :: neqminsp

integer, intent(out), optional              :: nminsp

integer, intent(out), optional              :: ngassp

integer, pointer, optional, dimension(:)    :: idbase

integer, pointer, optional, dimension(:)    :: ideqminsp

integer, pointer, optional, dimension(:)    :: idminsp 

integer, pointer, optional, dimension(:)    :: idgassp

character(len=100), pointer, optional       :: namesp(:)

character(len=100), pointer, optional       :: nameph(:)

character(len=100), pointer, optional       :: namesurf(:)

character(len=100), pointer, optional       :: namebase(:)

character(len=100), pointer, optional       :: nameminsp(:)

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                   :: &
 i, &
 ndim, &
 numgassp
logical                                   :: &
 havenumsp, &
 havenummobph, &
 havenumbase, &
 havenumsites, &
 havewcindex, &
 havehcindex, &
 haveecindex, &
 havenumph, &
 havewspindex, &
 havenumsurf, &
 havenamesp, &
 havenameph, &
 havenumreact, &
 havenumeqreact, &
 havenumkinreact, &
 havenamesurf, &
 haveidbase, &
 haveideqminsp, &
 haveneqminsp, &
 havenamebase, &
 haveidminsp, &
 havenminsp, &
 havengassp, &
 haveidgassp, &
 havenameminsp, &
 havenummobsp
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=" "
!%-------------------------------------------------------------
!% Check optional arguments
!%-------------------------------------------------------------
havenumsp=present(numsp)
havenummobph=present(nummobph)
havenumbase=present(numbase)
havenumsites=present(numsites)
havewcindex=present(wcindex)
havehcindex=present(hcindex)
haveecindex=present(ecindex)
havenumph=present(numph)
havenumreact=present(numreact)
havenumeqreact=present(numeqreact)
havenumkinreact=present(numkinreact)
havewspindex=present(wspindex)
havenumsurf=present(numsurf)
havenamesp=present(namesp)
havenameph=present(nameph)
havenamesurf=present(namesurf)
haveidbase=present(idbase)
haveideqminsp=present(ideqminsp)
haveneqminsp=present(neqminsp)
havenamebase=present(namebase)
haveidminsp=present(idminsp)
havenminsp=present(nminsp)
havengassp=present(ngassp)
haveidgassp=present(idgassp)
havenameminsp=present(nameminsp)
havenummobsp=present(nummobsp)
!%-------------------------------------------------------------
!%Return the number of sepcies
!%-------------------------------------------------------------
if (havenumsp) then
 numsp=this%numsp
end if
!%-------------------------------------------------------------
!% Return the number of reactions
!%-------------------------------------------------------------
if (havenumreact) then
 numreact=this%numreact
end if
!%-------------------------------------------------------------
!% Return the number of equilibrium reactions
!%-------------------------------------------------------------
if (havenumeqreact) then
 numeqreact=0 
 do i=1,this%numreact
  if (.not.this%iskinreact(i)) numeqreact=numeqreact+1
 end do 
end if
!%-------------------------------------------------------------
!% Return the number of kinetic reactions
!%-------------------------------------------------------------
if (havenumkinreact) then
 numkinreact=0 
 do i=1,this%numreact
  if (this%iskinreact(i)) numkinreact=numkinreact+1
 end do 
end if
!%-------------------------------------------------------------
!%Return the number of mobile phases
!%-------------------------------------------------------------
if (havenummobph) then
 nummobph=1+this%numgasph
end if
!%-------------------------------------------------------------
!%Return the number of primary species 
!%-------------------------------------------------------------
if (havenumbase) then
 numbase=this%numaqprisp
end if
!%-------------------------------------------------------------
!%Return the total number of sites 
!%-------------------------------------------------------------
if (havenumsites) then
 numsites=this%numsites
end if
!%-------------------------------------------------------------
!%Return the index of water component
!%-------------------------------------------------------------
if (havewcindex) then
 wcindex=this%wcindex
end if
!%-------------------------------------------------------------
!%Return the index of proton component
!%-------------------------------------------------------------
if (havehcindex) then
 hcindex=this%hcindex
end if
!%-------------------------------------------------------------
if (haveecindex) then
 ecindex=this%ecindex
end if
!%-------------------------------------------------------------
if (havenumph) then
 numph=this%numph
end if
!%-------------------------------------------------------------
if (havewspindex) then
 call get_sp_index_(this,'h2o',wspindex)
end if
!%-------------------------------------------------------------
if (havenumsurf) then
 numsurf=this%numsurf
end if
!%-------------------------------------------------------------
!%Return the name of the species
!%-------------------------------------------------------------
if (havenamesp) then
 call check_pointer_ (namesp,this%numsp,.true.)
 do i=1,this%numsp
  namesp(i)=this%pspecies(i)%ptr%name
 end do
end if
!%-------------------------------------------------------------
!%Return the name of the base
!%-------------------------------------------------------------
if (havenamebase) then
 call check_pointer_ (namebase,this%numaqprisp,.true.)
 do i=1,this%numaqprisp
  namebase(i)=this%pspecies(this%idaqprisp(i))%ptr%name
 end do
end if
!%-------------------------------------------------------------
!%Return the name of phases 
!%-------------------------------------------------------------
if (havenameph) then
 call check_pointer_ (nameph,this%numph,.true.)
 do i=1,this%numph
  call get_name_ (this%pphase(i)%ptr,nameph(i))
 end do
end if
!%-------------------------------------------------------------
!%Return the name of surfaces 
!%-------------------------------------------------------------
if (havenamesurf) then
 call check_pointer_ (namesurf,this%numsurf,.true.)
 do i=1,this%numph
  call get_name_ (this%psurface(i)%ptr,namesurf(i))
 end do
end if
!%-------------------------------------------------------------
!%Return the global index of primary species 
!%-------------------------------------------------------------
if (haveidbase) then
 call check_pointer_ (idbase,this%numaqprisp,.true.)
 idbase=this%idaqprisp
end if
!%-------------------------------------------------------------
if (haveideqminsp) then
 call get_eq_min_sp_index_ &
 (this, &
  ideqminsp, &
  ndim)
end if
!%-------------------------------------------------------------
if (haveneqminsp) then
 call get_num_eq_min_sp_ &
 (this, &
  neqminsp)
end if
!%-------------------------------------------------------------
if (havenminsp.or.haveidminsp.or.havenameminsp) then
 call get_min_sp_index_ &
 (this, &
  iserror, &
  nminsp=nminsp, &
  idminsp=idminsp, &
  nameminsp=nameminsp)
end if
!%-------------------------------------------------------------
if (havengassp.or.haveidgassp) then
 call get_gas_sp_index_ &
 (this, &
  iserror, &
  ngassp=ngassp, &
  idgassp=idgassp)
end if
!%-------------------------------------------------------------
!%Return the number of mobile species
!%-------------------------------------------------------------
if (havenummobsp) then

    !1st we get the number of aqueous species
    call get_numsp_(this%pphase(this%aqphindex)%ptr,nummobsp)

    !Then we get the gaseous spcecies
    do i=1,this%numgasph
        call get_numsp_(this%pphase(this%idgasph(i))%ptr,numgassp)
        nummobsp=nummobsp+numgassp
    enddo

 
if (havenummobph) then
 nummobph=1+this%numgasph
end if



end if
!%-------------------------------------------------------------





return
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_iumob_pchemsys &
   (this, &
    iumob, &
    naqcol, &
    ngascol, &
    nnonaqcol, &
    c, &
    numsp, &
    ithcomp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the icomp aqueous components in the mobile phases (aqueous, gas, and non-aqueous phases)
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)    :: this

real*8, pointer, dimension(:)                :: iumob 

integer, intent(in)                          :: numsp

integer, intent(in)                          :: ithcomp

integer, intent(out)                         :: naqcol

integer, intent(out)                         :: ngascol

integer, intent(out)                         :: nnonaqcol

real*8, intent(in), dimension(numsp)         :: c

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
integer                :: &
 i, &
 isp1, &
 isp2
character(len=100)     :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
!%------------------------------------------------------------
if (numsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------
if (ithcomp<=0.or.ithcomp>this%numaqprisp) then
 msg='Error in component index'
 goto 10
end if
!%------------------------------------------------------------
!%----------------------------------------in the aqueous phase
!%------------------------------------------------------------
call get_iposspsph_ (this,this%aqphindex,isp1,isp2)
naqcol=1
ngascol=this%numgasph
nnonaqcol=0
!%-------------------------------------------------------------
call check_pointer_ (iumob,naqcol+ngascol+nnonaqcol,.true.)
iumob(naqcol)=dot_product(this%ueq(ithcomp,isp1:isp2),c(isp1:isp2))
!%------------------------------------------------------------
!%-------------------------------------------in the gas phases
!%------------------------------------------------------------
do i=1,this%numgasph
 call get_iposspsph_ (this,this%idgasph(i),isp1,isp2)
 iumob(naqcol+i)=dot_product(this%ueq(ithcomp,isp1:isp2),c(isp1:isp2))
end do
!%------------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'Chemical System:'
print *,'Chemical System:',this%name
print *,'Service: compute_iumob_'
print *,msg
print *,'*******************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine specia_from_solution_type_pchemsys &
   (this, &
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
	cputime, & 
    dc, &
    dg, &
    sktrk, &
    dsktrk, &
    nchemiter)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Specia according type of solution
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)        :: this

integer, intent(in)                              :: numcomp

integer, intent(in)                              :: numsites

integer, intent(in)                              :: numsp            ! Number of species 

integer, intent(in)                              :: numsurf

integer, intent(out)                             :: nummin

integer, intent(in), dimension(numcomp)          :: icon

real*8, intent(out), dimension(numsp)            :: c                ! Concentration vector 

real*8, intent(out), dimension(numsp)            :: g

real*8, intent(in), dimension(numsp)             :: alpha

real*8, intent(in), dimension(numcomp)           :: cguess

real*8, intent(in), dimension(numcomp)           :: ctot

real*8, intent(in), dimension(numsites,numsurf)  :: txoh

real*8, intent(in), dimension(numsites,numsurf)  :: c1

real*8, intent(in), dimension(numsites,numsurf)  :: c2

real*8, intent(in), dimension(numsites,numsurf)  :: spsurfarea

real*8, pointer , dimension(:)                   :: simin

character(len=*), intent(in), dimension(numcomp) :: constraint

character(len=*), intent(in), dimension(numcomp) :: namecomp

character(len=100), pointer, dimension(:)        :: namemin 

real*8, intent(out)                              :: ionstr

logical, intent(out)                             :: iserror

logical, intent(out)                             :: isconvergence

real*8, pointer, optional, dimension(:,:)        :: dc

real*8, pointer, optional, dimension(:,:)        :: dg 

real*8, pointer, optional, dimension(:)          :: sktrk

real*8, pointer, optional, dimension(:,:)        :: dsktrk

integer, intent(out), optional                   :: nchemiter 

real*8, intent(out), optional                    :: cputime         ! CPU time consumed during speciation calculations 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                     :: &
 i, &
 j, &
 isps, &
 naqt, &
 iter, &
 ipos1, &
 ipos2, &
 npri, &
 idiverg, &
 icw, &
 ich, &
 ispsw, &
 info 
real*8                      :: &
 chgbal, &
 dd, &
 errmax, &
 resmax, &
 cputime1, &
 cputime2 
real*8, pointer             :: &
 siloc(:) => null (), &
 u(:) => null (), &
 du(:,:) => null (), &
 deltacunk(:) => null ()
real*8, pointer             :: &
 dcloc(:,:) => null (), &
 dgloc(:,:) => null (), &
 dsiloc(:,:) => null (), &
 dchgbal(:) => null () , &
 ctotloc(:) => null (), &
 cguessloc(:) => null ()
character(len=100), pointer :: &
 constraintloc(:) => null ()
integer, pointer            :: &
 iconloc(:) => null (), &
 idmin(:) => null ()
logical                     :: &
 isconv, &
 isupmxiter, &
 isdivergence, &
 isanomalous, &
 havedc, &
 havedg, &
 havenchemiter, &
 havecputime, &
 isderivatives, &
 isupmxitergam, &
 isbe, &
 isbepri 
real*8, pointer             :: &
 jacobian(:,:) => null (), &
 residual(:) => null (), &
 dionstr(:) => null ()
integer                     :: &
 indx(this%numaqprisp)
character(len=100)          :: &
 namepri, &
 nameoutput, &
 msg 
real*8, parameter           :: &
 zero=1.0d-20
logical, parameter          :: &
 isitergam=.false.     
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
!% Local parameters 
!%-------------------------------------------------------------
integer, parameter :: &
    constr_primary=0, &
    constr_uaq=1, &
    constr_chargebalance=2, &
    constr_activity=3, &
    constr_mineral=4, &
    constr_gas=5
real*8, parameter :: &
    r0=0.0d0, &
    r1=1.0d0
!%-------------------------------------------------------------
msg=''
iserror=.false.
!%-------------------------------------------------------------
!% Check the number of species
!%-------------------------------------------------------------
if (numsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%-------------------------------------------------------------
!% Check the number of surfaces 
!%-------------------------------------------------------------
if (numsurf/=this%numsurf) then
 msg='Error in number of surfaces'
 goto 10
end if
!%-------------------------------------------------------------
!% Check the number of primary aqueus species 
!%-------------------------------------------------------------
if (numcomp>this%numaqprisp) then
 msg='Error in number of components'
 goto 10
end if
!%-------------------------------------------------------------
!% Check optional arguments
!%-------------------------------------------------------------
havedc=present(dc)
havedg=present(dg)
havenchemiter=present(nchemiter)
havecputime=present(cputime)
!%-------------------------------------------------------------
!% Initialice variables
!%-------------------------------------------------------------
isderivatives=.true.
isconvergence=.true. 
!%-------------------------------------------------------------
!% Allocate and initialice local pointers 
!%-------------------------------------------------------------
call check_pointer_ (jacobian,this%numaqprisp,this%numaqprisp,.true.)
call check_pointer_ (residual,this%numaqprisp,.true.)
call check_pointer_ (siloc,this%numsp,.true.)
call check_pointer_ (dcloc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dgloc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dsiloc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (deltacunk,this%numaqprisp,.true.)
call check_pointer_ (u,this%numaqprisp,.true.)
call check_pointer_ (du,this%numaqprisp,this%numaqprisp,.true.)
call check_pointer_ (ctotloc,this%numaqprisp,.true.)
call check_pointer_ (cguessloc,this%numaqprisp,.true.)
call check_pointer_ (constraintloc,this%numaqprisp,.true.)
call check_pointer_ (iconloc,this%numaqprisp,.true.)
!%-------------------------------------------------------------
!% Get water components index and species index
!%-------------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wcindex=icw,wspindex=ispsw, &
                     hcindex=ich)
!%-------------------------------------------------------------
!% Initialice local variables 
!%-------------------------------------------------------------
ctotloc=zero
cguessloc=zero
constraintloc=' '
iconloc=constr_uaq
!%-------------------------------------------------------------
!% Ordering vectors according components
!%-------------------------------------------------------------
isbe=.false. 
do j=1,numcomp
 isbepri=.false.
 do i=1,this%numaqprisp
   namepri=this%pspecies(this%idaqprisp(i))%ptr%name
   if (namepri==namecomp(j)) then
   if (namepri=='h+'.or.namepri=='H+') isbe=.true. 
   isbepri=.true.
   if (ctot(j)>r0) then
    ctotloc(i)=ctot(j)
   end if
   if (cguess(j)>r0) then
    cguessloc(i)=cguess(j)
   end if
   constraintloc(i)=constraint(j)
   iconloc(i)=icon(j)
  end if
 end do
 if (.not.isbepri) then
  msg='Error, this species was not defined as primary:'
  call add_ (msg,namecomp(j))
  iserror=.true.
  goto 20  
 end if
end do
!%-----------------------------------------------------------
!% If the proton component is present in the chemical system
!% but is not defined for the user, then assign pH=7 
!% and the h+ component is constrained to activity 
!%-----------------------------------------------------------
if (.not.isbe.and.ich>0) then
    cguessloc(ich)=1.0d-7 
    iconloc(ich)=constr_activity
    ctotloc(ich)=1.0d-7 
end if
!%-----------------------------------------------------------
!% In the case of water component 
!%-----------------------------------------------------------
if (icw>0) then
    cguessloc(icw)=r1/kgwmol
    iconloc(icw)=constr_primary
end if
!%------------------------------------------------------------------
!% First and last position for aqueous species 
!%------------------------------------------------------------------
call get_iposspsph_ (this,this%aqphindex,ipos1,ipos2)
!%------------------------------------------------------------------
!% Copy the initial primary concentrations in concentration vector
!%------------------------------------------------------------------
c(this%idaqprisp)=cguessloc
!%------------------------------------------------------------------
!% Fill derivatives of concentration of primary species with 
!% respect to primary species 
!%------------------------------------------------------------------
do i=1,this%numaqprisp
 dcloc(this%idaqprisp(i),i)=r1
end do
!%------------------------------------------------------------------
!% Initialice variables 
!%------------------------------------------------------------------
errmax=r0
idiverg=0
iter=0 
isconv=.false. 
!%-------------------------------------------------------------
!% Start the iterative process 
!% Newton Raphson Method 
!%-------------------------------------------------------------
call cpu_time (cputime1)
do
iter=iter+1
jacobian = r0
residual = r0 
!%-------------------------------------------------------------
!% Compute aqueous complexes 
!%-------------------------------------------------------------
call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  isitergam, &       
      isanomalous, &
      isupmxitergam, &
      isderivatives, &
      msg, &
      iserror) 
    if (iserror.or.isupmxitergam.or.isanomalous) goto 20
!%-------------------------------------------------------------
!% Compute partial gas pressures of gases 
!%-------------------------------------------------------------
siloc=c
dsiloc=dcloc
do i=1,this%numgasph 
     call compute_secondaries_ &
                  (this, &
                   siloc, &
                   g, &
                   dsiloc, &
                   dgloc, &
                   this%aqphindex, &
                   this%idgasph(i), &
                   ionstr, &
                   dionstr, &
				   1.0d0, &
				   isitergam, &
                   isanomalous, &
                   isupmxitergam, &
                   isderivatives, &
                   msg, &
                   iserror) 
    if (iserror.or.isupmxitergam.or.isanomalous) goto 20
end do
!%-------------------------------------------------------------
!% Compute saturation of mineral phases 
!%-------------------------------------------------------------
do i=1,this%numminph 
    call compute_secondaries_ &
         (this, &
          siloc, &
          g, &
          dsiloc, &
          dgloc, &
          this%aqphindex, &
          this%idminph(i), &
          ionstr, &
          dionstr, &
		  1.0d0, &
		  isitergam, &
          isanomalous, &
          isupmxitergam, &
          isderivatives, &
          msg, &
          iserror) 
    if (iserror.or.isupmxitergam.or.isanomalous) goto 20 
end do
!%-----------------------------------------------
!% Build the jacobian and residual 
!%-----------------------------------------------
do i=1,this%numaqprisp
 
 select case(iconloc(i))
 
 case(constr_primary)
         residual(i)  = r0
         jacobian(i,i)= r1
 case(constr_uaq) 
         u=matmul(this%ueq(:,ipos1:ipos2),c(ipos1:ipos2))
         du=matmul(this%ueq(:,ipos1:ipos2),dcloc(ipos1:ipos2,:))
         residual (i)   = ctotloc(i) - u(i)
         jacobian (i,:) = du(i,:) 
 case(constr_chargebalance)
         call compute_charge_balance_(this%pphase(this%aqphindex)%ptr, &
		                              chgbal,c(ipos1:ipos2),iserror)
         call compute_dcharge_balance_(this%pphase(this%aqphindex)%ptr, &
		                              dchgbal,dcloc(ipos1:ipos2,:),iserror)
         if (iserror) goto 20
         residual (i)   = - chgbal
         jacobian (i,:) =  dchgbal
 case(constr_activity) 
         isps=this%idaqprisp(i)
         residual (i)   = ctotloc(i)-(c(isps)*g(isps))
         jacobian (i,:) = c(isps)*dgloc(isps,:)
         jacobian (i,i) = jacobian(i,i)+g(isps)
 case(constr_mineral)
         call get_sp_index_(this,constraintloc(i),isps)
         if(isps==0)  then
           iserror=.true. 
           msg='Error, not found in the chemical system the species:'
           call add_ (msg,constraintloc(i))
           goto 20
         end if
         residual (i)   = r1 - siloc(isps)
         jacobian (i,:) = dsiloc(isps,:) 
 case(constr_gas)
         call get_sp_index_(this,constraintloc(i),isps)
         if(isps==0)  then
            iserror=.true. 
            msg='Error, not found in the chemical systems the species:'
            call add_ (msg,constraint(i))
            goto 20
         end if
         residual (i)   = r1 - (siloc(isps)/ctotloc (i))
         jacobian (i,:) = dsiloc(isps,:) / ctotloc(i) 
 end select
 
end do
!%---------------------------------------------------------------------
!% Write information about Newton Raphson 
!%---------------------------------------------------------------------
if (this%iswriteinfo) then
 call write_newton_raphson_info_ &
   (jacobian, &
    residual, &
    this%numaqprisp, &
    this%numaqprisp, &
    iter, &
    iouwiinfo, &
    'specia_ (according solution type)')
end if
!%-------------------------------------------------------------------
!% Solve the lineal system
!%-------------------------------------------------------------------
deltacunk=residual
call ludcmp (jacobian,this%numaqprisp,this%numaqprisp,indx,dd,msg,iserror) 
if (iserror) goto 20 
call lubksb (jacobian,this%numaqprisp,this%numaqprisp,indx,deltacunk)
!%------------------------------------------------------------------
!% Check the convergence in unknowns and residual 
!%------------------------------------------------------------------
   call check_convergence_ &
   (this, &
    isconv, &
    isdivergence, &
    isupmxiter, &
    c(this%idaqprisp), &
    deltacunk, &
    residual, &
    this%numaqprisp, &
    this%numaqprisp, &
    errmax, &
    iter, &
    idiverg, &
    this%tolunknr, &
    this%tolresnr)
 
    c(this%idaqprisp) = c(this%idaqprisp) + deltacunk
 
    if (isupmxiter.or.isdivergence) then
        write(6,*) 'Error, convergence problems'
        goto 20
    end if
 
    if (isconv) exit
 
end do
!%------------------------------------------------------------------
!% 
!%------------------------------------------------------------------
call cpu_time (cputime2)
!%------------------------------------------------------------------
!% 
!%------------------------------------------------------------------
if (.not.havedc.or..not.havedg) isderivatives=.false.
!%------------------------------------------------------------------
!% Update the solution 
!%------------------------------------------------------------------
call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    this%aqphindex, &
    0, &
    ionstr, &
    dionstr, &
	1.0d0, &
	isitergam, &
    isanomalous, &
    isupmxitergam, &
    isderivatives, &
    msg, &
    iserror) 
    if (iserror.or.isupmxitergam.or.isanomalous) goto 20

!%------------------------------------------------------------------
!% Solve adsorption
!%------------------------------------------------------------------
do i=1,this%numsurf
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    txoh(:,i), &
    c1(:,i), &
    c2(:,i), &
    spsurfarea(:,i), &
    numsites, &
    ionstr, &
    i, &
    isderivatives, &
    msg, &
    iserror) 
    if (iserror) goto 20
end do
!%------------------------------------------------------------------
!% Compute kinetic
!%------------------------------------------------------------------
call compute_kinetic_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    alpha, &
    c, &
    0.0d0, &
    msg, &
    iserror, &
    sktrk=sktrk, &
    dsktrk=dsktrk)
    if (iserror) goto 20
!%------------------------------------------------------------
!% Optional 
!%------------------------------------------------------------
if (havedc) then
 call check_pointer_ (dc,this%numsp,this%numaqprisp,.true.)
 dc=dcloc
end if
!%------------------------------------------------------------
!% Optional 
!%------------------------------------------------------------
if (havedg) then
 call check_pointer_ (dg,this%numsp,this%numaqprisp,.true.)
 dg=dgloc
end if
!%------------------------------------------------------------
!% Optional 
!%------------------------------------------------------------
if (havecputime) then
  cputime=cputime2-cputime1
end if
!%------------------------------------------------------------
!% Compute saturation indices for minerals
!%------------------------------------------------------------
call compute_min_sat_index_ &
   (this, &
    simin, &
    idmin, &
    namemin, &
    nummin, &
    c, &
    g, &
    numsp, &
    iserror)
 if (iserror) then
  msg='Error when calling compute_min_sat_index_'
  goto 20
 end if
!%---------------------------------------------------------------
20 continue 
!%---------------------------------------------------------------
!% Deallocate local pointers
!%---------------------------------------------------------------
call check_pointer_ (jacobian,1,1,.false.)
call check_pointer_ (residual,1,.false.)
call check_pointer_ (dcloc,1,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
call check_pointer_ (dsiloc,1,1,.false.)
call check_pointer_ (siloc,1,.false.)
call check_pointer_ (deltacunk,1,.false.)
call check_pointer_ (dchgbal,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (u,1,.false.)
call check_pointer_ (du,1,1,.false.)
call check_pointer_ (idmin,1,.false.)
call check_pointer_ (ctotloc,1,.false.)
call check_pointer_ (cguessloc,1,.false.)
call check_pointer_ (constraintloc,1,.false.)
call check_pointer_ (iconloc,1,.false.)
!%--------------------------------------------------------------
if (havenchemiter) nchemiter=iter
!%--------------------------------------------------------------
if (.not.isconv) isconvergence=.false.
!%--------------------------------------------------------------
if (iserror) goto 10
!%--------------------------------------------------------------
!%--------------------------------------------------------------
return
 
10 continue 
print *,'******************'
print *,'Chemical System:'
print *,'Name:', this%name
print *,'Service: specia_'
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
subroutine specia_from_u_pchemsys &
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
    capint, &
    capext, &
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
!   $Description: Compute the chemical speciation from total of components 
! in all phases
!
!   Jabobian matrix presents the following structure 
! 
! 
!   J11    !   J12    !   J13   R1
!   -------------------------   --
!   J21    !   J22    !   J23   R2
!   -------------------------   --
!   J31    !   J32    !   J33   R3 
!
! where the subindice 3 corresponding to total mol of gas species 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in), target :: this                 ! Type parent chemical system variable 

logical, intent(out)                              :: iserror              ! iserror=true, then there was an error. 

character(len=*), intent(out)                     :: msg                  ! Massage error

logical, intent(out)                              :: isconvergence        ! If .true. there wasn`t convergence
   
real*8, intent(out)                               :: ionstr               ! Ionic strength 
 
real*8, intent(out)                               :: factoromgw           ! Mass of free water 

integer, intent(in)                               :: nsp                  ! Number of species 

integer, intent(in)                               :: npri                 ! Number of aqueous primary species. 

integer, intent(in)                               :: ntxoh                ! Total number of sites for sorption or cation exchange.

integer, intent(in)                               :: nsurf                ! Total number of surfaces 

real*8, intent(inout), dimension(nsp)             :: c                    ! Concentration vector. 

real*8, intent(inout), dimension(nsp)             :: g                    ! Activity coefficient vector. 

real*8, intent(in), dimension(nsp)                :: alpha                ! Reactive area of minerals in k+1 

real*8, intent (in), dimension(npri)              :: t                    ! Concentrations of total components 
 
real*8, intent (in), dimension(ntxoh,nsurf)       :: txoh                 ! Total of sites. 

real*8, intent (in), dimension(ntxoh,nsurf)       :: capint               ! Internal capacitance 

real*8, intent (in), dimension(ntxoh,nsurf)       :: capext               ! External capacitance 

real*8, intent (in), dimension(ntxoh,nsurf)       :: spsurfarea           ! 

real*8, intent (in)                               :: dtime                ! Time increment 

logical, intent(in)                               :: iseqmin              ! If true the equilibrate with mineral phases

logical, intent(in)                               :: iseqgas              ! If true the equilibrate with gases phases 

real*8, intent(in)                                :: volgas               ! Volume of gas

real*8, intent(in)                                :: pgas                 ! Gas preassure

real*8, intent(in)                                :: temp                 ! Temperature

real*8, pointer, optional, dimension(:,:)         :: dc                   ! Derivatives of the concentrations. 

real*8, pointer, optional, dimension(:)           :: sktrk                ! Kinetics changes in species due kinetic reactions [nsp]

real*8, pointer, optional, dimension(:,:)         :: dsktrk               ! Derivatives of the Sktrk

real*8, pointer, optional, dimension(:,:)         :: dg                   ! Derivatives of the activity coefficients

real*8, pointer, optional, dimension(:)           :: simin                ! Saturation indices for minerals 

integer, intent(out), optional                    :: nchemiter            ! Number of iterations

real*8, intent (in), optional, dimension(npri)    :: cguess               ! Initial concentrations for primary species 

real*8, intent(in), optional                      :: thetakin             ! Temporal weight for kinetic term 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License: Barcelona, 2007 
!
!-------------------------------------------------------------------------
integer                      :: &
 nmin, &
 nminglob, &
 i, &
 j, &
 k, &
 iter, &
 nunk, &
 itry, &
 idiverg, &
 ispsw, &
 ipos1, &
 ipos2, &
 isp, &
 info 
integer, pointer            :: &
 idminsp(:)  => null (), &
 idminglob(:) => null (), &
 indx(:)  => null ()
real*8                      :: &
 dd, &
 errmax, &
 icpu1, &
 icpu2, &
 dtcpuaqspc, &
 sum1 
real*8, pointer             :: &
 cunk(:)  => null (), &
 si(:)  => null (), &
 sktrk1(:)  => null (), &
 dsktrk1(:,:)  => null (), &
 cold(:)  => null (), &
 gold(:)  => null (), &
 tc(:)  => null (), &
 residual(:) => null (), &
 dcloc(:,:)  => null (), &
 dsi(:,:)  => null (), &
 dgloc(:,:)  => null (), &
 dc2(:,:)  => null (), &
 du2(:,:)  => null (), &
 dtc(:,:)  => null (), &
 jacobian(:,:)  => null (), &
 dionstr(:)  => null (), &
 delta(:)  => null (), &
 u(:,:) => null (), &
 cguess1(:) => null ()
real*8, pointer         :: &
 nh2o => null (), &
 ngas => null ()
logical                     :: &
 isbenegmin, &
 isbenewmin, &
 isneggas, &
 isdivergence, &
 isconv, &
 isupmxiter, &
 isanomalous, &
 havecguess, &
 havedc, &
 havedg, &
 havedsktrk, &
 havesimin, &
 havenchemiter, &
 havethetakin, &
 isderivate, &
 isupmxitergam, &
 iscompgas
character(len=100)          :: &
 nameoutput
character(len=100), pointer :: &
 name => null ()
real*8, parameter           :: &
 r0=0.0d0, &
 r1=1.0d0, &
 r5=5.0d0  
!-------------------------------------------------------------------------
!
!   $code
!
msg=' '
iserror=.false.
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
isconvergence=.true.
isderivate=.true.
factoromgw=r1 
isneggas=.false. 
if (this%numgasph>0) then
 iscompgas=.true.
else 
 iscompgas=.false.
end if 
!%-----------------------------------------------------------------
!%Check optional arguments
!%-----------------------------------------------------------------
havedc=present(dc)
havedg=present(dg)
havedsktrk=present(dsktrk)
havesimin=present(simin)
havenchemiter=present(nchemiter)
havecguess=present(cguess)
havethetakin=present(thetakin)
!%------------------------------------------------------------------
!% Check number of species
!%------------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 iserror=.true.
 goto 10
end if
!%------------------------------------------------------------------
!% Check number of surfaces
!%------------------------------------------------------------------
if (nsurf/=this%numsurf) then
 msg='Error in number of surfaces'
 iserror=.true. 
 goto 10
end if
!%-----------------------------------------------------------------
!% Check number of primary species
!%-----------------------------------------------------------------
if(npri/=this%numaqprisp) then
  msg='Error in number of primary species'
  iserror=.true.
  goto 10
end if
!%-----------------------------------------------------------------
!%Find the component and species index of water specie 
!%-----------------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wspindex=ispsw)
!%-----------------------------------------------------------------
!%Allocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (cold,this%numsp,.true.)
call check_pointer_ (si,this%numsp,.true.)
call check_pointer_ (gold,this%numsp,.true.)
call check_pointer_ (dcloc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dsi,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dgloc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (tc,this%numaqprisp,.true.)
call check_pointer_ (dtc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (cguess1,this%numaqprisp,.true.)
!%-----------------------------------------------------------------
!%Copy c and g vectors 
!%-----------------------------------------------------------------
cold=c
gold=g
!%--------------------------------------------------------------------
!%Return the indices of minerals defined in the chemical system
!%If dtime=0 then use the general components matrix (U) and considers 
!%all minerals (equilibrium and kinetic).
!%If the dtime>0 then use the equilibrium components matrix (Ueq) 
!%(only minerals in equilibrium)
!%--------------------------------------------------------------------
if (dtime>r0) then
 call get_chem_info_ &
   (this, &
    msg, &
    iserror, &
    ideqminsp=idminglob, &
    neqminsp=nminglob) 
 u => this%ueq
else
 call get_chem_info_ &
   (this, &
    msg, &
    iserror, &
    idminsp=idminglob, &
    nminsp=nminglob)
 u => this%u
end if
 if (iserror) goto 20
!%--------------------------------------------------------------------
! Determine the first set of minerals 
!%--------------------------------------------------------------------
call check_mineral_ &
   (this, &
    c, &
    idminglob, &
    nminglob, &
    idminsp, &
    nmin, &
	isbenegmin, &
    isbenewmin)
!%----------------------------------------------------------------
!%Determine the initial values of concentrations of primary 
!%species  
!%----------------------------------------------------------------
if (havecguess) then
 cguess1=cguess
 if (this%wcindex/=0) cguess1(this%wcindex)=r1/kgwmol
else
 cguess1=t/r5 
 if (this%wcindex/=0) cguess1(this%wcindex)=r1/kgwmol
 if (this%ecindex/=0) cguess1(this%ecindex)=1.0d-7 
 if (this%hcindex/=0) cguess1(this%hcindex)=1.0d-7 
end if 
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
c(this%idaqprisp)=cguess1
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%Fill dc
!%----------------------------------------------------------------
do i=1,this%numaqprisp
 dcloc(this%idaqprisp(i),i)=r1
end do
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
! Start loop for negative and try
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
itry=0
 
do
 
 errmax=r0
 idiverg=0
 g=gold
 itry=itry+1
!%--------------
 if (itry>mxtry) then
  msg='Error, number of trys to change mineral sets exceeded'
  isconvergence=.false.
  goto 20
 end if
!%----------------------------------------------------------------
! Compute total number of unknowns of the system equations 
!%----------------------------------------------------------------
 nunk=this%numaqprisp 
 if (iseqmin) nunk=nunk+nmin
 if (iscompgas.and.iseqgas) nunk=nunk+1 
!%----------------------------------------------------------------
!%Allocate local pointers 
!%----------------------------------------------------------------
 call check_pointer_ (jacobian,nunk,nunk,.true.)
 call check_pointer_ (residual,nunk,.true.)
 call check_pointer_ (dtc,this%numaqprisp,this%numaqprisp,.true.)
 call check_pointer_ (indx,nunk,.true.)
 call check_pointer_ (delta,nunk,.true.)
 call check_pointer_ (cunk,nunk,.true.)
!%----------------------------------------------------------------
!% Put the initial solution 
!%----------------------------------------------------------------
 cunk(1:this%numaqprisp) = cguess1
 if (this%wcindex>0) then
  nh2o => cunk(this%wcindex)
 else
  allocate (nh2o)
  nh2o = r1/kgwmol
 end if
!%----------------------------------------------------------------
! ngas pointer to 
!%----------------------------------------------------------------
 if (iscompgas.and.iseqgas) then
  ngas => cunk(nunk)
 end if 
!%----------------------------------------------------------------
!% Initialice the unknowns vector 
!%----------------------------------------------------------------
 if (iseqmin.and.nmin>0) then
  cunk(this%numaqprisp+1:this%numaqprisp+nmin)=c(idminsp)
 end if
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
! Loop for Newton Raphson
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
iter = 0

do
 
  
  iter=iter+1
!%----------------------------------------------------------------
!% Zeroing arrays
!%----------------------------------------------------------------
  jacobian=r0
  residual=r0
  tc=r0
  dtc=r0
!%----------------------------------------------------------------
!% Compute factoromgw
!%----------------------------------------------------------------
  factoromgw=nh2o*kgwmol
!%----------------------------------------------------------------
! Specia aqueous phase 
!%----------------------------------------------------------------
 if (this%iswriteinfo) call cpu_time (icpu1)
 
 call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .true., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
 
   if (this%iswriteinfo) then
     call cpu_time (icpu2)
     dtcpuaqspc=icpu2-icpu1
   end if
 
   if (iserror.or.isupmxitergam.or.isanomalous) exit
!%--------------------------------------------------------------
! Compute concentration for gases (only if number of gas phases>0 
!%-------------------------------------------------------------- 
if (iscompgas) then 
 if (iseqgas) then
  call check_pointer_ (dc2,this%numsp,1,.true.)
  call check_pointer_ (du2,this%numaqprisp,1,.true.)
 end if 
!%--------------------------------------------------------------
! Compute gas preassure for gas phases 
!%-------------------------------------------------------------- 
 do i=1,this%numgasph
   call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      this%idgasph(i), &
      ionstr, &
      dionstr, &
	  r1, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
   if (iserror.or.isupmxitergam.or.isanomalous) exit
   call get_iposspsph_(this,this%idgasph(i),ipos1,ipos2)
   if (iseqgas) then 
     residual(nunk)=residual(nunk)+sum(c(ipos1:ipos2))
	 do j=ipos1,ipos2 
	  jacobian(nunk,1:this%numaqprisp)= &
	  jacobian(nunk,1:this%numaqprisp) + &
	  dcloc(j,1:this%numaqprisp)
     end do
	 dc2(ipos1:ipos2,1)=c(ipos1:ipos2)/(factoromgw*pgas) 
	 c(ipos1:ipos2)=ngas*c(ipos1:ipos2)/(factoromgw*pgas)
     dcloc(ipos1:ipos2,1:this%numaqprisp)= &
	 ngas*dcloc(ipos1:ipos2,1:this%numaqprisp)/(factoromgw*pgas)
  else 
     c(ipos1:ipos2)=c(ipos1:ipos2)*volgas/(rgas*(temp+273.15d0))
     dcloc(ipos1:ipos2,1:this%numaqprisp)= &
	 dcloc(ipos1:ipos2,1:this%numaqprisp)*volgas/(rgas*(temp+273.15d0))
  end if 
 end do
 if (iseqgas) then
  du2=matmul(u,dc2)
  residual(nunk)=pgas-residual(nunk)
  jacobian(1:this%numaqprisp,nunk:nunk)=du2 
  call check_pointer_ (dc2,1,1,.false.)
  call check_pointer_ (du2,1,1,.false.)
 end if
end if 
!%--------------------------------------------------------------
! Compute saturation indices for minerals 
!%--------------------------------------------------------------
si=c
dsi=dcloc
do i=1,this%numminph
 
    call compute_secondaries_ &
     (this, &
      si, &
      g, &
      dsi, &
      dgloc, &
      this%aqphindex, &
      this%idminph(i), &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam.or.isanomalous) exit
 
end do
!%---------------------------------------------------------------
! Specia surface complexes
!%---------------------------------------------------------------
do i=1,this%numsurf
 
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    txoh(:,i)/factoromgw, &
    capint(:,i), &
    capext(:,i), &
    spsurfarea(:,i), &
    ntxoh, &
    ionstr, &
    i, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror) exit
 
end do
!%------------------------------------------------------------------
! Compute kinetic if have theta kin 
!%------------------------------------------------------------------
if (havethetakin) then  
 call compute_kinetic_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    alpha/factoromgw, &
    cold/factoromgw, &
    dtime, &
    msg, &
    iserror, &
    sktrk=sktrk1, &
    dsktrk=dsktrk1)
  if (iserror) goto 20
  tc = - thetakin * matmul(u,sktrk1) * dtime 
  dtc = - thetakin * matmul(u,dsktrk1) * dtime 
end if 
!%---------------------------------------------------------------
! Assembles the part corresponding to the equilibrium minerals
!%---------------------------------------------------------------
if (iseqmin.and.nmin>0) then
!%---------------------
!% Build J21
!%---------------------
     call check_pointer_ (dc2,this%numsp,nmin,.true.)
     call check_pointer_ (du2,this%numaqprisp,nmin,.true.)
     do j=1,nmin
        isp=idminsp(j)
        dc2(isp,j)=r1
        residual(this%numaqprisp+j)=r1-si(isp)
        jacobian(this%numaqprisp+j,1:this%numaqprisp)=dsi(isp,:)
     end do
     du2=matmul(u,dc2)
!%---------------------
!% Build J12
!%---------------------
     jacobian(1:this%numaqprisp,this%numaqprisp+1:this%numaqprisp+nmin)=factoromgw*du2
     call check_pointer_ (dc2,1,1,.false.)
     call check_pointer_ (du2,1,1,.false.)
 
end if
!%-------------------------------------------------------------------
! Build J11 y R1
!%-------------------------------------------------------------------
  tc = tc + matmul(u,c)
  dtc = dtc + matmul(u,dcloc)
  residual(1:this%numaqprisp) = t - factoromgw*tc
  jacobian(1:this%numaqprisp,1:this%numaqprisp) = factoromgw*dtc
!%-----------------------------------------------------------------
! If the water component mass balance is considered, the compute
! derivatives 
!%-----------------------------------------------------------------
  if (this%wcindex>0.and.this%iswcompbal) then
   jacobian(1:this%numaqprisp,this%wcindex)=kgwmol*tc
  else if (this%wcindex>0.and..not.this%iswcompbal) then
   jacobian(this%wcindex,:)=r0
   jacobian(:,this%wcindex)=r0
   jacobian(this%wcindex,this%wcindex)=r1
   residual(this%wcindex)=r0 
  end if
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!% Optional: Write Newton Raphson Information
!%--------------------------------------------------------------------
if (this%iswriteinfo) then
 
 call write_newton_raphson_info_(jacobian,residual,nunk,nunk,iter,iouwiinfo,'specia_')
 
 if (iseqmin.and.nmin>0) then
  write (iouwiinfo,*) '------------------'
  write (iouwiinfo,*) 'Saturation Indices' 
  write (iouwiinfo,*) '------------------'
  do i=1,nmin
    name => this%pspecies(idminsp(i))%ptr%name
    write (iouwiinfo,3),name,'=',si(idminsp(i))
  end do
  write (iouwiinfo,4),'Time spent in aqueous speciation:',dtcpuaqspc,'[s]'
 end if
 
2 format(<nunk>e10.3,a2,e10.3)
3 format(a20,a2,e10.3)
4 format(a33,e10.3,a8)
5 format(a20,2x,a10)
 
end if
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%----------------------------------------------------------------------
delta=residual
if (any(isnan(delta))) then
 msg='Error, NAN values in residual vector during Newton-Raphson process, iter='
 call add_ (msg,iter)
 iserror=.true.
 goto 20
end if
if (any(isnan(jacobian))) then
 msg='Error, NAN values in jacobian matrix during Newton-Raphson process, iter='
 call add_ (msg,iter)
 iserror=.true.
 goto 20
end if
!%----------------------------------------------------------------------
! Solve the lineal system
!%----------------------------------------------------------------------
call ludcmp (jacobian,nunk,nunk,indx,dd,msg,iserror)
if (iserror) goto 20 
call lubksb (jacobian,nunk,nunk,indx,delta)
!%----------------------------------------------------------------------
!% Write information about solution 
!%----------------------------------------------------------------------
if (this%iswriteinfo) then
 
  write (iouwiinfo,*) '-----------------------'
  write (iouwiinfo,*) 'Newton Raphson solution' 
  write (iouwiinfo,*) '-----------------------'
  write (iouwiinfo,*) '     cunk     delta    '
  write (iouwiinfo,*) '-----------------------' 
  do i=1,this%numaqprisp
    name => this%pspecies(this%idaqprisp(i))%ptr%name
    write (iouwiinfo,6),name,'=',cunk(i),delta(i)
  end do
  if(iseqmin) then
   do i=1,nmin
    name => this%pspecies(idminsp(i))%ptr%name
    write (iouwiinfo,6),name,'=',cunk(this%numaqprisp+i),delta(this%numaqprisp+i)
   end do
  end if 
6 format(a20,a2,2e10.3)
 
end if
!%----------------------------------------------------------------------
! Check convergence 
!%----------------------------------------------------------------------
call check_convergence_ &
   (this, &
    isconv, &
    isdivergence, &
    isupmxiter, &
    cunk(1:this%numaqprisp), &
    delta(1:this%numaqprisp), &
    residual, &
    nunk, &
    this%numaqprisp, &
    errmax, &
    iter, &
    idiverg, &
    this%tolunknr, &
    this%tolresnr)
!%----------------------------------------------------------------------
! Update the solution
!%----------------------------------------------------------------------
    cunk = cunk + delta 
    c(this%idaqprisp)=cunk(1:this%numaqprisp)
	if (ispsw>0) c(ispsw)=r1/kgwmol 
    if (iseqmin.and.nmin>0) then
       c(idminsp)=cunk(this%numaqprisp+1:this%numaqprisp+nmin)
    end if
    
    if (isconv.or.isdivergence.or.isupmxiter) then
	  if (isconv.and.this%iswriteinfo) then
	   call write_ (this,c,g,25.0d0,ionstr,this%numsp,iouwiinfo, &
	                nh2o,msg,iserror,simin=si)
	  end if 
	  exit
	end if
 
end do
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
! Finished loop for Newton Raphson
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
if (iserror) goto 20
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
if (isanomalous.or.isupmxitergam) then
 isconvergence=.false.
 goto 20
end if
!%-------------------------------------------------------------------
!% Check if the maximum number of divergences was e
!%-------------------------------------------------------------------
if (isdivergence) then
 msg='Convergence problems, divergence detected'
 isconvergence=.false.
 goto 20
end if
!%------------------------------------------------------------------
!% Check if the maximum number of iterations was exedded
!%------------------------------------------------------------------
if (isupmxiter) then
  msg='Convergence problems, maximum number iterations exceeded'
  isconvergence=.false.
  goto 20
end if
!%----------------------------------------------------------------------
! Check concentrations of equilibrium minerals and check the saturation 
! of another minerals
!%----------------------------------------------------------------------
if (iseqmin) then
 call check_mineral_ &
 (this, &
  c, &
  idminglob, &
  nminglob, &
  idminsp, &
  nmin, &
  isbenegmin, &
  isbenewmin, &
  si)
end if 
!%----------------------------------------------------------------------
! Check if the solution is subsaturated if the equilibrium of gases is 
! imposed 
!%----------------------------------------------------------------------
if (iscompgas) then
 call check_gas_ (this,c,isneggas)
 if (isneggas) iscompgas=.false.
end if 
!%----------------------------------------------------------------------
! if there aren't new saturated minerals and there aren't negative 
! minerals, then the problen are finished
!%---------------------------------------------------------------------- 
if (.not.isbenegmin .and..not.isbenewmin.and..not.isneggas) then
  exit
else
  if (this%wcindex>0) then
   cguess1(this%wcindex)=cunk(this%wcindex)
  end if 
  if (this%iswriteinfo) then
    write (iouwiinfo,*),'-----------------------------'
    write (iouwiinfo,*),'THE MINERAL SET WAS CHANGED!!'
    write (iouwiinfo,*),'New mineral set:'
    write (iouwiinfo,*),'-----------------------------'
    write (iouwiinfo,*),'try=', itry
    write (iouwiinfo,5),'Mineral','Sat.Indices'
    write (iouwiinfo,*),'-----------------------------'
    do i=1,nmin
     name => this%pspecies(idminsp(i))%ptr%name
     write (iouwiinfo,3),name,'SI',si(idminsp(i))
    end do
    write (iouwiinfo,*),'-----------------------------'
  end if
end if
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
end do
!%----------------------------------------------------------------
!%----------------------------------------------------------------
! Finished loop for negative and try
!%----------------------------------------------------------------
!%----------------------------------------------------------------
if (.not.havedc &
        .or. &
    .not.havedg &
        .or. &
   .not.havedsktrk) isderivate=.false.
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
! If there was convergence, then update the speciation 
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .true., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
 
if (iserror.or.isupmxitergam.or.isanomalous) goto 20
!%--------------------------------------------------------------
! Compute concentration for gases
!%--------------------------------------------------------------
if (iscompgas) then 
 do i=1,this%numgasph
   call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      this%idgasph(i), &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
   if (iserror.or.isupmxitergam.or.isanomalous) exit
   call get_iposspsph_(this,this%idgasph(i),ipos1,ipos2)
   if (iseqgas) then
     c(ipos1:ipos2)=ngas*c(ipos1:ipos2)/(factoromgw*pgas)
     dcloc(ipos1:ipos2,1:this%numaqprisp)= &
	 ngas*dcloc(ipos1:ipos2,1:this%numaqprisp)/(factoromgw*pgas)
  else 
     c(ipos1:ipos2)=c(ipos1:ipos2)*volgas/(rgas*(temp+273.15d0))
     dcloc(ipos1:ipos2,1:this%numaqprisp)= &
	 dcloc(ipos1:ipos2,1:this%numaqprisp)*volgas/(rgas*(temp+273.15d0))
  end if 
 end do
end if 
!%------------------------------------------------------------------
! Make the speciation of surface complexes
!%------------------------------------------------------------------
do i=1,this%numsurf
 
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    txoh(:,i)/factoromgw, &
    capint(:,i), &
    capext(:,i), &
    spsurfarea(:,i), &
    ntxoh, &
    ionstr, &
    i, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror) goto 20
 
end do
!%------------------------------------------------------------------
! Compute kinetic
!%------------------------------------------------------------------
call compute_kinetic_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    alpha/factoromgw, &
    cold/factoromgw, &
    dtime, &
    msg, &
    iserror, &
    sktrk=sktrk, &
    dsktrk=dsktrk)
if (iserror) goto 20
!%------------------------------------------------------------------
!%Optional
!%------------------------------------------------------------------
if (havedc) then
 call check_pointer_ (dc,this%numsp,this%numaqprisp,.true.)
 dc=dcloc
end if
!%------------------------------------------------------------------
!%Optional
!%------------------------------------------------------------------
if (havedg) then
 call check_pointer_ (dg,this%numsp,this%numaqprisp,.true.)
 dg=dgloc
end if
!%------------------------------------------------------------------
!% Optional
!%------------------------------------------------------------------
if (havesimin) then
 call check_pointer_ (simin,this%numsp,.true.)
 simin=si
end if
!%------------------------------------------------------------------
20 continue     
!%------------------------------------------------------------------
! Deallocate local pointers 
!%------------------------------------------------------------------
call check_pointer_ (cunk,1,.false.)
call check_pointer_ (cguess1,1,.false.)
call check_pointer_ (delta,1,.false.)
call check_pointer_ (indx,1,.false.)
call check_pointer_ (jacobian,1,1,.false.)
call check_pointer_ (residual,1,.false.)
call check_pointer_ (si,1,.false.)
call check_pointer_ (dsi,1,1,.false.)
call check_pointer_ (tc,1,.false.)
call check_pointer_ (dtc,1,1,.false.)
call check_pointer_ (idminsp,1,.false.)
call check_pointer_ (dcloc,1,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (idminglob,1,.false.)
call check_pointer_ (sktrk1,1,.false.)
call check_pointer_ (dsktrk1,1,1,.false.)
!%------------------------------------------------------------------
!% Nullify local pointers 
!%------------------------------------------------------------------
u => null ()
name => null () 
if (this%wcindex>0) then
 nh2o => null ()  
else
 deallocate (nh2o)
 nh2o => null ()
end if  
!%------------------------------------------------------------------
!% Nullify pointer to ngas unknown 
!%------------------------------------------------------------------
ngas => null ()
!%------------------------------------------------------------------
if (havenchemiter) nchemiter=iter 
!%------------------------------------------------------------------
if (iserror &
      .or. &
    isanomalous &
     .or. &
    isupmxiter &
      .or. &
   isupmxitergam &
    .or. &
   isdivergence) then
 c=cold
 g=gold
end if
!%------------------------------------------------------------------
call check_pointer_ (cold,1,.false.)
call check_pointer_ (gold,1,.false.)
if (iserror &
      .or. &
    isanomalous &
     .or. &
    isupmxiter &
      .or. &
   isupmxitergam &
    .or. &
   isdivergence) then
 goto 10
end if
!%------------------------------------------------------------------
!%------------------------------------------------------------------
!%------------------------------------------------------------------
return
 
10 continue 
return 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine add_and_specia_pchemsys &
   (this, &
    c, &
    g, &
    ionstr, &
    simin, &
    namespadd, &
    moladd, &
    mh2o, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    nsp, &
    nsurf, &
    ntxoh, &
    omgwfreek, &
    isconvergence, &
    iserror, &
    uadd, &
    dc, &
    sktrk, &
    dsktrk, &
    dg, &
    ioutput, &
    nchemiter)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Add mol of components (uad) or mass of water (mh2o) and make the speciation. 
! The equilibrium with other phases is not imposed. 
!
!   $Arguments:
!
 
 
type (t_parentchemicalsystem), intent(in) :: this

integer, intent(in)                                      :: nsp ! Number of species. 

real*8, intent(inout), dimension(nsp)                    :: c ! concentration vector

real*8, intent(inout), dimension(nsp)                    :: g ! activity coefficients vector

real*8, intent(out), dimension(nsp)                      :: simin 

real*8, intent(inout)                                    :: omgwfreek ! Mass of free water 

real*8, intent (out)                                     :: ionstr   ! Ionic strength

real*8, pointer, dimension(:,:), optional                :: dc

real*8, pointer, dimension(:), optional                  :: sktrk

real*8, pointer, dimension(:,:), optional                :: dsktrk

real*8, pointer, dimension(:,:), optional                :: dg 

logical, intent(out)                                     :: iserror

logical, intent(out)                                     :: isconvergence ! If .true. there was convergence. 

real*8, intent (in)                                      :: mh2o ! mass of water 

real*8, intent (in)                                      :: moladd ! Mass of specie added/sustracted (in mol) 

integer, intent(in)                                      :: nsurf

integer, intent(in)                                      :: ntxoh

character(len=*), intent(in)                             :: namespadd  ! Name of species added/sustracted

real*8, intent(in), dimension(ntxoh,nsurf)               :: txoh

real*8, intent(in), dimension(ntxoh,nsurf)               :: c1

real*8, intent(in), dimension(ntxoh,nsurf)               :: c2

real*8, intent(in), dimension(ntxoh,nsurf)               :: spsurfarea

real*8, intent(in), optional, dimension(this%numaqprisp) :: uadd 

integer, intent(in), optional                            :: ioutput 

integer, intent(out), optional                           :: nchemiter
 
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
 i, &
 j, &
 k, &
 iter, &
 nunk, &
 idiverg, &
 ispsw, &
 icw, &
 ipos1, &
 ipos2, &
 npri, &
 isp, &
 nph, &
 nstep, &
 istep, &
 itry
integer, pointer            :: &
 indx(:)
real*8                      :: &
 dd, &
 errmax, &
 dtomg, &
 omg, &
 omgwcrystk, &
 omgwfreek1, &
 omgwcrystk1, &
 mh2o1, &
 moladd1
real*8, pointer             :: &
 cunk(:) => null (), &
 cold(:) => null (), &
 tk1(:) => null (), &
 tk(:) => null (), &
 residual(:) => null (), &
 dcloc(:,:) => null (), &
 dsi(:,:) => null (), &
 dgloc(:,:) => null (), &
 dtk1(:,:) => null (), &
 dtm(:) => null (), &
 jacobian(:,:) => null (), &
 dionstr(:) => null (), &
 dtcunk(:) => null ()
logical                     :: &
 isdivergence, &
 isupmxiter, &
 isupmxitergam, &
 isanomalous, &
 haveioutput, &
 havedc, &
 havedg, &
 havedsktrk, &
 haveuadd, &
 havenchemiter, &
 isderivate, &
 be
character(len=100)         :: &
 msg, &
 nameoutput, &
 name 
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
msg=' '
iserror=.false.
!%-----------------------------------------------------------------
 
!%-----------------------------------------------------------------
isderivate=.true.
isconvergence=.false.
!%-----------------------------------------------------------------
! Check optional arguments
!%-----------------------------------------------------------------
havedc=present(dc)
havedg=present(dg)
havedsktrk=present(dsktrk)
haveioutput=present(ioutput)
haveuadd=present(uadd)
havenchemiter=present(nchemiter)
!%------------------------------------------------------------------
! Check number of species
!%------------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 iserror=.true.
 goto 10
end if
!%------------------------------------------------------------------
! Check number of surfaces
!%------------------------------------------------------------------
if (nsurf/=this%numsurf) then
 msg='Error in number of surfaces'
 goto 10
end if
!%-----------------------------------------------------------------
!% Optional: Open unit for to write
!%-----------------------------------------------------------------
if (this%iswriteinfo.and.haveioutput) then
 nameoutput=''
 call add_ (nameoutput,ioutput)
 nameoutput(6:30)='add_and_specia_info.dat'
 open(unit=ioutput,file=nameoutput(1:30),status='unknown')
end if
!%-----------------------------------------------------------------
!% Give the water component and species index
!%-----------------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wcindex=icw,wspindex=ispsw)
!%-----------------------------------------------------------------
! Check if the species added is defined in the chemical system
!%-----------------------------------------------------------------
call get_sp_index_  (this,namespadd,isp)
if (moladd>0.0d0.and.isp==0) then
 msg='Error, not found the species:'
 call add_ (msg,namespadd)
 goto 10
end if
!%-----------------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------------
nunk=this%numaqprisp
call check_pointer_ (cold,this%numsp,.true.)
call check_pointer_ (dcloc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dsi,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dgloc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (tk,this%numaqprisp,.true.)
call check_pointer_ (tk1,this%numaqprisp,.true.)
call check_pointer_ (dtm,this%numaqprisp,.true.)
call check_pointer_ (dtk1,nunk,nunk,.true.)
call check_pointer_ (jacobian,nunk,nunk,.true.)
call check_pointer_ (residual,nunk,.true.)
call check_pointer_ (indx,nunk,.true.)
call check_pointer_ (dtcunk,nunk,.true.)
call check_pointer_ (cunk,nunk,.true.)
call check_pointer_ (dsi,this%numsp,this%numaqprisp,.true.)
!%----------------------------------------------------------------
tk=matmul(this%ueq,c)
!%----------------------------------------------------------------
!%---------------------------------------------------Fill c and dc
!%----------------------------------------------------------------
cold=c
!%----------------------------------------------------------------
do i=1,this%numaqprisp
 dcloc(this%idaqprisp(i),i)=1.0d0
end do
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
nstep=0
itry=0
 
do
 
itry=itry+1
nstep=nstep+1
moladd1=moladd/nstep
mh2o1=mh2o/nstep
omgwfreek1=omgwfreek
c=cold
 
if (itry>mxtry) then
  msg='Error, convergence problems'
  goto 20 
end if
 
 
do istep=1,nstep
!%----------------------------------------------------------------
! Compute the increment of mass of water 
!%----------------------------------------------------------------
dtomg=0.0d0
!%----------------------------------------------------------------
if (ispsw>0.and.isp>0) then
 dtomg=this%ueq(icw,isp)*moladd1*kgwmol
end if
!%----------------------------------------------------------------
if (ispsw>0) then
 dtomg=dtomg+mh2o1
end if 
!%----------------------------------------------------------------
omgwfreek1=omgwfreek1+dtomg
!%----------------------------------------------------------------
! Compute the variation of mass of components
!%----------------------------------------------------------------
dtm=0.0d0
if (moladd>0.0d0.and.isp>0) then
 dtm=this%ueq(:,isp)*moladd1
end if
if (mh2o>0.0d0.and.haveuadd) then
 dtm=dtm+mh2o1*uadd
end if
!%----------------------------------------------------------------
! Compute the first guess approach
!%----------------------------------------------------------------
cunk=(c(this%idaqprisp)+dtm/omgwfreek1)/2.0d0
cunk=dabs(cunk)
!%----------------------------------------------------------------
! Loop for Newton Raphson
!%----------------------------------------------------------------
iter = 0
idiverg = 0
!%----------------------------------------------------------------
c(this%idaqprisp)=cunk
!%----------------------------------------------------------------
do
 
  iter=iter+1
!%----------------------------------------------------------------
! Initialice variables  
!%----------------------------------------------------------------
  jacobian=0.0d0
  residual=0.0d0
  tk1=0.0d0
  dtk1=0.0d0
!%----------------------------------------------------------------
! Specia aqueous phase and aqueous phase with another phases
!%----------------------------------------------------------------
  call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .true., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam) goto 20
!%---------------------------------------------------------------
! Specia surface complexes
!%---------------------------------------------------------------
 do i=1,this%numsurf
 
  call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    txoh(:,i), &
    c1(:,i), &
    c2(:,i), &
    spsurfarea(:,i), &
    ntxoh, &
    ionstr, &
    i, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror) goto 20
 
 end do
!%-------------------------------------------------------------------
  tk1=matmul(this%ueq,c)
  dtk1=matmul(this%ueq,dcloc)
  residual=tk1*omgwfreek1-tk*omgwfreek-dtm
  jacobian=omgwfreek1*dtk1
!%--------------------------------------------------------------------
! Optional: Write Newton Raphson Information
!%--------------------------------------------------------------------
 if (this%iswriteinfo.and.haveioutput) then
 
  call write_newton_raphson_info_ &
   (jacobian, &
    residual, &
    nunk, &
    nunk, &
    iter, &
    ioutput, &
    'specia_')

 end if
!%----------------------------------------------------------------------
! Solve the lineal system
!%----------------------------------------------------------------------
 residual=-residual
 dtcunk=residual
 call ludcmp (jacobian,nunk,nunk,indx,dd,msg,iserror)
 if (iserror) goto 20 
 call lubksb (jacobian,nunk,nunk,indx,dtcunk)
!%----------------------------------------------------------------------
! Check convergence and update the solution
!%----------------------------------------------------------------------
 call check_convergence_ &
   (this, &
    isconvergence, &
    isdivergence, &
    isupmxiter, &
    cunk(1:this%numaqprisp), &
    dtcunk(1:this%numaqprisp), &
    residual, &
    nunk, &
    this%numaqprisp, &
    errmax, &
    iter, &
    idiverg, &
    this%tolunknr, &
    this%tolresnr)
 
    cunk = cunk + dtcunk
 
    cunk=dabs(cunk)

    c(this%idaqprisp)=cunk(1:this%numaqprisp)

	if (isconvergence.or.isdivergence.or.isupmxiter) exit
 
end do
!%-------------------------------------------------------------------
!% Finished loop for Newton Raphson
!%-------------------------------------------------------------------
 
    if (.not.isconvergence) exit
 
end do
 
    if (isconvergence) exit
 
end do
!%-------------------------------------------------------------------
! Update mass of free water 
!%-------------------------------------------------------------------
omgwfreek=omgwfreek1
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
if (.not.havedc &
        .or. &
    .not.havedg &
        .or. &
   .not.havedsktrk) isderivate=.false.
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
! Update speciation
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .true., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam) goto 20
!%------------------------------------------------------------------
! Compute saturation indices for minerals
!%------------------------------------------------------------------
simin=c
do i=1,this%numminph
 
 call compute_secondaries_ &
     (this, &
      simin, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      this%idminph(i), &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      .false., &
      msg, &
      iserror)
 
    if (iserror) goto 20
 
end do
!%------------------------------------------------------------------
! Make the speciation in gas phases
!%------------------------------------------------------------------
do i=1,this%numgasph
 
  call compute_secondaries_ &
    (this, &
     c, &
     g, &
     dcloc, &
     dgloc, &
     this%aqphindex, &
     this%idgasph(i), &
     ionstr, &
     dionstr, &
	 1.0d0, &
	 .false., &
     isanomalous, &
     isupmxitergam, &
     isderivate, &
     msg, &
     iserror)
 
    if (iserror.or.isupmxitergam) goto 20
 
end do
!%------------------------------------------------------------------
! Make the speciation of surface complexes
!%------------------------------------------------------------------
do i=1,this%numsurf
 
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    txoh(:,i), &
    c1(:,i), &
    c2(:,i), &
    spsurfarea(:,i), &
    ntxoh, &
    ionstr, &
    i, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror) goto 20
 
end do
!%------------------------------------------------------------------
!%----------------------------------------------------------Optional
!%------------------------------------------------------------------
if (havedc) then
 call check_pointer_ (dc,this%numsp,this%numaqprisp,.true.)
 dc=dcloc
end if
!%------------------------------------------------------------------
!% Optional
!%------------------------------------------------------------------
if (havedg) then
 call check_pointer_ (dg,this%numsp,this%numaqprisp,.true.)
 dg=dgloc
end if
!%------------------------------------------------------------------
!% Deallocate vectors and arrays
!%------------------------------------------------------------------
20 continue 
call check_pointer_ (dsi,1,1,.false.)
call check_pointer_ (tk1,1,.false.)
call check_pointer_ (tk,1,.false.)
call check_pointer_ (dtk1,1,1,.false.)
call check_pointer_ (cold,1,.false.)
call check_pointer_ (dcloc,1,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (cunk,1,.false.)
call check_pointer_ (dtcunk,1,.false.)
call check_pointer_ (indx,1,.false.)
call check_pointer_ (jacobian,1,1,.false.)
call check_pointer_ (residual,1,.false.)
call check_pointer_ (dtm,1,.false.)
!%-----------------------------------------------------------------
if (havenchemiter) then
 nchemiter=iter
end if 
!%-----------------------------------------------------------------
if (iserror) goto 10 
!%-----------------------------------------------------------------
!%------------------------------------------------------------------
!%------------------------------------------------------------------
return
 
2 format(<nunk>e10.3,a2,e10.3)
3 format(a20,a2,e10.3) 
10 continue 
print *,'************************************'
print *,'Chemical System:'
print *,'Name:', this%name
print *,'Service: specia_'
print *, msg
print *,'Maximun error in unknown=', errmax
print *,'Iterations=', iter
print *,'Tolerance in unknown=', this%tolunknr
print *,'Tolerance in residual=', this%tolresnr
print *,'************************************'
iserror=.true.
c=cold
call check_pointer_ (cold,1,.false.)
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine equilibrate_pchemsys &
   (this, &
    c, &
    g, &
    si, &
    ionstr, &
    nameph, &
    siph, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    nsp, &
    nsurf, &
    ntxoh, &
    omgwfreek, &
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
!   $Description: Equilibrate solution with mineral phase. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)  :: this

real*8, intent (inout)                     :: omgwfreek         ! mass of water per unit of volume in k

real*8, intent (out)                       :: ionstr

real*8, intent (in)                        :: siph         ! Saturation indice to saturation

integer, intent(in)                        :: nsp

integer, intent(in)                        :: nsurf

integer, intent(in)                        :: ntxoh

character(len=*), intent(in)               :: nameph         ! Name of species added

real*8, intent(inout), dimension(nsp)      :: c ! concentration vector

real*8, intent(inout), dimension(nsp)      :: g ! activity coefficients vector

real*8, intent(out), dimension(nsp)        :: si 

real*8, intent(in), dimension(ntxoh,nsurf) :: txoh

real*8, intent(in), dimension(ntxoh,nsurf) :: c1

real*8, intent(in), dimension(ntxoh,nsurf) :: c2

real*8, intent(in), dimension(ntxoh,nsurf) :: spsurfarea

logical, intent(out)                       :: iserror

logical, intent(out)                       :: isconvergence

real*8, pointer, optional, dimension(:,:)  :: dc

real*8, pointer, optional, dimension(:)    :: sktrk

real*8, pointer, optional, dimension(:,:)  :: dsktrk

real*8, pointer, optional, dimension(:,:)  :: dg

integer, intent(out), optional             :: nchemiter 
 
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
 i, &
 j, &
 k, &
 iter, &
 nunk, &
 idiverg, &
 ispsw, &
 icw, &
 ipos1, &
 ipos2, &
 npri, &
 imin, &
 itry, &
 istep, &
 nstep
integer, pointer            :: &
 indx(:) => null (), &
 idsp(:) => null ()
real*8                      :: &
 dd, &
 errmax, &
 omg, &
 omgwcrystk, &
 omgwfreek1, &
 dtomg, &
 siphloc, &
 a, &
 a2, &
 domgwfreek1, &
 amountph
real*8, pointer             :: &
 cunk(:) => null (), &
 cold(:) => null (), &
 tk1(:) => null (), &
 tk(:) => null (), &
 residual(:) => null (), &
 dcloc(:,:) => null (), &
 dsi(:,:) => null (), &
 dgloc(:,:) => null (), &
 dtk1(:,:) => null (), &
 dtm(:) => null (), &
 jacobian(:,:) => null (), &
 dionstr(:) => null (), &
 dtcunk(:) => null (), &
 um(:) => null (), &
 dum(:) => null ()
logical                     :: &
 isdiv, &
 isconv, &
 isupmxiter, &
 isupmxitergam, &
 isanomalous, &
 havedc, &
 havedg, &
 havedsktrk, &
 havenchemiter, &
 isderivate, &
 be
character(len=100)         :: &
 msg, &
 nameoutput, &
 name 
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
msg=' '
iserror=.false.
isconvergence=.true.
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
isderivate=.true.
!%-----------------------------------------------------------------
!% Check optional arguments
!%-----------------------------------------------------------------
havedc=present(dc)
havedg=present(dg)
havedsktrk=present(dsktrk)
havenchemiter=present(nchemiter)
!%------------------------------------------------------------------
!%------------------------------------------------------------------
!%------------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 iserror=.true.
 goto 10
end if
!%------------------------------------------------------------------
!% Check number of surfaces
!%------------------------------------------------------------------
if (nsurf/=this%numsurf) then
 msg='Error in number of surfaces'
 iserror=.true. 
 goto 10
end if
!%-----------------------------------------------------------------
!% Check if the phase one mineral
!%-----------------------------------------------------------------
do i=1,this%numph
 call get_if_sp_is_present_ (this%pphase(i)%ptr,nameph,be)
 if (be) exit
end do
if (be) then
 call get_if_mineral_ (this%pphase(i)%ptr,be)
 if (.not.be) then
  msg='Error, the phase is not mineral: '
  call add_ (msg,nameph)
  goto 10
 end if
else
 msg='Error in name of phase: '
 call add_ (msg,nameph)
 goto 10
end if
!%-----------------------------------------------------------------
!% Optional: Check the phases for equilibrium
!%-----------------------------------------------------------------
call get_sp_index_  (this,nameph,imin)
if (imin==0) then
 msg='Error, not found the species: '
 call add_ (msg,nameph)
 goto 10
end if
!%-----------------------------------------------------------------
!% Give the water component and species index
!%-----------------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wcindex=icw,wspindex=ispsw)
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
nunk=this%numaqprisp+1
call check_pointer_ (cold,this%numsp,.true.)
call check_pointer_ (dcloc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dsi,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dgloc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (tk,this%numaqprisp,.true.)
call check_pointer_ (tk1,this%numaqprisp,.true.)
call check_pointer_ (dtm,this%numaqprisp,.true.)
call check_pointer_ (dtk1,nunk,nunk,.true.)
call check_pointer_ (jacobian,nunk,nunk,.true.)
call check_pointer_ (residual,nunk,.true.)
call check_pointer_ (indx,nunk,.true.)
call check_pointer_ (dtcunk,nunk,.true.)
call check_pointer_ (cunk,nunk,.true.)
call check_pointer_ (um,this%numaqprisp,.true.)
call check_pointer_ (dum,this%numaqprisp,.true.)
cold=c
!%----------------------------------------------------------------
!% Compute T in k
!%----------------------------------------------------------------
tk=matmul(this%ueq,cold)
!%----------------------------------------------------------------
!% Fill dc
!%----------------------------------------------------------------
do i=1,this%numaqprisp
 dcloc(this%idaqprisp(i),i)=1.0d0
end do
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
nstep=0
itry=0
siphloc=10.0d0**(siph)
amountph=1.0d0
!%----------------------------------------------------------------
!% Compute the water in hydrated minerals
!%----------------------------------------------------------------
if (ispsw>0) then
 call compute_omgwcryst_(this,omgwcrystk,omgwfreek,c,nsp,msg,iserror)
 if (iserror) goto 20
end if
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
 
do
 
 
itry=itry+1
nstep=nstep+1
amountph=amountph/real(nstep)
 
do istep=1,nstep
!%----------------------------------------------------------------
!% Compute increment of mass of water
!%----------------------------------------------------------------
if (ispsw>0) then
 dtomg=this%ueq(icw,imin)*amountph*kgwmol
 omg=omgwfreek+omgwcrystk+dtomg
end if
!%----------------------------------------------------------------
!% Compute the variation of mass of components
!%----------------------------------------------------------------
dtm=this%ueq(:,imin)*amountph
!%----------------------------------------------------------------
!% Initialice cunk 
!%----------------------------------------------------------------
cunk(1:this%numaqprisp) = c(this%idaqprisp)
cunk(nunk)=(amountph/(omgwfreek+dtomg))/2.0d0
!%----------------------------------------------------------------
c=0.0d0
!%----------------------------------------------------------------
c(this%idaqprisp)=dabs(cunk(1:this%numaqprisp))
!%----------------------------------------------------------------
c(imin)=cunk(nunk)
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!% Loop for Newton Raphson
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
iter = 0
idiverg = 0
 
do
 
  iter=iter+1

  jacobian=0.0d0
  residual=0.0d0
  tk1=0.0d0
  dtk1=0.0d0
!%----------------------------------------------------------------
!% Specia aqueous phase and aqueous phase with another phases
!%----------------------------------------------------------------
  call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .true., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam) goto 20
!%----------------------------------------------------------------
!% Specia aqueous phase and aqueous phase with another phases
!%----------------------------------------------------------------
si=c
dsi=dcloc
do i=1,this%numminph
  call compute_secondaries_ &
     (this, &
      si, &
      g, &
      dsi, &
      dgloc, &
      this%aqphindex, &
      this%idminph(i), &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      .true., &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam) goto 20
end do
!%---------------------------------------------------------------
!%---------------------------------------Specia surface complexes
!%---------------------------------------------------------------
do i=1,this%numsurf
 
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    txoh(:,i), &
    c1(:,i), &
    c2(:,i), &
    spsurfarea(:,i), &
    ntxoh, &
    ionstr, &
    i, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror) goto 20
 
end do
 
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
a=1.0d0+kgwmol*this%ueq(icw,imin)*c(imin)
a2=a*a
omgwfreek1=omg/a
domgwfreek1=-omg*kgwmol*this%ueq(icw,imin)/a2
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
um=this%ueq(:,imin)*c(imin)
dum=this%ueq(:,imin)
!%-------------------------------------------------------------------
!%------------------------------------Build the jacobian and residual
!%-------------------------------------------------------------------
tk1=matmul(this%ueq,c)
dtk1=matmul(this%ueq,dcloc)
!%-------------------------------------------------------------------
residual(1:this%numaqprisp)=omgwfreek1*tk1-omgwfreek*tk-dtm
jacobian(1:this%numaqprisp,1:this%numaqprisp)=omgwfreek1*dtk1
!%-------------------------------------------------------------------
jacobian(1:this%numaqprisp,nunk)=domgwfreek1*um+omgwfreek1*dum
!%-------------------------------------------------------------------
residual(nunk)=si(imin)-siphloc
jacobian(nunk,1:this%numaqprisp)=dsi(imin,:)
!%--------------------------------------------------------------------
!%--------------------------Optional: Write Newton Raphson Information
!%--------------------------------------------------------------------
 if (this%iswriteinfo) then
 
  call write_newton_raphson_info_ &
   (jacobian, &
    residual, &
    nunk, &
    nunk, &
    iter, &
    iouwiinfo, &
    'specia_')

 end if
!%----------------------------------------------------------------------
!% Solve the lineal system
!%----------------------------------------------------------------------
 residual=-residual
 dtcunk=residual
 call ludcmp (jacobian,nunk,nunk,indx,dd,msg,iserror)
 if (iserror) goto 20 
 call lubksb (jacobian,nunk,nunk,indx,dtcunk)
!%----------------------------------------------------------------------
!% Check convergence and update the solution
!%----------------------------------------------------------------------
 call check_convergence_ &
   (this, &
    isconv, &
    isdiv, &
    isupmxiter, &
    cunk(1:this%numaqprisp), &
    dtcunk(1:this%numaqprisp), &
    residual, &
    nunk, &
    this%numaqprisp, &
    errmax, &
    iter, &
    idiverg, &
    this%tolunknr, &
    this%tolresnr)
 
    cunk = cunk + dtcunk
 
    c(this%idaqprisp)=cunk(1:this%numaqprisp)

    c(imin)=cunk(nunk)
 
    if (isconv.or.isdiv.or.isupmxiter) exit
 
end do
!%-------------------------------------------------------------------
!%-----------------------------------Finished loop for Newton Raphson
!%-------------------------------------------------------------------
if (.not.isconv.or.isdiv.or.isupmxiter) exit
 
end do
 
if (isconv.and.c(imin)<0.0d0) then
 amountph = amountph + 1.0d0
else if (isconv.and.c(imin)>0.0d0) then
 exit
end if
 
if (itry>mxtry.and..not.isconv) then
 msg='Error, convergence problems'
 isconvergence=.false.
 goto 20
end if
 
end do
!%------------------------------------------------------------------
!%----------------------------------------------Finished loop of try
!%------------------------------------------------------------------
omgwfreek=omgwfreek1
!%------------------------------------------------------------------
!%------------------------------------------------------------------
!%------------------------------------------------------------------
c(imin)=0.0d0
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
if (.not.havedc &
        .or. &
    .not.havedg &
        .or. &
   .not.havedsktrk) isderivate=.false.
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%--------------------------------------------------Update speciation
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .true., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
 
   if (iserror.or.isupmxitergam) goto 20
!%------------------------------------------------------------------
!%---------------------------Compute saturation indexes for minerals
!%------------------------------------------------------------------
si=c
do i=1,this%numminph
 
 call compute_secondaries_ &
     (this, &
      si, &
      g, &
      dcloc, &
      dgloc, &
      this%aqphindex, &
      this%idminph(i), &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      .false., &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam) goto 20
 
end do
!%------------------------------------------------------------------
!%---------------------------------Make the speciation in gas phases
!%------------------------------------------------------------------
do i=1,this%numgasph
 
  call compute_secondaries_ &
    (this, &
     c, &
     g, &
     dcloc, &
     dgloc, &
     this%aqphindex, &
     this%idgasph(i), &
     ionstr, &
     dionstr, &
	 1.0d0, &
	 .false., &
     isanomalous, &
     isupmxitergam, &
     isderivate, &
     msg, &
     iserror)
 
    if (iserror.or.isupmxitergam) goto 20
 
end do
!%------------------------------------------------------------------
!%--------------------------Make the speciation of surface complexes
!%------------------------------------------------------------------
do i=1,this%numsurf
 
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    txoh(:,i), &
    c1(:,i), &
    c2(:,i), &
    spsurfarea(:,i), &
    ntxoh, &
    ionstr, &
    i, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror) goto 20
 
end do
!%------------------------------------------------------------------
!%----------------------------------------------------------Optional
!%------------------------------------------------------------------
if (havedc) then
 call check_pointer_ (dc,this%numsp,this%numaqprisp,.true.)
 dc=dcloc
end if
!%------------------------------------------------------------------
!%----------------------------------------------------------Optional
!%------------------------------------------------------------------
if (havedg) then
 call check_pointer_ (dg,this%numsp,this%numaqprisp,.true.)
 dg=dgloc
end if
!%------------------------------------------------------------------
!%-------------------------------------Deallocate vectors and arrays
!%------------------------------------------------------------------
20 continue 
call check_pointer_ (dsi,1,1,.false.)
call check_pointer_ (tk1,1,.false.)
call check_pointer_ (tk,1,.false.)
call check_pointer_ (dtk1,1,1,.false.)
call check_pointer_ (cold,1,.false.)
call check_pointer_ (dcloc,1,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (cunk,1,.false.)
call check_pointer_ (dtcunk,1,.false.)
call check_pointer_ (indx,1,.false.)
call check_pointer_ (idsp,1,.false.)
call check_pointer_ (jacobian,1,1,.false.)
call check_pointer_ (residual,1,.false.)
call check_pointer_ (um,1,.false.)
call check_pointer_ (dum,1,.false.)
call check_pointer_ (dtm,1,.false.)
!%-----------------------------------------------------------------
if (havenchemiter) then 
 nchemiter=iter
end if 
!%-----------------------------------------------------------------
if (iserror) goto 10
!%-----------------------------------------------------------------
!%------------------------------------------------------------------
!%------------------------------------------------------------------
return
2  format(<nunk>e10.3,a2,e10.3)
3  format(a20,a2,e10.3) 
10 continue
print *,'************************************'
print *,'Chemical System:'
print *,'Name:', this%name
print *,'Service: specia_'
print *, msg
print *,'Maximun error in unknown=', errmax
print *,'Iterations=', iter
print *,'Tolerance in unknown=', this%tolunknr
print *,'Tolerance in residual=', this%tolresnr
print *,'************************************'
if (associated(cold)) then
 c=cold
 call check_pointer_ (cold,1,.false.)
end if 
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine equilibratebis_pchemsys &
   (this, &
    c, &
    g, &
    simin, &
    ionstr, &
    nameph, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    nsp, &
    nsurf, &
    ntxoh, &
    omgwfreek, &
    omgwcrystk, &
    isnonconvergence, &
    iserror, &
    siph, &
    amountph, &
    dc_ext, &
    sktrk_ext, &
    dsktrk_ext, &
    dg_ext)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

real*8, intent(inout)                     :: omgwfreek ! mass of water per unit of volume in k

real*8, intent(inout)                     :: omgwcrystk             ! mass of water per unit of volume tetained

real*8, intent (out)                      :: ionstr

integer, intent(in)                       :: nsp

integer, intent(in)                       :: nsurf

integer, intent(in)                       :: ntxoh

character(len=*), intent(in)              :: nameph         ! Name of species added

real*8, intent(inout)                     :: c(nsp) ! concentration vector

real*8, intent(inout)                     :: g(nsp) ! activity coefficients vector

real*8, intent(out)                       :: simin(nsp)

real*8, intent(in)                        :: txoh(ntxoh,nsurf)

real*8, intent(in)                        :: c1(ntxoh,nsurf)

real*8, intent(in)                        :: c2(ntxoh,nsurf)

real*8, intent(in)                        :: spsurfarea(ntxoh,nsurf)

logical, intent(out)                      :: iserror

logical, intent(out)                      :: isnonconvergence

real*8, pointer, optional                 :: dc_ext(:,:)

real*8, pointer, optional                 :: sktrk_ext(:)

real*8, pointer, optional                 :: dsktrk_ext(:,:)

real*8, pointer, optional                 :: dg_ext(:,:)

real*8, optional                          :: siph

real*8, optional                          :: amountph 
 
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
 i, &
 j, &
 k, &
 iter, &
 nunk, &
 idiverg, &
 ispsw, &
 icw, &
 ipos1, &
 ipos2, &
 npri, &
 isp
integer, pointer            :: &
 indx(:), &
 idsp(:)
real*8                      :: &
 dd, &
 errmax, &
 omg, &
 omgwfreek1, &
 dtomg, &
 siphloc, &
 amountphloc, &
 a, &
 a2, &
 domgwfreek1
real*8, pointer             :: &
 cunk(:), &
 cold(:), &
 tk1(:), &
 tk(:), &
 residual(:), &
 dc(:,:), &
 dsi(:,:), &
 dg(:,:), &
 dtk1(:,:), &
 dtm(:), &
 jacobian(:,:), &
 dionstr(:), &
 dtcunk(:), &
 um(:), &
 dum(:)
logical                     :: &
 isdivergence, &
 isconvergence, &
 isupmxiter, &
 isupmxitergam, &
 isanomalous, &
 havedc, &
 havedg, &
 havedsktrk, &
 isderivate, &
 havesiph, &
 haveamountph, &
 be
character(len=100)         :: &
 msg, &
 nameoutput, &
 name 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
msg=' '
iserror=.false.
!%-----------------------------------------------------------------
!% Initialice variables 
!%-----------------------------------------------------------------
isnonconvergence=.false.
isderivate=.true.
!%-----------------------------------------------------------------
!% Check optional arguments
!%-----------------------------------------------------------------
havedc=present(dc_ext)
havedg=present(dg_ext)
havedsktrk=present(dsktrk_ext)
havesiph=present(siph)
haveamountph=present(amountph)
!%------------------------------------------------------------------
!% Check number of species 
!%------------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 iserror=.true.
 goto 10
end if
!%------------------------------------------------------------------
!% Check number of surfaces
!%------------------------------------------------------------------
if (nsurf/=this%numsurf) then
 msg='Error in number of surfaces'
 goto 10
end if
!%-----------------------------------------------------------------
!% Check if the phase one mineral
!%-----------------------------------------------------------------
do i=1,this%numph
 call get_if_sp_is_present_ (this%pphase(i)%ptr,nameph,be)
 if (be) exit
end do
if (be) then
 call get_if_mineral_ &
   (this%pphase(i)%ptr, &
    be)
 if (.not.be) then
  msg='Error, non mineral phase:'
  call add_ (msg,nameph)
  goto 10
 end if
else
 msg='Error, chemical species not defined in the chemical system:'
 call add_ (msg,nameph)
 goto 10
end if
!%-----------------------------------------------------------------
!% Optional: Check the phases for equilibrium
!%-----------------------------------------------------------------
siphloc=1.0d0
amountphloc=10.0d0
call get_sp_index_  (this,nameph,isp)
if (isp.eq.0) then
 msg='Error, not found the species:'
 call add_ (msg,nameph)
 goto 10
end if
if (havesiph) then
 siphloc=siph
end if
if (haveamountph) then
 amountphloc=amountph
end if
!%-----------------------------------------------------------------
!% Give the water component and species index
!%-----------------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wcindex=icw,wspindex=ispsw)
!%-----------------------------------------------------------------
!% Compute number of unknowns 
!%-----------------------------------------------------------------
nunk=this%numaqprisp+1
!%-------------------------------------------------------------------
!% Allocate local pointers 
!%-------------------------------------------------------------------
call check_pointer_ (cold,this%numsp,.true.)
call check_pointer_ (dc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dsi,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dg,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (tk,this%numaqprisp,.true.)
call check_pointer_ (tk1,this%numaqprisp,.true.)
call check_pointer_ (dtm,this%numaqprisp,.true.)
call check_pointer_ (dtk1,nunk,nunk,.true.)
call check_pointer_ (jacobian,nunk,nunk,.true.)
call check_pointer_ (residual,nunk,.true.)
call check_pointer_ (indx,nunk,.true.)
call check_pointer_ (dtcunk,nunk,.true.)
call check_pointer_ (cunk,nunk,.true.)
call check_pointer_ (um,this%numaqprisp,.true.)
call check_pointer_ (dum,this%numaqprisp,.true.)
cold=c
c=0.0d0
!%----------------------------------------------------------------
!% Compute T in k
!%----------------------------------------------------------------
tk=matmul(this%ueq,cold)
!%----------------------------------------------------------------
!% Fill c and dc
!%----------------------------------------------------------------
do i=1,this%numaqprisp
 dc(this%idaqprisp(i),i)=1.0d0
end do
!%----------------------------------------------------------------
!% Compute increment of mass of water
!%----------------------------------------------------------------
if (ispsw>0) then
 dtomg=this%ueq(icw,isp)*amountphloc*kgwmol
 omg=omgwfreek+omgwcrystk+dtomg
end if
!%----------------------------------------------------------------
!% Compute the variation of mass of components
!%----------------------------------------------------------------
dtm=this%ueq(:,isp)*amountphloc
!%----------------------------------------------------------------
!%----------------------------------------------------------------
c(ispsw)=1.0d0/kgwmol
c(this%idaqprisp)=cold(this%idaqprisp)
c(isp)=amountphloc/omgwfreek
cunk(1:this%numaqprisp)=cold(this%idaqprisp)
cunk(nunk)=c(isp)
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!% Loop for Newton Raphson
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
iter = 0
 
do
 
  iter=iter+1
 
  jacobian=0.0d0
  residual=0.0d0
  tk1=0.0d0
  dtk1=0.0d0
!%----------------------------------------------------------------
!% Specia aqueous phase and aqueous phase with another phases
!%----------------------------------------------------------------
  call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dc, &
      dg, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .true., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam) goto 20 
!%----------------------------------------------------------------
!% Specia aqueous phase and aqueous phase with another phases
!%----------------------------------------------------------------
simin=c
dsi=dc
do i=1,this%numminph
  call compute_secondaries_ &
     (this, &
      simin, &
      g, &
      dsi, &
      dg, &
      this%aqphindex, &
      this%idminph(i), &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      .true., &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam) goto 20
end do
!%---------------------------------------------------------------
!% Specia surface complexes
!%---------------------------------------------------------------
do i=1,this%numsurf
 
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dc, &
    dg, &
    txoh(:,i), &
    c1(:,i), &
    c2(:,i), &
    spsurfarea(:,i), &
    ntxoh, &
    ionstr, &
    i, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror) goto 20
 
end do
 
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
omgwfreek1=omgwfreek+cunk(icw)*kgwmol-1.0d0
domgwfreek1=kgwmol
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
um=this%ueq(:,isp)*c(isp)
dum=this%ueq(:,isp)
!%-------------------------------------------------------------------
!% Build the jacobian and residual
!%-------------------------------------------------------------------
tk1=matmul(this%ueq,c)
dtk1=matmul(this%ueq,dc)
!%-------------------------------------------------------------------
residual(1:this%numaqprisp)=omgwfreek1*tk1-omgwfreek*tk-dtm
jacobian(1:this%numaqprisp,1:this%numaqprisp)=omgwfreek1*dtk1
!%-------------------------------------------------------------------
jacobian(1:this%numaqprisp,nunk)=omgwfreek1*dum
!%-------------------------------------------------------------------
residual(nunk)=simin(isp)-siphloc
jacobian(nunk,1:this%numaqprisp)=dsi(isp,:)
jacobian(1:this%numaqprisp,icw)=jacobian(1:this%numaqprisp,icw)+domgwfreek1*tk1
!%--------------------------------------------------------------------
!% Optional: Write Newton Raphson Information
!%--------------------------------------------------------------------
 if (this%iswriteinfo) then
 
  call write_newton_raphson_info_ &
   (jacobian, &
    residual, &
    nunk, &
    nunk, &
    iter, &
    iouwiinfo, &
    'specia_')
 
 
 
2   format(<nunk>e10.3,a2,e10.3)
3   format(a20,a2,e10.3)
 end if
!%----------------------------------------------------------------------
!% Solve the lineal system
!%----------------------------------------------------------------------
 residual=-residual
 dtcunk=residual
 call ludcmp (jacobian,nunk,nunk,indx,dd,msg,iserror)
 if (iserror) goto 20 
 call lubksb (jacobian,nunk,nunk,indx,dtcunk)
!%----------------------------------------------------------------------
!% Check convergence and update the solution
!%----------------------------------------------------------------------
 call check_convergence_ &
   (this, &
    isconvergence, &
    isdivergence, &
    isupmxiter, &
    cunk(1:this%numaqprisp), &
    dtcunk(1:this%numaqprisp), &
    residual, &
    nunk, &
    this%numaqprisp, &
    errmax, &
    iter, &
    idiverg, &
    this%tolunknr, &
    this%tolresnr)
 
    cunk = cunk + dtcunk
 
      c(this%idaqprisp)=cunk(1:this%numaqprisp)
    c(isp)=cunk(nunk)
 
      if (isconvergence.or.isdivergence.or.isupmxiter) exit
 
end do
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!% Finished loop for Newton Raphson
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
 omgwfreek=omgwfreek1
omgwcrystk=omg-omgwfreek1
!%------------------------------------------------------------------
c(isp)=0.0d0
!%------------------------------------------------------------------
!%------------------------------------------------------------------
!%------------------------------------------------------------------
if (isdivergence) then
 msg='Convergence problems, divergence detected'
 isnonconvergence=.true.
 goto 10
end if
!%------------------------------------------------------------------
!%------------------------------------------------------------------
!%------------------------------------------------------------------
if (isupmxiter) then
 msg='Convergence problems, maximum number of iterations exceeded'
 isnonconvergence=.true.
 goto 10
end if
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
if (.not.havedc &
        .or. &
    .not.havedg &
        .or. &
   .not.havedsktrk) isderivate=.false.
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!% Update speciation
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
call compute_secondaries_ &
     (this, &
      c, &
      g, &
      dc, &
      dg, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .true., &
      isanomalous, &
      isupmxitergam, &
      isderivate, &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam) goto 20
!%------------------------------------------------------------------
!% Compute saturation indexes for minerals
!%------------------------------------------------------------------
simin=c
do i=1,this%numminph
 
 call compute_secondaries_ &
     (this, &
      simin, &
      g, &
      dc, &
      dg, &
      this%aqphindex, &
      this%idminph(i), &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      .false., &
      msg, &
      iserror)
 
    if (iserror.or.isupmxitergam) goto 20
 
end do
!%------------------------------------------------------------------
!% Make the speciation in gas phases
!%------------------------------------------------------------------
do i=1,this%numgasph
 
  call compute_secondaries_ &
    (this, &
     c, &
     g, &
     dc, &
     dg, &
     this%aqphindex, &
     this%idgasph(i), &
     ionstr, &
     dionstr, &
	 1.0d0, &
	 .false., &
     isanomalous, &
     isupmxitergam, &
     isderivate, &
     msg, &
     iserror)
 
    if (iserror.or.isupmxitergam) goto 20
 
end do
!%------------------------------------------------------------------
!% Make the speciation of surface complexes
!%------------------------------------------------------------------
do i=1,this%numsurf
 
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dc, &
    dg, &
    txoh(:,i), &
    c1(:,i), &
    c2(:,i), &
    spsurfarea(:,i), &
    ntxoh, &
    ionstr, &
    i, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror) goto 20
 
end do
!%------------------------------------------------------------------
!% Optional
!%------------------------------------------------------------------
if (havedc) then
 call check_pointer_ (dc_ext,this%numsp,this%numaqprisp,.true.)
 dc_ext=dc
end if
!%------------------------------------------------------------------
!% Optional
!%------------------------------------------------------------------
if (havedg) then
 call check_pointer_ (dg_ext,this%numsp,this%numaqprisp,.true.)
 dg_ext=dg
end if
!%------------------------------------------------------------------
20 continue 
!%------------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------------
call check_pointer_ (dsi,1,1,.false.)
call check_pointer_ (tk1,1,.false.)
call check_pointer_ (tk,1,.false.)
call check_pointer_ (dtk1,1,1,.false.)
call check_pointer_ (cold,1,.false.)
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dg,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (cunk,1,.false.)
call check_pointer_ (dtcunk,1,.false.)
call check_pointer_ (indx,1,.false.)
call check_pointer_ (jacobian,1,1,.false.)
call check_pointer_ (residual,1,.false.)
call check_pointer_ (um,1,.false.)
call check_pointer_ (dum,1,.false.)
if (iserror) goto 10 
!%-----------------------------------------------------------------
!%------------------------------------------------------------------
!%------------------------------------------------------------------
return
 
10 continue 
print *,'************************************'
print *,'Chemical System:'
print *,'Name:', this%name
print *,'Service: specia_'
print *, msg
print *,'Maximun error in unknown=', errmax
print *,'Iterations=', iter
print *,'Tolerance in unknown=', this%tolunknr
print *,'Tolerance in residual=', this%tolresnr
print *,'************************************'
if (associated(cold)) then 
 c=cold
 call check_pointer_ (cold,1,.false.)
end if 
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine specia_from_cpri_pchemsys &
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
    c1, &
    c2, &
    supsurfarea, &
    mxtxoh, &
    nsurf, &
    dtime, &
	volgas, &
    ionstr, &
	faccap, &
    isanomalous, &
    isupmxitergam, &
    reset, &
    iserror, &
    dc, &
    dg, &
    sktrk, &
    dsktrk)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Make the speciation from concentration of primary 
! aqueous species and compute derivates and kinetic terms (optional)
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)   :: this

real*8, intent(in)                          :: temp

integer, intent(in)                         :: nsp             ! Number of species 

integer, intent(in)                         :: npri            ! Number of primary species 

integer, intent(in)                         :: mxtxoh

integer, intent(in)                         :: nsurf

real*8, intent(inout), dimension(nsp)       :: c 

real*8, intent(inout), dimension(nsp)       :: g 

real*8, intent(in), dimension(nsp)          :: cold

real*8, intent(in), dimension(nsp)          :: alpha

real*8, intent(in), dimension(npri)         :: cpri

real*8, intent(in), dimension(mxtxoh,nsurf) :: txoh

real*8, intent(in), dimension(mxtxoh,nsurf) :: c1

real*8, intent(in), dimension(mxtxoh,nsurf) :: c2

real*8, intent(in), dimension(mxtxoh,nsurf) :: supsurfarea

real*8, intent(in)                          :: dtime           ! Time increment 

real*8, intent(in)                          :: volgas

logical, intent(in)                         :: reset

real*8, intent(out)                         :: ionstr

real*8, intent(in)                          :: faccap          ! Capillary correction for water activity 

logical, intent(out)                        :: iserror         ! iserror=true, then there was an error 

logical, intent(out)                        :: isanomalous     ! isanomalous, there was an anomalous concentration 

logical, intent(out)                        :: isupmxitergam

real*8, pointer, optional, dimension(:,:)   :: dc

real*8, pointer, optional, dimension(:,:)   :: dg

real*8, pointer, optional, dimension(:)     :: sktrk

real*8, pointer, optional, dimension(:,:)   :: dsktrk
 
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
 i, &
 j, &
 nmobph, &
 ipos1, &
 ipos2, &
 iminpri(1)
real*8                       :: &
 minpri
real*8, pointer              :: &
 dcloc(:,:) => null (), &
 dgloc(:,:) => null (), &
 dionstr(:) => null ()
logical                      :: &
 havedc, &
 havedg, &
 havesktrk, &
 havedsktrk, &
 isderivate
character(len=100)           :: &
 msg 
real*8, parameter            :: &
 r0=0.0d0, &
 r1=1.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------------------
!% Check optional arguments
!%-----------------------------------------------------------------------
havedc=present(dc)
havedg=present(dg)
havesktrk=present(sktrk)
havedsktrk=present(dsktrk)
!%-----------------------------------------------------------------------
!% Check number of species 
!%-----------------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%-----------------------------------------------------------------------
!% Check number of surfaces 
!%-----------------------------------------------------------------------
if (nsurf/=this%numsurf) then
 msg='Error in number of surfaces'
 goto 10
end if
!%-----------------------------------------------------------------------
!% Check the number of components 
!%-----------------------------------------------------------------------
if (npri/=this%numaqprisp) then
 msg='Error in number of components'
 goto 10
end if
!%-----------------------------------------------------------------------
!% Check negative primary concentrations
!%-----------------------------------------------------------------------
minpri=minval(cpri)
iminpri=minloc(cpri)
if (minpri<=r0) then
 msg='Error, anomalous concentrations in primary species:'
 call add_ (msg,this%pspecies(this%idaqprisp(iminpri(1)))%ptr%name)
 goto 10
end if
!%-----------------------------------------------------------------------
isderivate=.false.
isupmxitergam=.false.
!%-----------------------------------------------------------------------
if (havedc.or.havedg) then
 call check_pointer_ (dcloc,this%numsp,this%numaqprisp,.true.)
 call check_pointer_ (dgloc,this%numsp,this%numaqprisp,.true.)
 do i=1,this%numaqprisp
  dcloc(this%idaqprisp(i),i)=r1
 end do
 isderivate=.true.
end if
!%----------------------------------------------------------------------
!% Reset c
!%----------------------------------------------------------------------
if (reset) c=r0
!%----------------------------------------------------------------------
c(this%idaqprisp)=cpri
!%----------------------------------------------------------------------
!% For aqueous phase
!%----------------------------------------------------------------------
call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    this%aqphindex, &
    0, &
    ionstr, &
    dionstr, &
	faccap, &
	.true., &
    isanomalous, &
    isupmxitergam, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror.or.isupmxitergam.or.isanomalous) goto 20
!%----------------------------------------------------------------------- 
!% For gas phases
!%-----------------------------------------------------------------------
do i=1,this%numgasph
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    this%aqphindex, &
    this%idgasph(i), &
    ionstr, &
    dionstr, &
	1.0d0, &
	.false., &
    isanomalous, &
    isupmxitergam, &
    isderivate, &
    msg, &
    iserror)
 
  if (iserror.or.isupmxitergam.or.isanomalous) goto 20
 
  call get_iposspsph_ (this,this%idgasph(i),ipos1,ipos2)
 
  c(ipos1:ipos2)=c(ipos1:ipos2)*alpha(ipos1:ipos2)
  do j=1,this%numaqprisp
     dcloc(ipos1:ipos2,j)=dcloc(ipos1:ipos2,j)*alpha(ipos1:ipos2)
  end do
 
end do
!%-----------------------------------------------------------------------
!% For surfaces
!%-----------------------------------------------------------------------
do i=1,this%numsurf
 
 call compute_secondaries_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    txoh(:,i), &
    c1(:,i), &
    c2(:,i), &
    supsurfarea(:,i), &
    mxtxoh, &
    ionstr, &
    i, &
    isderivate, &
    msg, &
    iserror)
 
    if (iserror) goto 20
 
end do
!%-----------------------------------------------------------------------
!% Compute kinetic
!%-----------------------------------------------------------------------
call compute_kinetic_ &
   (this, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    alpha, &
    cold, &
    dtime, &
    msg, &
    iserror, &
    sktrk=sktrk, &
    dsktrk=dsktrk)
    if (iserror) goto 20
!%-----------------------------------------------------------------------
!% Check NAN numbers in concentration and activity coefficients vectors
!%-----------------------------------------------------------------------    
do i=1,this%numsp 
 if (isnan(c(i)).or.isnan(g(i))) then
  isanomalous=.true.
  goto 20
 end if
end do     
!%-----------------------------------------------------------------------
!% If have dg
!%-----------------------------------------------------------------------
if (havedg) then
 call check_pointer_ (dg,nsp,npri,.true.)
 dg=dgloc
end if
!%-----------------------------------------------------------------------
!% If have dc
!%-----------------------------------------------------------------------
if (havedc) then
 call check_pointer_ (dc,nsp,npri,.true.)
 dc=dcloc
end if
!%---------------------------------------------------------------
20 continue 
!%---------------------------------------------------------------
!% Deallocate local pointers 
!%---------------------------------------------------------------
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (dcloc,1,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
if (iserror) goto 10
!%-------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:',this%name
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
subroutine compute_dumob1_pchemsys &
   (this, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute Ujth*dcmob/dc1
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)         :: this

integer, intent(in)                               :: numsp

real*8, intent(in), dimension(numsp)              :: c 

real*8, pointer, dimension(:,:)                   :: dumob

integer, intent(out)                              :: nrow

integer, intent(out)                              :: ncol

integer, intent(out)                              :: nmobph

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
integer                               :: &
 i, &
 j, &
 ipos1, &
 ipos2, &
 ipri1, &
 ipri2
real*8, pointer                       :: &
 dc(:,:) => null (), &
 dg(:,:) => null (), &
 g(:) => null (), &
 cloc(:) => null (), &
 dionstr(:) => null ()
real*8                                :: &
 ionstr
character(len=100)                    :: &
 msg
logical                               :: &
 isanomalous, &
 isupmxitergam 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
if (numsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------
nrow=this%numaqprisp
ncol=this%numaqprisp
nmobph=1+this%numgasph
!%------------------------------------------------------
call check_pointer_ (dumob,nmobph*nrow,ncol,.true.)
!%------------------------------------------------------
call check_pointer_ (dc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dg,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (g,this%numsp,.true.)
call check_pointer_ (cloc,this%numsp,.true.)
cloc=c
g=1.0d0
!%-------------------------------
call get_iposspsph_ (this, this%aqphindex, ipos1, ipos2)
!%-----------------------------by aqueous phase
do i=1,this%numaqprisp
 dc(this%idaqprisp(i),i)=1.0d0
end do
!%-------------------------------------------------------
call compute_secondaries_ &
   (this, &
    cloc, &
    g, &
    dc, &
    dg, &
    this%aqphindex, &
    0, &
    ionstr, &
    dionstr, &
	1.0d0, &
	.true., &
    isanomalous, &
    isupmxitergam, &
    .true., &
    msg, &
    iserror)
 
if (iserror.or.isupmxitergam) goto 20
 
ipri1=1
ipri2=nrow
dumob(ipri1:ipri2,:)= &
matmul(this%ueq(:,ipos1:ipos2),dc(ipos1:ipos2,:))
 
 
!%----------------------------------------by gas phase
ipri1=nrow+1
do i=1,this%numgasph
 
 call compute_secondaries_ &
   (this, &
    cloc, &
    g, &
    dc, &
    dg, &
    this%aqphindex, &
    this%idgasph(i), &
    ionstr, &
    dionstr, &
	1.0d0, &
	.false., &
    isanomalous, &
    isupmxitergam, &
    .true., &
    msg, &
    iserror)
 
 if (iserror.or.isupmxitergam) goto 20
 
 call get_iposspsph_ (this, this%idgasph(i), ipos1, ipos2)
 ipri2=ipri2+nrow
 
 dumob(ipri1:ipri2,:)= matmul(this%ueq(:,ipos1:ipos2),dc(ipos1:ipos2,:))
 
 ipri1=ipri2+1
 
end do
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------------------
call check_pointer_ (g,1,.false.)
call check_pointer_ (cloc,1,.false.)
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dg,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_dumob_'
print *,msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dumob2_pchemsys &
   (this, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    ndimder, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute U*dcmob
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)         :: this

integer, intent(in)                               :: numsp

integer, intent(in)                               :: ndimder

real*8, intent(in), dimension(numsp,ndimder)      :: dc 

real*8, pointer, dimension(:,:)                   :: dumob

integer, intent(out)                              :: nrow

integer, intent(out)                              :: ncol

integer, intent(out)                              :: nmobph

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
integer                               :: &
 i, &
 j, &
 isp1, &
 isp2, &
 ipri1, &
 ipri2
real*8                                :: &
 ionstr
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
!% Check the number of species 
!%------------------------------------------------------
if (numsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------
nrow=this%numaqprisp
ncol=ndimder
nmobph=1+this%numgasph
!%------------------------------------------------------
call check_pointer_ (dumob,nmobph*nrow,ncol,.true.)
!%------------------------------------------------------
!% for aqueous phase
!%------------------------------------------------------
call get_iposspsph_ (this, this%aqphindex,isp1,isp2)
ipri1=1
ipri2=nrow
dumob(ipri1:ipri2,:)= matmul(this%ueq(:,isp1:isp2),dc(isp1:isp2,:))
!%------------------------------------------------------
!% by gas phase
!%------------------------------------------------------
ipri1=nrow+1
do i=1,this%numgasph
 
 call get_iposspsph_ (this, this%idgasph(i), isp1, isp2)
 ipri2=ipri2+nrow
 dumob(ipri1:ipri2,:)= matmul(this%ueq(:,isp1:isp2),dc(isp1:isp2,:))
 ipri1=ipri2+1
 
end do
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_dumob_'
print *,msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dcmob1_pchemsys &
   (this, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    iserror, &
    dg, &
    g)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute dcmob from c
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)   :: this

integer, intent(in)                         :: numsp

real*8, intent(in)                          :: c(numsp)

real*8, pointer                             :: dcmob(:,:)

integer, intent(out)                        :: nrow

integer, intent(out)                        :: ncol

integer, intent(out)                        :: nmobph

logical, intent(out)                        :: iserror

real*8, pointer, optional                   :: dg(:,:)

real*8, pointer, optional                   :: g(:) 
 
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
 i, &
 j, &
 ipos1, &
 ipos2, &
 isp1, &
 isp2
real*8, pointer                       :: &
 dc(:,:) => null (), &
 dgloc(:,:) => null (), &
 dionstr(:) => null (), &
 gloc(:) => null (), &
 cloc(:) => null () 
real*8                                :: &
 ionstr
character(len=100)                    :: &
 msg
logical                               :: &
 isanomalous, &
 havedg, &
 haveg, &
 isupmxitergam 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
havedg=present(dg)
haveg=present(g)
!%------------------------------------------------------
if (numsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------
nrow=this%numsp
ncol=this%numaqprisp
nmobph=1+this%numgasph
!%------------------------------------------------------
call check_pointer_ (dcmob,nmobph*nrow,ncol,.true.)
!%------------------------------------------------------
call check_pointer_ (gloc,this%numsp,.true.)
call check_pointer_ (cloc,this%numsp,.true.)
call check_pointer_ (dc,this%numsp,this%numaqprisp,.true.)
call check_pointer_ (dgloc,this%numsp,this%numaqprisp,.true.)
cloc=c
gloc=1.0d0
!%-------------------------------
call get_iposspsph_ (this, this%aqphindex, ipos1, ipos2)
!%-----------------------------by aqueous phase
do i=1,this%numaqprisp
 dc(this%idaqprisp(i),i)=1.0d0
end do
!%---------------------------------------------
 
call compute_secondaries_ &
   (this, &
    cloc, &
    gloc, &
    dc, &
    dgloc, &
    this%aqphindex, &
    0, &
    ionstr, &
    dionstr, &
	1.0d0, &
	.true., &
    isanomalous, &
    isupmxitergam, &
    .true., &
    msg, &
    iserror)
 
 if (iserror.or.isupmxitergam) goto 20
 
isp1=1
isp2=nrow
dcmob(isp1:isp2,:)=dc(ipos1:ipos2,:)
!%----------------------------------------by gas phase
isp1=nrow+1
do i=1,this%numgasph
 
 call compute_secondaries_ &
   (this, &
    cloc, &
    gloc, &
    dc, &
    dgloc, &
    this%aqphindex, &
    this%idgasph(i), &
    ionstr, &
    dionstr, &
	1.0d0, &
	.false., &
    isanomalous, &
    isupmxitergam, &
    .true., &
    msg, &
    iserror)
 
    if (iserror.or.isupmxitergam) goto 20
 
call get_iposspsph_ (this, this%idgasph(i), ipos1, ipos2)
 
 
isp2=isp2+nrow
 
dcmob(isp1:isp2,:)=dc(ipos1:ipos2,:)
 
isp1=isp2+1
 
end do
!%-----------------------------------------------------------
if (havedg) then
 call check_pointer_ (dg,this%numsp,this%numaqprisp,.true.)
 dg=dgloc
end if
!%-----------------------------------------------------------
if (haveg) then
 call check_pointer_ (g,this%numsp,.true.)
 g=gloc
end if
!%-----------------------------------------------------------
20 continue 
!%------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (gloc,1,.false.)
call check_pointer_ (cloc,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_dcmob_'
print *,msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_alkalinity_pchemsys &
   (this, &
    alk, &
    c, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:Compute alkalinity from concentrations vector
!
!   In general the alcalinity is defined as
!   alkalinity = 2*[co3-2] + [hco3-]
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in), target :: this

integer, intent(in)                               :: nsp

real*8, intent(in), dimension(nsp), target        :: c 

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
integer                                           :: &
 isps, &
 isps1, &
 isps2
character(len=100)                                :: &
 msg
character(len=100), pointer                       :: &
 namesps => null ()
real*8, pointer                                   :: &
 conc => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
alk=0.0d0
!%------------------------------------------------------
!% Check the number of species 
!%------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10 
end if 
!%------------------------------------------------------
!%------------------------------------------------------
!%------------------------------------------------------
if (this%aqphindex==0) return
!%------------------------------------------------------
call get_iposspsph_ (this,this%aqphindex,isps1,isps2)
!%------------------------------------------------------
do isps=isps1,isps2
 namesps => this%pspecies(isps)%ptr%name
 conc => c(isps)
 select case (namesps)
 case ('hco3-','cahco3+','mghco3+','nahco3','nahco3(aq)')
  alk = alk + conc
 case ('co3-2','co3-2(aq)','caco3','caco3(aq)','mgco3','mgco3(aq)','naco3-')
  alk = alk + 2.0d0 * conc
 case ('oh-')
  alk = alk + conc
 case ('h+')
  alk = alk - conc
 end select
end do
!%------------------------------------------------------
!% Nullify local pointers 
!%------------------------------------------------------
namesps => null ()
conc => null ()
!%------------------------------------------------------
return
 
10 continue 
print *,'****************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_alkalinity_'
print *,msg
print *,'****************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dcmob2_pchemsys &
   (this, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    nsp, &
    npri, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute dcmob/dc1 from dc
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

integer, intent(in)                       :: nsp

integer, intent(in)                       :: npri

real*8, intent(in), dimension(nsp,npri)   :: dc

real*8, pointer, dimension(:,:)           :: dcmob 

integer, intent(out)                      :: nrow

integer, intent(out)                      :: ncol

integer, intent(out)                      :: nmobph

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
 i, &
 j, &
 ipos1, &
 ipos2, &
 isp1, &
 isp2
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 
 

 

!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------
nrow=this%numsp
ncol=npri
nmobph=1+this%numgasph
!%------------------------------------------------------
call check_pointer_ (dcmob,nmobph*nrow,ncol,.true.)
!%------------------------------------------------------
call get_iposspsph_ (this, this%aqphindex, ipos1, ipos2)
!%-------------------------------------------by aqueous phase
isp1=ipos1
isp2=ipos2
dcmob(isp1:isp2,:)=dc(ipos1:ipos2,:)
!%----------------------------------------------by gas phases
do i=1,this%numgasph
 call get_iposspsph_ (this,this%idgasph(i),ipos1,ipos2)
 isp1=i*nrow+ipos1
 isp2=i*nrow+ipos2
 dcmob(isp1:isp2,:)=dc(ipos1:ipos2,:)
end do
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_dcmob_'
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
subroutine compute_dcads_pchemsys &
   (this, &
    dcads, &
    nrow, &
    ncol, &
    dc, &
    nsp, &
    npri, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute dcmob/dc1
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

integer, intent(in)                       :: nsp

integer, intent(in)                       :: npri

real*8, intent(in), dimension(nsp,npri)   :: dc

real*8, pointer, dimension(:,:)           :: dcads 

integer, intent(out)                      :: nrow

integer, intent(out)                      :: ncol

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
 isurf, &
 ipos1, &
 ipos2
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------
nrow=this%numsp
ncol=npri
!%------------------------------------------------------
call check_pointer_ (dcads,nrow,ncol,.true.)
!%------------------------------------------------------
do isurf=1,this%numsurf
 call get_iposspsurf_ (this,isurf,ipos1,ipos2)
 dcads(ipos1:ipos2,:)=dc(ipos1:ipos2,:)
end do
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_dcads_'
print *,msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_secondaries_pchemsys &
   (this, &
    c, &
    g, &
    dc, &
    dg, &
    iph1, &
    iph2, &
    ionstr, &
    dionstr, &
	faccap, &
	isiterative, &
	isanomalous, &
    isupmxiter, &
    isderivates, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Solve the action mass law for homogeneous and heterogeneous reactions. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)            :: this

real*8, intent(inout), dimension(this%numsp), target :: c            ! Vector of concentration of species

real*8, intent(inout), dimension(this%numsp), target :: g            ! Activity coefficients vector

real*8, intent(inout), dimension(:,:)                :: dc           ! Derivatives of the concentration of species with respect to primary species

real*8, intent(inout), dimension(:,:)                :: dg           ! Derivatives of the activity coefficients with respect to primary species

integer, intent(in)                                  :: iph1         ! Global index of the phase 1

integer, intent(in)                                  :: iph2         ! Global index of the phase 2

real*8, intent(out)                                  :: ionstr       ! Ionic strength 

real*8, pointer, dimension(:)                        :: dionstr      ! Derivatives of the ionic strength 

real*8, intent(in)                                   :: faccap       ! Capillary correction for water activity 

logical, intent(in)                                  :: isiterative  ! isiterative=true, is applied the Picard method

logical, intent(out)                                 :: isanomalous  ! isanomalous=true, then there was an anomalous concentration 
 
logical, intent(out)                                 :: isupmxiter   ! 

logical, intent(out)                                 :: iserror      ! iserror=true, then there was an error
 
logical, intent(in)                                  :: isderivates  ! isderivates=true, then compute the derivatives. 

character(len=*), intent(out)                        :: msg          ! Message error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer &
 i, &
 j, &
 iph, &
 numreact, &
 ireact, &
 iter, &
 isps1, &
 isps2, &
 isps, &
 icw, &
 ispsw
real*8, pointer             :: &
 c2, &
 g2
real*8                      :: &
 mxerror, &
 cw, &
 error, &
 c2old 
logical                     :: &
 isconvergence, &
 isbeaqph
real*8, pointer             :: &
 gloc(:) => null (), &
 dcloc(:) => null (), &
 dgloc(:,:) => null ()
integer, pointer            :: &
 idreaction(:) => null () 
type(t_reaction), pointer   :: &
 reaction => null ()
type(t_phase), pointer      :: &
 phase => null ()
real*8, parameter           :: &
 tolrel=1.0d-5, &    ! Relative tolerance 
 tolabs=1.0d-10, &   ! Absolute tolerance
 r0=0.0d0, &
 r1=1.0d0
logical, parameter    :: &
 isnumeric=.false.
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------------
!% Initialice variables
!%--------------------------------------------------------------------
isanomalous=.false.
isupmxiter=.false. 
isbeaqph=.false.
!%--------------------------------------------------------------------
!% 1) Determine the reactions between phases iph1 and iph2
!%--------------------------------------------------------------------
call get_idreaction_(this,idreaction,numreact,iph1,iph2,.false.)
!%--------------------------------------------------------------------
!% If iph2=0, then 
!%--------------------------------------------------------------------
if (iph2==0) then
  iph=iph1
else
  iph=iph2
end if
!%--------------------------------------------------------------------
!% Get general information in the parent chemical system 
!%--------------------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wcindex=icw,wspindex=ispsw)
if (ispsw>0) then
 cw=c(ispsw)
 c(ispsw)=r1
end if
!%--------------------------------------------------------------------
!% Get the global position of species defined in the chemical system
!%--------------------------------------------------------------------
call get_iposspsph_ (this,iph,isps1,isps2)
if(iph==this%aqphindex) isbeaqph=.true.
!%--------------------------------------------------------------------
phase => this%pphase(iph)%ptr
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!% Start Picard to compute secondary species
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
iter=0
isconvergence=.false. 
do
 
 
 mxerror=r0
 iter=iter+1

!%-------------------------------------------------------
 do i=1,numreact
 
  ireact = idreaction (i)
  reaction => this%preaction(ireact)%ptr
  isps = this%idreactsp(ireact)
  g2 => g(isps)
  c2 => c(isps)
  c2old = c2
 
  call compute_x_ (reaction,c2,c(this%idaqprisp),g(this%idaqprisp),g2,iserror) 
  if (iserror) goto 10

 
!%--------------------------------------------------------------------
!% Check anomalous concentrations (only for aqueous phase)
!%--------------------------------------------------------------------
  if ((c2>maxanomc.or.c2<minanomc).and.isbeaqph) then
     isanomalous=.true.
     msg='Anomalous concentrations in the species:'
     call add_ (msg,this%pspecies(isps)%ptr%name)
     write(6,*) msg
     goto 10
  end if


if (isderivates) then 
      call compute_dx_ &
	  (reaction, &
	   dcloc, &
	   c(this%idaqprisp), &
       g(this%idaqprisp), &
       g2, &
       dc(this%idaqprisp,:), &
       dg(this%idaqprisp,:), &
       dg(isps,:), &
       iserror, &
       c2)
  
      if (iserror) goto 20
   
      if (icw>0) dcloc(icw)=r0
 
      dc(isps,:) = dcloc
end if 


 
  error=dabs((c2-c2old)/c2) 
  if (error>mxerror) mxerror=error
  
 end do
 
!%---------------------------------------------------------------------
!% Compute activity coefficients
!%---------------------------------------------------------------------
  call compute_act_coeff_(phase,gloc,c(isps1:isps2),iserror,param=ionstr,factor=faccap)
  
  if (iserror) goto 10 
 
  g(isps1:isps2) = gloc

  
if (isderivates) then 
     
	 if (isnumeric) then   ! Compute numeric derivatives of activity coefficients 
		
		call compute_dact_coeff_ &
        (phase, &
         dgloc, &
         c(isps1:isps2), &
         dc(isps1:isps2,:), &
         1.0d-2, &
         iserror, &
         dparam=dionstr, &
		 factor=faccap)
      
	  else		
		
		call compute_dact_coeff_ &   ! Compute analitic derivatives 
        (phase, &
         dgloc, &
         c(isps1:isps2), &
         dc(isps1:isps2,:), &
         iserror, &
         g=g(isps1:isps2), &
         dparam=dionstr, &
		 factor=faccap)
        
      end if 	
      
	  if (iserror) goto 10
      
	  if (icw>0) then
          dgloc(:,icw)=r0
          dionstr(icw)=r0
      end if
 
      dg(isps1:isps2,:) = dgloc 


end if 

!%---------------------------------------------------------------------
!% Check convergence 
!%---------------------------------------------------------------------   
  if (isiterative.and.iter>1) then
   isconvergence=(mxerror<=tolrel)
  else if (.not.isiterative) then
   isconvergence=.true. 
  end if

  if (isconvergence) exit
  
  isupmxiter=(iter>mxitergam)
 
  if (isupmxiter) then
     msg='Non convergence in gamma iterations'
     write(6,*) msg
     if (this%iswriteinfo) write(iouwiinfo,*) msg
     goto 10
  end if

  
 
end do
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
if (ispsw>0) c(ispsw)=cw
!%--------------------------------------------------------------------
10 continue 
!%--------------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------------
call check_pointer_ (idreaction,1,.false.)
call check_pointer_ (gloc,1,.false.)
call check_pointer_ (dcloc,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
reaction => null ()
phase => null ()
if (iserror) goto 20
!%-------------------------------------------------------------------
return
 
20 iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_adsorption_pchemsys &
   (this, &
    c, &
    g, &
    dc, &
    dg, &
    txoh, &
    capint, &
    capext, &
    spsurfarea, &
    mxnsite, &
    ionstr, &
    ithsurf, &
    isderivates, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Solve adsorption from activity of primary species and make the mass balance of xoh
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)               :: this

integer, intent(in)                                     :: ithsurf

integer, intent(in)                                     :: mxnsite

real*8, intent(inout), dimension(this%numsp), target    :: c             ! Concentration vector. 

real*8, intent(in), dimension(this%numsp), target       :: g             ! Activity coefficient vector. 

real*8, intent(in), dimension(mxnsite), target          :: txoh

real*8, intent(in), dimension(mxnsite), target          :: capint

real*8, intent(in), dimension(mxnsite), target          :: capext

real*8, intent(in), dimension(mxnsite), target          :: spsurfarea

real*8, intent(inout), dimension(:,:)                   :: dc

real*8, intent(in), dimension(:,:)                      :: dg

real*8, intent(in)                                      :: ionstr        ! Ionic strength 

logical, intent(in)                                     :: isderivates

logical, intent(out)                                    :: iserror       ! iserror=true, then there was an error

character(len=*), intent(out)                           :: msg           ! Error message. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                      :: &
 isps1, &
 isps2, &
 nsite, &
 nreact, &
 i, &
 j, &
 ireact, &
 iter, &
 isps, &
 nsk, &
 nsp, &
 icw, &
 ispsw, &
 ifail, &
 npri, &
 irank  
real*8, pointer                             :: &
 cd, &
 gd 
real*8                                      :: &
 cw, &
 sigma
integer, pointer                            :: &
 idreact(:) => null (), &
 idxoh(:) => null ()
real*8, pointer                             :: &
 sk(:) => null (), &
 dcloc(:) => null (), &
 cloc(:) => null (), &
 txohloc(:) => null (), &
 capintloc(:) => null (), &
 capextloc(:) => null (), &
 spsurfarealoc(:) => null (), &
 dcdsk(:,:) => null (), &
 b(:,:) => null (), &
 a(:,:) => null (), &
 work(:) => null ()
logical                                     :: &
 isbeactive, &
 isconvergence, &
 isupmxiter
type(t_reaction), pointer                   :: &
 reaction => null ()
type(t_surface), pointer                    :: &
 surface => null ()
real*8, parameter                           :: &
 r0=0.0d0, &
 r1=1.0d0, &
 tol=1.0d-6
!-------------------------------------------------------------------------
!
!   $code
!

!%---------------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------------
if (ithsurf>this%numsurf) then
 msg='Error in surface index, index='
 call add_ (msg,ithsurf) 
 goto 10
end if
!%--------------------------------------------------------------
surface => this%psurface(ithsurf)%ptr
!%--------------------------------------------------------------
!% Check if the surfaces is active
!%--------------------------------------------------------------
call get_numsites_ (surface,nsite)
!%--------------------------------------------------------------
txohloc => txoh(1:nsite)
capintloc => capint(1:nsite)
capextloc => capext(1:nsite)
spsurfarealoc => spsurfarea(1:nsite)
!%--------------------------------------------------------------
call get_if_active_ (surface,txohloc,isbeactive)
if (.not.isbeactive) return
!%--------------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wcindex=icw,wspindex=ispsw)
if (ispsw>0) then
 cw=c(ispsw)
 c(ispsw)=r1
end if
!%--------------------------------------------------------------
!% Get reaction indexes
!%--------------------------------------------------------------
call get_idreaction_ (this,idreact,nreact,this%aqphindex,ithsurf,.true.)
!%-------------------------------------------------------------
call get_iposspsurf_ (this,ithsurf,isps1,isps2)
!%-------------------------------------------------------------
call get_numsp_ (surface,nsp)
!%-------------------------------------------------------------
call get_num_tot_sk_ (surface,nsk)
!%-------------------------------------------------------------
call check_pointer_ (dcdsk,this%numsp,nsk,.true.)
!%-------------------------------------------------------------
call get_idxoh_ (surface,idxoh,nsite,iserror)
if (iserror) goto 20
do i=1,nsite
 isps=isps1-1+idxoh(i)
 dcdsk(isps,r1)=r1 
end do 
!%-------------------------------------------------------------
cloc => c(isps1:isps2)
!%-------------------------------------------------------------
!% Init value for unknowns for adsorption 
!%-------------------------------------------------------------
call init_ (surface,sk,cloc,txohloc,capintloc,spsurfarealoc,ionstr)
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Solve adsorption using Newton Raphson
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
iter=0
!%--------------------------------------------------------------
do
 iter=iter+1
 do i=1,nreact
  ireact=idreact(i)
  reaction => this%preaction(ireact)%ptr
  isps = this%idreactsp(ireact)
  gd => g(isps)
  cd => c(isps)
  call compute_x_ (reaction,cd,c(this%idaqprisp),g(this%idaqprisp),gd,iserror,sk)
  if (iserror) goto 20
  call compute_dx_ (reaction,dcloc,cd,sk,iserror)
  if (iserror) goto 20
  dcdsk(isps,:)=dcloc
 end do
!%-------------------------------------------------------------
!% Compute sk
!%-------------------------------------------------------------
 call update_sk_(surface,sk,cloc,dcdsk(isps1:isps2,1:nsk),txohloc, &
                 capintloc,capextloc,spsurfarealoc,ionstr,isconvergence, &
                 isupmxiter,iter,iserror)
 if (isupmxiter) then
   msg='Error, over number of iterations in adsorption'
   iserror=.true.
   goto 20
  end if
 if (iserror) goto 20
!%-------------------------------------------------------------
!% If there was convergence
!%-------------------------------------------------------------
if (isconvergence) then
 
  do i=1,nreact
      ireact=idreact(i)
	  reaction => this%preaction(ireact)%ptr
      isps = this%idreactsp(ireact)
      gd => g(isps)
      cd => c(isps) 
      call compute_x_ (reaction,cd,c(this%idaqprisp),g(this%idaqprisp),gd,iserror,sk)
      if (iserror) goto 20
      if (isderivates) then   
         call compute_dx_ (reaction,dcloc,cd,sk,iserror)
         if (iserror) goto 20
         dcdsk(isps,:)=dcloc
      end if
  end do
 
  exit
 
end if
!%------------------------------------------------------------
end do
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
! Change to molality only valid for cation exchange model
!%------------------------------------------------------------
call change_to_mol_(surface,cloc,txohloc,iserror)
!%------------------------------------------------------------
!% Compute derivates
!%------------------------------------------------------------
if (isderivates) then
 !npri=this%numaqprisp
 !call check_pointer_(b,nsp,npri,.true.) 
 !call check_pointer_(a,nsp,nsk,.true.) 
 !call check_pointer_(work,4*nsk,.true.)
 
 !do j=1,30

  !do i=1,npri 
   ! ifail=0 
   ! a=dcdsk(isps1:isps2,1:nsk)
   ! b(1:nsp,i)= dc(isps1:isps2,i)
   ! call f04jaf(nsp,nsk,a,nsp,b(1:nsp,i),tol,sigma,irank,work,4*nsk,ifail)
  !end do  
    
  
  do i=1,nreact
 
    ireact = idreact(i)
	reaction => this%preaction(ireact)%ptr
    isps = this%idreactsp(ireact)
    cd => c(isps)
    gd => g(isps)
 
    call compute_dx_ &
     (reaction, &
     dcloc, &
     c(this%idaqprisp), &
     g(this%idaqprisp), &
     gd, &
     dc(this%idaqprisp,:), &
     dg(this%idaqprisp,:), &
     dg(isps,:), &
     iserror, &
     xj=cd)  !  , &
     !sk=sk, &
     !dsk=b(1:nsk,1:npri))
 
    if (iserror) goto 20
 
    dc(isps,:) = dcloc
 
  end do
  
!end do 
 
end if
!%------------------------------------------------------------
if (ispsw>0) c(ispsw)=cw
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate and nullify local pointers
!%------------------------------------------------------------
call check_pointer_ (idreact,1,.false.)
call check_pointer_ (dcloc,1,.false.)
call check_pointer_ (dcdsk,1,1,.false.)
call check_pointer_ (sk,1,.false.)
!call check_pointer_(b,1,1,.false.) 
!call check_pointer_(a,1,1,.false.) 
!call check_pointer_(work,1,.false.) 
call check_pointer_(idxoh,1,.false.) 
reaction => null ()
surface => null ()
cloc => null () 
txohloc => null ()
capintloc => null ()
capextloc => null ()
spsurfarealoc => null ()
if (iserror) goto 10
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
subroutine compute_min_sat_index_pchemsys &
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
!   $Description: Compute saturation indices of mineral species. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)    :: this

real*8, pointer, dimension(:)                :: simin

integer, intent(in)                          :: nsp

real*8, intent(in), dimension(nsp)           :: c

real*8, intent(in), dimension(nsp)           :: g

integer, pointer, dimension(:)               :: idminsp

character(len=*), pointer, dimension(:)      :: nameminsp

integer, intent(out)                         :: nminsp

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
integer                          :: &
 i, &
 j, &
 isp, &
 sum, &
 ipos1, &
 ipos2, &
 numprisp, &
 isp2
real*8,           pointer       :: &
 dc(:,:) => null (), &
 dg(:,:) => null (), &
 dionstr(:) => null (), &
 si(:) => null (), &
 gloc(:) => null () 
real*8                          :: &
 ionstr
character(len=100), pointer     :: &
 namesp(:) => null ()
logical                         :: &
 isanomalous, &
 isupmxitergam
character(len=100)              :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%-------------------------------------------------------
if (nsp/=this%numsp) then
  msg='Error in number of species'
  goto 10
end if
!%-------------------------------------------------------
if (this%numminph==0) return
!%-------------------------------------------------------
call check_pointer_ (si,this%numsp,.true.)
call check_pointer_ (gloc,this%numsp,.true.)
!%-------------------------------------------------------
!% Compute saturation indices for minerals
!%-------------------------------------------------------
sum=0
si=c
gloc=g 
do i=1,this%numminph
 
 call compute_secondaries_ &
   (this, &
    si, &
    gloc, &
    dc, &
    dg, &
    this%aqphindex, &
    this%idminph(i), &
    ionstr, &
    dionstr, &
	1.0d0, &
	.false., &
    isanomalous, &
    isupmxitergam, &
    .false., &
    msg, &
    iserror)
 
 if (iserror.or.isupmxitergam) goto 20
 
 call get_numsp_ (this%pphase(this%idminph(i))%ptr,isp)
 sum=sum+isp
 
end do
!%------------------------------------------------------ 
nminsp=sum
!%------------------------------------------------------
!% Allocate local pointers 
!%-------------------------------------------------------
call check_pointer_ (simin,nminsp,.true.)
call check_pointer_ (nameminsp,nminsp,.true.)
call check_pointer_ (idminsp,nminsp,.true.)
!%------------------------------------------------------
call get_chem_info_ (this,msg,iserror,namesp=namesp)
isp=0
isp2=0
do i=1,this%numminph
 call get_iposspsph_ (this, this%idminph(i), ipos1, ipos2)
 do j=ipos1,ipos2
  isp2=isp2+1
  idminsp (isp2) = j
 end do
 simin (isp+1:isp + (ipos2-ipos1)+1) = si (ipos1:ipos2)
 nameminsp (isp+1:isp + (ipos2-ipos1+1)) = namesp (ipos1:ipos2)
 isp = isp + ipos2 - ipos1 + 1
end do
 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
20 continue 
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (si,1,.false.)
call check_pointer_ (gloc,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dg,1,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_min_sat_index_'
print *, msg
print*, '*******************************'
iserror=.true.
return
 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_pchemsys &
   (this, &
    temp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update parameters depending of the temperature in the 
! chemical system
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout)    :: this     ! Parent chemical system variable

real*8, intent(in)                              :: temp     ! Temperature in celcius.

logical, intent(out)                            :: iserror  ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                          :: &
 i
character(len=100)               :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
! Update the temperature depending parameters in the phases
!%-----------------------------------------------------------
do i=1,this%numph
 call update_ (this%pphase(i)%ptr,temp,iserror)
 if (iserror) goto 10
end do
!%-----------------------------------------------------------
! Update the temperature depending parameters in the surfaces
!%-----------------------------------------------------------
do i=1,this%numsurf
 call update_ (this%psurface(i)%ptr,temp,iserror)
 if (iserror) goto 10
end do
!%----------------------------------------------------------
! Update the temperature depending parameters for reactions
!%----------------------------------------------------------
do i=1,this%numreact
 call update_ (this%preaction(i)%ptr,temp,iserror)
 if (iserror) goto 10 
end do
!%---------------------------------------------------------
 
return
 
10 continue 
print *,'******************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: update_ (temperature)'
print *, msg
print*, '******************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_porosity_pchemsys &
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
!   $Description: Update the porosity according the mineral concentrations
!%
!%
!%    por = porosity (in/out)
!%    c = vector of concentrations in k+1
!%    cold = vector of concentrations in k
!%    nsp = number of species
!%    omgwfreek = mass of water free per unit of volume in k
!%    omgwfreek1 = mass of water free per unit of volume in k+1
!%    iserror = .true. there was error
!%
!%   delta(m3min)= sum( (cminj(k+1)*omgwfree(k+1) -
!%                       cminj(k+1)*omgwfree(k))* VMj*1.0d-3 )
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
 
type (t_parentchemicalsystem), intent(in), target    :: this        ! Type chemical system variable

integer, intent(in)                                  :: nsp         ! Number of species

real*8, intent(in), dimension(nsp)                   :: ck1         ! Concentration of species in k+1

real*8, intent(in), dimension(nsp)                   :: ck          ! Concentration of species in k

real*8, intent(inout)                                :: porosity    ! Porosity 

real*8, intent(in)                                   :: omgwfreek1  ! Mass of free water in k+1

real*8, intent(in)                                   :: omgwfreek   ! Mass of free water in k

logical, intent(out)                                 :: iserror     ! iserror=true, then there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                      :: &
 i, &
 j, &
 isps1, &
 isps2
real*8                                       :: &
 sum, &
 dvol
character(len=100)                           :: &
 msg 
character(len=100), pointer                  :: &
 name => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%----------------------------------------------------------
iserror=.false.
msg=' '
!%----------------------------------------------------------
!% Check the number of species
!%----------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%----------------------------------------------------------
if (this%numminph==0) return
!%----------------------------------------------------------
sum=0.0d0
do i=1,this%numminph
 call get_iposspsph_ (this, this%idminph(i),isps1,isps2)
 do j=isps1,isps2
   name => this%pspecies(j)%ptr%name
   dvol = ck1 (j) * omgwfreek1 - ck (j) * omgwfreek
   call change_chem_unit_(this,dvol,name,'mol','m3',iserror)
   if (iserror) goto 10
 end do
 sum = sum + dvol
end do
!%----------------------------------------------------------
!% Update the porosity
!%----------------------------------------------------------
porosity = porosity - sum
!%----------------------------------------------------------
!% Nullify local pointers 
!%----------------------------------------------------------
name => null ()
!%----------------------------------------------------------
return
 
10 continue 
print *,'*****************************'
print *,'Chemical System:'
print *,'Name:', this%name
print *,'Service: update_ (porosity)'
print *, msg
print *,'*****************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_min_area_pchemsys &
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
!   $Description: Update the reactive surface of kinetic minerals
!
!   $Arguments:
!
 
 
type (t_parentchemicalsystem), intent(in)    :: this     ! Type parent chemical system variable. 

real*8, intent(inout), dimension(this%numsp) :: area     ! Area 

real*8, intent(in), dimension(this%numsp)    :: area0    ! Initial area 

real*8, intent(in), dimension(this%numsp)    :: c        ! Concentration 

real*8, intent(in), dimension(this%numsp)    :: c0       ! Initial concentration 

character(len=*), intent(out)                :: msg      ! Error message. 

logical, intent(out)                         :: iserror  ! iserror=true, then there was an error. 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8              :: &
 vm, &     ! Molar volume 
 area01, &
 area00 
integer             :: &
 ireact, &
 isps 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
msg='' 
iserror=.false. 
!%------------------------------------------------------------
!% Compute minimum default area (according sphere of 
!% r=1.0d-8 [m]
!%------------------------------------------------------------
area00=4.0d0*pi*r0*r0
!%------------------------------------------------------------
!% Only update the area of kinetic minerals 
!%------------------------------------------------------------
do ireact=1,this%numreact
 if (this%iskinreact(ireact)) then ! If the reaction is kinetic
  isps=this%idreactsp(ireact)
  if (isps>0) then ! If the reaction is associated with a species 
      if (area0(isps)>0.0d0) then
	      area01=area0(isps)
	  else
          area01=area00 
	  end if 
	  if (c0(isps)>0.0d0) then 
          area(isps)=area01*(c(isps)/c0(isps))**(2.0d0/3.0d0)
		  if (area(isps)<area00) area(isps)=area00
      else 
          call get_prop_ (this%pspecies(isps)%ptr,vm,'molvol',msg,iserror)
          if (iserror) return 
          vm=vm*1.0d-6 ! cm3_min/mol  => m3_min/mol  
          area(isps)=(((9.0d0*area01)/(r0*r0)) * (vm * c(isps))**2.0d0)**(1.0d0/3.0d0)
          if (area(isps)<area01) area(isps)=area01
      end if
  end if 
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
subroutine update_txoh_pchemsys &
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
 
 
type (t_parentchemicalsystem), intent(in)                    :: this     ! Type parent chemical system variable. 

real*8, intent(inout), dimension(this%numsites,this%numsurf) :: txoh     ! Total sites 

real*8, intent(in), dimension(this%numsites,this%numsurf)    :: txoh0    ! Initial total sites

integer, intent(in), dimension(this%numsurf)                 :: idtxohmin ! Initial total sites

real*8, intent(in), dimension(this%numsp)                    :: c        ! Concentration 

real*8, intent(in), dimension(this%numsp)                    :: c0       ! Initial concentration 
 
character(len=*), intent(out)                                :: msg      ! Error message. 

logical, intent(out)                                         :: iserror  ! iserror=true, then there was an error. 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8              :: &
 vm, &     ! Molar volume 
 txoh01, &
 txoh00 
integer             :: &
 isurf, &
 imin, &
 nsite, &
 isite 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
msg='' 
iserror=.false. 
!%------------------------------------------------------------
!% Compute minimum default area (according sphere of 
!% r=1.0d-8 [m]
!%------------------------------------------------------------
txoh00=4.0d0*pi*r0*r0
!%------------------------------------------------------------
!% Only update the area of kinetic minerals 
!%------------------------------------------------------------
do isurf=1,this%numsurf
  imin=idtxohmin(isurf)
  if (imin>0) then  
      call get_numsites_ (this%psurface(isurf)%ptr,nsite)
	  do isite=1,nsite
	     if (txoh0(isite,isurf)>0.0d0) then
	       txoh01=txoh0(isite,isurf)
	     else
           txoh01=txoh00 
	     end if 
	     if (c0(imin)>0.0d0) then 
            txoh(isite,isurf)=txoh01*(c(imin)/c0(imin))**(2.0d0/3.0d0)
		    if (txoh(isite,isurf)<txoh00) txoh(isite,isurf)=txoh00
         else 
          call get_prop_ (this%pspecies(imin)%ptr,vm,'molvol',msg,iserror)
          if (iserror) return 
          vm=vm*1.0d-6 ! cm3_min/mol  => m3_min/mol  
          txoh(isite,isurf)=(((9.0d0*txoh01)/(r0*r0)) * (vm * c(imin))**2.0d0)**(1.0d0/3.0d0)
          if (txoh(isite,isurf)<txoh01) txoh(isite,isurf)=txoh01
         end if
      end do
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
subroutine change_chem_unit_pchemsys &
   (this, &
    c, &
    namesp, &
    unitin, &
    unitout, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Change units used for chemical system. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)    :: this

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
integer                          :: &
 isps
real*8                           :: &
 value
character(len=100)               :: &
 msg
type(t_species), pointer         :: &
 pspecies => null () 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%--------------------------------------------------------
!% Find the global species indice 
!%--------------------------------------------------------
call get_sp_index_(this,namesp,isps)
if (isps==0) then
 msg='Error not found the species:' 
 call add_ (msg,namesp)
 goto 10
end if
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
pspecies => this%pspecies(isps)%ptr
!%--------------------------------------------------------
 select case (unitout)
!%--------------------------------------------------------
!% 
!%--------------------------------------------------------
 case ('m3','m3/m3rock')
      select case (unitin)
      case ('m','mol','mol/m3rock')
         call get_prop_ (pspecies,value,'molvol',msg,iserror)
         if (iserror) goto 10
         c = c * value * 1.0d-6 
      case default
         msg='Error, not recognized unit(in):'
         call add_ (msg,unitin)
         goto 10
      end select
!%--------------------------------------------------------
!% 
!%--------------------------------------------------------
 case ('ppm','ppm/m3rock')
      select case (unitin)
      case ('m','mol','mol/m3rock')
         call get_prop_ (pspecies,value,'molweight',msg,iserror)
         if (iserror) goto 10
         c = c * value * 1.0d3
      case default
         msg='Error, not recognized unit(in):'
         call add_ (msg,unitin)
         goto 10
      end select
!%--------------------------------------------------------
!% 
!%--------------------------------------------------------
 case ('gr','gr/m3rock')
      select case (unitin)
      case ('m','mol','mol/m3rock')
         call get_prop_ (pspecies,value,'molweight',msg,iserror)
         if (iserror) goto 10
         c = c * value
      case default
         msg='Error, not recognized unit(in):'
         call add_ (msg,unitin)
         goto 10
      end select
!%--------------------------------------------------------
!% 
!%--------------------------------------------------------
 case ('kgr','kgr/m3rock')
      select case (unitin)
      case ('m','mol','mol/m3rock')
         call get_prop_ (pspecies,value,'molweight',msg,iserror)
         if (iserror) goto 10
         c = c * value * 1.0d-3 
      case default
         msg='Error, not recognized unit(in):'
         call add_ (msg,unitin)
         goto 10
      end select
!%--------------------------------------------------------
!% 
!%--------------------------------------------------------
 case ('bq','bq/m3rock')
      select case (unitin)
      case ('m','mol','mol/m3rock')
         call get_prop_ (pspecies,value,'t1_2',msg,iserror)
         if (iserror) goto 10
         c = c * value * 3.1536d6 * avogadro
      case default
         msg='Error, not recognized unit(in):'
         call add_ (msg,unitin)
         goto 10
      end select
!%--------------------------------------------------------
!% 
!%--------------------------------------------------------
 case ('m','mol','mol/m3rock')
      select case (unitin)
      case ('m3','m3/m3rock')
         call get_prop_ (pspecies,value,'molvol',msg,iserror)
         if (iserror) goto 10
         ! value = molar volume in cm3_min/mol 
         c = c * 1.0d6 / value
      case ('gr','gr/m3rock')
         call get_prop_ (pspecies,value,'molweight',msg,iserror)
         if (iserror) goto 10
         c = c / value
      case ('ppm','ppm/m3rock')
         call get_prop_ (pspecies,value,'molweight',msg,iserror)
         if (iserror) goto 10
         c = c / (1.0d3*value)
      case ('m','mol','mol/m3rock')
         return
      case default
         msg='Error, not recognized unit(in):'
         call add_ (msg,unitin)
         goto 10
      end select
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
 case default
      msg='Error, not recognized unit(out):'
      call add_ (msg,unitout)
      goto 10
 end select
!%-----------------------------------------------------
!% Nullify local pointers 
!%-----------------------------------------------------
pspecies => null ()
!%-----------------------------------------------------
return
 
10 continue 
print *,'*****************************'
print *,'Chemical System:'
print *,'Name:', this%name
print *,'Service: change_chem_unit_'
print *, msg
print *,'*****************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_usktrk1_pchemsys &
   (this, &
    usktrk, &
    ncomp, &
    c, &
    g, &
    alpha, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute U*Skt*rk 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)         :: this

integer, intent(out)                              :: ncomp

integer, intent(in)                               :: nsp

real*8, intent(in), dimension(nsp)                :: c

real*8, intent(in), dimension(nsp)                :: g

real*8, intent(in), dimension(nsp)                :: alpha

real*8, pointer, dimension(:)                     :: usktrk

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
real*8                                :: &
 rk
real*8, pointer                       :: &
 sktrk(:) => null (), &
 stqi(:) => null (), &
 cloc(:) => null (), &
 alphaloc(:) => null ()
character(len=100), pointer           :: &
 namesp(:) => null ()
integer                               :: &
 i, &
 j, &
 isp
character(len=100)                    :: &
 msg 
logical                               :: &
 isbe 
!-------------------------------------------------------------------------
!
!   $code
!

!%---------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------
if (nsp.ne.this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%---------------------------------------------------------
ncomp=this%numaqprisp
call check_pointer_ (usktrk,ncomp,.true.)
call check_pointer_ (sktrk,this%numsp,.true.)
call check_pointer_ (cloc,this%numsp,.true.)
call check_pointer_ (alphaloc,this%numsp,.true.)
cloc=c
alphaloc=alpha
!%------------------------
call get_chem_info_ (this,msg,iserror,namesp=namesp)
!%-------------------------------------------------------
do i=1,this%numreact

   if (this%iskinreact(i)) then
     isp=this%idreactsp(i)
     stqi => this%stq(i,:)
     call compute_rk_ &
      (this%preaction(i)%ptr, &
       rk, &
       cloc(this%idaqprisp), &
       g(this%idaqprisp), &
       cloc, &
       g, &
       ncomp, &
       nsp, &
       namesp, &
       cloc(isp), &
       alphaloc(isp), &
       iserror)
 
      if (iserror) then
        msg='Error when calling compute_rk_'
        goto 20
      end if
 
      sktrk=sktrk+stqi*rk
 
    end if
end do
!%----------------------------------------------------------
usktrk=matmul(this%ueq,sktrk)
!%----------------------------------------------------------
20 continue 
!%----------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------
call check_pointer_ (sktrk,1,.false.)
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (cloc,1,.false.)
call check_pointer_ (alphaloc,1,.false.)
stqi => null ()
if (iserror) goto 10
!%----------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Chemical system:'
print *,'Name:',this%name
print *,'Service: compute_usktrk_'
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
subroutine compute_usktrk2_pchemsys &
   (this, &
    usktrk, &
    ncomp, &
    sktrk, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute U*Skt*rk
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)   :: this

integer, intent(out)                        :: ncomp

integer, intent(in)                         :: nsp

real*8, intent(in), dimension(nsp)          :: sktrk

real*8, pointer, dimension(:)               :: usktrk

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
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%---------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------
if (nsp.ne.this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%---------------------------------------------------------
ncomp=this%numaqprisp
call check_pointer_ (usktrk,ncomp,.true.)
!%----------------------------------------------------------
usktrk=matmul(this%ueq,sktrk)
!%----------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Chemical system:'
print *,'Name:',this%name
print *,'Service: compute_usktrk_'
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
subroutine compute_umob_pchemsys &
   (this, &
    umob, &
    npri, &
    nmobph, &
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
 
type (t_parentchemicalsystem), intent(in)    :: this

integer, intent(out)                         :: nmobph

integer, intent(out)                         :: npri

integer, intent(in)                          :: nsp

real*8, intent(in), dimension(nsp)           :: c

real*8, pointer, dimension(:,:)              :: umob

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
integer                               :: &
 i, &
 ipos1, &
 ipos2
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
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------
nmobph=1+this%numgasph
npri=this%numaqprisp
call check_pointer_ (umob,npri,nmobph,.true.)
!%------------------------------------------------------------
call get_iposspsph_ (this, this%aqphindex, ipos1, ipos2)
umob(:,1)=matmul(this%ueq(:,ipos1:ipos2),c(ipos1:ipos2))
!%--------------------------------------------------------------
do i=1,this%numgasph
  call get_iposspsph_ (this,this%idgasph(i),ipos1,ipos2)
  umob(:,1+i)= matmul(this%ueq(:,ipos1:ipos2),c(ipos1:ipos2))
end do
!%------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_umob_'
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
subroutine compute_uads_pchemsys &
   (this, &
    uads, &
    npri, &
    c, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute Uads*cads
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)         :: this

integer, intent(out)                              :: npri

integer, intent(in)                               :: nsp

real*8, intent(in), dimension(nsp)                :: c

real*8, pointer, dimension(:)                     :: uads

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
integer                               :: &
 isps1, &
 isps2, &
 i
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%--------------------------------------------------------------
msg=''
iserror=.false.
!%--------------------------------------------------------------
!% Check the number of species
!%--------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%--------------------------------------------------------------
call check_pointer_ (uads,this%numaqprisp,.true.)
npri=this%numaqprisp
!%--------------------------------------------------------------
if (this%numsurf/=0) then
 do i=1,this%numsurf
  call get_iposspsurf_ (this,i,isps1, isps2)
  uads=uads+matmul(this%ueq(:,isps1:isps2),c(isps1:isps2))
 end do
end if
!%--------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Chemical System:'
print *,'Name:',this%name
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
subroutine compute_ith_total_pchemsys &
   (this, &
    itotal, &
    nametotal, &
    c, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute totals according name of components
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)    :: this

integer, intent(in)                          :: nsp

real*8, intent(in), dimension(nsp)           :: c 

real*8, intent(out)                          :: itotal

logical, intent(out)                         :: iserror

character(len=*), intent(in)                 :: nametotal 
 
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
logical                                      :: &
 isbe
integer                                      :: &
 isp, &
 ipos1, &
 ipos2, &
 ipri, &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
msg=' '
iserror=.false.
!%------------------------------------------------------------
isbe=.true.
itotal=0.0d0
!%------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------
do i=1,this%numaqprisp
 isp=this%idaqprisp(i)
 if (nametotal==this%pspecies(isp)%ptr%name) then
  isbe=.true.
  ipri=i
  exit
 end if
end do
!%------------------------------------------------------------
if (.not.isbe) then
 msg='Error finding aqueous component'
 goto 10
end if
!%------------------------------------------------------------
call get_iposspsph_ (this, this%aqphindex, ipos1, ipos2)
itotal=dot_product (this%ueq(ipri,ipos1:ipos2),c(ipos1:ipos2))
!%------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_iuaq_'
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
subroutine copy_pchemsys &
   (copied, &
    this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy parent chemical system object 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

type (t_parentchemicalsystem), intent(out):: copied 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer        :: &
 i, &
 ndim1, &
 ndim2 
type(t_phase), pointer :: &
 phase(:) => null ()
type(t_reaction), pointer :: &
 reaction(:) => null ()
type(t_surface), pointer :: &
 surface(:) => null ()
character(len=100)       :: &
 namecomp
logical                  :: &
 iserror 
!-------------------------------------------------------------------------
!
!   $code
!
if (this%numph>0) then
 allocate (phase(this%numph))
 do i=1,this%numph
  call create_ (phase(i))
  phase(i)=this%pphase(i)%ptr
 end do
end if

if (this%numreact>0) then
 allocate (reaction(this%numreact))
 do i=1,this%numreact
  call create_ (reaction(i))
  reaction(i)=this%preaction(i)%ptr
 end do
end if

if (this%numsurf>0) then
 allocate (surface(this%numsurf))
 do i=1,this%numsurf
  call create_ (surface(i))
  surface(i)=this%psurface(i)%ptr
 end do
end if

select case (this%isgaussjordan)
case (.true.) 
 namecomp='gj'
case (.false.) 
 namecomp='sv'
end select


call set_ &
   (copied, &
    this%name, &
    this%tempref, &
    phase, &
    surface, &
    reaction, &
    this%numph, &
    this%numsurf, &
    this%numreact, &
    namecomp, &
    this%tolunknr, &
    this%tolresnr, &
    this%deltasatmin, &
    iserror)

if (this%numph>0) then
 do i=1,this%numph
  call destroy_ (phase(i))
 end do
 deallocate (phase)  
end if

if (this%numreact>0) then
 do i=1,this%numreact
  call destroy_ (reaction(i))
 end do
 deallocate (reaction)
end if

if (this%numsurf>0) then
 do i=1,this%numsurf
  call destroy_ (surface(i))
 end do
 deallocate(surface)
end if

copied%iswriteinfo=this%iswriteinfo
copied%iswcompbal=this%iswcompbal
copied%iswriteinfo=this%iswriteinfo

return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_r_from_stqtr_pchemsys &
   (this, &
    stqtr, &
    namesp, &
    nsp, &
    idreact, &
    nreact, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute r solving the following system equations.
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)      :: this

integer, intent(in)                            :: nsp  ! Number of species

integer, intent(in)                            :: nreact     ! Number of reactions

real*8, intent(inout), dimension(nsp)          :: stqtr             ! St r (in) / r (out)

integer, intent(in), dimension(nreact)         :: idreact ! Global indices of reactions

logical, intent(out)                           :: iserror

character(len=*), intent(out)                  :: msg

character(len=*), intent(in), dimension(nsp)   :: namesp
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer               :: &
 a(:,:) => null (), &
 work(:) => null (), &
 b(:) => null ()
real*8                        :: &
 sigma, &
 coeff
integer                       :: &
 ifail, &
 irank, &
 ireact, &
 i, &
 isp
character(len=100)            :: &
 name
real*8, parameter             :: &
 tol=1.0d-6, &
 r0=0.0d0  
!-------------------------------------------------------------------------
!
!   $code
!

!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check the number of equations 
!%-------------------------------------------------------------
if (nreact>nsp) then
 msg='Error in number of reactions'
 goto 10
end if
!%-------------------------------------------------------------
!% Allocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (a,nsp,nreact,.true.)
call check_pointer_ (b,nsp,.true.)
!%--------------------------------------------------------------
!% Find the equilibrium reactions 
!%--------------------------------------------------------------
do i=1,nreact
 
    ireact=idreact(i)
    do isp=1,nsp
     name=namesp(isp)
     call get_coeff_stq_(this%preaction(ireact)%ptr,name,coeff)
     a(isp,i)=coeff
    end do
 
end do
!%-------------------------------------------------------------
call check_pointer_ (work,4*nreact,.true.)
!%-------------------------------------------------------------
!%Solve using 
!%-------------------------------------------------------------
ifail=0
b=stqtr
call f04jaf(nsp,nreact,a,nsp,b,tol,sigma,irank,work,4*nreact,ifail)
stqtr=r0
stqtr(1:nreact)=b(1:nreact)
!%------------------------------------------------------------
!%Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (a,1,1,.false.)
call check_pointer_ (b,1,.false.)
call check_pointer_ (work,1,.false.)
!%-------------------------------------------------------------
return
 
10 iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_speciation_pchemsys &
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
!   $Description: Write in txt the speciation information. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)    :: this

integer, intent(in)                          :: nsp

integer, intent(in)                          :: ioutput

real*8, intent(in)                           :: ionstr

real*8, intent(in)                           :: temp

real*8, intent(in)                           :: omgwfree   ! Mass of free water 

real*8, intent(in), dimension(nsp)           :: c          ! Concentration vector 

real*8, intent(in), dimension(nsp)           :: g          ! Activity coefficients vector 

logical, intent(out)                         :: iserror    ! iserror=true, there was an error. 

character(len=*), intent(out)                :: msg        ! Massage error. 

real*8, intent(in), dimension(nsp), optional :: simin

real*8, intent(in), dimension(nsp), optional :: sktrk
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                     :: &
 ih, &
 i, &
 j, &
 icw, &
 ispsw, &
 isps1, &
 isps2, &
 naqsp
real*8                                      :: &
 gamma, &
 conc, &
 act, &
 mol, &
 ionstrloc, &
 alk, &
 tot, &
 chgbal, &
 tds, &
 density  
logical                                     :: &
 havesimin, &
 havesktrk, &
 isupmxitergam, &
 isanomalous, &
 iswritekin, &
 ismolality 
character(len=40)                           :: &
 name, &
 namecomp 
real*8, pointer                             :: &
 u(:) => null (), &
 t(:) => null (), &
 siloc(:) => null (), &
 gloc(:) => null (), &
 dc(:,:) => null (), &
 dg(:,:) => null (), &
 dionstrloc(:) => null () 
real*8, parameter                      :: &
 zero=1.0d-20, &
 r0=0.0d0, &
 r1=1.0d0  
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Initialice variables
!%-----------------------------------------------------------
isupmxitergam=.false.
isanomalous=.false. 
ismolality=.false.
!%-----------------------------------------------------------
!% Check optional argument
!%-----------------------------------------------------------
havesimin=present(simin)
havesktrk=present(sktrk)
!%-----------------------------------------------------------
!% Determine if There are changes due kinetic reactions
!%-----------------------------------------------------------
iswritekin=.false.
if (havesktrk) then 
 if (sum(dabs(sktrk))/=r0) iswritekin=.true. 
end if 
!%-----------------------------------------------------------
!% Get the proton index 
!%-----------------------------------------------------------
name='h+'
call get_sp_index_ (this,name,ih)
!%-----------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wcindex=icw,wspindex=ispsw)
if (iserror) goto 10
!%-----------------------------------------------------------
!% Compute totals 
!%-----------------------------------------------------------
call check_pointer_ (u,this%numaqprisp,.true.)
call check_pointer_ (t,this%numaqprisp,.true.)
call get_iposspsph_ (this,this%aqphindex,isps1,isps2)
u=matmul(this%u(:,isps1:isps2),c(isps1:isps2))
t=matmul(this%u,c)
!%-----------------------------------------------------------
!% Compute the Total Solid Dissolved of the solution 
!%-----------------------------------------------------------
call get_numsp_ (this%pphase(this%aqphindex)%ptr,naqsp)
call compute_mass_ (this%pphase(this%aqphindex)%ptr,tds,c(isps1:isps2),naqsp,iserror)
if (iserror) goto 20 
!%-----------------------------------------------------------
!% Compute density 
!%-----------------------------------------------------------
!call compute_density_ (this%pphase(this%aqphindex)%ptr,r0,temp,c(isps1:isps2),&
!                       density,ismolality,iserror)
!%-----------------------------------------------------------
!% Write temperature, ionic strength, mass of free water, 
!% water activity, pH and alkalinity
!%-----------------------------------------------------------
write(ioutput,9) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '-------------------------------'
write(ioutput,3) 'Temperature [C] = ',temp 
write(ioutput,4) 'Ionic strength [mol/kg water] = ',ionstr
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
if (tds/=r0) then 
     write(ioutput,3) 'TDS [kg/kg water] = ',tds 
end if
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!if (density/=r0) then 
!     write(ioutput,3) 'Density[kg/m3 water] = ',density
!end if
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
write(ioutput,4) 'Mass of water [kg] = ',omgwfree
if (ispsw>0) then
 write(ioutput,3) 'Water activity = ',g(ispsw)
end if
if (ih>0) then
 write(ioutput,3) 'pH = ',-dlog10(c(ih)*g(ih))
end if
call compute_alkalinity_ (this,alk,c,nsp,iserror)
if (iserror) goto 10
if (alk/=r0) then
 write(ioutput,4) 'Alkalinity [mol/kg water] = ',alk
end if
!%-----------------------------------------------------------
!% Compute and write charge balance in the aqueous phase. 
!%-----------------------------------------------------------
call compute_charge_balance_ (this%pphase(this%aqphindex)%ptr,chgbal,c(isps1:isps2),iserror)
if (iserror) goto 20 
write(ioutput,3) 'Charge balance = ',chgbal 
!%------------------------------------------------------------------
!%Write speciation information fo species organized by components
!%------------------------------------------------------------------
write(ioutput,9) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '-------------------------------'
write(ioutput,5) 'Species','Mol','Molality', 'Gamma', 'Activity', &
                 'log Gamma','log Activity'
write(ioutput,9) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '-------------------------------'
do j=1,this%numaqprisp
 namecomp = this%pspecies(this%idaqprisp(j))%ptr%name 
 conc = u(j) * omgwfree
 tot = t(j) * omgwfree
 if (dabs(tot)>zero) then  
   write(ioutput,11) namecomp,tot,conc,u(j) 
   do i=1,this%numsp
     if (dabs(this%u(j,i))>zero) then 
       name=this%pspecies(i)%ptr%name
       if (c(i)>zero) then
        conc=c(i)
        gamma=g(i)
        if (i/=ispsw) then
          act=conc*gamma
          mol=conc*omgwfree
        else
          mol=omgwfree/kgwmol
          act=gamma
        end if
        write(ioutput,2) name,mol,conc,gamma,act,dlog10(gamma),dlog10(act)
       end if
     end if 
   end do
 end if 
end do 
!%------------------------------------------------------------------
!% Compute saturation indices (optional)
!%------------------------------------------------------------------
call check_pointer_ (siloc,this%numsp,.true.)
if (havesimin.and.this%numminph>0) then
 isanomalous=.false.
 isupmxitergam=.false. 
 siloc=simin
else
 call check_pointer_ (gloc,this%numsp,.true.)
 siloc=c
 gloc=g
 do i=1,this%numminph
  call compute_secondaries_ &
     (this, &
      siloc, &
      gloc, &
      dc, &
      dg, &
      this%aqphindex, &
      this%idminph(i), &
      ionstrloc, &
      dionstrloc, &
	  r1, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      .false., &
      msg, &
      iserror)
  if (iserror.or.isanomalous.or.isupmxitergam) goto 20
 end do
 
end if
!%------------------------------------------------------------------
!% Write saturation indices 
!%------------------------------------------------------------------
if (this%numminph>0) then
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 write(ioutput,*) '             Saturation Indices         '
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 write(ioutput,*) 'Mineral species      IAP/K    log(IAP/K)      '
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 do i=1,this%numminph
  call get_iposspsph_ (this,this%idminph(i),isps1,isps2)
  do j=isps1,isps2
   name=this%pspecies(j)%ptr%name
   conc=siloc(j)
   write(ioutput,8) name,conc,dlog10(conc)
  end do
 end do
 
end if
!%------------------------------------------------------------------
!% Write infomation about surface complexes
!%------------------------------------------------------------------
if (this%numsurf>0) then
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 write(ioutput,*) '             Surface complexes          '
 do i=1,this%numsurf
  name='Surface'
  call add_ (name,i)
  call get_iposspsurf_(this,i,isps1,isps2)
  write(ioutput,9) '----------------------------------------'// &
                   '----------------------------------------'// &
                   '-------------------------------'
  write(ioutput,*) name
  write(ioutput,9) '----------------------------------------'// &
                   '----------------------------------------'// &
                   '-------------------------------'
 write(ioutput,*) 'Species                  Molality              '
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
  do j=isps1,isps2
   name=this%pspecies(j)%ptr%name
   conc=c(j)
   write(ioutput,7) name,conc
  end do
 end do 
end if 
!%------------------------------------------------------------------
!% Optional: Write kinetic information
!%------------------------------------------------------------------
if (iswritekin) then
!%-------------------------------------------------------
!% Consumption/production of concentrations due
!% kinetic reactions
!%-------------------------------------------------------
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 write(ioutput,*) ' Consumption/production of concentrations'// &
                  ' due kinetic reactions '
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 write(ioutput,*) 'species          dc/dt [mol/s]          '
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 do i=1,this%numsp
   name=this%pspecies(i)%ptr%name
   conc=sktrk(i)*omgwfree
   write(ioutput,7) name,conc
 end do
!%-------------------------------------------------------
!% Consumption/production of components due kinetic 
!% reactions
!%-------------------------------------------------------
 call get_iposspsph_ (this,this%aqphindex,isps1,isps2)
 u=matmul(this%ueq(:,isps1:isps2),sktrk(isps1:isps2))
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 write(ioutput,*) ' Consumption/production of components'// &
                  ' due kinetic reactions '
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 write(ioutput,*) 'Component          du/dt [mol/s]        '
 write(ioutput,9) '----------------------------------------'// &
                  '----------------------------------------'// &
                  '-------------------------------'
 do i=1,this%numaqprisp
  name=this%pspecies(this%idaqprisp(i))%ptr%name
  conc=u(i)*omgwfree
  write(ioutput,7) name,conc
 end do
 
end if
!%---------------------------------------------------------------
write(ioutput,9) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '-------------------------------'
20 continue
!%---------------------------------------------------------------
!% Deallocate local pointers
!%--------------------------------------------------------------- 
call check_pointer_ (u,1,.false.)
call check_pointer_ (t,1,.false.)
call check_pointer_ (siloc,1,.false.)
call check_pointer_ (gloc,1,.false.)
if (iserror.or.isanomalous.or.isupmxitergam) goto 10
!%---------------------------------------------------------------
return
2 format(5x,a17,4e15.4,5x,f10.4,5x,f10.4)
3 format(a32,f10.4)
4 format(a32,e10.4)
5 format(5x,a7,a15,a22,a16,a15,3x,a15,a15)
6 format(a5,a1,2e15.3)
7 format(a15,e15.3)
8 format(a13,e15.3,2x,f8.2)
9 format(a113)
11 format(a17,3e15.4)
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public Subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_idreaction_pchemsys &
   (this, &
    idreact, &
    nreact, &
    ithph, &
    jthph, &
    isbesurf)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the global reaction indices between phases 
! (ithph and jthph) or between phase and surface (isbesurf=.true.). 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)         :: this

integer, intent(in)                               :: ithph

integer, intent(in)                               :: jthph

integer, pointer, dimension(:)                    :: idreact 

integer, intent(out)                              :: nreact

logical, intent(in)                               :: isbesurf 
 
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
 index, &
 i
integer, pointer                      :: &
 idreactloc(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
index=0
nreact=0
!%-----------------------------------------------------------
call check_pointer_ (idreactloc,this%numreact,.true.)
!%-----------------------------------------------------------
if (isbesurf) index=this%numph
!%-----------------------------------------------------------
if (jthph==0) then
!%---------------
 do i=1,this%numreact
  if (this%ishomreact(i,ithph)) then
          nreact=nreact+1
		  idreactloc(nreact)=i
  end if
 end do

else
!%--------------- 
 do i=1,this%numreact
  if (this%ishetreact(i,ithph) &
               .and. &
      this%ishetreact(i,index+jthph)) then
          nreact=nreact+1
		  idreactloc(nreact)=i
  end if
 end do

end if
!%-----------------------------------------------------------
if (nreact>0) then
 call check_pointer_ (idreact,nreact,.true.)
 idreact=idreactloc(1:nreact)
end if
!%-----------------------------------------------------------
!% Deallocate local pointers   
!%-----------------------------------------------------------
call check_pointer_ (idreactloc,1,.false.)
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public Subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_min_sp_index_pchemsys &
   (this, &
    iserror, &
	nameminsp, &
    idminsp, &
    nminsp)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return global index corresponding to mineral species 
! and the total number of mineral species. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)             :: this      ! Type parent chemical system variable. 

integer, pointer, optional, dimension(:)              :: idminsp   ! Global index of the mineral species. 

character(len=100), pointer, optional, dimension(:)   :: nameminsp ! Name of the mineral species. 

integer, intent(out), optional                        :: nminsp    ! Number of mineral species. 

logical, intent(out)                                  :: iserror   ! iserror=true, there was an error. 
 
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
 ipos1, &
 ipos2, &
 isploc, &
 isp, &
 i
logical                               :: &
 havenminsp, &
 haveidminsp, &
 havenameminsp
integer, pointer                      :: &
 ivector(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
!%-----------------------------------------------------------
!% Check optional arguments
!%-----------------------------------------------------------
havenminsp=present(nminsp)
haveidminsp=present(idminsp)
havenameminsp=present(nameminsp)
!%-----------------------------------------------------------
if (havenminsp) nminsp=0 
if (this%numminph==0) return
if (.not.haveidminsp.and..not.havenminsp.and..not.havenameminsp) return
!%-----------------------------------------------------------
 call check_pointer_ (ivector,this%numsp,.true.)
 isploc=0
!%-----------------------------------------------------------
 do i=1,this%numminph
  call get_iposspsph_ (this,this%idminph(i),ipos1,ipos2)
  do isp=ipos1,ipos2
   isploc=isploc+1
   if (haveidminsp.or.havenameminsp) ivector(isploc)=isp
  end do
 end do
!%-----------------------------------------------------------
if (haveidminsp) then
 call check_pointer_ (idminsp,isploc,.true.)
 idminsp=ivector(1:isploc)
end if
!%-----------------------------------------------------------
if (havenameminsp) then
 call check_pointer_ (nameminsp,isploc,.true.)
 do i=1,isploc
  nameminsp(i)=this%pspecies(ivector(i))%ptr%name
 end do
end if
!%-----------------------------------------------------------
if (havenminsp) nminsp=isploc
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (ivector,1,.false.)
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public Subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_idreaction_from_sps_pchemsys &
   (this, &
    idreact, &
    nreact, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the global index of reactions from global index 
! of secundary species 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)    :: this

integer, intent(inout)                       :: nreact

integer, intent(inout), dimension(nreact)    :: idreact

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
integer                :: &
 id(this%numreact), &
 idreact1(nreact)
integer                :: &
 isps, &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
!%-----------------------------------------------------------
idreact1=idreact
!%-----------------------------------------------------------
idreact=0 
!%-----------------------------------------------------------
do i=1,nreact
 isps=idreact1(i) 
 id=0 
 where (this%idreactsp==isps)
  id=1
 end where 
 idreact(i:i)=maxloc(id) 
end do 
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public Subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_gas_sp_index_pchemsys &
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
 
type (t_parentchemicalsystem), intent(in)         :: this

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
integer                               :: &
 ipos1, &
 ipos2, &
 isploc, &
 isp, &
 i
logical                               :: &
 havengassp, &
 haveidgassp
integer, pointer                      :: &
 ivector(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!
 
 


!%-----------------------------------------------------------
iserror=.false.
!%-----------------------------------------------------------
ngassp=0
!%-----------------------------------------------------------
!% Check optional arguments
!%-----------------------------------------------------------
havengassp=present(ngassp)
haveidgassp=present(idgassp)
!%-----------------------------------------------------------
if (this%numgasph==0) return
if (.not.haveidgassp.and..not.havengassp) return
!%-----------------------------------------------------------
 call check_pointer_ (ivector,this%numsp,.true.)
 isploc=0
!%-----------------------------------------------------------
 do i=1,this%numgasph
  call get_iposspsph_ (this,this%idgasph(i),ipos1,ipos2)
  do isp=ipos1,ipos2
   isploc=isploc+1
   if (haveidgassp) ivector(isploc)=isp
  end do
 end do
!%-----------------------------------------------------------
if (haveidgassp) then
 call check_pointer_ (idgassp,isploc,.true.)
 idgassp=ivector(1:isploc)
end if
!%-----------------------------------------------------------
if (havengassp) ngassp=isploc
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (ivector,1,.false.)
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_mineral_pchemsys &
   (this, &
    c, &
    idminglob, &
    nminglob, &
    idmin, &
    nmin, &
    isbenegmin, &
    isbenewmin, &
    si)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Check the minerals present into the chemical system if
! the saturation indexes is present as argument,
! check the saturated minerals and build the new
! equilibrium mineral set.
!
!   $Arguments:
!
 
 
type (t_parentchemicalsystem), intent(in)            :: this

real*8, intent(inout), dimension(this%numsp)         :: c

integer, intent(in)                                  :: nminglob    !Global number of minerals

integer, intent(out)                                 :: nmin

integer, intent(in), dimension(nminglob)             :: idminglob

integer, pointer, dimension(:)                       :: idmin 

logical, intent(out)                                 :: isbenegmin  !Indicates if there are negative minerals

logical, intent(out)                                 :: isbenewmin  !Indicates if there are new minerals

real*8, intent(in), dimension(this%numsp), optional  :: si
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                     :: &
 i, &
 nnewmin, &
 j, &
 imaxval, &
 nnegmin, &
 nposmin, &
 imin
logical                     :: &
 isbemin, &
 isbesatmin, &
 havesi
real*8                      :: &
 mxsat
real*8, parameter           :: &
 r0=0.0d0, &
 r1=1.0d0
integer, pointer            :: &
 idnegmin(:) => null (), &
 idposmin(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%---------------------------------------------------------------
nnegmin=0
nposmin=0
nmin=0 
isbenegmin=.false.
isbenewmin=.false.
!%---------------------------------------------------------------
!%If there are not mineral then return 
!%---------------------------------------------------------------
if (nminglob==0) return
!%---------------------------------------------------------------
!% Check optional arguments 
!%---------------------------------------------------------------
havesi=present(si)
!%---------------------------------------------------------------
call check_pointer_ (idnegmin,nminglob,.true.)
call check_pointer_ (idposmin,nminglob,.true.)
!%---------------------------------------------------------------
!% Check negative minerals and minerals present
!%---------------------------------------------------------------
do i=1,nminglob
     imin=idminglob(i)
     if (c(imin)<r0) then
        c(imin)=r0
        nnegmin=nnegmin+1
        idnegmin(nnegmin)=imin
     else if (c(imin)>r0) then
        nposmin=nposmin+1
        idposmin(nposmin)=imin
     end if
end do
!%---------------------------------------------------------------
!% Check saturated minerals
!%---------------------------------------------------------------
   nnewmin=0
   isbesatmin=.false.
   isbemin=.false.
!%---------------------------------------------------------------
!% Check saturated minerals if si argument is present
!%---------------------------------------------------------------
if (havesi) then
   mxsat=maxval(si(idminglob))
   if (mxsat>=r1+this%deltasatmin) then
    do i=1,nminglob
      if(si(idminglob(i))==mxsat) then
          imaxval=i
          isbesatmin=.true.
        exit
      end if
    end do
   end if
!%---------------------------------------------------------------
!% Check if the saturated mineral was considered in the set
!% and add a new set if not
!%---------------------------------------------------------------
   if (isbesatmin) then
         do j=1,nposmin
          if (idposmin(j)==idminglob(imaxval)) then
            isbemin=.true.
            exit
          end if
         end do
   end if
!%------------------------------
   if (isbesatmin.and..not.isbemin) then
         do j=1,nnegmin
          if (idnegmin(j)==idminglob(imaxval)) then
            isbemin=.true.
            exit
          end if
         end do
   end if
!%------------------------------   
   if (isbesatmin.and..not.isbemin) then
        nnewmin=nnewmin+1
        nposmin=nposmin+1
        idposmin(nposmin)=idminglob(imaxval)
   end if
!%------------------------------   
end if
!%---------------------------------------------------------------
!% Allocate new set of minerals
!%---------------------------------------------------------------
nmin=nposmin
if (nmin>0) then
call check_pointer_ (idmin,nmin,.true.)
idmin=idposmin(1:nmin)
end if
!%--------------------------------------------------------------
!% If there are negative mineral then isbenegmin=true 
!%--------------------------------------------------------------
if (nnegmin>0) isbenegmin=.true.
!%--------------------------------------------------------------
!% If there are new minerals then isbenewmin=true 
!%--------------------------------------------------------------
if (nnewmin>0) isbenewmin=.true.
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (idposmin,1,.false.)
call check_pointer_ (idnegmin,1,.false.)
!%--------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_gas_pchemsys &
   (this, &
    c, &
    isneggas)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Check the gas species are negative (the problem is
! when the equilibrium with gas phases are imposed 
!
!   $Arguments:
!
 
 
type (t_parentchemicalsystem), intent(in)            :: this ! Type chemical system variable. 

real*8, intent(inout), dimension(this%numsp)         :: c ! Concentration vector. 

logical, intent(out)                                 :: isneggas ! isneggas=true, the gas concentration is negative. 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                     :: &
 i, &
 iph, &
 isps, &
 isps1, &
 isps2
!-------------------------------------------------------------------------
!
!   $code
!
!%---------------------------------------------------------------
isneggas=.false.
!%---------------------------------------------------------------
do iph=1,this%numgasph
  call get_iposspsph_ (this,this%idgasph(iph),isps1,isps2)
  do isps=isps1,isps2
   if (c(isps)<0.0d0) then
    c(isps)=0.0d0 
    isneggas=.true.
   end if 
  end do 
end do 
!%--------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_aq_sps_pos_pchemsys &
   (this, &
    isps1, &
    isps2, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the first and last global index of aqueous species.
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)        :: this ! Type parent chemical system variable. 

integer, intent(out)                             :: isps1 ! Global index of the first aqueous species. 

integer, intent(out)                             :: isps2 ! Global index of the last aqueous species. 

logical, intent(out)                             :: iserror ! iserror=true, there was an error. 
 
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
isps1=0 
isps2=0 
!%------------------------------------------------------------
if (this%aqphindex==0) return
!%------------------------------------------------------------
call get_iposspsph_ (this,this%aqphindex,isps1,isps2)
!%------------------------------------------------------------
 
return
 
10 continue 
print *, '**************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: get_aq_sps_pos_'
print *, msg
print *, '**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_kinetic_pchemsys &
   (this, &
    ck1, &
    g, &
    dck1, &
    dg, &
    areak1, &
    ck, &
    dtime, &
    msg, &
    iserror, &
    sktrk,   & 
    dsktrk)    
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the kinetic calculations. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)                 :: this      ! Type parent chemical system variable. 

real*8, intent(inout), dimension(this%numsp)              :: ck1       ! Concentration vector in k+1.

real*8, intent(in), dimension(this%numsp)                 :: areak1    ! Area of minerals.

real*8, intent(in), dimension(this%numsp)                 :: ck        ! Concentrations vector in k

real*8, intent(in), dimension(this%numsp)                 :: g         ! Activity coefficients vector 

real*8, intent(in), dimension(this%numsp,this%numaqprisp) :: dck1      ! Derivatives of the concentrations. 

real*8, intent(in), dimension(this%numsp,this%numaqprisp) :: dg        ! Derivatives of the activity coefficients. 

real*8, intent(in)                                        :: dtime     ! Time increment [s] 

logical, intent(out)                                      :: iserror   ! iserror=true, then there was an error 

character(len=*), intent(out)                             :: msg       ! Message error. 

real*8, pointer, dimension(:), optional                   :: sktrk     ! Changes in species due kinetic reactions [numsp]

real*8, pointer, dimension(:,:), optional                 :: dsktrk    ! Derivatives of the changes in the species due kinetic reactions. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License: CSIC-UPC
!
!-------------------------------------------------------------------------
integer                              :: &
 i, &
 j, &
 isps, &
 iaqsps1, &
 iaqsps2, &
 naqsps, &
 npri, &
 ireact, &
 nkinreact, &
 ispsw
real*8                               :: &
 rk, &
 cw, &                   
 cold, &
 area 
real*8, pointer                      :: &
 drk(:)  => null (), &
 stqi(:) => null (), &
 uaq(:) => null (), &
 duaq(:,:) => null ()
character(len=100), pointer          :: &
 namesp(:) => null ()
logical                              :: &
 havesktrk, &
 havedsktrk, &
 isbe 
type(t_reaction), pointer            :: &
 preaction => null ()
real*8, parameter                    :: &
 r0=0.0d0, &
 r1=1.0d0
!-------------------------------------------------------------------------
!
!   $code
!
!%---------------------------------------------------------
msg=''
iserror=.false.
!%---------------------------------------------------------
!% Check optional arguments
!%---------------------------------------------------------
havesktrk=present(sktrk)
havedsktrk=present(dsktrk)
!%---------------------------------------------------------
!% If dsktrk array is present then allocate and zeoring 
!%---------------------------------------------------------
if (havesktrk) then
  call check_pointer_ (sktrk,this%numsp,.true.)
end if
!%---------------------------------------------------------
!% If sktrk vector is present then allocate and zeoring 
!%---------------------------------------------------------
if (havedsktrk) then
  call check_pointer_ (dsktrk,this%numsp,this%numaqprisp,.true.)
end if
!%---------------------------------------------------------
!% Get information (the name of the species) 
!%---------------------------------------------------------
call get_chem_info_ (this,msg,iserror,namesp=namesp)
if (iserror) goto 20 
!%---------------------------------------------------------
!% Compute u y du/dc1
!%---------------------------------------------------------
call get_aq_sps_pos_ (this,iaqsps1,iaqsps2,iserror)
call get_numsp_ (this%pphase(this%aqphindex)%ptr,naqsps)
if (iserror) goto 20 
call check_pointer_ (uaq,this%numaqprisp,.true.)
uaq=matmul(this%ueq(:,iaqsps1:iaqsps2),ck1(iaqsps1:iaqsps2))
if (havedsktrk) then
 call check_pointer_ (duaq,this%numaqprisp,this%numaqprisp,.true.)
 duaq=matmul(this%ueq(:,iaqsps1:iaqsps2),dck1(iaqsps1:iaqsps2,:))
end if
!%---------------------------------------------------------
!% Get the index of water species 
!%---------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wspindex=ispsw)
if (ispsw>0) then
 cw=ck1(ispsw)
 ck1(ispsw)=r1
end if
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------  
do i=1,this%numreact
 if (this%iskinreact(i)) then   
   preaction => this%preaction(i)%ptr
   cold=r0
   area=r0
   isps=this%idreactsp(i)
   stqi => this%stq(i,:) 
   if (isps>0) then 
    cold = ck(isps)
	area = areak1(isps)
   end if    
!%------------------------------------------------------------
!% Compute reaction rate
!%------------------------------------------------------------
   call compute_rk_ &
   (preaction, &
    rk, &
	ck1(this%idaqprisp), &
    g(this%idaqprisp), &
    ck1(iaqsps1:iaqsps2), &
    g(iaqsps1:iaqsps2), &
    this%numaqprisp, &
    naqsps, &
    namesp(iaqsps1:iaqsps2), &
    cold, &
    area, &
    iserror)
    if (iserror) goto 20       
!%-----------------------------------------------------------
!% If the mineral concentration is negative, then 
!% rk = ck/dtime
!%-----------------------------------------------------------
    if (isps>0) then 
        ck1(isps)=ck(isps) + stqi(isps) * rk * dtime
        if (ck1(isps)<r0) then
            ck1(isps)=r0
	        rk=ck(isps)/dtime  
	        isbe=.false.
        else if (ck1(isps)==r0) then
            isbe=.false. 
	    else
	        isbe=.true. 
        end if 
	end if 
!%-----------------------------------------------------------
!% Optional: compute sktrk
!%-----------------------------------------------------------
    if (havesktrk) then
      sktrk = sktrk + stqi * rk
    end if
!%-----------------------------------------------------------
!% Optional Compute dsktrk/dc1
!%-----------------------------------------------------------
    if (havedsktrk) then
         call compute_drk_ &
         (preaction, &
          drk, &
          ck1(this%idaqprisp), &
          g(this%idaqprisp), &
          dck1(this%idaqprisp,:), &
          dg(this%idaqprisp,:), &
          ck1(iaqsps1:iaqsps2), &
          dck1(iaqsps1:iaqsps2,:), &
          g(iaqsps1:iaqsps2), &
          dg(iaqsps1:iaqsps2,:), &
          namesp(iaqsps1:iaqsps2), &
          this%numaqprisp, &
          naqsps, &
          this%numaqprisp, &
          area, &
          isbe, &
          iserror)
          if (iserror) goto 20
          do j=1,this%numsp
              dsktrk(j,:) = dsktrk(j,:) + stqi(j) * drk
          end do
    end if 
 end if
end do
!%------------------------------------------------------------
!% Zeroing the derivatives with respect to water
!%------------------------------------------------------------
if (havedsktrk.and.this%wcindex>0) then
  dsktrk(1:this%numsp,this%wcindex)=r0
end if
!%------------------------------------------------------------
!% Restart the water concentration in the vector 
!%------------------------------------------------------------
if (ispsw>0) ck1(ispsw)=cw
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (drk,1,.false.)
call check_pointer_ (uaq,1,.false.)
call check_pointer_ (duaq,1,1,.false.)
stqi => null ()
preaction => null ()
if (iserror) goto 10
!%------------------------------------------------------------
return
 
10 iserror=.true.
return
end subroutine
!%************************************************************
!%****************Public subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dusktrk1_pchemsys &
   (this, &
    dusktrk, &
    npri, &
    c, &
    alpha, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute derivative of changes in components due kinetic 
! reactions with respect to primary aqueous species. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)        :: this

integer, intent(in)                              :: nsp

integer, intent(out)                             :: npri

real*8, intent(in), dimension(nsp)               :: alpha

real*8, intent(in), dimension(nsp)               :: c

real*8, pointer, dimension(:,:)                  :: dusktrk

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
integer                              :: &
 i, &
 j, &
 isps
real*8                               :: &
 ionstr
real*8, pointer                      :: &
 dsktrk(:,:) => null (), &
 drk(:) => null (), &
 g(:) => null (), &
 dc(:,:) => null (), &
 dg(:,:) => null (), &
 dionstr(:) => null (), &
 stqi(:) => null (), &
 cloc(:) => null ()
character(len=100), pointer          :: &
 namesp(:) => null ()
logical                              :: &
 isanomalous, &
 isupmxitergam
character(len=100)                   :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%----------------------------------------------------------------
npri=this%numaqprisp
!%----------------------------------------------------------------
!% Allocate local pointers 
!%----------------------------------------------------------------
call check_pointer_ (dsktrk,this%numsp,npri,.true.)
call check_pointer_ (g,this%numsp,.true.)
call check_pointer_ (cloc,this%numsp,.true.)
call check_pointer_ (dg,this%numsp,npri,.true.)
call check_pointer_ (dc,this%numsp,npri,.true.)
cloc=c
!%----------------------------------------------------------------
!% Allocate dusktrk 
!%----------------------------------------------------------------
call check_pointer_ (dusktrk,npri,npri,.true.)
!%----------------------------------------------------------------
!
!%----------------------------------------------------------------
do i=1,this%numaqprisp
 dc(this%idaqprisp(i),i)=1.0d0
end do
!%----------------------------------------------------------------
!% Compute the aqueous speciation 
!%----------------------------------------------------------------
call compute_secondaries_ &
     (this, &
      cloc, &
      g, &
      dc, &
      dg, &
      this%aqphindex, &
      0, &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .true., &
      isanomalous, &
      isupmxitergam, &
      .true., &
      msg, &
      iserror)
 
if (iserror.or.isupmxitergam) goto 20
!%----------------------------------------------------
call get_chem_info_ (this,msg,iserror,namesp=namesp)
!%----------------------------------------------------
do i=1,this%numreact
 if (this%iskinreact(i))  then
  isps=this%idreactsp(i)
  stqi => this%stq(i,:)
  call compute_drk_ &
    (this%preaction(i)%ptr, &
     drk, &
     cloc(this%idaqprisp), &
     g(this%idaqprisp), &
     dc(this%idaqprisp,:), &
     dg(this%idaqprisp,:), &
     cloc, &
     dc, &
     g, &
     dg, &
     namesp, &
     npri, &
     nsp, &
     npri, &
     alpha(isps), &
     .true., &
     iserror)
  if (iserror) goto 20
  do j=1,this%numsp
    dsktrk(j,:) = dsktrk (j,:) + stqi(j) * drk
  end do
 end if
end do
!%------------------------------------------------------------
!% Multiply by equilibrium components matrix 
!%------------------------------------------------------------
dusktrk=matmul(this%ueq,dsktrk)
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dg,1,1,.false.)
call check_pointer_ (drk,1,.false.)
call check_pointer_ (dsktrk,1,1,.false.)
call check_pointer_ (g,1,.false.)
call check_pointer_ (cloc,1,.false.)
stqi => null ()
if (iserror) goto 10
!%------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_dusktrk_'
print *, msg
print *,'**************************'
iserror=.true.
return
 
 
end subroutine
!%************************************************************
!%****************Public subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dusktrk2_pchemsys &
   (this, &
    dusktrk, &
    npri, &
    dsktrk, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)        :: this

integer, intent(in)                              :: nsp

integer, intent(in)                              :: npri

real*8, intent(in), dimension(nsp,npri)          :: dsktrk

real*8, pointer, dimension(:,:)                  :: dusktrk

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

!%----------------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%----------------------------------------------------------------
if (npri/=this%numaqprisp) then
 msg='Error in number of primary species'
 goto 10
end if
!%----------------------------------------------------------------
call check_pointer_ (dusktrk,npri,npri,.true.)
!%------------------------------------------------------------
dusktrk=matmul(this%ueq,dsktrk)
!%------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: compute_dusktrk_'
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
subroutine get_iposspsph_pchemsys &
   (this, &
    ithph, &
    isps1, &
    isps2)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return first and last global indices of species of the 
! ith phase
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

integer, intent(in)                       :: ithph

integer, intent(out)                      :: isps1

integer, intent(out)                      :: isps2 
 
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
 i, &
 nsp 
!-------------------------------------------------------------------------
!
!   $code
!
isps1=0
isps2=0
!%------------------------------------------------------------
if (ithph>this%numph) return
!%------------------------------------------------------------
do i=1,ithph
 isps1=isps2+1 
 call get_numsp_ (this%pphase(i)%ptr,nsp)
 isps2=isps2+nsp
end do
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_iposspsurf_pchemsys &
   (this, &
    ithsurf, &
    isps1, &
    isps2)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return first and last global indices of species of the 
! ith surface
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in):: this

integer, intent(in)                      :: ithsurf

integer, intent(out)                     :: isps1

integer, intent(out)                     :: isps2 
 
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
 i, &
 numsp, &
 sum 
!-------------------------------------------------------------------------
!
!   $code
!
isps1=0
isps2=0
numsp=0
sum=0
!%------------------------------------------------------------
if (ithsurf>this%numsurf) return
!%------------------------------------------------------------
do i=1,this%numph
 call get_numsp_ (this%pphase(i)%ptr,numsp)
 sum=sum+numsp
end do
!%------------------------------------------------------------
isps1=sum
isps2=sum
numsp=0
!%------------------------------------------------------------
do i=1,ithsurf
 isps1=isps1+numsp
 call get_numsp_ (this%psurface(i)%ptr,numsp)
 isps2=isps2+numsp
end do
!%------------------------------------------------------------
isps1=isps1+1
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_molsalt_pchemsys &
   (this, &
    molsalt, &
    c)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)           :: this

real*8, intent(out)                                 :: molsalt

real*8, intent(in), dimension(this%numsp), target   :: c
 
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
 i, &
 isps1, &
 isps2
type(t_phase), pointer                :: &
 phase => null ()
real*8, pointer                       :: &
 cloc(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
!%------------------------------------------------------------
call get_iposspsph_(this,this%aqphindex,isps1,isps2)
!%------------------------------------------------------------
phase => this%pphase(this%aqphindex)%ptr
cloc => c(isps1:isps2)
!%------------------------------------------------------------
call compute_sum_c_ (phase,molsalt,cloc)
!%------------------------------------------------------------
!% Nullify local pointers 
!%------------------------------------------------------------
phase => null ()
cloc => null ()
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_mass_salt_pchemsys &
   (this, &
    mass, &
    mol, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the mass of salt in the aqueous phase
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this

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
!   $License:
!
!-------------------------------------------------------------------------
integer                       :: &
 i, &
 isps1, &
 isps2, &
 nsp1
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
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------
!% Get the positions of aqueous species
!%------------------------------------------------------------
call get_iposspsph_(this,this%aqphindex,isps1,isps2)
!%------------------------------------------------------------
call get_numsp_ (this%pphase(this%aqphindex)%ptr,nsp1)
!%------------------------------------------------------------
call compute_mass_ &
   (this%pphase(this%aqphindex)%ptr, &
    mass, &
    mol(isps1:isps2), &
    nsp1, &
    iserror)
!%------------------------------------------------------------
return
10 continue 
print *, '**************************'
print *, 'Chemical System:'
print *, 'Name:',this%name
print *, 'Service: compute_mass_salt_'
print *,  msg
print *, '**************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dmassSalt_dc_pchemsys &
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
 
type (t_parentchemicalsystem), intent(in) :: this

real*8,dimension(:,:),pointer             :: dc

real*8, dimension(:),pointer              :: dmassdc

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
integer                       :: &
 ipos1, &
 ipos2
 real*8,pointer, dimension(:,:)::aqdc
!-------------------------------------------------------------------------
!
!   $code
!

!get auqous index
call get_iposspsph_ (this, this%aqphindex, ipos1, ipos2)

!build dcaq/dpri
call check_pointer_ (aqdc,ipos2-ipos1+1,size(dc,2),.true.)
aqdc=dc(ipos1:ipos2,:)

!get dmass/dpri
call compute_dmass_dc_(this%pphase(this%aqphindex)%ptr, dmassdc, aqdc,iserror)
if (IsError) goto 10

!deallocate aqdc
call check_pointer_ (aqdc,ipos2-ipos1+1,size(dc,2),.false.)

return
10 continue 
print *, '**************************'
print *, 'Chemical System:'
print *, 'Name:',this%name
print *, 'Service: compute_dmassSalt_dc_'
print *, '**************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_new_chemsys_pchemsys &
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
!   $Description: Return new chemical system according name of phases and surfaces.  
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)     :: this

type (t_parentchemicalsystem), pointer        :: pnewchemsys

integer, intent(in)                           :: numph
 
integer, intent(in)                           :: numsurf

character(len=*), intent(in)                  :: nameph(numph)

character(len=*), intent(in)                  :: namesurf(numsurf)

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
type(t_phase), pointer                :: &
 phase(:) => null ()
type(t_surface), pointer           :: &
 surface(:) => null ()
type(t_reaction), pointer             :: &
 reaction(:) => null ()
integer, pointer                      :: &
 idph(:) => null (), &
 idsurf(:) => null ()
integer, pointer                      :: &
 idreaction(:) => null (), &
 idreactionloc(:) => null ()
integer                               :: &
 nreactloc, &
 i, &
 j, &
 ireact, &
 nsites
character(len=100)                    :: &
 name, &
 msg, &
 comp
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------
if (associated(pnewchemsys)) then
 nullify (pnewchemsys)
end if
!%--------------------------------------------------------
if (numph>0) then
 call check_pointer_ (idph,numph,.true.)
end if
if (numsurf>0) then
 call check_pointer_ (idsurf,numsurf,.true.)
end if
!%--------------------------------------------------------
do i=1,numph
 call get_ph_index_ (this,idph(i),nameph(i))
end do
!%-------------------------------------------------------
do i=1,numsurf
 call get_surf_index_ (this,idsurf(i),nsites, namesurf(i))
end do
!%-------------------------------------------------------
call check_pointer_ (idreactionloc,this%numreact,.true.)
ireact=0
!%-------------------------------------------------------
!%-----------------------------Find homogeneous reactions
!%-------------------------------------------------------
do i=1,numph
 call get_idreaction_ (this,idreaction,nreactloc,idph(i),0,.false.)
 if (nreactloc>0) then
  idreactionloc(ireact+1:ireact+nreactloc)=idreaction
  ireact=ireact+nreactloc
 end if
end do
!%-------------------------------------------------------
!% Find heterogeneous reactions between aqueous phase and
!% another phases
!%-------------------------------------------------------
 do j=1,numph
 
   if(idph(j).ne.this%aqphindex) then
 
     call get_idreaction_ (this,idreaction,nreactloc,this%aqphindex,idph(j),.false.)
 
   if (nreactloc>0) then
     idreactionloc(ireact+1:ireact+nreactloc)=idreaction
     ireact=ireact+nreactloc
   end if
 
 
   end if
 
  end do
!%-------------------------------------------------------
!% Find heterogeneous reactons between aqueous phase and 
!% surfaces
!%-------------------------------------------------------
do i=1,numsurf
 
  call get_idreaction_ (this,idreaction,nreactloc,this%aqphindex,idsurf(i),.true.)
  if (nreactloc>0) then
    idreactionloc(ireact+1:ireact+nreactloc)=idreaction
    ireact=ireact+nreactloc
  end if
end do
!%-------------------------------------------------------
nreactloc=ireact
!%-------------------------------------------------------
if (numph>0) allocate (phase(numph))
do i=1,numph
 call create_ (phase(i))
 phase(i) = this%pphase(idph(i))%ptr
end do
!%---------------------------------------------------------
if (numsurf>0) allocate (surface(numsurf))
do i=1,numsurf
 call create_ (surface(i))
 surface(i)=this%psurface(idsurf(i))%ptr 
end do
!%---------------------------------------------------------
if (nreactloc>0) allocate (reaction(nreactloc))
do i=1,nreactloc
  call create_ (reaction(i))
  reaction(i)=this%preaction(idreactionloc(i))%ptr
end do
!%-------------------------------------------------------
if (this%isgaussjordan) then
 comp='gj'
else
 comp='sv'
end if
!%------------------------------------------------------- 
name='Chemical System Derived'
!%-------------------------------------------------------
!% Allocate the new chemical system 
!%-------------------------------------------------------
allocate (pnewchemsys)
!%-------------------------------------------------------
!% Create the new chemical system 
!%-------------------------------------------------------
call create_ (pnewchemsys)
!%-------------------------------------------------------
!% Set the new chemical system 
!%-------------------------------------------------------
call set_ &
   (pnewchemsys, &
    name, &
    this%tempref, &
    phase, &
    surface, &
    reaction, &
    numph, &
    numsurf, &
    nreactloc, &
    comp, &
    this%tolunknr, &
    this%tolresnr, &
    this%deltasatmin, &
    iserror) 
if (iserror) then
 msg='Error when calling set_'
 goto 20
end if
!%--------------------------------------------------------
20 continue      
if (numph>0) then
 do i=1,numph
  call destroy_ (phase(i))
 end do
 deallocate (phase)
end if
!%--------------------------------------------------------
if (nreactloc>0) then
 do i=1,nreactloc
  call destroy_ (reaction(i))
 end do
 deallocate (reaction)
end if
!%--------------------------------------------------------
if (numsurf>0) then
 do i=1,numsurf
  call destroy_ (surface(i))
 end do
 deallocate (surface)
end if
!%--------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------
call check_pointer_ (idph,1,.false.)
call check_pointer_ (idreactionloc,1,.false.)
call check_pointer_ (idreaction,1,.false.)
call check_pointer_ (idsurf,1,.false.)
if (iserror) goto 10
!%------------------------------------------------------------
return
 
10 continue 
print *, '*************************'
print *, 'Chemical System:'
print *, 'Name:',this%name
print *, 'Service: get_new_chemsys_'
print *,  msg
print *, '*************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_convergence_pchemsys &
   (this, &
    isconvergence, &
    isdivergence, &
    isupmxiter, &
    cunk, &
    dtcunk, &
    residual, &
    nres, &
    nunk, &
    mxerrorunk, &
    iter, &
    idiverg, &
    tolunknr, &
    tolresnr)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Check convergence computing relative error and error in 
! residual for several speciation algorithm performed using newton raphson 
! method
!
!   $Arguments:
!
 
type(t_parentchemicalsystem), intent(in), target  :: this

integer, intent(in)                               :: nunk

integer, intent(in)                               :: iter

integer, intent(in)                               :: nres

real*8, intent(in), dimension(nunk)               :: cunk

real*8, intent(in), dimension(nres)               :: residual

integer, intent(inout)                            :: idiverg

real*8, intent(inout)                             :: mxerrorunk

real*8, intent(in)                                :: tolunknr

real*8, intent(in)                                :: tolresnr

logical, intent(out)                              :: isconvergence

logical, intent(out)                              :: isdivergence

logical, intent(out)                              :: isupmxiter

real*8, intent(inout), dimension(nunk)            :: dtcunk 
 
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
 error(:) => null ()
real*8                                    :: &
 mxerrorunkold, &
 mxerrorres
integer                                   :: &
 i, &
 ipri 
character(len=100), pointer               :: &
 name => null () 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (error,nunk,.true.)
!%-----------------------------------------------------------
isconvergence=.false.
isdivergence=.false.
isupmxiter=.false.
mxerrorunkold=mxerrorunk
!%-----------------------------------------------------------
!% Compute maximum error in residual
!%-----------------------------------------------------------
mxerrorres=maxval(dabs(residual))
!%-----------------------------------------------------------
!% Compute maximum error in unknowns 
!%-----------------------------------------------------------
where(dabs(cunk)>zero)
   error = dtcunk/cunk
end where
mxerrorunk = maxval (dabs(error))
!%----------------------------------------------------------
!% If the maximum error>facmax, then scale the increment in 
!% unknowns 
!%----------------------------------------------------------
if (mxerrorunk>facmax) then
  dtcunk = dtcunk * facmax/mxerrorunk
end if
!%-----------------------------------------------------------
!% Check convergence
!%-----------------------------------------------------------
isconvergence=(mxerrorunk<=tolunknr.and.mxerrorres<=tolresnr)
!%----------------------------------------------------------
!% Check if maximum number of iterations was exceed
!%----------------------------------------------------------
isupmxiter=(iter>mxiter)
!%----------------------------------------------------------
!% Check if there is divergence
!%----------------------------------------------------------
if (mxerrorunk>mxerrorunkold) then
  idiverg=idiverg+1
  isdivergence=(idiverg>mxdiverg)
end if
!%-----------------------------------------------------------
!% Write convergence information (optional)
!%-----------------------------------------------------------
if (this%iswriteinfo) then
 write (iouwiinfo,*) '-----------------------------------------'
 write (iouwiinfo,*) '     Relative error                      '
 write (iouwiinfo,5) 'unknown tol.=',tolunknr
 write (iouwiinfo,5) 'residual tol.=',tolresnr
 write (iouwiinfo,20)'-----------------------------------------'
 write (iouwiinfo,20)'          Delta unk.   unk.     error    '
 write (iouwiinfo,20)'-----------------------------------------'
 do i=1,this%numaqprisp
   ipri = this%idaqprisp(i)
   name => this%pspecies(ipri)%ptr%name
   write (iouwiinfo,10) name,'=',dtcunk(i),cunk(i),error(i)
 end do
 if (isconvergence) then
   write (iouwiinfo,*) '  CONVERGENCE ARRIVED!!!'
 end if
 write (iouwiinfo,*) '-------------------------------------'
end if
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (error,1,.false.)
name => null ()
!%-----------------------------------------------------------
return
5 format(a15,e10.3)
10 format(a8,a2,3e10.3)
20 format(a41)
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_ph_index_pchemsys &
   (this, &
    ithph, &
    nameph)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the global index of phase from name of the phase 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in):: this

integer, intent(out)                     :: ithph

character(len=*), intent(in)             :: nameph 
 
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
 namephloc
integer                                 :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
ithph=0
do i=1,this%numph
 
 call get_name_ (this%pphase(i)%ptr,namephloc)
 if (namephloc==nameph) then
  ithph=i
  exit  
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
subroutine begin_element_handler (name,attributes)
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
 
 
 
 call read_xml_loc_ (name,attributes)
 
 
 
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc_pchemsys (name,attributes,this,iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
character(len=*), intent(in)                          :: name

type(dictionary_t), intent(in)                        :: attributes

type(t_parentchemicalsystem), intent(inout), optional :: this

logical, intent(out), optional                        :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: Read the chemical system according xml type file
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer        :: status, n
real*8         :: reallocal(1)
integer        :: integerlocal(1)
 
character(len=100), save:: &
 namechemsys, &
 components
character(len=100), save:: &
 pitzerdatabase
character(len=100), pointer, save:: &
 nameph(:) => null (), &
 isconvph(:) => null (), &
 namereact(:) => null (), &
 typereact(:) => null (), &
 namesp(:) => null (), &
 namespreact(:,:) => null (), &
 nameprop(:,:) => null (), &
 unitprop(:,:) => null (), &
 namepropph(:,:) => null (), &
 modelpropph(:,:) => null (), &
 typeph(:) => null (), &
 actmodel(:) => null (), &
 namecatsp(:,:,:) => null (), &
 namerrlaw(:) => null (), &
 namepropsite(:,:,:) => null (), &
 nameintf(:) => null (), &
 modelintf(:) => null (), &
 typetermrrl(:,:) => null ()
real*8, save            :: &
 temp, &
 tolunknrloc, &
 tolresnrloc, &
 deltasatminloc 
real*8, pointer, save   :: &
 valueprop(:,:) => null (), &
 valuepropph(:,:) => null (), &
 coeffreact(:,:) => null (), &
 logk(:,:) => null (), &
 templogk(:,:) => null (), &
 ea(:) => null (), &
 attrrrlaw1(:) => null (), &
 attrrrlaw2(:) => null (), &
 attrrrlaw3(:) => null (), &
 attrrrlaw4(:) => null (), &
 valuepropsite(:,:,:) => null (), &
 attrtermrrl(:,:,:) => null (), &
 attrsprrl(:,:,:) => null ()
logical, save           :: &
 beph, &
 isbesurf, &
 bereact, &
 bekinetic, &
 bepitzer, &
 besp, &
 iswriteinfo
integer, save           :: &
 iph, &
 isurf, &
 isite, &
 ispsite, &
 isp, &
 ispglob, &
 ireact, &
 ilogk, &
 irrlaw, &
 iprop, &
 iterm
integer                 :: &
 ntotsp, &
 isp1, &
 isp2, &
 i, &
 j
integer, pointer, save           :: &
 numspph(:) => null (), &
 numpropph(:) => null (), &
 numspreact(:) => null (), &
 numprop(:) => null (), &
 numsite(:) => null (), &
 numspsite(:,:)=> null (), &
 numlogk(:) => null (), &
 numterm(:) => null (), &
 numspterm(:,:)=> null (), &
 idreactrrlaw(:) => null (), &
 numpropsite(:,:) => null (), &
 idprisite(:,:) => null (), &
 typerrlaw(:) => null (), &
 nattrtermrrl(:) => null ()
character(len=100)      :: &
 id, &
 msg
logical        :: &
 havethis, &
 betemp, &
 haveiserror
logical, pointer     :: &
 iskinreact(:)=> null ()
logical, pointer, save :: &
 isareadep(:) => null ()
type(t_species), allocatable:: &
 species(:) 
 type(t_phase), allocatable:: &
 phase(:) 
type(t_surface), allocatable:: &
 surface(:) 
type(t_reaction), allocatable:: &
 reaction(:) 
type(t_reactionratelaw), pointer:: &
 rrlaw => null ()
integer, parameter :: &
 mxdim=100 
!-------------------------------------------------------------------------
!
!   $code
!


!%-----------------------------------------------------------------
havethis = present(this)
haveiserror=present(iserror)
call lowercase (name)
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% If this is present like argument, then storage the information. 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
if (havethis) then
!%-----------------
  msg=''
  iserror=.false.
!%-----------------
  this%iswriteinfo=iswriteinfo
!%-----------------
  if (iph>0) then
   allocate(phase(iph))
  end if
!%-----------------
  if (isurf>0) then
   allocate (surface(isurf))
  end if
!%-----------------
  if (ireact>0) then
   allocate(reaction(ireact))
  end if
!%-----------------
  allocate(iskinreact(ireact))
  iskinreact=.false.
!%----------------------------------------
!% Create and set species
!%----------------------------------------
  allocate(species(ispglob))
  do i=1,ispglob
    call create_(species(i))
    call set_name_(species(i),namesp(i))
    call set_prop_(species(i),valueprop(i,1:numprop(i)), &
	               nameprop(i,1:numprop(i)),unitprop(i,1:numprop(i)),numprop(i))
  end do
!%----------------------------------------
isp1=0
isp2=0
!%-----------------------------------------------
!% Create and set phases
!%-----------------------------------------------
  do i=1,iph
       isp1=isp1+1
       isp2=isp2+numspph(i)
 
     call create_ (phase(i))
     
     call set_ &
     (phase(i), &
      species(isp1:isp2), &
      numspph(i), &
      numpropph(i), &
      nameph(i), &
      typeph(i), &
      actmodel(i), &
      namepropph(i,1:numpropph(i)), &
      modelpropph(i,1:numpropph(i)), &
      valuepropph(i,1:numpropph(i)), &
      isconvph(i), &
      iserror, &
      pitzerdatabase)
      
      if (iserror) then
         msg='Error when calling set_'
         goto 20
      end if
 
 
       isp1=isp2
 
  end do
!%-----------------------------------------------
!% Create and set surfaces
!%-----------------------------------------------
 do i=1,isurf
      isp1=isp1+1
      isp2=isp2+sum(numspsite(i,1:numsite(i)))
 
      call create_ (surface(i))

      call set_ &
      (surface(i), &
       nameintf(i), &
       modelintf(i), &
       species(isp1:isp2), &
       sum(numspsite(i,1:numsite(i))), &
       numsite(i), &
       numspsite(i,1:numsite(i)), &
       namepropsite(i,1,1:numpropsite(i,1)), &
       valuepropsite(i,1:numsite(i),1:numpropsite(i,1)), &
       numpropsite(i,1), &
       idprisite(i,1:numsite(i)), &
       iserror)
 
      if (iserror) then
        msg='Error when calling set_'
        goto 20
      end if
 
     isp1=isp2
 end do
 
!%-----------------------------------------------
!% Reaction and reaction rate law
!%-----------------------------------------------
 
  do i=1,ireact
 
   if (typereact(i)=='kinetic') then
     iskinreact(i)=.true.
   end if
 
   if(iskinreact(i)) then
        irrlaw=idreactrrlaw(i)
        if (irrlaw==0) then
          msg='Error, not defined reaction rate law for reaction'
        call add_ (msg,namereact(i))
        goto 20
      else
        ntotsp=sum(numspterm(irrlaw,1:numterm(irrlaw)))
		allocate (rrlaw)
        call create_ (rrlaw)
        call set_ & ! Set reaction rate law
        (rrlaw, &
         typerrlaw(irrlaw), &
         namerrlaw(irrlaw), &
         numspterm(irrlaw,1:numterm(irrlaw)), &
         attrtermrrl(irrlaw,1:nattrtermrrl(irrlaw),1:numterm(irrlaw)), &
         attrsprrl(irrlaw,1:ntotsp,1:numterm(irrlaw)), &
         namecatsp(irrlaw,1:ntotsp,1:numterm(irrlaw)), &
         typetermrrl(irrlaw,1:numterm(irrlaw)), &
         numterm(irrlaw), &
         nattrtermrrl(irrlaw), &
         ntotsp, &
         ea(irrlaw), &
         attrrrlaw1(irrlaw), &
         attrrrlaw2(irrlaw), &
         attrrrlaw3(irrlaw), &
		 attrrrlaw4(irrlaw), &
         isareadep(irrlaw), &
         iserror)
 
 
      end if
 
   end if
 
    call create_ (reaction(i))
 
    call set_ &
    (reaction(i), &
     namereact(i), &
     iskinreact(i), &
     coeffreact(i,1:numspreact(i)), &
     namespreact(i,1:numspreact(i)), &
     logk(i,1:numlogk(i)), &
     templogk(i,1:numlogk(i)), &
     numlogk(i), &
     5, &
     numspreact(i), &
     iserror, &
     rrlaw=rrlaw)
 
      if (iskinreact(i)) then
      call destroy_ (rrlaw)
      deallocate (rrlaw)
    end if
 
 
    if (iserror) goto 20
 
  end do
 
!%-----------------------------------------------
!% Set in chemical system
!%-----------------------------------------------
  call set_ &
   (this, &
    namechemsys, &
    temp, &
    phase, &
    surface, &
    reaction, &
    iph, &
    isurf, &
    ireact, &
    components, &
    tolunknrloc, &
    tolresnrloc, &
    deltasatminloc, &
    iserror)
  if (iserror) then
     msg='Error when calling set_'
     goto 20
  end if
!%--------------------------------------------------------
!% Deallocate pointer to species
!%--------------------------------------------------------
if (ispglob>0) then
 do i=1,iph
  call destroy_ (species(i))
 end do
 deallocate(species)
end if
!%--------------------------------------------------------
!% Deallocate to phase pointer
!%--------------------------------------------------------
if (iph>0) then
 do i=1,iph
  call destroy_ (phase(i))
 end do
 deallocate(phase)
end if
!%--------------------------------------------------------
!% Deallocate to reaction pointer
!%--------------------------------------------------------
if (ireact>0) then
 do i=1,ireact
  call destroy_ (reaction(i))
 end do
 deallocate(reaction)
end if
!%--------------------------------------------------------
!% Deallocate to surface pointer
!%--------------------------------------------------------
if (isurf>0) then
 do i=1,isurf
  call destroy_ (surface(i))
 end do
 deallocate(surface)
end if
!%--------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------
deallocate(iskinreact)
deallocate(nameph)
deallocate(isconvph)
deallocate(namereact)
deallocate(typereact)
deallocate(namesp)
deallocate(namespreact)
deallocate(nameprop)
deallocate(unitprop)
deallocate(modelpropph)
deallocate(namepropph)
deallocate(typeph)
deallocate(namecatsp)
deallocate(namerrlaw)
deallocate(namepropsite)
deallocate(nameintf)
deallocate(modelintf)
deallocate(valueprop)
deallocate(valuepropph)
deallocate(coeffreact)
deallocate(logk)
deallocate(templogk)
deallocate(ea)
deallocate(valuepropsite)
deallocate(numspph)
deallocate(numpropph)
deallocate(numspreact)
deallocate(numprop)
deallocate(numsite)
deallocate(numspsite)
deallocate(numlogk)
deallocate(numterm)
deallocate(numspterm)
deallocate(idreactrrlaw)
deallocate(numpropsite)
deallocate(idprisite)
deallocate(typerrlaw)
deallocate(nattrtermrrl)
deallocate(typetermrrl)
deallocate(attrtermrrl)
deallocate(attrsprrl)
deallocate(attrrrlaw1)
deallocate(attrrrlaw2)
deallocate(attrrrlaw3)
deallocate(attrrrlaw4)
deallocate(actmodel)
deallocate(isareadep)
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% If this is not present like argument, then read 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
else
 
 select case (name)
!%-------------------------------------
  case ('chemicalsystem')
 
allocate(nameph(mxdim))
allocate(isconvph(mxdim))
allocate(namereact(mxdim))
allocate(typereact(mxdim))
allocate(namesp(mxdim))
allocate(namespreact(mxdim,mxdim))
allocate(nameprop(mxdim,mxdim))
allocate(unitprop(mxdim,mxdim))
allocate(namepropph(mxdim,mxdim))
allocate(modelpropph(mxdim,mxdim))
allocate(typeph(mxdim))
allocate(actmodel(mxdim))
allocate(namecatsp(mxdim,mxdim,mxdim))
allocate(namerrlaw(mxdim))
allocate(namepropsite(mxdim,mxdim,mxdim))
allocate(nameintf(mxdim))
allocate(modelintf(mxdim))
allocate(valueprop(mxdim,mxdim))
allocate(valuepropph(mxdim,mxdim))
allocate(coeffreact(mxdim,mxdim))
allocate(logk(mxdim,mxdim))
allocate(templogk(mxdim,mxdim))
allocate(ea(mxdim))
allocate(valuepropsite(mxdim,mxdim,mxdim))
allocate(numspph(mxdim))
allocate(numpropph(mxdim))
allocate(numspreact(mxdim))
allocate(numprop(mxdim))
allocate(numsite(mxdim))
allocate(numspsite(mxdim,mxdim))
allocate(numlogk(mxdim))
allocate(numterm(mxdim))
allocate(numspterm(mxdim,mxdim))
allocate(idreactrrlaw(mxdim))
allocate(numpropsite(mxdim,mxdim))
allocate(idprisite(mxdim,mxdim))
allocate(typerrlaw(mxdim))
allocate(nattrtermrrl(mxdim))
allocate(typetermrrl(mxdim,mxdim))
allocate(attrtermrrl(mxdim,mxdim,mxdim))
allocate(attrsprrl(mxdim,mxdim,mxdim))
allocate(attrrrlaw1(mxdim))
allocate(attrrrlaw2(mxdim))
allocate(attrrrlaw3(mxdim))
allocate(attrrrlaw4(mxdim))
allocate(isareadep(mxdim))
        isareadep=.true. 
        attrtermrrl=0.0d0
        attrsprrl=0.0d0
		numspterm=0 
        namecatsp=''
        nameprop=' ' 
        unitprop=' ' 
        typetermrrl=''
        namepropph=' ' 
        modelpropph=' ' 
        numterm=0
        nattrtermrrl=0
        ea=0.0d0
		tolunknrloc=0.0d0
        tolresnrloc=0.0d0
        idreactrrlaw=0
        attrrrlaw1=0.0d0
        attrrrlaw2=0.0d0
		attrrrlaw3=0.0d0
        attrrrlaw4=0.0d0
        numpropph=0
        numprop=0
        numpropsite=0
        numsite=0
        pitzerdatabase=' '
        ireact=0
        ispglob=0
        isp=0
        iph=0
        isurf=0
        irrlaw=0
        ilogk=0
        temp=0.0d0
        id=' '
        components=''
        call get_value (attributes,"name", id, status)
        call lowercase (id)
        namechemsys=id
        id=' '
        call get_value (attributes,"temp", id, status)
        if (status<0) goto 10
        n=0
        reallocal(1)=0.0d0
        call build_data_array (id,reallocal,n)
        temp = reallocal(1)
        id='' 
        call get_value (attributes,"writeinfo", id, status)
        if (status<0) then
          iswriteinfo=.false.
        else
          select case (id)
          case ('yes','si','y')
           iswriteinfo=.true.
          case default
           iswriteinfo=.false.
          end select
        end if
        id=' '
        call get_value (attributes,"components", id, status)
        components=id
!%------------------------------------------------------------
!% Tolerancia in unknown newton raphson
!%------------------------------------------------------------
         id=' '
         reallocal(1)=0.0d0
         call get_value (attributes,"tolunknr", id, status)
         n=0
         call build_data_array (id,reallocal,n)
         tolunknrloc = reallocal(1)
!%------------------------------------------------------------
!% Tolerancia in residual newton raphson
!%------------------------------------------------------------
        id=' '
        reallocal(1)=0.0d0
        call get_value (attributes,"tolresnr", id, status)
        n=0
        call build_data_array (id,reallocal,n)
        tolresnrloc = reallocal(1)
!%------------------------------------------------------------
!% Tolerancia in residual newton raphson
!%------------------------------------------------------------
        id=' '
        reallocal(1)=0.0d0
        call get_value (attributes,"deltasatmin", id, status)
        n=0
        call build_data_array (id,reallocal,n)
        deltasatminloc = reallocal(1)
!%------------------------------------
 case ('phases')
!%------------------------------------
    iph=0
    beph=.true.
    isbesurf=.false.
    bereact=.false.
    besp=.false.
!%------------------------------------
 case ('phase')
!%------------------------------------
    iph=iph+1
    iprop=0
    isp=0
    besp=.false.
    beph=.true.
    isbesurf=.false.
    bepitzer=.false.
    bereact=.false.
    id=' '
    call get_value (attributes,"name", id, status)
    call lowercase (id)
    nameph(iph) = id
    id=' '
    call get_value (attributes,"type", id, status)
    if (status.lt.0) goto 10
    call lowercase (id)
    typeph(iph) = id
    id=' '
    call get_value (attributes,"model", id, status)
    call lowercase (id)
    actmodel(iph) = id
    if (actmodel(iph).eq.' ') actmodel(iph)='ideal'
    id=' '
    call get_value (attributes,"convention", id, status)
    call lowercase (id)
    isconvph(iph) = id
!%------------------------------------
 case ('surfaces')
 !%------------------------------------
    idprisite=0
    isurf=0
    beph=.false.
    isbesurf=.true.
    bereact=.false.
    besp=.false.
!%------------------------------------
 case ('surface')
!%------------------------------------
    isurf=isurf+1
    isite=0
    isp=0
    id=' '
    call get_value (attributes,"name", id, status)
    call lowercase (id)
    nameintf(isurf) = id
    id=' '
    call get_value (attributes,"model", id, status)
    if (status.lt.0) goto 10
    call lowercase (id)
    modelintf(isurf) = id
    id=' '
!%------------------------------------
 case ('site')
!%------------------------------------
    isite=isite+1
    numsite(isurf)=isite
    ispsite=0
    iprop=0
    isbesurf=.true.
    beph=.false.
    besp=.false.
    bereact=.false.
!%------------------------------------
 case ('filedatabase')
!%------------------------------------
    select case (actmodel(iph))
    case ('pitzer','pz')
     bepitzer=.true.
    end select
    if (bepitzer) then
     id=' '
     call get_value (attributes,"name", id, status)
     pitzerdatabase=id
    end if
!%------------------------------------
 case ('species')
!%------------------------------------
      besp=.true.
      iprop=0
      isp=isp+1
!%--------------
      if(beph) then
      ispglob=ispglob+1
      numspph(iph)=isp
      id=' '
      call get_value (attributes,"name", id, status)
      if (status/=0) goto 10
      call lowercase (id)
      namesp(ispglob) = id
!%--------------
    else if (isbesurf) then
      ispglob=ispglob+1
      ispsite=ispsite+1
      numspsite(isurf,isite)=ispsite
      id=' '
      call get_value (attributes,"name", id, status)
      if (status/=0) goto 10
      namesp(ispglob)=id
      call get_value (attributes,"type", id, status)
      if(status.ge.0.and.id.eq.'primary') then
      idprisite(isurf,isite)=isp
      end if
!%--------------
    else if (bereact) then    ! For reaction
      numspreact(ireact)=isp
      call get_value (attributes,"name", id, status)
      if (status/=0) goto 10
      call lowercase (id)
      namespreact(ireact,isp) = id
      call get_value (attributes,"coeff", id, status)
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
      coeffreact(ireact,isp) = reallocal(1)
   else if (bekinetic) then        ! For reaction rate law
      numspterm(irrlaw,iterm)=isp
      call get_value (attributes,"name", id, status)
      call lowercase (id)
      namecatsp(irrlaw,isp,iterm) = id
      call get_value (attributes,"attr1", id, status)
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
      attrsprrl(irrlaw,isp,iterm) = reallocal(1)
  end if
!%------------------------------------
 case ('property')
!%------------------------------------
     iprop=iprop+1
!%-------------
   if (besp) then
     numprop(ispglob)=iprop
     id=' '
     call get_value (attributes,"name", id, status)
     call lowercase (id)
     nameprop(ispglob,iprop) = id
     id=' '
     call get_value (attributes,"unit", id, status)
     call lowercase (id)
     unitprop(ispglob,iprop) = id
     id=' ' 
     call get_value (attributes,"value", id, status)
     n=0
     reallocal(1)=0.0d0
     call build_data_array (id,reallocal,n)
     valueprop(ispglob,iprop) = reallocal(1)
!%-------------
   else if (beph) then
     id=' ' 
     numpropph(iph)=iprop
     call get_value (attributes,"name", id, status)
     call lowercase (id)
     namepropph(iph,iprop) = id
     id=' ' 
     call get_value (attributes,"model", id, status)
     call lowercase (id)
     modelpropph(iph,iprop) = id
     id=' '
     call get_value (attributes,"value", id, status)
     n=0
     reallocal(1)=0.0d0
     call build_data_array (id,reallocal,n)
       valuepropph(iph,iprop) = reallocal(1)
!%-------------
   else if (isbesurf) then
     numpropsite(isurf,isite)=iprop
     id=' '
       call get_value (attributes,"name", id, status)
     namepropsite(isurf,isite,iprop)=id
     id=' '
       call get_value (attributes,"value", id, status)
     n=0
     reallocal(1)=0.0d0
     call build_data_array (id,reallocal,n)
       valuepropsite(isurf,isite,iprop) = reallocal(1)
   end if
!%------------------------------------
 case ('logk')
!%------------------------------------
    ilogk=ilogk+1
    numlogk(ireact)=ilogk
    call get_value (attributes,"temp", id, status)
    if (status/=0) goto 10
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    templogk(ireact,ilogk) = reallocal(1)
    call get_value (attributes,"value", id, status)
    if (status/=0) goto 10
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    logk(ireact,ilogk) = reallocal(1)
!%------------------------------------
 case ('reactions')
!%------------------------------------
    ireact=0
    irrlaw=0
    beph=.false.
    isbesurf=.false.
!%------------------------------------
 case ('reaction')
!%------------------------------------
    beph=.false.
	isbesurf=.false. 
	bekinetic=.false.
    bereact=.true.
    ireact=ireact+1
    ilogk=0
    isp=0
    id=' '
    call get_value (attributes,"name", id, status)
    call lowercase (id)
    namereact(ireact) = id
    call get_value (attributes,"type", id, status)
    call lowercase (id)
    typereact(ireact) = id
    if (typereact(ireact)=='kinetic') bekinetic=.true.
    isp=isp+1
    numspreact(ireact)=isp
    namespreact(ireact,isp) = namereact(ireact)
    coeffreact(ireact,isp) = -1.0d0
!%------------------------------------
 case ('reactionratelaw')
!%------------------------------------
  bereact=.false.
  if (bekinetic) then
    irrlaw=irrlaw+1
    idreactrrlaw(ireact)=irrlaw
    isp=0
    iterm=0
    id=''
    call get_value (attributes,"name", id, status)
    call lowercase (id)
    namerrlaw(irrlaw)=id
    id=''
    call get_value (attributes,"type", id, status)
    call lowercase (id)
!%-------------
    select case (id)
    case ('lasaga')
     typerrlaw(irrlaw)=1
     nattrtermrrl(irrlaw)=3
    case ('monod')
     typerrlaw(irrlaw)=2
     nattrtermrrl(irrlaw)=0
    case default
     typerrlaw(irrlaw)=0
     nattrtermrrl(irrlaw)=0
    end select
!%-------------
    call get_value (attributes,"ea", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    ea(irrlaw) = reallocal(1)
!%--------
    id=''
    call get_value (attributes,"attr1", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    attrrrlaw1(irrlaw) = reallocal(1)
!%--------
    id=''
    call get_value (attributes,"attr2", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    attrrrlaw2(irrlaw) = reallocal(1)
!%--------
    id=''
    call get_value (attributes,"attr3", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    attrrrlaw3(irrlaw) = reallocal(1)
!%--------
    id=''
    call get_value (attributes,"attr4", id, status)
    n=0
    reallocal(1)=0.0d0
    call build_data_array (id,reallocal,n)
    attrrrlaw4(irrlaw) = reallocal(1)
!%--------
    id=''
    call get_value (attributes,"areadep", id, status)
    select case (id)
    case ('NOT','not','n','N')
     isareadep(irrlaw)=.false. 
    end select 
  end if
!%------------------------------------
 case ('term')
!%------------------------------------
  if (bekinetic) then
      iterm=iterm+1
      numterm(irrlaw)=iterm
      isp=0
      if (typerrlaw(irrlaw)==1) then    ! For reaction rate law type
       typetermrrl(irrlaw,iterm)='catalyst'
       id=' '
       call get_value (attributes,"attr1", id, status)  ! rate constan
       if (status/=0) goto 10
       n=0
       reallocal(1)=0.0d0
       call build_data_array (id,reallocal,n)
       attrtermrrl(irrlaw,1,iterm) = reallocal(1)
       call get_value (attributes,"attr2", id, status)  ! theta
       n=0
       reallocal(1)=0.0d0
       call build_data_array (id,reallocal,n)
       attrtermrrl(irrlaw,2,iterm) = reallocal(1)
       call get_value (attributes,"attr3", id, status)   ! eta
       n=0
       reallocal(1)=0.0d0
       call build_data_array (id,reallocal,n)
       attrtermrrl(irrlaw,3,iterm) = reallocal(1)
      else if (typerrlaw(irrlaw)==2) then  ! For reaction rate law typ
       id=' '
       call get_value (attributes,"type", id, status)
       if (status/=0) goto 10
       typetermrrl(irrlaw,iterm)=id
      end if
  end if
!%------------------------------------
 case ('stoichiometry')
!%------------------------------------ 
 case default
 
   goto 10
 
 end select
 
 
 
end if
!%------------------------------------------------------------
return
 
10 continue 
print *,'*******************'
print *,'Chemical System:'
if (havethis) print *,'Name:',this%name
print *,'Service: read_xml_'
print *,'error in tag', name
print *,'*******************'
stop 
 
20 continue 
print *,'*******************'
print *,'Chemical System:'
if (havethis) print *,'Name:',this%name
print *,'Service: read_xml_'
print *,msg
iserror=.true.
print *,'*******************'

 
end subroutine

!%************************************************************
!%***************Private Subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_idreactionb_pchemsys &
   (this, &
    idreaction, &
    numreact, &
    iph1, &
    iph2, &
    isbesurf)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Give global reactions pointer of homogeneous and
!%   heterogeneous reactions
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in)         :: this

integer, intent(in)                               :: iph1

integer, intent(in)                               :: iph2

integer, pointer, dimension(:)                    :: idreaction 

integer, intent(out)                              :: numreact

logical, intent(in)                               :: isbesurf 
 
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
 i, &
 j, &
 isp, &
 isp1, &
 isp2, &
 tot, &
 ireact, &
 numspreact, &
 ndim
character(len=100)                   :: &
 namesp(this%numsp)
character(len=100), pointer          :: &
 namespreact(:) => null()
integer                              :: &
 idreactionlocal(this%numreact)
logical                              :: &
 beph1, &
 beph2 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%--------------------------------------------------------------
tot=0
numreact=0
ndim=1
if (iph2>0) ndim=2
 
select case (ndim)
 
case (1) ! Homogeneous reactions
 
 do i=1,this%numreact
 
      isp1=0
      
	  call get_namesp_ (this%preaction(i)%ptr, namespreact,numspreact)
      do j=1,numspreact
       call get_if_sp_is_present_ (this%pphase(iph1)%ptr,namespreact(j),beph1)
       if (beph1) then
        isp1=isp1+1
       end if
      end do
      if(isp1==numspreact) then
         numreact=numreact+1
         idreactionlocal(numreact)=i
      end if 
 end do
 
case (2) ! Heterogeneous reactions
 
   do i=1,this%numreact
      isp1=0
      isp2=0
       
      call get_namesp_ (this%preaction(i)%ptr, namespreact,numspreact)
 
      do j=1,numspreact
 
          call get_if_sp_is_present_ &
            (this%pphase(iph1)%ptr, &
             namespreact(j), &
             beph1)
 
        if (isbesurf) then
          call get_if_sp_is_present_ &
            (this%psurface(iph2)%ptr, &
             namespreact(j), &
             beph2)
            else
              call get_if_sp_is_present_ &
            (this%pphase(iph2)%ptr, &
             namespreact(j), &
             beph2)
          end if
 
          if (beph1) then
          isp1=isp1+1
        end if
        if (beph2) then
          isp2=isp2+1
        end if
 
        end do
 
        isp=isp1+isp2
        if(isp1<numspreact.and.isp2<numspreact &
         .and.isp==numspreact) then
        numreact=numreact+1
        idreactionlocal(numreact)=i
        end if
 
 
   end do
 
 
end select
!%-----------------------------------------------------------
if (numreact>0) then
  call check_pointer_ (idreaction,numreact,.true.)
  idreaction = idreactionlocal(1:numreact)
end if
!%-----------------------------------------------------------
call check_pointer_ (namespreact,1,.false.)
!%-----------------------------------------------------------
return
end subroutine

!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine build_components_matrix_gauss_jordan_pchemsys &
   (this, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build components matrix according gauss-jordan elimination
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout)   :: this

character(len=*), intent(out)                  :: msg

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
 
integer                         :: &
 i, &
 j
real*8                          :: &
 coeff
character(len=100), pointer     :: &
 namepri(:) => null()
logical                         :: &
 be 
!-------------------------------------------------------------------------
!
!   $code
!
 
 


!%---------------------------------------------------------
msg=''
iserror=.false.
!%---------------------------------------------------------
call check_pointer_ (namepri,this%numaqprisp,.true.)
call check_pointer_ (this%ueq,this%numaqprisp,this%numsp,.true.)
!%----------------------------------------------------------
call check_pointer_ (this%u,this%numaqprisp,this%numsp,.true.)
!%------------------------------------------------------------
!% Assign name of primary species
!%------------------------------------------------------------
do i=1,this%numaqprisp
 
 this%ueq(i,this%idaqprisp(i)) = 1.0d0
!%--------- 
 this%u(i,this%idaqprisp(i)) = 1.0d0
!%---------
 namepri(i)=this%pspecies(this%idaqprisp(i))%ptr%name
 
end do
!%---------------------------------------------------------
do i=1,this%numreact
 
   do j=1,this%numaqprisp
 
    call get_coeff_stq_ (this%preaction(i)%ptr,namepri(j),coeff)
 
    if (.not.this%iskinreact(i)) then
     this%ueq(j,this%idreactsp(i)) = coeff
	end if 
     
	this%u(j,this%idreactsp(i)) = coeff 
 
   end do
 
 
end do
!%------------------------------------------------------------
!% Eliminate of the column u the electron species
!%------------------------------------------------------------
if (this%ecindex.ne.0) then
 this%ueq(:,this%ecindex)=0.0d0
 this%u(:,this%ecindex)=0.0d0
end if 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (namepri,1,.false.)
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
subroutine build_stq_pchemsys &
   (this, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build stoichiometric matrix STQ 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout), target   :: this

character(len=*), intent(out)                          :: msg

logical, intent(out)                                   :: iserror 
 
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
 j 
real*8, pointer                 :: &
 coeff => null ()
character(len=100), pointer     :: &
 name => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
!%---------------------------------------------------------
msg=''
iserror=.false.
!%---------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------
call check_pointer_ (this%stq,this%numreact,this%numsp,.true.)
!%---------------------------------------------------------
!% Assign name of primary species
!%---------------------------------------------------------
do i=1,this%numreact
 do j=1,this%numsp
  coeff => this%stq(i,j)
  name => this%pspecies(j)%ptr%name
  call get_coeff_stq_ (this%preaction(i)%ptr,name,coeff)
 end do
end do
!%------------------------------------------------------------
!% Nullify local pointers 
!%------------------------------------------------------------
coeff => null ()
name => null ()
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
subroutine switch_base2_pchemsys &
   (this, &
    naqprisp, &
	nadsprisp, &
    nameaqprisp, &
	nameadsprisp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Switch the set of primary species of the chemical system 
! according its names. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(inout), target   :: this           ! Type parent chemical system variable

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
integer                         :: &
 i, &
 j, &
 ithntemp, &
 isps, &
 ipri, &
 isec, &
 nph, &
 nsurf, &
 nreact 
real*8                          :: &
 tempref, &
 tolunknr, &
 tolresnr, &
 deltasatmin 
real*8, pointer                 :: &
 s1(:,:) => null (), &
 s2(:,:) => null (), &
 logktemp(:) => null (), &
 templogk(:) => null (), &
 logk(:,:) => null (), &
 stq(:,:) => null (), &
 s2lu(:,:) => null (), &
 stqloc(:) => null ()
integer, pointer                :: &
 indx(:) => null (), &
 ntemp(:) => null (), &             ! Number of temperatures per reaction 
 idaqprisp(:) => null (), &
 idreactsp(:) => null ()
type(t_reaction), pointer         :: &
 reaction(:) => null ()
type(t_reactionratelaw), pointer         :: &
 rrlaw => null ()
type(t_species), pointer          :: &
 species => null ()
type(t_surface), pointer          :: &
 surface => null ()
character(len=100)                :: &
 msg, &
 name  
character(len=100), pointer       :: &
 namesp(:) => null ()
logical                           :: &
 ispri  
integer, parameter                :: &
 mxntemp=100, &    ! Maximum number of temperatures 
 ncoeff=5
!-------------------------------------------------------------------------
!
!   $code
!
 
!%---------------------------------------------------------
msg=''
iserror=.false.
!%---------------------------------------------------------
!% Check the number of primary species
!%---------------------------------------------------------
if (naqprisp/=this%numaqprisp) then
 msg='Error in number of primary species'
 goto 10
end if
!%---------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------
call check_pointer_ (s1,this%numreact,this%numaqprisp,.true.)
call check_pointer_ (s2,this%numreact,this%numreact,.true.)
call check_pointer_ (logk,this%numreact,mxntemp,.true.)
call check_pointer_ (s2lu,this%numreact,this%numreact,.true.)
call check_pointer_ (stq,this%numreact,this%numsp,.true.)
call check_pointer_ (ntemp,this%numreact,.true.)
call check_pointer_ (indx,this%numsp,.true.)
call check_pointer_ (namesp,this%numsp,.true.)
call check_pointer_ (stqloc,this%numsp,.true.)
call check_pointer_ (idaqprisp,this%numaqprisp,.true.)
call check_pointer_ (idreactsp,this%numreact,.true.)
!%------------------------------------------------------------
!% Allocate local reaction objects 
!%------------------------------------------------------------
if (this%numreact>0) then
 allocate(reaction(this%numreact))
end if
!%------------------------------------------------------------
!% Create and copy previous reaction objects
!% Destroy the prevous reaction objects 
!%------------------------------------------------------------
do i=1,this%numreact
 call create_ (reaction(i))
 reaction(i)=this%preaction(i)%ptr
 call get_logk_ (this%preaction(i)%ptr,logktemp,templogk,ithntemp,iserror)
 ntemp(i)=ithntemp
 logk(i,1:ithntemp)=logktemp
 do j=1,this%numsp
      species => this%pspecies(j)%ptr
      call set_pspecies_ (reaction(i),species)
 end do
 call destroy_ (this%preaction(i)%ptr)
end do 
!%------------------------------------------------------------
!% Build S1 and S2 
!%------------------------------------------------------------
ipri=0
isec=0 
do i=1,this%numsp
 ispri=.false. 
 name=this%pspecies(i)%ptr%name
 do j=1,this%numaqprisp
   if (name==nameaqprisp(j)) then
     ipri=ipri+1
     idaqprisp(ipri)=i
	 s1(:,ipri)=this%stq(:,i)
	 ispri=.true.
	 exit
   end if
 end do 
 if (.not.ispri) then
    isec=isec+1
	idreactsp(isec)=i
	s2(:,isec)=this%stq(:,i)
 end if
end do 
!%------------------------------------------------------------
!% Check the number of species 
!%------------------------------------------------------------
if ((ipri+isec)/=this%numsp) then
 iserror=.true. 
 msg='Error in the new set of primary species'
 goto 20
end if
!%------------------------------------------------------------
!% Swith the chemical base 
!%------------------------------------------------------------
call switch_base_  &
  (this%numaqprisp, &
   this%numreact, &
   mxntemp, &
   this%numsp, &
   s2, &
   s1, &
   logk, &
   s2lu, &
   stq, &
   indx) 
!%-------------------------------------------------------------------------
!% Create and set the reaction objects 
!%-------------------------------------------------------------------------
do i=1,this%numreact
  isps=1
  isec=idreactsp(i)
  ithntemp=ntemp(i)
  if (this%iskinreact(i)) then
    call get_prrlaw_ (reaction(i),rrlaw)
  end if 
  name=this%pspecies(isec)%ptr%name  
  namesp(isps)=name
  stqloc(isps)=-1.0d0
  do j=1,this%numaqprisp
    ipri=idaqprisp(j)
	if (dabs(stq(i,j))>zero) then
     isps=isps+1
     namesp(isps)=this%pspecies(ipri)%ptr%name
     stqloc(isps)=stq(i,j)
    end if
  end do
!%-------------------------------------------------------------------------
!% Set in the reaction object 
!%------------------------------------------------------------------------- 
  call create_ (this%preaction(i)%ptr)
  call set_ &
   (this%preaction(i)%ptr, &
    name, &
    this%iskinreact(i), &
    stqloc(1:isps), &
    namesp(1:isps), &
    logk(i,1:ithntemp), &
    templogk(1:ithntemp), &
    ithntemp, &
    ncoeff, &
    isps, &
    iserror, &
	rrlaw=rrlaw)
  call destroy_ (reaction(i))
end do
!%--------------------------------------------------------------
!% Build the chemical base
!%--------------------------------------------------------------
call build_chemical_base_(this,msg,iserror)
if (iserror) goto 20
!%--------------------------------------------------------------
!% Build the stoichiometric matrix
!%--------------------------------------------------------------
call build_stq_ (this,msg,iserror)
if (iserror) goto 20
!%----------------------------------------------------------------
!% Build the components matrix
!%----------------------------------------------------------------
if (this%isgaussjordan) then
  call build_components_matrix_gauss_jordan_(this,msg,iserror)
else
  call build_components_matrix_sing_values_(this,msg,iserror)
end if
if (iserror) goto 20
!%----------------------------------------------------------------
!% If the problem is isoterm update parameters that depend of the 
!% temperature
!%----------------------------------------------------------------
call update_ (this, this%tempref,iserror)
if (iserror) then
   msg='Error when calling update_'
   goto 20
end if
!%----------------------------------------------------------------
!%----------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate and nulllify local pointers 
!%------------------------------------------------------------
call check_pointer_ (s1,1,1,.false.)
call check_pointer_ (s2,1,1,.false.)
call check_pointer_ (logk,1,1,.false.)
call check_pointer_ (logktemp,1,.false.)
call check_pointer_ (templogk,1,.false.)
call check_pointer_ (stq,1,1,.false.)
call check_pointer_ (s2lu,1,1,.false.)
call check_pointer_ (indx,1,.false.)
call check_pointer_ (ntemp,1,.false.)
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (stqloc,1,.false.)
call check_pointer_ (idaqprisp,1,.false.)
call check_pointer_ (idreactsp,1,.false.)
rrlaw => null ()
species => null ()
surface => null ()
if (this%numreact>0) then 
 deallocate(reaction)
 reaction => null ()
end if
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: set_new_chemical_base_'
print *,'********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_new_chemical_base_pchemsys &
   (this, &
    nameaqprisp, &
	nameadsprisp, &
    c, &
	nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Get the name of the species of the new chemical base. 
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in), target      :: this           ! Type parent chemical system variable

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
 j, &
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
call get_iposspsph_ (this,this%aqphindex,isps1,isps2)
!%------------------------------------------------------------
!% Allocate pointers 
!%------------------------------------------------------------
call check_pointer_ (nameaqprisp,this%numaqprisp,.true.)
!%------------------------------------------------------------
!% Determine the new set of primary aqueous species
!%------------------------------------------------------------
do i=1,this%numaqprisp
  cmax=0.0d0
  do j=isps1,isps2
     name=this%pspecies(j)%ptr%name
     if (this%ueq(i,j)/=0.0d0) then
       if (c(j)>cmax) then
	     nameaqprisp(i)=name
		 cmax=c(j)
	   end if
	 end if 
  end do
end do
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: get_new_chemical_base_'
print *,'********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine build_components_matrix_sing_values_pchemsys &
   (this, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build components matrix according sigular values descomposition
!
!   $Arguments:
!
 
type (t_parentchemicalsystem)   :: &
 this
character(len=*)                :: &
 msg
logical                         :: &
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
real*8, pointer                 :: &
 aa(:,:), &
 sv(:), &
 vv(:,:), &
 work(:), &
 uu(:,:)
integer, pointer                :: &
 iord(:), &
 idreact(:)
integer                         :: &
 ireact, &
 ierr, &
 irank, &
 i, &
 j, &
 ii
 real*8                          :: &
 coeff
character(len=100)              :: &
 name 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%---------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------
call check_pointer_ (this%ueq,this%numaqprisp,this%numsp,.true.)
call check_pointer_ (idreact,this%numreact,.true.)
!%---------------------------------------------------
ireact=0
do i=1,this%numreact
 
 if (.not.this%iskinreact(i)) then
  ireact=ireact+1
  idreact(ireact)=i
 end if
 
end do
!%---------------------------------------------------
if (ireact>0) then
!%---------------------------------------------------
call check_pointer_ (aa,this%numsp,this%numsp,.true.)
call check_pointer_ (work,this%numsp,.true.)
call check_pointer_ (vv,this%numsp,this%numsp,.true.)
call check_pointer_ (sv,this%numsp,.true.)
call check_pointer_ (iord,this%numsp,.true.)
call check_pointer_ (uu,this%numsp,this%numsp,.true.)
!%--------------------------------------------------
do i=1,ireact
 do j=1,this%numsp
    name=this%pspecies(j)%ptr%name
    call get_coeff_stq_ (this%preaction(idreact(i))%ptr,name,coeff)
    aa(i,j)=coeff
  end do
end do
!%-------------------------------------------Descomposition
ierr=0
call svd &
    (this%numsp, &
     ireact, &
     this%numsp, &
     aa, &
     sv, &
     .false., &
     uu, &
     .true., &
     vv, &
     ierr, &
     work)
!%----------------------------------------------------------
irank = 0
do i=1,this%numsp
  if (dabs(sv(i))<1.0d-5) then
      irank=irank+1
      iord(irank) = i
  end if
end do
!%------------------------------------------------------------
if (irank.ne.this%numaqprisp) then
 msg='Error in rank of components matrix, should be:'
 call add_ (msg,this%numaqprisp)
 goto 10
end if
!%------------------------------------------------------------
do i=1,this%numaqprisp
    do j=1,this%numsp
      ii = iord(i)
      this%ueq(i,j)= vv(j,ii)
    end do
end do
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (aa,1,1,.false.)
call check_pointer_ (work,1,.false.)
call check_pointer_ (vv,1,1,.false.)
call check_pointer_ (sv,1,.false.)
call check_pointer_ (iord,1,.false.)
call check_pointer_ (uu,1,1,.false.)
 
else
 
   do i=1,this%numaqprisp
      this%ueq(i,this%idaqprisp(i)) = 1.0d0
   end do
 
end if
!%------------------------------------------------------------
if (this%ecindex/=0) this%ueq(:,this%ecindex)=0.0d0
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (idreact,1,.false.)
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
SUBROUTINE READ_PRIM_master25_pchemsys &
 (THIS,NPRI,NAQX,NADS,MXTOT,IPRINT ,IUQDB  ,IOUTPUT,NAPRI,NAAQX,NAADS,&
 PROPADS,ZT,A0T, WMolT  ,BdotTot,IOptXhl,iact,iserror)
 
IMPLICIT REAL*8 (A-H,O-Z)
!-------------------------------------------------------------------------
!
!   $Description: Read primary species of chemical database (master25)
!%
!%DESCRIPTION
!%   Read primary species of chemical database (master25) ans store their
!%   ionic radii and charges
!%
!%ARGUMENTS : SCALARS
!%   NPRI            Number of primary species of chemical system
!%   NAQX           Number of secondary species of chemical system
!%   MXTOT            Maximum total number of species
!%   IPRINT           Controls printing of messagges
!%   IUQDB            Chemical Database (master25) unit number
!%   IOUTPUT          Output file unit number
!%
!%ARGUMENTS : ARRAYS
!%   NAPRI             Names of primary species of chemical system
!%   NAAQX            Names of secondary species of chemical system
!%   ZT               Ionic charge of all species of chemical system
!%   A0T              Ionic radius of all species of chemical system
!%
!%INTERNAL VARIABLES: SCALARS
!%    CARG            master25 value of ionic charge
!%    NAMSP           master25 name of primary species
!%    RAD             master25 value of ionic radius
!%
!%HISTORY
!%   Created by J.Carrera and C.Ayora (Jan,1998)
!
!   $Arguments:
!
 
CHARACTER*20 NAMSP
CHARACTER(LEN=*) NAPRI,NAAQX,NAADS
 
DIMENSION NAPRI(NPRI), NAAQX(NAQX), ZT(NPRI+NAQX), A0T(NPRI+NAQX), &
          WMolT(NPRI+NAQX), BdotTot(NPRI+NAQX),PROPADS(1,NADS),NAADS(NADS)
Integer*4 IOptXhl 
 
logical :: iserror

type(t_parentchemicalsystem) :: this

character(len=100) msg
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
INTEGER, PARAMETER :: BDOTMOD=3, FRACMASSOPTION=2, ACTIVE=2,r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
! 
iserror=.false. 
msg=''
!%--------------------------------------- Begin reading primary species fr
1 If (IOptXhl==FRACMASSOPTION.or.IACT==BDOTMOD) Then
  READ (IUQDB, *, ERR=20) NAMSP, RAD, CARG, BDot, WMol
Else
  READ (IUQDB, *, ERR=10) NAMSP, RAD, CARG
EndIf
 
IF (IPRINT==ACTIVE) WRITE(IOUTPUT,*) NAMSP
 
IF (NAMSP=='null') RETURN
 
DO J=1,NPRI       ! Check if NAMSP is primary species of chemical sy
   IF (NAMSP==NAPRI(J)) THEN
      ZT(J)=CARG
      A0T(J)=RAD
      If (IOptXhl==FRACMASSOPTION) WMolT(J) = WMol
      If (IACT==BDOTMOD) Then
        BdotTot(J) = Bdot
        If (Bdot<r0) Then
           msg='Error, negative BDOT for species'
           call add_ (msg,NamSp)
           goto 30
        EndIf
      EndIf
      GOTO 1
   END IF
END DO
 
DO J=1, NAQX   ! Check if NAMSP is secondary species of chemical sys
   IF (NAMSP.EQ.NAAQX(J)) THEN
       ZT(J+NPRI)=CARG
       A0T(J+NPRI)=RAD
       If (IOptXhl==FRACMASSOPTION) WMolT(J+NPRI) = WMol
       If (IAct==BDOTMOD) Then
         BdotTot(J+NPRI) = Bdot
         If (Bdot<r0) Then
           msg='Error, negative BDOT for species'
           call add_ (msg,NamSp)
           goto 30
         EndIf
       EndIf
       GOTO 1
   END IF
END DO

DO J=1, NADS   ! Check if any adsorbed specie is primary 
   IF (NAMSP==NAADS(J)) THEN
       PROPADS(1,J)=CARG
       GOTO 1
   END IF
END DO

GOTO 1                             ! NAMSP is not part of chemical s
 
10 continue  

msg=' ERROR when reading prim. species from chem. database'
GOTO 30

return

20 continue
if (IOptXhl==FRACMASSOPTION) then
 msg='Error, molecular weight not defined in the data base'
 goto 30
end if
if (IACT==BDOTMOD) then
 msg='Error, BDOT parameter not defined in the data base'
 goto 30
end if
return

30 continue  
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: select_cmob_'
print *, msg
print *,'***********************'
iserror=.true.
return

END SUBROUTINE
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
SUBROUTINE READ_SCND_master25_pchemsys &
 (NPRI  ,NAQX ,NGAS   ,NGASIN ,NTEMP  ,MXBASIS,NBASIS &
 ,MXTOT  ,NREAC  ,IPRINT ,IOUTPUT,IUQDB  ,ALOGK &
 ,NAPRI   ,INDXPRI,NAAQX  ,INDXSEC,NAGAS   ,INDXGAS &
 ,INDX   ,NAMSEC ,STCPR  ,STCRC   ,ZT     ,A0T   ,WMolT &
 ,BDotTot,IOptXhl,IAct)
 
 
 
IMPLICIT REAL*8 (A-H,O-Z)
!-------------------------------------------------------------------------
!
!   $Description: Read secondary species of chemical database (QDB)
!%
!%DESCRIPTION
!%   Read secondary species of chemical database (QDB) ans store their
!%   ionic radii and charges
!%
!%ARGUMENTS : SCALARS
!%   NPRI            Number of primary species of chemical system
!%   NAQX           Number of secondary species of chemical system
!%   NGAS             Number of gasses in chemical system
!%   NGASIN           Location of first gas reaction
!%   MXTOT            Maximum total number of species
!%   IPRINT           Controls printing of messagges
!%   IUQDB            Chemical Database (QDB) unit number
!%   IOUTPUT          Output file unit number
!%
!%ARGUMENTS : ARRAYS
!%   NAPRI             Names of primary species of chemical system
!%   NAAQX            Names of secondary species of chemical system
!%   NAGAS             Names of gasses of chemical system
!%   INDXPRI(I)       Init 0, becomes 1 when Ith primary species is found
!%   INDXSEC(I)       Init 0, becomes 1 when Ith secndary spcies is found
!%   INDXGAS(I)       Init 0, becomes 1 when Ith gas is found in QDB
!%   ZT               Ionic charge of all species of chemical system
!%   A0T              Ionic radius of all species of chemical system
!%
!%INTERNAL VARIABLES: SCALARS
!%   IFAIL            Error indicator from  LEAST_SQUARES_FIT
!%   NCT              Num QDB species matched by chemical system species
!%   NSP              Number of QDB secondary species participating in rea
!%   CARG             QDB value of ionic charge
!%   RAD              QDB value of ionic radius
!%   IERR             number of incomplete reactions (missing species)
!%
!%INTERNAL VARIABLES: ARRAYS
!%   ALGEK            Equilib. constant at NTEMP temp. points (read from Q
!%   NAM              Names of species involved in reaction (read from QDB
!%   STC              Stoichiometric coefficients of reaction (read from Q
!%
!%HISTORY
!%   Created by J.Carrera and C.Ayora (Jan,1998)
!
!   $Arguments:
!
 
PARAMETER (MXLOC=100)
 
CHARACTER(LEN=*)  NAPRI   ,NAAQX  ,NAGAS ,NAMSEC 
CHARACTER(LEN=LEN(NAPRI)) NAM
 
DIMENSION ALOGK (MXTOT,NTEMP)  , NAPRI (NPRI)    , INDXPRI (NPRI) &
        , NAAQX (NAQX) , INDXSEC (NAQX), NAGAS (NGAS) &
        , INDXGAS (NGAS) &
        , INDX (MXBASIS),NAMSEC (MXTOT) &
        , STCPR (MXTOT,NPRI) , STCRC (MXTOT,MXTOT) &
        , ZT(NPRI+NAQX), A0T(NPRI+NAQX) &
        , WMolT (NPRI+NAQX), BdotTot (NPRI+NAQX) &
        , ALGEK (MXLOC)  , NAM (MXLOC)     , STC (MXLOC)
Integer*4 IOptXhl 
 
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

 
!%--------------------------------------------------------------- Initiali
 
NREAC=0                  ! Beginning of counter for this kind of rea
IERR=0
IF (MXLOC.LT.MXTOT) THEN
   WRITE (IOUTPUT,8001) MXTOT
8001    FORMAT(' ERROR: MXLOC in subroutine READ_SCND_master25_ has ', &
          'to be increased to',I3)
   STOP
END IF
IF (MXLOC.LT.NTEMP) THEN
   WRITE (IOUTPUT,8001) MXTOT
   STOP
END IF
 
!%******************************* Begin loop over all secondary species fr
 
1 NREAC=NREAC+1
If (IOptXhl==2 .or. IAct==3) Then
  READ(IUQDB,*,ERR=1001) NAM(1) , NSP , (STC(I+1) &
               ,NAM(I+1),I=1,NSP) ,(ALGEK(L),L=1,NTEMP),RAD, CARG &
               ,BDot, WMol
Else
  READ(IUQDB,*,ERR=1001) NAM(1) , NSP , (STC(I+1) &
               ,NAM(I+1),I=1,NSP) ,(ALGEK(L),L=1,NTEMP),RAD, CARG
EndIf
 
IF (IPRINT.EQ.2) WRITE(IOUTPUT,'(A)') NAM(1)
 
IF (NAM(1).EQ.'null') THEN       ! ends reading secondary species fr
   NREAC=NREAC-1
   IF(IERR.GT.0) THEN
      WRITE(IOUTPUT,8005) IERR
8005       FORMAT (' ERROR: The chemical database definition of',I2, &
      ' secondary species'/ &
      ' contains species missing in the chemical system')
      STOP
   END IF

   
   DO I=1,NAQX
    IF (INDXSEC(I)==0) THEN
     WRITE(IOUTPUT,*)'ERROR, NOT FOUND IN THE DATA BASE THE SPECIES:',NAAQX(I)
     STOP
    END IF 
   END DO
   
   
   RETURN
END IF
 
STC(1)=-1
NCT=0
 
!%---------------------------- Find if species NAM(1) is part of chemical
 
DO 100 K=1,NSP+1
 
!%------------- First, check if K-th species is a primary species of chem
   DO J = 1,NPRI
      IF (NAM(K).EQ.NAPRI(J)) THEN  ! It is J-th "    "     "   "
         NCT = NCT+1
         INDXPRI(J) = NREAC
         STCPR(NREAC,J) = STC(K)
         IF (K.EQ.1) THEN
            ZT(J)=CARG
            A0T(J)=RAD
            If (IOptXhl .eq. 2) WMolT(J) = WMol
            If (IAct .eq. 3) Then
              BdotTot(J) = Bdot
              If (Bdot .lt. 0.0) Then
                Write(IOutput,8002) Nam(K)(1:Index(Nam(K),' ')-1)
                Write(*,8002) Nam(K)(1:Index(Nam(K),' ')-1)
              EndIf
            EndIf
         END IF
         GOTO 100
      END IF
   END DO
 
!%------------- Now, check if K-th species is a secondary species of chem
   DO J=1, NAQX
      IF (NAM(K).EQ.NAAQX(J)) THEN  ! It is the J-th "    "     "
         NCT=NCT+1
         INDXSEC(J)=NREAC
         STCRC(NREAC,J)=STC(K)
         ALOGK(NREAC,1:NTEMP)=ALGEK(1:NTEMP)
         IF (K.EQ.1) THEN
            ZT(J+NPRI)=CARG
            A0T(J+NPRI)=RAD
            If (IOptXhl .eq. 2) WMolT(J+NPRI) = WMol
            If (IAct .eq. 3) THEN
              BdotTot(J+NPRI) = Bdot
              If (Bdot .lt. 0.0) Then
                Write(IOutput,8002) Nam(K)(1:Index(Nam(K),' ')-1)
                Write(*,8002) Nam(K)(1:Index(Nam(K),' ')-1)
              EndIf
            EndIf
         END IF
         GOTO 100
      END IF
   END DO
 
!%---------------------------------------------- Now check other (Type2) s
 
   DO J=1, NGAS
      IF (NAM(K).EQ.NAGAS(J)) THEN
         NCT=NCT+1
         INDXGAS(J)=NREAC
         STCRC(NREAC, NGASIN+J)=STC(K)
         ALOGK(NREAC,1:NTEMP)=ALGEK(1:NTEMP)
         GOTO 100
      END IF
   END DO
 
!%-------------------------------------------------------------- Check for
 
   IF (NAM(K).EQ.'H2O' .OR. NAM(K) .EQ. 'h2o') THEN
      NCT=NCT+1
      GOTO 100
   END IF
 
   IF (K.EQ.1) THEN      ! This species is not part of the chemical
      NREAC=NREAC-1
      GOTO 1
   END IF
 
!%-------------- If arrived here, it means that it has not found any speci
!%-------------- the chemical system matching the K-th reactant in this re
 
   WRITE (IOUTPUT, 8030) NAM(K), NAM(1)
8030    FORMAT('ERROR: species', A, ',which is part of reaction', A, &
          ',is missing in the chemical system')
 
100 CONTINUE
 
!%--------------------- Check if number of identified species equals that
IF (NCT .NE. NSP+1) THEN
   WRITE (IOUTPUT, 8040) NAM(1), NSP+1, NCT
   WRITE (*, 8040) NAM(1), NSP+1, NCT
8040    FORMAT(' Counting ERROR in secondary species ', A,/ &
          ' the number of species involved is', I2,/ &
          ' only', I2, ' have been found in the chemical system')
   IERR=IERR+1
END IF
 
NAMSEC (NREAC) = NAM(1)
 
 
!%----------------------------------------------- Go on to a new reaction
GOTO 1
 
1001 WRITE (IOUTPUT, 8010)
WRITE (*,8010)
8010 FORMAT (' ERROR when reading secondary species ', &
        'from chem. database')
8002 FORMAT (' WARNING: No Truesdell-Jones Bdot data for species ',A, &
        '. Instead, Bdot of Debye-Hueckel is used')
STOP
END SUBROUTINE
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
SUBROUTINE read_genr_master25_pchemsys &
 (NAMTYP ,NTP1   ,NPROP  ,NPRI  ,NAQX ,NTP2   ,NTP2IN ,NTEMP &
 ,MXBASIS,NBASIS ,MXTOT  ,NREAC  ,IPRINT ,IOUTPUT,IUQDB &
 ,NAMTP1 ,INDXTP1,VARTP1 ,ALOGK  ,NAPRI   ,INDXPRI,NAAQX &
 ,INDXSEC,NAMTP2 ,INDXTP2,NAMSEC ,STCPR &
 ,STCRC  ,IOPGASTR,VARTP2)
 
IMPLICIT NONE
!-------------------------------------------------------------------------
!
!   $Description: Read generic reactions from chemical database (QDB)
!%
!%DESCRIPTION
!%   Read  Mineral reactions (ITYPRC=1), Gases (ITYPRC=1) or adsorption
!%   reactions (ITYPRC=2); compute lnK, lnK-temp. coefficients and stores
!%   stoichiometric coefficients.
!%   When called for reading minerals, then type 1 reactions refer to mine
!%                                     and type 2 reactions refer to gasse
!%   When called for reading gasses, then type 1 reactions refer to gasses
!%                                     and type 2 reactions can be anythin
!%   When called for reading adsorbed species, then type 1 reactions refer
!%                     to adsorbed species, and type 2 reactions refer to
!%
!%ARGUMENTS : SCALARS
!%
!%   NAMTYP           Name of type of reactions to be read
!%   NTP1             Number of Type 1 reactions (minerals, gasses or adso
!%   NPROP              Number of VARTP parameters for each species
!%   NPRI            Number of primary species of chemical system
!%   NAQX           Number of secondary species of chemical system
!%   IPRINT           Controls printing of messagges
!%   IUQDB            Chemical Database (QDB) unit number
!%   IOUTPUT          Output file unit number
!%   NTP2             Number of Type 2 reactions (minerals, gasses or adso
!%   NTP2IN           Initial position of Type 2 reactions
!%   IOPGASTR         Indicates if gas is transported (adding dominant gas
!%
!%ARGUMENTS : ARRAYS
!%   NAMTP1           Name of Type 1 reactions
!%   INDXTP1(I)       Initially zero, becomes nonzero whenever the I-th
!%                    type 1 reaction has been identified in the QDB
!%   VARTP1           Real variable associated to type 1 reactions
!%                    (molar volume, charge, etc, depending on reaction)
!%   NAPRI             Names of primary species of chemical system
!%   NAAQX            Names of secondary species of chemical system
!%   INDXPRI(I)      Init 0, becomes 1 when Ith primary species is found i
!%   INDXSEC(I)      Init 0, becomes 1 when Ith secndary spcies is found i
!%   NAMTP2           Name of Type 2 reactions
!%   INDXTP2(I)       Initially zero, becomes nonzero whenever the I-th
!%                    type 2 reaction has been identified in the QDB
!%   VARTP2           Real variable associated to type 2 reactions
!%                    (molar volume, charge, etc, depending on reaction)
!%
!%INTERNAL VARIABLES: SCALARS
!%   IERR             Counter of number of incomplete reactions, missing s
!%   IFAIL            Error indicator from  LEAST_SQUARES_FIT
!%   NCT              Num QDB species matched by chemical system species
!%   NFILA            Row index of reaction
!%   NFSTIN           Initial position of Type 1 reactions
!%   NSP              Number of QDB secondary species participating in rea
!%   VMOL             QDB value of variable associated to type 1 reactions
!%                    (molar volume, charge, etc, depending on reaction)
!%
!%INTERNAL VARIABLES: ARRAYS
!%   ALGEK            Equilib. constant at NTEMP temp. points (read from Q
!%   NAM              Names of species involved in reaction (read from QDB
!%   STC              Stoichiometric coefficients of reaction (read from Q
!%
!%HISTORY
!%   Created by J.Carrera and C.Ayora (Jan,1998)
!
!   $Arguments:
!
 
INTEGER, PARAMETER::MXLOC=100
INTEGER NPRI,NTP1,NAQX,NTP2,MXTOT,MXBASIS,NTEMP,NPROP,NTP2IN, &
       NBASIS,IPRINT,IOUTPUT,IUQDB,IOPGASTR,NREAC
CHARACTER(LEN=*) NAMTP1(NTP1),NAPRI(NPRI),NAAQX(NAQX),NAMTP2(NTP2), &
                 NAMSEC(MXTOT),NAMTYP
INTEGER INDXTP1 (NTP1),INDXPRI (NPRI), INDXSEC (NAQX), &
        INDXTP2 (NTP2)
REAL*8  VARTP1 (NPROP, NTP1),ALOGK(MXTOT,NTEMP), &
        STCPR (MXTOT,NPRI),STCRC(MXTOT,MXTOT), &
        ALGEK (MXLOC),STC (MXLOC),VARTP2(NPRI)
 
! LOCAL VARIABLES
INTEGER NFSTIN,IERR,ITYPRC,NFILA,I,K,NSP,L,NCT,J
REAL*8  VPROP(NPROP),VBARGAS1,VBARGAS2,VMOL,ZD
CHARACTER(LEN=100) NAM(MXLOC) 
INTEGER, PARAMETER :: MINERAL=1, &
 SURFACE=2, &
 GAS=3
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
IOPGASTR=0 
!%------------------------------------------------------------------------
!% IF NTP1=0 THE FIND THE NEXT NULL AND RETURN  
!%------------------------------------------------------------------------
IF (NTP1==0) THEN 
 DO  
  READ(IUQDB,*,ERR=1001) NAM(1)
  IF (NAM(1)=='null') EXIT
 END DO 
 RETURN 
end if 
!%--------------------------------------------------------------- Initiali 
NFSTIN=NREAC           ! Beginning of counter for this kind of react
IERR=0
IF (MXLOC.LT.MXTOT) THEN
   WRITE (IOUTPUT,8001) MXTOT
8001    FORMAT(' ERROR: MXLOC in subroutine READ_GENR_master25_ has ', &
          'to be increased to',I3)
   STOP
END IF
 
IF (MXLOC.LT.NTEMP) THEN
   WRITE (IOUTPUT,8001) NTEMP
   STOP
END IF
 
IF (NAMTYP.EQ.'minerals') THEN    ! CAI&FBP
   ITYPRC=MINERAL
ELSE IF (NAMTYP.EQ.'adsorbed species') THEN
   ITYPRC=SURFACE
ELSE IF (NAMTYP.EQ.'gasses') THEN ! CAI&FBP
   ITYPRC=GAS                       ! CAI&FBP
ELSE
   WRITE (IOUTPUT,8002) NAMTYP
8002    FORMAT(' ERROR: NAMTYP in subroutine READ_GENR_master25_ is ',A20/ &
          ' and it must be: minerals, gasses or adsorbed species')
   STOP
END IF
 
!%**************************** Begin loop over all reactions of this type
 
1 NREAC=NREAC+1
!%---------------------------------------------------------------- Read fr
 
IF (ITYPRC.EQ.MINERAL) THEN                       ! Mineral
 
   READ(IUQDB,*,ERR=1001) NAM(1) , VMOL , NSP &
               ,(STC(I+1),NAM(I+1),I=1,NSP) ,(ALGEK(L),L=1,NTEMP)

   VPROP(1) = VMOL
  
ELSE IF (ITYPRC.EQ.SURFACE) THEN                  ! Adsorbed species type
   NAM='' 
   READ(IUQDB,*,ERR=1001) NAM(1) , NSP , (STC(I+1) &
               ,NAM(I+1),I=1,NSP) ,(ALGEK(L),L=1,NTEMP) ,ZD
   VPROP(1) = ZD
 
ELSE IF (ITYPRC.EQ.GAS) THEN                  ! Gases
 
   READ(IUQDB,*,ERR=1001) NAM(1) , VBARGAS1, VBARGAS2 , NSP &
               ,(STC(I+1),NAM(I+1),I=1,NSP) ,(ALGEK(L),L=1,NTEMP)
   VPROP(1) = VBARGAS1
   VPROP(2) = VBARGAS2
 
END IF

 
IF (IPRINT.EQ.2) WRITE(IOUTPUT,'(A)') NAM(1)
 
IF (NAM(1).EQ.'null') THEN         ! ends reading this reaction type
   NREAC=NREAC-1
   IF (IERR.GT.0) THEN
      WRITE(IOUTPUT,8011) IERR,NAMTYP
8011       FORMAT (' ERROR: The chemical database stoichiometric', &
              ' definition of',I2,' ',A/' contains species', &
              ' missing in you chemical system')
      STOP
   END IF
   
   DO I=1,NTP1
     IF(INDXTP1(I)==0) THEN
        WRITE (IOUTPUT,*) 'NOT FOUND IN THE DATA BASE THE SPECIES:',NAMTP1(I)
        STOP
     END IF 
   END DO 
   
   
   
   RETURN
END IF
 
!%--------------------------- Find if reaction NAM(1) is part of chemical
 
DO I=1,NTP1

         IF (NAM(1).EQ.NAMTP1(I)) THEN  ! NAM(1) is part of our chemical system!

!* ----------------------------------------------- Store everything OF NAM(1) !!!            
            NCT = 1    ! Initalize counter of NAM(1) species in our chem. system
            INDXTP1(I) = NREAC
            NFILA = NFSTIN+I
            NAMSEC (NFILA) = NAM(1)
            STCRC (NFILA,NFILA) = -1D0
            ALOGK(NFILA,1:NTEMP)=ALGEK(1:NTEMP)
	        DO J=1,NPROP
              VARTP1(J, I) = VPROP(J)
	        END DO



          IF (ITYPRC.EQ.GAS .AND. I.EQ.NTP1 .and. iopgastr.gt.0) THEN
            GOTO 1
          ELSE
            GOTO 100 ! We have found a reaction, go on to store involved species
          ENDIF
       
       END IF

END DO

DO I=1,NPRI

         IF (NAM(1).EQ.NAPRI(I)) THEN  ! NAM(1) is part of our chemical system!
           DO K=2,NSP+1
              DO J = 1,NTP1
                IF (NAM(K).EQ.NAMTP1(J)) THEN  
                  NCT = 1    ! Initalize counter of NAM(1) species in our chem. system
                  INDXTP1(J) = NREAC
                  NFILA = NFSTIN + J
                  NAMSEC (NFILA) = NAM(1)
                  STCPR (NFILA,I) = -1D0
                  ALOGK(NFILA,1:NTEMP)=ALGEK(1:NTEMP)
	              VARTP2(I) = VPROP(NPROP)
	              GOTO 100 ! We have found a reaction, go on to store involved species
                 ENDIF
              END DO  
            END DO
          END IF 

END DO 

   

 
!%NAM(1) is not part of our chemical system, go on to check next reaction
NREAC = NREAC-1
GOTO 1
 
 
 
!%***************************** Store stoichiometric coefficients of NSP s
100 CONTINUE
 
DO 200 K=2,NSP+1
 
!%------------- First, check if K-th species is a primary species of chem
   DO J = 1,NPRI
      IF (NAM(K).EQ.NAPRI(J)) THEN  ! It is J-th "    "     "   "
         NCT = NCT+1
         INDXPRI(J) = NREAC
         STCPR(NFILA,J) = STC(K)
         GOTO 200
      END IF
   END DO
 
!%------------- Now, check if K-th species is a secondary species of chem
   DO J=1, NAQX
      IF (NAM(K).EQ.NAAQX(J)) THEN  ! It is the J-th "    "     "
          NCT=NCT+1
          INDXSEC(J)=NREAC
          STCRC(NFILA,J)=STC(K)
          GOTO 200
      END IF
   END DO
 
!%---------------------------------------------- Now check other (Type2) s
 
   DO J=1, NTP2
      IF (NAM(K).EQ.NAMTP2(J)) THEN
         NCT=NCT+1
         INDXTP2(J)=NREAC
         STCRC(NFILA,NTP2IN+J)=STC(K)
         GOTO 200
      END IF
   END DO
 
!%-------------------------------------------------------------- Check for
 
   IF (NAM(K)=='H2O'.OR.NAM(K)=='h2o') THEN
      NCT=NCT+1
      GOTO 200
   END IF
 
!%-------------- If arrived here, it means that it has not found any speci
!%-------------- the chemical system matching the K-th reactant in this re
 
   WRITE (IOUTPUT, 8030) NAM(K), NAM(1)
8030    FORMAT('ERROR: species', A, ',which is part of reaction', A, &
          ',is missing in the chemical system')
 
200 CONTINUE


 
!%--------------------- Check if number of identified species equals that
IF (NCT .NE. NSP+1) THEN
   WRITE (IOUTPUT, 8040) NAM(1), NSP+1, NCT
8040    FORMAT('Counting ERROR in reaction ', A,/ &
          ' the number of species involved is',I2,/ &
          ' only', I2,1x,'were found in the chemical system ')
   IERR=IERR+1
END IF
 
!%----------------------------------------------- Go on to a new reaction
GOTO 1
 
1001 WRITE(IOUTPUT, 8008) NAMTYP, NAM(1)
8008 FORMAT('ERROR: format error when reading ',A,'line ',A,' in thermodynamic data base')
NREAC = NREAC - 1
STOP 
END SUBROUTINE
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%PURPOSE
!%   Read temperature points and build ATBSFN
!%
!%DESCRIPTION
!%   This routine reads temperature points (first line of chemical data ba
!%   builds ATBSFN, matrix containing the values of log K polynomials
!%   basis functions at temp. points. It also checks MXBASIS
!%
!%ARGUMENTS : SCALARS
!%   NTEMP            Output, Number of temperature points
!%   MXBASIS          Input, Maximum number of basis functions as specifie
!%                           in the calling routine
!%   NBASIS:          Output, number of "basis" functions of log K depend.
!%   IUQDB            Chemical Database (QDB) unit number
!%   IOUTPUT          Output file unit number
!%
!%ARGUMENTS : ARRAYS
!%   ATBSFN           Matrix of basis functions evaluated at temperature p
!%
!%INTERNAL VARIABLES: SCALARS
!%   NAME             Text in first line of data base
!%   TK               Constant for changing deg Celsius to deg Kelvin(TK=2
!%
!%INTERNAL VARIABLES: ARRAYS
!%   TEMPC            Temperature data points
!%
!%HISTORY
!%   Created by J.Carrera and C.Ayora (Jan,1998)
!%************************************************************
SUBROUTINE read_temp_master25_pchemsys &
 (NTEMP,TEMPC,NBASIS,MXDIM,IUQDB,IOUTPUT)
 
IMPLICIT REAL * 8 (A-H,O-Z)
 
DIMENSION TEMPC(MXDIM)
 
CHARACTER * 20 NAME
 
!%------------------------------------------------------ Read temperature
 
READ (IUQDB, *, ERR=6001) NAME, NTEMP, (TEMPC(I), I=1,NTEMP)
 
!%------------------------------------------------------------- Check dime
 
IF (NTEMP.GT.MXDIM) THEN
   WRITE (*,8002) NTEMP
   WRITE (IOUTPUT,8002) NTEMP
8002    FORMAT ('ERROR:  increase MXDIM to',I3, &
                   'in subroutine READ_TEMP_master25_')
   STOP
END IF
!%---------------------------------------------------------- Build Matrix
IF (NTEMP>1) THEN
   NBASIS=5
ELSE
   NBASIS=1
END IF
 
RETURN
 
6001 WRITE (*,8010)
WRITE (IOUTPUT,8010)
8010 FORMAT (' ERROR when reading temp. points (first line) in QDB')
STOP
 
END SUBROUTINE
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
SUBROUTINE READ_AQSP_SYS_PCHEMSYS &
  (ITEMP,IACT,ICONV,NAQT,NPRI,TC2,LABEL,NAAQT,MXTOT,MXLABEL,ISERROR)
  
implicit double precision (a-h,o-z)  
implicit integer (i-n)
LOGICAL, INTENT(OUT) :: ISERROR
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
  
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
CHARACTER(LEN=*) LABEL(MXLABEL)
CHARACTER(LEN=*) NAAQT(MXTOT) 
 
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
 
ISERROR=.FALSE. 
 
  
 
!%      read control parameters of the chemical system
 
  READ (1,'(A)',ERR=9001) LABEL(1)      !title of the problem
 
!%      read aqueous species of the system
  READ (1,'(A)',ERR=9002) LABEL(2)      !---------------------------
 
  READ (1,'(A)',ERR=9002) LABEL(3)      !aqueous species
 
  READ (1,'(f10.0,2i5)',ERR=9002) TC2,IACT,ICONV        !temp. to calculate
 
  IF (TC2.EQ.25.0D0) THEN
    ITEMP=0
  ELSE
    ITEMP=1
  ENDIF
 
 
        READ (1,'(A)',ERR=9003) LABEL(4)      !aqueous primary species
            
        I=0
130     CONTINUE
        I=I+1
        READ (1,*,ERR=9003) NAAQT(I)
        IF (NAAQT(I).EQ.'*') THEN
          NPRI=I-1
          GO TO 120
        ENDIF
        GO TO 130
120     CONTINUE

!%       READ AQUEOUS COMPLEXES
        READ (1,'(A)',ERR=9004) LABEL(5)             !aqueous complexes
           
        J=NPRI
132     CONTINUE
        J=J+1
        READ (1,*,ERR=9004) NAAQT(J)
        IF (NAAQT(J).EQ.'*') THEN
          NAQT=J-1
          GO TO 122
        END IF
        GO TO 132
122     CONTINUE

        RETURN

9001    WRITE(*,*) 'ERROR READING TITLE'
        ISERROR=.TRUE.
        RETURN
9002    WRITE(*,*) 'ERROR READING TEMPERATURE AND THERMODYNAMIC MODEL'
        ISERROR=.TRUE.
        RETURN
9003    WRITE(*,*) 'ERROR READING PRIMARY AQUEOUS species'
        ISERROR=.TRUE.
        RETURN
9004    WRITE(*,*) 'ERROR READING AQUEOUS COMPLEXES'
        ISERROR=.TRUE.
        RETURN
 
END SUBROUTINE
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
SUBROUTINE READ_MIGA_SYS_PCHEMSYS &
  (NGAS,NMIN, NMINEQ, NMINKIN, NAQT, NAAQT, &
    IMEQ, LABELM, NAGAS,NAMIN,MXLABEL,MXSP,iserror)
 
implicit double precision (a-h,o-z)  
implicit integer (i-n)
logical, intent(out)   :: iserror
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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


 
   CHARACTER(LEN=*) LABELM(MXLABEL)
   CHARACTER(LEN=*) NAMIN(MXSP),NAGAS(MXSP), NAAQT(MXSP)
   DIMENSION IMEQ(MXSP)
 
iserror=.false. 
!FBP ini
!%        read minerals
         READ (1,'(A)',ERR=9001) LABELM(1)             !minerals
              
         M=0
         NMINEQ=0
         NMINKIN=0
140      CONTINUE
         M=M+1
         READ (1,*,ERR=9002) NAMIN(M),IMEQ(M)
         IF (NAMIN(M).EQ.'*') GOTO 150
         IF (IMEQ(M).EQ.0) THEN 
            NMINEQ=NMINEQ+1
         ELSE
           NMINKIN=NMINKIN+1
         ENDIF
         GO TO 140
150      CONTINUE
         NMIN=NMINEQ+NMINKIN


!%       read kinetic data from kinetics.dat file 
!%        IF (NMINKIN.GT.0) CALL DATAKIN &
!%        (NAQT,  NMIN, NMINKIN, NAAQT, IMEQ, NAMIN,       &
!%           EA,  NKIN, THRESH, DINS,   FORA,    RK, PCAT, &
!%         NCAT, NACAT, filebase, kin_filename)  !FBP added kin_filename

!%       read gases
        READ (1,'(A)',ERR=9003) LABELM(3) 		!gases
        K=0
170     CONTINUE
        K=K+1
        READ (1,*,ERR=9004) NAGAS(K)

        IF (NAGAS(K).EQ.'*') THEN
          NGAS=K-1

          GO TO 190
        END IF
        GO TO 170
190     CONTINUE

        RETURN

9001    WRITE(*,*) 'ERROR READING HEADING OF MINERALS OF THE SYSTEM'
        iserror=.true.
        return
9002    WRITE(*,*) 'ERROR READING NAMES OF THE MINERALS OF THE SYSTEM'
        iserror=.true.
        return
9003    WRITE(*,*) 'ERROR READING HEADING OF GASES OF THE SYSTEM'
        iserror=.true.
        return
9004    WRITE(*,*) 'ERROR READING NAMES OF THE GASES OF THE SYSTEM'
        iserror=.true.
        return
 
END SUBROUTINE
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
SUBROUTINE READ_FROM_DATAKIN_PCHEMSYS &
 (NPRI,  & 
  NAQX, &
  NMIN,  &
  NMINKIN, &
  MXTERM, &
  NAPRI, & 
  NAAQX, & 
  IMEQ, & 
  NAMIN, &
  EA, & 
  NKIN, & 
  THRESH, & 
  DINS, &   
  FORA, &    
  RK, & 
  PCAT, &
  NCAT, &  
  NACAT, & 
  FILEBASE, &
  NAMEFILE, &
  ISERROR) 
 
IMPLICIT NONE 
!-------------------------------------------------------------------------
!
!   $Description: Read kinetic data base
!
!   $Arguments:
!

INTEGER, INTENT(IN)                                                      :: NMIN

INTEGER, INTENT(IN)                                                      :: MXTERM

INTEGER, INTENT(IN)                                                      :: NPRI

INTEGER, INTENT(IN)                                                      :: NAQX

INTEGER, INTENT(IN)                                                      :: NMINKIN

CHARACTER(LEN=*), INTENT(IN), DIMENSION(NMIN)                            :: NAMIN

CHARACTER(LEN=*), INTENT(IN), DIMENSION(NPRI)                            :: NAPRI

CHARACTER(LEN=*), INTENT(IN), DIMENSION(NAQX)                            :: NAAQX

CHARACTER(LEN=*), INTENT(OUT), DIMENSION(NMIN,MXTERM,NPRI+NAQX+NMIN)     :: NACAT

INTEGER, INTENT(IN), DIMENSION(NMIN)                                     :: IMEQ

INTEGER, INTENT(OUT), DIMENSION(NMIN)                                    :: NKIN

REAL*8, INTENT(OUT), DIMENSION(NMIN)                                     :: EA

REAL*8, INTENT(OUT), DIMENSION(NMIN)                                     :: THRESH

REAL*8, INTENT(OUT), DIMENSION(NMIN,5)                                   :: DINS

REAL*8, INTENT(OUT), DIMENSION(NMIN,MXTERM)                              :: FORA

REAL*8, INTENT(OUT), DIMENSION(NMIN,MXTERM)                              :: RK

REAL*8, INTENT(OUT), DIMENSION(NMIN,MXTERM,NPRI+NAQX+NMIN)               :: PCAT

INTEGER, INTENT(OUT), DIMENSION(NMIN,MXTERM)                             :: NCAT
   
CHARACTER(LEN=*), INTENT(IN)                                             :: NAMEFILE 

CHARACTER(LEN=*), INTENT(IN)                                             :: FILEBASE

LOGICAL, INTENT(OUT)                                                     :: ISERROR
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
CHARACTER(LEN=20)  :: &
NAMEDUM, &
NACATDUM(MXTERM,MXTERM)
CHARACTER(LEN=80)  :: &
REFERENCE, &
LABEL, &
LBAUX, &
NAMEFILELOC 
INTEGER     :: &
IMKIN(50), &
NCATDUM(MXTERM), &
M, &
K, &
J, &
KK, &
MINCOUNT, &
NKINDUM, &
LBASE
REAL*8         :: &
PCATDUM(MXTERM,MXTERM), &
RKDUM(MXTERM), &
DINSDUM(MXTERM), &
FORADUM(MXTERM)
REAL*8         :: &
EADUM, &
THRESHDUM
!-------------------------------------------------------------------------
!
!   $code
!
 
ISERROR=.FALSE. 
 
IF(NMINKIN==0) RETURN 
IF (NAMEFILE==' ') THEN
	NAMEFILELOC='kinetics.dat'
ELSE
    NAMEFILELOC=NAMEFILE 
ENDIF
    call lastletter_ (lbase,filebase)
    lbAux=filebase(1:lbase)//namefileloc
	open(unit=26,file=lbAux,status='old',err=5000)
    write(6,*) '=======> Reading database:',namefileloc 
!%---------------------------------------
!%      initialize flag
!%---------------------------------------
       DO M=1,NMIN
         IMKIN(M)=0
       END DO
!%---------------------------------------      
!%     read database
!%---------------------------------------
       READ (26,*)

       MINCOUNT=0
20     CONTINUE
       IF (MINCOUNT.GE.NMINKIN) GO TO 4000
       READ (26,*,ERR=5002) NAMEDUM,NKINDUM,EADUM,THRESHDUM,REFERENCE
       IF (NAMEDUM.EQ.'NULL'.OR.NAMEDUM.EQ.'null') GO TO 4000
       DO K=1,NKINDUM
         READ (26,*,ERR=5002) RKDUM(K),DINSDUM(K),FORADUM(K),NCATDUM(K), &
         (NACATDUM(K,KK),PCATDUM(K,KK),KK=1,NCATDUM(K))
       END DO
!%---------------------------------------
!%      assign kinetic values to minerals
!%---------------------------------------  

       DO 100 M=1,NMIN
         IF (IMEQ(M).EQ.0) GO TO 100
         IF (NAMEDUM.EQ.NAMIN(M)) THEN
           MINCOUNT=MINCOUNT+1
           IMKIN(M)=1
           NKIN(M)=NKINDUM
           EA(M)=EADUM
           THRESH(M)=THRESHDUM
           DO K=1,NKINDUM
             RK(M,K)=RKDUM(K)
             
             DINS(M,K)=DINSDUM(K)
             FORA(M,K)=FORADUM(K)

             NCAT(M,K)=NCATDUM(K)

             DO 50 KK=1,NCATDUM(K) 
               DO J=1,NPRI
                 IF (NACATDUM(K,KK).EQ.NAPRI(J)) THEN
                   PCAT(M,K,J)=PCATDUM(K,KK)
                   NACAT(M,K,J)=NACATDUM(K,KK)
                   GO TO 50
                 END IF
               END DO
               DO J=1,NAQX
                 IF (NACATDUM(K,KK).EQ.NAAQX(J)) THEN
                   PCAT(M,K,NPRI+J)=PCATDUM(K,KK)
                   NACAT(M,K,NPRI+J)=NACATDUM(K,KK)
                   GO TO 50
                 END IF
               END DO
               WRITE (*,*) 'CATALYST IS NOT AQ.SPECIES: ',NACATDUM(K,KK)
               ISERROR=.TRUE. 
			   RETURN 
50           CONTINUE
            END DO  
          END IF   
100     CONTINUE

        GO TO 20
4000    CONTINUE
!%--------------------------------------------------------------
!%       check if all minerals are in the database kinetics.dat
!%--------------------------------------------------------------
        DO M=1,NMIN
         IF (IMEQ(M).GT.0.AND.IMKIN(M).EQ.0) THEN
           WRITE (*,*) ' MINERAL MISSING IN KINETICS.DAT:',NAMIN(M)        
           WRITE (2,*)
           WRITE (2,*) ' MINERAL MISSING IN KINETICS.DAT:',NAMIN(M)        
           ISERROR=.TRUE.
		   RETURN
         END IF
        END DO

        CLOSE (UNIT=26)
        RETURN

5000    WRITE (*,*) ' ERROR OPENING',NAMEFILELOC
        ISERROR=.TRUE.
		RETURN
5001    WRITE (*,*) ' ERROR READING',NAMEFILELOC
        ISERROR=.TRUE.
		RETURN
5002    WRITE (*,*) 'ERROR READING DATA OF MINERAL ',NAMEDUM
        ISERROR=.TRUE.
		RETURN
 
END SUBROUTINE
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
SUBROUTINE READ_ADS_SYS_PCHEMSYS(NADS,NMOD,NPRI,NAPRI,NAADS, &
                           NAMOD,NADSMOD,MXDIM,MXLABEL,iserror) 

IMPLICIT NONE
!******************************************************************
!  INTERNAL VARIABLES:SCALARS
!    NAADSDUM 
!  EXTERNAL VARIABLES:SCALARS
!    NADS         number of total adsorbed species
!    NSTT         number of total primary species of adsorption
!    NMOD         number of adsorption models                           
!  EXTERNAL VARIABLES:ARRAYS
!    NADSMOD         adsorption model name
!    NST          number of primary species of adsorption for each
!                 model
!    NADSMOD      number of adsorbed species for each model
!    NAADS        name of total adsorbed species
!    LABELS
!  EXTERNAL VARIABLES:WORKSPACE
!    NAST         name of primary species for each model
!    NAADSMOD     name of adsorbed species for each model
!********************************************************************          
        INTEGER MXLABEL,MXDIM,NADS,NMOD,NPRI
        CHARACTER(LEN=*) NAADS(MXDIM),NAMOD(MXDIM),NAPRI(MXDIM)
        CHARACTER(LEN=20)JADSDUM(MXDIM),LABELS(MXLABEL),NAADSDUM
        INTEGER NADSMOD(MXDIM)
        INTEGER I,K,J  
        LOGICAL ISPRI,iserror
        
        iserror=.false.  

!       read surface complexes
        READ (1,'(A)',ERR=9001) LABELS(1)              !surface complexes
   
        NMOD=0
        NADS=0
        I=0 
 
210     READ (1,*,ERR=9002) NAADSDUM   
        IF (NAADSDUM.EQ.'*') GO TO 250
215     I=I+1
        NMOD=I
        JADSDUM(I)= NAADSDUM
        K=0
        J=0

!       Read kinetics model
        IF (NAADSDUM.EQ.'EX+KIN') THEN
           JADSDUM(I)='EX'
        ENDIF
        IF (NAADSDUM.EQ.'TL+KIN') THEN
           JADSDUM(I)='TL'

        ENDIF
        IF (NAADSDUM.EQ.'DL+KIN') THEN
           JADSDUM(I)='DL'
          
	  ENDIF
        IF (NAADSDUM.EQ.'CC+KIN') THEN
           JADSDUM(I)='CC'
	   
        ENDIF
        IF (NAADSDUM.EQ.'EQ+KIN') THEN
           JADSDUM(I)='EQ'
	   
        ENDIF
        IF (NAADSDUM.EQ.'LA+KIN') THEN
           JADSDUM(I)='LA'
	    
        ENDIF
        IF (NAADSDUM.EQ.'FR+KIN') THEN
           JADSDUM(I)='FR'
           
	  ENDIF
        IF (NAADSDUM.EQ.'KD+KIN') THEN
           JADSDUM(I)='KD'
          
	  ENDIF

      
       NAMOD(I)=JADSDUM(I)
      
!       Read adsorbate names + kinetic constand + kinetic order
       
       K=0
       ISPRI=.TRUE.
220    READ (1,*,ERR=9003) NAADSDUM
	   	
		
		  IF (NAADSDUM.EQ.'EX'.OR.NAADSDUM.EQ.'TL'.OR. &
           NAADSDUM.EQ.'DL'.OR.NAADSDUM.EQ.'CC'.OR. &
           NAADSDUM.EQ.'FR'.OR.NAADSDUM.EQ.'LA'.OR. &
           NAADSDUM.EQ.'EQ'.OR.NAADSDUM.EQ.'KD'.OR. &
           NAADSDUM.EQ.'EX+KIN'.OR.NAADSDUM.EQ.'TL+KIN'.OR. &
           NAADSDUM.EQ.'DL+KIN'.OR.NAADSDUM.EQ.'CC+KIN'.OR. &
           NAADSDUM.EQ.'FR+KIN'.OR.NAADSDUM.EQ.'LA+KIN'.OR. &
           NAADSDUM.EQ.'EQ+KIN'.OR.NAADSDUM.EQ.'KD+KIN') THEN

          GO TO 215
            
		  END IF
     
          IF (NAADSDUM .EQ. '*') GOTO 250
	      IF (ISPRI) THEN 
	       ISPRI=.FALSE.
	       SELECT CASE (NAMOD(I))
	       CASE('EX','TL','DL','CC')
	        NPRI=NPRI+1
		    NAPRI(NPRI)=NAADSDUM
		   END SELECT 
	      ELSE
	       NADS=NADS+1
	       K=K+1
		   NAADS(NADS)=NAADSDUM   
		   NADSMOD(I)=K 
		  END IF 
 
          GO TO 220

        
 
 
250    CONTINUE
       IF (NADS==0.AND.NMOD/=0) GOTO 9003


       RETURN

9001   WRITE(*,*) 'ERROR READING HEADING SURF. COMPL. OF THE SYSTEM'
       iserror=.true.
       return

9002   WRITE(*,*) 'ERROR READING MODEL OF ADSORPTION OF THE SYSTEM'
       iserror=.true.
       return

9003   WRITE(*,*) 'ERROR READING NAMES OF SURF, COMPL. OF THE SYSTEM'
       iserror=.true.
       return

END SUBROUTINE 
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
SUBROUTINE switch_base1_pchemsys &
  (NPRI, &
   NREAC, &
   NTEMP, &
   NSP, &
   STCRC, &
   STCPR, &
   ALOGK, &
   S2LU, &
   STQT, &
   INDX)

!********************************************************************************
!*
!* PURPOSE
!*    Compute STQT, EKT, and EKTCOEF
!* DESCRIPTION
!*    Change basis (set of primary species) from that of
!*    database (QDB) to that of the chemical system (input). Equilibrium
!*    constants and logK-temp. coefficients are changed accordingly
!*
!*    Theory:
!*
!*      Given  S1 (STCPR) stoichiometric coefficients of primary variables
!*             S2 (STCRC)       "             "       "  secondary variables
!*      Mass action law can be written as:
!*
!*             S1*log(a1)+S2*log(a2)=log(K)
!*
!*      We wish to write it as:
!*
!*             log(a2)=S1'*log(a1)+log(K')
!*
!*      where S1' is the solution of (-S2)*S1'=S1
!*            log(K') is the solution of (-S2)*log(K')=log(K)
!*
!*      The output of this routine is STQT=S1', EKT=log(K'), 
!*      together with ETKCOEF (coefficients of logK-temp. polynomials),
!*      which  is the solution of (-S2)*EKTCOEF=-COEF
!*
!* ARGUMENTS : SCALARS
!*    NPRI            Number of primary species of chemical system
!*    NAQX           Number of secondary species of chemical system
!*    MXTOT            Maximum total number of species
!*    IOUTPUT          Output file unit number
!*    IPRINT           Controls printing of messagges 
!*    NBASIS           Number of functions in temperature interpolation
!*    NTEMP            Number of temperature data points for interpolation
!*    NREAC            Total number of reactions 
!*
!* ARGUMENTS : ARRAYS 
!*    STQT             Stoichiometric matrix expressed in the primary species of
!*                     the chemical system (S1' above)
!*    EKTCOEF          Temperature Interpolation Coefficients of all reactions
!*    EKT              Equilibrium constants of all reactions
!*    STCPR(I,J)       Stoichiometric coefficient of J-th primary species of QDB
!*                        in J-th reaction (S1 above)
!*    STCRC(I,J)       Stoichiometric coefficient of J-th secondary species 
!*                        of QDB (including mineral, gases,etc) in J-th reaction 
!*                        (S2 above)
!*    INDX             Integer Workspace
!*    EKTCOEF2(I)      Dummy vector with the values of EKTCOEF for one sec. species
!*
!* HISTORY
!*    Created by J.Carrera and C.Ayora (Jan,1998)
!********************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
     


      DIMENSION  STCPR(NREAC,NPRI) , STCRC(NREAC,NREAC), &
                 S2LU(NREAC,NREAC) , STQT(NREAC,NSP), &
                 ALOGK(NREAC,NTEMP), &
                 INDX(NREAC)
	  CHARACTER*100 MSG
	  LOGICAL ISERROR

!-------------------------------------Change sign of STCRC

      DO I=1,NREAC
         DO J=1,NREAC
            S2LU(J,I)=-STCRC(J,I)
            STCRC(J,I)=-STCRC(J,I)
         END DO
      END DO

! -------------------------------------L U Decomposition of S2LU

      CALL LUDCMP (S2LU,NREAC,NREAC,INDX,DW,MSG,ISERROR)

!* -------------------------------------Solve S2LU*STQT=STCPR

      DO I=1,NPRI
         DO J=1,NREAC
            STQT(J,I)=STCPR(J,I)
         END DO
      END DO
      DO I=1,NPRI
         CALL LUBKSB (S2LU,NREAC,NREAC,INDX,STQT(1:NREAC,I))
		 CALL MPROVE (STCRC,S2LU,NREAC,NREAC,INDX,STCPR(1:NREAC,I),STQT(1:NREAC,I))
		 CALL MPROVE (STCRC,S2LU,NREAC,NREAC,INDX,STCPR(1:NREAC,I),STQT(1:NREAC,I))
		 CALL MPROVE (STCRC,S2LU,NREAC,NREAC,INDX,STCPR(1:NREAC,I),STQT(1:NREAC,I))		  
      END DO

!* -------------------------------------Solve S2LU*EKT=ALOGK
   
      DLN=DLOG(10D0)
      DO J=1,NTEMP 
	   CALL LUBKSB (S2LU,NREAC,NREAC,INDX,ALOGK(1:NREAC,J))
      END DO 
      
      EKT = DLN*EKT

      RETURN
END SUBROUTINE
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_inv_points_pchemsys &
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
! the water activity. 
! Also compute the constrained water activity.
! The method impremented was presented by Risacher and Clement (2001).  
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in) :: this     ! Type parent chemical system variable

integer, intent(in)                       :: nreact   ! Number of reactions 

integer, intent(in), dimension(nreact)    :: idreact  ! Global index of reactions thta must be checked 

logical, intent(out), dimension(nreact)   :: isconstr ! isconstr=true, then the reaction constraint the water activity  

real*8, intent(out)                       :: awconstr ! Value of the constrained water activity 

logical, intent(out)                      :: iserror  ! iserror=true, then there was an error

character(len=*), intent(out)             :: msg      ! Message error  
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                   :: &
 i, &
 j, &
 irow, &
 icw, &
 ispsw, &
 ireact, &
 nindp, &
 nindp1, &
 nnindp, &
 ifail, &
 irank
logical                                   :: &
 isnonsingular 
real*8, pointer                           :: &
 mat1(:,:) => null (), &
 mat(:,:) => null (), &
 mat2(:,:) => null (), &
 stq(:,:) => null (), &
 stq1(:,:) => null (), &
 identity(:,:) => null (), &
 work(:), &
 logk(:) => null (), &
 mat3(:,:) => null (), &
 b3(:) => null ()
real*8              :: &
 det, &
 sigma   
integer, pointer             :: &
 nonindp(:) => null (), &
 indp(:) => null (), &
 indp1(:) => null ()
real*8, parameter            :: &
 singular=0.0d0, &
 tol=1.0d-6 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Zeroing the constrained water activity 
!%------------------------------------------------------------
awconstr=0.0d0 
!%------------------------------------------------------------
!% If there are not reactions then return 
!%------------------------------------------------------------
if (nreact==0) return 
!%------------------------------------------------------------
!%Find water component index
!%------------------------------------------------------------
call get_chem_info_ (this,msg,iserror,wcindex=icw,wspindex=ispsw)
if (iserror) goto 10 
if (ispsw==0) return 
!%------------------------------------------------------------
!%Allocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (mat1,nreact,nreact,.true.)
call check_pointer_ (mat2,nreact,nreact,.true.)
call check_pointer_ (mat3,this%numsp,this%numsp,.true.)
call check_pointer_ (b3,this%numsp,.true.)
call check_pointer_ (mat,nreact,nreact,.true.)
call check_pointer_ (stq,nreact,this%numsp,.true.)
call check_pointer_ (stq1,nreact,this%numsp,.true.)
call check_pointer_ (logk,nreact,.true.)
call check_pointer_ (identity,nreact,nreact,.true.)
call check_pointer_ (work,4*this%numsp,.true.)
call check_pointer_ (indp,nreact,.true.)
call check_pointer_ (indp1,nreact,.true.)
call check_pointer_ (nonindp,nreact,.true.)
!%------------------------------------------------------------
!%Build the stoichiometric matrix 
!%------------------------------------------------------------
do i=1,nreact 
 ireact=idreact(i)
 call get_lnk_ (this%preaction(ireact)%ptr,logk(i))
 stq(i,:)=this%stq(ireact,:)
 stq(i,this%idreactsp(ireact))=0.0d0
 identity(i,i)=1.0d0
end do
stq1=stq
!%------------------------------------------------------------
!%Zeroing the water comlumn 
!%------------------------------------------------------------
stq(:,ispsw)=0.0d0
!%------------------------------------------------------------
!%Compute matrix m 
!%------------------------------------------------------------
mat=matmul(stq,transpose(stq))
!%------------------------------------------------------------
!%Deterine the independent and dependent set of reactions. 
!%------------------------------------------------------------
nindp=0
nindp1=0
nnindp=0
mat1=mat
ifail=1
!%------------------------------------------------------------
!%Compute determinant
!%------------------------------------------------------------
call f03aaf(mat1,nreact,nreact,det,work(1:nreact),ifail) 
if (det/=singular) then 
  nindp1=nreact
  do i=1,nindp1
   indp1(i)=i 
  end do
else if(det==singular.and.nreact>2) then  
  mat1=mat
  nindp=0
  nnindp=0
 do irow=1,nreact
  mat1(irow,:)=identity(irow,:)
  mat2=mat1
  ifail=1
  !%------------------------------------------------------------
  !%Compute determinant
  !%------------------------------------------------------------
  call f03aaf(mat2,nreact,nreact,det,work(1:nreact),ifail)
  if (det==singular) then
   nindp=nindp+1
   indp(nindp)=irow
  else
   nnindp=nnindp+1
   nonindp(nnindp)=irow
   mat1(irow,:)=mat(irow,:)
  end if
 end do
!%------------------------------------------------------------
!%------------------------------------------------------------
!%Check the independent set of reactions 
!%------------------------------------------------------------
!%------------------------------------------------------------
 do i=1,nindp
   isnonsingular=.true.
   mat1=identity
   mat1(indp(i),:)=mat(indp(i),:)
   do j=1,nnindp
     mat2=mat1
     mat2(nonindp(j),:)=mat(nonindp(j),:)
     ifail=1
     !%------------------------------------------------------------
     !%Compute determinant
     !%------------------------------------------------------------
     call f03aaf(mat2,nreact,nreact,det,work(1:nreact),ifail)
     if (det==singular) then 
       isnonsingular=.false.
       exit
     end if 
   end do 
   if (isnonsingular) then
     nindp1=nindp1+1
     indp1(nindp1)=indp(i)  
   end if  
 end do
end if
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
!%Storage the information 
!%------------------------------------------------------------
isconstr=.true.
do i=1,nindp1 
 ireact=indp1(i)
 isconstr(ireact)=.false. 
 stq1(ireact,:)=0.0d0
end do
!%-------------------------------------------------------------
!%Compute constrained water activity 
!%-------------------------------------------------------------
if (nindp1<nreact) then 
 mat3=matmul(transpose(stq1),stq1)
 b3=matmul(transpose(stq1),logk)
 ifail=0
 call f04jaf(this%numsp,this%numsp,mat3,this%numsp,b3,tol, &
             sigma,irank,work,4*this%numsp,ifail)
 awconstr=dexp(b3(ispsw))
end if 
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!%Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (mat1,1,1,.false.)
call check_pointer_ (mat2,1,1,.false.)
call check_pointer_ (mat,1,1,.false.)
call check_pointer_ (stq,1,1,.false.)
call check_pointer_ (stq1,1,1,.false.)
call check_pointer_ (logk,1,.false.)
call check_pointer_ (identity,1,1,.false.)
call check_pointer_ (indp,1,.false.)
call check_pointer_ (indp1,1,.false.)
call check_pointer_ (nonindp,1,.false.)
call check_pointer_ (b3,1,.false.)
call check_pointer_ (mat3,1,1,.false.)
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
subroutine select_cmineq_pchemsys &
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
 
type (t_parentchemicalsystem), intent(in) :: this      ! Type parent chemical system variable

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
integer                      :: &
 isps1, &
 isps2, &
 iph,   &
 i     
 integer,save:: ndim
 integer, dimension(:),pointer,save::ideqminsp
character(len=100)           :: &
 msg 
 logical,save::&
 haveideqminsp=.false.
!-------------------------------------------------------------------------
!
!   $code
!
 
iserror=.false.
msg=''
!%----------------------------------------------------------------
!% Check the number of species 
!%----------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%----------------------------------------------------------------
call check_pointer_ (cmineq,this%numsp,.true.)
!%----------------------------------------------------------------
call check_pointer_ (cmineq,this%numsp,.true.)
!%----------------------------------------------------------------
!% We get the index of the equilirium minerals
!%----------------------------------------------------------------
if (.not.haveideqminsp) then
    call get_eq_min_sp_index_(this,ideqminsp,ndim)
    haveideqminsp=.true.
endif

cmineq=0

!we save the values
do i=1,ndim
    cmineq(ideqminsp(i))=c(ideqminsp(i))
enddo


!%----------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Chemical System:'
print *,'Name:',this%name
print *,'Service: select_cmineq_'
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
subroutine compute_density_pchemsys &
    (this       ,&
    phasecode   ,&
    pliq        ,&
    temp        ,&
    c           ,&
    density     ,&
    iserror     ,&
    dc          ,&
    ddensitydc)


    

type(t_parentchemicalsystem), intent(in)      :: this         !Type parent chemical system. 

logical,intent(out)                     :: IsError

integer,intent(in)                      :: phasecode    !Phase for which the density will be calculated

real*8,intent(in)                       :: pliq         !liquid pressure

real*8,intent(in)                       :: temp         !temperature

real*8,pointer,dimension(:)             :: c         !concentration array

real*8,dimension(:,:),optional          :: dc  !concentration derivatives (needed to calculate density derivatives)

real*8,intent(out)                      :: density      !Density     

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
real*8,pointer,dimension(:) :: &
caq
real*8,pointer,dimension(:,:) :: &
dcaq

!-------------------------------------------------------------------------
!
!   $code
!


if (size(c)/=this%numsp) then
 msg='Error in number of species'
 goto 10
end if


if (phasecode.eq.ip_liquid_phase) then

    !%------------------------------------------------------------
    !% Get the positions of aqueous species
    !%------------------------------------------------------------
    call get_iposspsph_(this,this%aqphindex,isps1,isps2)
    !%------------------------------------------------------------
    call get_numsp_ (this%pphase(this%aqphindex)%ptr,nsp1)

    !% Fix c and dc for aqueous species only
    call check_pointer_ (caq,isps2-isps1+1,.true.)
    caq(isps1:isps2)=c(isps1:isps2)
    
    if (present(dc)) then
        call check_pointer_ (dcaq,isps2-isps1+1,size(dc,2),.true.)
        dcaq(isps1:isps2,1:size(dc,2))=dc(isps1:isps2,1:size(dc,2))
    endif

    if  (present(ddensitydc) ) then
        if (.not.present(dc)) then
            msg="dc are needed for calcutate density derivatives"
            isError=.true.
            goto 10
        endif

         call compute_density_(this%pphase(this%aqphindex)%ptr,pliq,temp,caq,density,.false.,iserror,dcaq,ddensitydc)
        if (isError) goto 10
    else
        call compute_density_(this%pphase(this%aqphindex)%ptr,pliq,temp,caq,density,.false.,iserror)
        if (isError) goto 10
    endif

else

    msg='Density function only supported for liquid phase'
    isError=.true.
    goto 10
endif

call check_pointer_ (caq,isps2-isps1+1,.false.)
if (present(dc)) call check_pointer_ (dcaq,isps2-isps1+1,size(dc,2),.false.)

!%----------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Parent Chemical system:'
print *,'Name:',this%name
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
subroutine compute_from_setre_pchemsys &
   (this, &
    ck,   & 
    setre, &
	nreact, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_parentchemicalsystem), intent(in):: this

integer, intent(in)                          :: nsp            ! Number of species 

real*8, intent(in), dimension(nsp)           :: ck             ! concentration in k.       

real*8, intent(inout), dimension(nsp)        :: setre          ! Changes due equilibrium reactions [mol/s] 

logical, intent(out)                         :: iserror        ! iserror=true, there was an error. 

Integer, intent(out)                         :: nreact

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
character(len=100), pointer :: &
 namesp(:) => null (), &
 namesp1(:) => null ()
integer                     :: &
 i, &
 j, &
 isps, &
 ireact, &
 ireact1, &
 isps1, &
 isps2
real*8, pointer             :: &
 b(:) => null ()
Integer, dimension(:), pointer :: idreact
!-------------------------------------------------------------------------
!
!   $code
!
 

!%---------------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------------
!% Check the number of species 
!%---------------------------------------------------------------
if (nsp/=this%numsp) then
 msg='Error in number of species' 
 goto 10 
end if 

call get_chem_info_ (this,msg,iserror,namesp=namesp)
!%---------------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------------
call check_pointer_ (b,nsp,.true.)
call check_pointer_ (namesp1,nsp,.true.)
call check_pointer_ (idreact,this%numreact,.true.)
if (iserror) goto 20
!%---------------------------------------------------------------
!% Determine and storage indices of equilibrium reactions 
!%---------------------------------------------------------------
nreact=0 
do i=1,this%numreact
 isps=this%idreactsp(i)
 if (.not.this%iskinreact(i).and.ck(isps)>0.0d0) then
  nreact=nreact+1
  idreact(nreact)=i
 end if
end do

!%---------------------------------------------------------------
!% If there aren't equilibrium reactions return
!%---------------------------------------------------------------
if (nreact>0) then 
!%---------------------------------------------------------------
!% Select the changes in the aqueous species 
!%---------------------------------------------------------------
   call get_iposspsph_ (this,this%aqphindex,isps1,isps2)
   isps=isps2-isps1+1
   b(1:isps)=setre(isps1:isps2)
   namesp1(1:isps)=namesp(isps1:isps2)
!%---------------------------------------------------------------
!% Select the equations corresponding to sorption complexes
!%---------------------------------------------------------------
   do i=1,this%numsurf 
      call get_iposspsurf_ (this,i,isps1,isps2)
      do j=isps1,isps2
	   isps=isps+1
       b(isps)=setre(j)
       namesp1(isps)=namesp(j)
	  end do
   end do
!%---------------------------------------------------------------
!% Not consider changes in the water species !!!!!
!%---------------------------------------------------------------
 where (namesp1=='h2o')
   namesp1=''
   b=0.0d0
 end where
!%----------------------------------------------------------------
!% Compute re from Set re
!%----------------------------------------------------------------
 call compute_r_from_stqtr_ &
   (this, &
    b, &
    namesp1(1:isps), &
    isps, &
    idreact(1:nreact), &
    nreact, &
    msg, &
    iserror)
 if (iserror) goto 20
!%----------------------------------------------------------------
!% Copy b to Setre
 Setre = 0.0
 do i=1,nreact
  Setre(idreact(i)) = b(i)
 enddo
 nreact=this%numreact
!%----------------------------------------------------------------
end if 
!%---------------------------------------------------------------
20 continue 
!%---------------------------------------------------------------
!% Deallocate local pointers 
!%---------------------------------------------------------------
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (namesp1,1,.false.)
call check_pointer_ (b,1,.false.)
call check_pointer_ (idreact,1,.false.)
if (iserror) goto 10
!%---------------------------------------------------------------
return
 
10 continue 
print *,'*********************************'
print *,'Chemical System:'
print *,'Name:', this%name
print *,'Service: specia_eqmin_from_setre_'
print *, msg
print *,'*********************************'
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
end module m_parentchemicalsystem
