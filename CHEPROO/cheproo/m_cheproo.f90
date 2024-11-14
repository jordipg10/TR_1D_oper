module m_cheproo
!-------------------------------------------------------------------------
!
!   $Description: The chemical processes are represented by a class called 
! CHEPROO (CHemical PRocesses Object-Oriented). 
! This class encapsulates two different classes associated between them. 
! On one hand the nodal chemistry class and on the other hand the chemical 
! system class. 
! Nodal chemistry class encapsulates geochemical state variables 
! (e.g. concentrations) and chemical system class encapsulates the 
! thermodynamic-kinetic models. 
! Both classes are associated by means of pointer link. 
! Thus, several objects of the nodal chemistry type could be associated to 
! one same chemical system object (see figure \ref{fig:cheproo}).
! Chemical system class offers functions mainly used by the nodal chemistry 
! class (e.g. chemical speciation functions).
!
!   $Use: m_chemicalsystem.
! m_nodalchemistry.
! flib_xpath.
! flib_sax.
! m_general_tools_cheproo. 
! m_constants
!
!   $Author: Sergio Andrés Bea Jofré. 
!
!   $License: UPC-CSIC, 2007. 
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_chemicalsystem
use m_nodalchemistry
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
private                        ::
!%-----------------------------------------------------------------
!% Public services
!%-----------------------------------------------------------------
public                         :: &
create_ &                    ! Create cheproo object 
,create_nodalchem_ &         ! Create nodal chemistries 
,create_rmrmtbox_ &          ! Create nodal chemistries for reactive multi-rate mass transfer box
,destroy_ &                  ! Destroy cheproo object  
,read_xml_ &                 ! Read the CHEPROO object from xml file. 
,read_txt_ &                 ! Read CHEPROO attributes from ascii file. 
,init_ &
,init_nodalchem_ &           ! Initialice nodal chemistries 
,init_rmrmtbox_ &            ! Initialice nodal chemistries (reactive multirate mass transfer) 
,update_r_sia_ith_ &         ! Update reaction term for SIA for ith nodal chemistry. 
,set_  &                     ! Set general attributes in the CHEPROO object. 
,set_from_uaq_ith_  & 
,set_from_cpri_ith_  &       ! Set ith nodal chemistry from cpri
,set_from_cpri_ &            ! Set all nodal chemistry objects from cpri
,save_nodalchem_  &
,get_chem_info_ &            ! Return general chemical information for all the nodal chemistry
,get_chem_info_ith_ &        ! Return general chemical information for a single nodal chemistry 
,get_cpri_ith_ &             ! Get primary species from the ith nodalchemistry
,get_cpri_water_ &           ! Get primary species from a water for the component definition of some nodalchemistry
,get_caq_water_ &            ! Return the concentration for all the aqueous species defined by the water "ithw"
,get_cgas_water_ &           ! Return the concentration for all the gaseous species defined by the water "ithw"
,get_ipos_unk_ &
,get_umob_water_ &           ! Return mobile components ofr ith water. 
,get_uaq_nch_sia_ &          ! return the aqueous components for all nodal chemistry objects 
,get_uaq_ &                  ! Returns the uaq=Uaqith*caqjth 
,get_uads_ &                 ! Returns the uads=Uadsith*cadsjth 
,get_uaq_water_ &            ! Return the aqueous components for ith water defined in CHEPROO
,get_uaq_nch_ &              ! Return the aqueous components for ith nodal chemistry 
,get_ph_nch_ &               ! Return the ph for ith nodal chemistry 
,get_numsp_ &                ! Return the number of species
,get_namesp_ &               ! Return the name of species
,get_name_comp_ith_  &       ! Return the name of component for the it nodal chemistry
,get_aq_sps_pos_ &           ! Return the first and last indices of the aqueous spcecies. 
,get_ncomp_dsa_ &            ! Return the number of components for each nodal chemistry object and their corresponding hash index.       
,get_c_ithnch_ &             ! Return the concentrations of ith nodal chemistry objectname of species
,get_ujth_cmobith_ &         ! Return mobile components of ith nodal chemistry 
,get_ujth_cadsith_ &         ! Return adsorbed omponents of ith nodal chemistry 
,get_iumob_ith_ &
,get_iumob_water_ &
,get_usktrk_ith_ &           ! Return the change in components for ith nodal chemistry  
,get_dumob_ith_ &            ! Return derivate of mobile components. 
,get_duaq_ &                 ! Return duaq[npri]=Uaqith*dcaqjth, where ith and jth are referring to nodal chemistry index for U aqand dcaq respectively.   
,get_dusktrk_ &              ! Return dusktrk[ncomp]=theta*U*dusktrk for ith nodal chemistry object 
,get_duads_ &                ! Return duads[npri]=Uith*duadsjth
,get_umin_  &                ! Return umin[npri]=Uith*cminjth
,get_umineq_  &              ! Return umineq[npri]=Uith*cmineqjth
,get_duads_ith_ & 
,get_numunk_sia_ &           ! Return the number of unknowns for SIA
,get_numunk_dsa_ &          
,get_mxnumunk_dsa_ &         ! Return maximum number of primary species defined in cheproo object. 
,mix_ &                      ! Mixing waters
,solve_react_trp_step_ &     ! Solve one reactive transport step (for SIA or DSA)
,build_indep_term_dsa_ &     ! Build independent term for DSA
,build_indep_term_sia_ &     ! Build independent term for SIA
,build_jacob_resid_dsa_ &    ! Build the jacobian and residual for DSA (Direct Substitution Approach)
,update_and_check_ &         ! Update and check covergence for SIA and DSA 
,check_compzone_ &           ! Check components zone in all nodal chemistry objects 
,compute_total_mol_ &        ! Compute the total mol of components in all nodal chemistries defined in CHEPROO
,compute_f_df_rmrmt_ &       ! Compute f and df for reactive multi-rate mass transfer box 
,compute_eq_min_from_trp_ &  ! Compute concentrations of equilibrium minerals building Sere
,compute_eq_min_from_SeRe_ & ! Compute concentrations of equilibrium minerals processing Sere
,compute_aqdensity_ith_ &    ! Compute aqueous density (and derivatives if asked) for the ith nodal chemistry
,equilibrate_ith_ &	         ! equilibrate ith nodal chemistry with ith mineral set.
,reaction_path_ith_ &        ! Reaction path operations on nodal chemistry obejcts. 
,copy_compdef_ &             ! Copy the component definition from one of the nodal chemistry group to the other
!%,operator(*) &             ! Multiply the CHEPROO object by a scalar.  
,assignment(=) &             ! Copy a cheproo object in other cheproo object. 
,update_ &                   ! Update different chemical/physical parameters (e.g. porosity). 
,write_backup_ &             ! Write backup file. 
,write_signature_ &          ! Write CHEPROO signature in an output file. 
,write_ &                    ! Write in ascii the attributes encapsulated in CHEPROO. 
,test_
!%-----------------------------------------------------------------
!% Private services
!%-----------------------------------------------------------------
private :: &
make_chem_info_sia_ &
,make_preview_info_dsa_ &
,make_preview_info_dsa_cheproo &
,make_zonal_chem_info_ithnch_ &
,get_jacobpos_band_ &
,build_trp_eq_sia_ & 
,build_jacob_resid_dsa_ithnch_ &
,build_trp_eq_sia_ithnch_ &
,solve_react_trp_step_sia_ &
,solve_react_trp_step_dsa_ &        
,add_jacobian_dsa_ &
,add_f_df_dsa_ith_ &                ! Add f and df in the jacobian and residula of DSA
,compute_in_out_mol_ &              ! Compute the input/output mol of components when solving reactive transport step. This service is performed by mol balance
,begin_element_handler &
,read_xml_loc &
,compute_kin_mol_ &
,compute_f_df_rmrmt_cheproo &
,init_nodalchem_cheproo &
,init_nodalchem_from_xml_cheproo &
,init_nodalchem_from_output_cheproo &
,init_rmrmtbox_cheproo &
,get_umob_water1_cheproo &
,get_umob_water2_cheproo &
,init_from_backup_cheproo &
,write_backup_cheproo
!%--------------------------------------------------------
!% Private parameters
!%--------------------------------------------------------
logical, parameter :: iswchsys=.true.  ! If .true. the chemical system object will be written in an outtpt file 

integer, parameter :: sniasolver=1     ! Use SNIA approach for to solve reactive transport step 

integer, parameter :: siasolver=2      ! Use SIA approach for to solve reactive transport step 

integer, parameter :: dsasolver=3      ! Use DSA approach for to solve reactive transport step 

real*8, parameter  :: zerosia=1.0d-15  ! Zero for unknowns in SIA solver

real*8, parameter  :: zerodsa=1.0d-100 ! Zero for unknowns in DSA solver
!%-----------------------------------------------------------------
!% Type mineral set
!%-----------------------------------------------------------------
type, private:: t_minset

character(len=100)                        :: name        ! Name of mineral set 

character(len=100), pointer, dimension(:) :: namemin     ! Name of minerals in set 

real*8, pointer, dimension(:)             :: cmin        ! Concentration of minerals

character(len=100), pointer, dimension(:) :: unitcmin    ! Unit for mineral concentrations 

real*8, pointer, dimension(:)             :: areamin     ! Initial area for minerals 

character(len=100), pointer, dimension(:) :: unitareamin ! Unit of mineral surface area

integer                                   :: nummin      ! Number of minerals 

end type
!%-----------------------------------------------------------------
!% Type gas set
!%-----------------------------------------------------------------
type, private:: t_gasset

character(len=100)                        :: name     ! Name of gas set 

character(len=100), pointer, dimension(:) :: namegas  ! Name of gas species

character(len=100), pointer, dimension(:) :: unitcgas ! Unit for concentrations

real*8, pointer, dimension(:)             :: cgas     ! Partial preassure of gas 

integer                                   :: numgas   ! Number of gas species 

end type
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Type surface
!%-----------------------------------------------------------------
type, private:: t_surf

character(len=100)            :: name           ! Name of surface 

logical                       :: isequilibrate  ! If true, the solution is equilibrated

real*8, pointer, dimension(:) :: txoh           ! Concentrations of sites in the surface [numsite]

real*8, pointer, dimension(:) :: capint         ! Capacitance 1

real*8, pointer, dimension(:) :: capext         ! Capacitance 2

real*8, pointer, dimension(:) :: spsurfarea     ! Specific surface area [numsite]

integer                       :: numsites       ! Number of sites in surface

character(len=100)            :: namemin        ! Name of mineral associated to surface   

end type
!%-----------------------------------------------------------------
!% Type surface set
!%-----------------------------------------------------------------
type, private:: t_surfset

character(len=100)                  :: name    ! Name of surface set 

type(t_surf), pointer, dimension(:) :: surf    ! t_surf pointer 

integer                             :: numsurf ! Number of surfaces 

end type
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Type vector 
!%-----------------------------------------------------------------
type, private:: t_vector

real*8, pointer, dimension(:) :: vector   ! Vector

integer                       :: nrow     ! Number of rows 

end type
!%-----------------------------------------------------------------
!% Type array 
!%-----------------------------------------------------------------
type, private:: t_array

real*8, pointer, dimension(:,:) :: array   ! Array[nrow,ncol]

integer                         :: nrow    ! Number of rows

integer                         :: ncol    ! Number of columns 

end type
!%-----------------------------------------------------------------
!% Type output
!%-----------------------------------------------------------------
type, private:: t_output

character(len=100), pointer, dimension(:) :: namecomp   ! Name of components

character(len=100), pointer, dimension(:) :: namesp     ! Name of species 

character(len=100), pointer, dimension(:) :: unitsp     ! Unit for concentrations 

integer                                   :: numcomp    ! Number of components

integer                                   :: numsp      ! Number of species

end type
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!% Type pointer to CHEPROO object 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
type, public:: t_pcheproo

type(t_cheproo), pointer:: ptr  

end type
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Type definition
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
type, public::t_cheproo
 
private                                        ::
 
character(len=100)                             :: name          ! Name of CHEPROO
 
type(t_pnodalchemistry), pointer, dimension(:) :: pnodalchem0   ! Initial Nodal Chemistry list

type(t_pnodalchemistry), pointer, dimension(:) :: pnodalchemk   ! Nodal Chemistry list in k

type(t_pnodalchemistry), pointer, dimension(:) :: pnodalchemk1  ! Nodal Chemistry list in k+1

type(t_pnodalchemistry), pointer, dimension(:) :: pwater        ! List of nodal chemistry pointers (waters)
 
type(t_pnodalchemistry), pointer, dimension(:) :: prmrmtboxk1   ! List of nodal chemistry objects for reactive multi-rate mass transfer (in k+1)

type(t_pnodalchemistry), pointer, dimension(:) :: prmrmtboxk    ! List of nodal chemistry objects for reactive multi-rate mass transfer (in k) 
 
type(t_pchemsys), pointer, dimension(:)        :: pchemsys      ! List of chemical system pointers
 
type(t_minset), pointer, dimension(:)          :: minset        ! Private mineral set information
 
type(t_gasset), pointer, dimension(:)          :: gasset        ! Private gas set information
 
type(t_surfset), pointer, dimension(:)         :: surfset       ! Private surface set information
 
type(t_output)                                 :: output        ! Private output set information
 
integer                                        :: numnch        ! Number of nodal chemistry objects

integer                                        :: numw          ! Number of water objects

integer                                        :: numrmrmtbox   ! Number of nodal chemistry objects for reactive multi-rate mass transfer 

integer                                        :: numchemsys    ! Number of chemical system objects 

integer                                        :: numminset     ! Number of mineral sets

integer                                        :: numgasset     ! Number of gas sets

integer                                        :: numsurfset    ! Number of surfaces sets
 
logical                                        :: isderstored   ! If the derivates must be stored in the nodal chemistries

logical                                        :: iswriteinfo   ! If the CHEPROO information must be write

logical, pointer, dimension(:)                 :: islockchemsys ! [numchemsys]
 
integer                                        :: itypesolver   ! 1=SNIA, 2=SIA, 3=DSA
 
real*8                                         :: zero          ! zero for unknowns
 
end type t_cheproo
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!% Public interfaces
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
interface create_
 
module procedure create_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface create_nodalchem_
 
module procedure create_nodalchem_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface create_rmrmtbox_
 
module procedure create_rmrmtbox_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_f_df_rmrmt_ 
 
module procedure compute_f_df_rmrmt_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface destroy_
 
module procedure destroy_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface save_nodalchem_
 
module procedure save_nodalchem_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_xml_
 
module procedure read_xml_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_txt_
 
module procedure read_txt_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface init_nodalchem_
 
module procedure init_nodalchem_cheproo
module procedure init_nodalchem_from_xml_cheproo
module procedure init_nodalchem_from_output_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface init_
 
module procedure init_from_backup_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface init_rmrmtbox_
 
module procedure init_rmrmtbox_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_
 
module procedure set_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_from_uaq_ith_
 
module procedure set_from_uaq_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_from_cpri_ith_
 
module procedure set_from_cpri_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_from_cpri_
 
module procedure set_from_cpri_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_cpri_ith_
 
module procedure get_cpri_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_caq_water_
 
module procedure get_caq_water_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_ipos_unk_
 
module procedure get_ipos_unk_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_mxnumunk_dsa_
 
module procedure get_mxnumunk_dsa_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_umob_water_
 
module procedure get_umob_water1_cheproo
module procedure get_umob_water2_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_uaq_water_
 
module procedure get_uaq_water_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_ph_nch_
 
module procedure get_ph_nch_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_uaq_nch_
 
module procedure get_uaq_nch_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_uaq_
 
module procedure get_uaq_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_uads_
 
module procedure get_uads_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_duaq_
 
module procedure get_duaq_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_duads_
 
module procedure get_duads_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_usktrk_ith_
 
module procedure get_usktrk_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_dusktrk_
 
module procedure get_dusktrk_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_uaq_nch_sia_
 
module procedure get_uaq_nch_sia_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_namesp_
 
module procedure get_namesp_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_aq_sps_pos_
 
module procedure get_aq_sps_pos_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_c_ithnch_
 
module procedure get_c_ithnch_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_numsp_
 
module procedure get_numsp_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_ujth_cmobith_
 
module procedure get_ujth_cmobith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_dumob_ith_jth_
 
module procedure get_dumob_ith_jth_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_ujth_cadsith_
 
module procedure get_ujth_cadsith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_dumob_ith_
 
module procedure get_dumob_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_duads_ith_
 
module procedure get_duads_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_chem_info_
 
module procedure get_chem_info_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_chem_info_ith_
 
module procedure get_chem_info_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_numunk_sia_
 
module procedure get_numunk_sia_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_numunk_dsa_
 
module procedure get_numunk_dsa_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_iumob_ith_
 
module procedure get_iumob_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_iumob_water_
 
module procedure get_iumob_water_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface mix_
 
module procedure mix_nodalchem_cheproo
module procedure mix_desimoni_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface update_r_sia_ith_
 
module procedure update_r_sia_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface update_and_check_
 
module procedure update_and_check_convergence_dsa_cheproo
module procedure update_r_and_check_convergence_sia_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface solve_react_trp_step_
 
module procedure solve_react_trp_step_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_jacob_resid_dsa_
 
module procedure build_jacob_resid_dsa_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_indep_term_dsa_
 
module procedure build_indep_term_dsa_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_indep_term_sia_
 
module procedure build_indep_term_sia_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_total_mol_
 
module procedure compute_total_mol_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface check_compzone_
 
module procedure check_compzone_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface reaction_path_ith_
 
module procedure reaction_path_ith_cheproo
module procedure reaction_path_kinetic_ith_cheproo
module procedure reaction_path_ev_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface equilibrate_ith_
 
module procedure equilibrate_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface assignment (=)
 
module procedure copy_cheproo
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface operator (*)
 
module procedure scalarbycheproo_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface update_
 
module procedure update_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface write_
 
module procedure write_output1_cheproo
module procedure write_output2_cheproo
module procedure write_attributes_cheproo
module procedure write_ith_nch1_cheproo
module procedure write_ith_nch2_cheproo
module procedure write_from_namesp1_cheproo
module procedure write_from_namesp2_cheproo
module procedure write_mol_balance_cheproo
module procedure write_output_tecplot_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface write_backup_
 
module procedure write_backup_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface write_signature_
 
module procedure write_signature_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_ncomp_dsa_
module procedure get_ncomp_dsa_cheproo
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface copy_compdef_
module procedure copy_compdef_cheproo
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_cpri_water_
module procedure get_cpri_water_cheproo
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_umin_
module procedure get_umin_cheproo
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_umineq_
module procedure get_umineq_cheproo
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_name_comp_ith_
module procedure get_name_comp_ith_cheproo
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------


!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%----------------Private Interfaces----------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
interface make_chem_info_sia_
 
module procedure make_chem_info_sia_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface make_preview_info_dsa_
 
module procedure make_preview_info_dsa_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface make_zonal_chem_info_ithnch_
 
module procedure make_zonal_chem_info_ithnch_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_jacobpos_band_
 
module procedure get_jacobpos_band_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_trp_eq_sia_
 
module procedure build_trp_eq_sia_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_jacob_resid_dsa_ithnch_
 
module procedure build_jacob_resid_dsa_ithnch_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface build_trp_eq_sia_ithnch_
 
module procedure build_trp_eq_sia_ithnch_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_eq_min_from_trp_
 
module procedure compute_eq_min_from_trp_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_eq_min_from_SeRe_
 
module procedure compute_eq_min_from_SeRe_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface solve_react_trp_step_sia_
 
module procedure solve_react_trp_step_sia_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface solve_react_trp_step_dsa_
 
module procedure solve_react_trp_step_dsa_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface add_jacobian_dsa_
 
module procedure add_jacobian_dsa_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface add_f_df_dsa_ith_
 
module procedure add_f_df_dsa_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_in_out_mol_
 
module procedure compute_in_out_mol_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_kin_mol_
 
module procedure compute_kin_mol_cheproo
 
end interface

!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_aqdensity_ith_
 
module procedure compute_aqdensity_ith_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface test_
 
module procedure test_cheproo
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface begin_element_handler
 
module procedure begin_element_handler
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_xml_loc
 
module procedure read_xml_loc
 
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
subroutine create_cheproo &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create CHEPROO object.
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout):: this    ! Type CHEPROO object 
 
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
 i 
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
!% Write signature 
!%-----------------------------------------------------------
call write_signature_ (this,6)
!%----------------------------------------------------------
this%pnodalchem0 => null ()
this%pnodalchemk => null ()
this%pnodalchemk1 => null ()
this%prmrmtboxk1 => null ()
this%prmrmtboxk => null ()
this%pwater => null ()
this%minset => null ()
this%gasset => null ()
this%surfset => null ()
this%pchemsys => null ()
!%-----------------------------------------------------------
this%output%namecomp => null ()
this%output%namesp => null ()
this%output%unitsp => null ()
this%output%numcomp=0
this%output%numsp=0
!%-----------------------------------------------------------
this%name=' '
this%numw=0
this%numrmrmtbox=0 
this%numminset=0
this%numgasset=0
this%numsurfset=0
this%numnch = 0
this%isderstored=.false.
this%itypesolver=0
this%iswriteinfo=.false.
this%islockchemsys => null ()
this%zero=0.0d0
!%-----------------------------------------------------------
this%numchemsys=0
this%pchemsys => null ()
!%-----------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_cheproo &
   (this)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy CHEPROO object
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)        :: this 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 i, &
 j 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%--------------------------------------------------------------
this%name=' '
!%--------------------------------------------------------------
this%isderstored=.false.
!%--------------------------------------------------------------
! Destroy and deallocate nodal chemistry
!%--------------------------------------------------------------
if (this%numnch>0) then
 do i=1,this%numnch
  call destroy_ (this%pnodalchem0(i)%ptr)
  call destroy_ (this%pnodalchemk(i)%ptr)
  call destroy_ (this%pnodalchemk1(i)%ptr)
  deallocate (this%pnodalchem0(i)%ptr)
  deallocate (this%pnodalchemk(i)%ptr)
  deallocate (this%pnodalchemk1(i)%ptr)
 end do
 deallocate (this%pnodalchem0)
 deallocate (this%pnodalchemk)
 deallocate (this%pnodalchemk1)
 this%pnodalchem0 => null ()
 this%pnodalchemk => null ()
 this%pnodalchemk1 => null ()
 this%numnch=0
end if
!%--------------------------------------------------------------
! Destroy and deallocate nodal chemistry (RMRMT)
!%--------------------------------------------------------------
if (this%numrmrmtbox>0) then
 do i=1,this%numrmrmtbox
  call destroy_ (this%prmrmtboxk1(i)%ptr)
  call destroy_ (this%prmrmtboxk(i)%ptr)
  deallocate (this%prmrmtboxk1(i)%ptr)
  deallocate (this%prmrmtboxk(i)%ptr)
 end do
 deallocate (this%prmrmtboxk1)
 deallocate (this%prmrmtboxk)
 this%prmrmtboxk1 => null ()
 this%prmrmtboxk => null ()
 this%numrmrmtbox=0
end if
!%--------------------------------------------------------------
! Destroy and deallocate nodal chemistry (water)
!%--------------------------------------------------------------
if (this%numw>0) then
 do i=1,this%numw
  call destroy_ (this%pwater(i)%ptr)
  deallocate (this%pwater(i)%ptr)
 end do
 deallocate (this%pwater)
 this%pwater => null ()
 this%numw=0
end if
!%--------------------------------------------------------------
! Destroy mineral sets 
!%--------------------------------------------------------------
if (this%numminset>0) then
 do i=1,this%numminset
  this%minset(i)%name=''
  if (this%minset(i)%nummin>0) then
   call check_pointer_ (this%minset(i)%namemin,1,.false.)
   call check_pointer_ (this%minset(i)%unitcmin,1,.false.)
   call check_pointer_ (this%minset(i)%unitareamin,1,.false.)
   call check_pointer_ (this%minset(i)%cmin,1,.false.)
   call check_pointer_ (this%minset(i)%areamin,1,.false.)
  end if
 end do
 deallocate (this%minset)
end if
this%numminset=0
this%minset => null()
!%--------------------------------------------------------------
! Destroy gas sets 
!%--------------------------------------------------------------
if (this%numgasset>0) then
 do i=1,this%numgasset
  this%gasset(i)%name=''
  if (this%gasset(i)%numgas>0) then
   call check_pointer_ (this%gasset(i)%namegas,1,.false.)
   call check_pointer_ (this%gasset(i)%unitcgas,1,.false.)
   call check_pointer_ (this%gasset(i)%cgas,1,.false.)
  end if
 end do
 deallocate (this%gasset)
end if
this%numgasset=0
this%gasset => null()
!%--------------------------------------------------------------
! Destroy surface sets 
!%--------------------------------------------------------------
if (this%numsurfset>0) then
 do i=1,this%numsurfset
  this%surfset(i)%name=''
  do j=1,this%surfset(i)%numsurf
   call check_pointer_ (this%surfset(i)%surf(j)%txoh,1,.false.)
  end do
  if (this%surfset(i)%numsurf>0) then
   deallocate(this%surfset(i)%surf)
  end if
 end do
 deallocate (this%surfset)
end if
this%numsurfset=0
this%surfset => null()
!%--------------------------------------------------------------
call check_pointer_ (this%output%namecomp,1,.false.)
call check_pointer_ (this%output%namesp,1,.false.)
call check_pointer_ (this%output%unitsp,1,.false.)
this%output%numcomp=0
this%output%numsp=0
!%--------------------------------------------------------------
if (this%numchemsys>0) then
    do i=1,this%numchemsys
      if (this%islockchemsys(i)) then  
		call destroy_ (this%pchemsys(i)%ptr)
        deallocate (this%pchemsys(i)%ptr)
	  end if 
	  this%pchemsys(i)%ptr => null ()
    end do
    deallocate(this%pchemsys)
    this%pchemsys => null ()
    this%numchemsys=0
end if 
!%--------------------------------------------------------------
this%itypesolver=0
!%--------------------------------------------------------------
this%zero=0.0d0
!%--------------------------------------------------------------
this%iswriteinfo=.false.
!%--------------------------------------------------------------
this%islockchemsys=.false.
!%--------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine test_cheproo &
   (this, &
    ntime, &
    iproof)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: General test for CHEPROO. 
! At the moment for test allocation of memory.
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout):: this    ! Type CHEPROO variable

integer, intent(in)            :: ntime   ! Times

integer, intent(in)            :: iproof  ! iproof=1
                                          ! iproof=0 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer                :: &
 c(:) => null(), &
 cmob(:,:)=> null()
integer                        :: &
 i, &
 j
type(t_nodalchemistry), pointer:: &
 nch(:)
logical iserror 
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------
if (iproof==1) then
 if (this%numnch>0) then
 allocate (nch(this%numnch))
 do i=1,this%numnch
  call create_ (nch(i))
 end do
  do j=1,ntime
   do i=1,this%numnch
    nch(i)=this%pnodalchemk1(i)%ptr
   end do
  end do
 
  do i=1,this%numnch
    call destroy_ (this%pnodalchemk1(i)%ptr)
  end do
 
 deallocate (nch)
 
 end if
end if
!%----------------------------------------------------------
if (iproof==2) then
 do j=1,ntime
  do i=1,this%numnch
   call get_chem_info_ (this%pnodalchemk1(i)%ptr,iserror,c=c,cmob=cmob)
  end do
 end do
 call check_pointer_ (c,1,.false.)
 call check_pointer_ (cmob,1,1,.false.)
end if
!%-----------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_cheproo &
   (this, &
    namefile, &
    iserror)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read the CHEPROO object from xml file. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout):: this      ! Type CHEPROO variable

character(len=*), intent(in)   :: namefile  ! Name of the file 

logical, intent(out)           :: iserror   ! iserror=true, then there was an error
 
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
 name
type(dictionary_t) :: &
 attributes
integer:: &
 iostat
type(xml_t):: &
 fxml
character(len=100) :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
! Open xml file 
!%----------------------------------------------------------
call open_xmlfile(namefile, fxml, iostat)
if (iostat/=0) then
 msg='Error when open file'
 call add_ (msg,namefile)
 goto 10
end if
!%-----------------------------------------------------------
write(6,*) '---------------------------------------------------'
write(6,*) '=======> Reading CHEPROO object'
!%-----------------------------------------------------------
call xml_parse(fxml,begin_element_handler = begin_element_handler)
!%-----------------------------------------------------------
! End and close xml file
!%-----------------------------------------------------------
call endfile_xmlfile (fxml)
call close_xmlfile(fxml)
!%-----------------------------------------------------------
! Read information about nodal chemistry objects 
! (waters and mineral, gas and surfaces sets)
!%-----------------------------------------------------------
call read_xml_loc(name,attributes,this=this,msg=msg,iserror=iserror)
if (iserror) goto 10
!%-----------------------------------------------------------
if (iserror) goto 10
!%-----------------------------------------------------------
write(6,*) '=======> Reading CHEPROO object finished'
write(6,*) '---------------------------------------------------'
!%-----------------------------------------------------------
return
 
10 continue 
print *,'*******************'
print *,'CHEPROO:  '
print *,'Name:',this%name
print *,'Service: read_xml_ '
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
subroutine read_txt_cheproo &
   (this, &
    namefile, &
    namethdb, &
    namekindb, &
    filebase, &
    itypesolver, &
    ioptxhl, &
    iserror, &
    zero)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read CHEPROO attributes from ascii file. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout) :: this         ! Type CHEPROO variable 

character(len=*), intent(in)    :: namefile     ! Name of CHEPROO file 

character(len=*), intent(in)    :: namethdb     ! Name of thermodynamic data base

character(len=*), intent(in)    :: namekindb    ! Name of kinetic data base

character(len=*), intent(in)    :: filebase     ! Path 

integer, intent(in)             :: itypesolver  ! Indice pf type solver 

integer, intent(in)             :: ioptxhl      ! Option if to compute ioptxhl 

logical, intent(out)            :: iserror      ! iserror=true, then there was an error 

real*8, intent(in), optional    :: zero         ! zero for concentrations 
 
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
 msg, &
 name, &
 name1, &
 name2, &
 inamecomp, &
 iconstr, &
 label, &
 namemin, &
 nametask 
integer                        :: &
 lbase, &
 niwtype, &
 nbwtype, &
 nrwtype, &
 ncomp, &
 iicon, &
 ilabel, &
 isps, &
 iwater, &
 i, &
 j, &
 nwater, &
 itypechemsys, &
 nsurf, &
 iequil, &
 ntask, &
 nmix, &
 iminset 
logical                        :: &
 isconvergence, &
 isboundary, &
 isrecharge, &
 isbe, &
 havezero, &
 isspecia
real*8                         :: &
 temp, &
 icguess, &
 ictot, &
 cputime, &
 prop, &
 conc 
integer, parameter             :: &
 mxdim=100 
real*8, pointer                :: &
 value1(:) => null (), &
 value2(:) => null ()
integer, pointer               :: &
 icon(:) => null ()
character(len=100), pointer    :: &
 namesps(:) => null (), &
 unit(:) => null (), &
 constraint(:) => null ()
character(len=100), parameter  :: &
 label1='INITIAL AND BOUNDARY WATER TYPES', &
 label2='INITIAL MINERAL ZONES', &
 label3='INITIAL SURFACE ADSORPTION ZONES', &
 label4='GAS BOUNDARY ZONES', &
 label5='OUTPUT', &
 label6='TASK'
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Check optional arguments
!%-----------------------------------------------------------
havezero=present(zero)
!%-----------------------------------------------------------
!% Check the reactive transport solver and assign zero 
!% concentrations 
!%-----------------------------------------------------------
select case (itypesolver)
case (dsasolver)
 this%itypesolver=dsasolver
 itypechemsys=2
 if (havezero) then 
    this%zero=zero
 else
    this%zero=zerodsa
 end if
case (siasolver)
 this%itypesolver=siasolver
 itypechemsys=1
 if (havezero) then 
    this%zero=zero
 else
    this%zero=zerosia
 end if
case (sniasolver)
 this%itypesolver=sniasolver
 itypechemsys=1 
 if (havezero) then 
    this%zero=zero
 else
    this%zero=zerosia
 end if
case default 
 msg='Error, not recognized reactive transport solver:'
 call add_ (msg,itypesolver)
 goto 10 
end select  
!%-----------------------------------------------------------
!% Assign the name to CHEPROO object
!%-----------------------------------------------------------
this%name='CHEPROO object'
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
write(6,*) '---------------------------------------------------'
write(6,*) '=======> Reading CHEPROO object'
!%-----------------------------------------------------------
!% Read the chemical system object 
!%-----------------------------------------------------------
this%numchemsys=1
allocate (this%pchemsys(this%numchemsys))
call check_pointer_ (this%islockchemsys,this%numchemsys,.true.)
do i=1,this%numchemsys
 this%islockchemsys(i)=.true. 
 allocate (this%pchemsys(i)%ptr)
 call create_ (this%pchemsys(i)%ptr)
end do
call read_txt_ (this%pchemsys(1)%ptr,namefile,namethdb,namekindb, &
                filebase,itypechemsys,ioptxhl,iserror)
if (iserror) then
 msg='Error when calling read_txt_'
 goto 10
end if
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!% Read waters 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!% Open unit 
!%-----------------------------------------------------------
name=namefile
call lastletter_ (lbase,filebase)
name=filebase(1:lbase)//name
open(unit=1,file=name,status='old',err=30)
!%-----------------------------------------------------------
!% Find label corresponding to waters definition 
!%-----------------------------------------------------------
call find_label_ (label1,1,isbe)
if (.not.isbe) then
  msg='Error, not found the label INITIAL AND BOUNDARY WATER TYPES'
  goto 10
end if 
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (namesps,mxdim,.true.)
call check_pointer_ (icon,mxdim,.true.)
call check_pointer_ (unit,mxdim,.true.)
call check_pointer_ (value1,mxdim,.true.)
call check_pointer_ (value2,mxdim,.true.)
call check_pointer_ (constraint,mxdim,.true.)
!%-----------------------------------------------------------
!% Initialice local variables 
!%-----------------------------------------------------------
unit='m'
niwtype=0
nbwtype=0
nrwtype=0 
isboundary=.true.
isrecharge=.true. 
!%-----------------------------------------------------------
!% Read number of waters 
!%-----------------------------------------------------------
read (1,*,err=40) niwtype, nbwtype,nrwtype
this%numw=niwtype+nbwtype+nrwtype
if (this%numw>0) then
 print *,'==========> Starting speciation of water types'
 allocate (this%pwater(this%numw))
  nwater=niwtype
  iwater=0
!%-----------------------------------------------------------
!% Read water definition 
!%-----------------------------------------------------------
5 do i=1,nwater
     ncomp=0 
     iwater=iwater+1
     allocate (this%pwater(iwater)%ptr)
     call create_ (this%pwater(iwater)%ptr)
     call set_(this%pwater(iwater)%ptr,iserror, &
               chemsys=this%pchemsys(1)%ptr,isderstored=.true.)
     if (iserror) goto 10
     name=' '
     read (1,*,err=40) ilabel,temp,name
     if (name==' ') then
      name='Water'
     end if  
     call add_ (name,iwater)
     call set_(this%pwater(iwater)%ptr,iserror,name=name)
     if (iserror) goto 10 
     if (i/=ilabel) goto 40
     read (1,*)
     do 
        read (1,*) inamecomp,iicon,icguess,ictot,iconstr
        if (inamecomp=='*'.or.ncomp>=mxdim) exit 
        ncomp=ncomp+1
        namesps(ncomp)=inamecomp 
        icon(ncomp)=iicon
        value1(ncomp)=icguess
        value2(ncomp)=ictot
        constraint(ncomp)=iconstr 
     end do 
     print *,'=======> Starting speciation',iwater  
     call set_from_solution_type_ &
            (this%pwater(iwater)%ptr, &
             namesps(1:ncomp), &
             value2(1:ncomp), &
             value1(1:ncomp), &
             unit(1:ncomp), &
             icon(1:ncomp), &
             constraint(1:ncomp), &
             ncomp, &
             temp, &
             isconvergence, & 
			 cputime, &
             iserror)
     if (iserror) goto 20 
     if (.not.isconvergence) then 
        msg='Convergence problems when specia water:'
        call add_ (msg,i)
		iserror=.true. 
        goto 20 
     end if  
     call write_ (this%pwater(iwater)%ptr,6,iserror)
     if (iserror) goto 20 
	 print *,'=======> Finishing speciation water',iwater 
 end do
 if (isboundary) then
  nwater=nbwtype
  isboundary=.false.
  goto 5 
 end if 
 if (isrecharge) then
  nwater=nrwtype
  isrecharge=.false.
  goto 5 
 end if  
 write(6,*) '==========> Speciation of water types finished'
else
 msg='Error, the number of waters is 0'
 iserror=.true. 
 goto 20 
end if 
!%-----------------------------------------------------------
!% Find label corresponding to initial mineral zones
!%-----------------------------------------------------------
call find_label_ (label2,1,isbe)
if (isbe) then 
!%-----------------------------------------------------------
!% Read the initial mineral zones 
!%-----------------------------------------------------------
read (1,*) this%numminset
if(this%numminset>0) then
  write(6,*) '==========> Reading mineral zones'
  allocate(this%minset(this%numminset))
  do i=1,this%numminset
    read (1,*) ilabel
    if (i/=ilabel) goto 50
    write(6,*) '====> Reading mineral zone:',ilabel
    read (1,*)
    isps=0 
    name='mineral set'
    call add_ (name,i) 
    this%minset(i)%name=name
    do 
     isps=isps+1
     read (1,*,err=60) namesps(isps),value1(isps),value2(isps)
     if (namesps(isps)=='*') then
      isps=isps-1
      exit
     end if 
    end do
    this%minset(i)%nummin=isps
    if (isps>0) then
     allocate(this%minset(i)%namemin(isps))
     allocate(this%minset(i)%cmin(isps))
     allocate(this%minset(i)%areamin(isps))
     allocate(this%minset(i)%unitcmin(isps))
     allocate(this%minset(i)%unitareamin(isps))
     this%minset(i)%namemin=namesps(1:isps)
     this%minset(i)%cmin=value1(1:isps)
     this%minset(i)%areamin=value2(1:isps)
     this%minset(i)%unitcmin='m3/m3rock'
     this%minset(i)%unitareamin='m2/m3rock'
    end if 
    do j=1,this%minset(i)%nummin
     name=this%minset(i)%namemin(j)
     call get_sp_index_ (this%pchemsys(1)%ptr,name,isps)
     if (isps==0) then
      msg='Error, not defined in the chemical system the species:'
      call add_ (msg,name)
      iserror=.true. 
      goto 20 
     end if 
    end do 
  end do  
  write(6,*) '==========> Reading mineral zones finished'
end if
end if
!%-----------------------------------------------------------
!% Find label corresponding to surface zones
!%-----------------------------------------------------------
call find_label_ (label3,1,isbe)
if (isbe) then
!%-----------------------------------------------------------
!% Get the total numbers of surfaces 
!%-----------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numsurf=nsurf)
if (iserror) goto 20 
!%-----------------------------------------------------------
!% Read the initial output 
!%-----------------------------------------------------------
read (1,*) this%numsurfset
if(this%numsurfset>0) then
  write(6,*) '==========> Reading surface zones'
  allocate(this%surfset(this%numsurfset))
  do i=1,this%numsurfset
    read (1,*,err=70) ilabel,iequil
    if (i/=ilabel) goto 50
    write(6,*) '====> Reading surface zone:',ilabel
    name='surface set'
    call add_ (name,i) 
    this%surfset(i)%name=name
    this%surfset(i)%numsurf=nsurf 
    this%surfset(i)%surf => null ()
    allocate(this%surfset(i)%surf(nsurf))
    do j=1,this%surfset(i)%numsurf
     if (iequil/=0) then
       this%surfset(i)%surf(j)%isequilibrate=.true. 
     else
       this%surfset(i)%surf(j)%isequilibrate=.false.  
     end if 
     this%surfset(i)%surf(j)%numsites=1
     this%surfset(i)%surf(j)%txoh => null ()
     this%surfset(i)%surf(j)%capint => null ()
     this%surfset(i)%surf(j)%capext => null ()
     this%surfset(i)%surf(j)%spsurfarea => null ()
     call check_pointer_ (this%surfset(i)%surf(j)%txoh,this%surfset(i)%surf(j)%numsites,.true.)
     call check_pointer_ (this%surfset(i)%surf(j)%capint,this%surfset(i)%surf(j)%numsites,.true.)
     call check_pointer_ (this%surfset(i)%surf(j)%capext,this%surfset(i)%surf(j)%numsites,.true.)
     call check_pointer_ (this%surfset(i)%surf(j)%spsurfarea,this%surfset(i)%surf(j)%numsites,.true.)
     name='surface'
     call add_ (name,j) 
     this%surfset(i)%surf(j)%name=name 
     read (1,*) name,ictot
     this%surfset(i)%surf(j)%capint(1)=ictot
     read (1,*) name,ictot
     this%surfset(i)%surf(j)%capext(1)=ictot
     read (1,*) name,ictot
     this%surfset(i)%surf(j)%spsurfarea(1)=ictot
     read (1,*) name,ictot,namemin
     this%surfset(i)%surf(j)%txoh(1)=ictot
     this%surfset(i)%surf(j)%namemin=namemin 
    end do
    write(6,*) '==========> Reading surface zones finished'
  end do 
end if
end if
!%-----------------------------------------------------------
!% Find label corresponding to output species
!%-----------------------------------------------------------
call find_label_ (label5,1,isbe)
if (isbe) then
!%-----------------------------------------------------------
!% Read the initial output 
!%-----------------------------------------------------------
isps=0 
write(6,*) '==========> Reading output'
do 
 isps=isps+1
 read (1,*) namesps(isps)
 if (namesps(isps)=='*') then
    isps=isps-1
    exit
 end if 
end do
call find_repeated_(namesps(1:isps),isbe,isps,name)
if (isbe) then
  msg='Error, the species is repeated:'
  call add_ (msg,name)
  iserror=.true. 
  goto 20
end if
if (isps>0) then 
 this%output%numsp=isps
 call check_pointer_ (this%output%namesp,this%output%numsp,.true.)
 this%output%namesp=namesps(1:this%output%numsp)
 call check_pointer_ (this%output%unitsp,this%output%numsp,.true.)
 this%output%unitsp='m'
 do i=1,this%output%numsp 
   call get_sp_index_ (this%pchemsys(1)%ptr,namesps(i),isps)
   if (isps==0) then
      msg='Error, not found in the chemical system the species:'
      call add_ (msg,namesps(i))
      iserror=.true. 
      goto 20
   end if
  end do
end if
write(6,*) '==========> Reading output finished'
end if


!%-----------------------------------------------------------
!% Find label corresponding to task
!%-----------------------------------------------------------
call find_label_ (label6,1,isbe)
if (isbe) then
!%-----------------------------------------------------------
!% Read the initial output 
!%-----------------------------------------------------------
isps=0 
write(6,*) '==========> Reading task'
read (1,*) ntask
!%-----------------------------------------------------------
!% Create nodal chemistry objects. Each one by task 
!%-----------------------------------------------------------
call create_nodalchem_ (this,ntask,iserror)
if (iserror) then
  msg='Error when call create_nodalchem_'
  goto 20
end if
do i=1,ntask
  read (1,*) nametask
  select case (nametask)
  
  case ('MIX')
    
    read (1,*) name
    read (1,*) nmix
    write(6,*) '===========> Mixing waters'
    do j=1,nmix
       read (1,*,err=80) prop,ilabel 
       call get_chem_info_ (this%pwater(ilabel)%ptr,iserror,name=name1)
       write(6,*) '====> Water:',j 
       write(6,*) '====> Name:',name1
       write(6,*) '====> Proportion:',prop
       if (iserror) goto 20
       if (j==1) then
         this%pnodalchemk1(i)%ptr=prop*this%pwater(ilabel)%ptr
       else
         this%pnodalchemk1(i)%ptr = this%pnodalchemk1(i)%ptr + prop*this%pwater(ilabel)%ptr  
       end if
    end do  
    
    call set_(this%pnodalchemk1(i)%ptr,iserror,name=name)
    if (iserror) then
       msg='Error when call set_'
       goto 20
    end if
    
   case ('ADD MINERAL')
    
    read (1,*) name
    write(6,*) '===========> Adding minerals'
    read (1,*,err=80) iminset,ilabel,isspecia 
    call get_chem_info_ (this%pwater(ilabel)%ptr,iserror,name=name1)
    if (iserror) goto 20
    write(6,*) '====> Water:',ilabel 
    write(6,*) '====> Mineral:'
    this%pnodalchemk1(i)%ptr=this%pnodalchemk1(ilabel)%ptr
    do j=1,this%minset(iminset)%nummin
    
      write(6,*) '====>',this%minset(iminset)%namemin(j)
      call add_(this%pnodalchemk1(i)%ptr, &
                this%minset(iminset)%namemin(j), &
                this%minset(iminset)%cmin(j), &
                this%minset(iminset)%unitcmin(j), &
                isspecia, &
                iserror)    
      if (iserror) goto 20           
       
    end do  
    
    call set_(this%pnodalchemk1(i)%ptr,iserror,name=name)
    if (iserror) then
       msg='Error when call set_'
       goto 20
    end if   
  
  case ('ADD SPECIES')
    
    read (1,*) name
    write(6,*) '===========> add species'
    read (1,*,err=80) ilabel,name2,conc 
    call get_chem_info_ (this%pwater(ilabel)%ptr,iserror,name=name1)
    if (iserror) goto 20
    write(6,*) '====> Water:',ilabel 
    write(6,*) '====> Specie:'
    this%pnodalchemk1(i)%ptr=this%pwater(ilabel)%ptr
    isspecia=.true. 
    call add_(this%pnodalchemk1(i)%ptr, &
              name2, &
              conc, &
              'mol/kg water', &
              isspecia, &
              iserror)   
    if (iserror) then
       msg='Error when call add_'
       goto 20
    end if 
    call set_(this%pnodalchemk1(i)%ptr,iserror,name=name)
    if (iserror) then
       msg='Error when call set_'
       goto 20
    end if   
    
    
  case default
      
      msg='Error, task not implemented'
      iserror=.true. 
      goto 20
  end select
end do

write(6,*) '==========> Reading task finished'
end if


!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
write(6,*) '=======> Reading CHEPROO object finished'
write(6,*) '---------------------------------------------------'
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Close unit 
!%-----------------------------------------------------------
close (unit=1)
!%-----------------------------------------------------------
!% Deallocate local pointers
!%-----------------------------------------------------------
call check_pointer_ (namesps,1,.false.)
call check_pointer_ (icon,1,.false.)
call check_pointer_ (unit,1,.false.)
call check_pointer_ (value1,1,.false.)
call check_pointer_ (value2,1,.false.)
call check_pointer_ (constraint,1,.false.)
if (iserror) goto 10 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
return
 
10 continue 
print *,'*******************'
print *,'CHEPROO:  '
print *,'Name:',this%name
print *,'Service: read_txt_ '
print *, msg
print *,'*******************'
iserror=.true.
return
30  msg='Error when open file:'
call add_ (msg,namefile)  
goto 10 
40  msg='Error reading number of water zones' 
goto 10
50  msg='Error reading index of mineral zone' 
goto 10 
60  msg='Error reading infomation of the mineral set' 
goto 10 
70 msg='Error reading label and equilibrate indices in surface set'
goto 10
80 msg='Error reading proportion and water index for task'
goto 10
return 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine init_from_backup_cheproo &
   (this, &
    iobackup, &
	iserror, &
	isopen)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Init nodal chemistry objects from output files generated
! for CHEPROO. This service is usefull for to init a CHEPROO object from 
! the results of other simulation. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)             :: this              ! Type CHEPROO variable 

integer, intent(in)                         :: iobackup          ! Backup unit 

logical, intent(out)                        :: iserror           ! iserror=true, then there was an error

logical, intent(in), optional               :: isopen            ! 
 
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
integer            :: &
 inch
logical            :: &
 haveisopen
 
iserror=.false.
msg=' ' 
haveisopen=present(isopen)
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
if (.not.haveisopen.or.(haveisopen.and..not.isopen)) then 
 open(iobackup,file='restart.dat',status='old',form='formatted')
end if 
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
do inch=1,this%numnch
   call read_backup_ (this%pnodalchemk1(inch)%ptr,iobackup,iserror)
   if (iserror) goto 20 
end do
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
do inch=1,this%numrmrmtbox
   call read_backup_ (this%prmrmtboxk1(inch)%ptr,iobackup,iserror)
   if (iserror) goto 20
end do
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
if (.not.haveisopen.or.(haveisopen.and..not.isopen)) then 
 close(iobackup)
end if 
20 continue
if (iserror) goto 10 
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
return

 
10 continue 
print *,'****************************'
print *,'CHEPROO:  '
print *,'Name:',this%name
print *,'Service: init_ '
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
subroutine write_backup_cheproo &
   (this, &
    iobackup, &
	iserror, &
	isopen)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write backup. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                :: this              ! Type CHEPROO variable 

integer, intent(in)                         :: iobackup          ! Backup unit 

logical, intent(out)                        :: iserror           ! iserror=true, then there was an error

logical, intent(in), optional               :: isopen            !  
 
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
integer            :: &
 inch
logical            :: &
 haveisopen
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
iserror=.false.
msg=' ' 
haveisopen=present(isopen)
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
if (.not.haveisopen.or.(haveisopen.and..not.isopen)) then 
 open(iobackup,file='backup.dat',status='unknown',form='formatted')
end if 
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
do inch=1,this%numnch
   call write_backup_ (this%pnodalchemk1(inch)%ptr,iobackup,iserror)
   if (iserror) goto 20
end do
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
do inch=1,this%numrmrmtbox
   call write_backup_ (this%prmrmtboxk1(inch)%ptr,iobackup,iserror)
   if (iserror) goto 20 
end do
!%--------------------------------------------------------
!% Close backup unit 
!%--------------------------------------------------------
if (.not.haveisopen.or.(haveisopen.and..not.isopen)) then 
 close(iobackup)
end if 
20 continue
if (iserror) goto 10 
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
return

 
10 continue 
print *,'****************************'
print *,'CHEPROO:  '
print *,'Name:',this%name
print *,'Service: init_ '
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
subroutine init_rmrmtbox_cheproo &
   (this, &
    iserror, &
    idwaterset, &
    idminset, &
    idsurfset, &
    dtime)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Init the nodal chemistries from :
! - water types defined in CHEPROO
! - set of mineral performed in CHEPROO
! - surface performed in cheproo 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)                             :: this         ! Type CHEPROO variable 

logical, intent(out)                                        :: iserror      ! iserror=true, then there was an error 

integer, intent(in), dimension(this%numrmrmtbox), optional  :: idwaterset   ! Indice of water [numnch]

integer, intent(in), dimension(this%numrmrmtbox), optional  :: idminset     ! Indice of mineral set [numnch]

integer, intent(in), dimension(this%numrmrmtbox), optional  :: idsurfset    ! Indice of surface set [numnch]

real*8, intent(in), optional                                :: dtime        ! Time increment [s]
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 ithnch, &
 i, &
 ndim, &
 iset
logical                                :: &
 haveidwaterset, &
 haveidminset, &
 haveidsurfset, &
 isspecia, &
 isconvergence
character(len=100)                     :: &
 name, &
 msg 
real*8                                 :: &
 volnch,  &
 volgas,  &
 omgwfree
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
msg=' '
iserror=.false.
!%------------------------------------------------------------
!% Check optional arguments
!%------------------------------------------------------------
haveidwaterset=present(idwaterset)
haveidminset=present(idminset)
haveidsurfset=present(idsurfset)
!%------------------------------------------------------------
! Start loop to nodal chemistry objects
!%------------------------------------------------------------
if (haveidwaterset.or.haveidminset.or.haveidsurfset) then
 
 isconvergence=.true. 
 
 do ithnch=1,this%numrmrmtbox
  
  pnch  => this%prmrmtboxk1(ithnch)%ptr

  
!%--------------------------------
!% Add water (optional)
!%--------------------------------
 
   if (haveidwaterset) then
    iset=idwaterset(ithnch)

    !We keep the omegafree and the node volume from the original set
    call get_chem_info_ (pnch,iserror,volnch=volnch)
	if (iserror) goto 20
	call get_chem_info_ (pnch,iserror,volgas=volgas)
	if (iserror) goto 20
    call get_chem_info_ (pnch,iserror,omgwfree=omgwfree)
    if (iserror) goto 20

    if (iset>0.and.iset<=this%numw) then
       
	   pnch=this%pwater(iset)%ptr
       name='Batch Experiment'
       call add_ (name,ithnch)
       call set_ (pnch,iserror,name=name)
	   if (iserror) goto 20
	   call set_ (pnch,iserror,volnch=volnch)
	   if (iserror) goto 20
	   call set_ (pnch,iserror,volgas=volgas)
	   if (iserror) goto 20
	   call set_ (pnch,iserror,omgwfree=omgwfree)
	   if (iserror) goto 20

    else
       msg='Error in water index:'
	   call add_ (msg,iset) 
	   goto 20
    end if
   end if
 
!%--------------------------------------
!% Add mineral set (optional)
!%--------------------------------------
   if (haveidminset) then
    
	isspecia=.false.
    
	iset=idminset(ithnch)
    
	if (iset>0.and.iset<=this%numminset) then
       
	   do i=1,this%minset(iset)%nummin

		 call add_ &
         (pnch, &
          this%minset(iset)%namemin(i), &
          this%minset(iset)%cmin(i), &
          this%minset(iset)%unitcmin(i), &
          isspecia, &
          iserror, &
          alpha=this%minset(iset)%areamin(i), &
          isconvergence=isconvergence, &
          dtime=dtime)
         
		 if (iserror) goto 10 
		 
		 if (.not.isconvergence) then
		   iserror=.true.
		   msg='Convergence problems when adding mineral:'
		   call add_ (msg,this%minset(iset)%namemin(i)) 
		   goto 20  
		 end if 
		 
		 
	   end do

    end if
   end if
 
!%--------------------------------------
!% Add surface set (optional)
!%--------------------------------------
 
   if (haveidsurfset) then
      isspecia=.false.
      iset=idsurfset(ithnch)
      if (iset>0.and.iset<=this%numsurfset) then
       do i=1,this%surfset(iset)%numsurf
         
		 ndim=this%surfset(iset)%surf(i)%numsites 
		 
		 
		 call add_ &
         (pnch, &
          this%surfset(iset)%surf(i)%name, &
          this%surfset(iset)%surf(i)%txoh(1:ndim), &
          this%surfset(iset)%surf(i)%capint(1:ndim), &
          this%surfset(iset)%surf(i)%capext(1:ndim), &
          this%surfset(iset)%surf(i)%spsurfarea(1:ndim), &
          ndim, &
          isspecia, &
          iserror, &
          isconvergence=isconvergence, &
          dtime=dtime, &
		  namemin=this%surfset(iset)%surf(i)%namemin)
		  if (iserror) goto 10 

		  if (.not.isconvergence) then
		   iserror=.true.
		   msg='Convergence problems when adding surface:'
		   call add_ (msg,this%surfset(iset)%surf(i)%name) 
		   goto 20  
		  end if 
       end do
      end if
   end if
 
!%-----------------------
   if (iserror) goto 20
!%------------------------------------------------------------
!%  
!%------------------------------------------------------------ 
 end do
20 continue
!%------------------------------------------------------------
!% Finish loop to nodal chemistry objects
!%------------------------------------------------------------
pnch => null ()
if (iserror) goto 10 
end if
!%-------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: init_rmrmtbox_ '
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
subroutine compute_f_df_rmrmt_cheproo &
  (this, &
   f, &
   df, &
   nunk, &
   npri, &
   imethod, &
   iconnbox, &
   posunk, &
   mxnbox, &
   theta, &
   alpha, &
   phi, &
   dtime, &
   isconvergence, &
   tolunk, &
   tolres, &
   mxiter, &
   mxiterminchange, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  
!
!   $Arguments:
!
 
type(t_cheproo), intent(inout)                         :: this            ! Immobile box at k+1. Type nodal chemistry variable. 

integer, intent(in)                                    :: nunk            ! dimensions of F

integer, intent(in)                                    :: npri            ! dimensions of F

real*8, dimension(nunk)                                :: f               ! F [ndim]

real*8, dimension(nunk,npri)                           :: df              ! dF [ndim,ndim]

integer, intent(in)                                    :: mxnbox          ! Maximum number of boxes per node

integer, intent(in), dimension(mxnbox+1,this%numnch)   :: iconnbox        ! Correspondence between boxes and nodes

integer, intent(in), dimension(2,this%numnch)          :: posunk          ! Position of first and last unknowns for each nodal chemistry of the mobile zone

integer, intent(in)                                    :: imethod         ! Index of reactive multi-rate mass transfer method 

real*8, intent(in)                                     :: dtime           ! Time increment (s)

real*8, intent(in)                                     :: theta           ! Temporal weight

real*8, intent(in), dimension(mxnbox,this%numnch)      :: alpha           ! Alpha for boxes 

real*8, intent(in), dimension(mxnbox,this%numnch)      :: phi             ! Porosity for boxes 

real*8, intent(in)                                     :: tolunk          ! Tolerance in the unknowns for Newton-Raphson

real*8, intent(in)                                     :: tolres          ! Tolerance in the residual for Newton-Raphson

integer, intent(in)                                    :: mxiter          ! Maximum number of iterations for Newton-Raphson

integer, intent(in)                                    :: mxiterminchange ! Maximum number of iterations for changes in mineral presence

logical, intent(out)                                   :: isconvergence   ! isconvergence=true, then there was convergence

logical, intent(out)                                   :: iserror         ! iserror=true, then there was an error
 
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
real*8, pointer              :: &
 f1(:) => null (), &
 df1(:,:) => null () 
type(t_nodalchemistry), pointer :: &
 pmobk => null (), &
 pmobk1 => null (), &
 pboxk => null (), &
 pboxk1 => null ()
integer                         :: &
 inch, &
 nbox, &
 ibox, &
 ndim1, &
 ipos1, &
 ipos2 
 
!%------------------------------------------------------------
iserror=.false. 
msg=''
!%------------------------------------------------------------
!% Initialice variables 
!%------------------------------------------------------------
f  = 0.0d0
df = 0.0d0
!%------------------------------------------------------------
!% Loop for nodal chemistry objects 
!%------------------------------------------------------------
do inch=1,this%numnch

   nbox = iconnbox(1,inch)
   pmobk1 => this%pnodalchemk1(inch)%ptr
   pmobk  => this%pnodalchemk(inch)%ptr
   ipos1=posunk(1,inch)
   ipos2=posunk(2,inch)
   
   do ibox=1,nbox
      
	  pboxk1 => this%prmrmtboxk1(iconnbox(ibox+1,inch))%ptr
      pboxk  => this%prmrmtboxk(iconnbox(ibox+1,inch))%ptr

      call compute_f_df_rmrmt_ &
          (pboxk1, &
           f1, &
           df1, &
           ndim1, &
           pboxk, &
           pmobk1, &
           pmobk, &
           theta, &
           alpha(ibox,inch), &
           phi(ibox,inch), &
           dtime, &
           isconvergence, &
           tolunk, &
           tolres, &
           mxiter, &
           mxiterminchange, &
		   imethod, &
           iserror)

	   if (iserror.or..not.isconvergence) goto 20 
       
       f(ipos1:ipos2) = f(ipos1:ipos2) + f1
       df(ipos1:ipos2,1:ndim1) = df(ipos1:ipos2,1:ndim1) + df1
       

   end do


end do 
!%-------------------------------------------------------
ndim1=0 
!%-------------------------------------------------------
20 continue 
!%-------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------
call check_pointer_ (f1,1,.false.)
call check_pointer_ (df1,1,1,.false.)
pmobk => null ()
pmobk1 => null ()
pboxk => null ()
pboxk1 => null ()
!%-------------------------------------------------------
if (iserror) goto 10 
!%-------------------------------------------------------
return
 
10 continue 
print *,'****************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: compute_f_df_rmrmt_'
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
subroutine init_nodalchem_cheproo &
   (this, &
    iserror, &
    idwaterset, &
    idminset, &
    idsurfset, &
    dtime)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Init the nodal chemistries from :
! - water types defined in CHEPROO
! - set of mineral performed in CHEPROO
! - surface performed in cheproo 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)                        :: this         ! Type CHEPROO variable 

logical, intent(out)                                   :: iserror      ! iserror=true, then there was an error 

integer, intent(in), dimension(this%numnch), optional  :: idwaterset   ! Indice of water [numnch]

integer, intent(in), dimension(this%numnch), optional  :: idminset     ! Indice of mineral set [numnch]

integer, intent(in), dimension(this%numnch), optional  :: idsurfset    ! Indice of surface set [numnch]

real*8, intent(in), optional                           :: dtime        ! Time increment [s]
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 ithnch, &
 i, &
 ndim, &
 iset
logical                                :: &
 haveidwaterset, &
 haveidminset, &
 haveidsurfset, &
 isequilibrate, &
 isconvergence
character(len=100)                     :: &
 name, &
 msg 
real*8                                 :: &
 volnch,  &
 volgas,  &
 omgwfree
type(t_nodalchemistry), pointer        :: &
 pnch => null (), &
 pnch0 => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
msg=' '
iserror=.false.
!%------------------------------------------------------------
!% Check optional arguments
!%------------------------------------------------------------
haveidwaterset=present(idwaterset)
haveidminset=present(idminset)
haveidsurfset=present(idsurfset)
!%------------------------------------------------------------
! Start loop to nodal chemistry objects
!%------------------------------------------------------------
if (haveidwaterset.or.haveidminset.or.haveidsurfset) then
 
 isconvergence=.true. 
 
 do ithnch=1,this%numnch
  
  pnch0 => this%pnodalchem0(ithnch)%ptr
  pnch  => this%pnodalchemk1(ithnch)%ptr

  
!%--------------------------------
!% Add water (optional)
!%--------------------------------
 
   if (haveidwaterset) then
    iset=idwaterset(ithnch)

    !We keep the omegafree and the node volume from the original set
    call get_chem_info_ (pnch,iserror,volnch=volnch)
	if (iserror) goto 20
	call get_chem_info_ (pnch,iserror,volgas=volgas)
	if (iserror) goto 20
    call get_chem_info_ (pnch,iserror,omgwfree=omgwfree)
    if (iserror) goto 20

    if (iset>0.and.iset<=this%numw) then
       pnch=this%pwater(iset)%ptr
       pnch0=this%pwater(iset)%ptr
       name='Batch Experiment'
       call add_ (name,ithnch)
       call set_ (pnch0,iserror,name=name)
	   if (iserror) goto 20
	   call set_ (pnch,iserror,name=name)
	   if (iserror) goto 20
	
	   call set_ (pnch0,iserror,volnch=volnch)
	   if (iserror) goto 20
	   call set_ (pnch0,iserror,volgas=volgas)
	   if (iserror) goto 20
	   call set_ (pnch,iserror,volnch=volnch)
	   if (iserror) goto 20
	   call set_ (pnch,iserror,volgas=volgas)
	   if (iserror) goto 20
	 
	   call set_ (pnch0,iserror,omgwfree=omgwfree)
	   if (iserror) goto 20
	   call set_ (pnch,iserror,omgwfree=omgwfree)
	   if (iserror) goto 20

    else
       msg='Error in water index:'
	   call add_ (msg,iset) 
	   goto 20
    end if
   end if
 
!%--------------------------------------
!% Add mineral set (optional)
!%--------------------------------------
   if (haveidminset) then
    
	isequilibrate=.false.
    
	iset=idminset(ithnch)
    
	if (iset>0.and.iset<=this%numminset) then
       
	   do i=1,this%minset(iset)%nummin

		 call add_ &
         (pnch0, &
          this%minset(iset)%namemin(i), &
          this%minset(iset)%cmin(i), &
          this%minset(iset)%unitcmin(i), &
          isequilibrate, &
          iserror, &
          alpha=this%minset(iset)%areamin(i), &
          isconvergence=isconvergence, &
          dtime=dtime)
         if (iserror) goto 10 

		 call add_ &
         (pnch, &
          this%minset(iset)%namemin(i), &
          this%minset(iset)%cmin(i), &
          this%minset(iset)%unitcmin(i), &
          isequilibrate, &
          iserror, &
          alpha=this%minset(iset)%areamin(i), &
          isconvergence=isconvergence, &
          dtime=dtime)
         
		 if (iserror) goto 10 
		 
		 if (.not.isconvergence) then
		   iserror=.true.
		   msg='Convergence problems when adding mineral:'
		   call add_ (msg,this%minset(iset)%namemin(i)) 
		   goto 20  
		 end if  		 
	   
	   end do

    end if
   end if
 
!%--------------------------------------
!% Add surface set (optional)
!%--------------------------------------
 
   if (haveidsurfset) then
      iset=idsurfset(ithnch)
      if (iset>0.and.iset<=this%numsurfset) then
       do i=1,this%surfset(iset)%numsurf
         ndim=this%surfset(iset)%surf(i)%numsites 
		 
		 call add_ &
         (pnch0, &
          this%surfset(iset)%surf(i)%name, &
          this%surfset(iset)%surf(i)%txoh(1:ndim), &
          this%surfset(iset)%surf(i)%capint(1:ndim), &
          this%surfset(iset)%surf(i)%capext(1:ndim), &
          this%surfset(iset)%surf(i)%spsurfarea(1:ndim), &
          ndim, &
          this%surfset(iset)%surf(i)%isequilibrate, &
          iserror, &
          isconvergence=isconvergence, &
          dtime=dtime, &
		  namemin=this%surfset(iset)%surf(i)%namemin)
		 if (iserror) goto 10 
		 
		 call add_ &
         (pnch, &
          this%surfset(iset)%surf(i)%name, &
          this%surfset(iset)%surf(i)%txoh(1:ndim), &
          this%surfset(iset)%surf(i)%capint(1:ndim), &
          this%surfset(iset)%surf(i)%capext(1:ndim), &
          this%surfset(iset)%surf(i)%spsurfarea(1:ndim), &
          ndim, &
          this%surfset(iset)%surf(i)%isequilibrate, &
          iserror, &
          isconvergence=isconvergence, &
          dtime=dtime, &
		  namemin=this%surfset(iset)%surf(i)%namemin)
		  if (iserror) goto 10 

		  if (.not.isconvergence) then
		   iserror=.true.
		   msg='Convergence problems when adding surface:'
		   call add_ (msg,this%surfset(iset)%surf(i)%name) 
		   goto 20  
		  end if 
       end do
      end if
   end if
 
!%-----------------------
   if (iserror) goto 20
!%------------------------------------------------------------
!%  
!%------------------------------------------------------------ 
 end do
20 continue
!%------------------------------------------------------------
!% Finish loop to nodal chemistry objects
!%------------------------------------------------------------
pnch0 => null ()
pnch => null ()
if (iserror) goto 10 
end if
!%-------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: init_chem_'
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
subroutine init_nodalchem_from_xml_cheproo &
   (this, &
    namefile, &
    iserror)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Init nodal chemistry objects from xml file 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout):: this       ! Type CHEPROO variable 

character(len=*), intent(in)   :: namefile   ! Name of the file 

logical, intent(out)           :: iserror    ! iserror=true, then there was an error
 
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
 name
type(dictionary_t) :: &
 attributes
integer:: &
 iostat
type(xml_t):: &
 fxml
character(len=100) :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
! Open xml file 
!%----------------------------------------------------------
call open_xmlfile(namefile, fxml, iostat)
if (iostat/=0) then
 msg='Error when open file'
 call add_ (msg,namefile)
 goto 10
end if
!%-----------------------------------------------------------
call xml_parse &
 (fxml, &
  begin_element_handler = begin_element_handler)
!%--------------------------------------------------------
! End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%--------------------------------------------------------
!%--------------------------------------------------------
return
 
10 continue 
print *,'****************************'
print *,'CHEPROO:  '
print *,'Name:',this%name
print *,'Service: init_nodalchem_ '
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
subroutine init_nodalchem_from_output_cheproo &
   (this, &
    namesptotfile, &
	namespconcfile, &
	timelabel, &
	isconvergence, &
	iserror, &
	namespminfile, &
	namespparamfile, &
	dtime)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Init nodal chemistry objects from output files generated
! for CHEPROO. This service is usefull for to init a CHEPROO object from 
! the results of other simulation. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)             :: this              ! Type CHEPROO variable 

character(len=*), intent(in)                :: namesptotfile     ! spatial_tot_cheproo.dat 

character(len=*), intent(in)                :: namespconcfile    ! spatial_conc_cheproo.dat 

character(len=*), intent(in)                :: timelabel         ! Time label 

logical, intent(out)                        :: isconvergence     ! isconvergence=true, then there was convergence 

logical, intent(out)                        :: iserror           ! iserror=true, then there was an error

character(len=*), intent(in), optional      :: namespminfile     ! spatial_tot_cheproo.dat 

character(len=*), intent(in), optional      :: namespparamfile   ! spatial_tot_cheproo.dat 

real*8, intent(in), optional                :: dtime             ! time increment [s]
 
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
 msg, &
 label 
logical            :: &
 havedtime, &
 haveiouspvolfracmin, &
 havenamespparamfile, &
 havenamespminfile, &
 isbelabel, &
 isbe 
integer            :: &
 ncomp, &
 i, &
 j, &
 k, &
 inch, &
 ipos1, &
 ipos2, &
 ipos3, &
 nmin, &
 nsp 
character(len=100), pointer :: &
 namelabel(:) => null (), &
 namelabel1(:) => null (), &
 constraint(:) => null (), &
 unit(:) => null (), &
 namemin(:) => null ()
integer, pointer            :: &
 icon(:) => null (), &
 id(:) => null ()  
real*8, pointer             :: &
 value(:) => null (), &
 value1(:) => null ()  
real*8                      :: &
 cputime, &
 temp  
type(t_nodalchemistry), pointer :: &
 pnch => null () 
integer, parameter              :: &
 nparam=10 
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
!% Check optional arguments 
!%----------------------------------------------------------
havedtime=present(dtime)
havenamespparamfile=present(namespparamfile)
havenamespminfile=present(namespminfile)
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
isconvergence=.false. 
!%--------------------------------------------------------
!% Open files
!%--------------------------------------------------------
open(unit=1,file=namesptotfile,status='old',err=20)
open(unit=3,file=namespconcfile,status='old',err=20)
if(havenamespparamfile) then 
   open(unit=2,file=namespparamfile,status='old',err=20)
   !%--------------------------------------------------------
   !% Find the time label 
   !%--------------------------------------------------------
   isbelabel=.false. 
    do 
      read(1,40) label
      if (label==' ') then
        msg='Error, time label not found' 
        goto 10 
      end if 
      call find(label,'Time',isbe) 
      if (isbe) then
        call find(label,timelabel,isbe)
        if (isbe) then   
          isbelabel=.true. 
          exit
        end if 
      end if 
    end do 	 
    call check_pointer_ (namelabel,nparam,.true.)
    call check_pointer_ (value,nparam,.true.)
    if (isbelabel) then
      read(1,*)  
      read(1,*)  (namelabel(j), j=1,ncomp)
	  call findInArray_(namelabel,'volume',ipos1)
	  call findInArray_(namelabel,'temperature',ipos2)
	  call findInArray_(namelabel,'mass of water',ipos3)
      do i=1,this%numnch
         pnch => this%pnodalchemk1(i)%ptr 
         read(1,*)  inch,(value(j), j=1,nparam) 
	     if (ipos1>0) then   
		   call set_ (pnch,iserror,volnch=value(ipos1))
		   if (iserror) goto 30
		 end if
		 if (ipos2>0) then   
		   call set_ (pnch,iserror,temp=value(ipos2))
		   if (iserror) goto 30
		 end if
		 if (ipos3>0) then   
		   call set_ (pnch,iserror,omgwfree=value(ipos3))
		   if (iserror) goto 30
		 end if 
      end do
    end if 

end if 
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!% Find the time label in spatial_tot_cheproo.dat
!%--------------------------------------------------------
isbelabel=.false. 
do 
 read(1,40) label
 if (label==' ') then
   msg='Error, time label not found' 
   goto 10 
 end if 
 call find(label,'Time',isbe) 
 if (isbe) then
   call find(label,timelabel,isbe)
   if (isbe) then   
      isbelabel=.true. 
      exit
   end if 
 end if 
end do 
!%--------------------------------------------------------
!% Find the time label in spatial_conc_cheproo.dat
!%--------------------------------------------------------
isbelabel=.false. 
do 
 read(3,40) label
 if (label==' ') then
   msg='Error, time label not found' 
   goto 10 
 end if 
 call find(label,'Time',isbe) 
 if (isbe) then
   call find(label,timelabel,isbe)
   if (isbe) then   
      isbelabel=.true. 
      exit
   end if 
 end if 
end do 
!%--------------------------------------------------------
!% 
!%--------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numbase=ncomp,numsp=nsp)
if (iserror) goto 30 
call check_pointer_ (namelabel,ncomp,.true.)
call check_pointer_ (id,ncomp,.true.)
call check_pointer_ (namelabel1,nsp,.true.)
call check_pointer_ (value,ncomp,.true.)
call check_pointer_ (value1,nsp,.true.)
call check_pointer_ (unit,ncomp,.true.)
call check_pointer_ (icon,ncomp,.true.)
call check_pointer_ (constraint,ncomp,.true.)
icon=1
unit='m'
if (isbelabel) then
  read(1,*)  
  read(3,*)  
  read(1,*)  (namelabel(j), j=1,ncomp)
  read(3,*)  (namelabel1(j), j=1,nsp)
  do i=1,ncomp
    call findInArray_(namelabel1,namelabel(i),ipos1)
    if (ipos1>0) then
       id(i)=ipos1
	end if
  end do 
  do i=1,this%numnch
    pnch => this%pnodalchemk1(i)%ptr
	call set_ (pnch,iserror,isderstored=.true.)
    if (iserror) goto 30
	call get_chem_info_ (pnch,iserror,temp=temp)
    if (iserror) goto 30 
    read(1,*)  inch,(value(j), j=1,ncomp)
	read(3,*)  inch,(value1(j), j=1,nsp)
	call set_from_solution_type_ (pnch,namelabel,value1(id),value,unit,icon,constraint,ncomp,temp,isconvergence,cputime,iserror) 
    if (iserror) then 
	   msg='Error when call set_from_solution_type_'
	   goto 30 
	end if
	if (.not.isconvergence) return 
     
  end do
end if 
!%--------------------------------------------------------
!% Read the volumetric fraction of minerals 
!%--------------------------------------------------------
if(havenamespminfile) then 
   open(unit=3,file=namespminfile,status='old',err=20)
   !%--------------------------------------------------------
   !% Find the time label 
   !%--------------------------------------------------------
   isbelabel=.false. 
   do 
      read(1,40) label
      if (label==' ') then
        msg='Error, time label not found' 
        goto 10 
      end if 
      call find(label,'Time',isbe) 
      if (isbe) then
        call find(label,timelabel,isbe)
        if (isbe) then   
          isbelabel=.true. 
          exit
        end if 
      end if 
    end do 	
	call get_chem_info_ (this%pchemsys(1)%ptr,iserror,nminsp=nmin,nameminsp=namemin)
    if (iserror) goto 30  
    call check_pointer_ (namelabel,nmin,.true.)
    call check_pointer_ (value,nmin,.true.)
    if (isbelabel.and.nmin>0) then
      read(1,*)  
      read(1,*)  (namelabel(j), j=1,nmin)
      do j=1,nmin 
	    call findInArray_(namelabel,namemin(j),ipos1)  
        do i=1,this%numnch
           pnch => this%pnodalchemk1(i)%ptr 
           read(1,*)  inch,(value(k), k=1,nmin) 
	       if (ipos1>0) then   
		      call add_ (pnch,namelabel(ipos1),value(ipos1),'m3min/m3rock',.false.,iserror,1.0d-5,isconvergence,dtime=dtime)
		      if (iserror) goto 30
		   end if
        end do
      end do 
	end if 

end if
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
30 continue
call check_pointer_ (namelabel,1,.false.)
call check_pointer_ (namelabel1,1,.false.)
call check_pointer_ (value,1,.false.)
call check_pointer_ (value1,1,.false.)
call check_pointer_ (id,1,.false.)
call check_pointer_ (unit,1,.false.)
call check_pointer_ (icon,1,.false.)
call check_pointer_ (constraint,1,.false.)
call check_pointer_ (namemin,1,.false.)
pnch => null ()
if (iserror) goto 10 
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
return

20  msg='Error when open file'
goto 10

40 format(a50)
 
10 continue 
print *,'****************************'
print *,'CHEPROO:  '
print *,'Name:',this%name
print *,'Service: init_nodalchem_ '
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
subroutine set_cheproo &
   (this, &
    iserror, &
    name, &
    temp, &
    omgwfree, &
    volgas, &
	pgas, &
	pliq, &
    volnch, &
	liqdens, &
	iswriteinfo, &
	tempbox, &
	omgwfreebox, &
	volbox)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set general attributes in the CHEPROO object. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)                            :: this              ! Type CHEPROO variable. 
 
 character(len=*), intent(in), optional                     :: name             ! Name of the CHEPROO object. 

logical, intent(out)                                       :: iserror           ! iserror=true, there was an error. 

real*8, intent(in), dimension(this%numnch), optional       :: temp              ! Temperature [number of nodal chemistry] (in celcius). 

real*8, intent(in), dimension(this%numnch), optional       :: omgwfree          ! Mass of free water  [number of nodal chemistry] [kgw]. 

real*8, intent(in), dimension(this%numnch), optional       :: volgas            ! Gas volume [number of nodal chemistry] [m3/m3rock]

real*8, intent(in), dimension(this%numnch), optional       :: pgas              ! Gas preassure [number of nodal chemistry] [Mpa]

real*8, intent(in), dimension(this%numnch), optional       :: pliq              ! liquid preassure [number of nodal chemistry] [Mpa]

real*8, intent(in), dimension(this%numnch), optional       :: liqdens           ! liquid density [number of nodal chemistry] [kg/m3]

real*8, intent(in), dimension(this%numnch), optional       :: volnch            ! Volume [number of nodal chemistry] [m3]

real*8, intent(in), dimension(this%numrmrmtbox), optional  :: omgwfreebox       ! Mass of free water  [number of nodal chemistry] [kgw]. 

real*8, intent(in), dimension(this%numrmrmtbox), optional  :: tempbox           ! Temperature [number of rmrmt boxes] (in celcius). 

real*8, intent(in), dimension(this%numrmrmtbox), optional  :: volbox            ! Volume [number of nodal chemistry] [m3]

logical, intent(in), optional                              :: iswriteinfo       ! iswriteinfo=true, write general information.  
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 ndim1, &
 i
logical                                :: &
 havename, &
 havetemp, &
 havetempbox, &
 haveomgwfree, &
 haveomgwfreebox, &
 havevolgas, &
 havepgas, &
 havepliq, &
 haveliqdens, &
 havevolnch, &
 havevolbox, &
 haveiswriteinfo
character(len=100)                     :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check optional arguments
!%------------------------------------------------------------
havename=present(name)
havetemp=present(temp)
havetempbox=present(tempbox)
haveomgwfree=present(omgwfree)
havevolgas=present(volgas)
havevolnch=present(volnch)
havepgas=present(pgas)
havepliq=present(pliq)
haveliqdens=present(liqdens)
haveiswriteinfo=present(iswriteinfo)
haveomgwfreebox=present(omgwfreebox)
havevolbox=present(volbox)
!%------------------------------------------------------------
!% Set the name of CHEPROO object 
!%------------------------------------------------------------
if (havename) then
 this%name=name
end if
!%------------------------------------------------------------
!% Set the temperature 
!%------------------------------------------------------------
if (havetemp) then
 do i=1,this%numnch
  call set_ (this%pnodalchemk1(i)%ptr,iserror,temp=temp(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set the temperature in nodal chemistry boxes 
!%------------------------------------------------------------
if (havetempbox) then
 do i=1,this%numrmrmtbox
  call set_ (this%prmrmtboxk1(i)%ptr,iserror,temp=tempbox(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set the liquid density 
!%------------------------------------------------------------
if (haveliqdens) then
 do i=1,this%numnch
  call set_ (this%pnodalchemk1(i)%ptr,iserror,liqdens=liqdens(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set the mass of free water 
!%------------------------------------------------------------
if (haveomgwfree) then
 do i=1,this%numnch
  call set_ (this%pnodalchemk1(i)%ptr,iserror,omgwfree=omgwfree(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set the mass of free water in nodal chemistry boxes
!%------------------------------------------------------------
if (haveomgwfreebox) then
 do i=1,this%numrmrmtbox
  call set_ (this%prmrmtboxk1(i)%ptr,iserror,omgwfree=omgwfreebox(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set the nodal chemistry volume 
!%------------------------------------------------------------
if (havevolnch) then
 do i=1,this%numnch
  call set_ (this%pnodalchemk1(i)%ptr,iserror,volnch=volnch(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set the nodal chemistry volume in nodal chemistry boxes
!%------------------------------------------------------------
if (havevolbox) then
 do i=1,this%numnch
  call set_ (this%pnodalchemk1(i)%ptr,iserror,volnch=volbox(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set the gas volume 
!%------------------------------------------------------------
if (havevolgas) then
 do i=1,this%numnch
  call set_ (this%pnodalchemk1(i)%ptr,iserror,volgas=volgas(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set the gas preassure
!%------------------------------------------------------------
if (havepgas) then
 do i=1,this%numnch
  call set_ (this%pnodalchemk1(i)%ptr,iserror,pgas=pgas(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set the liquid preassure
!%------------------------------------------------------------
if (havepliq) then
 do i=1,this%numnch
  call set_ (this%pnodalchemk1(i)%ptr,iserror,pliq=pliq(i))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Set iswriteinfo
!%------------------------------------------------------------
if (haveiswriteinfo) then
 this%iswriteinfo=iswriteinfo
end if
!%------------------------------------------------------------
return
10 continue
print *,'*******************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: set_'
print *,msg
print *,'*******************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_chem_info_cheproo &
   (this, &
    iserror, &
    xsalt, &
    xwater, &
    actwater, &
    omgwcryst, &
	faccap, &
    ispsw, &
    icw, &
    namebase, &
    numw, &
    nummobph,&
	nummobsp, &
	numbase, &
	numnch)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return general chemical information 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                         :: this           ! Type CHEPROO variable 

logical, intent(out)                                 :: iserror        ! iserror=true, then there was an error 

real*8, pointer, optional, dimension(:)              :: xsalt          ! Mass fraction of salt in liquid. 

real*8, pointer, optional, dimension(:)              :: xwater         ! Fraction of water in liquid. 
 
real*8, pointer, optional, dimension(:)              :: actwater       ! Water activity [numnch]

real*8, pointer, optional, dimension(:)              :: faccap         ! Capillary correction for water activity 

real*8, pointer, optional, dimension(:)              :: omgwcryst      ! Mass of water fixed in minerals

character(len=*), pointer, optional, dimension(:)    :: namebase       ! Name list of components

integer, intent(out), optional                       :: numw           ! Number of waters 

integer, intent(out), optional                       :: ispsw          ! Index of water species

integer, intent(out), optional                       :: icw            ! Global index of water component 

integer, intent(out), optional                       :: nummobph       ! Number of mobile phases. 

integer, intent(out), optional                       :: numbase        ! Number of chemical components 

integer, intent(out), optional                       :: nummobsp       ! Number of mobile species

integer, intent(out), optional                       :: numnch         ! Number of nodal chemistry objects
 
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
 inch
character(len=100)                      :: &
 msg
logical                                 :: &
 havexsalt, &
 havexwater, &
 haveactwater, &
 haveomgwcryst, &
 haveispsw, &
 haveicw, &
 havenummobph, &
 havenummobsp, &
 havefaccap, &
 havenumw, &
 havenamebase, &
 havenumbase, &
 havenumnch
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check the number of chemical system objects 
!%------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not associated chemical system'
 goto 10
end if
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
havexsalt=present(xsalt)
havexwater=present(xwater)
haveactwater=present(actwater)
haveomgwcryst=present(omgwcryst)
haveispsw=present(ispsw)
haveicw=present(icw)
havenummobph=present(nummobph)
havenummobsp=present(nummobsp)
havefaccap=present(faccap)
havenumw=present(numw)
havenamebase=present(namebase)
havenumbase=present(numbase)
havenumnch=present(numnch)
!%------------------------------------------------------------
!% Return the salt fraction in liquid vector
!%------------------------------------------------------------
if (havexsalt) then
 call check_pointer_ (xsalt,this%numnch,.true.)
 do inch=1,this%numnch
   call get_chem_info_ (this%pnodalchemk1(inch)%ptr,iserror,xsalt=xsalt(inch))
   if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Return the capillary correction for water activity 
!%------------------------------------------------------------
if (havefaccap) then
 call check_pointer_ (faccap,this%numnch,.true.)
 do inch=1,this%numnch
   call get_chem_info_ (this%pnodalchemk1(inch)%ptr,iserror,faccap=faccap(inch))
   if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Return the water fraction in liquid vector 
!%------------------------------------------------------------
if (havexwater) then
 call check_pointer_ (xwater,this%numnch,.true.)
 do inch=1,this%numnch
  call get_chem_info_ (this%pnodalchemk1(inch)%ptr,iserror,xwater=xwater(inch))
  if (iserror) goto 10 
 end do
end if
!%------------------------------------------------------------
!% Return the water activity vector 
!%------------------------------------------------------------
if (haveactwater) then
 call check_pointer_ (actwater,this%numnch,.true.)
 do inch=1,this%numnch
  call get_chem_info_ (this%pnodalchemk1(inch)%ptr,iserror,actwater=actwater(inch))
  if (iserror) goto 10
 end do
end if
!%------------------------------------------------------------
!% Return the vector of mass of water in minerals 
!%------------------------------------------------------------
if (haveomgwcryst) then
 call check_pointer_ (omgwcryst,this%numnch,.true.)
 do inch=1,this%numnch
  call get_chem_info_ (this%pnodalchemk1(inch)%ptr,iserror,omgwcryst=omgwcryst(inch))
  if (iserror) goto 10
 end do
end if
!%-----------------------------------------------------------
!% Return the global index of water species 
!%-----------------------------------------------------------
if (haveispsw) then
 call get_chem_info_(this%pchemsys(1)%ptr,iserror,wspindex=ispsw)
 if (iserror) goto 10
end if
!%------------------------------------------------------------
!% Return the global index of water component 
!%------------------------------------------------------------
if (haveicw) then
 call get_chem_info_ (this%pchemsys(1)%ptr,iserror,wcindex=icw)
 if (iserror) goto 10
end if
!%------------------------------------------------------------
!% Return the number of mobile phases 
!%------------------------------------------------------------
if (havenummobph) then
 call get_chem_info_ (this%pchemsys(1)%ptr,iserror,nummobph=nummobph)
 if (iserror) goto 10
end if
!%------------------------------------------------------------
!% Return the number of mobile speceies 
!%------------------------------------------------------------
if (havenummobsp) then
 call get_chem_info_ (this%pchemsys(1)%ptr,iserror,nummobsp=nummobsp)
 if (iserror) goto 10
end if
!%------------------------------------------------------------
!% Return the number of waters 
!%------------------------------------------------------------
if (havenumw) then
   numw=this%numw
end if
!%------------------------------------------------------------
!% Return the name list of chemical components
!%------------------------------------------------------------
if (havenamebase) then
   call get_chem_info_ (this%pchemsys(1)%ptr,iserror,namebase=namebase)
   if (iserror) goto 10 
end if
!%------------------------------------------------------------
!% Return the number of components 
!%------------------------------------------------------------
if (havenumbase) then
   call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numbase=numbase)
   if (iserror) goto 10 
end if
!%------------------------------------------------------------
!% Return the number of nodal chemistry objects
!%------------------------------------------------------------
if (havenumnch) then
   numnch=this%numnch
end if
!%------------------------------------------------------------
return
 
10 continue  
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_chem_info_'
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
subroutine get_chem_info_ith_cheproo &
   (this, &
    ith, &
    iserror, &
    xsalt, &
    xwater, &
    dxwaterdc, &
    actwater, &
    omgwcryst, &
    ispsw, &
    icw, &
	nummobph,&
	nummobsp)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return general chemical information about a specific nodal chemistry
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)               :: this           ! Type CHEPROO variable 

integer,intent(in)                         :: ith            ! Index of the nodal chemistry

logical, intent(out)                       :: iserror        ! iserror=true, then there was an error 

real*8,optional                            :: xsalt          ! Mass fraction of salt in liquid. 

real*8, optional                           :: xwater         ! Fraction of water in liquid. 
 
real*8, optional,pointer,dimension(:)      :: dxwaterdc      ! derivative of water mass fraction with respect to primary species

real*8, optional                           :: actwater       ! Water activity

real*8, optional                           :: omgwcryst      ! Mass of water fixed in minerals

integer, intent(out), optional             :: ispsw          ! Global index of water species 

integer, intent(out), optional             :: icw            ! Global index of water component 

integer, intent(out), optional             :: nummobph       ! Number of mobile phases. 

integer, intent(out), optional             :: nummobsp       ! Number of mobile species

 
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
 inch
character(len=100)                      :: &
 msg
logical                                 :: &
 havexsalt, &
 havexwater, &
 havedxwaterdc, &
 haveactwater, &
 haveomgwcryst, &
 haveispsw, &
 haveicw, &
 havenummobph, &
 havenummobsp
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check the number of chemical system objects 
!%------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not associated chemical system'
 goto 10
end if
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
havexsalt=present(xsalt)
havexwater=present(xwater)
havedxwaterdc=present(dxwaterdc)
haveactwater=present(actwater)
haveomgwcryst=present(omgwcryst)
haveispsw=present(ispsw)
haveicw=present(icw)
havenummobph=present(nummobph)
havenummobsp=present(nummobsp)
!%------------------------------------------------------------
!% Return the salt fraction in liquid vector
!%------------------------------------------------------------
if (havexsalt) then

   call get_chem_info_ (this%pnodalchemk1(ith)%ptr,iserror,xsalt=xsalt)
   if (iserror) goto 10

end if
!%------------------------------------------------------------
!% Return the water fraction in liquid vector 
!%------------------------------------------------------------
if (havexwater) then

  call get_chem_info_ (this%pnodalchemk1(ith)%ptr,iserror,xwater=xwater)
  if (iserror) goto 10 

end if
!%------------------------------------------------------------
!% Return the derivative from water fraction in liquid vector 
!%------------------------------------------------------------
if (havedxwaterdc) then
 
  call get_chem_info_ (this%pnodalchemk1(ith)%ptr,iserror,dxwaterdc=dxwaterdc)
  if (iserror) goto 10 
 
end if
!%------------------------------------------------------------
!% Return the water activity vector 
!%------------------------------------------------------------
if (haveactwater) then
  call get_chem_info_ (this%pnodalchemk1(ith)%ptr,iserror,actwater=actwater)
  if (iserror) goto 10

end if
!%------------------------------------------------------------
!% Return the vector of mass of water in minerals 
!%------------------------------------------------------------
if (haveomgwcryst) then
  call get_chem_info_ (this%pnodalchemk1(ith)%ptr,iserror,omgwcryst=omgwcryst)
  if (iserror) goto 10
 end if
!%-----------------------------------------------------------
!% Return the global index of water species 
!%-----------------------------------------------------------
if (haveispsw) then
 call get_chem_info_(this%pchemsys(1)%ptr,iserror,wspindex=ispsw)
 if (iserror) goto 10
end if
!%------------------------------------------------------------
!% Return the global index of water component 
!%------------------------------------------------------------
if (haveicw) then
 call get_chem_info_ (this%pchemsys(1)%ptr,iserror,wcindex=icw)
 if (iserror) goto 10
end if
!%------------------------------------------------------------
!% Return the number of mobile phases 
!%------------------------------------------------------------
if (havenummobph) then
 call get_chem_info_ (this%pchemsys(1)%ptr,iserror,nummobph=nummobph)
 if (iserror) goto 10
end if
!%------------------------------------------------------------
!% Return the number of mobile speceies 
!%------------------------------------------------------------
if (havenummobsp) then
 call get_chem_info_ (this%pchemsys(1)%ptr,iserror,nummobsp=nummobsp)
 if (iserror) goto 10
end if
!%------------------------------------------------------------



return
 
10 continue  
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_chem_info_ith'
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
subroutine get_iumob_ith_cheproo &
  (this, &
   iumob, &
   naqcol, &
   ngascol, &
   nnonaqcol, &
   ithcomp, &
   ithnch, &
   isbackwards, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Return the ithcomp aqueous component in the mobile phases defined in the chemical system for ith nodal chemistry. 
!
!   $Arguments:
!
 
 
 
type(t_cheproo), intent(in) :: this        ! Type CHEPROO variable

real*8, pointer             :: iumob(:)

integer, intent(out)        :: naqcol

integer, intent(out)        :: ngascol

integer, intent(out)        :: nnonaqcol

integer, intent(in)         :: ithcomp     ! Index of ith components 

integer, intent(in)         :: ithnch      ! Index of nodal chemistry object 

logical, intent(in)         :: isbackwards ! If .true. then get in k list 

logical, intent(out)        :: iserror     ! iserror=true, then there was an error 
 
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
 numsp
character(len=100)     :: &
 msg
type(t_nodalchemistry), pointer :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
if (ithnch<=0.or.ithnch>this%numnch) then
  msg='Error in nodal chemistry index, index='
  call add_ (msg,ithnch) 
  goto 10
end if
!%--------------------------------------------------------------
if (isbackwards) then 
 pnch => this%pnodalchemk(ithnch)%ptr 
else
 pnch => this%pnodalchemk1(ithnch)%ptr 
end if
!%--------------------------------------------------------------
call get_iumob_ &
  (pnch, &
   iumob, &
   naqcol, &
   ngascol, &
   nnonaqcol, &
   ithcomp, &
   iserror)
 
if (iserror) then
 msg='Error when call public service (get_iumob_)'
 goto 10
end if
!%-------------------------------------------------------------- 
nullify(pnch)
!%--------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_iumob_ith_'
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
subroutine get_iumob_water_cheproo &
  (this, &
   iumob, &
   naqcol, &
   ngascol, &
   nnonaqcol, &
   ithcomp, &
   ithw, &
   iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Reaturn the ithcomp aqueous component in the mobile phases defined in the chemical system for ithw water type. 
!
!   $Arguments:
!
 
type(t_cheproo), intent(in)  :: this

real*8, pointer              :: iumob(:)

integer, intent(out)         :: naqcol

integer, intent(out)         :: ngascol

integer, intent(out)         :: nnonaqcol

integer, intent(in)          :: ithcomp

integer, intent(in)          :: ithw

logical, intent(out)         :: iserror      ! iserror=true, then there was an error 
 
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
 numsp
character(len=100)     :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%--------------------------------------------------------------
!% Check index 
!%--------------------------------------------------------------
if (ithw<=0.or.ithw>this%numw) then
  msg='Error in water index, index='
  call add_ (msg,ithw)
  goto 10
end if
!%--------------------------------------------------------------
!%  
!%--------------------------------------------------------------
call get_iumob_ &
  (this%pwater(ithw)%ptr, &
   iumob, &
   naqcol, &
   ngascol, &
   nnonaqcol, &
   ithcomp, &
   iserror)
if (iserror) then
 msg='Error when calling get_iumob_'
 goto 10
end if
!%--------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_iumob_water_'
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
subroutine set_from_uaq_ith_cheproo &
   (this, &
    unew, &
    dtime, &
    isconvergence, &
    ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set from u the ithnch nodal chemistry object. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)           :: this           ! Type CHEPROO variable 

real*8, intent(in), dimension(:)          :: unew 

real*8, intent(in)                        :: dtime          ! Time increment [s]

integer, intent(in)                       :: ithnch         ! Index of nodal chemistry 

logical, intent(out)                      :: iserror        ! iserror=true, then there was an error 

logical, intent(out)                      :: isconvergence  ! isconvergence=true, then there was convergence
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
character(len=100)                        :: &
 msg
integer                                   :: &
 npri
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=' '
!%--------------------------------------------------------------
!% Check nodal chemistry index
!%--------------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
  msg='Error in nodal chmistry index, index='
  call add_ (msg,ithnch)
  goto 10
end if
!%--------------------------------------------------------------
npri=size(unew)
call set_from_uaq_ &
  (this%pnodalchemk1(ithnch)%ptr, &
   unew, &
   npri, &
   dtime, &
   isconvergence, &
   iserror)
if (iserror) then
 msg='Error when calling set_from_uaq_'
 call add_ (msg, ithnch)
 goto 10
end if
!%--------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: set_from_uaq_ith_'
print *,msg
print *,'*************************'
iserror=.true.
return
 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_nodalchem_cheproo &
   (this, &
    numnch, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create nodal cheistry objects. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)       :: this     ! Type CHEPROO variable 

integer, intent(in)                   :: numnch   ! Number of nodal chemistry objects 

logical, intent(out)                  :: iserror  ! iserror=true, then there was an error 
 
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
 inch
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%------------------------------------------------------------
!% If the nodal chemistry objects were previously created then 
!% destroy and deallocate the memeory space 
!%------------------------------------------------------------
if (this%numnch>0) then
 do inch=1,this%numnch
  call destroy_ (this%pnodalchem0(inch)%ptr)
  call destroy_ (this%pnodalchemk(inch)%ptr)
  call destroy_ (this%pnodalchemk1(inch)%ptr)
  deallocate (this%pnodalchem0(inch)%ptr)
  deallocate (this%pnodalchemk(inch)%ptr)
  deallocate (this%pnodalchemk1(inch)%ptr)
  this%pnodalchem0(inch)%ptr => null ()
  this%pnodalchemk(inch)%ptr => null ()
  this%pnodalchemk1(inch)%ptr => null ()
 end do
 deallocate (this%pnodalchemk)
 deallocate (this%pnodalchemk1)
end if
!%------------------------------------------------------------
this%numnch=numnch
!%------------------------------------------------------------
if (this%numnch>0) then
 allocate (this%pnodalchem0(this%numnch))
 allocate (this%pnodalchemk(this%numnch))
 allocate (this%pnodalchemk1(this%numnch))
 do inch=1,this%numnch
   allocate (this%pnodalchem0(inch)%ptr)
   allocate (this%pnodalchemk(inch)%ptr)
   allocate (this%pnodalchemk1(inch)%ptr)
   call create_ (this%pnodalchem0(inch)%ptr)
   call create_ (this%pnodalchemk(inch)%ptr)
   call create_ (this%pnodalchemk1(inch)%ptr)
   call set_ (this%pnodalchem0(inch)%ptr,iserror,chemsys=this%pchemsys(1)%ptr)
   if (iserror) goto 10 
   call set_ (this%pnodalchemk(inch)%ptr,iserror,chemsys=this%pchemsys(1)%ptr)
   if (iserror) goto 10 
   call set_ (this%pnodalchemk1(inch)%ptr,iserror,chemsys=this%pchemsys(1)%ptr)
   if (iserror) goto 10 
 end do
end if
!%------------------------------------------------------------
return
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: create_nodalchem_'
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
subroutine create_rmrmtbox_cheproo &
   (this, &
    numrmrmtbox, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create nodal chemistry objects for multirate mass transfer. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)       :: this         ! Type CHEPROO variable 

integer, intent(in)                   :: numrmrmtbox  ! Number of nodal chemistry objects (box)

logical, intent(out)                  :: iserror      ! iserror=true, then there was an error 
 
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
 inch
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%------------------------------------------------------------
!% If the nodal chemistry objects were previously created then 
!% destroy and deallocate the memeory space 
!%------------------------------------------------------------
if (this%numrmrmtbox>0) then
 do inch=1,this%numrmrmtbox
  call destroy_ (this%prmrmtboxk1(inch)%ptr)
  call destroy_ (this%prmrmtboxk(inch)%ptr)
  deallocate (this%prmrmtboxk1(inch)%ptr)
  deallocate (this%prmrmtboxk(inch)%ptr)
  this%prmrmtboxk1(inch)%ptr => null ()
  this%prmrmtboxk(inch)%ptr => null ()
 end do
 deallocate (this%prmrmtboxk1)
 deallocate (this%prmrmtboxk)
end if
!%------------------------------------------------------------
this%numrmrmtbox=numrmrmtbox
!%------------------------------------------------------------
if (this%numrmrmtbox>0) then
 allocate (this%prmrmtboxk1(this%numrmrmtbox))
 allocate (this%prmrmtboxk(this%numrmrmtbox))
 do inch=1,this%numrmrmtbox
   allocate (this%prmrmtboxk1(inch)%ptr)
   allocate (this%prmrmtboxk(inch)%ptr)
   call create_ (this%prmrmtboxk1(inch)%ptr)
   call create_ (this%prmrmtboxk(inch)%ptr)
   call set_ (this%prmrmtboxk1(inch)%ptr,iserror,chemsys=this%pchemsys(1)%ptr)
   if (iserror) goto 10
   call set_ (this%prmrmtboxk(inch)%ptr,iserror,chemsys=this%pchemsys(1)%ptr)
   if (iserror) goto 10 
 end do
end if
!%------------------------------------------------------------
return
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: create_rmrmtbox_ '
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
subroutine get_numunk_dsa_cheproo &
   (this, &
    nunkdsa, &
    iserror, &
	nunkdsanch)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the total number of unknowns for DSA approach.
! Also return the component per nodal chemistry  
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                           :: this       ! Type CHEPROO variable 

integer, intent(out)                                   :: nunkdsa    ! Total number of unknowns for DSA

logical, intent(out)                                   :: iserror    ! iserror=true, then there was an error 

integer, pointer, dimension(:),optional                :: nunkdsanch
 
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
 ithnch, &
 ncomp
character(len=100)                    :: &
 msg 
integer, pointer                      :: &
 id(:) => null ()
logical                               :: &
 havenunkdsanch
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Check type of reactive transport solver 
!%--------------------------------------------------------------
if (this%itypesolver/=dsasolver) then
  msg='Error, DSA solver is not performed in CHEPROO'
  goto 10
end if
!%------------------------------------------------------------
!% Check optional arguments
!%------------------------------------------------------------
havenunkdsanch=present(nunkdsanch)
!%------------------------------------------------------------
nunkdsa=0
call check_pointer_ (id,this%numnch,.true.)
!%------------------------------------------------------------
do ithnch=1,this%numnch
 call get_chem_info_ (this%pnodalchemk1(ithnch)%ptr,iserror,npri=ncomp)
 if (iserror) goto 10
 nunkdsa=nunkdsa+ncomp
 id(ithnch)=ncomp
end do
!%------------------------------------------------------------
!%If nunkdsanch is present then store the number of unknown 
!%in each nodal chemistry object. 
!%------------------------------------------------------------
if (havenunkdsanch) then
 call check_pointer_ (nunkdsanch,this%numnch,.true.)
 nunkdsanch=id 
end if  
!%------------------------------------------------------------
!%Deallocate pointers 
!%------------------------------------------------------------
call check_pointer_ (id,1,.false.)
!%------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_numunk_dsa_'
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
subroutine get_ncomp_dsa_cheproo &
   (this, &
    ncomp, &
	hashindex, &
	nnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the number of components for each nodal chemistry
! object and their corresponding hash index. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                           :: this       ! Type CHEPROO variable 

integer, intent(in)                                    :: nnch       ! Number of nodal chemistry objects 

integer, intent(inout), dimension(nnch)                :: ncomp      ! Number of components

integer, intent(inout), dimension(nnch)                :: hashindex  ! Hash indices of the components zone

logical, intent(out)                                   :: iserror    ! iserror=true, then there was an error 
 
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
 ithnch, &
 npri, &
 hash
character(len=100)                    :: &
 msg 
logical                               :: &
 havenunkdsanch
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Check the reactive transport solver
!%--------------------------------------------------------------
if (this%itypesolver/=dsasolver) then
  msg='Error, DSA solver is not performed in CHEPROO'
  goto 10
end if
!%--------------------------------------------------------------
!% Check the number of nodal chemistry objects 
!%--------------------------------------------------------------
if (nnch/=this%numnch) then
  msg='Error in number of nodal chemistry objects'
  goto 10
end if
!%------------------------------------------------------------
!% Loop on nodal chemistry objects 
!%------------------------------------------------------------
do ithnch=1,this%numnch
 call get_chem_info_ (this%pnodalchemk1(ithnch)%ptr,iserror,npri=npri,hashcompz=hash)
 if (iserror) goto 10
 ncomp(ithnch)=npri
 hashindex(ithnch)=hash
end do
!%------------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_ncomp_dsa_'
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
subroutine save_nodalchem_cheproo &
   (this, &
    isbackwards)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Save nodal chemistry 
! isbackwards=true => save nodal chemistry k+1 in nodal chemistry in k.
! isbackwards=false => save nodal chemistry k in nodal chemistry in k+1. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)       :: this

logical, intent(in)                   :: isbackwards 
 
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
 inch 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
if (this%numnch>0) then
!%------------------------------------------------------------
!% isbackward=true   then nodalchem(k)=nodalchem(k+1) 
!% isbackward=false  then rmrmtbox(k)=rmrmtbox(k+1)
!%------------------------------------------------------------
 if (isbackwards) then 
  do inch=1,this%numnch
   this%pnodalchemk(inch)%ptr=this%pnodalchemk1(inch)%ptr
  end do
  do inch=1,this%numrmrmtbox
   this%prmrmtboxk(inch)%ptr=this%prmrmtboxk1(inch)%ptr
  end do
 else
!%------------------------------------------------------------
!% isbackward=true   then nodalchem(k+1)=nodalchem(k) 
!% isbackward=false  then rmrmtbox(k+1)=rmrmtbox(k)
!%------------------------------------------------------------
  do inch=1,this%numnch
   this%pnodalchemk1(inch)%ptr=this%pnodalchemk(inch)%ptr
  end do
  do inch=1,this%numrmrmtbox
   this%prmrmtboxk1(inch)%ptr=this%prmrmtboxk(inch)%ptr
  end do
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
subroutine get_ujth_cmobith_cheproo &
   (this, &
    umob, &
    ncomp, &
    nmobph, &
    jthnch, &
	ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the mobile aqueous components in the aqueous, gas phases.
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this      ! Type CHEPROO variable 

real*8, pointer, dimension(:,:)        :: umob 

integer, intent(in)                    :: jthnch

integer, intent(in)                    :: ithnch

integer, intent(out)                   :: nmobph

integer, intent(out)                   :: ncomp

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
  
character(len=100)                    :: &
 msg
integer                               :: &
 hash, &
 i 
type(t_chemicalsystem), pointer       :: &
 jthchemsys => null ()
type(t_nodalchemistry), pointer       :: &
 ithnodalchem => null (), &
 jthnodalchem => null ()
real*8, pointer                       :: &
 vector(:) => null (), &
 cmob (:,:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
if ((ithnch<=0.or.ithnch>this%numnch) &
            .or. &
	(jthnch<=0.or.jthnch>this%numnch)) then
 msg='Error in nodal chemistry index'
 goto 10
end if
!%------------------------------------------------------------
ithnodalchem => this%pnodalchemk1(ithnch)%ptr
jthnodalchem => this%pnodalchemk1(jthnch)%ptr
!%------------------------------------------------------------
call get_chem_info_ (ithnodalchem,iserror,cmob=cmob)
if (iserror) goto 10
nmobph=size(cmob,2)
!%------------------------------------------------------------- 
call get_chem_info_(jthnodalchem,iserror,chemsys=jthchemsys, &
                    hashcompz=hash,npri=ncomp)
if (iserror) goto 10
!%------------------------------------------------------------- 
call check_pointer_ (umob,ncomp,nmobph,.true.)
!%------------------------------------------------------------- 
do i=1,nmobph 
 call make_lin_trf_ (jthchemsys,vector,cmob(:,i),hash,iserror)
 if (iserror) goto 20
 umob(:,i)=vector
end do
!%------------------------------------------------------------- 
20 continue 
!%------------------------------------------------------------- 
!% Deallocate local pointers 
!%------------------------------------------------------------- 
call check_pointer_ (vector,1,.false.)
call check_pointer_ (cmob,1,1,.false.)
jthchemsys => null ()
ithnodalchem => null ()
jthnodalchem => null ()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_ujth_cmobith_'
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
subroutine get_uaq_nch_sia_cheproo &
  (this, &
   uaq, &
   ndim, &
   iserror)              
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the aqueous components for all nodal chemistry 
! objects. Components are stored in a vector
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this      ! Type CHEPROO variable 

integer, intent(in)                   :: ndim

real*8, intent(out), dimension(ndim)  :: uaq

logical, intent(out)                  :: iserror   ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)            :: &
 msg
integer                       :: &
 ithnch, &
 ipos1, &
 ipos2, &
 npri, &
 ndim1
real*8, pointer               :: &
 caq(:) => null (), &
 uaq1(:) => null ()
type(t_nodalchemistry), pointer :: &
 pnch => null ()
type(t_chemicalsystem), pointer :: &
 pchemsys => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%---------------------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------------------
!% Check if there were created chemcial system objects 
!%---------------------------------------------------------------------
if (this%numchemsys==0) then 
 msg='Error, not defined chemical system objects'
 goto 10
end if 
!%---------------------------------------------------------------------
!% Get the number of primary species (used in sia)
!%---------------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numbase=npri)
if (iserror) goto 10 
ndim1=npri*this%numnch
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
if (ndim1/=ndim) then
 msg='Error in dimension of uaq vector'
 goto 10  
end if
!%---------------------------------------------------------------------
!% Initialice variables 
!%---------------------------------------------------------------------
ipos1=1
ipos2=0
uaq=0.0d0 
do ithnch=1,this%numnch
 
  pnch => this%pnodalchemk1(ithnch)%ptr 
  
  call get_chem_info_(pnch,iserror,caq=caq,chemsys=pchemsys)

  if (iserror) goto 20 
  
  call make_lin_trf_ (pchemsys,uaq1,caq,0,iserror)

  if (iserror) goto 20
  
  ipos2=ipos2+npri
  
  uaq(ipos1:ipos2)=uaq1
  
  ipos1=ipos2+1
 
end do
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------------
!% Nullify local pointers 
!%-------------------------------------------------------------------
pnch => null ()
pchemsys => null ()
call check_pointer_ (caq,1,.false.)
call check_pointer_ (uaq1,1,.false.)
if (iserror) goto 10 
!%-------------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_uaq_nch_sia_'
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
subroutine get_umob_water1_cheproo &
   (this, &
    umob, &
    ncomp, &
    nmobph, &
    namew, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return mobile components ofr ith water. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this

real*8, pointer, dimension(:,:)       :: umob 

character(len=*), intent(in)          :: namew

integer, intent(out)                  :: ncomp

integer, intent(out)                  :: nmobph

logical, intent(out)                  :: iserror    ! iserror=true, then there was an error 
 
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
 ithw, &
 i
character(len=100)                    :: &
 name, &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
ithw=0
do i=1,this%numw
 
 call get_chem_info_ (this%pwater(i)%ptr,iserror,name=name)
 if (namew==name) then
  ithw=i
  exit
 end if
 
end do
!%------------------------------------------------------------
if (ithw==0) then
 msg='Error, not defined in CHEPROO de water type'
 call add_ (msg,namew)
 goto 10
end if
!%------------------------------------------------------------
call get_chem_info_ &
   (this%pwater(ithw)%ptr, &
    iserror, &
    umob=umob, &
    npri=ncomp, &
    nummobph=nmobph)
 
if (iserror) then
 msg='Error when calling get_chem_info_'
 goto 10
end if
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_umob_water_'
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
subroutine get_umob_water2_cheproo &
   (this, &
    umob, &
    ncomp, &
    nmobph, &
    ithw, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return uaq components for any water defined in the chemical module
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this

real*8, pointer, dimension(:,:)       :: umob 

integer, intent(out)                  :: ncomp

integer, intent(out)                  :: nmobph

integer, intent(in)                   :: ithw

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
integer                               :: &
 i
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
if (ithw>=0.or.ithw>this%numw) then
 msg='Error in water index, index='
 call add_ (msg,ithw) 
 goto 10
end if
!%------------------------------------------------------------
call get_chem_info_ &
   (this%pwater(ithw)%ptr, &
    iserror, &
    umob=umob, &
    npri=ncomp, &
    nummobph=nmobph)
 
if (iserror) then
 msg='Error when calling get_chem_info_'
 goto 10
end if
 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_umob_water_'
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
subroutine get_uaq_water_cheproo &
   (this, &
    uaq, &
    npri, &
    ithw, &
    iserror, &
    ithnch)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return uaq components for any water defined in CHEPROO
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this

integer, intent(in)                   :: npri

real*8, intent(out), dimension(npri)  :: uaq

integer, intent(in)                   :: ithw

logical, intent(out)                  :: iserror    ! iserror=true, then there was an error 

integer, intent(in), optional         :: ithnch     ! Index of the nodal chemistry object that the user want to the U
 
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
 hash
real*8, pointer                       :: &
 vector(:) => null ()
 real*8, pointer                       :: &
 c(:) => null ()
character(len=100)                    :: &
 msg 
logical                               :: &
haveithnch
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check optional arguments
!%------------------------------------------------------------
haveithnch=present(ithnch)
!%-------------------------------------------------------------
if (ithw<=0.or.ithw>this%numw) then
 msg='Error in water index, index='
 call add_ (msg,ithw) 
 goto 10
end if
!%-------------------------------------------------------------
if (haveithnch) then
    if (ithnch<=0.or.ithnch>this%numnch) then
        msg='Error in nodal chemistry index, index='
        call add_ (msg,ithnch) 
        goto 10
    endif
end if 
!%------------------------------------------------------------
call get_chem_info_ (this%pwater(ithw)%ptr,iserror,c=c)
if (iserror) goto 20 
if (haveithnch) then
 call get_chem_info_ (this%pnodalchemk1(ithnch)%ptr,iserror,hashcompz=hash) !We get the hash for the U of the node
 if (iserror) goto 20 
else
 hash=0
end if 
call make_lin_trf_ (this%pchemsys(1)%ptr,vector,c,hash,iserror)  !We multiply the concentrations by U to get uaq
if (iserror) goto 20 
!%-------------------------------------------------------------
uaq=vector
!%-------------------------------------------------------------
20 continue 
call check_pointer_ (vector,1,.false.)
call check_pointer_ (c,1,.false.)
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_uaq_water_'
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
subroutine get_uaq_nch_cheproo &
   (this, &
    uaq, &
    npri, &
    ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return uaq components for any water defined in CHEPROO
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this

integer, intent(in)                   :: npri

real*8, intent(out), dimension(npri)  :: uaq

integer, intent(in)                   :: ithnch

logical, intent(out)                  :: iserror    ! iserror=true, then there was an error  
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer                       :: &
 vector(:) => null ()
 real*8, pointer                       :: &
 c(:) => null ()
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nch index, index='
 call add_ (msg,ithnch) 
 goto 10
end if
!%------------------------------------------------------------
call get_chem_info_ (this%pnodalchemk1(ithnch)%ptr,iserror,c=c)
if (iserror) goto 20 
call make_lin_trf_ (this%pchemsys(1)%ptr,vector,c,0,iserror)  !We multiply the concentrations by U to get uaq
if (iserror) goto 20 
!%-------------------------------------------------------------
uaq=vector
!%-------------------------------------------------------------
20 continue 
call check_pointer_ (vector,1,.false.)
call check_pointer_ (c,1,.false.)
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_uaq_nch_'
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
subroutine get_ph_nch_cheproo &
   (this, &
    ph, &
    ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return ph for ith nodal chemistry object 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this

real*8, intent(out)                   :: ph

integer, intent(in)                   :: ithnch

logical, intent(out)                  :: iserror    ! iserror=true, then there was an error  
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
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nch index, index='
 call add_ (msg,ithnch) 
 goto 10
end if
!%------------------------------------------------------------
call get_chem_info_ (this%pnodalchemk1(ithnch)%ptr,iserror,ph=ph)
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_ph_nch_'
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
subroutine get_uaq_cheproo &
   (this, &
    uaq, &
    npri, &
    ith, &
	jth, &
	isbackward, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return uaq[npri]=Uaqith*caqjth, where ith and jth are referring
! to nodal chemistry index for U aqand caq respectively. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this       ! Type CHEPROO variable

integer, intent(in)                   :: npri       ! Number of components

real*8, intent(inout), dimension(npri):: uaq        ! Aqueous components 

integer, intent(in)                   :: ith        ! Nodal chemistry index for U

integer, intent(in)                   :: jth        ! Nodal chemistry index for c

logical, intent(in)                   :: isbackward ! if true, then use the nodalchemistry list in k
                                                    ! if false, then use the nodalchemistry list in k+1

logical, intent(out)                  :: iserror    ! iserror=true, then there was an error 
 
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
 npri1, &
 hash
real*8, pointer                       :: &
 caq(:) => null (), &
 uaq1(:) => null ()
character(len=100)                    :: &
 msg 
type(t_nodalchemistry), pointer       :: &
 pnchith => null (), &
 pnchjth => null ()
type(t_chemicalsystem), pointer       :: &
 pchemsysith => null()
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check the ith and jth nodal chemistry indices
!%-------------------------------------------------------------
if (ith<=0.or.ith>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,ith) 
 goto 10
end if
if (jth<=0.or.jth>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,jth) 
 goto 10
end if
!%------------------------------------------------------------
!% Assign the corresponding nodal chemistries pointers
!%------------------------------------------------------------
if (isbackward) then
 pnchith => this%pnodalchemk(ith)%ptr
 pnchjth => this%pnodalchemk(jth)%ptr
else
 pnchith => this%pnodalchemk1(ith)%ptr
 pnchjth => this%pnodalchemk1(jth)%ptr
end if 
!%------------------------------------------------------------
!% Get components hash index of ith nodal chemistry object
!%------------------------------------------------------------
call get_chem_info_ (pnchith,iserror,hashcompz=hash,npri=npri1,chemsys=pchemsysith)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Check the number of components of the ith components zone
!%------------------------------------------------------------
if (npri/=npri1) then
 msg='Error in number of components'
 goto 20 
end if
!%------------------------------------------------------------
!% Get concentrations vector of jth nodal chemistry object 
!%------------------------------------------------------------
call get_chem_info_ (pnchjth,iserror,caq=caq)
if (iserror) goto 20 

call make_lin_trf_ (pchemsysith,uaq1,caq,hash,iserror)
if (iserror) goto 20 
!%-------------------------------------------------------------
!% Copy 
!%-------------------------------------------------------------
uaq=uaq1
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!% Deallocate local pointers
!%-------------------------------------------------------------
call check_pointer_ (caq,1,.false.)
call check_pointer_ (uaq1,1,.false.)
pnchith => null()
pnchjth => null ()
pchemsysith => null()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_uaq_'
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
subroutine get_uads_cheproo &
   (this, &
    uads, &
    npri, &
    ith, &
	jth, &
	isbackward, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return uaq[npri]=Uaqith*caqjth, where ith and jth are referring
! to nodal chemistry index for U aqand caq respectively. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this       ! Type CHEPROO variable

integer, intent(in)                   :: npri       ! Number of components

real*8, intent(inout), dimension(npri):: uads       ! Aqueous components on the surfaces 

integer, intent(in)                   :: ith        ! Nodal chemistry index for U

integer, intent(in)                   :: jth        ! Nodal chemistry index for c

logical, intent(in)                   :: isbackward ! if true, then use the nodalchemistry list in k
                                                    ! if false, then use the nodalchemistry list in k+1

logical, intent(out)                  :: iserror    ! iserror=true, then there was an error 
 
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
 npri1, &
 hashcompz
real*8, pointer                       :: &
 cads(:) => null (), &
 uads1(:) => null ()
character(len=100)                    :: &
 msg 
type(t_nodalchemistry), pointer       :: &
 pnchith => null (), &
 pnchjth => null ()
type(t_chemicalsystem), pointer       :: &
 pchemsysith => null()
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check the ith and jth nodal chemistry indices
!%-------------------------------------------------------------
if (ith<=0.or.ith>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,ith) 
 goto 10
end if
if (jth<=0.or.jth>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,jth) 
 goto 10
end if
!%------------------------------------------------------------
!% Assign the corresponding nodal chemistries pointers
!%------------------------------------------------------------
if (isbackward) then
 pnchith => this%pnodalchemk(ith)%ptr
 pnchjth => this%pnodalchemk(jth)%ptr
else
 pnchith => this%pnodalchemk1(ith)%ptr
 pnchjth => this%pnodalchemk1(jth)%ptr
end if 
!%------------------------------------------------------------
!% Get components hash index of ith nodal chemistry object
!%------------------------------------------------------------
call get_chem_info_ (pnchith,iserror,hashcompz=hashcompz,npri=npri1,chemsys=pchemsysith)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Check the number of components of the ith components zone
!%------------------------------------------------------------
if (npri/=npri1) then
 msg='Error in number of components'
 goto 20 
end if
!%------------------------------------------------------------
!% Get concentrations vector of jth nodal chemistry object 
!%------------------------------------------------------------
call get_chem_info_ (pnchjth,iserror,cads=cads)
if (iserror) goto 20 

call make_lin_trf_ (pchemsysith,uads1,cads,hashcompz,iserror)
if (iserror) goto 20 
!%-------------------------------------------------------------
!% Copy 
!%-------------------------------------------------------------
uads=uads1
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!% Deallocate local pointers
!%-------------------------------------------------------------
call check_pointer_ (cads,1,.false.)
call check_pointer_ (uads1,1,.false.)
pnchith => null()
pnchjth => null ()
pchemsysith => null()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_uaq_'
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
subroutine get_duaq_cheproo &
   (this, &
    duaq, &
    npriith, &
	nprijth, &
    ith, &
	jth, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return duaq[npri]=Uaqith*dcaqjth, where ith and jth are referring
! to nodal chemistry index for U aqand dcaq respectively. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                        :: this       ! Type CHEPROO variable

integer, intent(in)                                 :: npriith    ! Number of components

integer, intent(in)                                 :: nprijth    ! Number of components

real*8, intent(inout), dimension(npriith,nprijth)   :: duaq        ! Aqueous components 

integer, intent(in)                                 :: ith        ! Nodal chemistry index for U

integer, intent(in)                                 :: jth        ! Nodal chemistry index for c

logical, intent(out)                                :: iserror    ! iserror=true, then there was an error 
 
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
 npri, &
 msp, &
 hash, &
 nsp
real*8, pointer                       :: &
 dcaq(:,:) => null (), &
 duaq1(:,:) => null ()
character(len=100)                    :: &
 msg 
type(t_nodalchemistry), pointer       :: &
 pnchith => null (), &
 pnchjth => null ()
type(t_chemicalsystem), pointer       :: &
 pchemsysith => null()
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check the ith and jth nodal chemistry indices
!%-------------------------------------------------------------
if (ith<=0.or.ith>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,ith) 
 goto 10
end if
if (jth<=0.or.jth>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,jth) 
 goto 10
end if
!%------------------------------------------------------------
!% Assign the corresponding nodal chemistries pointers
!%------------------------------------------------------------
pnchith => this%pnodalchemk1(ith)%ptr
pnchjth => this%pnodalchemk1(jth)%ptr
!%------------------------------------------------------------
!% Get components hash index of ith nodal chemistry object
!%------------------------------------------------------------
call get_chem_info_ (pnchith,iserror,hashcompz=hash,npri=npri,chemsys=pchemsysith)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Check the number of components of the ith components zone
!%------------------------------------------------------------
if (npriith/=npri) then
 msg='Error in number of components'
 goto 20 
end if
!%------------------------------------------------------------
!% Get concentrations vector of jth nodal chemistry object 
!%------------------------------------------------------------
call get_chem_info_ (pnchjth,iserror,dcmob=dcaq,npri=npri)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Check the number of components of the ith components zone
!%------------------------------------------------------------
if (nprijth/=npri) then
 msg='Error in number of components'
 goto 20 
end if

!%------------------------------------------------------------
!% Ask for the number of species of the chemical system
!%------------------------------------------------------------
call get_chem_info_(pchemsysith,iserror,numsp=nsp)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Compute duaq=U*dcaq
!%------------------------------------------------------------
call make_lin_trf_ (pchemsysith,duaq1,dcaq(1:nsp,:),hash,iserror)
if (iserror) goto 20 
!%-------------------------------------------------------------
!% Copy 
!%-------------------------------------------------------------
duaq=duaq1
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!% Deallocate local pointers
!%-------------------------------------------------------------
call check_pointer_ (dcaq,1,1,.false.)
call check_pointer_ (duaq1,1,1,.false.)
pnchith => null()
pnchjth => null ()
pchemsysith => null()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_duaq_'
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
subroutine get_duads_cheproo &
   (this, &
    duads, &
    npriith, &
	nprijth, &
    ith, &
	jth, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return duaq[npri]=Uaqith*dcaqjth, where ith and jth are referring
! to nodal chemistry index for U aqand dcaq respectively. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                        :: this       ! Type CHEPROO variable

integer, intent(in)                                 :: npriith    ! Number of components

integer, intent(in)                                 :: nprijth    ! Number of components

real*8, intent(inout), dimension(npriith,nprijth)   :: duads      ! Derivatives of the aqueous components

integer, intent(in)                                 :: ith        ! Nodal chemistry index for U

integer, intent(in)                                 :: jth        ! Nodal chemistry index for c

logical, intent(out)                                :: iserror    ! iserror=true, then there was an error 
 
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
 npri, &
 hashcompz
real*8, pointer                       :: &
 dcads(:,:) => null (), &
 duads1(:,:) => null ()
character(len=100)                    :: &
 msg 
type(t_nodalchemistry), pointer       :: &
 pnchith => null (), &
 pnchjth => null ()
type(t_chemicalsystem), pointer       :: &
 pchemsysith => null()
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check the ith and jth nodal chemistry indices
!%-------------------------------------------------------------
if (ith<=0.or.ith>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,ith) 
 goto 10
end if
if (jth<=0.or.jth>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,jth) 
 goto 10
end if
!%------------------------------------------------------------
!% Assign the corresponding nodal chemistry pointers
!%------------------------------------------------------------
pnchith => this%pnodalchemk1(ith)%ptr
pnchjth => this%pnodalchemk1(jth)%ptr
!%------------------------------------------------------------
!% Get components hash index of ith nodal chemistry object
!%------------------------------------------------------------
call get_chem_info_ (pnchith,iserror,hashcompz=hashcompz,npri=npri,chemsys=pchemsysith)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Check the number of components of the ith components zone
!%------------------------------------------------------------
if (npriith/=npri) then
 msg='Error in number of components'
 goto 20 
end if
!%------------------------------------------------------------
!% Get concentrations vector of jth nodal chemistry object 
!%------------------------------------------------------------
call get_chem_info_ (pnchjth,iserror,dcads=dcads,npri=npri)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Check the number of components of the ith components zone
!%------------------------------------------------------------
if (nprijth/=npri) then
 msg='Error in number of components'
 goto 20 
end if 
!%------------------------------------------------------------
!% Compute duads=U*dcads
!%------------------------------------------------------------
call make_lin_trf_ (pchemsysith,duads1,dcads,hashcompz,iserror)
if (iserror) goto 20 
!%-------------------------------------------------------------
!% Copy 
!%-------------------------------------------------------------
duads=duads1
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!% Deallocate local pointers
!%-------------------------------------------------------------
call check_pointer_ (dcads,1,1,.false.)
call check_pointer_ (duads1,1,1,.false.)
pnchith => null()
pnchjth => null ()
pchemsysith => null()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_duads_'
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
subroutine get_dusktrk_cheproo &
   (this, &
    dusktrk, &
    npri, &
    theta, &
    ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return dusktrk[ncomp]=theta*U*dusktrk for ith nodal chemistry
! object 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                        :: this       ! Type CHEPROO variable

integer, intent(in)                                 :: npri       ! Number of components

real*8, intent(inout), dimension(npri,npri)         :: dusktrk    ! Aqueous components 

integer, intent(in)                                 :: ithnch     ! Nodal chemistry index  

real*8, intent(in)                                  :: theta      ! Temporal weight 

logical, intent(out)                                :: iserror    ! iserror=true, then there was an error 
 
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
 npri1, &
 hash
real*8, pointer                       :: &
 array1(:,:) => null (), &
 array2(:,:) => null ()
character(len=100)                    :: &
 msg 
type(t_nodalchemistry), pointer       :: &
 pnchith => null ()
type(t_chemicalsystem), pointer       :: &
 pchemsysith => null()
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check the nodal chemmistry index
!%-------------------------------------------------------------
if (ithnch<0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,ithnch) 
 goto 10
end if
!%------------------------------------------------------------
!% Assign the corresponding nodal chemistry pointer
!%------------------------------------------------------------
pnchith => this%pnodalchemk1(ithnch)%ptr
!%------------------------------------------------------------
!% Get components hash index of ith nodal chemistry object
!%------------------------------------------------------------
call get_chem_info_ (pnchith,iserror,hashcompz=hash,dsktrk=array1,chemsys=pchemsysith,npri=npri1)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Check the number of components of the ith components zone
!%------------------------------------------------------------
if (npri/=npri1) then
 msg='Error in number of components'
 goto 20 
end if
!%------------------------------------------------------------
!% Get concentrations vector of jth nodal chemistry object 
!%------------------------------------------------------------ 
call make_lin_trf_ (pchemsysith,array2,array1,hash,iserror)
if (iserror) goto 20 
!%-------------------------------------------------------------
!% Copy 
!%-------------------------------------------------------------
dusktrk=theta*array2
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!% Deallocate local pointers
!%-------------------------------------------------------------
call check_pointer_ (array1,1,1,.false.)
call check_pointer_ (array2,1,1,.false.)
pnchith => null()
pchemsysith => null()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_dusktrk_'
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
subroutine write_output1_cheproo &
   (this, &
    time, &
	isnch0, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Generate information about species in file (txt) 
! according namesp list
!
!   - Spatial distribution fo species
!   - Spatial distribtion of components
!   - Spatial distribution of saturation indices
!   - Spatial distribution of consumption/production rates
!   - Spatial distribution of mol of species 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this     ! Type CHEPROO variable. 

real*8, intent(in)                     :: time     ! Time increment 

logical, intent(in)                    :: isnch0   ! isnch0=true, write the list in 0  

logical, intent(out)                   :: iserror  ! iserror=true, there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 i, &
 ipos, &
 last
character(len=100)                     :: &
 raiz, &
 msg
logical                                :: &
 iswritename
character(len=100), pointer            :: &
 nameparam(:) => null ()
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
integer, parameter                     :: &
 nparam=10, &
 ioutconc=1000, &
 iouttot=1001, &
 ioutsat=1002, &
 ioutkin=1003, &
 ioutmol=1004, &
 ioutminarea=1005, &
 ioutact=1006, &
 ioutequil=1007, &
 ioutchemparam=1008, &
 ioutvolfracmin=1009 , &
 ioutkgr=1010
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=' '
!%----------------------------------------------------------
call check_pointer_ (nameparam,nparam,.true.)
nameparam(1)='pH'
nameparam(2)='ionic strength'
nameparam(3)='water activity'
nameparam(4)='mass of water'
nameparam(5)='temperature'
!nameparam(6)='xsalt'
nameparam(6)='volume'
nameparam(7)='liquid density'
nameparam(8)='liquid pressure'
nameparam(9)='gas pressure'
nameparam(10)='capilar factor'
!%----------------------------------------------------------
!% Write headers 
!%----------------------------------------------------------
if (this%output%numsp>0) then
 raiz='spatial_conc'
 call lastletter_ (ipos,raiz)
 open(unit=ioutconc,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
 write(ioutconc,*)'------------------------------------------------'
 write(ioutconc,*)'       Spatial distribution of species          '
 write(ioutconc,*)'                  [mol/kgw]                     '
 write(ioutconc,*)'------------------------------------------------'
 write(ioutconc,1)'Time:',time
 write(ioutconc,*)'------------------------------------------------'
end if 
!%----------------------------------------------------------
raiz='spatial_tot'
call lastletter_ (ipos,raiz)
open(unit=iouttot,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(iouttot,*)'-----------------------------------------------'
write(iouttot,*)'       Spatial distribution of components      '
write(iouttot,*)'                   [mol/kgw]                   '
write(iouttot,*)'-----------------------------------------------'
write(iouttot,1)'Time:',time
write(iouttot,*)'-----------------------------------------------'
!%----------------------------------------------------------
raiz='spatial_sat'
call lastletter_ (ipos,raiz)
open(unit=ioutsat,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(ioutsat,*)'------------------------------------------------'
write(ioutsat,*)'   Spatial distribution of saturation indices   '
write(ioutsat,*)'                 log(IAP/Ke)                    '
write(ioutsat,*)'------------------------------------------------'
write(ioutsat,1)'Time:',time
write(ioutsat,*)'------------------------------------------------' 
!%----------------------------------------------------------------------
raiz='spatial_kin'
call lastletter_ (ipos,raiz)
open(unit=ioutkin,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(ioutkin,*)'--------------------------------------------------'
write(ioutkin,*)'Consumption/production rates due kinetic reactions'
write(ioutkin,*)'                [mol/s]                           '
write(ioutkin,*)'--------------------------------------------------'
write(ioutkin,1)'Time:',time
write(ioutkin,*)'------------------------------------------------'
!%----------------------------------------------------------------------
raiz='spatial_mol'
call lastletter_ (ipos,raiz)
open(unit=ioutmol,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(ioutmol,*)'--------------------------------------------------'
write(ioutmol,*)'       Spatial distribution of mol of species     '
write(ioutmol,*)'                       [mol]                      '
write(ioutmol,*)'--------------------------------------------------'
write(ioutmol,1)'Time:',time
write(ioutmol,*)'------------------------------------------------'
!%----------------------------------------------------------------------
raiz='spatial_kgr'
call lastletter_ (ipos,raiz)
open(unit=ioutkgr,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(ioutkgr,*)'--------------------------------------------------'
write(ioutkgr,*)'       Spatial distribution of kgr of species     '
write(ioutkgr,*)'                       [kgr]                      '
write(ioutkgr,*)'--------------------------------------------------'
write(ioutkgr,1)'Time:',time
write(ioutkgr,*)'------------------------------------------------'
!%----------------------------------------------------------
raiz='spatial_min_area'
call lastletter_ (ipos,raiz)
open(unit=ioutminarea,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(ioutminarea,*)'---------------------------------------------------------'
write(ioutminarea,*)'   Spatial distribution of reactive surface of mineral   '
write(ioutminarea,*)'                  [m2/m3 nodal chemistry]                '
write(ioutminarea,*)'---------------------------------------------------------'
write(ioutminarea,1)'Time:',time
write(ioutminarea,*)'---------------------------------------------------------'
!%----------------------------------------------------------------------
raiz='spatial_act'
call lastletter_ (ipos,raiz)
open(unit=ioutact,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(ioutact,*)'--------------------------------------------------'
write(ioutact,*)'   Spatial distribution of species activities     '
write(ioutact,*)'                   [mol/kgw]                      '
write(ioutact,*)'--------------------------------------------------'
write(ioutact,1)'Time:',time
write(ioutact,*)'------------------------------------------------'
!%----------------------------------------------------------------------
raiz='spatial_equil'
call lastletter_ (ipos,raiz)
open(unit=ioutequil,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(ioutequil,*)'------------------------------------------------------'
write(ioutequil,*)'Consumption/production rates due equilibrium reactions'
write(ioutequil,*)'                [mol/s]                               '
write(ioutequil,*)'------------------------------------------------------'
write(ioutequil,1)'Time:',time
write(ioutequil,*)'----------------------------------------------------'
!%----------------------------------------------------------------------
raiz='spatial_chem_param'
call lastletter_ (ipos,raiz)
open(unit=ioutchemparam,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(ioutchemparam,*)'------------------------------------------------------'
write(ioutchemparam,*)'      Spatial distribution of chemical parameters     '
write(ioutchemparam,*)'------------------------------------------------------'
write(ioutchemparam,1)'Time:',time
write(ioutchemparam,*)'----------------------------------------------------'
!%----------------------------------------------------------------------
raiz='spatial_volfracmin'
call lastletter_ (ipos,raiz)
open(unit=ioutvolfracmin,file=raiz(1:ipos)//'_cheproo.dat',status='unknown')
write(ioutvolfracmin,*)'-------------------------------------------------------------'
write(ioutvolfracmin,*)'     Spatial distribution of volume fraction of minerals     '
write(ioutvolfracmin,*)'                       [m3 min/m3 rock]                      '
write(ioutvolfracmin,*)'-------------------------------------------------------------'
write(ioutvolfracmin,1)'Time:',time
write(ioutvolfracmin,*)'----------------------------------------------------'
!%----------------------------------------------------------------------
!% Start loop on nodal chemistry obejects  
!%----------------------------------------------------------------------
do i=1,this%numnch
 
  if (isnch0) then
    pnch => this%pnodalchem0(i)%ptr
  else
    pnch => this%pnodalchemk1(i)%ptr
  end if 
  
  if (i==1) then
   iswritename=.true.
  else
   iswritename=.false.
  end if
!%----------------------------------------------------------------------
!% Write concentrations 
!%----------------------------------------------------------------------
  if (this%output%numsp>0) then
   call write_sps_(pnch,ioutconc,'concentration',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
     msg='Error when calling write_sps_'
     call add_ (msg,i)
     goto 20
   end if
  end if 
!%----------------------------------------------------------
!% Write components 
!%----------------------------------------------------------
  call write_tot_ (pnch,iouttot,iswritename,iserror,integeritem=i)
  if (iserror) then
     msg='Error when calling write_tot_'
     call add_ (msg,i)
     goto 20
  end if
!%----------------------------------------------------------
!% Write saturation indices 
!%----------------------------------------------------------
   call write_sat_(pnch,ioutsat,i,iswritename,iserror)
   if (iserror) then
     msg='Error when calling write_sat_'
     call add_ (msg,i)
     goto 20
   end if
!%----------------------------------------------------------
!% Write changes in species due kinetic reactions
!%----------------------------------------------------------
   call write_sps_(pnch,ioutkin,'kinetic changes',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
      msg='Error when calling write_sps_'
      call add_ (msg,i)
      goto 20
   end if
!%----------------------------------------------------------
!% Write mol 
!%----------------------------------------------------------
   call write_sps_(pnch,ioutmol,'mol',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
     msg='Error when calling write_sps_'
     call add_ (msg,i)
     goto 20
   end if
!%----------------------------------------------------------
!% Write kgr  
!%----------------------------------------------------------
   !call write_sps_(pnch,ioutkgr,'kgr',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   !if (iserror) then
   !  msg='Error when calling write_sps_'
   !  call add_ (msg,i)
   !  goto 20
   !end if
!%----------------------------------------------------------
!% Write reactive surface of minerals
!%----------------------------------------------------------
   call write_min_area_ (pnch,ioutminarea,'',i,iswritename,.true.,iserror)
   if (iserror) then
     msg='Error when calling write_min_area_'
     call add_ (msg,i)
     goto 20
   end if
!%----------------------------------------------------------
!% Write the activities of the species
!%----------------------------------------------------------
   call write_sps_(pnch,ioutact,'activity',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
     msg='Error when calling write_sps_'
     call add_ (msg,i)
     goto 20
   end if
!%----------------------------------------------------------
!% Write changes in species due equilibrium reactions
!%----------------------------------------------------------
   call write_sps_(pnch,ioutequil,'equilibrium changes',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
      msg='Error when calling write_sps_'
      call add_ (msg,i)
      goto 20
   end if
!%----------------------------------------------------------
!% Write chemical parameters 
!%----------------------------------------------------------
   call write_chemical_parameters_ (pnch,ioutchemparam,nameparam,nparam,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
      msg='Error when calling write_chemical_parameters_'
      call add_ (msg,i)
      goto 20
   end if
!%----------------------------------------------------------
!% Write volumetric fraction of minerals
!%----------------------------------------------------------
   call write_volfracmin_ (pnch,ioutvolfracmin,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
      msg='Error when calling write_volfracmin_'
      call add_ (msg,i)
      goto 20
   end if
end do 
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
if (this%numrmrmtbox>0) then

if (this%output%numsp>0) then
 write(ioutconc,*)'---------------------------------------------------------'
 write(ioutconc,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
 write(ioutconc,*)'---------------------------------------------------------'
end if 
!%----------------------------------------------------------
write(iouttot,*)'---------------------------------------------------------'
write(iouttot,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
write(iouttot,*)'---------------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutkin,*)'---------------------------------------------------------'
write(ioutkin,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
write(ioutkin,*)'---------------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutmol,*)'---------------------------------------------------------'
write(ioutmol,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
write(ioutmol,*)'---------------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutkgr,*)'---------------------------------------------------------'
write(ioutkgr,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
write(ioutkgr,*)'---------------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutminarea,*)'---------------------------------------------------------'
write(ioutminarea,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
write(ioutminarea,*)'---------------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutact,*)'---------------------------------------------------------'
write(ioutact,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
write(ioutact,*)'---------------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutequil,*)'---------------------------------------------------------'
write(ioutequil,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
write(ioutequil,*)'---------------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutchemparam,*)'---------------------------------------------------------'
write(ioutchemparam,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
write(ioutchemparam,*)'---------------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutvolfracmin,*)'---------------------------------------------------------'
write(ioutvolfracmin,*)'Reactive Multi-Rate Mass transfer nodal chemistry objects'
write(ioutvolfracmin,*)'---------------------------------------------------------'

end if 
!%----------------------------------------------------------------------
!% Start loop on reactive multi-rate mass transfer nodal chemistry  
!%----------------------------------------------------------------------
do i=1,this%numrmrmtbox
 
  
  pnch => this%prmrmtboxk1(i)%ptr
  
  if (i==1) then
   iswritename=.true.
  else
   iswritename=.false.
  end if
!%----------------------------------------------------------------------
!% Write concentrations 
!%----------------------------------------------------------------------
  if (this%output%numsp>0) then
   call write_sps_(pnch,ioutconc,'concentration',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
     msg='Error when calling write_sps_'
     call add_ (msg,i)
     goto 20
   end if
  end if 
!%----------------------------------------------------------
!% Write components 
!%----------------------------------------------------------
  call write_tot_ (pnch,iouttot,iswritename,iserror,integeritem=i)
  if (iserror) then
     msg='Error when calling write_tot_'
     call add_ (msg,i)
     goto 20
  end if
!%----------------------------------------------------------
!% Write saturation indices 
!%----------------------------------------------------------
   call write_sat_(pnch,ioutsat,i,iswritename,iserror)
   if (iserror) then
     msg='Error when calling write_sat_'
     call add_ (msg,i)
     goto 20
   end if
!%----------------------------------------------------------
!% Write changes in species due kinetic reactions
!%----------------------------------------------------------
   call write_sps_(pnch,ioutkin,'kinetic changes',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
      msg='Error when calling write_sps_'
      call add_ (msg,i)
      goto 20
   end if
!%----------------------------------------------------------
!% Write mol 
!%----------------------------------------------------------
   call write_sps_(pnch,ioutmol,'mol',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
     msg='Error when calling write_sps_'
     call add_ (msg,i)
     goto 20
   end if
!%----------------------------------------------------------
!% Write kgr  
!%----------------------------------------------------------
   !call write_sps_(pnch,ioutkgr,'kgr',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   !if (iserror) then
   !  msg='Error when calling write_sps_'
   !  call add_ (msg,i)
   !  goto 20
   !end if
!%----------------------------------------------------------
!% Write reactive surface of minerals
!%----------------------------------------------------------
   call write_min_area_ (pnch,ioutminarea,'',i,iswritename,.true.,iserror)
   if (iserror) then
     msg='Error when calling write_min_area_'
     call add_ (msg,i)
     goto 20
   end if
!%----------------------------------------------------------
!% Write the activities of the species
!%----------------------------------------------------------
   call write_sps_(pnch,ioutact,'activity',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
     msg='Error when calling write_sps_'
     call add_ (msg,i)
     goto 20
   end if
!%----------------------------------------------------------
!% Write changes in species due equilibrium reactions
!%----------------------------------------------------------
   call write_sps_(pnch,ioutequil,'equilibrium changes',this%output%namesp,this%output%numsp,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
      msg='Error when calling write_sps_'
      call add_ (msg,i)
      goto 20
   end if
!%----------------------------------------------------------
!% Write chemical parameters 
!%----------------------------------------------------------
   call write_chemical_parameters_ (pnch,ioutchemparam,nameparam,nparam,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
      msg='Error when calling write_chemical_parameters_'
      call add_ (msg,i)
      goto 20
   end if
!%----------------------------------------------------------
!% Write volumetric fraction of minerals
!%----------------------------------------------------------
   call write_volfracmin_ (pnch,ioutvolfracmin,iserror,iswritename=iswritename,integeritem=i)
   if (iserror) then
      msg='Error when calling write_volfracmin_'
      call add_ (msg,i)
      goto 20
   end if

end do 
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
if (this%output%numsp>0) then
 write(ioutconc,*)'------------------------------------------------'
end if 
!%----------------------------------------------------------
write(iouttot,*)'-----------------------------------------------'
!%----------------------------------------------------------
write(ioutsat,*)'------------------------------------------------' 
!%----------------------------------------------------------------------
write(ioutkin,*)'------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutmol,*)'------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutkgr,*)'------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutminarea,*)'------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutact,*)'------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutequil,*)'------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutchemparam,*)'------------------------------------------------'
!%----------------------------------------------------------------------
write(ioutvolfracmin,*)'------------------------------------------------'
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
20 continue 
!%---------------------------------------------------------------
!% Nullify local pointers 
!%---------------------------------------------------------------
pnch => null()
call check_pointer_ (nameparam,1,.false.)
!%---------------------------------------------------------------
if (iserror) goto 10 
!%---------------------------------------------------------------
return
1 format(a6,e10.4)
2 format(a6,3x,<this%output%numsp>a12)
3 format(a6,3x,<this%output%numcomp>a12)
 
10 continue 
print *,'*********************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: write_'
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
subroutine write_output2_cheproo &
   (this, &
    ioutput, &
    istep, &
    itype, &
    iswritetitle, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write on output file in gid or meshplot format
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this

integer, intent(in)                    :: istep

integer, intent(in)                    :: ioutput

integer, intent(in)                    :: itype

logical, intent(out)                   :: iserror       ! iserror=true, then there was an error 

logical, intent(in)                    :: iswritetitle 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 inch, &
 isp, &
 iph, &
 imin, &
 nminsp
real*8                                 :: &
 omgwfree
real*8, pointer                        :: &
 res(:,:) => null (), &
 value1(:) => null (), &
 value2(:) => null ()
integer, pointer                       :: &
 idsp(:) => null (), &
 idminsp(:) => null ()
character(len=100)                     :: &
 msg
character(len=100), pointer            :: &
 name(:) => null ()
logical, pointer                       :: &
 beminsp(:) => null ()
integer, parameter                     :: &
 ascii_output=1, &
 gid_output=2
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=' '
!%----------------------------------------------------------
iph=0
!%----------------------------------------------------------
!% Check if there are chemical system objects 
!%----------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system objects in CHEPROO'
 goto 10
end if
!%----------------------------------------------------------
if (this%output%numsp>0 &
        .and. &
     this%numchemsys>0 &
        .and. &
      this%numnch>0) then
!%----------------------------------------------------------
 select case (itype)
!%----------------------------------------------------------
!% For GID output
!%----------------------------------------------------------
 case (gid_output)
    call get_chem_info_(this%pchemsys(1)%ptr,iserror,idminsp=idminsp)
    if (iserror) goto 20
    nminsp=size(idminsp)
    call check_pointer_ (res,this%numnch,this%output%numsp,.true.)
    call check_pointer_ (beminsp,this%output%numsp,.true.)
    call check_pointer_ (name,this%output%numsp,.true.)
    call check_pointer_(idsp,this%output%numsp,.true.)
    name=this%output%namesp
!%---------------
    do isp=1,this%output%numsp
      if (name(isp)=='h+'.or.name(isp)=='H+') then
          iph=isp
      end if
        call get_sp_index_(this%pchemsys(1)%ptr,name(isp),idsp(isp))
        do imin=1,nminsp
         if (idsp(isp)==idminsp(imin)) then
          beminsp(isp)=.true.
          exit
         end if
        end do
    end do
!%---------------
    do inch=1,this%numnch
      call get_chem_info_(this%pnodalchemk1(inch)%ptr,iserror, &
                          c=value1,g=value2,omgwfree=omgwfree)
      if (iserror) goto 20
      where (beminsp)
            res(inch,1:this%output%numsp)=value1(idsp)*omgwfree
      elsewhere
            res(inch,1:this%output%numsp)=value1(idsp)
      end where
      if (iph>0) then ! only por H+
        isp=idsp(iph)
        res(inch,iph)=-dlog10(value1(isp)*value2(isp))
      end if
    end do
!%---------------
    if (iph>0) name(iph)='pH'
!%---------------
    call write_resultsgid_ (istep,ioutput,name,res,this%numnch, &
                            this%output%numsp,iswritetitle)
20  continue 
!%----------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------
    call check_pointer_ (res,1,1,.false.)
    call check_pointer_ (idsp,1,.false.)
    call check_pointer_ (value1,1,.false.)
    call check_pointer_ (name,1,.false.)
    call check_pointer_ (value2,1,.false.)
    call check_pointer_ (idminsp,1,.false.)
    call check_pointer_ (beminsp,1,.false.)
    if (iserror) goto 10
 
 case default
   msg='Error, not implemented public service (only for GID)'
   goto 10
 end select
!%----------------------------------------------------------
end if
!%----------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'CHEPROO:'
print *,'Name:',this%name
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
subroutine write_output_tecplot_cheproo &
   (this, &
    ioutput, &
    time, &
	root, &
	x, &
	y, &
	z, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write on output file in tecplot format
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this

integer, intent(in)                    :: ioutput

real*8, intent(in)                     :: time 

character(len=*), intent(in)           :: root 

real*8, intent(in)                     :: x(this%numnch)

real*8, intent(in)                     :: y(this%numnch)

real*8, intent(in)                     :: z(this%numnch)

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
integer                                :: &
 inch, &
 nsp 
real*8, pointer                        :: &
 val(:,:) => null (), &
 c(:) => null ()
character(len=100)                     :: &
 msg, &
 head 
character(len=100), pointer            :: &
 namesp(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=' '
!%----------------------------------------------------------
!% Check if there are chemical system objects 
!%----------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system objects in CHEPROO'
 goto 10
end if
!%----------------------------------------------------------
call get_chem_info_(this%pchemsys(1)%ptr,iserror,namesp=namesp,numsp=nsp)
call check_pointer_ (val,this%numnch,nsp,.true.)
do inch=1,this%numnch
 call get_chem_info_(this%pnodalchemk1(inch)%ptr,iserror,c=c)
 if (iserror) goto 20 
 val(1:nsp,inch)=c  
end do
head='conentrations [mol/kgw]'
call write_resultstecplot_ (time,ioutput,head,root,namesp, &
                            val,this%numnch,nsp,x,y,z)

20  continue 
!%----------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------
call check_pointer_ (val,1,1,.false.)
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (c,1,.false.)
if (iserror) goto 10
!%----------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'CHEPROO:'
print *,'Name:',this%name
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
subroutine write_attributes_cheproo &
   (this, &
    ioutput, &
    time, &
	isopenfile, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write attributes of the CHEPROO object. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)        :: this        ! Type CHEPROO variable. 

integer, intent(in)                    :: ioutput     ! Output unit 

real*8, intent(in)                     :: time        ! Time increment  

logical, intent(in)                    :: isopenfile  ! If .true. then cheproo object open the ioutput file. 

logical, intent(out)                   :: iserror     ! iserror=true, there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License: Barcelona, 2007 
!
!-------------------------------------------------------------------------
integer                                :: &
 i
character(len=100)                     :: &
 msg, &
 file
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%----------------------------------------------------------
!% Open output file (optional)
!%----------------------------------------------------------
if (isopenfile) then
 file=''
 call add_(file,ioutput)
 file(6:23)='_cheproo_info.dat'
 open(unit=ioutput,file=file(1:23),status='unknown')
end if 
!%----------------------------------------------------------
!% Write CHEPROO signature
!%----------------------------------------------------------
call write_signature_ (this,ioutput)
!%----------------------------------------------------------
write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
write(ioutput,*) ' Time:', time
write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
write(ioutput,*) ' Reactive transport solver implemented: '
if (this%itypesolver==1) then
 write(ioutput,*) '  SNIA (Sequential Non-Iteration Approach)  '
else if(this%itypesolver==2) then
 write(ioutput,*) '  SIA (Sequential Iteration Approach)  '
 write(ioutput,1) 'Zero for unknown in SIA:',this%zero
else
 write(ioutput,*) '  DSA (Direct Substitution Approach)   '
 write(ioutput,1) 'Zero for unknown in DSA:',this%zero
end if
!%-----------------------------------------------------------
write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
!%----------------------------------------------------------
! 2) Write water objects
!%----------------------------------------------------------
if (this%numw>0) then
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*                 WATERS OBJECTS                *'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
do i=1, this%numw
 call write_ (this%pwater(i)%ptr,ioutput,iserror)
  if (iserror) then
       msg='Error when calling write_, water'
     call add_ (msg,i)
     goto 10
  end if
end do
write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
end if 
!%----------------------------------------------------------
! 3) Write nodal chemistry objects 
!%----------------------------------------------------------
if (this%numnch>0) then
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*               NODAL CHEMISTRY OBJECTS         *'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
do i=1, this%numnch
 call write_ (this%pnodalchemk1(i)%ptr,ioutput,iserror)
  if (iserror) then
       msg='Error when calling write_, nodal chemistry'
     call add_ (msg,i)
     goto 10
  end if
end do
write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
end if 
!%----------------------------------------------------------
! 4) Write nodal chemistry objects 
!%----------------------------------------------------------
if (this%numrmrmtbox>0) then
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*       Reactive Multi-Rate Mass Transfer Box   *'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
write(ioutput,*) '*************************************************'
do i=1, this%numrmrmtbox
 call write_ (this%prmrmtboxk1(i)%ptr,ioutput,iserror)
  if (iserror) then
       msg='Error when calling write_, nodal chemistry'
     call add_ (msg,i)
     goto 10
  end if
end do
write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
end if 
!%----------------------------------------------------------
! 5) Write chemical system objects (optional)
!%----------------------------------------------------------
if (iswchsys) then
 write(ioutput,*) '*************************************************'
 write(ioutput,*) '*************************************************'
 write(ioutput,*) '*************************************************'
 write(ioutput,*) '*          CHEMICAL SYSTEMS OBJECTS             *'
 write(ioutput,*) '*************************************************'
 write(ioutput,*) '*************************************************'
 write(ioutput,*) '*************************************************'
 call write_ (this%pchemsys(1)%ptr,ioutput,iserror)
 if (iserror) goto 10 
end if
!%-----------------------------------------------------------
!% Close output file 
!%----------------------------------------------------------
if (isopenfile) then
 close(unit=ioutput)
end if
!%----------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: write_'
print *, msg
print *,'**********************'
return
 
1 format(a25,e10.4)
2 format(a108)

end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_ith_nch1_cheproo &
   (this, &
    ithnch, &
    ioutput, &
	time, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write chemical information about ith nodal chemistry.
! If time argument is present then write concentrations vs. time. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)   :: this      ! Type CHEPROO variable.  

integer, intent(in)            :: ioutput   ! Output unit 

integer, intent(in)            :: ithnch    ! Index of nodal chemistry object. 

real*8, intent(in)             :: time      ! Time increment 

logical, intent(out)           :: iserror   ! iserror=true, there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 i
character(len=100)                     :: &
 msg 
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%----------------------------------------------------------
!% Check the nodal chemistry index.
!%----------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
  msg='Error in nodal chemistry index'
  goto 10
end if
!%----------------------------------------------------------
!% Assign nodal chemistry pointer 
!%----------------------------------------------------------
pnch => this%pnodalchemk1(ithnch)%ptr
!%----------------------------------------------------------
! Write complete chemical information for ith nodal chemistry.  
!%----------------------------------------------------------
write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
write(ioutput,*) 'CHEPROO Information'
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
write(ioutput,*) 'Time', time 


write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
write(ioutput,1) '==>', ithnch
write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
call write_ (pnch,ioutput,iserror)
if (iserror) then
   msg='Error when calling public service write_'
   call add_ (msg,ithnch)
   goto 20
end if
 
write(ioutput,2) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------' 
!%--------------------------------------------------------
20 continue 
pnch => null ()
if (iserror) goto 10 
!%--------------------------------------------------------
return
1 format(a3,i5)
2 format(a108)

10 continue 
print *,'**********************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: write_'
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
subroutine write_ith_nch2_cheproo &
   (this, &
    ithnch, &
    ioutput, &
    typeoutdata, &
	iswritehead, &
	time, &
    iserror, &
	isbackwards)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write chemical information about ith nodal chemistry.
! If time argument is present then write concentrations vs. time. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)   :: this        ! Type CHEPROO variable.  

integer, intent(in)            :: ioutput     ! Output unit 

character(len=*)               :: typeoutdata ! Type of output data

integer, intent(in)            :: ithnch      ! Index of nodal chemistry object. 

real*8, intent(in)             :: time        ! Time increment 

logical, intent(in)            :: iswritehead ! iswritehead=true, then write head

logical, intent(out)           :: iserror     ! iserror=true, there was an error

logical, intent(in), optional  :: isbackwards ! if true, then considering the list in k 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 ipos
character(len=100)                     :: &
 msg,  &
 head, &
 root
logical                                :: &
 haveisbackwards 
integer, parameter                     :: &
 nparam=9
character(len=100), pointer            :: &
 nameparam(:) => null ()
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%----------------------------------------------------------
!% Check the nodal chemistry index.
!%----------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
  msg='Error in nodal chemistry index'
  goto 10
end if
!%----------------------------------------------------------
!% Check optional arguments 
!%----------------------------------------------------------
haveisbackwards=present(isbackwards)
!%----------------------------------------------------------
!% Assign nodal chemistry pointer 
!%----------------------------------------------------------
if (haveisbackwards.and.isbackwards)  then
 pnch => this%pnodalchemk(ithnch)%ptr
else
 pnch => this%pnodalchemk1(ithnch)%ptr
end if
!%----------------------------------------------------------
!% Write concentrations 
!%----------------------------------------------------------
if (this%output%numsp>0) then
 
  select case (typeoutdata)
  
  case ('concentration')
  	
	root='temporal_conc'
    call lastletter_ (ipos,root)
    open(unit=ioutput,file=root(1:ipos)//'_cheproo.dat',status='unknown')
    if (iswritehead) then 
       write(ioutput,*)'----------------------------------------------------------'
	   head='Concentration of species vs. time in nodal chemistry'
	   call add_ (head,ithnch)
	   write(ioutput,*) head
	   write(ioutput,*)'                  [mol/kgw]                     '
       write(ioutput,*)'----------------------------------------------------------'
    end if 
	call write_sps_(pnch,ioutput,'concentration',this%output%namesp,this%output%numsp,iserror,iswritename=iswritehead,realitem=time)    
	if (iserror) goto 20

  case ('mol')
  	
	root='temporal_mol'
    call lastletter_ (ipos,root)
    open(unit=ioutput,file=root(1:ipos)//'_cheproo.dat',status='unknown')
    if (iswritehead) then 
       write(ioutput,*)'----------------------------------------------------------'
	   head='Mol of species vs. time in nodal chemistry'
	   call add_ (head,ithnch)
	   write(ioutput,*) head
	   write(ioutput,*)'----------------------------------------------------------'
    end if 
	call write_sps_(pnch,ioutput,'mol',this%output%namesp,this%output%numsp,iserror,iswritename=iswritehead,realitem=time)    
	if (iserror) goto 20

  case ('chemical parameters')
  	
	call check_pointer_ (nameparam,nparam,.true.)
    nameparam(1)='pH'
    nameparam(2)='ionic strength'
    nameparam(3)='water activity'
    nameparam(4)='mass of water'
    nameparam(5)='temperature'
	nameparam(6)='liquid density'
	nameparam(7)='liquid pressure'
    nameparam(8)='gas pressure'
	nameparam(9)='capilar factor' 
    
	
	root='temporal_chemparam'
    call lastletter_ (ipos,root)
    open(unit=ioutput,file=root(1:ipos)//'_cheproo.dat',status='unknown')
    if (iswritehead) then 
       write(ioutput,*)'----------------------------------------------------------'
	   head='Chemical paramaters vs. time in nodal chemistry'
	   call add_ (head,ithnch)
	   write(ioutput,*) head
	   write(ioutput,*)'----------------------------------------------------------'
    end if 
	call write_chemical_parameters_ (pnch,ioutput,nameparam,nparam,iserror,iswritename=iswritehead,realitem=time)    
	call check_pointer_ (nameparam,1,.false.)
	if (iserror) goto 20
    
  case ('activity') 

    root='temporal_act'
    call lastletter_ (ipos,root)
    open(unit=ioutput,file=root(1:ipos)//'_cheproo.dat',status='unknown')
    if (iswritehead) then 
       write(ioutput,*)'----------------------------------------------------------'
	   head='Activity of species vs. time in nodal chemistry'
	   call add_ (head,ithnch)
	   write(ioutput,*) head
	   write(ioutput,*)'                  [mol/kgw]                     '
       write(ioutput,*)'----------------------------------------------------------'
    end if
    call write_sps_(pnch,ioutput,'activity',this%output%namesp,this%output%numsp,iserror,iswritename=iswritehead,realitem=time)    
	if (iserror) goto 20
    

 case ('saturation') 

    root='temporal_sat'
    call lastletter_ (ipos,root)
    open(unit=ioutput,file=root(1:ipos)//'_cheproo.dat',status='unknown')
    if (iswritehead) then 
       write(ioutput,*)'----------------------------------------------------------'
	   head='Saturation of mineral species vs. time in nodal chemistry'
	   call add_ (head,ithnch)
	   write(ioutput,*) head
	   write(ioutput,*)'----------------------------------------------------------'
    end if 
	call write_sat_(pnch,ioutput,1,iswritehead,iserror)    
	if (iserror) goto 20

  case default
 
    head='Output data type not recognized:'
	call add_ (head,typeoutdata)
	print *,head  
    
  end select 
    
end if 
!%--------------------------------------------------------
20 continue 
pnch => null ()
if (iserror) goto 10 
!%--------------------------------------------------------
return
1 format(a3,i5)
2 format(a108)

10 continue 
print *,'**********************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: write_'
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
subroutine write_from_namesp1_cheproo &
   (this, &
    ioutput, &
    time, &
    namesp, &
    nsp, &
    iserror)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in output file the concentration of all nodes for numsp species.
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                    :: this      ! Type CHEPROO variable 

integer, intent(in)                             :: ioutput   ! Output unit 

integer, intent(in)                             :: nsp       ! Number of species. 

real*8, intent(in)                              :: time      ! Time increment 

character(len=*), intent(in), dimension(nsp)    :: namesp    ! Name of the species 

logical, intent(out)                            :: iserror   ! iserror=true, then there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 inch, &
 i 
logical                                :: &
 iswritename
character(len=100)                     :: &
 msg 
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%---------------------------------------------------------
if (nsp==0) return
!%---------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------
write(ioutput,3) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
write(ioutput,*) 'Concentration of species vs. node'
write(ioutput,3) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
write(ioutput,2) 'Time=',time
write(ioutput,3) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
write(ioutput,1) 'nch',(namesp(i),i=1,nsp)
!%---------------------------------------------------------
! Loop to nodal chemistry objects
!%---------------------------------------------------------
do inch=1, this%numnch
 pnch => this%pnodalchemk1(inch)%ptr
 if(inch==1) then
  iswritename=.true.
 else
  iswritename=.false.
 end if
!%---------------------------------------------------------
! 1) Write in txt ith nodal chemistry 
!%---------------------------------------------------------
 call write_sps_(pnch,ioutput,'concentration',namesp,nsp,iserror,iswritename=iswritename,integeritem=i)
 if (iserror) then
   msg='Error when calling write_sps_'
   call add_ (msg,inch)
   goto 20
 end if
end do
!%---------------------------------------------------------
! End loop to nodal chemistry objects
!%---------------------------------------------------------
write(ioutput,3) '----------------------------------------'// &
                 '----------------------------------------'// &
                 '--------------------------'
!%----------------------------------------------------------
20 continue
pnch => null ()
if (iserror) goto 10 
!%----------------------------------------------------------
return
 
10 continue 
print *,'*********************'
print *,'CHEPROO:'
print *,'Namë:',this%name
print *,'Service:  write_'
print *, msg
print *,'*********************'
iserror=.true.
return
 
1 format(a5,<nsp>a15)
2 format(a5,f10.1)
3 format(a108)
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_from_namesp2_cheproo &
   (this, &
    ioutput, &
    time, &
    namesp, &
    nsp, &
    ithnch, &
    iswritehead, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in txt file the concentrations of ith nodal chemistry object 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this

integer, intent(in)                    :: ioutput

integer, intent(in)                    :: nsp

integer, intent(in)                    :: ithnch

real*8, intent(in)                     :: time          ! Time increment    

character(len=*), intent(in)           :: namesp(nsp)

logical, intent(in)                    :: iswritehead 

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
character(len=100)                     :: &
 msg, &
 head  
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
  msg='Error in nodal chemistry index, index='
  call add_ (msg,ithnch)
  goto 10
end if
!%----------------------------------------------------------
if (nsp==0) return
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
head='Concentration of species vs. in nodal chemistry'
call add_ (head,ithnch)
!%----------------------------------------------------------
pnch => this%pnodalchemk1(ithnch)%ptr
 
call write_sps_(pnch,ioutput,'concentration',namesp,nsp,iserror,iswritename=iswritehead,realitem=time)
nullify(pnch)

if (iserror) then
  msg='Error when calling write_sps_'
  call add_ (msg,ithnch)
  goto 10
end if
!%-------------------------------------------------------------
 
return
 
 
10 continue 
print *,'*********************'
print *,'CHEPROO:'
print *,'Namë:',this%name
print *,'Service:  write_'
print *, msg
print *,'********************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_ujth_cadsith_cheproo &
   (this, &
    uads, &
    ncomp, &
    jthnch, &
	ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the mobile aqueous components in the aqueous, gas phases.
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this

real*8, pointer, dimension(:)          :: uads

integer, intent(in)                    :: jthnch

integer, intent(in)                    :: ithnch

integer, intent(out)                   :: ncomp

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
  
character(len=100)                    :: &
 msg
integer                               :: &
 hash, &
 i 
type(t_chemicalsystem), pointer       :: &
 jthchemsys => null ()
type(t_nodalchemistry), pointer       :: &
 ithnodalchem => null (), &
 jthnodalchem => null ()
real*8, pointer                       :: &
 cads(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
if ((ithnch<=0.or.ithnch>this%numnch) &
            .or. &
	(jthnch<=0.or.jthnch>this%numnch)) then
 msg='Error in nodal chemistry index'
 goto 10
end if
!%------------------------------------------------------------
ithnodalchem => this%pnodalchemk1(ithnch)%ptr
jthnodalchem => this%pnodalchemk1(jthnch)%ptr
!%------------------------------------------------------------
call get_chem_info_ (ithnodalchem,iserror,cads=cads)
if (iserror) goto 10
!%------------------------------------------------------------- 
call get_chem_info_(jthnodalchem,iserror,chemsys=jthchemsys, &
                    hashcompz=hash,npri=ncomp)
if (iserror) goto 10
!%------------------------------------------------------------- 
!% Compute ud=Ud*cd
!%------------------------------------------------------------- 
call make_lin_trf_ (jthchemsys,uads,cads,hash,iserror)
!%------------------------------------------------------------- 
20 continue 
!%------------------------------------------------------------- 
!% Deallocate and nullify local pointers 
!%------------------------------------------------------------- 
call check_pointer_ (cads,1,.false.)
jthchemsys => null ()
ithnodalchem => null ()
jthnodalchem => null ()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'*****************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_ujth_cadsith_'
print *,msg
print *,'*****************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_ipos_unk_cheproo &
   (this, &
    iposunk, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the global indice for unknowns for DSA
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this      ! Type CHEPROO variable. 

integer, pointer, dimension(:,:)       :: iposunk   ! Ïndices of the unknowns for DSA [numnch,2]

logical, intent(out)                   :: iserror   ! iserror=true, then there was an error.
 
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
 inch, &
 ncomp
character(len=100)                     :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
call check_pointer_ (iposunk,this%numnch,2,.true.)
!%---------------------------------------------------------
!% Initialice variable
!%---------------------------------------------------------
iposunk(1,1)=1
iposunk(1,2)=0
!%---------------------------------------------------------
do inch=1,this%numnch
 call get_chem_info_ (this%pnodalchemk1(inch)%ptr,iserror,npri=ncomp)
 if (iserror) goto 10
 iposunk(inch,2) = iposunk(inch,2) + ncomp
 iposunk(inch,1) = iposunk(inch,2) + 1
end do
!%------------------------------------------------------------
return
 
10 continue 
print *,'***************************'
print *,'CHEPROO:'
print *,'Namë:',this%name
print *,'Service: get_ipos_unk_'
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
subroutine get_usktrk_cheproo &
   (this, &
    usktrk, &
    npri, &
	theta, &
    ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return usktrk[ncomp]=Uith*sktrkjth
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this     ! Type CHEPROO variable 

integer, intent(in)                    :: npri     ! Number of primary species  

real*8, intent(inout), dimension(npri) :: usktrk   ! usktrk[npri]

integer, intent(in)                    :: ithnch   ! Nodal chemistry index

real*8, intent(in)                     :: theta    ! Temporal weight

logical, intent(out)                   :: iserror  ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
 msg 
integer                                :: &
 npri1
type(t_nodalchemistry), pointer        :: &
 pnchk1 => null (),&
 pnchk => null ()
real*8, pointer,dimension(:)                     :: &
 vector => null()
!-------------------------------------------------------------------------
!
!   $code
!
 
!%-----------------------------------------------------------
msg=''
iserror=.false.
!%-----------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%------------------------------------------------------------
!% Assign nodal chemistry pointer
!%------------------------------------------------------------
pnchk1 => this%pnodalchemk1(ithnch)%ptr
pnchk => this%pnodalchemk(ithnch)%ptr
call get_chem_info_ (pnchk1,iserror,usktrk=vector,npri=npri1)
if (iserror) goto 10
!%------------------------------------------------------------
!% Check the number of components 
!%------------------------------------------------------------
if (npri/=npri1) then
 msg='Error in number of components'
 goto 20
end if
!%------------------------------------------------------------
!% usktrktheta=theta*usktrk1  
!%------------------------------------------------------------
usktrk=theta*vector
!%------------------------------------------------------------
!% Now in k 
!%------------------------------------------------------------
call get_chem_info_ (pnchk,iserror,usktrk=vector,npri=npri1)
if (iserror) goto 10
!%------------------------------------------------------------
!% Check the number of components 
!%------------------------------------------------------------
if (npri/=npri1) then
 msg='Error in number of components'
 goto 20
end if
!%------------------------------------------------------------
!% usktrktheta = usktrktheta + (1-theta) usktrk
!%------------------------------------------------------------
usktrk=usktrk+(1.0d0-theta)*vector
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointer 
!%------------------------------------------------------------
call check_pointer_ (vector,1,.false.)
pnchk1 => null()
pnchk => null()
if (iserror) goto 10
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_usktrk_ith_'
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
subroutine get_dumob_ith_cheproo &
   (this, &
    dumob, &
    nrow, &
    ncol, &
    nmobph, &
    ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return derivate of mobile components. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                :: this

real*8, pointer                             :: dumob (:,:)

integer, intent(in)                         :: ithnch

integer, intent(out)                        :: nrow

integer, intent(out)                        :: ncol

integer, intent(out)                        :: nmobph

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
character(len=100)                     :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%----------------------------------------------------------
call get_dumob_ &
  (this%pnodalchemk1(ithnch)%ptr, &
   dumob, &
   nrow, &
   ncol, &
   nmobph, &
   iserror)
 
if (iserror) then
 msg='Error when calling get_dumob_'
 call add_ (msg,ithnch)
 goto 10
end if
 
!%----------------------------------------------------------
return
 
10 continue 
print *,'***************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_dumob_ith_'
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
subroutine get_dumob_ith_jth_cheproo &
   (this, &
    dumob, &
    nrow, &
    ncol, &
    nmobph, &
    ithnch, &
    jthnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: This subroutine compute U_jth d cmob_ith
!                                               ---------
!                                               d c11_ith
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                       :: this

real*8, pointer, dimension(:,:)                    :: dumob 

integer, intent(in)                                :: ithnch

integer, intent(in)                                :: jthnch

integer, intent(out)                               :: nrow

integer, intent(out)                               :: ncol

integer, intent(out)                               :: nmobph

logical, intent(out)                               :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
 msg 
type(t_chemicalsystem), pointer        :: &
 pchemsys

!-------------------------------------------------------------------------
!
!   $code
!
!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
nrow=0
ncol=0
nmobph=0
!%----------------------------------------------------------
!% Check the ith nodal chemistry index 
!%----------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index:'
 call add_ (msg,ithnch)
 goto 10
end if
!%----------------------------------------------------------
!% Check the jth nodal chemistry index 
!%----------------------------------------------------------
if (jthnch<=0.or.jthnch>this%numnch) then
 msg='Error in nodal chemistry index:'
 call add_ (msg,jthnch)
 goto 10
end if
!%----------------------------------------------------------
call get_chem_info_ (this%pnodalchemk1(ithnch)%ptr,iserror, &
                     chemsys=pchemsys)
if (iserror) goto 10 
!%----------------------------------------------------------

!%----------------------------------------------------------
pchemsys => null ()
!%----------------------------------------------------------
return
 
10 continue 
print *,'***************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_dumob_ith_'
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
subroutine equilibrate_ith_cheproo &
   (this, &
    ithnch, &
    iserror, &
    ioutput, &
    ithminset, &
    ithgasset)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: isequilibrate ith nodal chemistry with ith mineral set. 
!
!   $Arguments:
! 
 
type (t_cheproo), intent(inout)        :: this

integer, intent(in)                    :: ithnch

integer, intent(in), optional          :: ioutput

integer, intent(in), optional          :: ithminset

integer, intent(in), optional          :: ithgasset

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
integer                                :: &
 imin
character(len=100)                     :: &
 msg, &
 nameout
type(t_nodalchemistry), pointer        :: &
 pnch
logical                                :: &
 haveithminset, &
 haveithgasset, &
 haveioutput 
character(len=100), pointer            :: &
 name => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=' '
!%--------------------------------------------------------------
!% Check optional arguments 
!%--------------------------------------------------------------
haveithminset=present(ithminset)
haveithgasset=present(ithgasset)
haveioutput=present(ioutput)
!%---------------------------------------------------------------
! If not present optional arguments izmin or izads the error
!%---------------------------------------------------------------
if (.not.haveithminset.and..not.haveithgasset) then
 msg='Error, not defined set'
 goto 10
end if
!%----------------------------------------------------------Check
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%---------------------------------------------------------------
! Open the ioutput unit 
!%---------------------------------------------------------------
if (haveioutput) then
 nameout=''
 call add_ (nameout,ioutput)
 nameout(6:21)='_equilibrate.dat'
 open(unit=ioutput,file=nameout(1:21),status='unknown')
 call write_signature_ (this,ioutput)
end if 
!%---------------------------------------------------------------
! Pointer nodal chemitry ith
!%---------------------------------------------------------------
pnch => this%pnodalchemk1(ithnch)%ptr
!%---------------------------------------------------------------
! Add mineral set
!%---------------------------------------------------------------
if (haveithminset) then
!%-------
 if (ithminset>this%numminset) then
  msg='Error in index of mineral set'
  goto 10
 end if
!%-------
 
!%--------------------------------------------------
! Use equilibrate service in nodalchemistry object 
!%--------------------------------------------------
 
 do imin=1,this%minset(ithminset)%nummin
 
  name =>  this%minset(ithminset)%namemin(imin)
  
  call equilibrate_(pnch,name,iserror,ioutput)
  
  if (iserror) then
    msg='Error when calling equilibrate_'
    goto 10
  end if
 
 end do

end if
!%--------------------------------------------------
if (haveioutput) then
 close (unit=ioutput)
end if
!%--------------------------------------------------
nullify(pnch)
name => null ()
!%--------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: equilibrate_ith_'
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
subroutine reaction_path_ith_cheproo &
   (this, &
    ithnch, &
    nstep, &
    iserror, &
    ioutput, &
    ithminset)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)        :: this        ! Type CHEPROO variable 

integer, intent(in)                    :: nstep       ! Number of steps 

integer, intent(in)                    :: ithnch      ! Index of nodal chemistry object 

integer, intent(in), optional          :: ioutput     ! Output unit 

integer, intent(in), optional          :: ithminset   ! Mineral set index 

logical, intent(out)                   :: iserror     ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
 msg, &
 nameout
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
logical                                :: &
 haveithminset, &
 haveioutput
integer                                :: &
 imin 
real*8, pointer                        :: &
 c => null ()
character(len=100), pointer            :: &
 namespout(:) => null ()
character(len=100), pointer            :: &
 name => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%--------------------------------------------------------------
iserror=.false.
msg=' '
!%--------------------------------------------------------------
haveithminset=present(ithminset)
haveioutput=present(ioutput)
!%--------------------------------------------------------------
!%--------------If not present optional arguments izmin or izads
!%--------------the error
!%--------------------------------------------------------------
if (.not.haveithminset) then
 msg='Error, not defined sets for to make the reaction path'
 goto 10
end if
!%----------------------------------------------------------Check
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%----------------------------------------------------------------
if (haveioutput) then
 nameout=''
 call add_ (nameout,ioutput)
 nameout(6:24)='_reaction_path.dat'
 open(unit=ioutput,file=nameout(1:24),status='unknown')
 call write_signature_ (this,ioutput)
end if
!%-------------------------------------Pointer nodal chemitry ith
pnch => this%pnodalchemk1(ithnch)%ptr
!%---------------------------------------------------------------
!% Add mineral set
!%---------------------------------------------------------------
if (haveithminset) then
 
 if (ithminset>this%numminset) then
  msg='Error, mineral set index greater than total number of sets'
  goto 10
 end if
 
 
 do imin=1,this%minset(ithminset)%nummin
 
  name => this%minset(ithminset)%namemin(imin)
  c => this%minset(ithminset)%cmin(imin)
  namespout => this%output%namesp
 
  call reaction_path_ (pnch,name,c,1.0d0,nstep,iserror, &
                       ioutput=ioutput,namespout=namespout)
 
  if (iserror) then
    msg='Error when calling reaction_path_'
    goto 10
  end if
 
 end do
end if
!%------------------------------------------------------------
if (haveioutput) then
 close (unit=ioutput)
end if
!%------------------------------------------------------------
pnch => null ()
name => null ()
c => null ()
namespout => null ()
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: reaction_path_ith_'
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
subroutine reaction_path_kinetic_ith_cheproo &
   (this, &
    ithnch, &
    nstep, &
	time, &
	theta, &
	isconvergence, &
    iserror, &
    ioutput)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Make the kinetic reaction path for ith nodal chemistry object
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)        :: this             ! Type CHEPROO variable 

integer, intent(in)                    :: nstep            ! Number of steps 

real*8, intent(in)                     :: time             ! Total time 

integer, intent(in)                    :: ithnch           ! Nodal chemistry index 

real*8, intent(in)                     :: theta            ! Temporal weight 

integer, intent(in), optional          :: ioutput          ! Output unit 

logical, intent(out)                   :: iserror          ! iserror=true, then there was an error

logical, intent(out)                   :: isconvergence    ! isconvergence=true, there was convergence in mol balance
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
 msg, &
 nameout 
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
logical                                :: &
 haveioutput
!-------------------------------------------------------------------------
!
!   $code
!
!%--------------------------------------------------------------
iserror=.false.
msg=' '
!%--------------------------------------------------------------
!% Check optional arguments 
!%--------------------------------------------------------------
haveioutput=present(ioutput)
!%--------------------------------------------------------------
!% Check the nodal chemistry indice 
!%--------------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%----------------------------------------------------------------
!% Open output unit 
!%----------------------------------------------------------------
if (haveioutput) then
 nameout=''
 call add_ (nameout,ioutput)
 nameout(6:27)='_reaction_kinetic.dat'
 open(unit=ioutput,file=nameout(1:27),status='unknown')
 call write_signature_ (this,ioutput)
end if
!%----------------------------------------------------------------
!% Pointer nodal chemitry ith
!%----------------------------------------------------------------
pnch => this%pnodalchemk1(ithnch)%ptr
!%---------------------------------------------------------------
!% Call the corresponding service in the nodal chemistry object 
!%---------------------------------------------------------------
call reaction_path_ (pnch,time,nstep,theta,isconvergence,iserror, &
                     ioutput=ioutput,namespout=this%output%namesp)
if (iserror) goto 10 
!%---------------------------------------------------------------
!% Close the output unit 
!%---------------------------------------------------------------
if (haveioutput) then 
 close(unit=ioutput)
end if 
!%---------------------------------------------------------------
!% Nullify local pointers 
!%---------------------------------------------------------------
nullify(pnch)
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: reaction_path_ith_'
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
subroutine reaction_path_ev_ith_cheproo &
   (this, &
    ithnch, &
    ntime, &
    nstep, &
    isequilibrate, &
    fraccryst, &
    iswcompbal, &
    iserror, &
    ioutput, &
    temprank)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Evaporate/dilute the ith nodal chemistry object. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)        :: this             ! Type CHEPROO variable 

integer, intent(in)                    :: nstep            ! Number of steps 

integer, intent(in)                    :: ithnch           ! Index of nodal chemistry object 

logical, intent(out)                   :: iserror          ! iserror=true, then there was an error

logical, intent(in)                    :: isequilibrate    ! isequilibrate=true, the solution is equilibrated according to 
                                                           ! mineral and gas phases. 

real*8, intent(in)                     :: fraccryst        ! Crystalline fractionation factor (0-1)  

real*8, intent(in)                     :: ntime            ! Evaporation/dilution times 

logical, intent(in)                    :: iswcompbal       ! iswcompbal=true, the mass balance of water component is calculated 

integer, intent(in), optional          :: ioutput          ! Unit for output 

real*8, intent(in), optional           :: temprank(2)      ! Temperature range  
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
type (t_nodalchemistry), pointer       :: &
 pnch => null ()
character(len=100)                     :: &
 msg, &
 nameout
integer                                :: &
 nsp
logical                                :: &
 haveioutput 
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=' '
!%--------------------------------------------------------------
!% Check optional arguments
!%--------------------------------------------------------------
haveioutput=present(ioutput)
!%--------------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%--------------------------------------------------------------
!% Open write unit 
!%--------------------------------------------------------------
if (haveioutput) then
 nameout=''
 call add_ (nameout,ioutput)
 nameout(6:31)='_reaction_evaporation.dat'
 open(unit=ioutput,file=nameout(1:31),status='unknown')
 call write_signature_ (this,ioutput)
end if 
!%--------------------------------------------------------------
pnch => this%pnodalchemk1(ithnch)%ptr
!%--------------------------------------------------------------
!call reaction_path_ &
!  (pnch, &
!   ntime, &
!   nstep, &
!   isequilibrate, &
!   fraccryst, &
!   iswcompbal, &
!   iserror, &
!   ioutput=ioutput, &
!   temprank=temprank, &
!   namespout=this%output%namesp)
!%------------------------------------------------------------
!% Close write unit 
!%------------------------------------------------------------
if (haveioutput) then
 close (unit=ioutput)
end if
!%------------------------------------------------------------
!% Nullify local nodal chemistry pointer 
!%------------------------------------------------------------
pnch => null ()
!%------------------------------------------------------------
if (iserror) goto 10
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: reaction_path_ith_'
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
subroutine get_numsp_cheproo &
   (this, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the number of species 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this

integer, intent(out)                   :: nsp

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
character(len=100)                     :: &
 msg 
integer                                :: &
 nsp1 
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=' '
!%--------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numsp=nsp)
if (iserror) goto 10 
!%--------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_numsp_'
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
subroutine get_namesp_cheproo &
   (this, &
    namesp, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of the species 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                   :: this

integer, intent(in)                            :: nsp

character(len=*), intent(out), dimension(nsp)  :: namesp

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
character(len=100)                     :: &
 msg 
integer                                :: &
 nsp1 
character(len=100), pointer            :: &
 namesp1(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=' '
!%--------------------------------------------------------------
!% Initializing variables
!%--------------------------------------------------------------
namesp='' 
!%--------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,namesp=namesp1)
if (iserror) goto 10 
nsp1=size(namesp1)
!%--------------------------------------------------------------
!% Check the number of species 
!%--------------------------------------------------------------
if (nsp1/=nsp) then
 msg='Error in number of species'
 goto 20
end if 
namesp=namesp1
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
20 continue 
call check_pointer_ (namesp1,1,.false.)
if (iserror) goto 10 
!%--------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_namesp_'
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
subroutine get_c_ithnch_cheproo &
   (this, &
    c, &
    nsp, &
    ithnch, &
    iserror,&
    isbackward)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute reaction term for SIA (Sequiential iteration Approach) for ith nodal chemistry. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this

integer, intent(in)                    :: nsp

integer, intent(in)                    :: ithnch ! Index of nodal chemistry object 

real*8, intent(out), dimension(nsp)    :: c

logical, intent(out)                   :: iserror ! If .true. there was error. 

logical, intent(in),optional           :: isbackward !if true we give the concentration of nodal
                                                     !chemistry in time k  
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
 msg 
integer                                :: &
 nsp1 
real*8, pointer                        :: &
 c1(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=' '
!%--------------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%---------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numsp=nsp1)
if (iserror) goto 10 
!%---------------------------------------------------------------
if (nsp1/=nsp) then
 msg='Error in number of species'
 goto 10
end if 
!%---------------------------------------------------------------
if (present(isbackward)) then
    if (isbackward) then
        call get_chem_info_ (this%pnodalchemk(ithnch)%ptr,iserror,c=c1)
    else
        call get_chem_info_ (this%pnodalchemk1(ithnch)%ptr,iserror,c=c1)
    endif
else
    call get_chem_info_ (this%pnodalchemk1(ithnch)%ptr,iserror,c=c1)
endif
c=c1
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
20 continue 
call check_pointer_ (c1,1,.false.)
if (iserror) goto 10 
!%--------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_c_ithnch_'
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
subroutine update_r_sia_ith_cheproo &
   (this, &
    utra, &
    ncomp, &
    dtime, &
    theta, &
    isconvchem, &
    ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute reaction term for SIA (Sequiential iteration Approach) for ith nodal chemistry. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)        :: this

integer, intent(in)                    :: ncomp

real*8, intent(inout)                  :: utra(ncomp)

real*8, intent(in)                     :: dtime ! Time increment 

real*8, intent(in)                     :: theta ! Temporal weight 

integer, intent(in)                    :: ithnch

logical, intent(out)                   :: iserror ! If .true. there was error. 

logical, intent(out)                   :: isconvchem ! If .true. there was convergence in chemistry step. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
 msg 
type(t_nodalchemistry), pointer        :: &
 pnchk => null (), &
 pnchk1 => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=' '
!%--------------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%---------------------------------------------------------------
pnchk => this%pnodalchemk(ithnch)%ptr
pnchk1 => this%pnodalchemk1(ithnch)%ptr
!%---------------------------------------------------------------
call update_r_sia_ (pnchk1,pnchk,utra,ncomp,dtime,theta,isconvchem,iserror)
pnchk => null ()
pnchk1 => null ()
if (iserror) then
  msg='Error when calling update_r_sia_'
  call add_ (msg,ithnch)
  goto 10
end if
!%--------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: update_r_sia_ith_'
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
subroutine build_trp_eq_sia_cheproo &
   (this, &
    trpmatrix, &
    b, &
    nnch, &
    nmobph, &
    mxconn, &
    nband, &
    typesto, &
    theta, &
    dtime, &
    iconn, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    istrpgas, &
    utindep, &
    rsia, &
    ncomp, &
    ithcomp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build the transport equation for all nodes for SIA 
!
!   $Arguments:
!
 
type(t_cheproo), intent(in)                   :: this

integer, intent(in)                           :: nmobph ! Total number of mobile phases

integer, intent(in)                           :: nnch 

integer, intent(in)                           :: mxconn

integer, intent(in)                           :: nband

integer, intent(in)                           :: typesto

integer, intent(in)                           :: ithcomp

integer, intent(in)                           :: ncomp

real*8, intent(in)                            :: theta

real*8, intent(in)                            :: dtime

real*8, intent(in),dimension(nnch,2*nband+1)  :: efloww_ij

real*8, intent(in),dimension(nnch,2*nband+1)  :: eflowg1_ij

real*8, intent(in),dimension(nnch,2*nband+1)  :: eflowg2_ij

real*8, intent(in), dimension(nnch,ncomp)     :: rsia

type(t_vector), intent(in), dimension(nnch)   :: utindep

integer, intent(in), dimension(mxconn,nnch)   :: iconn

logical, intent(in)                           :: istrpgas

real*8, pointer, dimension(:)                 :: b

real*8, pointer, dimension(:,:)               :: trpmatrix

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
integer                               :: &
 i, &
 ithnch, &
 ic, &
 ncon, &
 nbandt
real*8                                :: &
 irsia, &
 iutindep
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
if (this%itypesolver/=siasolver.and.this%itypesolver/=sniasolver) then
 msg='Error, SIA (or SNIA) solver is not performed in CHEPROO'
 goto 10
end if 
!%------------------------------------------------------------
if (nnch/=this%numnch) then
 msg='Error in nodal chemistry objects'
 goto 10
end if 
!%------------------------------------------------------------
nbandt=3*nband+1
!%------------------------------------------------------------
!% Allocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (trpmatrix,nbandt,nnch,.true.)
call check_pointer_ (b,nnch,.true.)
!%------------------------------------------------------------
do ithnch=1,nnch
 
 irsia=rsia(ithnch,ithcomp)
 iutindep=utindep(ithnch)%vector(ithcomp)
 
 call build_trp_eq_sia_ithnch_ &
   (this, &
    trpmatrix, &
    b, &
    nnch, &
    nmobph, &
    mxconn, &
    nband, &
    theta, &
    dtime, &
    iconn(1:mxconn,ithnch), &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    istrpgas, &
    iutindep, &
    irsia, &
    ithnch, &
    msg, &
    iserror)
 
 if (iserror) goto 10
end do
!%-------------------------------------------------------------
return
 
10 continue 
print *,'****************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: build_trp_eq_sia_'
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
subroutine update_r_and_check_convergence_sia_cheproo &
   (this, &
    rsiak1, &
    utra, &
    nnch, &
    ncomp, &
    inchmxerrsia, &
    mxerrsia, &
    namemxerrsia, &
    dtime, &
    theta, &
    tolsia, &
	factor, &
    isconvsia, &
    isconvchem, &
    iter, &
    isopensystem, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute reaction term for SIA (Sequential Iteration Approach) 
! for all nodal chemistry objects defined in the CHEPROO and check convergence. 
! inchmxerrsia=ith nodal chemistry where the error is maximun 
! nammxerrsia=name of ith component where de error in rsia is maximum
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)                    :: this

integer, intent(in)                                :: ncomp ! Number of components 

integer, intent(in)                                :: nnch ! Number of nodal chemistry objects

integer, intent(in)                                :: iter ! Iteration 

real*8, pointer, dimension(:,:)                    :: rsiak1 ! Reaction term for SIA in k+1 

real*8, intent(out)                                :: mxerrsia    ! Maximum error in SIA

character(len=*), intent(out)                      :: namemxerrsia ! name of ith component where de error in rsia

logical, intent(out)                               :: iserror      ! If true there was error

logical, intent(out)                               :: isconvsia    ! If true, there was convergence in SIA

logical, intent(out)                               :: isconvchem   ! If true there was convergence in chemistry 

integer, intent(out)                               :: inchmxerrsia    ! O  ith nodal chemistry where the error is maximum

real*8, intent(in), dimension(nnch,ncomp), target  :: utra ! Transport solution 

real*8, intent(in)                                 :: dtime ! Time increment 

real*8, intent(in)                                 :: theta ! Temporal weight 

real*8, intent(in)                                 :: tolsia ! Tolerance in SIA

real*8, intent(in)                                 :: factor  ! Factor to scale solution 

logical, intent(in)                                :: isopensystem ! If true the system is considered open 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
character(len=100)                     :: &
 msg, &
 name
integer                                :: &
 ithnch
real*8                                 :: &
 mxerror
logical                                :: &
 isconvsiaith
type(t_nodalchemistry), pointer        :: &
 pnchk => null (), &
 pnchk1 => null ()
real*8, pointer                        :: &
 u(:) => null (), &
 rk1(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%--------------------------------------------------------------
iserror=.false.
msg=" "
!%--------------------------------------------------------------
!% Initialice variables 
!%--------------------------------------------------------------
isconvsia=.true.
mxerrsia=0.0d0
namemxerrsia=" "
!%--------------------------------------------------------------
!% Check the number of nodal chemistry objects
!%--------------------------------------------------------------
if (nnch/=this%numnch) then
 msg='Error in number of nodal chemistry objects'
 goto 10
end if
!%--------------------------------------------------------------
!% Allocate rsia in k+1
!%--------------------------------------------------------------
call check_pointer_ (rsiak1,nnch,ncomp,.true.)
!%-----------------------------------------------------------------
! Compute r sia for all nodes
!%-----------------------------------------------------------------
do ithnch=1,this%numnch
 
 pnchk => this%pnodalchemk(ithnch)%ptr
 pnchk1 => this%pnodalchemk1(ithnch)%ptr
 u => utra(ithnch,1:ncomp)
 rk1 => rsiak1(ithnch,1:ncomp)

 call update_and_check_ &
   (pnchk1, &
    pnchk, &
    rk1, &
    u, &
    ncomp, &
    mxerror, &
    name, &
    dtime, &
    theta, &
    tolsia, &
	factor, &
    isconvsiaith, &
    isconvchem, &
    iter, &
    isopensystem, &
    this%zero, &
    iserror)
 
 if (iserror) goto 20 

 
 if (.not.isconvchem) then
  isconvsia=.false.
  exit
 end if
 
 if (.not.isconvsiaith) isconvsia=.false.
 
 if (mxerror>mxerrsia) then
  inchmxerrsia=ithnch
  mxerrsia=mxerror
  namemxerrsia=name
 end if
 
end do
!%--------------------------------------------------------------
!% Nullify local pointers 
!%--------------------------------------------------------------
20 continue 
pnchk => null ()
pnchk1 => null ()
rk1 => null ()
u => null ()
if (iserror) goto 10 
!%--------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: update_and_check_'
print *, msg
print *,'********************************'
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine solve_react_trp_step_cheproo &
   (this, &
    iconn, &
    mxconn, &
    nnch, &
    typesto, &
    mxitertr, &
    mxdiverg, &
    mxchcompz, &
    theta, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    nband, &
    ngas, &
    dtime, &
    istrpgas, &
    isupdatejac, &
    caudalboundw, &
    caudalrechw, &
    caudalrechg, &
    caudalboundg, &
    itypeboundw, &
    itypeboundg, &
    gammaboundw, &
    isetboundw, &
    isetrechw, &
    isconvergence, &
	iter, &
    factor, &
    tolunk, &
    tolres, &
    isopensystem, &
    time, &
	iswcompbal, &
	mxerror, &
	issolvermrmt, &
	imethod, &
	mxnbox, &
	iconnbox, &
	alpha, &
	phi, &
	iserror, & 
	iouwinfo, &
    molinout, &
	molkin)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Solve one reactive transport step (for DSA or for SIA). 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)                                   :: this             ! Type cheproo variable 

integer, intent(in)                                               :: mxconn           ! Maximum number of connected nodal chemistries
 
integer, intent(in)                                               :: nnch             ! Number of nodal chemistries objects

integer, intent(in)                                               :: typesto          ! Type matrix storage of jacobian (only for DSA)

integer, intent(in)                                               :: mxitertr         ! Maximum number of transport iterarions

integer, intent(in)                                               :: mxdiverg         ! Maximum number of divergence

integer, intent(in)                                               :: mxchcompz        ! Maximum number of changes of components zones

integer, intent(in)                                               :: nband            ! Half-band

integer, intent(in)                                               :: ngas             ! Number of gases

real*8, intent(in)                                                :: theta            ! Temporal weight 

real*8, intent(in)                                                :: dtime            ! Time increment [s]

real*8, intent(in)                                                :: factor           ! Maximum factor for updating concentrations after a Newton-Raphson iteration (only for DSA) 

real*8, intent(in)                                                :: tolunk           ! Tolerance for unknowns 

real*8, intent(in)                                                :: tolres           ! Tolerance in residual (only for DSA)

real*8, intent(in)                                                :: time             ! Absolute time  

real*8, intent(in), dimension(nnch,2*nband+1)                     :: efloww_ij        ! Water flow field

real*8, intent(in), dimension(nnch,2*nband+1)                     :: eflowg1_ij       ! Gas flow field (1)

real*8, intent(in), dimension(nnch,2*nband+1)                     :: eflowg2_ij       ! Gas flow field (2) 

real*8, intent(in), dimension(nnch)                               :: caudalboundw     ! Water flow in the boundary 

real*8, intent(in), dimension(nnch)                               :: caudalrechw      ! Rechafge water flow

real*8, intent(in), dimension(nnch)                               :: caudalboundg     ! Gas flow in the boundary 

real*8, intent(in), dimension(nnch)                               :: caudalrechg      ! Recharge gas flow 

real*8, intent(in), dimension(nnch)                               :: gammaboundw      ! Coefficient for the mixed boundary condition (only necessary if itypeboundw = 4)

integer, intent(in), dimension(mxconn,nnch)                       :: iconn            ! Conectivity between nodal chemistries 

integer, intent(in), dimension(nnch)                              :: isetboundw       ! Number of the water corresponding to the boundary specified in CHEPROO input file. 

integer, intent(in), dimension(nnch)                              :: isetrechw        ! Number of the water corresponding to recharge specified in CHEPROO input file.  

integer, intent(in), dimension(nnch)                              :: itypeboundw      ! Type of boundary condition: 
                                                                                      ! 0: natural 
															                          ! 1: Fixed concentration
																                      ! 2: Fixed mass flux equal to liquid flux times concentration
																                      ! 3: Fixed mass flux
																                      ! 4: Mixed condition 

integer, intent(in), dimension(nnch)                              :: itypeboundg 
 
logical, intent(in)                                               :: istrpgas         ! If .true. then gas is transported 

logical, intent(in)                                               :: isupdatejac      ! If .true. then update jacobian in each iteration (only for DSA)

logical, intent(in)                                               :: isopensystem     ! If .true. the system is considered open 

logical, intent(out)                                              :: isconvergence    ! If .true. there was convergence in the reactive transport step 

integer, intent(out)                                              :: iter             ! Number of iterations for to solve the reactive transport step 

logical, intent(out)                                              :: iserror          ! If .true. there was error 

logical, intent(in)                                               :: iswcompbal       ! If true, then water components mass balance is evaluated 

real*8, intent(in)                                                :: mxerror          ! Maximum absolute error in residual (only for DSA)
 
logical, intent(in)                                               :: issolvermrmt      ! .true. if there is error

integer, intent(in)                                               :: imethod 

integer, intent(in)                                               :: mxnbox

real*8, intent(in), dimension(mxnbox,this%numnch)                 :: alpha           ! Alpha for boxes 

real*8, intent(in), dimension(mxnbox,this%numnch)                 :: phi             ! Porosity for boxes 

integer, intent(in), dimension(mxnbox+1,this%numnch)              :: iconnbox

integer, intent(in), optional                                     :: iouwinfo         ! Output unit for to write the Newton-Raphson infomation 

real*8, pointer, dimension(:), optional                           :: molinout         ! Mol of components that input/output in the system. 

real*8, pointer, dimension(:), optional                           :: molkin           ! Mol of components consumed/produced due kinetic reactions
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License: Barcelona, 2007 
!
!-------------------------------------------------------------------------
character(len=100)                   :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Check the time increment 
!%--------------------------------------------------------------
if (dtime==0.0d0) then
 msg='Error, time increment is equal to 0.0d0'
 goto 10
end if
!%--------------------------------------------------------------
select case (this%itypesolver)
case (siasolver,sniasolver)
 call solve_react_trp_step_sia_ &
   (this, &
    iconn, &
    mxconn, &
    nnch, &
    typesto, &
    mxitertr, &
    mxdiverg, &
    theta, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    nband, &
    ngas, &
    dtime, &
    istrpgas, &
    caudalboundw, &
    caudalrechw, &
    caudalrechg, &
    caudalboundg, &
    itypeboundw, &
    itypeboundg, &
    gammaboundw, &
    isetboundw, &
    isetrechw, &
    isconvergence, &
	iter, &
    factor, &
    tolunk, &
    isopensystem, &
    time, &
	iswcompbal, &
	mxerror, & 
    iserror, &
    molinout=molinout, &
	molkin=molkin, &
	iouwinfo=iouwinfo)
 
case (dsasolver)
 call solve_react_trp_step_dsa_ &
   (this, &
    iconn, &
    mxconn, &
    nnch, &
    typesto, &
    mxitertr, &
    mxdiverg, &
    mxchcompz, &
    theta, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    nband, &
    ngas, &
    dtime, &
    istrpgas, &
    isupdatejac, &
    caudalboundw, &
    caudalrechw, &
    caudalrechg, &
    caudalboundg, &
    itypeboundw, &
    itypeboundg, &
    gammaboundw, &
    isetboundw, &
    isetrechw, &
    isconvergence, &
	iter, &
    factor, &
    tolunk, &
    tolres, &
    isopensystem, &
    time, &
	mxerror, &
	issolvermrmt, &
	imethod, &
	mxnbox, &
	iconnbox, &
	alpha, &
	phi, &
	iserror, &
	molinout=molinout, &
	molkin=molkin, &
	iouwinfo=iouwinfo)
end select
!%--------------------------------------------------------------
return
 
10 continue 
print *,'******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: solve_react_trp_step_'
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
subroutine build_indep_term_dsa_cheproo &
   (this, &
    tindep, &
    nnch, &
    nsp, &
    theta, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    nband, &
    ngas, &
    dtime, &
    istrpgas, &
    caudalboundw, &
    caudalrechw, &
    caudalrechg, &
    caudalboundg, &
    itypeboundw, &
    itypeboundg, &
    isetboundw, &
    isetrechw, &
    gammaboundw, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build the independent term for DSA 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                    :: this

real*8, pointer, dimension(:,:)                 :: tindep        ! Independent term 

integer, intent(in)                             :: nnch          ! Number of nodal chemistry objects 

integer, intent(in)                             :: nband

integer, intent(in)                             :: ngas

integer, intent(out)                            :: nsp           ! Number of species 

real*8, intent(in)                              :: theta

real*8, intent(in)                              :: dtime         ! Time increment 

real*8, intent(in), dimension(nnch,2*nband+1)   :: efloww_ij

real*8, intent(in), dimension(nnch,2*nband+1)   :: eflowg1_ij

real*8, intent(in), dimension(nnch,2*nband+1)   :: eflowg2_ij

real*8, intent(in), dimension(nnch)             :: caudalboundw

real*8, intent(in), dimension(nnch)             :: caudalrechw

real*8, intent(in), dimension(nnch)             :: caudalboundg

real*8, intent(in), dimension(nnch)             :: caudalrechg

integer, intent(in), dimension(nnch)            :: isetboundw

integer, intent(in), dimension(nnch)            :: isetrechw

integer, intent(in), dimension(nnch)            :: itypeboundw

integer, intent(in), dimension(nnch)            :: itypeboundg

real*8, intent(in), dimension(nnch)             :: gammaboundw

logical, intent(in)                             :: istrpgas

logical, intent(out)                            :: iserror 
 
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
real*8                               :: &
 alpha, &
 fii, &
 eij, &
 caudal, &
 gamma 
integer                              :: &
 iband, &
 ithnch, &
 jthnch, &
 nmobph, &
 ib, &
 ir, &
 tbound, &
 nbands, &
 nbandt
real*8, pointer                      :: &
 c(:) => null (), &
 cmob(:,:) => null (), &
 sktrk(:) => null (), &
 ptindep(:) => null ()
type(t_nodalchemistry), pointer      :: &
 pnchk => null (), &
 pnchk1 => null (), &
 pwater => null () 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Check the solver type
!%--------------------------------------------------------------
if (this%itypesolver/=dsasolver) then
  msg='Error, DSA solver is not performed in CHEPROO'
  goto 10
end if
!%--------------------------------------------------------------
!% Check the number of nodal chemistry objects 
!%--------------------------------------------------------------
if (nnch/=this%numnch) then
  msg='Error in number of nodal chemistry objects'
  goto 10
end if
!%--------------------------------------------------------------
!% Check the number of chemical system objects 
!%--------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system'
 goto 10
end if
!%--------------------------------------------------------------
!% Check the time increment 
!%--------------------------------------------------------------
if (dtime==0.0d0) then
  msg='Error, time increment is equal to 0.0d0'
  goto 10
end if
!%--------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numsp=nsp)
call check_pointer_ (tindep,nsp,nnch,.true.)
!%--------------------------------------------------------------
nbands=nband+1
nbandt=2*nband+1
!%--------------------------------------------------------------
do ithnch=1,nnch
    
	pnchk => this%pnodalchemk(ithnch)%ptr
	pnchk1 => this%pnodalchemk1(ithnch)%ptr
!%--------------------------------------------------------------
    ptindep => tindep(1:nsp,ithnch)	
!%--------------------------------------------------------------
!% Get information to nodal chemistry in k + 1
!% Get the mass of free water   
!%--------------------------------------------------------------
    call get_chem_info_ (pnchk1,iserror,omgwfree=fii)
	if (iserror) then
      msg='Error when calling get_chem_info_'
      call add_ (msg,jthnch)
      goto 20
    end if
!%--------------------------------------------------------------
	ib=isetboundw(ithnch)
    ir=isetrechw(ithnch)
    tbound=itypeboundw(ithnch)
    caudal=caudalboundw(ithnch)
    gamma=gammaboundw(ithnch)
!%--------------------------------------------------------------
! 1) Add boundary term
!%--------------------------------------------------------------
    if ((ib>0.and.tbound==1).or. &
	    (ib>0.and.tbound==4).or. &
        (ib>0.and.tbound==2.and.caudal>0.0d0)) then
	     pwater => this%pwater(ib)%ptr	 
         call get_chem_info_(pwater,iserror,cmob=cmob,nummobph=nmobph)
         if (iserror) then
          msg='Error when calling get_chem_info_'
          call add_ (msg,ib)
          goto 20
         end if
     if (tbound==1) then
       alpha=1.0d0
      else if(tbound==2) then
       alpha=caudal
      else if(tbound==4) then 
       alpha=gamma
      end if
      call add_(ptindep,alpha,cmob(1:nsp,1),nsp)
    end if
!%--------------------------------------------------------------
! 2) Add recharge term
!%--------------------------------------------------------------
    if (ir>0.and.caudalrechw(ithnch)>0.0d0) then
         pwater => this%pwater(ir)%ptr
	     call get_chem_info_(pwater,iserror,cmob=cmob,nummobph=nmobph)
         if (iserror) then
          msg='Error when calling get_chem_info_'
          call add_ (msg,ir)
          goto 20
         end if
         alpha=caudalrechw(ithnch)
         call add_(ptindep,alpha,cmob(1:nsp,1),nsp)
    end if
!%--------------------------------------------------------------
! 3) Add storage term
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Get information to nodal chemistry in k 
!%--------------------------------------------------------------
    call get_chem_info_ (pnchk,iserror,c=c,sktrk=sktrk)
    if (iserror) then
      msg='Error when calling get_chem_info_'
      call add_ (msg,ithnch)
      goto 20
    end if
!%--------------------------------------------------------------    
    alpha=fii/dtime
    call add_ (ptindep,alpha,c,nsp)
!%cprovi--------------------------------------------------------
!%cprovi--------------------------------------------------------------
!%cprovi--------------------------------------------------------------
!% call get_chem_info_(pnchk,iserror,caq=c)
!% alpha=fii/dtime
!%  call add_(ptindep,alpha,c,nsp)
!%  call get_chem_info_ (pnchk,iserror,omgwfree=fii,cmineq=c)    
!%  alpha=fii/dtime
!%  call add_ (ptindep,alpha,c,nsp)
!%  call get_chem_info_ (pnchk,iserror,cads=c)
!%  call add_ (ptindep,alpha,c,nsp)
!%cprovi--------------------------------------------------------------
!%cprovi--------------------------------------------------------------
!%cprovi--------------------------------------------------------------
    if (theta/=1.0d0) then
!%--------------------------------------------------------------
! 4) Add kinetic term (theta-1)*Sktrk
!%--------------------------------------------------------------
    alpha=(theta-1.0d0)*fii
    call add_ (ptindep,alpha,sktrk,nsp)
!%--------------------------------------------------------------
! 4) Add contribution of k
!%--------------------------------------------------------------
   do iband=1,nbandt
      jthnch=iband+ithnch-nbands
      eij=efloww_ij(ithnch,iband)
      if (jthnch>=1.and.jthnch<=nnch) then
        pnchk => this%pnodalchemk(jthnch)%ptr 
	    call get_chem_info_ (pnchk,iserror,cmob=cmob,nummobph=nmobph)
        if (iserror) then
          msg='Error when calling get_chem_info_'
          call add_ (msg,jthnch)
          goto 20
        end if
        alpha=(theta-1.0d0)*eij
        call add_ (ptindep,alpha,cmob(1:nsp,1),nsp)
      end if
   end do
end if
!%--------------------------------------------------------------
 
end do
!%--------------------------------------------------------------
20 continue 
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (sktrk,1,.false.)
call check_pointer_ (cmob,1,1,.false.)
call check_pointer_ (c,1,.false.)
!%--------------------------------------------------------------
!% Nullify local pointers 
!%--------------------------------------------------------------
pnchk => null ()
pnchk1 => null ()
pwater => null ()
ptindep => null ()
if (iserror) goto 10
!%--------------------------------------------------------------
return
 
10 continue 
print *,'************************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: build_indep_term_dsa_'
print *, msg
print *,'************************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine build_indep_term_sia_cheproo &
   (this, &
    tindep, &
    nnch, &
    nsp, &
    theta, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    nband, &
    ngas, &
    dtime, &
    istrpgas, &
    caudalboundw, &
    caudalrechw, &
    caudalrechg, &
    caudalboundg, &
	gammaboundw, &
    itypeboundw, &
    itypeboundg, &
    isetboundw, &
    isetrechw, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build the independent term for SIA (or SNIA). 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                    :: this

real*8, pointer, dimension(:,:)                 :: tindep          ! Independent term for SIA [nsp x nnch]

integer, intent(out)                            :: nsp             ! Number of species

logical, intent(out)                            :: iserror         ! iserror=true, then there was an error

integer, intent(in)                             :: nnch            ! Number of nodal chemistry objects 

integer, intent(in)                             :: nband            

integer, intent(in)                             :: ngas            ! Number of gas phases 

real*8, intent(in)                              :: theta           ! Tempral weight 

real*8, intent(in)                              :: dtime           ! Time increment 

real*8, intent(in), dimension(nnch,2*nband+1)   :: efloww_ij

real*8, intent(in), dimension(nnch,2*nband+1)   :: eflowg1_ij

real*8, intent(in), dimension(nnch,2*nband+1)   :: eflowg2_ij

real*8, intent(in), dimension(nnch)             :: caudalboundw

real*8, intent(in), dimension(nnch)             :: caudalrechw

real*8, intent(in), dimension(nnch)             :: caudalboundg

real*8, intent(in), dimension(nnch)             :: caudalrechg

real*8, intent(in), dimension(nnch)             :: gammaboundw

integer, intent(in), dimension(nnch)            :: isetboundw

integer, intent(in), dimension(nnch)            :: isetrechw

integer, intent(in), dimension(nnch)            :: itypeboundw

integer, intent(in), dimension(nnch)            :: itypeboundg

logical, intent(in)                             :: istrpgas
 
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
real*8                               :: &
 alpha, &
 fii, &
 eij, &
 gamma, &
 caudal
integer                              :: &
 iband, &
 ithnch, &
 jthnch, &
 nmobph, &
 ib, &
 ir, &
 tbound, &
 nbands, &
 nbandt
real*8, pointer                      :: &
 c(:) => null (), &
 cmob(:,:) => null (), &
 ptindep(:) => null ()
type(t_nodalchemistry), pointer      :: &
 pnchk => null (), &
 pnchk1 => null (), &
 pwater => null ()
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------
if (this%itypesolver/=siasolver.and.this%itypesolver/=sniasolver) then
  msg='Error, SIA (or SNIA) solver is not performed in CHEPROO'
  goto 10
end if
!%--------------------------------------------------------------
!% Check the number of nodal chemistry objects 
!%--------------------------------------------------------------
if (nnch/=this%numnch) then
  msg='Error in number of nodal chemistry objects'
  goto 10
end if
!%--------------------------------------------------------------
!% Check the number of chemical system objects. 
!%--------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system'
 goto 10
end if
!%--------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numsp=nsp)
call check_pointer_ (tindep,nsp,nnch,.true.)
!%--------------------------------------------------------------
nbands=nband+1
nbandt=2*nband+1
!%--------------------------------------------------------------
! Start loop by nodal chemistry objects 
!%--------------------------------------------------------------
do ithnch=1,nnch
    
	pnchk => this%pnodalchemk(ithnch)%ptr
	pnchk1 => this%pnodalchemk1(ithnch)%ptr
	
	ptindep => tindep(1:nsp,ithnch)
	
    call get_chem_info_ (pnchk1,iserror,omgwfree=fii)
    if (iserror) then
       msg='Error when calling get_chem_info_ in nodal chemistry:'
       call add_ (msg,ithnch)
       goto 20
    end if
!%--------------------------------------------------------------    
    ib=isetboundw(ithnch)
    ir=isetrechw(ithnch)
    tbound=itypeboundw(ithnch)
	caudal=caudalboundw(ithnch)
	gamma=gammaboundw(ithnch)
!%--------------------------------------------------------------
! Add boundary term
! theta* g 
!%--------------------------------------------------------------
    if ((ib>0.and.tbound==1).or. &
	    (ib>0.and.tbound==4).or. &
        (ib>0.and.tbound==2.and.caudal>0.0d0)) then
        pwater => this%pwater(ib)%ptr
		call get_chem_info_(pwater,iserror,cmob=cmob,nummobph=nmobph)
        if (iserror) then
         msg='Error when calling get_chem_info_ in water:'
         call add_ (msg,ib)
         goto 20
        end if
        if (tbound==1) then
         alpha=1.0d0
        else if(tbound==2) then
         alpha=caudal
        else if(tbound==4) then
         alpha=gamma
        end if
        call add_(ptindep,alpha,cmob(1:nsp,1),nsp)
    end if
!%--------------------------------------------------------------
! Add recharge term
!%--------------------------------------------------------------
      if (ir>0.and.caudalrechw(ithnch)>0.0d0) then
         pwater => this%pwater(ir)%ptr
		 call get_chem_info_(pwater,iserror,cmob=cmob,nummobph=nmobph)
         if (iserror) then
         msg='Error when calling get_chem_info_ in water:'
         call add_ (msg,ir)
         goto 20
         end if
        alpha=caudalrechw(ithnch)
        call add_(ptindep,alpha,cmob(:,1),nsp)
      end if
!%--------------------------------------------------------------
! Add contribution in k
!%--------------------------------------------------------------
    call get_chem_info_(pnchk,iserror,cmob=cmob,nummobph=nmobph)
    if (iserror) then
     msg='Error when calling get_chem_info_ in nodal chemistry:'
     call add_ (msg,ithnch)
     goto 20
    end if
!%--------------------------------------------------------------
! Add storage term
! fii
! --- cmob 
! dt 
!%--------------------------------------------------------------
 alpha=fii/dtime
 call add_ (ptindep,alpha,cmob(1:nsp,1),nsp)
!%--------------------------------------------------------------
if (theta/=1.0d0) then
!%--------------------------------------------------------------
! Add contribution in k
!%--------------------------------------------------------------
   do iband=1,nbandt
          jthnch=iband+ithnch-nbands
          eij=efloww_ij(ithnch,iband)
          if (jthnch>=1.and.jthnch<=nnch) then
               pnchk => this%pnodalchemk(jthnch)%ptr
			   call get_chem_info_(pnchk,iserror,cmob=cmob,nummobph=nmobph)
               if (iserror) then
                 msg='Error when calling get_chem_info_ in nodal chemistry:'
                 call add_ (msg,ithnch)
                 goto 20
               end if
               alpha=(theta-1.0d0)*eij
               call add_ (ptindep,alpha,cmob(1:nsp,1),nsp)
          end if
   end do
end if
!%--------------------------------------------------------------
end do
!%--------------------------------------------------------------
! Finish loop by nodal chemistry objects 
!%--------------------------------------------------------------
20 continue 
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (cmob,1,1,.false.)
call check_pointer_ (c,1,.false.)
!%--------------------------------------------------------------
!% Nullify local pointers 
!%--------------------------------------------------------------
pnchk => null ()
pnchk1 => null ()
pwater => null ()
ptindep => null ()
if (iserror) goto 10
!%--------------------------------------------------------------
return
 
10 continue 
print *,'******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: build_indep_term_sia_'
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
subroutine build_jacob_resid_dsa_cheproo &
   (this, &
    jacobian, &
    residual, &
    nunk, &
    nbandjac, &
    iconn, &
	itypeboundw, &
    nnch, &
    mxconn, &
    nband, &
    istrpgas, &
    typesto, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    ngas, &
    theta, &
    tindep, &
    nsp, &
    dtime, &
	issolvermrmt, &
	imethod, &
	mxnbox, &
	iconnbox, &
	alpha, &
	phi, &
	isconvrmrmt, &
	iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build the jacobian and residual for DSA (Direct Substitution Approach)
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                                     :: this

integer, intent(in)                                              :: nnch

integer, intent(out)                                             :: nunk

integer, intent(in)                                              :: mxconn

integer, intent(in)                                              :: nband

integer, intent(in)                                              :: typesto

integer, intent(in)                                              :: nsp

integer, intent(in)                                              :: ngas

integer, intent(out)                                             :: nbandjac

real*8, intent(in)                                               :: theta                     ! Temporal weight 

real*8, intent(in)                                               :: dtime                     ! Time increment

integer, intent(in)                                              :: itypeboundw(nnch)

real*8, intent(in)                                               :: efloww_ij(nnch,2*nband+1)

real*8, intent(in)                                               :: eflowg1_ij(nnch,2*nband+1)

real*8, intent(in)                                               :: eflowg2_ij(nnch,2*nband+1)

real*8, intent(in)                                               :: tindep(nsp,nnch)           ! Independent term for D
 
real*8, pointer                                                  :: jacobian(:,:)

real*8, pointer                                                  :: residual(:)

integer, intent(in)                                              :: iconn(mxconn,nnch)                         ! Nodes connectivity

logical, intent(in)                                              :: istrpgas          ! .true. the gas phases is transported

logical, intent(out)                                             :: iserror         ! .true. if there is error 

logical, intent(in)                                              :: issolvermrmt      ! .true. reactive multi-rate mass transfer is solved

logical, intent(out)                                             :: isconvrmrmt     

integer, intent(in)                                              :: imethod 

integer, intent(in)                                              :: mxnbox

real*8, intent(in), dimension(mxnbox,this%numnch)                :: alpha           ! Alpha for boxes 

real*8, intent(in), dimension(mxnbox,this%numnch)                :: phi             ! Porosity for boxes 

integer, intent(in), dimension(mxnbox,this%numnch)               :: iconnbox 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 ithnch, &
 jthnch, &
 icz, &
 mxnunkz,      & ! Maximun nodal unknowns
 indx,         & ! Index of components set
 icn, &
 nsploc, &
 ipos1, &
 ipos2, &
 ndim, &
 nczjthnch, &
 nmobph, &
 i, &
 j, &
 nbox
type(t_array), pointer                 :: &
 dumob(:,:) => null (), &
 dusktrk(:) => null (), &
 du(:) => null (), &
 duads(:) => null ()
type(t_vector), pointer                :: &
 umob(:,:) => null (), &
 usktrk(:) => null (), &
 u(:) => null (), &
 uads(:) => null (), &
 utindep(:) => null ()
integer, pointer                       :: &
 indexconn(:,:) => null (), &
 czlocconn(:,:) => null (), &
 posunk(:,:) => null (), &
 nunkconn(:,:) => null (), &
 ncznod(:) => null ()
character(len=100)                     :: &
 msg
logical                                :: &
 isconvergence         
integer, parameter                     :: &
 bandstored=0, &
 sparsstored=1 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
isconvrmrmt=.true. 
!%--------------------------------------------------------------
!% Check if the reactive transport solver is DSA
!%--------------------------------------------------------------
if (this%itypesolver/=dsasolver) then
  msg='Error, DSA solver is not performed in CHEPROO'
  goto 10
end if
!%---------------------------------------------------------------
!% Check if there are chemical system objects in the CHEPROO 
!% object 
!%---------------------------------------------------------------
if (this%numchemsys==0) then
  msg='Error, not defined chemical system'
  goto 10
end if
!%---------------------------------------------------------------
!% Check the number of species
!%---------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numsp=nsploc)
if (nsploc/=nsp) then
  msg='Error in number of species'
  goto 10
end if
!%----------------------------------------------------------------
!% Check the number of nodal chemistry objects 
!%----------------------------------------------------------------
if (nnch/=this%numnch) then
  msg='Error in number of nodal chemistry objects'
  goto 10
end if
!%----------------------------------------------------------------
!% Allocate local pointers 
!%----------------------------------------------------------------
call check_pointer_ (indexconn,mxconn,nnch,.true.)
call check_pointer_ (czlocconn,mxconn,nnch,.true.)
call check_pointer_ (nunkconn,mxconn,nnch,.true.)
call check_pointer_ (ncznod,nnch,.true.)
call check_pointer_ (posunk,2,nnch,.true.)
!%---------------------------------------------------------------
!% Make the component sets associated each nodal chemistry
!%---------------------------------------------------------------
call make_preview_info_dsa_ &
   (this, &
    czlocconn, &
    indexconn, &
    nunkconn, &
    ncznod, &
    posunk, &
    nunk, &
    mxnunkz, &
    nnch, &
    iconn, &
    mxconn, &
    msg, &
    iserror)
 
 if (iserror) goto 20
!%-----------------------------------------------------------------
! Preview tasks, Allocate and initialice jacobian and residual
!%-----------------------------------------------------------------
!%----------------------------------------------------------------
nbandjac=(nband+1)*mxnunkz-1
!%----------------------------------------------------------------
!% Allocate and zeroing jacobian and residual 
!%----------------------------------------------------------------
call check_pointer_ (jacobian,3*nbandjac+1,nunk,.true.)
call check_pointer_ (residual,nunk,.true.) 
!%----------------------------------------------------------------
!% End preview tasks
!% Build the jacobian and residual DSA
!%----------------------------------------------------------------
do jthnch=1,nnch
 
 nczjthnch=ncznod(jthnch)
 
 call make_zonal_chem_info_ithnch_ &
   (this, &
    nczjthnch, &
    indexconn(1:nczjthnch,jthnch), &
    nunkconn(1:nczjthnch,jthnch), &
    tindep(1:nsp,jthnch), &
    utindep, &
    u, &
    umob, &
    uads, &
    usktrk, &
    du, &
    dumob, &
    duads, &
    dusktrk, &
    nmobph, &
    jthnch, &
    nsp, &
    .true., &
    msg, &
    iserror)
 
 if (iserror) goto 20
 
 
 call build_jacob_resid_dsa_ithnch_ &
   (this, &
    jacobian, &
    residual, &
    nnch, &
    nmobph, &
    nczjthnch, &
    nbandjac, &
    mxconn, &
    nband, &
    nunk, &
    ngas, &
    theta, &
    dtime, &
	itypeboundw, &
    iconn(1:mxconn,jthnch), &
	efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    posunk, &
    czlocconn(1:mxconn,jthnch), &
    istrpgas, &
    utindep, &
    u, &
    umob, &
    uads, &
    usktrk, &
    du, &
    dumob, &
    duads, &
    dusktrk, &
    jthnch, &
    msg, &
    iserror)
 
 if (iserror) goto 20
 
 
!%------------------------------------------------------------
!% Solve reactive multi-rate mass transfer 
!%------------------------------------------------------------
if (issolvermrmt) then
 
 nbox=iconnbox(1,jthnch)
 call add_f_df_dsa_ith_ &
   (this, &
    jacobian, &
    residual, &
    nunk, &
    nbandjac, &
	nsp, &
    typesto, &
	nbox, &
    theta, &
    dtime, &
	posunk, &
	imethod, &
	iconnbox(2:nbox+1,jthnch), &
	alpha(1:nbox,jthnch), &
	phi(1:nbox,jthnch), &
	1.0d-4, &        ! Tolerance in unknowns 
	1.0d-10, &       ! Tolerance in residual 
	100, &           ! Maximum number of iterations 
	10, &            ! Maximum number of changes in component zones 
	isconvrmrmt, &
	jthnch, &
    iserror)

 if (iserror.or..not.isconvrmrmt) goto 20

end if 



end do
!%------------------------------------------------------------
!% Change the sign of the residual 
!%------------------------------------------------------------
residual=-residual
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
do i=1,nczjthnch
  deallocate(u(i)%vector)
  deallocate(uads(i)%vector)
  deallocate(utindep(i)%vector)
  deallocate(usktrk(i)%vector)
  deallocate(dusktrk(i)%array)
  deallocate(du(i)%array)
  deallocate(duads(i)%array)
 do j=1,nmobph
  deallocate(umob(i,j)%vector)
  deallocate(dumob(i,j)%array)
 end do
end do
deallocate (u)
deallocate (utindep)
deallocate (du)
deallocate (usktrk)
deallocate (umob)
deallocate (uads)
deallocate (dumob)
deallocate (duads)
deallocate (dusktrk)
call check_pointer_ (ncznod,1,.false.)
call check_pointer_ (czlocconn,1,1,.false.)
call check_pointer_ (nunkconn,1,1,.false.)
call check_pointer_ (posunk,1,1,.false.)
call check_pointer_ (indexconn,1,1,.false.)
if (iserror) goto 10
!%------------------------------------------------------------
!%------------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: build_jacob_resid_dsa_'
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
subroutine compute_total_mol_cheproo &
   (this, &
    molaq, &
    molprec, &
    molads, &
    ncomp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the total mol of components in all nodal chemistries defined in CHEPROO
! molaq[npri] moles of components in solution
! molprec[npri] moles of components precipitated
! molads[npri] moles of components adsorbed
!
!   $Arguments:
!
 
type(t_cheproo), intent(in)           :: this       ! Type CHEPROO variable 

integer, intent(out)                  :: ncomp      ! Number of components 

real*8, pointer, dimension(:)         :: molaq      

real*8, pointer, dimension(:)         :: molprec

real*8, pointer, dimension(:)         :: molads

logical, intent(out)                  :: iserror    ! iserror=true, then there was an error 
 
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
 inch
real*8, pointer                       :: &
 vector1(:) => null (), &
 vector2(:) => null (), &
 vector3(:) => null ()
character(len=100)                    :: &
 msg 
type(t_nodalchemistry), pointer       :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Check if there are chemical system objects inside CHEPROO 
!% object 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys(1)%ptr)) then
 msg='Error, not associated chemical system in CHEPROO'
 goto 10
end if
!%------------------------------------------------------------
call get_chem_info_(this%pchemsys(1)%ptr,iserror,numbase=ncomp)
!%------------------------------------------------------------
!% Allocate and initiallize pointers 
!%------------------------------------------------------------
call check_pointer_ (molaq,ncomp,.true.)
call check_pointer_ (molprec,ncomp,.true.)
call check_pointer_ (molads,ncomp,.true.)
!%------------------------------------------------------------
!% Loop for nodal chemistry objects 
!%------------------------------------------------------------
do inch=1,this%numnch
 pnch => this%pnodalchemk1(inch)%ptr
 call get_chem_info_ (pnch,iserror,molaq=vector1,molads=vector2,molprec=vector3)
 if (iserror) goto 20
 molaq=molaq+vector1
 molads=molads+vector2
 molprec=molprec+vector3
end do
!%------------------------------------------------------------
!% Loop for reactive multi-rate mass transfer objects 
!%------------------------------------------------------------
do inch=1,this%numrmrmtbox
 pnch => this%prmrmtboxk1(inch)%ptr
 call get_chem_info_ (pnch,iserror,molaq=vector1,molads=vector2,molprec=vector3)
 if (iserror) goto 20
 molaq=molaq+vector1
 molads=molads+vector2
 molprec=molprec+vector3
end do
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
! Deallocate local pointers
!%-------------------------------------------------------------
call check_pointer_ (vector1,1,.false.)
call check_pointer_ (vector2,1,.false.)
call check_pointer_ (vector3,1,.false.)
pnch => null ()
if (iserror) goto 10
!%-------------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: compute_total_mol_'
print *, msg
print *,'*******************************'
return
 
 end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_mol_balance_cheproo &
   (this, &
    molinout, &
	molkin, &
    molaq_ini, &
    molprec_ini, &
    molads_ini, &
    molaq_fin, &
    molprec_fin, &
    molads_fin, &
    ncomp, &
    ioutput, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write mol balance information for reactive transport 
!
!   $Arguments:
!

type (t_cheproo), intent(in)          :: this

integer, intent(in)                   :: ncomp 

integer, intent(in)                   :: ioutput

real*8, intent(in), dimension(ncomp)  :: molinout

real*8, intent(in), dimension(ncomp)  :: molkin

real*8, intent(in), dimension(ncomp)  :: molaq_ini

real*8, intent(in), dimension(ncomp)  :: molprec_ini

real*8, intent(in), dimension(ncomp)  :: molads_ini
 
real*8, intent(in), dimension(ncomp)  :: molaq_fin

real*8, intent(in), dimension(ncomp)  :: molprec_fin

real*8, intent(in), dimension(ncomp)  :: molads_fin

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
character(len=100)                   :: &
 msg
integer                              :: &
 ncomp1, &
 i
character(len=100), pointer          :: &
 namecomp(:) => null ()
real*8                               :: &
 increment, &
 inout, &
 kinchange, &
 error, &
 total_ini, &
 total_fin   
real*8, parameter                    :: &
 zero=1.0d-4
!-------------------------------------------------------------------------
!
!   $code
!
 
!%--------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system'
 goto 10
end if
!%--------------------------------------------------------------
!% Get chemical information 
!%--------------------------------------------------------------
call get_chem_info_ &
 (this%pchemsys(1)%ptr, &
  iserror,numbase=ncomp1, &
  namebase=namecomp)
if(iserror) goto 20   
!%--------------------------------------------------------------
!% Check the number of components 
!%--------------------------------------------------------------
if (ncomp1/=ncomp) then
 msg='Error in number of components'
 goto 20
end if
!%--------------------------------------------------------------
!% Open unit 
!%--------------------------------------------------------------
open(unit=ioutput,file='mol_balance.dat',status='unknown')
!%--------------------------------------------------------------
!% Write the signature in output unit 
!%--------------------------------------------------------------
call write_signature_ (this,ioutput)
!%--------------------------------------------------------------
write (ioutput,*) '--------------------------------------------'// &
                  '-------------------'
write (ioutput,*) 'CHEPROO name:',this%name
write (ioutput,*) '--------------------------------------------'// &
                  '-------------------'
!%--------------------------------------------------------------
write (ioutput,*) '--------------------------------------------'// &
                  '-------------------'
write (ioutput,*) '           Mol balance  of components              '
write (ioutput,*) '--------------------------------------------'// &
                  '-------------------'
do i=1,ncomp 
 
 increment=molaq_fin(i)+molprec_fin(i)+molads_fin(i)
 increment=increment-molaq_ini(i)-molprec_ini(i)-molads_ini(i)
 inout=molinout(i)
 kinchange=molkin(i)
!%--------------------------------------------------------------
!% Compute absolute error in components 
!%--------------------------------------------------------------
 error=dabs(inout-increment)
!%--------------------------------------------------------------
 write (ioutput,*) '--------------------------------------------'// &
                   '-------------------'
 write (ioutput,*) '--------------------------------------------'// &
                   '-------------------'
 write (ioutput,1) 'Component:',namecomp(i),'[mol]'
 write (ioutput,*) '--------------------------------------------'// &
                   '-------------------'
 write (ioutput,2) 'Initial','Final','Difference'
 write (ioutput,*) '--------------------------------------------'// &
                    '-------------------'
 write (ioutput,3) 'Disolved:',molaq_ini(i),molaq_fin(i),molaq_fin(i)-molaq_ini(i)
 write (ioutput,3) 'Precipitated:',molprec_ini(i),molprec_fin(i),molprec_fin(i)-molprec_ini(i)
 write (ioutput,3) 'Adsorbed:',molads_ini(i),molads_fin(i),molads_fin(i)-molads_ini(i)
 total_ini=molaq_ini(i)+molprec_ini(i)+molads_ini(i)
 total_fin=molaq_fin(i)+molprec_fin(i)+molads_fin(i)
 write (ioutput,3) 'Total:',total_ini,total_fin,total_fin-total_ini
 write (ioutput,*) '--------------------------------------------'// &
                   '-------------------'
 write (ioutput,4) 'Total increment:', increment 
 write (ioutput,4) 'Input-Output:', inout
 write (ioutput,4) 'Kinetic Changes:', kinchange
 write (ioutput,4) 'Absolute Error:', error
 error=error/dabs(inout)
 write (ioutput,4) 'Relative Error:', error 
 write (ioutput,*) '--------------------------------------------'// &
                   '-------------------'
 write (ioutput,*) '--------------------------------------------'// &
                   '-------------------'
end do
write (ioutput,*) '--------------------------------------------'// &
                  '-------------------'
!%--------------------------------------------------------------
20 continue
!%--------------------------------------------------------------
!% Close unit 
!%--------------------------------------------------------------
close(ioutput) 
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (namecomp,1,.false.)
if (iserror) goto 10
!%--------------------------------------------------------------
return
 
1 format (2a10,a5)
2 format (15x,3a15)
3 format (a15,3e15.7)
4 format (a16,e15.7)
 
 
10 continue 
print *,'******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: write_mol_balance_'
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
subroutine update_and_check_convergence_dsa_cheproo &
  (this,                 & 
   isconvergence,          & 
   isdivergence, &
   isanomalous, &
   isupmxitergam, &
   mxerrorunk,           & 
   mxerrorres,           & 
   namemxunk, &
   namemxres, &
   inchmxunk, &
   inchmxres, &
   deltacunk,            & 
   nunk,                 & 
   residual,             & 
   tolunk,               & 
   tolres,               & 
   dtime,                & 
   factor,               & 
   iserror)                
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Check the convergence computing the relative error and
!  the residual error. After update all nodal chemistries from
!  cpri (for DSA)
!
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)                 :: this            ! Type CHEPROO variable.  

integer, intent(in)                             :: nunk            ! Number of unknowns 

real*8, intent(in), dimension(nunk), target     :: deltacunk       ! Delta solution 

real*8, intent(in), dimension(nunk), target     :: residual        ! Residual 

real*8, intent(in)                              :: tolres          ! Tolerance in residual

real*8, intent(in)                              :: tolunk          ! Tolerance in unknown 

real*8, intent(in)                              :: dtime           ! Time increment 

real*8, intent(in)                              :: factor

logical, intent(out)                            :: isconvergence   

logical, intent(out)                            :: isdivergence

logical, intent(out)                            :: iserror

logical, intent(out)                            :: isanomalous

logical, intent(out)                            :: isupmxitergam

integer, intent(out)                            :: inchmxunk

integer, intent(out)                            :: inchmxres

real*8, intent(out)                             :: mxerrorres

real*8, intent(out)                             :: mxerrorunk

character(len=*), intent(out)                   :: namemxunk

character(len=*), intent(out)                   :: namemxres 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)            :: &
 msg
integer                       :: &
 ithnch, &
 ipos1, &
 ipos2, &
 npri
logical                       :: &
 isconv, &
 isupmxiter
real*8                        :: &
 mxerrorunkold, &
 mxunk, &
 mxres
real*8, pointer               :: &
 delta(:) => null (), &
 res(:) => null ()
character(len=100)            :: &
 name1, &
 name2 
type(t_nodalchemistry), pointer :: &
 pnchk1 => null (), &
 pnchk => null ()
real*8, parameter             :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
!%---------------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------
if (this%itypesolver/=dsasolver) then
  msg='Error, DSA solver is not performed in CHEPROO'
  goto 10
end if
!%---------------------------------------------------------------------
!% Initialice variables. 
!%---------------------------------------------------------------------
isconvergence=.true.
isdivergence=.false.
isanomalous=.false.
isupmxitergam=.false.
mxerrorunkold=mxerrorunk
mxerrorunk=r0
mxerrorres=r0
ipos1=1
ipos2=0
!%---------------------------------------------------------------------
do ithnch=1,this%numnch
 
  pnchk1 => this%pnodalchemk1(ithnch)%ptr 
  pnchk => this%pnodalchemk(ithnch)%ptr 
  
  call get_chem_info_ (pnchk1,iserror,npri=npri)
  if (iserror) goto 20
 
  ipos2=ipos2+npri
  
  delta => deltacunk(ipos1:ipos2)
  res => residual(ipos1:ipos2)
 
  call update_and_check_ &
  (pnchk1, &
   pnchk, &
   isconv, &
   isanomalous, &
   isupmxiter, &
   mxunk, &
   mxres, &
   name1, &
   name2, &
   delta, &
   res, &
   npri, &
   tolunk, &
   tolres, &
   factor, &
   dtime, &
   .true., &
   this%zero, &
   50.0d0, &
   iserror)
 
 if (iserror) then
    isconvergence=.false.
    msg='Error when call update_and_check_ in nodal chemistry:'
    call add_ (msg,ithnch)
    goto 20
 end if
!%-------------------------------------------------------------------
!% Check anomalous concentrations 
!%------------------------------------------------------------------- 
 if (isanomalous) then
    isconvergence=.false.
	goto 20 
 end if
!%-------------------------------------------------------------------
!% Check if there was convergence in the node 
!%-------------------------------------------------------------------
 if (.not.isconv) isconvergence=.false.
!%-------------------------------------------------------------------
!% Check the maximum iterations in Picard (gamma) was exeeded 
!%-------------------------------------------------------------------
 if (isupmxiter) then
    isupmxitergam=.true.
    exit
 end if 
!%-------------------------------------------------------------------
!% Check the maximum error in unknowns 
!%-------------------------------------------------------------------
 if (mxunk>mxerrorunk) then
    mxerrorunk=mxunk
    namemxunk=name1
    inchmxunk=ithnch
 end if
!%-------------------------------------------------------------------
!% Check the maximum error in residual 
!%-------------------------------------------------------------------
 if (mxres>mxerrorres) then
    mxerrorres=mxres
    namemxres=name2
    inchmxres=ithnch
 end if
 
 
 ipos1=ipos2+1
 
end do
!%-------------------------------------------------------------------
!% Check divergence 
!%-------------------------------------------------------------------
if (mxerrorunk>mxerrorunkold.and..not.isconvergence) then
 isdivergence=.true.
end if
!%-------------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------------
!% Nullify local pointers 
!%-------------------------------------------------------------------
pnchk1 => null ()
pnchk => null ()
res => null ()
delta => null ()
if (iserror) goto 10 
!%-------------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: update_and_check_'
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
subroutine get_numunk_sia_cheproo &
   (this, &
    nunksia, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the number of unknowns for SIA (or SNIA).
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this      ! Type CHEPROO variable 

integer, intent(out)                   :: nunksia   ! 

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
character(len=100)                     :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
nunksia=0
!%------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system'
 goto 10
end if
!%------------------------------------------------------------
!% Get the number of components for SIA
!%------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numbase=nunksia)
if (iserror) goto 10 
!%------------------------------------------------------------
return
 
10 continue 
print *,'****************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_numunk_sia_'
print *, msg
print *,'****************************'
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine mix_nodalchem_cheproo &
   (this, &
    kth, &
    ith, &
    jth, &
    nstep, &
    iserror, &
    ioutput, &
	ithminset, &
	isopenfile)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Mix the ith and jth nodal chemistry and store the results in the kth nodal chemistry object. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)        :: this       ! Type cheproo variable 

integer, intent(in)                    :: ith        ! Global index of nodal chemistry 1

integer, intent(in)                    :: jth        ! Global index of nodal chemistry 2

integer, intent(in)                    :: kth        ! Global index of nodal chemistry for storage

integer, intent(in)                    :: nstep      ! Number of mixing steps

logical, intent(out)                   :: iserror    ! If .true. there was error

integer, intent(in), optional          :: ioutput    ! Unit of output 

integer, intent(in), optional          :: ithminset  ! Index of mineral set 

logical, intent(in), optional          :: isopenfile ! If true, open the file 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
msg
logical                                :: &
haveioutput, &
haveithminset, &
haveisopenfile, &
isconvergence, &
iswritename 
real*8                                 :: &
p1, &
p2, &
itime1, &
itime2, &
istep 
real*8, pointer                        :: &
alpha(:) => null ()
integer                                :: &
i, &
j
type(t_nodalchemistry), pointer        :: &
nch1 => null (), &
nch2 => null (), &
nch3 => null ()
character(len=100)                     :: &
name 
type(t_pnodalchemistry), pointer       :: &
ptrnch(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system'
 goto 10
end if
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
haveioutput=present(ioutput)
haveithminset=present(ithminset)
haveisopenfile=present(isopenfile)
!%------------------------------------------------------------
if (haveithminset) then
 if (ithminset<=0.or.ithminset>this%numminset) then
  msg='Error in mineral set index'
  goto 10  
 end if 
end if 
!%------------------------------------------------------------
!% Check the nodal chemistry indices 
!%------------------------------------------------------------
if ((ith<=0.or.ith>this%numnch) &
             .or. &
    (jth<=0.or.jth>this%numnch) &
	         .or. &
    (kth<=0.or.kth>this%numnch)) then
	msg='Error in nodal chemistry index'
	goto 10
end if 	 
!%------------------------------------------------------------
!% Allocate and create local nodal chemistry objects 
!%------------------------------------------------------------
allocate (nch1) 
allocate (nch2)
allocate (nch3)
call create_ (nch1)
call create_ (nch2)
call create_ (nch3)
nch1 = this%pnodalchemk1(ith)%ptr 
nch2 = this%pnodalchemk1(jth)%ptr
call set_ (nch1,iserror,name='End member 1')
if (iserror) goto 20
call set_ (nch2,iserror,name='End member 2')
if (iserror) goto 20
!%------------------------------------------------------------ 
if (haveioutput) then
  allocate(ptrnch(nstep+1))
  call check_pointer_ (alpha,nstep+1,.true.)
  if ((.not.haveisopenfile).or.(haveisopenfile.and.isopenfile)) then  ! If true then open the file 
     name=''
     call add_ (name,ioutput)
     name(6:14)='_mix.dat'
     open(unit=ioutput,file=name(1:14),status='unknown')
  end if 
  call write_signature_ (this,ioutput)
  write(ioutput,1) '*****************************************************************'
  write(ioutput,*) 'Beginning mixing end members calculations'
  write(ioutput,*) 'Number of mixing steps:',nstep 
!%------------------------------------------------------------   
  if (haveithminset) then 
   write(ioutput,1) '*****************************************************************'
   write(ioutput,*) 'WARNING!!!: EACH MIXING STEP IS EQUILIBRATED WITH MINERAL PHASES:'
   do i=1,this%minset(ithminset)%nummin
    write(ioutput,*) this%minset(ithminset)%namemin(i)
   end do 
   write(ioutput,1) '*****************************************************************'
  end if 
!%------------------------------------------------------------   
!% Write end members 1 and 2
!%------------------------------------------------------------   
  call write_ (nch1,ioutput,iserror)
  if (iserror) goto 20
  call write_ (nch2,ioutput,iserror)
  if (iserror) goto 20
  if (haveithminset) then 
    do j=1,this%minset(ithminset)%nummin
     name = this%minset(ithminset)%namemin(j)
     call equilibrate_(nch1,name,iserror,isconvergence=isconvergence)
     if (.not.isconvergence.or.iserror) goto 20 
    end do 
  end if
  p1=1.0d0
  p2=0.0d0
  write(ioutput,1) '*****************************************************************'
  write(ioutput,2) 'MIX =',p1,'MEMBER1','+',p2,'MEMBER2'
  write(ioutput,1) '*****************************************************************'
  name='Mixing step'
  call add_ (name,0)
  call set_ (nch1,iserror,name=name)
  if (iserror) goto 20
  call write_ (nch1,ioutput,iserror)
  if (iserror) goto 20
end if 
!%------------------------------------------------------------
!% Start cpu time 
!%------------------------------------------------------------
call cpu_time (itime1)
!%------------------------------------------------------------
!% Start mixing steps 
!%------------------------------------------------------------
istep=0
do i=0,nstep 
 istep=istep+1
 
 p1 = (real(nstep)-real(i))/real(nstep)
 p2 = real(i)/real(nstep)
 alpha(istep)=p1
 nch1 = this%pnodalchemk1(ith)%ptr 
 nch2 = this%pnodalchemk1(jth)%ptr
 
 
 nch3 = p1 * nch1  +  p2 * nch2

!%------------------------------------------------------------ 
!% Equilibrate with mineral species (optional)
!%------------------------------------------------------------ 
 if (haveithminset) then 
  do j=1,this%minset(ithminset)%nummin
   name = this%minset(ithminset)%namemin(j)
   call equilibrate_(nch3,name,iserror,isconvergence=isconvergence)
   if (.not.isconvergence.or.iserror) goto 20 
  end do 
 end if 
!%------------------------------------------------------------ 
 if (haveioutput) then
  write(ioutput,1) '*****************************************************************'
  write(ioutput,2) 'MIX =',p1,'MEMBER1','+',p2,'MEMBER2'
  write(ioutput,1) '*****************************************************************'
  name='Mixing step'
  call add_ (name,i)
  call set_ (nch3,iserror,name=name)
  if (iserror) goto 20
  call write_ (nch3,ioutput,iserror)
  if (iserror) goto 20
  allocate(ptrnch(istep)%ptr)
  call create_ (ptrnch(istep)%ptr) 
  ptrnch(istep)%ptr=nch3
 end if 
!%------------------------------------------------------------ 
end do 
!%------------------------------------------------------------
!% Finishing cpu time 
!%------------------------------------------------------------
call cpu_time (itime2)
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
this%pnodalchemk1(kth)%ptr=nch3
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
if (haveioutput) then
  write(ioutput,1) '*****************************************************************'
  write(ioutput,*) 'Time spent:',itime2-itime1,'[s]'
  write(ioutput,1) '*****************************************************************'
  write(ioutput,1) '*****************************************************************'
  write(ioutput,*) 'Concentrations vs. mixing steps'
  write(ioutput,1) '-----------------------------------------------------------------'
  do i=1,nstep+1
   if (i==1) then
      iswritename=.true.
   else
      iswritename=.false.
   end if 
   call write_sps_(ptrnch(i)%ptr,ioutput,'concentration',this%output%namesp,this%output%numsp,iserror,realitem=alpha(i),iswritename=iswritename)
   call destroy_ (ptrnch(i)%ptr)
   deallocate(ptrnch(i)%ptr)
   if (iserror) goto 20 
  end do
  deallocate (ptrnch)
!%------------------------------------------------------------
!% Close the unit file 
!%------------------------------------------------------------
  if ((.not.haveisopenfile).or. &
    (haveisopenfile.and.isopenfile)) then 
    close (unit=ioutput)  
  end if
end if 
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Destroy and deallocate local variables and pointers 
!%------------------------------------------------------------
call destroy_ (nch1)
call destroy_ (nch2)
call destroy_ (nch3)
deallocate (nch1) 
deallocate (nch2)
deallocate (nch3)
call check_pointer_ (alpha,1,.false.)
if (iserror) goto 10
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
return

10 continue 
print *,'****************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: mix_nodalchem_'
print *, msg
print *,'****************************'
iserror=.true.
return
1 format (a65)  
2 format (a5,1x,f6.2,1x,a7,1x,a1,1x,f6.2,1x,a7)  
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine mix_desimoni_cheproo &
   (this, &
    ith, &
    jth, &
    nstep, &
	nsp, &
	nreact, &
	cmix, &
	secder, &
    iserror, &
    ioutput, &
	isopenfile)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Calculates the concentrations of all concentrations as a
!                 function of the mixing ratio (alpha) of two end members
!                 as well as the second derivatives (secder) with respect
!                 to the mixing ratio: 
!                     Set * secder = d^2c / dalpha^2
!                     where Set is the transpose of the equilibrium 
!                     stoichiometric matrix

!   $Global variables 

type (t_cheproo), intent(inout)        :: this

integer, intent(in)                    :: ith, jth    ! identifiers of end member 1 and 2

integer, intent(in)                    :: nstep       ! Number of mixing ratio steps

integer, intent(out)                   :: nsp          ! Number of species

integer, intent(out)                   :: nreact         ! Number of species

real*8, pointer                        :: cmix(:,:)

real*8, pointer                        :: secder(:,:)

logical, intent(out)                   :: iserror

integer, intent(in), optional          :: ioutput ! Unit of output 

logical, intent(in), optional          :: isopenfile ! If true, open the file 

!  $Local variables 

character(len=100)                     :: &
msg,name

logical                                :: &
haveioutput, &
haveisopenfile

type(t_nodalchemistry), pointer        :: &
nch1 => null (), &
nch2 => null ()

type(t_nodalchemistry), pointer, dimension(:) :: &
nchmix => null ()

real*8, pointer        :: cminus(:),c(:),cplus(:),d2c_dalpha2(:)

integer               :: i, j

real*8                :: alpha, itime1, itime2

!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%------------------------------------------------------------
haveioutput=present(ioutput)
haveisopenfile=present(isopenfile)
!%------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system'
 goto 10
end if
!%------------------------------------------------------------
if ((ith<=0.or.ith>this%numnch) &
             .or. &
    (jth<=0.or.jth>this%numnch)) then
	msg='Error in nodal chemistry index'
	goto 10
end if 	 
!%------------------------------------------------------------
call get_chem_info_(this%pchemsys(1)%ptr,iserror,numreact=nreact,numsp=nsp)
allocate(nch1)
allocate(nch2)
call create_(nch1)
call create_(nch2)
nch1 = this%pnodalchemk1(ith)%ptr
nch2 = this%pnodalchemk1(jth)%ptr
allocate(nchmix(0:nstep))
allocate(cmix(1:nstep-1,nsp))
allocate(secder(1:nstep-1,nreact))
allocate(d2c_dalpha2(nsp))
!%------------------------------------------------------------
if (haveioutput) then
  if ((.not.haveisopenfile).or.(haveisopenfile.and.isopenfile)) then  ! If true then open the file 
     name=''
     call add_ (name,ioutput)
     name(6:14)='_mix.dat'
     open(unit=ioutput,file=name(1:14),status='unknown')
  end if 
  call write_signature_ (this,ioutput)
  write(ioutput,1) '*****************************************************************'
  write(ioutput,*) 'Beginning mixing end members calculations'
  write(ioutput,*) 'Number of mixing steps:',nstep 
!%------------------------------------------------------------   
!% Write end members 1 and 2
!%------------------------------------------------------------   
  call write_ (nch1,ioutput,iserror)
  if (iserror) goto 20
  call write_ (nch2,ioutput,iserror)
  if (iserror) goto 20
end if
!%------------------------------------------------------------
!% Start cpu time 
!%------------------------------------------------------------
call cpu_time (itime1)

!%------------------------------------------------------------
!% Mix end members as a function of mixing ratio alpha
!%------------------------------------------------------------
do i=0,nstep
  alpha = real(i)/real(nstep)
  call create_(nchmix(i))
  nchmix(i) = alpha*nch1 + (1-alpha)*nch2
  if (haveioutput) then
    write(ioutput,1) '*****************************************************************'
    write(ioutput,'(a,e10.5)') 'Mixing ratio =', alpha
    write(ioutput,1) '*****************************************************************'
    if (iserror) goto 20
    call write_ (nchmix(i),ioutput,iserror)
    if (iserror) goto 20
  end if 
enddo

! Calculate second derivatives
do i=1,nstep-1

! Calculate second derivative of concentrations with respect to alpha
  call get_chem_info_(nchmix(i-1),iserror,c=cminus)
  call get_chem_info_(nchmix(i)  ,iserror,c=c     )
  call get_chem_info_(nchmix(i+1),iserror,c=cplus )
  d2c_dalpha2 = (cminus-2*c+cplus)*nstep*nstep
  call get_from_setre_(nchmix(i),d2c_dalpha2,nreact,nsp,iserror)

! Copy
  cmix(i,:)=c
  secder(i,:)=d2c_dalpha2

enddo
!%------------------------------------------------------------
!% Finishing cpu time 
!%------------------------------------------------------------
call cpu_time (itime2)

if (haveioutput) then
  write(ioutput,1) '*****************************************************************'
  write(ioutput,*) 'Time spent:',itime2-itime1,'[s]'
  write(ioutput,1) '*****************************************************************'
!%------------------------------------------------------------
!% Close the unit file 
!%------------------------------------------------------------
  if ((.not.haveisopenfile).or. &
    (haveisopenfile.and.isopenfile)) then 
    close (unit=ioutput)  
  end if
end if 

20 continue 
call destroy_ (nch1)
call destroy_ (nch2)
deallocate (nch1) 
deallocate (nch2)
call check_pointer_(cminus,1,.false.)
call check_pointer_(c,1,.false.)
call check_pointer_(d2c_dalpha2,1,.false.)
call check_pointer_(cplus,1,.false.)
return

10 continue 
print *,'****************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: mix_desimoni_'
print *, msg
print *,'****************************'
iserror=.true.
return

1 format (a65)  
2 format (a5,1x,f6.2,1x,a7,1x,a1,1x,f6.2,1x,a7)  

end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_mxnumunk_dsa_cheproo &
   (this, &
    mxnunk, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return maximum number of primary species defined in cheproo object. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)           :: this

integer, intent(out)                   :: mxnunk

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
 
integer                                :: &
 nunk, &
 inch
character(len=100)                     :: &
 msg 
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------
mxnunk=0
do inch=1,this%numnch
 pnch => this%pnodalchemk1(inch)%ptr
 call get_chem_info_(pnch,iserror,npri=nunk)
 if (iserror) goto 20
 if (nunk>mxnunk) mxnunk=nunk
end do
!%--------------------------------------------------------------
20 continue 
pnch => null ()
if (iserror) goto 10 
!%--------------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_mxnumunk_dsa_'
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
subroutine set_from_cpri_ith_cheproo &
   (this, &
    cpri, &
    dtime, &
    ithnch, &
    isreset, &
    isanomalous, &
    isupmxitergam, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set ith nodal chemistry from cpri 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)        :: this

real*8, intent(in), dimension(:)       :: cpri 

integer, intent(in)                    :: ithnch

real*8, intent(in)                     :: dtime

logical, intent(in)                    :: isreset 

logical, intent(out)                   :: iserror

logical, intent(out)                   :: isanomalous

logical, intent(out)                   :: isupmxitergam
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
 msg 
type(t_nodalchemistry), pointer        :: &
 pnchk1 => null (), &
 pnchk => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
iserror=.false.
msg=''
isupmxitergam=.false.
isanomalous=.false.
!%------------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch) 
 goto 10
end if
!%------------------------------------------------------------
pnchk1 => this%pnodalchemk1(ithnch)%ptr
pnchk => this%pnodalchemk(ithnch)%ptr
call set_from_cpri_ (pnchk1,pnchk,cpri,dtime,isreset,isanomalous,isupmxitergam,iserror)
pnchk1 => null ()
pnchk => null ()
if (iserror) goto 10
!%------------------------------------------------------------
return
 
10 continue 
print *,'**********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: set_from_cpri_ith_'
print *, msg
print *,'**********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_from_cpri_cheproo &
   (this, &
    cpri, &
    nunk, &
    dtime, &
    isreset, &
    isanomalous, &
    isupmxitergam, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set all nodal chemistry objects from cpri (only if 
! DSA solver is performed in the CHEPROO object). 
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)                :: this           ! Type CHEPROO variable 

integer, intent(in)                            :: nunk           ! Number of unknowns

real*8, intent(in), dimension(nunk), target    :: cpri           ! Concentration of primary species

real*8, intent(in)                             :: dtime          ! Time increment 

logical, intent(in)                            :: isreset 

logical, intent(out)                           :: iserror

logical, intent(out)                           :: isanomalous

logical, intent(out)                           :: isupmxitergam
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
 msg 
integer                                :: &
 ithnch, &
 ncomp, &
 ipos1, &
 ipos2, &
 nunkdsa  
type(t_nodalchemistry), pointer        :: &
 pnchk1 => null (), &
 pnchk => null ()
real*8, pointer                        :: &
 pcpri(:) => null ()
integer, pointer                       :: &
 id(:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
isupmxitergam=.false.
isanomalous=.false.
!%------------------------------------------------------------
!%Check if the reactive transport solver performed in the 
!%CHEPROO object 
!%------------------------------------------------------------
if (this%itypesolver/=dsasolver) then
 msg='Error not performed DSA solver'
 goto 10
end if 
!%------------------------------------------------------------
!%Check the total number of unknowns in DSA
!%------------------------------------------------------------
call get_numunk_dsa_ (this,nunkdsa,iserror,nunkdsanch=id)
if (iserror) goto 20
if (nunk/=nunkdsa) then
 msg='Error in number of unknowns for DSA'
 goto 10
end if 
!%------------------------------------------------------------
ncomp=0 
ipos1=1
ipos2=0
do ithnch=1,this%numnch
 pnchk1 => this%pnodalchemk1(ithnch)%ptr
 pnchk => this%pnodalchemk(ithnch)%ptr
 ipos2 = ipos2 + id(ithnch)
 pcpri => cpri(ipos1:ipos2)
 call set_from_cpri_ (pnchk1,pnchk,pcpri,dtime,isreset,isanomalous, &
                      isupmxitergam,iserror)
 if (iserror) goto 20
 ipos1 = ipos2 + 1 
end do
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
pnchk1 => null ()
pnchk => null ()
pcpri => null () 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (id,1,.false.)
if (iserror) goto 10 
!%------------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: set_from_cpri_'
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
subroutine get_duads_ith_cheproo &
   (this, &
    duads, &
    n1, &
    n2, &
    ithnch, &
	isbackwards, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_cheproo), intent(in), target   :: this

real*8, pointer, dimension(:,:)        :: duads

integer, intent(in)                    :: ithnch

integer, intent(out)                   :: n1

integer, intent(out)                   :: n2

logical, intent(in)                    :: isbackwards

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
character(len=100)                     :: &
 msg 
type(t_nodalchemistry), pointer        :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%----------------------------------------------------------
!% If isbackwards=true, then use nodal chemistry list in k 
!%----------------------------------------------------------
if (isbackwards) then
 pnch => this%pnodalchemk(ithnch)%ptr
else
 pnch => this%pnodalchemk1(ithnch)%ptr
end if 
!%----------------------------------------------------------
call get_duads_(pnch,duads,n1,n2,iserror)
if (iserror) goto 10
!%----------------------------------------------------------
nullify (pnch)  
!%----------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_duads_ith_'
print *,msg
print *,'*************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_cheproo &
   (targetobj, &
    sourceobj)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy CHEPROO object.
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)         :: sourceobj

type (t_cheproo), intent(out)        :: targetobj 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: The  target object must be previously created 
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                              :: &
 i, &
 j, &
 ndim  
!-------------------------------------------------------------------------
!
!   $code
!
targetobj%name=sourceobj%name
targetobj%numminset=sourceobj%numminset
targetobj%numgasset=sourceobj%numgasset
targetobj%numsurfset=sourceobj%numsurfset
targetobj%isderstored=sourceobj%isderstored
targetobj%itypesolver=sourceobj%itypesolver
targetobj%iswriteinfo=sourceobj%iswriteinfo
targetobj%zero=sourceobj%zero
!%-----------------------------------------------------------
! Copy chemical system list
!%-----------------------------------------------------------
targetobj%numchemsys=sourceobj%numchemsys
if (targetobj%numchemsys>0) then
  allocate (targetobj%pchemsys(targetobj%numchemsys))
  call check_pointer_ (targetobj%islockchemsys,targetobj%numchemsys,.true.)
  do i=1,targetobj%numchemsys
     targetobj%islockchemsys(i)=.false.    
 	 targetobj%pchemsys(i)%ptr => sourceobj%pchemsys(i)%ptr
  end do
end if
!%-----------------------------------------------------------
! Copy nodal chemistry lists
!%-----------------------------------------------------------
if (targetobj%numnch>0) then
  do i=1,targetobj%numnch
   call destroy_ (targetobj%pnodalchemk(i)%ptr)
   call destroy_ (targetobj%pnodalchemk1(i)%ptr)
   deallocate (targetobj%pnodalchemk(i)%ptr)
   deallocate (targetobj%pnodalchemk1(i)%ptr)
  end do
  deallocate(targetobj%pnodalchemk)
  deallocate(targetobj%pnodalchemk1)
end if
targetobj%numnch=sourceobj%numnch
if (targetobj%numnch>0) then
 allocate(targetobj%pnodalchemk(targetobj%numnch))
 allocate(targetobj%pnodalchemk1(targetobj%numnch))
 do i=1,targetobj%numnch
   allocate (targetobj%pnodalchemk(i)%ptr)
   allocate (targetobj%pnodalchemk1(i)%ptr)
   call create_ (targetobj%pnodalchemk(i)%ptr)
   call create_ (targetobj%pnodalchemk1(i)%ptr)
   targetobj%pnodalchemk(i)%ptr=sourceobj%pnodalchemk(i)%ptr
   targetobj%pnodalchemk1(i)%ptr=sourceobj%pnodalchemk1(i)%ptr
 end do
end if
!%------------------------------------------------------------
if (targetobj%numw>0) then
  do i=1,targetobj%numw
   call destroy_ (targetobj%pwater(i)%ptr)
   deallocate (targetobj%pwater(i)%ptr)
  end do
  deallocate(targetobj%pwater)
end if
targetobj%numw=sourceobj%numw
if (targetobj%numw>0) then
 allocate(targetobj%pwater(targetobj%numw))
 do i=1,targetobj%numw
   allocate (targetobj%pwater(i)%ptr)
   call create_ (targetobj%pwater(i)%ptr)
   targetobj%pwater(i)%ptr=sourceobj%pwater(i)%ptr
 end do
end if
!%------------------------------------------------------------
!% Copy reactive multi-rate mass transfer objects 
!%------------------------------------------------------------
if (targetobj%numrmrmtbox>0) then
  do i=1,targetobj%numrmrmtbox
   call destroy_ (targetobj%prmrmtboxk1(i)%ptr)
   call destroy_ (targetobj%prmrmtboxk(i)%ptr)
   deallocate (targetobj%prmrmtboxk1(i)%ptr)
   deallocate (targetobj%prmrmtboxk(i)%ptr)
  end do
  deallocate(targetobj%prmrmtboxk1)
  deallocate(targetobj%prmrmtboxk)
end if
targetobj%numw=sourceobj%numrmrmtbox
if (targetobj%numrmrmtbox>0) then
 allocate(targetobj%prmrmtboxk1(targetobj%numrmrmtbox))
 allocate(targetobj%prmrmtboxk(targetobj%numrmrmtbox))
 do i=1,targetobj%numrmrmtbox
   allocate (targetobj%prmrmtboxk1(i)%ptr)
   allocate (targetobj%prmrmtboxk(i)%ptr)
   call create_ (targetobj%prmrmtboxk1(i)%ptr)
   call create_ (targetobj%prmrmtboxk(i)%ptr)
   targetobj%prmrmtboxk1(i)%ptr=sourceobj%prmrmtboxk1(i)%ptr
   targetobj%prmrmtboxk(i)%ptr=sourceobj%prmrmtboxk(i)%ptr
 end do
end if
!%------------------------------------------------------------
!% Copy mineral set 
!%------------------------------------------------------------
if(sourceobj%numminset>0) then
 if(associated(targetobj%minset)) deallocate(targetobj%minset)
 allocate(targetobj%minset(targetobj%numminset))
 do i=1,targetobj%numminset
   targetobj%minset(i)%name=sourceobj%minset(i)%name
   targetobj%minset(i)%nummin=sourceobj%minset(i)%nummin
   targetobj%minset(i)%namemin => null ()
   allocate(targetobj%minset(i)%namemin(targetobj%minset(i)%nummin))
   targetobj%minset(i)%unitcmin => null ()
   allocate(targetobj%minset(i)%unitcmin(targetobj%minset(i)%nummin))
   targetobj%minset(i)%unitareamin => null ()
   allocate(targetobj%minset(i)%unitareamin(targetobj%minset(i)%nummin))
   targetobj%minset(i)%cmin => null ()
   allocate(targetobj%minset(i)%cmin(targetobj%minset(i)%nummin))
   targetobj%minset(i)%areamin => null ()
   allocate(targetobj%minset(i)%areamin(targetobj%minset(i)%nummin))
   targetobj%minset(i)%cmin=sourceobj%minset(i)%cmin
   targetobj%minset(i)%namemin=sourceobj%minset(i)%namemin
   targetobj%minset(i)%unitcmin=sourceobj%minset(i)%unitcmin
   targetobj%minset(i)%unitareamin=sourceobj%minset(i)%unitareamin
   targetobj%minset(i)%areamin=sourceobj%minset(i)%areamin
  end do
end if
!%------------------------------------------------------------
!% Copy gas set 
!%------------------------------------------------------------
if(sourceobj%numgasset>0) then
 if(associated(targetobj%gasset)) deallocate(targetobj%gasset)
 allocate(targetobj%gasset(targetobj%numgasset))
 do i=1,targetobj%numgasset
   targetobj%gasset(i)%name=sourceobj%gasset(i)%name
   targetobj%gasset(i)%numgas=sourceobj%gasset(i)%numgas
   targetobj%gasset(i)%namegas => null ()
   allocate(targetobj%gasset(i)%namegas(targetobj%gasset(i)%numgas))
   targetobj%gasset(i)%unitcgas => null ()
   allocate(targetobj%gasset(i)%unitcgas(targetobj%gasset(i)%numgas))
   targetobj%gasset(i)%cgas => null ()
   allocate(targetobj%gasset(i)%cgas(targetobj%gasset(i)%numgas))
   targetobj%gasset(i)%cgas=sourceobj%gasset(i)%cgas
   targetobj%gasset(i)%namegas=sourceobj%gasset(i)%namegas
   targetobj%gasset(i)%unitcgas=sourceobj%gasset(i)%unitcgas
 end do
end if
!%------------------------------------------------------------
!% Copy surface set 
!%------------------------------------------------------------
if(sourceobj%numsurfset>0) then
 if(associated(targetobj%surfset)) deallocate(targetobj%surfset)
 allocate(targetobj%surfset(targetobj%numsurfset))
 do i=1,targetobj%numsurfset
   targetobj%surfset(i)%name=sourceobj%surfset(i)%name
   targetobj%surfset(i)%numsurf=sourceobj%surfset(i)%numsurf
   if (targetobj%surfset(i)%numsurf>0) then
    allocate (targetobj%surfset(i)%surf(targetobj%surfset(i)%numsurf))
    do j=1,targetobj%surfset(i)%numsurf
     targetobj%surfset(i)%surf(j)%name=sourceobj%surfset(i)%surf(j)%name
	 targetobj%surfset(i)%surf(j)%numsites=sourceobj%surfset(i)%surf(j)%numsites
	 ndim=size(sourceobj%surfset(i)%surf(j)%txoh)
	 allocate(targetobj%surfset(i)%surf(j)%txoh(ndim))
	 allocate(targetobj%surfset(i)%surf(j)%capint(ndim))
	 allocate(targetobj%surfset(i)%surf(j)%capext(ndim))
	 allocate(targetobj%surfset(i)%surf(j)%spsurfarea(ndim))
     targetobj%surfset(i)%surf(j)%txoh=sourceobj%surfset(i)%surf(j)%txoh
     targetobj%surfset(i)%surf(j)%capint=sourceobj%surfset(i)%surf(j)%capint
     targetobj%surfset(i)%surf(j)%capext=sourceobj%surfset(i)%surf(j)%capext
     targetobj%surfset(i)%surf(j)%spsurfarea=sourceobj%surfset(i)%surf(j)%spsurfarea
    end do
   end if 
 end do
end if
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_cpri_ith_cheproo &
   (this, &
    cpri, &
    npri, &
    ithnch, &
	isbackwards, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_cheproo), intent(in), target       :: this

real*8, pointer, dimension(:)              :: cpri 

integer, intent(out)                       :: npri

logical, intent(in)                        :: isbackwards

logical, intent(out)                       :: iserror 

integer, intent(in)                        :: ithnch
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
type(t_nodalchemistry), pointer    :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%----------------------------------------------------------
if (ithnch<=0.or.ithnch>this%numnch) then
 msg='Error in nodal chemistry index, index='
 call add_ (msg,ithnch)
 goto 10
end if
!%----------------------------------------------------------
if (isbackwards) then
 pnch => this%pnodalchemk(ithnch)%ptr
else
 pnch => this%pnodalchemk1(ithnch)%ptr
end if 
!%----------------------------------------------------------
call get_chem_info_ (pnch,iserror,cpri=cpri,npri=npri)
if (iserror) goto 10
!%----------------------------------------------------------
pnch => null ()
!%----------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_cpri_ith_'
print *,msg
print *,'**********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_caq_water_cheproo &
   (this, &
    caq, &
    ithw, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the concentration for all the aqueous species defined by the water "ithw"
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)            :: this

real*8, pointer, dimension(:)           :: caq 

integer, intent(in)                     :: ithw         ! water index
                                        
logical, intent(out)                    :: iserror      ! iserror=true, then there was an error 
 
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
 numsp
character(len=100)                    :: &
 msg 
type(t_chemicalsystem),  pointer    :: pchemsys
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
if (ithw<=0.or.ithw>this%numw) then
 msg='Error in water index, index='
 call add_ (msg,ithw) 
 goto 10
end if
!%------------------------------------------------------------

call get_chem_info_ &
   (this%pwater(ithw)%ptr, &
    iserror, &
    numsp=numsp)
if (iserror) goto 10

call get_chem_info_ &
   (this%pwater(ithw)%ptr, &
    iserror, &
    caq=caq,chemsys=pchemsys)
if (iserror) goto 10

 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_cmob_aq_water_'
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
function scalarbycheproo_cheproo &
  (scalar, &
   this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Multiply concentrations in nodal chemistry for scalar value. 
!
!   $Arguments:
!
 
type(t_cheproo), intent(in)           :: this                    ! Type cheproo variable. 

type(t_cheproo)                       :: scalarbycheproo_cheproo ! Type cheproo variable. 

real*8, intent(in)                    :: scalar                  ! Scalar. 
 
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
 iserror 
character(len=100)       :: &
 msg
integer                  :: &
 ithnch
type(t_nodalchemistry), pointer :: &
 pnch1 => null (), &
 pnch2 => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!
!%------------------------------------------------------------------------
msg='' 
!%------------------------------------------------------------------------
!% Create the cheproo object  
!%------------------------------------------------------------------------
call create_ (scalarbycheproo_cheproo)
!%------------------------------------------------------------------------
!% Copy the cheproo object 
!%------------------------------------------------------------------------
scalarbycheproo_cheproo = this
!%------------------------------------------------------------------------
!% Multiply the scalar by the nodal chemistry objects in k+1
!%------------------------------------------------------------------------
do ithnch=1,this%numnch

  pnch1 => scalarbycheproo_cheproo%pnodalchemk1(ithnch)%ptr
  pnch2 => this%pnodalchemk1(ithnch)%ptr
  
  pnch1 = scalar * pnch2

end do
!%------------------------------------------------------------------------
!% Nullify local pointers
!%------------------------------------------------------------------------
pnch1 => null ()
pnch2 => null ()
!%------------------------------------------------------------------------
return 
!%------------------------------------------------------------------------
20 continue 
print *,'*************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: (*)'
print *, msg
print *,'*************************'
iserror=.true.
return
end function
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_compzone_cheproo &
   (this, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Check the definition of the chemical system for each
! nodal chemistry object.
!
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)     :: this       ! Type CHEPROO variable.

logical, intent(out)                :: iserror    ! iserror=true, then there was an error.
 
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
 i
character(len=100)                  :: &
 msg
logical                             :: &
 iscompzchanged 
type(t_nodalchemistry), pointer     :: &
 pnchk => null (), &
 pnchk1 => null () 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not associated chemical system'
 goto 10
end if
!%-------------------------------------------------------------
!% Check the number of nodal chemistry objects. 
!%-------------------------------------------------------------
if (this%numnch==0) then
 msg='Warning, number of nodal chemistry objects = 0'
 goto 10
end if
!%-------------------------------------------------------------
do i=1,this%numnch
   pnchk1 => this%pnodalchemk1(i)%ptr
   pnchk => this%pnodalchemk(i)%ptr 
   call check_compzone_(pnchk1,iserror,pnchk)
   if (iserror) goto 20
end do
!%-------------------------------------------------------------
20 continue
!%-------------------------------------------------------------
!% Nullify local pointers.
!%-------------------------------------------------------------
pnchk1 => null ()
pnchk => null ()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: check_pointer_'
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
subroutine update_cheproo &
   (this, &
    iserror, &
    porosity)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update different chemical/physical parameters (e.g. porosity). 
!
!   $Arguments:
!
 
type(t_cheproo), intent(in)                                      :: this

logical, intent(out)                                             :: iserror 

real*8, intent(inout), dimension(this%numnch), optional, target  :: porosity
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)               :: &
 msg
integer                          :: &
 ithnch
logical                          :: &
 haveporosity 
type(t_nodalchemistry), pointer  :: &
 pnchk1 => null (), &
 pnchk => null ()
real*8, pointer                  :: &
 por => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
!%Check optional arguments
!%------------------------------------------------------------
haveporosity=present(porosity) 
!%------------------------------------------------------------
!% Upcate the porosity 
!%------------------------------------------------------------
if (haveporosity) then
  do ithnch=1,this%numnch
     pnchk1 => this%pnodalchemk1(ithnch)%ptr
     pnchk => this%pnodalchemk(ithnch)%ptr
     por => porosity(ithnch)
     call update_(pnchk1,por,pnchk,iserror)
     if (iserror) goto 10
  end do
  pnchk1 => null ()
  pnchk => null ()
  por => null ()
end if
!%------------------------------------------------------------
return
 
10 continue 
print *,'*****************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: update_'
print *, msg
print *,'*****************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%***************PRIVATE SUBROUTINES**************************
!%************************************************************
 
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
integer        :: status
 
logical                        :: iserror
 
 
call read_xml_loc (name,attributes,iserror=iserror)
 
 
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc (name,attributes,this,msg,iserror)
implicit none 
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
character(len=*), intent(in)   :: name

type(dictionary_t), intent(in) :: attributes

type(t_cheproo), optional      :: this

logical, optional              ::iserror

character(len=*), optional     ::msg 
 
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
 n, &
 status, &
 i, &
 sum
real*8, save   :: &
 zero
logical, save  :: &
 bewz, &
 beminz, &
 benode, &
 bepri, &
 besurfz, &
 begasz, &
 beoutsp, &
 beoutcomp, &
 bechemsys, &
 isderstored, &
 iswriteinfo
real*8, pointer, save:: &
 ctotwz(:,:), &
 cguesswz(:,:), &
 concminz(:,:), &
 areaminz(:,:), &
 tempwz(:), &
 txoh(:,:,:), &
 capint(:,:,:), &
 capext(:,:,:), &
 spsurfarea(:,:,:), &
 concgasz(:,:), &
 temptask(:,:), &
 ntimetask(:)
character(len=100), pointer, save:: &
 namew(:), &
 namespwz(:,:), &
 namespminz(:,:), &
 constraintwz(:,:), &
 unitwz(:,:), &
 unitconcminz(:,:), &
 unitareaminz(:,:), &
 namesurf(:,:), &
 namespgasz(:,:), &
 unitgasz(:,:), &
 nameminset(:), &
 namegasset(:), &
 namesurfset(:), &
 nameoutsp(:), &
 nameouttot(:), &
 unitoutsp(:), &
 nametask(:), &
 equiltask(:)
character(len=100), save:: &
 namefilechemsys, &
 namechp, &
 typesolver
integer, pointer, save:: &
 iconwz(:,:), &
 numsurf(:), &
 numspminz(:), &
 numspwz(:), &
 numsite(:,:), &
 idwz(:), &
 idminz(:), &
 idsurfz(:), &
 numspgasz(:), &
 idgasz(:), &
 nchtask(:), &
 wtask(:), &
 nsteptask(:), &
 outputtask(:), &
 settask(:)
integer, save :: &
 itask, &
 isp, &
 iwz, &
 iminset, &
 isurfset, &
 isurf, &
 igasset, &
 noutsp, &
 nouttot, &
 nsp, &
 isite
integer                  :: &
 numaqcomp, &
 numminsp, &
 numgassp, &
 numsurfprisp, &
 j
integer:: &
 integerlocal(1)
real*8 :: &
 reallocal(1), &
 cputime 
character(len=100), pointer:: &
 nameaqcomp(:), &
 nameminsp(:), &
 namegassp(:), &
 namesurfprisp(:), &
 namesp(:)
character(len=100)             :: &
 id
logical                        :: &
 havethis, &
 be, &
 ok, &
 isrepeated, &
 besp, &
 isconvergence 
integer, pointer               :: &
 idaqprisp(:)
integer, parameter             :: &
 mxdim=100
character(len=100)             :: &
 msgloc 
!-------------------------------------------------------------------------
!
!   $code
!

 

!%---------------------------------------------------------
havethis = present (this)
call lowercase (name)
!%---------------------------------------------------------
if (havethis) then
!%---------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------
!% Read chemical system object 
!%---------------------------------------------------------
if (bechemsys) then
    this%numchemsys=1
	call check_pointer_ (this%islockchemsys,this%numchemsys,.true.)
	this%islockchemsys(1)=.true. 
	allocate(this%pchemsys(this%numchemsys))
	allocate(this%pchemsys(this%numchemsys)%ptr)
	call create_ (this%pchemsys(this%numchemsys)%ptr)
	call read_xml_ (this%pchemsys(this%numchemsys)%ptr,namefilechemsys, iserror)
    if (iserror) then
       msg='Error when calling read_xml_'
       goto 10
    end if
else
    msg='Error, not defined chemical system in cheproo file'
    goto 10
end if
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
select case (typesolver)
case ('sia','SIA','Sia')
  this%itypesolver=2
  if (zero==0.0d0) then
   this%zero=zerosia
  else
   this%zero=zero
  end if
case ('snia','SNIA','Snia')
  this%itypesolver=1
case ('dsa','DSA','Dsa')
  this%itypesolver=3
  if (zero==0.0d0) then
   this%zero=zerodsa
  else
   this%zero=zero
  end if
case default
  msg='Error in reactive transport solver'
  call add_ (msg,typesolver)
  goto 10
end select
!%--------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
this%name=namechp
this%isderstored=isderstored
this%iswriteinfo=iswriteinfo
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
  if (iwz>0) then ! If there are water set
    this%numw=iwz
      allocate (this%pwater(this%numw))
  else
    msg='Error, not defined waters type'
      goto 10
  end if
 
 
   print *,'==========> Starting speciation of water types'
 
   call check_pointer_ (nameminsp,mxdim,.true.)
   call check_pointer_ (namegassp,mxdim,.true.)
   
   do i=1,iwz
      allocate (this%pwater(i)%ptr)
      call create_ (this%pwater(i)%ptr)
      numaqcomp=numspwz(i)
      call check_pointer_ (nameaqcomp,numaqcomp,.true.)
      nameaqcomp(1:numaqcomp)=namespwz(i,1:numspwz(i))
      numminsp=0
      numgassp=0
      numsurfprisp=0
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
      do j=1,numspwz(i)
        select case (iconwz(i,j))
        case (4)
          numminsp=numminsp+1
          nameminsp(numminsp)=constraintwz(i,j)
        case (5)
          numgassp=numgassp+1
          namegassp(numgassp)=constraintwz(i,j)
        end select
      end do
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
      call set_(this%pwater(i)%ptr, &
      iserror, &
      name=namew(i), &
      chemsys=this%pchemsys(1)%ptr, &
      isderstored=isderstored)
      
      if (iserror) then
       msg='Error when calling set_'
       goto 10
      end if
      
      id=namew(i)
	  call lastletter_(j,id)
	  print *,'====> Specia water:',id(1:j)
	  
      call set_from_solution_type_ &
       (this%pwater(i)%ptr, &
        namespwz(i,1:numspwz(i)), &
        ctotwz(i,1:numspwz(i)), &
        cguesswz(i,1:numspwz(i)), &
        unitwz(i,1:numspwz(i)), &
        iconwz(i,1:numspwz(i)), &
        constraintwz(i,1:numspwz(i)), &
        numspwz(i), &
        tempwz(i), &
	    isconvergence, &
		cputime, &
        iserror)
    if (iserror) then
      msg='Error when call set_from_solution_type_'
      goto 10
    end if
	if (.not.isconvergence) then
      msg='Convergence problems in water:'
      id=namew(i)
	  call lastletter_(j,id)
	  call add_ (msg,id(1:j))
      goto 10
    end if

    call write_ (this%pwater(i)%ptr,6,iserror)
    
    if (iserror) then
      msg='Error when call write_'
      goto 10
    end if
 
 
   end do
 
   print *,'==========> Finishing speciation of water types'
 
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
   if(iminset>0) then
    this%numminset=iminset
    allocate(this%minset(iminset))
    do i=1,this%numminset
     this%minset(i)%nummin=numspminz(i)
     this%minset(i)%name=nameminset(i)
     allocate(this%minset(i)%namemin(numspminz(i)))
     allocate(this%minset(i)%cmin(numspminz(i)))
     allocate(this%minset(i)%areamin(numspminz(i)))
     allocate(this%minset(i)%unitcmin(numspminz(i)))
     allocate(this%minset(i)%unitareamin(numspminz(i)))
     this%minset(i)%namemin=namespminz(i,1:numspminz(i))
     this%minset(i)%cmin=concminz(i,1:numspminz(i))
     this%minset(i)%areamin=areaminz(i,1:numspminz(i))
     this%minset(i)%unitcmin=unitconcminz(i,1:numspminz(i))
     this%minset(i)%unitareamin=unitareaminz(i,1:numspminz(i))
    end do
   end if
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
   if(igasset>0) then
    this%numgasset=igasset
    allocate(this%gasset(igasset))
    do i=1,this%numgasset
     this%gasset(i)%numgas=numspgasz(i)
     this%gasset(i)%name=namegasset(i)
     allocate(this%gasset(i)%namegas(numspgasz(i)))
     allocate(this%gasset(i)%cgas(numspgasz(i)))
     allocate(this%gasset(i)%unitcgas(numspgasz(i)))
     this%gasset(i)%namegas=namespgasz(i,1:numspgasz(i))
     this%gasset(i)%cgas=concgasz(i,1:numspgasz(i))
     this%gasset(i)%unitcgas=unitgasz(i,1:numspgasz(i))
    end do
   end if
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
   if(isurfset>0) then
    this%numsurfset=isurfset
    allocate(this%surfset(this%numsurfset))
    do i=1,this%numsurfset
     isurf=numsurf(i)
     this%surfset(i)%numsurf=isurf
     this%surfset(i)%name=namesurfset(i)
     allocate(this%surfset(i)%surf(isurf))
      do j=1,isurf
        isite=numsite(i,j)
        this%surfset(i)%surf(j)%name=namesurf(i,j)
        this%surfset(i)%surf(j)%numsites=isite
        allocate (this%surfset(i)%surf(j)%txoh(isite))
        allocate (this%surfset(i)%surf(j)%capint(isite))
        allocate (this%surfset(i)%surf(j)%capext(isite))
        allocate (this%surfset(i)%surf(j)%spsurfarea(isite))
        this%surfset(i)%surf(j)%txoh=txoh(i,j,1:isite)
        this%surfset(i)%surf(j)%capint=capint(i,j,1:isite)
        this%surfset(i)%surf(j)%capext=capext(i,j,1:isite)
        this%surfset(i)%surf(j)%spsurfarea=spsurfarea(i,j,1:isite)
      end do
    end do
   end if
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
  if (nouttot>0) then
    call find_repeated_(nameouttot(1:nouttot),isrepeated,nouttot,id)
    if (isrepeated) then
     msg='Error, component repeated in output tag:'
     call add_ (msg,id)
     goto 10
    end if
    this%output%numcomp=nouttot
    allocate(this%output%namecomp(this%output%numcomp))
    this%output%namecomp=nameouttot(1:nouttot)
    call get_chem_info_ &
     (this%pchemsys(1)%ptr, &
      iserror, &
      idbase=idaqprisp)
 
     if (iserror) then
      msg='Error when calling get_chem_info_'
      goto 10
     end if
 
     do i=1,nouttot
       call get_sp_index_ (this%pchemsys(1)%ptr,nameouttot(i), isp)
       if (isp.eq.0) then
            msg='Error, in output tag, not defined in the chemical system the species'
            call add_ (msg,nameouttot(i))
            goto 10
       end if
       be=.false.
       do j=1,size(idaqprisp)
         if (idaqprisp(j).eq.isp) then
          be=.true.
         end if
       end do
       if (.not.be) then
           call check_pointer_ (idaqprisp,1,.false.)
                msg='Error in output tag, not defined in the chemical system the component'
           call add_ (msg,nameouttot(i))
           goto 10
       end if
     end do
     deallocate (idaqprisp)
  end if
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
  if (noutsp>0) then
    call find_repeated_(nameoutsp(1:noutsp),isrepeated,noutsp,id)
    if (isrepeated) then
     msg='Error, species repeated in output tag:'
     call add_ (msg,id)
     goto 10
    end if
    this%output%numsp=noutsp
    allocate(this%output%namesp(this%output%numsp))
    this%output%namesp=nameoutsp(1:noutsp)
    allocate(this%output%unitsp(this%output%numsp))
    this%output%unitsp=unitoutsp(1:noutsp)
    do i=1,noutsp
      call get_sp_index_ (this%pchemsys(1)%ptr, nameoutsp(i), isp)
      if (isp==0) then
          msg='Error, in output tag not defined in the chemical system the species'
          call add_ (msg,nameoutsp(i))
          goto 10
      end if
    end do
  end if
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
 if (itask>0) then
 
      call create_nodalchem_ (this,itask,iserror)
      if (iserror) goto 10 
 
      do i=1,itask
       print *,'=======> Solving task',i
 
       call init_nodalchem_ &
      (this, &
       iserror, &
       idwaterset=wtask(i:i))
 
       select case (nametask(i))
     case ('add','ADD')
      call reaction_path_ith_ &
      (this, &
       nchtask(i), &
       nsteptask(i), &
       iserror, &
       outputtask(i), &
       settask(i))
     case ('ev/dil','EV/DIL')
 
 
     end select
     print *,'=======> Finishing task',i
    end do
 
end if
!%--------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
 call check_pointer_ (ctotwz,1,1,.false.)
 call check_pointer_ (cguesswz,1,1,.false.)
 call check_pointer_ (concminz,1,1,.false.)
 call check_pointer_ (areaminz,1,1,.false.)
 call check_pointer_ (tempwz,1,.false.)
 call check_pointer_ (txoh,1,1,1,.false.)
 call check_pointer_ (capint,1,1,1,.false.)
 call check_pointer_ (capext,1,1,1,.false.)
 call check_pointer_ (spsurfarea,1,1,1,.false.)
 call check_pointer_ (concgasz,1,1,.false.)
 call check_pointer_ (namespwz,1,1,.false.)
 call check_pointer_ (namespminz,1,1,.false.)
 call check_pointer_ (constraintwz,1,1,.false.)
 call check_pointer_ (unitwz,1,1,.false.)
 call check_pointer_ (unitconcminz,1,1,.false.)
 call check_pointer_ (unitareaminz,1,1,.false.)
 call check_pointer_ (namesurf,1,1,.false.)
 call check_pointer_ (namespgasz,1,1,.false.)
 call check_pointer_ (unitgasz,1,1,.false.)
 call check_pointer_ (nameminset,1,.false.)
 call check_pointer_ (namegasset,1,.false.)
 call check_pointer_ (namesurfset,1,.false.)
 call check_pointer_ (nameoutsp,1,.false.)
 call check_pointer_ (unitoutsp,1,.false.)
 call check_pointer_ (nameouttot,1,.false.)
 call check_pointer_ (iconwz,1,1,.false.)
 call check_pointer_ (numsurf,1,.false.)
 call check_pointer_ (numspminz,1,.false.)
 call check_pointer_ (numspwz,1,.false.)
 call check_pointer_ (numsite,1,1,.false.)
 call check_pointer_ (idwz,1,.false.)
 call check_pointer_ (idminz,1,.false.)
 call check_pointer_ (idsurfz,1,.false.)
 call check_pointer_ (numspgasz,1,.false.)
 call check_pointer_ (idgasz,1,.false.)
 call check_pointer_ (nameaqcomp,1,.false.)
 call check_pointer_ (nameminsp,1,.false.)
 call check_pointer_ (namegassp,1,.false.)
 call check_pointer_ (namesurfprisp,1,.false.)
 call check_pointer_ (namew,1,.false.)
 call check_pointer_ (nametask,1,.false.)
 call check_pointer_ (nchtask,1,.false.)
 call check_pointer_ (wtask,1,.false.)
 call check_pointer_ (ntimetask,1,.false.)
 call check_pointer_ (nsteptask,1,.false.)
 call check_pointer_ (temptask,1,1,.false.)
 call check_pointer_ (equiltask,1,.false.)
 call check_pointer_ (settask,1,.false.)
 call check_pointer_ (outputtask,1,.false.)
 
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
else
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
  select case(name)
!%-----------------------------
  case('cheproo')
!%----------------------------- 
   allocate(ctotwz(mxdim,mxdim))
   allocate(cguesswz(mxdim,mxdim))
   allocate(concminz(mxdim,mxdim))
   allocate(areaminz(mxdim,mxdim))
   allocate(tempwz(mxdim))
   allocate(txoh(mxdim,mxdim,mxdim))
   allocate(capint(mxdim,mxdim,mxdim))
   allocate(capext(mxdim,mxdim,mxdim))
   allocate(spsurfarea(mxdim,mxdim,mxdim))
   allocate(concgasz(mxdim,mxdim))
   allocate(namespwz(mxdim,mxdim))
   allocate(namespminz(mxdim,mxdim))
   allocate(constraintwz(mxdim,mxdim))
   allocate(unitwz(mxdim,mxdim))
   allocate(unitconcminz(mxdim,mxdim))
   allocate(unitareaminz(mxdim,mxdim))
   allocate(namesurf(mxdim,mxdim))
   allocate(namespgasz(mxdim,mxdim))
   allocate(unitgasz(mxdim,mxdim))
   allocate(nameminset(mxdim))
   allocate(namegasset(mxdim))
   allocate(namesurfset(mxdim))
   allocate(nameoutsp(mxdim))
   allocate(unitoutsp(mxdim))
   allocate(nameouttot(mxdim))
   allocate(iconwz(mxdim,mxdim))
   allocate(numsurf(mxdim))
   allocate(numspminz(mxdim))
   allocate(numspwz(mxdim))
   allocate(numsite(mxdim,mxdim))
   allocate(idwz(mxdim))
   allocate(idminz(mxdim))
   allocate(idsurfz(mxdim))
   allocate(numspgasz(mxdim))
   allocate(idgasz(mxdim))
   allocate(nameaqcomp(mxdim))
   allocate(nameminsp(mxdim))
   allocate(namegassp(mxdim))
   allocate(namesurfprisp(mxdim))
   allocate(namew(mxdim))
   allocate(nametask(mxdim))
   allocate(nchtask(mxdim))
   allocate(wtask(mxdim))
   allocate(ntimetask(mxdim))
   allocate(nsteptask(mxdim))
   allocate(temptask(2,mxdim))
   allocate(equiltask(mxdim))
   allocate(settask(mxdim))
   allocate(outputtask(mxdim))
 
   txoh=0.0d0
   capint=0.0d0
   capext=0.0d0
   spsurfarea=0.0d0
   isurfset=0 
   itask=0
   iwz=0
   iminset=0
   igasset=0
   isurf=0
   beoutsp=.false.
   beoutcomp=.false.
   nouttot=0
   noutsp=0
   bechemsys=.false.
 
!%-------------------
   id=' '
   call get_value (attributes,"name", id, status)
   if (status<0) then
    namechp='CHEPROO'
   else
    namechp=id
   end if
!%------------------
   id=''
   call get_value (attributes,"derivatestored", id, status)
   if (status<0) then
    isderstored=.false.
   else
    call lowercase (id)
         if (id.eq.'yes') then
          isderstored=.true.
    else
     isderstored=.false.
    end if
   end if
!%------------------
   id=''
   call get_value (attributes,"writeinfo", id, status)
   if (status<0) then
    iswriteinfo=.false.
   else
    call lowercase (id)
         if (id.eq.'yes') then
          iswriteinfo=.true.
    else
     iswriteinfo=.false.
    end if
   end if
!%-------------------
   id=' '
   call get_value (attributes,"zero", id, status)
   n=0
   reallocal(1)=0.0d0
   call build_data_array (id,reallocal,n)
   zero=reallocal(1)
!%-------------------
   id=' '
   call get_value (attributes,"typesolver", id, status)
   typesolver=id
!%------------------------------
  case('chemicalsystem')
   bechemsys=.true.
!%------------------------------
  case('file')
   if (bechemsys) then
    id=' '
    call get_value (attributes,"name", id, status)
    if (status/=0) then
      msgloc='Error, not defined the chemical system file name'
      goto 20
    end if
    namefilechemsys=id
   end if
!%----------------------------
  case('waters')
!%-----------------------------
   bewz=.true.
!%-----------------------------   
  case('water')
!%-----------------------------
   iwz=iwz+1
   isp=0
   id=' '
   call get_value (attributes,"name", id, status)
   namew(iwz)=id
   id=' '
   call get_value (attributes,"id", id, status)
   n=0
   integerlocal(1)=0
   call build_data_array (id,integerlocal,n)
   idwz(iwz)=integerlocal(1)
   id=' '
   call get_value (attributes,"temp", id, status)
   n=0
   reallocal(1)=0.0d0
   call build_data_array (id,reallocal,n)
   tempwz(iwz)=reallocal(1)
!%-----------------------------
  case('sets')
!%-----------------------------  
    bewz=.false.
!%-----------------------------
  case('set')
    bewz=.false.
    beminz=.false.
    besurfz=.false.
    begasz=.false.
    beoutsp=.false.
    beoutcomp=.false.
    isurf=0
    isp=0
	id='' 
    call get_value (attributes,"type", id, status)
    select case (id)
    case ('mineral')
      beminz=.true.
      iminset=iminset+1
    case ('surface')
      besurfz=.true.
      isurfset=isurfset+1
    case ('gas')
      begasz=.true.
      igasset=igasset+1
    end select
    id=' '
    call get_value (attributes,"id", id, status)
    n=0
    integerlocal(1)=0
    call build_data_array (id,integerlocal,n)
    id=' '
    call get_value (attributes,"name", id, status)
    if (beminz) then
      nameminset(iminset)=id
    else if (besurfz) then
      namesurfset(isurfset)=id
    else if (begasz) then
      namegasset(igasset)=id
    end if
!%---------------
       if (beminz) then
         idminz(iminset)=iminset
!%---------------
     else if (besurfz) then
       idsurfz(isurfset)=isurfset
!%---------------
       else if (begasz) then
       idgasz(igasset)=igasset
     end if
!%---------------
  case('species')
     isp = isp + 1
       id=' '
       call get_value (attributes,"name", id, status)
       call lowercase (id)
!%--------------------
       if (bewz) then
      numspwz(iwz)=isp
      namespwz(iwz,isp)=id
      id=' '
        call get_value (attributes,"icon", id, status)
!%---------------------
          if (status<0) then
             iconwz(iwz,isp)=0
          else
		    select case (id)
            case('conc')
             iconwz(iwz,isp)=0
			case('tot')
             iconwz(iwz,isp)=1
            case('chgbal')
             iconwz(iwz,isp)=2
            case('act')
             iconwz(iwz,isp)=3
            case('eqgas')
             iconwz(iwz,isp)=5
            case('eqmin')
             iconwz(iwz,isp)=4
            case default
             msgloc='Error, not defined constraint type'
            call add_ (msg,id)
             goto 20
           end select
          end if
!%---------------------          
		  id=' '
          call get_value (attributes,"ctot", id, status)
          if (id.eq.' ') then
              msgloc='Error, not defined ctot by component'
              call add_ (msg,namespwz(iwz,isp))
              goto 20
          end if
          n=0
          reallocal(1)=0.0d0
          call build_data_array (id,reallocal,n)
          ctotwz(iwz,isp)=reallocal(1)
          id=' '
          call get_value (attributes,"cguess", id, status)
          if (id.eq.' ') then
             msgloc='Error, not defined cguess by component'
          call add_ (msg,namespwz(iwz,isp))
          goto 20
      end if
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
      cguesswz(iwz,isp)=reallocal(1)
      id=' '
      call get_value (attributes,"constraint", id, status)
      call lowercase (id)
      constraintwz(iwz,isp)=id
      id=' '
      call get_value (attributes,"unit", id, status)
      call lowercase (id)
      unitwz(iwz,isp)=id
      if (id==' ') unitwz(iwz,isp)='m'
!%--------------------
 
     else if (beminz) then
 
        numspminz(iminset)=isp
      namespminz(iminset,isp)=id
      id=' '
        call get_value (attributes,"conc", id, status)
      if (id.eq.' ') then
             msgloc='Error, not defined conc by mineral species'
        call add_ (msg,namespminz(iminset,isp))
                  goto 20
           end if
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
        concminz(iminset,isp)=reallocal(1)
      id=' '
      call get_value (attributes,"area", id, status)
      n=0
      reallocal(1)=0.0d0
      call  build_data_array (id,reallocal,n)
        areaminz(iminset,isp)=reallocal(1)
        if (id.eq.' ')  areaminz(iminset,isp)=1.0d0
        id=' '
        call get_value (attributes,"unitconc", id, status)
      call lowercase (id)
        unitconcminz(iminset,isp)=id
        if (id.eq.' ') unitconcminz(iminset,isp)='m'
        id=' '
        call get_value (attributes,"unitarea", id, status)
      call lowercase (id)
      unitareaminz(iminset,isp)=id
      if (id.eq.' ') unitareaminz(iminset,isp)='m2/m3rock'
!%--------------------
       else if (begasz) then
 
           numspgasz(igasset)=isp
      namespgasz(igasset,isp)=id
      id=' '
           call get_value (attributes,"conc", id, status)
      if (id.eq.' ') then
            msgloc='Error, not defined conc for species'
       call add_ (msg,namespgasz(igasset,isp))
            goto 20
      end if
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
           concgasz(igasset,isp)=reallocal(1)
      id=' '
      call get_value (attributes,"unitconc", id, status)
      unitgasz(igasset,isp)=id
!%--------------------
      else if (beoutsp) then
        noutsp=isp
        id=' '
        call get_value (attributes,"name", id, status)
        if (status<0) then
           msgloc='Error, not defined name by output species'
           goto 20
        end if
        nameoutsp(isp)=id
        id=''
        call get_value (attributes,"unit", id, status)
        if (status<0) then
         unitoutsp(isp)='m'
      else
         unitoutsp(isp)=id
      end if
!%--------------------
      else if (beoutcomp) then
           nouttot=isp
      id=' '
      call get_value (attributes,"name", id, status)
      if(status<0) then
            msgloc='Error, not defined name by output components'
            goto 20
      end if
      nameouttot(isp)=id
 
      end if
!%--------------------------
  case ('surface')
      isite=0
      isurf=isurf+1
      numsurf(isurfset)=isurf
      id=' '
      call get_value (attributes,"name", id, status)
        namesurf(isurfset,isurf)=id
!%--------------------------
  case ('site')
      isite=isite+1
      numsite(isurfset,isurf)=isite
      id=' '
      call get_value (attributes,"txoh", id, status)
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
      txoh(isurfset,isurf,isite)=reallocal(1)
      id=' '
      call get_value (attributes,"capint", id, status)
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
      capint(isurfset,isurf,isite)=reallocal(1)
      id=' '
      call get_value (attributes,"capext", id, status)
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
      capext(isurfset,isurf,isite)=reallocal(1)
      id=' '
      call get_value (attributes,"spsurfarea", id, status)
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
      spsurfarea(isurfset,isurf,isite)=reallocal(1)
!%-------------------------
  case ('output')
      nouttot=0
      noutsp=0
!%-------------------------
  case ('totals')
      beoutsp=.false.
      beoutcomp=.true.
      beminz=.false.
      besurfz=.false.
      begasz=.false.
      besp=.false.
      bewz=.false.
      isp=0
!%-------------------------
  case ('concentrations')
      beoutsp=.true.
      beoutcomp=.false.
      beminz=.false.
      besurfz=.false.
      begasz=.false.
      besp=.false.
      bewz=.false.
      isp=0
!%-------------------------
  case ('task')
      itask=itask+1
      id=''
      call get_value (attributes,"name", id, status)
      nametask(itask)=id
!%---------
      id=''
      call get_value (attributes,"nch", id, status)
      n=0
      integerlocal(1)=1
      call build_data_array (id,integerlocal,n)
      nchtask(itask)=integerlocal(1)
!%---------
      id=''
      call get_value (attributes,"water", id, status)
      n=0
      integerlocal(1)=1
      call build_data_array (id,integerlocal,n)
      wtask(itask)=integerlocal(1)
!%---------
      id=''
      call get_value (attributes,"ntime", id, status)
      n=0
      reallocal(1)=1.0d0
      call build_data_array (id,reallocal,n)
      ntimetask(itask)=reallocal(1)
!%---------
      id=''
      call get_value (attributes,"nstep", id, status)
      n=0
      integerlocal(1)=1
      call build_data_array (id,integerlocal,n)
      nsteptask(itask)=integerlocal(1)
!%---------
      id=''
      call get_value (attributes,"temp1",id, status)
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
      temptask(1,itask)=reallocal(1)
!%---------
      id=''
      call get_value (attributes,"temp2",id, status)
      n=0
      reallocal(1)=0.0d0
      call build_data_array (id,reallocal,n)
      temptask(2,itask)=reallocal(1)
!%---------
      id=''
      call get_value (attributes,"isequilibrate",id, status)
      equiltask(itask)=id
!%---------
      id=''
      call get_value (attributes,"set",id, status)
      n=0
      integerlocal(1)=1
      call build_data_array (id,integerlocal,n)
      settask(itask)=integerlocal(1)
!%---------
      id=''
      call get_value (attributes,"outputunit",id, status)
      n=0
      integerlocal(1)=1
      call build_data_array (id,integerlocal,n)
      outputtask(itask)=integerlocal(1)
!%------------------------
  case default
      msgloc='Error, not recongnized tag'
      call add_ (msg,name)
      goto 20
  end select
 
end if
!%---------------------------------------------------------
return
 
10 iserror=.true.
return
 
20 continue 
print *,'**********************'
print *,'CHEPROO:'
if (havethis) print *,'Name:',this%name
print *,'Service: read_xml_'
print *, msgloc
print *,'**********************'
stop
 
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine solve_react_trp_step_dsa_cheproo &
   (this, &
    iconn, &
    mxconn, &
    nnch, &
    typesto, &
    mxitertr, &
    mxdiverg, &
    mxchcompz, &
    theta, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    nband, &
    ngas, &
    dtime, &
    istrpgas, &
    isupdatejac, &
    caudalboundw, &
    caudalrechw, &
    caudalrechg, &
    caudalboundg, &
    itypeboundw, &
    itypeboundg, &
    gammaboundw, &
    isetboundw, &
    isetrechw, &
    isconvergence, &
	iter, &
    factor, &
    tolunk, &
    tolres, &
    isopensystem, &
    time, &
	mxerror, &
	issolvermrmt, &
	imethod, &
	mxnbox, &
	iconnbox, &
	alpha, &
	phi, &
	iserror, &
	molinout, &
	molkin, &
	iouwinfo)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Solve one reactive transport time step (for DSA)
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)                                :: this            ! Type CHEPROO variable. 

integer, intent(in)                                            :: mxconn          ! Maximum number of connected nodal chemistries. 

integer, intent(in)                                            :: nnch            ! Number of nodal chemistry objects 

integer, intent(in)                                            :: typesto         ! Type of storage in the jacobian matrix. 
  
integer, intent(in)                                            :: mxitertr        ! Maximum number of iterations. 

integer, intent(in)                                            :: mxdiverg        ! Maximum number of divergences. 

integer, intent(in)                                            :: mxchcompz       ! Maximum number of allowed changes in components zones. 

integer, intent(in)                                            :: nband

integer, intent(in)                                            :: ngas            ! Number of gases. 

real*8, intent(in)                                             :: theta           ! Temporal weight  

real*8, intent(in)                                             :: dtime           ! Time increment 

real*8, intent(in)                                             :: factor          ! Factor to update the unknowns. 

real*8, intent(in)                                             :: tolunk          ! Tolerance in unknowns. 

real*8, intent(in)                                             :: tolres          ! Tolerance in the residual 

real*8, intent(in)                                             :: time            ! Time increment 

real*8, intent(in), dimension(nnch,2*nband+1)                  :: efloww_ij       ! Water conductancy matrix (in spars)

real*8, intent(in), dimension(nnch,2*nband+1)                  :: eflowg1_ij

real*8, intent(in), dimension(nnch,2*nband+1)                  :: eflowg2_ij

real*8, intent(in), dimension(nnch)                            :: caudalboundw

real*8, intent(in), dimension(nnch)                            :: caudalrechw

real*8, intent(in), dimension(nnch)                            :: caudalboundg

real*8, intent(in), dimension(nnch)                            :: caudalrechg

real*8, intent(in), dimension(nnch)                            :: gammaboundw

integer, intent(in), dimension(mxconn,nnch)                    :: iconn          ! Indices of connected nodal chemistries. 

integer, intent(in), dimension(nnch)                           :: isetboundw

integer, intent(in), dimension(nnch)                           :: isetrechw

integer, intent(in), dimension(nnch)                           :: itypeboundw

integer, intent(in), dimension(nnch)                           :: itypeboundg

logical, intent(in)                                            :: isupdatejac   ! isupdatejac=true, the jacobian matrix is updated in each Newton-Raphson iteration. 

logical, intent(in)                                            :: istrpgas      ! istrpgas=true, then gasses is transported. 

logical, intent(in)                                            :: isopensystem  ! isopensystem=true, the system is considered open.  

logical, intent(out)                                           :: isconvergence ! isconvergence=true, then there was convergence in the reactive transport calculations. 

integer, intent(out)                                           :: iter          ! Number of iteration for newton raphson 

logical, intent(out)                                           :: iserror       ! iserror=true, then there was an error. 

real*8, intent(in)                                             :: mxerror       ! Maximum error in residual 

real*8, pointer, optional, dimension(:)                        :: molinout

real*8, pointer, optional, dimension(:)                        :: molkin

integer, intent(in), optional                                  :: iouwinfo 

logical, intent(in)                                            :: issolvermrmt      ! .true. if there is error

integer, intent(in)                                            :: imethod         ! Method for to solve reactive multi-rate mass transfer 

integer, intent(in)                                            :: mxnbox

real*8, intent(in), dimension(mxnbox,this%numnch)              :: alpha           ! Alpha for boxes 

real*8, intent(in), dimension(mxnbox,this%numnch)              :: phi             ! Porosity for boxes 

integer, intent(in), dimension(mxnbox+1,this%numnch)           :: iconnbox        ! Connectivity between boxes 
 
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
 msg, &
 namemxunk, &
 namemxres
integer                              :: &
 nunk, &
 nsp, &
 nbandjac, &
 info, &
 ichange, &
 idiverg, &
 inchmxunk, &
 inchmxres, &
 ndim, &
 inch
real*8, pointer                      :: &
 jacobian(:,:) => null (), &
 residual(:) => null (), &
 deltacunk(:) => null (), &
 tindep(:,:) => null ()
integer, pointer                     :: &
 icolt(:) => null ()
real*8                               :: &
 mxerrorunk, &
 mxerrorres, &
 itime1, &
 itime2, &
 dtindep,           & ! Cpu time spending building the independent t
 dtjacres,          & ! Cpu time spending building the jacobian/resi
 dtsolve,           & ! Cpu time spending solving the lineal system
 dtupdate,          & ! Cpu time spending updatting the solution
 dteqmin              ! Cpu time spending computing equilibrium mine
logical                               :: &
 iscompzchanged, &
 upmxchcompz,       & ! .true. if maximum number of component zones
 upmxiter,          & ! .true. if maximum number of iteration is arr
 divergence, &
 upmxdiverg,        & ! .true. if maximum number of divergence is ar
 havemolinout, &
 havemolkin, &
 haveiouwinfo, &
 isanomalous,       & ! .true. if there was anomalous concentrations
 upmxitergam, &
 upmxerrorres, &
 isconvrmrmt    
integer, parameter                   :: &
 siatype=1, &
 dsatype=2
real*8, parameter                    :: &
 r0=0.0d0 
!-------------------------------------------------------------------------
!
!   $code
!

iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Check optional arguments 
!%--------------------------------------------------------------
havemolinout=present(molinout)
havemolkin=present(molkin)
haveiouwinfo=present(iouwinfo)
!%--------------------------------------------------------------
!% Check the number of nodal chemistry objects
!%--------------------------------------------------------------
if (nnch/=this%numnch) then
  msg='Error in number of nodal chemistry objects'
  goto 10
end if
!%--------------------------------------------------------------
!% Check the number of chemical system objects 
!%--------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system in CHEPROO'
 goto 10
end if
!%--------------------------------------------------------------
!% Check the maximum error ion residual
!%--------------------------------------------------------------
if (mxerror<=r0) then
  msg='Error in maximum error in residual, mxerror='
  call add_ (msg,mxerror) 
  goto 10
end if
!%--------------------------------------------------------------
!% Initialice variables
!%--------------------------------------------------------------
upmxdiverg=.false.
upmxiter=.false.
upmxitergam=.false.
isanomalous=.false.
upmxerrorres=.false. 
!%--------------------------------------------------------------
!% Write convergence parameters information
!%--------------------------------------------------------------
if (this%iswriteinfo.and.haveiouwinfo) then
    write (iouwinfo,*) '------------------------------------'// &
                       '------------------------------------'
    write (iouwinfo,*) 'Service: solve_react_trp_step_'
    write (iouwinfo,*) '------------------------------------'// &
                       '------------------------------------'
    write (iouwinfo,*) '   Newton Raphson Information  '
    write (iouwinfo,*) 'Tolerance in unknowns:',tolunk
    write (iouwinfo,*) 'Tolerance in residual:',tolres
	write (iouwinfo,*) 'Mx. error in residual',mxerror
    write (iouwinfo,*) 'Mx. number of divergence',mxdiverg
    write (iouwinfo,*) 'Mx. number of trp. iter.',mxitertr
    write (iouwinfo,*) 'Time:',time
    write (iouwinfo,*) 'Time increment:',dtime
    write (iouwinfo,*) '------------------------------------'// &
                         '------------------------------------'
end if
!%----------------------------------------------------------------
! 1) Build independent term
!%----------------------------------------------------------------
call cpu_time (itime1)
 
call build_indep_term_dsa_ &
   (this, &
    tindep, &
    nnch, &
    nsp, &
    theta, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    nband, &
    ngas, &
    dtime, &
    istrpgas, &
    caudalboundw, &
    caudalrechw, &
    caudalrechg, &
    caudalboundg, &
    itypeboundw, &
    itypeboundg, &
    isetboundw, &
    isetrechw, &
    gammaboundw, &
    iserror)
 
 if (iserror) goto 20
 
call cpu_time (itime2)
dtindep=itime2-itime1
 
!%--------------------------------------------------------------
! 2) Start Newton Raphson iterative process
!%--------------------------------------------------------------
ichange=0
upmxchcompz=.false.
upmxdiverg=.false.
!%--------------------------------------------------------------
do
 
 idiverg=0
 isconvergence=.false.
 iscompzchanged=.false.
 isconvrmrmt=.true. 
 iter=0
!%--------------------------------------------------------------
!%--------------------------------------------------------------
 if (this%iswriteinfo.and.haveiouwinfo) then
   
   write (iouwinfo,*) '------------------------------------'// &
                      '------------------------------------'
   write (iouwinfo,3) 'ichange:',ichange
   write (iouwinfo,1) 'iter','converge','mx.errunk', &
                      'node','comp','mx.errres','node', &
                      'comp','dt(ind)','dt(jac/res)', &
                      'dt(solve)','dt(upd)'
   write (iouwinfo,*) '------------------------------------'// &
                      '------------------------------------'
 
 end if
!%--------------------------------------------------------------
!%--------------------------------------------------------------
do while (.not.isconvergence &
              .and. &
          .not.upmxiter &
              .and. &
          .not.upmxdiverg &
              .and. &
          .not.isanomalous &
              .and. &
          .not.upmxitergam &
              .and. &
          .not.upmxerrorres &
		      .and. &
           .not.iserror)
!%--------------------------------------------------------------
  iter=iter+1
!%--------------------------------------------------------------
! 2.1) Build the jacobian and residual for DSA
!%--------------------------------------------------------------
  call cpu_time (itime1)
 
  call build_jacob_resid_dsa_ &
   (this, &
    jacobian, &
    residual, &
    nunk, &
    nbandjac, &
    iconn, &
	itypeboundw, &
    nnch, &
    mxconn, &
    nband, &
    istrpgas, &
    typesto, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    ngas, &
    theta, &
    tindep, &
    nsp, &
    dtime, &
	issolvermrmt, &
	imethod, &
	mxnbox, &
	iconnbox, &
	alpha, &
	phi, &
	isconvrmrmt, &
	iserror)
 
  call cpu_time (itime2)
  dtjacres=itime2-itime1

  if (iserror) goto 20
!%--------------------------------------------------------------
!% If reactive multi-rate mass transfer is solved and there was
!% not convergence, then exit 
!%--------------------------------------------------------------
  if (issolvermrmt.and..not.isconvrmrmt) exit 
!%--------------------------------------------------------------
!% Allocate solution vector
!%--------------------------------------------------------------
 call check_pointer_ (deltacunk,nunk,.true.)
 call check_pointer_ (icolt,nunk,.true.)
 deltacunk=residual
!%--------------------------------------------------------------
! 2.2) Solve the lineal system
!%--------------------------------------------------------------
 if (isupdatejac) then
 
  call cpu_time (itime1)
 
  call f07bdf(nunk,nunk,nbandjac,nbandjac,jacobian,3*nbandjac+1, &
              icolt,info)
  
  if (info/=0) then 
    msg='Problems in LU factorization' 
    print *,msg
    goto 20 
  end if 
 
  call f07bef('n',nunk,nbandjac,nbandjac,1,jacobian,3*nbandjac+1, &
              icolt,deltacunk,nunk,info)

  if (info/=0) then 
    msg='Problems in LU factorization' 
    print *,msg
    goto 20 
  end if 

  call cpu_time (itime2)
  dtsolve=itime2-itime1
 
 end if
!%--------------------------------------------------------------
! 2.2) Check isconvergence and update the solution
!%--------------------------------------------------------------
  call cpu_time (itime1)
!%--------------------------------------------------------------
!% Determine the maximum error in residual
!% If the maximum error is > than maxerror then there was 
!% convergence in DSA
!%--------------------------------------------------------------
  mxerrorres=maxval(dabs(residual(1:nunk)))
  upmxerrorres=(mxerrorres>mxerror)
!%--------------------------------------------------------------
!% Update the solution 
!% WARNING!!!: The solution is expressed as
!% ln delta c1
!% We must change to delta c1
!%--------------------------------------------------------------
  call update_and_check_ &
  (this, &
   isconvergence, &
   divergence, &
   isanomalous, &
   upmxitergam, &
   mxerrorunk, &
   mxerrorres, &
   namemxunk, &
   namemxres, &
   inchmxunk, &
   inchmxres, &
   deltacunk, &
   nunk, &
   residual(1:nunk), &
   tolunk, &
   tolres, &
   dtime, &
   factor, &
   iserror)
 
 
   call cpu_time (itime2)
   dtupdate=itime2-itime1
!%--------------------------------------------------------------
!% Write information about time spent in diferent tasks
!%--------------------------------------------------------------
   if (this%iswriteinfo.and.haveiouwinfo) then
    write (iouwinfo,2) iter,isconvergence,mxerrorunk,inchmxunk, &
                        namemxunk,mxerrorres, inchmxres, &
                        namemxres,dtindep,dtjacres,dtsolve, &
                        dtupdate
   end if
!%-------------------------------------------------------------
   if (divergence) idiverg=idiverg+1
   upmxdiverg=(idiverg>mxdiverg)
   upmxiter=(iter>mxitertr)
   

end do
!%-------------------------------------------------------------
! 3) End fo newton raphson iterative process 
!%-------------------------------------------------------------
if (iserror) goto 20
!%-------------------------------------------------------------
! If there was convergence then
! 4) Compute equilibrium minerals concentrations from transport 
! equations
!%-------------------------------------------------------------
 if (isconvergence) then
   call cpu_time (itime1)
   call compute_eq_min_from_trp_ &
   (this, &
    iscompzchanged, &
    nnch, &
    nsp, &
    nband, &
    mxconn, &
    ngas, &
    theta, &
    dtime, &
    iconn, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    istrpgas, &
    tindep, &
	issolvermrmt, &
	imethod, &
	mxnbox, &
	iconnbox, &
	alpha, &
	phi, &
	msg, &
    iserror, &
	iouwinfo=iouwinfo)

    if (iserror) goto 20
    call cpu_time (itime2)
    dteqmin=itime2-itime1
!%-----------------------------------------------------------
    if (iscompzchanged) then
        ichange=ichange+1
    else
!%--------------------------------------------------------------
! Update the reactive surface of kinetic minerals 
! Update the txoh 
!%--------------------------------------------------------------
        do inch=1,this%numnch
           call update_min_area_ (this%pnodalchemk1(inch)%ptr,iserror) 
           if (iserror) goto 20
		   call update_txoh_ (this%pnodalchemk1(inch)%ptr,iserror) 
           if (iserror) goto 20 
        end do
        exit
   end if
   upmxchcompz=(ichange>mxchcompz)
 end if
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
 if (upmxdiverg) then
   if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) 'OVER MAXIMUM NUMBER OF ALLOWED DIVERGENCE'
   end if
   exit
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
 else if (upmxiter) then
   if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) 'OVER MAXIMUM NUMBER OF ALLOWED ITERATIONS'
   end if
   exit
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
 else if (upmxchcompz) then
   isconvergence=.false. 
   if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) 'OVER MAXIMUM NUMBER OF ALLOWED COMPONENTS CHANGED'
   end if
   exit
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
 else if (isanomalous) then
   if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) 'ANOMALOUS CONCENTRATIONS DETECTED'
   end if
   exit
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
 else if (upmxitergam) then
   if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) 'OVER MAXIMUM NUMBER OF ALLOWED GAMMA ITERATIONS'
   end if
   exit
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
 else if (upmxerrorres) then
   if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) 'OVER MAXIMUM ERROR IN RESIDUAL'
   end if
   exit
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
 else if (issolvermrmt.and..not.isconvrmrmt) then
   if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) 'NOT CONVERGENCE IN REACTIVE MULTI-RATE MASS TRANSFER CALCULATIONS'
   end if
   exit
 end if
 
end do
!%--------------------------------------------------------------
!% Write time spent computing mineral equilibrium
!%--------------------------------------------------------------
if (this%iswriteinfo.and.haveiouwinfo) then
  write (iouwinfo,*) '------------------------------------'// &
                     '------------------------------------'
  write (iouwinfo,4) 'dt eq.min:',dteqmin
  write (iouwinfo,*) '------------------------------------'// &
                     '------------------------------------'
  write (iouwinfo,*) '------------------------------------'// &
                     '------------------------------------'
end if 
!%--------------------------------------------------------------
! (Optionl task): Compute mass of components that input/output 
! of the system. Only if there was convergence 
!%--------------------------------------------------------------
if (havemolinout.and.isconvergence) then
          call compute_in_out_mol_ &
            (this, &
             molinout, &
             ndim, &
             nnch, &
             dtime, &
             caudalboundw, &
             caudalrechw, &
             itypeboundw, &
             isetboundw, &
             isetrechw, &
             gammaboundw, &
             iserror)
		  if (iserror) goto 20 
end if
!%--------------------------------------------------------------
! (Optional task): Compute mass of components preduced/consumed
! due kinetic reactions 
! Only if there was convergence 
!%--------------------------------------------------------------
if (havemolkin.and.isconvergence) then   
         call compute_kin_mol_ &
           (this, &
            molkin, &
            ndim, &
            dtime, &
            iserror)
         if (iserror) goto 20 
end if
!%--------------------------------------------------------------
20 continue 
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (tindep,1,1,.false.)
call check_pointer_ (jacobian,1,1,.false.)
call check_pointer_ (residual,1,.false.)
call check_pointer_ (deltacunk,1,.false.)
call check_pointer_ (icolt,1,.false.)         
if (iserror) goto 10
!%--------------------------------------------------------------
return
 
1 format(12a10)
2 format(i10,l5,e15.3,i9,a10,e10.3,i10,a11,4e10.2)
3 format(a10,i5)
4 format(a10,e10.2)
 
10 continue 
print *,'******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: solve_react_trp_step_'
print *, msg
print *,'******************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine solve_react_trp_step_sia_cheproo &
   (this, &
    iconn, &
    mxconn, &
    nnch, &
    typesto, &
    mxitertr, &
    mxdiverg, &
    theta, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    nband, &
    ngas, &
    dtime, &
    istrpgas, &
    caudalboundw, &
    caudalrechw, &
    caudalrechg, &
    caudalboundg, &
    itypeboundw, &
    itypeboundg, &
    gammaboundw, &
    isetboundw, &
    isetrechw, &
    isconvergence, &
	iter, &
    factor, &
    tolunk, &
    isopensystem, &
    time, &
	iswcompbal, &
	mxerror, &
    iserror, &
    molinout, &
	molkin, &
	iouwinfo)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Solve one time step of reactive transport (for SIA or SNIA)
!
!   $Arguments:
!
 
type (t_cheproo), intent(inout)      :: this  ! Type cheproo variable. 

integer, intent(in)                  :: mxconn

integer, intent(in)                  :: nnch

integer, intent(in)                  :: mxitertr

integer, intent(in)                  :: mxdiverg 

integer, intent(in)                  :: nband ! Half band 

integer, intent(in)                  :: ngas ! Number of gas phases 

integer, intent(in)                  :: typesto ! Type of storage 

real*8, intent(in)                   :: theta ! Temporal weight 

real*8, intent(in)                   :: dtime ! Time increment [s]

real*8, intent(in)                   :: factor

real*8, intent(in)                   :: tolunk

real*8, intent(in)                   :: time

real*8, intent(in)                   :: efloww_ij(nnch,2*nband+1)

real*8, intent(in)                   :: eflowg1_ij(nnch,2*nband+1)

real*8, intent(in)                   :: eflowg2_ij(nnch,2*nband+1)

real*8, intent(in)                   :: caudalboundw(nnch)

real*8, intent(in)                   :: caudalrechw(nnch)

real*8, intent(in)                   :: caudalboundg(nnch)

real*8, intent(in)                   :: caudalrechg(nnch)

real*8, intent(in)                   :: gammaboundw(nnch)

integer, intent(in)                  :: iconn(mxconn,nnch)

integer, intent(in)                  :: isetboundw(nnch)

integer, intent(in)                  :: isetrechw(nnch)

integer, intent(in)                  :: itypeboundw(nnch)

integer, intent(in)                  :: itypeboundg(nnch)

logical, intent(in)                  :: istrpgas

logical, intent(in)                  :: isopensystem

logical, intent(out)                 :: isconvergence

integer, intent(out)                 :: iter 

logical, intent(out)                 :: iserror

logical, intent(in)                  :: iswcompbal

real*8, intent(in)                   :: mxerror ! Maximum error 

real*8, pointer, optional            :: molinout(:) 

real*8, pointer, optional            :: molkin(:) 

integer, intent(in), optional        :: iouwinfo 
 
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
 msg, &
 namemxerrsia
logical                              :: &
 upmxiter, &
 isconvchem, &
 havemolinout, &
 havemolkin, &
 haveiouwinfo, &
 isiter 
real*8                               :: &
 dd, &
 mxerrorsia
integer                              :: &
 ithcomp, &
 ncomp, &
 info, &
 nsp, &
 nmobph, &
 nbandt, &
 ndim, &
 imxnch, &
 i, &
 nlock
integer, pointer                     :: &
 icolt(:) => null ()
type(t_vector), pointer              :: &
 utindep(:)  => null ()
real*8, pointer                      :: &
 rsia(:,:) => null (), &
 utra(:,:) => null (), &
 trpmatrix(:,:) => null (), &
 b(:) => null (), &
 tindep(:,:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Check optional arguments 
!%--------------------------------------------------------------
havemolinout=present(molinout)
havemolkin=present(molkin)
haveiouwinfo=present(iouwinfo)
!%--------------------------------------------------------------
!% Check the number of nodal chemistry objects 
!%--------------------------------------------------------------
if (nnch/=this%numnch) then
  msg='Error in number of nodal chemistry objects'
  goto 10
end if
!%--------------------------------------------------------------
!% Check the number of chemical system objects 
!%--------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system'
 goto 10
end if
!%--------------------------------------------------------------
!% Check the maximum error 
!%--------------------------------------------------------------
if (mxerror<=0.0d0) then 
 msg='Error in maximum error, mxerror='
 call add_ (msg,mxerror) 
 goto 10
end if 
!%--------------------------------------------------------------
!% Set if the water mass balance must be evaluated
!%--------------------------------------------------------------
call set_iswcompbal_ (this%pchemsys(1)%ptr,iswcompbal,iserror)
if (iserror) goto 10 
!%--------------------------------------------------------------
nmobph=ngas+1
nbandt=3*nband+1
!%--------------------------------------------------------------
!% 1) Build independent term for SIA (or SNIA)
!%--------------------------------------------------------------
call build_indep_term_sia_ &
   (this, &
    tindep, &
    nnch, &
    nsp, &
    theta, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    nband, &
    ngas, &
    dtime, &
    istrpgas, &
    caudalboundw, &
    caudalrechw, &
    caudalrechg, &
    caudalboundg, &
	gammaboundw, &
    itypeboundw, &
    itypeboundg, &
    isetboundw, &
    isetrechw, &
    iserror)
 if (iserror) then
   msg='Error when calling build_indep_term_sia_'
   goto 20
 end if
!%--------------------------------------------------------------
!% Write information about convergence (only for SIA)
!%--------------------------------------------------------------
if (this%iswriteinfo.and.haveiouwinfo.and.this%itypesolver==siasolver) then
      write (iouwinfo,*) '------------------------------------------------------------------------'
      write (iouwinfo,*) 'Service: solve_react_trp_step_'
      write (iouwinfo,*) '------------------------------------------------------------------------'
      write (iouwinfo,*) '      SIA information   '
      write (iouwinfo,*) 'Tolerance in unknowns:',tolunk
      write (iouwinfo,*) 'Mx. number of divergence',mxdiverg
      write (iouwinfo,*) 'Mx. number of trp. iter.',mxitertr
      write (iouwinfo,*) 'Time:',time
      write (iouwinfo,*) 'Time increment:',dtime
      write (iouwinfo,*) '------------------------------------------------------------------------'
      write (iouwinfo,*) '------------------------------------------------------------------------'
      write (iouwinfo,1) 'iter','converge','mx.errsia','node','comp'
      write (iouwinfo,*) '------------------------------------------------------------------------'
end if
!%--------------------------------------------------------------
!% Prepare preview information 
!%--------------------------------------------------------------
call make_chem_info_sia_ &
   (this, &
    tindep, &
    utindep, &
    nsp, &
    nnch, &
    ncomp, &
    msg, &
    iserror)
if (iserror) goto 20
!%--------------------------------------------------------------
!% Allocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (rsia,nnch,ncomp,.true.)
call check_pointer_ (utra,nnch,ncomp,.true.)
call check_pointer_ (icolt,nnch,.true.)
!%--------------------------------------------------------------
!%--------------------------------------------------------------
iter=0
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Start transport iterative process 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
do
!%--------------------------------------------------------------
!%--------------------------------------------------------------
 iter=iter+1
!%--------------------------------------------------------------
!% Loop for components 
!%-------------------------------------------------------------- 
 do ithcomp=1,ncomp
 
 
!%--------------------------------------------------------------
!% 2) Build the transport equation for sia for ith component 
!%--------------------------------------------------------------
   call build_trp_eq_sia_ &
   (this, &
    trpmatrix, &
    b, &
    nnch, &
    nmobph, &
    mxconn, &
    nband, &
    typesto, &
    theta, &
    dtime, &
    iconn, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    istrpgas, &
    utindep, &
    rsia, &
    ncomp, &
    ithcomp, &
    iserror)
 
  if (iserror) goto 20
 
!%--------------------------------------------------------------
!% 3) Solve the lineal system for each component 
!%--------------------------------------------------------------
  utra(:,ithcomp)=b
 
  call f07bdf(nnch,nnch,nband,nband,trpmatrix,nbandt,icolt,info)
 
  call f07bef('n',nnch,nband,nband,1,trpmatrix,nbandt,icolt,utra(:,ithcomp),nnch,info)
 
 end do
!%--------------------------------------------------------------
!% End loop for components 
!%--------------------------------------------------------------
 
!%--------------------------------------------------------------
!% 4) Update nodal chemistries and rsia and check convergence
!%--------------------------------------------------------------
 call update_and_check_ &
   (this, &
    rsia, &
    utra, &
    nnch, &
    ncomp, &
    imxnch, &
    mxerrorsia, &
    namemxerrsia, &
    dtime, &
    theta, &
    tolunk, &
	factor, &
    isconvergence, &
    isconvchem, &
    iter, &
    isopensystem, &
    iserror)
 
 if (iserror) goto 20
!%--------------------------------------------------------------
 upmxiter=(iter>mxitertr)
!%--------------------------------------------------------------
!% Write convergence information 
!%--------------------------------------------------------------
 if (this%iswriteinfo.and.haveiouwinfo) then
    write (iouwinfo,2) iter,isconvergence,mxerrorsia,imxnch,namemxerrsia
 end if
!%--------------------------------------------------------------
!% If not there was convergence in chemistry step, then return 
!%--------------------------------------------------------------
 if (.not.isconvchem) then
  if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) 'CONVERGENCE PROBLEMS IN CHEMISTRY STEP IN ITERATION',iter
  end if
  exit
 end if
!%--------------------------------------------------------------
!% If the solver is SNIA, then isconvsia=true
!%--------------------------------------------------------------
 if (this%itypesolver==sniasolver) isconvergence=.true.  
!%--------------------------------------------------------------
 if (isconvergence.and.iter>1) then
  if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) 'CONVERGENCE ARRIVED'
  end if
  exit
 end if 
!%--------------------------------------------------------------
 if (upmxiter) then
  if (this%iswriteinfo.and.haveiouwinfo) then
     write (iouwinfo,*) '*********************************'
     write (iouwinfo,*) 'Over mx. number or SIA iterations'
     write (iouwinfo,*) '*********************************'
  end if
  exit
 end if
!%--------------------------------------------------------------
end do
!%--------------------------------------------------------------
!% End transport iterative process
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% (Optional task): Compute mass of components that 
!% input/output of the system 
!%--------------------------------------------------------------
if (havemolinout.and.isconvergence) then
 
  call compute_in_out_mol_ &
        (this, &
         molinout, &
         ndim , &
         nnch, &
         dtime, &
         caudalboundw, &
         caudalrechw, &
         itypeboundw, &
         isetboundw, &
         isetrechw, &
         gammaboundw, &
         iserror)
  if (iserror) goto 20 
 
end if
!%--------------------------------------------------------------
!% (Optional task): Compute mass of components preduced/consumed
!% due kinetic reactions 
!%--------------------------------------------------------------
if (havemolkin.and.isconvergence) then   
   call compute_kin_mol_(this,molkin,ndim,dtime,iserror)
   if (iserror) goto 20 
end if 
!%--------------------------------------------------------------
if (this%iswriteinfo.and.haveiouwinfo) then
 write (iouwinfo,*) '------------------------------------------------------------------------'
end if
!%--------------------------------------------------------------
20 continue 
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (rsia,1,1,.false.)
call check_pointer_ (utra,1,1,.false.)
call check_pointer_ (trpmatrix,1,1,.false.)
call check_pointer_ (b,1,.false.)
call check_pointer_ (tindep,1,1,.false.)
call check_pointer_ (icolt,1,.false.)
!%--------------------------------------------------------------
do i=1,nnch
 call check_pointer_ (utindep(i)%vector,1,.false.)
end do
!%--------------------------------------------------------------
deallocate (utindep)
!%--------------------------------------------------------------
if (iserror) goto 10
!%--------------------------------------------------------------
return
1 format(5a10)
2 format(i10,l5,e15.3,i9,a10)
 
10 continue 
print *,'******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: solve_react_trp_step_'
print *, msg
print *,'******************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine make_preview_info_dsa_cheproo &
   (this, &
    czlocconn, &
    indexconn, &
    nunkconn, &
    ncznod, &
    posunk, &
    nunk, &
    mxnunkz, &
    nnch, &
    iconn, &
    mxconn, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Make the preview zonal information per each node (DSA)
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                    :: this

integer, intent(in)                             :: nnch                        ! Number of nodes

integer, intent(in)                             :: mxconn                      ! Maximun number of connected nodes

integer, intent(out)                            :: nunk                        ! Total number of unknowns

integer, intent(out)                            :: mxnunkz                     ! Maximun number of unknowns of all nodes

integer, intent(in)                             :: iconn(mxconn,nnch)

integer, intent(out)                            :: posunk(2,nnch)

integer, intent(out)                            :: czlocconn(mxconn,nnch)

integer, intent(out)                            :: indexconn(mxconn,nnch)

integer, intent(out)                            :: ncznod(nnch)

integer, intent(out)                            :: nunkconn(mxconn,nnch)

logical, intent(out)                            :: iserror

character(len=*), intent(out)                   :: msg 
 
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
 ithnch, &
 jthnch, &
 ipos1, &
 ipos2, &
 icz, &
 index, &
 indexithnch, &
 ic, &
 nunkith, &
 ncon, &
 nczjthnch, &
 nunkithnch
logical                            :: &
 be
integer, pointer                   :: &
 indexnod(:) => null (), &
 nunknod(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
call check_pointer_ (indexnod,nnch,.true.)
call check_pointer_ (nunknod,nnch,.true.)
!%----------------------------------------------------------
ipos1=1
ipos2=0
nunk=0
posunk=0
indexconn=0
czlocconn=0
ncznod=0
!%----------------------------------------------------------
do jthnch=1,nnch
 
   call get_chem_info_ (this%pnodalchemk1(jthnch)%ptr,iserror,hashcompz=index)
 
   call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numbase=nunkith,hashcompz=index)
 
   nunknod(jthnch)=nunkith
   indexnod(jthnch)=index
 
   nunk=nunk+nunkith
   ipos2=ipos2+nunkith
   posunk(1,jthnch)=ipos1
   posunk(2,jthnch)=ipos2
   ipos1=ipos1+nunkith
 
 
end do
!%-------------------------------------------------------
do jthnch=1,nnch
 ncon=iconn(1,jthnch)
 nunkconn(1,jthnch)=nunknod(jthnch)
 indexconn(1,jthnch)=indexnod(jthnch)
 czlocconn(1,jthnch)=ncon
 ncznod(jthnch)=1
 nczjthnch=1
  do ic=2,ncon
   be=.false.
   ithnch=iconn(ic,jthnch)
   indexithnch=indexnod(ithnch)
   nunkithnch=nunknod(ithnch)
     do icz=1,nczjthnch
      index=indexconn(icz,jthnch)
       if(index==indexithnch) then
          be=.true.
        exit
         end if
     end do
         if (be) then
         czlocconn(ic,jthnch)=icz
         be=.false.
       else
         nczjthnch=nczjthnch+1
         ncznod(jthnch)=nczjthnch
         czlocconn(ic,jthnch)=nczjthnch
         indexconn(nczjthnch,jthnch)=indexithnch
         nunkconn(nczjthnch,jthnch)=nunkithnch
       end if
 
  end do
 
 
end do
!%----------------------------------------------------------
mxnunkz=maxval(nunknod)
!%----------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------
call check_pointer_ (indexnod,nnch,.true.)
call check_pointer_ (nunknod,nnch,.true.)
!%----------------------------------------------------------
return
 
10 iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine make_zonal_chem_info_ithnch_cheproo &
   (this, &
    ncz, &
    czindx, &
    nc, &
    tindep, &
    utindep, &
    u, &
    umob, &
    uads, &
    usktrk, &
    du, &
    dumob, &
    duads, &
    dusktrk, &
    nmobph, &
    ithnch, &
    numsp, &
    isderivates, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the compositional information of the ith nodal
! chemistry object. 
!
!   $Arguments:
!
 
type(t_cheproo), intent(in)           :: this            ! Type CHEPROO object

integer, intent(in)                   :: ithnch          ! Index of the nodal chemistry object 

integer, intent(in)                   :: ncz             ! Total number of components sets associated to ithnch

integer, intent(out)                  :: nmobph          ! Total number of mobile phases 

integer, intent(in)                   :: numsp           ! Number of species 

integer, intent(in)                   :: czindx(ncz)     ! Components set index

integer, intent(in)                   :: nc(ncz)         ! Number of components per set of components

real*8, intent(in)                    :: tindep(numsp)

character(len=*), intent(out)         :: msg

type(t_vector) , pointer              :: utindep(:)

type(t_vector) , pointer              :: u(:)

type(t_vector) , pointer              :: usktrk(:)

type(t_vector) , pointer              :: umob(:,:)

type(t_vector) , pointer              :: uads(:)

type(t_array), pointer                :: du(:)

type(t_array), pointer                :: duads(:)

type(t_array), pointer                :: dumob(:,:)

type(t_array), pointer                :: dusktrk(:)

logical, intent(out)                  :: iserror

logical, intent(in)                   :: isderivates
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
real*8, pointer                       :: &
 c(:) => null (), &
 cpri(:) => null (), &
 cads(:) => null (), &
 cmob(:,:) => null (), &
 dc(:,:) => null (), &
 dcmob(:,:) => null (), &
 dcads(:,:) => null (), &
 sktrk(:) => null (), &
 dsktrk(:,:) => null (), &
 array(:,:) => null (), &
 vector(:) => null ()
integer                              :: &
 ipos1, &
 ipos2, &
 i, &
 j, &
 k, &
 l, &
 ndim1, &
 ndim2, &
 ncomp, &
 numsploc 
type(t_nodalchemistry), pointer     :: &
 pnch => null () 
!-------------------------------------------------------------------------
!
!   $code
!
 


!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Pointer to nodal chemistry object 
!%------------------------------------------------------------
pnch => this%pnodalchemk1(ithnch)%ptr
!%------------------------------------------------------------
if (associated(utindep)) then
 ndim1=size(utindep)
 do i=1,ndim1
  deallocate (utindep(i)%vector)
 end do
 deallocate (utindep)
end if
!%------------------------------------------------------------
if (associated(u)) then
 ndim1=size(u)
 do i=1,ndim1
  deallocate (u(i)%vector)
 end do
 deallocate (u)
end if
!%------------------------------------------------------------
if (associated(usktrk)) then
 ndim1=size(usktrk)
 do i=1,ndim1
  deallocate (usktrk(i)%vector)
 end do
 deallocate (usktrk)
end if
!%------------------------------------------------------------
if (associated(uads)) then
 ndim1=size(uads)
 do i=1,ndim1
  deallocate (uads(i)%vector)
 end do
 deallocate (uads)
end if
!%------------------------------------------------------------
if (associated(umob)) then
 ndim1=size(umob,1)
 ndim2=size(umob,2)
 do i=1,ndim1
  do j=1,ndim2
   deallocate (umob(i,j)%vector)
  end do
 end do
 deallocate (umob)
end if
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
if (isderivates) then
 
 if (associated(du)) then
  ndim1=size(du)
  do i=1,ndim1
   deallocate (du(i)%array)
  end do
  deallocate (du)
 end if
!%------------------------------------------------------------
 if (associated(duads)) then
  ndim1=size(duads)
  do i=1,ndim1
   deallocate (duads(i)%array)
  end do
  deallocate (duads)
 end if
!%------------------------------------------------------------
 if (associated(dusktrk)) then
  ndim1=size(dusktrk)
  do i=1,ndim1
   deallocate (dusktrk(i)%array)
  end do
  deallocate (dusktrk)
 end if
!%------------------------------------------------------------
 if (associated(dumob)) then
  ndim1=size(dumob,1)
  ndim2=size(dumob,2)
  do i=1,ndim1
   do j=1,ndim2
    deallocate (dumob(i,j)%array)
   end do
  end do
  deallocate (dumob)
 end if
 
end if
!%------------------------------------------------------------
!% Get chemical information 
!%------------------------------------------------------------
if (isderivates) then
 call get_chem_info_ &
 (pnch, &
  iserror, &
  c=c, &
  cmob=cmob, &
  cads=cads, &
  cpri=cpri, &
  sktrk=sktrk, &
  dc=dc, &
  dcmob=dcmob, &
  dcads=dcads, &
  dsktrk=dsktrk, &
  numsp=numsploc, &
  nummobph=nmobph)
  ndim1=size(dc,2)
else
 call get_chem_info_ &
 (pnch, &
  iserror, &
  c=c, &
  cmob=cmob, &
  cads=cads, &
  sktrk=sktrk, &
  numsp=numsploc, &
  nummobph=nmobph)
end if
if (iserror) then
  msg='Error when calling get_chem_info_'
  goto 20
end if
!%------------------------------------------------------------
!% Check the number of species 
!%------------------------------------------------------------
if (numsp/=numsploc) then
 iserror=.true. 
 msg='Error in number of species'
 goto 20
end if 
!%------------------------------------------------------------
!% Allocate pointers 
!%------------------------------------------------------------
allocate (utindep(ncz))
allocate (u(ncz))
allocate (umob(ncz,nmobph))
allocate (uads(ncz))
allocate (usktrk(ncz))
if (isderivates) then
 allocate (dumob(ncz,nmobph))
 allocate (du(ncz))
 allocate (duads(ncz))
 allocate (dusktrk(ncz))
end if
!%------------------------------------------------------------
do i=1,ncz
 
 ncomp=nc(i)
 nullify(utindep(i)%vector)
 nullify(u(i)%vector)
 nullify(uads(i)%vector)
 nullify(usktrk(i)%vector)
 call check_pointer_ (utindep(i)%vector,ncomp,.true.)
 call check_pointer_ (u(i)%vector,ncomp,.true.)
 call check_pointer_ (uads(i)%vector,ncomp,.true.)
 call check_pointer_ (usktrk(i)%vector,ncomp,.true.)
 
 utindep(i)%nrow=ncomp
 u(i)%nrow=ncomp
 uads(i)%nrow=ncomp
 usktrk(i)%nrow=ncomp
!%------------------------------------------------------------
!% utindep
!%------------------------------------------------------------
 call make_lin_trf_ (this%pchemsys(1)%ptr,vector,tindep,czindx(i),iserror)
 
 if (iserror) then
  msg='Error when calling make_lin_trf_'
  goto 20
 end if
 
 utindep(i)%vector=vector
!%------------------------------------------------------------
!% u
!%------------------------------------------------------------
 call make_lin_trf_ (this%pchemsys(1)%ptr,vector,c,czindx(i),iserror)
 
 if (iserror) then
  msg='Error when calling make_lin_trf_'
  goto 20
 end if
 
 u(i)%vector=vector
!%------------------------------------------------------------
!% usktrk
!%------------------------------------------------------------
 call make_lin_trf_ (this%pchemsys(1)%ptr,vector,sktrk,czindx(i),iserror)
 
 if (iserror) then
  msg='Error when calling make_lin_trf_'
  goto 20
 end if
 
 usktrk(i)%vector=vector
!%------------------------------------------------------------ 
!% uads
!%------------------------------------------------------------
 call make_lin_trf_ (this%pchemsys(1)%ptr,vector,cads,czindx(i),iserror)
 
 if (iserror) then
  msg='Error when calling make_lin_trf_'
  goto 20
 end if
 
 
 uads(i)%vector=vector
!%------------------------------------------------------------ 
!% umob
!%------------------------------------------------------------
  do j=1,nmobph
   nullify(umob(i,j)%vector)
   call check_pointer_ (umob(i,j)%vector,ncomp,.true.)
   umob(i,j)%nrow=ncomp
   call make_lin_trf_ (this%pchemsys(1)%ptr,vector,cmob(:,j),czindx(i),iserror)
   if (iserror) then
      msg='Error when calling make_lin_trf_'
      goto 20
   end if
   umob(i,j)%vector=vector
  end do
!%------------------------------------------------------------ 
!% If derivatives must be computed 
!%------------------------------------------------------------
  if (isderivates) then
 
   nullify(du(i)%array)
   nullify(duads(i)%array)
   nullify(dusktrk(i)%array)
 
   call check_pointer_ (du(i)%array,ncomp,ndim1,.true.)
   call check_pointer_ (duads(i)%array,ncomp,ndim1,.true.)
   call check_pointer_ (dusktrk(i)%array,ncomp,ndim1,.true.)
 
   du(i)%nrow=ncomp
   du(i)%ncol=ndim1
 
   duads(i)%nrow=ncomp
   duads(i)%ncol=ndim1
 
   dusktrk(i)%nrow=ncomp
   dusktrk(i)%ncol=ndim1
!%------------------------------------------------------------ 
!% du 
!%------------------------------------------------------------
   call make_lin_trf_(this%pchemsys(1)%ptr,array,dc,czindx(i),iserror)
 
   if (iserror) then
     msg='Error when calling make_lin_trf_'
     goto 20
   end if
!%------------------------------------------------------------ 
!% Change 
!% 
!%    du     du 
!%   ---- = ---- c1
!% dln c1    dc1
!%------------------------------------------------------------
 !  do k=1,ncomp  
  ! 	 do l=1,ndim1
 !	   array(k,l)=array(k,l)*cpri(l)
 !	 end do 
 !  end do 
   du(i)%array=array
!%------------------------------------------------------------ 
!% duads 
!%------------------------------------------------------------
    call make_lin_trf_(this%pchemsys(1)%ptr,array,dcads,czindx(i),iserror)
 
    if (iserror) then
      msg='Error when calling make_lin_trf_'
      goto 20
    end if
!%------------------------------------------------------------ 
!% Change 
!% 
!%    du     du 
!%   ---- = ---- c1
!% dln c1    dc1
!%------------------------------------------------------------
 !   do k=1,ncomp  
 !  	 do l=1,ndim1
 !	   array(k,l)=array(k,l)*cpri(l)
 !	 end do 
  !  end do
    
	duads(i)%array=array
!%------------------------------------------------------------ 
!% dusktrk
!%------------------------------------------------------------
    call make_lin_trf_ (this%pchemsys(1)%ptr,array,dsktrk,czindx(i),iserror)
 
    if (iserror) then
       msg='Error when calling make_lin_trf_'
       goto 20
    end if
!%------------------------------------------------------------ 
!% Change 
!% 
!%    du     du 
!%   ---- = ---- c1
!% dln c1    dc1
!%------------------------------------------------------------
 !  do k=1,ncomp  
  ! 	 do l=1,ndim1
!	   array(k,l)=array(k,l)*cpri(l)
!	 end do 
  ! end do
   
   dusktrk(i)%array=array
!%------------------------------------------------------------ 
!% dumob
!%------------------------------------------------------------
    ipos1=0
    ipos2=0
      do j=1,nmobph
        nullify(dumob(i,j)%array)
        call check_pointer_ (dumob(i,j)%array,ncomp,ndim1,.true.)
        dumob(i,j)%nrow=ncomp
        dumob(i,j)%ncol=ndim1
        ipos1=ipos2+1
        ipos2=ipos2+numsp
!%------------------------------------------------------------
!% Compute Uc 
!%------------------------------------------------------------
        call make_lin_trf_ (this%pchemsys(1)%ptr,array,dcmob(ipos1:ipos2,:),czindx(i),iserror)
 
         if (iserror) then
         msg='Error when calling make_lin_trf_'
         goto 20
       end if
!%------------------------------------------------------------ 
!% Change 
!% 
!%    du     du 
!%   ---- = ---- c1
!% dln c1    dc1
!%------------------------------------------------------------
     !  do k=1,ncomp  
   	 !    do l=1,ndim1
	 !      array(k,l)=array(k,l)*cpri(l)
	 !    end do 
     !  end do

       dumob(i,j)%array=array
!%--------------
 
    end do
 
 
 
  end if
   
end do
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (c,1,.false.)
call check_pointer_ (cpri,1,.false.)
call check_pointer_ (cmob,1,1,.false.)
call check_pointer_ (cads,1,.false.)
call check_pointer_ (sktrk,1,.false.)
call check_pointer_ (array,1,1,.false.)
call check_pointer_ (vector,1,.false.)
if (isderivates) then
 call check_pointer_ (dcmob,1,1,.false.)
 call check_pointer_ (dcads,1,1,.false.)
 call check_pointer_ (dsktrk,1,1,.false.)
 call check_pointer_ (dc,1,1,.false.)
end if
pnch => null ()
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
subroutine make_chem_info_sia_cheproo &
   (this, &
    tindep, &
    utindep, &
    numsp, &
    nnch, &
    ncomp, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Make the chemical information for all nodal chemistries 
!
!   $Arguments:
!
 
type(t_cheproo), intent(in)           :: this

integer, intent(in)                   :: nnch

integer, intent(in)                   :: numsp

integer, intent(out)                  :: ncomp

real*8, intent(in)                    :: tindep(numsp,nnch)

character(len=*), intent(out)         :: msg

type(t_vector), pointer               :: utindep(:)

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
real*8, pointer                       :: &
 vector(:) => null ()
integer                              :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%------------------------------------------------------------
if (associated(utindep)) deallocate (utindep)
!%------------------------------------------------------------
allocate (utindep(nnch))
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numbase=ncomp)
!%------------------------------------------------------------
do i=1,nnch
 
 nullify(utindep(i)%vector)
 
 call check_pointer_ (utindep(i)%vector,ncomp,.true.)
 
 utindep(i)%nrow=ncomp
 
 call make_lin_trf_ (this%pchemsys(1)%ptr,vector,tindep(:,i),0,iserror)
 
 if (iserror) then
  msg='Error when calling make_lin_trf_'
  goto 10
 end if
 
 utindep(i)%vector=vector

end do
!%------------------------------------------------------------
call check_pointer_ (vector,1,.false.)
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
subroutine build_jacob_resid_dsa_ithnch_cheproo &
   (this, &
    jacobian, &
    residual, &
    nnch, &
    nmobph, &
    ncz, &
    nbandjac, &
    mxconn, &
    nband, &
    nunk, &
    ngas, &
    theta, &
    dtime, &
	itypeboundw, &
    iconn, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    posunk, &
    czlocconn, &
    istrpgas, &
    d, &
    u, &
    umob, &
    uads, &
    usktrk, &
    du, &
    dumob, &
    duads, &
    dusktrk, &
    jthnch, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build the jacobian and residual for DSA for ith transport equation. 
!
!   $Arguments:
!
 
type(t_cheproo), intent(in)           :: this

integer, intent(in)                   :: jthnch ! Index of nodal chemistry 

integer, intent(in)                   :: ncz  ! Total number of components sets associated to ithnch

integer, intent(in)                   :: nmobph  ! Total number of mobile phases

integer, intent(in)                   :: nunk

integer, intent(in)                   :: nnch                       ! Number of nodal chemistries

integer, intent(in)                   :: ngas

integer, intent(in)                   :: mxconn                     ! Maximum number of connected nodes

integer, intent(in)                   :: nband

integer, intent(in)                   :: nbandjac 

real*8, intent(in)                    :: theta                       ! Temporal weight 

real*8, intent(in)                    :: dtime                       ! Time increment 

integer, intent(in)                   :: itypeboundw(nnch) 

real*8, intent(in)                    :: efloww_ij(nnch,2*nband+1)

real*8, intent(in)                    :: eflowg1_ij(nnch,2*nband+1)

real*8, intent(in)                    :: eflowg2_ij(nnch,2*nband+1)

real*8, pointer                       :: residual(:)                  ! Residual vector 

real*8, pointer                       :: jacobian(:,:)                ! Jacobian matrix

integer, intent(in)                   :: iconn(mxconn)                ! Connectivity between nodal chemistries 

integer, intent(in)                   :: posunk(2,nnch)

integer, intent(in)                   :: czlocconn(mxconn)

character(len=*), intent(out)         :: msg

type(t_vector), intent(in)            :: d(ncz)

type(t_vector), intent(in)            :: u(ncz)

type(t_vector), intent(in)            :: uads(ncz)

type(t_vector), intent(in)            :: usktrk(ncz)

type(t_vector), intent(in)            :: umob(ncz,nmobph)

type(t_array), intent(in)             :: du(ncz)

type(t_array), intent(in)             :: dumob(ncz,nmobph)

type(t_array), intent(in)             :: duads(ncz)

type(t_array), intent(in)             :: dusktrk(ncz)

logical, intent(in)                   :: istrpgas ! If .true. the gas is transported 

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
integer                               :: &
 i, &
 ipos1jthnch, &
 ipos2jthnch, &
 ipos1ithnch, &
 ipos2ithnch, &
 ijac1, &
 ijac2, &
 iczithnch, &
 iczjthnch, &
 ithnch, &
 ic, &
 ncon, &
 iband
real*8                                :: &
 alpha, &
 fii, &
 eij 
type(t_nodalchemistry), pointer       :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%------------------------------------------------------------
ipos1jthnch=posunk(1,jthnch)
ipos2jthnch=posunk(2,jthnch)
iczjthnch=1
ncon=iconn(1)
iband=nband+1
pnch => this%pnodalchemk1(jthnch)%ptr
call get_chem_info_ (pnch,iserror,omgwfree=fii)
if (iserror) then
 msg='Error when calling get_chem_info_ in nodal chemistry:'
 call add_ (msg,jthnch)
 goto 10
end if
eij=efloww_ij(jthnch,iband)
!%-------------------------------------------------------------
! Add in residual
! -u * tindep
!%-------------------------------------------------------------
alpha=-1.0d0
call add_ (residual(ipos1jthnch:ipos2jthnch),alpha,d(iczjthnch)%vector,d(iczjthnch)%nrow)
!%-------------------------------------------------------------
! Add in residual
! fii
! --- umob(k+1)
! dt
!%-------------------------------------------------------------
alpha=fii/dtime
do i=1,nmobph
 call add_ (residual(ipos1jthnch:ipos2jthnch),alpha,umob(iczjthnch,i)%vector,umob(iczjthnch,i)%nrow)
end do
!%-------------------------------------------------------------
! Add in residual
! fii
! --- uads(k+1)
! dt
!%-------------------------------------------------------------
call add_ (residual(ipos1jthnch:ipos2jthnch),alpha,uads(iczjthnch)%vector,uads(iczjthnch)%nrow)
!%-------------------------------------------------------------
! Add in jacobian
! fii
! --- dumob(k+1)
! dt
!%-------------------------------------------------------------
do i=1,nmobph
 call add_jacobian_dsa_ (jacobian,nunk,ipos1jthnch,ipos1jthnch,nbandjac,dumob(iczjthnch,i),alpha,msg,iserror)
 if (iserror) goto 10
end do
!%-------------------------------------------------------------
! Add in jacobian
! fii
! --- duads(k+1)
! dt
!%-------------------------------------------------------------
call add_jacobian_dsa_(jacobian,nunk,ipos1jthnch,ipos1jthnch,nbandjac,duads(iczjthnch),alpha,msg,iserror)
if (iserror) goto 10
!%-------------------------------------------------------------
! Add in residual,
! -fii*theta*u*skt*rk(k+1)
!%-------------------------------------------------------------
alpha=-fii*theta 
call add_ (residual(ipos1jthnch:ipos2jthnch),alpha,usktrk(iczjthnch)%vector,usktrk(iczjthnch)%nrow)
!%-------------------------------------------------------------
! Add in jacobian,
! -fii*theta*u*skt*drk(k+1)
!%-------------------------------------------------------------
call add_jacobian_dsa_(jacobian,nunk,ipos1jthnch,ipos1jthnch,nbandjac,dusktrk(iczjthnch),alpha,msg,iserror)
if (iserror) goto 10
!%-------------------------------------------------------------
! Add in residual,
! theta*Eij*umob(k+1)
!%-------------------------------------------------------------
alpha=theta*eij
call add_ (residual(ipos1jthnch:ipos2jthnch),alpha,umob(iczjthnch,1)%vector,umob(iczjthnch,1)%nrow)
!%-------------------------------------------------------------
! Add in jacobian,
! theta*Eij*dumob(k+1)
!%-------------------------------------------------------------
call add_jacobian_dsa_(jacobian,nunk,ipos1jthnch,ipos1jthnch,nbandjac,dumob(iczjthnch,1),alpha,msg,iserror)
if (iserror) goto 10
!%----------------------------------
! Add another mobile phase
!%----------------------------------
if (istrpgas) then
 do i=2,nmobph
  call add_ (residual(ipos1jthnch:ipos2jthnch),alpha,umob(iczjthnch,i)%vector,umob(iczjthnch,i)%nrow)
 end do
end if
!%-------------------------------------------------------------
! Now add connected nodes 
!%-------------------------------------------------------------
do ic=2,ncon
 ithnch=iconn(ic)
 iczithnch=czlocconn(ic)
 ipos1ithnch=posunk(1,ithnch)
 ipos2ithnch=posunk(2,ithnch)
 iband=nband+1-ithnch+jthnch
 eij=efloww_ij(ithnch,iband)
!%-------------------------------------------------------------
! Add in residual,
! theta*Eij*umob(k+1)
!%-------------------------------------------------------------
 alpha=theta*eij
 call add_ (residual(ipos1ithnch:ipos2ithnch),alpha,umob(iczithnch,1)%vector,umob(iczithnch,1)%nrow)
!%-------------------------------------------------------------
! Add in jacobian,
! theta*Eij*dumob(k+1)
!%-------------------------------------------------------------
 call add_jacobian_dsa_(jacobian,nunk,ipos1jthnch,ipos1ithnch,nbandjac,dumob(iczithnch,1),alpha,msg,iserror)
 if (iserror) goto 10
!%----------------------------------
! Add another mobile phase
!%----------------------------------
 if (istrpgas) then
  do i=2,nmobph
   call add_(residual(ipos1ithnch:ipos2ithnch),alpha,umob(iczithnch,i)%vector,umob(iczithnch,i)%nrow)
  end do
 end if
!%-------------------------------------------------------------
end do
!%-------------------------------------------------------------
pnch => null ()
!%-------------------------------------------------------------
return
 
10 continue 
print *,'***************************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: build_jacob_resid_dsa_ithnch_'
print *, msg
print *,'***************************************'
iserror=.true.
return
 
 end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine build_trp_eq_sia_ithnch_cheproo &
   (this, &
    trpmatrix, &
    b, &
    nnch, &
    nmobph, &
    mxconn, &
    nband, &
    theta, &
    dtime, &
    iconn, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    istrpgas, &
    iutindep, &
    irsia, &
    jthnch, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Build the transport equation for ith nodal chemistry 
!
!   $Arguments:
!
 
 
type(t_cheproo), intent(in)           :: this

integer, intent(in)                   :: jthnch

integer, intent(in)                   :: nmobph  ! Total number of mobile phases

integer, intent(in)                   :: nnch ! Number of nodes

integer, intent(in)                   :: mxconn ! Maximum number of connected nodes 

integer, intent(in)                   :: nband

real*8, intent(in)                    :: theta ! Temporal weight 

real*8, intent(in)                    :: dtime ! Time increment 

real*8, intent(in)                    :: iutindep

real*8, intent(in)                    :: irsia ! Reaction term of SIA 

real*8, intent(in)                    :: efloww_ij(nnch,2*nband+1)

real*8, intent(in)                    :: eflowg1_ij(nnch,2*nband+1)

real*8, intent(in)                    :: eflowg2_ij(nnch,2*nband+1)

integer, intent(in)                   :: iconn(mxconn) ! Connectivity between nodal chemistry objects 

character(len=*), intent(out)         :: msg

logical, intent(in)                   :: istrpgas ! If true, the gases are trasportated 

real*8, intent(inout)                 :: b(nnch)

real*8, intent(inout)                 :: trpmatrix(3*nband+1,nnch) ! Transport matrix 

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
integer                               :: &
 i, &
 ithnch, &
 ic, &
 ncon, &
 iband, &
 ipos
real*8                                :: &
 alpha, &
 fii, &
 eij 
type(t_nodalchemistry), pointer       :: &
 pnch => null ()
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false. 
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
pnch => this%pnodalchemk1(jthnch)%ptr
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
ncon=iconn(1)
iband=nband+1
ipos=2*nband+1
call get_chem_info_ (pnch,iserror,omgwfree=fii)
if (iserror) then
 msg='Error when calling get_chem_info_ in nodal chemistry:'
 call add_ (msg,jthnch)
 goto 10
end if
eij=efloww_ij(jthnch,iband)
!%-------------------------------------------------------------
! Add in independent term
! u * tindep
!%-------------------------------------------------------------
alpha=1.0d0
b(jthnch) = b(jthnch) + alpha*iutindep
!%-------------------------------------------------------------
! Add in independent term
! fii*rsia
!%-------------------------------------------------------------
alpha=fii
b(jthnch) = b(jthnch) + alpha*irsia
!%-------------------------------------------------------------
! Add in transport matrix
! fii
! --- umob(k+1)
! dt
!%-------------------------------------------------------------
alpha=fii/dtime
trpmatrix(ipos,jthnch) = trpmatrix(ipos,jthnch) + alpha
!%-------------------------------------------------------------
! Add in transport matrix
! theta*Eij*umob(k+1)
!%-------------------------------------------------------------
alpha=theta*eij
trpmatrix(ipos,jthnch) = trpmatrix(ipos,jthnch) + alpha
!%-------------------------------------------------------------
! Add in transport matrix information about connected nodal 
! chemistry objects 
!%-------------------------------------------------------------
do ic=2,ncon
 ithnch=iconn(ic)
 iband=nband+1-ithnch+jthnch
 ipos=2*nband+1+ithnch-jthnch
 eij=efloww_ij(ithnch,iband)
!%-------------------------------------------------------------
! Add in transport matrix
! theta*Eij*umob(k+1)
!%-------------------------------------------------------------
 alpha=theta*eij
 trpmatrix(ipos,jthnch) = trpmatrix(ipos,jthnch) + alpha
!%-------------------------------------------------------------
end do
!%-------------------------------------------------------------
!% Nullify pointers 
!%-------------------------------------------------------------
pnch => null ()
!%-------------------------------------------------------------
return
 
10 continue 
print *,'*************************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: build_trp_eq_sia_ithnch_'
print *, msg
print *,'*************************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_eq_min_from_trp_cheproo &
   (this, &
    iscompzchanged, &
    nnch, &
	nsp, &
    nband, &
    mxconn, &
    ngas, &
    theta, &
    dtime, &
    iconn, &
    efloww_ij, &
    eflowg1_ij, &
    eflowg2_ij, &
    istrpgas, &
    tindep, &
	issolvermrmt, &
	imethod, &
	mxnbox, &
	iconnbox, &
	alphabox, &
	phibox, &
	msg, &
    iserror, &
	iouwinfo)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:Compute the concentrations of equilibrium minerals 
! writting the transport equations.
!
!   $Arguments:
!
 
type(t_cheproo), intent(inout)                       :: this

integer, intent(in)                                  :: nsp

integer, intent(in)                                  :: nnch

integer, intent(in)                                  :: ngas

integer, intent(in)                                  :: mxconn

integer, intent(in)                                  :: nband

integer, intent(in)                                  :: mxnbox

integer, intent(in), dimension(mxnbox,this%numnch)   :: iconnbox        ! Correspondence between boxes and nodes

real*8, intent(in), dimension(mxnbox,this%numnch)    :: alphabox        ! Alpha for boxes 

real*8, intent(in), dimension(mxnbox,this%numnch)    :: phibox          ! Porosity for boxes 

integer, intent(in)                                  :: imethod         ! Index of reactive multi-rate mass transfer method 

logical, intent(in)                                  :: issolvermrmt 

real*8, intent(in)                                   :: theta

real*8, intent(in)                                   :: dtime

integer, intent(in),dimension(mxconn,nnch)           :: iconn

real*8, intent(in),dimension(nsp,nnch)               :: tindep 

real*8, intent(in),dimension(nnch,2*nband+1)         :: efloww_ij

real*8, intent(in),dimension(nnch,2*nband+1)         :: eflowg1_ij

real*8, intent(in),dimension(nnch,2*nband+1)         :: eflowg2_ij

logical, intent(in)                                  :: istrpgas

character(len=*), intent(out)                        :: msg

logical, intent(out)                                 :: iserror ! If .true. there was error 

logical, intent(out)                                 :: iscompzchanged ! If .true. the components zone was changed 

integer, intent(in), optional                        :: iouwinfo     ! Output unit
 
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
 jthnch, &
 ithnch, &
 ic, &
 ncon, &
 iband, &
 nbox, &
 ibox, &
 ndim 
real*8                        :: &
 alpha, &
 fii, &
 eij, &
 volii
real*8, pointer               :: &
 setre(:,:) => null (), &
 sktrk(:) => null (), &
 c(:) => null (), &
 cmob(:,:) => null (), &
 indep(:) => null (), &
 psetre(:) => null (), &
 f(:) => null (), &
 f1(:) => null (), &
 df1(:,:) => null ()
logical                       :: &
 ischanged, &
 isconvergence, &
 haveiouwinfo 
type(t_nodalchemistry), pointer        :: &
 pnchk => null (), &
 pnchk1 => null (), &
 pboxk1 => null (), &
 pboxk => null ()
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Check the number of nodal chemistry objects 
!%------------------------------------------------------------
if(nnch/=this%numnch) then
 msg='Error in number of nodal chemistry objects'
 goto 10 
end if
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
haveiouwinfo=present(iouwinfo)
!%------------------------------------------------------------
iscompzchanged=.false.
!%------------------------------------------------------------
!% Allocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (setre,nsp,nnch,.true.)
call check_pointer_ (indep,nsp,.true.)
!%------------------------------------------------------------
!% Loop for nodal chemistry objects 
!%------------------------------------------------------------
do jthnch=1,nnch
  
  ncon=iconn(1,jthnch)
  
  iband=nband+1
  eij=efloww_ij(jthnch,iband)
  psetre => setre(1:nsp,jthnch) 
  pnchk1 => this%pnodalchemk1(jthnch)%ptr
  pnchk => this%pnodalchemk(jthnch)%ptr

!%------------------------------------------------------------
!% Get chemical info
!%------------------------------------------------------------
  call get_chem_info_(pnchk1,iserror,c=c,cmob=cmob,sktrk=sktrk,omgwfree=fii,volnch=volii)
  if (iserror) goto 20
!%------------------------------------------------------------
!%  Add f 
!%------------------------------------------------------------
  if (issolvermrmt) then

     nbox=iconnbox(1,jthnch)
     
	 do ibox=1,nbox
      
	     pboxk1 => this%prmrmtboxk1(iconnbox(ibox+1,jthnch))%ptr
         pboxk  => this%prmrmtboxk(iconnbox(ibox+1,jthnch))%ptr

         call compute_f_df_rmrmt_ &
          (pboxk1, &
           f1, &
           df1, &
           ndim, &
           pboxk, &
           pnchk1, &
           pnchk, &
           theta, &
           alphabox(ibox,jthnch), &
           phibox(ibox,jthnch), &
           dtime, &
           isconvergence, &
           1.0d-4, &
           1.0d-10, &
           100, &
           10, &
		   imethod, &
           iserror, &
		   fsp=f)

	       if (iserror) goto 20 
       
!%------------------------------------------------------------
! Add f
!%------------------------------------------------------------
           alpha=volii
		   call add_ (psetre,alpha,f,nsp)
       

     end do

  end if 
!%------------------------------------------------------------
! Add independent term
!
!   -ti
!
!%------------------------------------------------------------
  alpha=-1.0d0
  call add_ (psetre,alpha,tindep(1:nsp,jthnch),nsp)
!%------------------------------------------------------------
! Add storage term
!
!   ck+1
!
!%------------------------------------------------------------
  alpha=fii/dtime
  call add_ (psetre,alpha,c,nsp)
!%------------------------------------------------------------
! Add kinetic term
!
! -theta*sktrk
!
!%------------------------------------------------------------
  alpha=-fii  
  call add_ (psetre,alpha,sktrk,nsp)
!%------------------------------------------------------------
! Add himself contribution
!%------------------------------------------------------------
  alpha=eij
  call add_ (psetre,alpha,cmob(1:nsp,1),nsp)
!%------------------------------------------------------------
! Add contribution of connected nodes
!%------------------------------------------------------------
  do ic=2,ncon
   ithnch=iconn(ic,jthnch)
   iband=nband+1-ithnch+jthnch
   eij=efloww_ij(ithnch,iband)
   alpha=eij 
   psetre => setre(1:nsp,ithnch) 
   call add_ (psetre,alpha,cmob(1:nsp,1),nsp)
  end do
 
end do
!%------------------------------------------------------------
! Set each nodal chemistry from Set re computed from transport 
! equations
!%------------------------------------------------------------
do ithnch=1,nnch
 pnchk1 => this%pnodalchemk1(ithnch)%ptr 
 pnchk => this%pnodalchemk(ithnch)%ptr 
 psetre => setre(1:nsp,ithnch)
 call set_from_setre_ (pnchk1,pnchk,ischanged,psetre,nsp,dtime,iserror)
 if (iserror) goto 20
!%------------------------------------------------------------
! If change de components zone, then update derivates
!%------------------------------------------------------------
 if (ischanged) then
    if (this%iswriteinfo.and.haveiouwinfo) then
      write (iouwinfo,*) '-------------------------------------------------------'
	  write (iouwinfo,*) 'COMPONENTS MATRIX WAS CHANGED IN NODAL CHEMISTRY',ithnch
	  write (iouwinfo,*) '-------------------------------------------------------'
    end if  
	call update_derivatives_ (pnchk1,iserror)
    if (iserror) goto 20
    iscompzchanged=.true.
 end if
!%-----------------------------------------------------------
end do
!%------------------------------------------------------------ 
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers
!%------------------------------------------------------------
call check_pointer_ (setre,1,1,.false.)
call check_pointer_ (c,1,.false.)
call check_pointer_ (sktrk,1,.false.)
call check_pointer_ (cmob,1,1,.false.)
call check_pointer_ (indep,1,.false.)
call check_pointer_ (f,1,.false.)
call check_pointer_ (f1,1,.false.)
call check_pointer_ (df1,1,1,.false.)
!%------------------------------------------------------------
!% Nullify local pointers 
!%------------------------------------------------------------
pnchk => null ()
pnchk1 => null ()
psetre => null ()
pboxk1 => null ()
pboxk  => null ()
if (iserror) goto 10
!%-------------------------------------------------------------
return
 
10 continue 
print *,'***************************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: compute_eq_min_from_trp_'
print *, msg
print *,'***************************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine add_jacobian_dsa_cheproo &
   (jacobian, &
    nunk, &
    iposjthnch, &
    iposithnch, &
    nbandjac, &
    der, &
    alpha, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Add information in the jacobian matrix for DSA 
! formulation. 
!
!   $Arguments:
!
 

integer, intent(in)                   :: iposjthnch

integer, intent(in)                   :: iposithnch

integer, intent(in)                   :: nbandjac

integer, intent(in)                   :: nunk

real*8, pointer, dimension(:,:)       :: jacobian 

type(t_array), intent(in)             :: der

real*8, intent(in)                    :: alpha

character(len=*), intent(out)         :: msg

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
integer                               :: &
 ijac1, &
 ijac2, &
 icoljthnch, &
 icol 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
do icoljthnch=1,der%ncol
 
  icol=iposjthnch+icoljthnch-1
 
 call get_jacobpos_band_ &
   (ijac1, &
    ijac2, &
    icol, &
    iposithnch, &
    der%nrow, &
    nbandjac, &
    msg, &
    iserror)
 
 if (iserror) goto 10
 
 call add_ &
 (jacobian(ijac1:ijac2,icol), &
  alpha, &
  der%array(:,icoljthnch), &
  der%nrow)
 
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
subroutine get_jacobpos_band_cheproo &
   (irowjac1, &
    irowjac2, &
    icolfull, &
    irowfull, &
    nrowfull, &
    nband, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Give the position in the jacobian for badn stored
!
!   $Arguments:
!
 
integer, intent(out)                  :: irowjac1

integer, intent(out)                  :: irowjac2

integer, intent(in)                   :: icolfull

integer, intent(in)                   :: irowfull

integer, intent(in)                   :: nrowfull

integer, intent(in)                   :: nband

character(len=*), intent(out)         :: msg

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
msg=''
iserror=.false.
!%------------------------------------------------------------
irowjac1=2*nband+1+irowfull-icolfull
irowjac2=irowjac1+nrowfull-1
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
subroutine compute_in_out_mol_cheproo &
   (this, &
    molinout, &
    ncomp, &
    nnch, &
    dtime, &
    caudalw, &
    caudalrechw, &
    itypeboundw, &
    isetboundw, &
    isetrechw, &
    gammaboundw, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the input/output mol of components when solving 
! reactive transport step. 
! This service is performed by mol balance of components. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)             :: this         ! Type CHEPROO variable 

real*8, pointer                          :: molinout(:)  ! Vector of in/out of mol 

integer, intent(in)                      :: nnch         ! Number of nodal chemistries

integer, intent(out)                     :: ncomp        ! Number of components

real*8, intent(in)                       :: dtime        ! Time increment

real*8, intent(in), dimension(nnch)      :: caudalw

real*8, intent(in), dimension(nnch)      :: caudalrechw

real*8, intent(in), dimension(nnch)      :: gammaboundw

integer, intent(in), dimension(nnch)     :: isetboundw

integer, intent(in), dimension(nnch)     :: isetrechw

integer, intent(in), dimension(nnch)     :: itypeboundw

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
character(len=100)                   :: &
 msg
integer                              :: &
 type, &
 inch, &
 nsp, &
 iwater
real*8, pointer                      :: &
 cmob1(:,:) => null(), &
 cmob2(:,:) => null(), &
 vector(:) => null()
real*8                               :: &
 floww, &
 flowrech, &
 gamma 
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Check the number of nodal chemistry objects 
!%--------------------------------------------------------------
if (nnch/=this%numnch) then
  msg='Error in number of nodal chemistry objects'
  goto 10
end if
!%--------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system'
 goto 10
end if
!%--------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numbase=ncomp,numsp=nsp)
if (iserror) goto 10
!%--------------------------------------------------------------
call check_pointer_ (vector,nsp,.true.)
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Loop to nodal chemistry objects 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
do inch=1,nnch
 
    floww=caudalw(inch)
    type=itypeboundw(inch)
    iwater=isetboundw(inch)
    gamma=gammaboundw(inch)
 
    if (iwater>0) then
     call get_chem_info_ (this%pwater(iwater)%ptr,iserror,cmob=cmob1)
     if (iserror) goto 20 
    end if
    
	call get_chem_info_ (this%pnodalchemk1(inch)%ptr,iserror,cmob=cmob2)
	if (iserror) goto 20 
!%--------------------------
!%--------------------------
!%--------------------------
    if (type==1.and.iwater>0) then
       vector=vector+cmob1(:,1)*dtime
    else if (type==2.and.iwater>0) then
       if (floww>0.0d0) then
         vector=vector+floww*cmob1(:,1)*dtime
       else
         vector=vector+floww*cmob2(:,1)*dtime
       end if
    else if (type==3.and.iwater>0) then
        vector=vector+cmob1(:,1)*dtime
    else if (type==4.and.iwater>0) then
        vector=vector+gamma*(cmob1(:,1)-cmob2(:,1))*dtime
    else if (type==5.and.iwater>0.and.floww>0.0d0) then
        vector=vector+floww*cmob1(:,1)*dtime
    else if (type==0.and.floww<0.0d0) then
        vector=vector+floww*cmob2(:,1)*dtime
    end if
!%--------------------------
!%--------------------------
!%--------------------------
    iwater=isetrechw(inch)
    floww=caudalrechw(inch)
 
    if (iwater>0) then
     call get_chem_info_(this%pwater(iwater)%ptr,iserror,cmob=cmob1)
     if (iserror) goto 20 
     vector=vector+floww*cmob1(:,1)*dtime
    end if
 
 
 
end do
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Finishing loop to nodal chemistry objects 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Make the lineal transformation multiplying by components 
!% matrix
!%--------------------------------------------------------------
call make_lin_trf_ (this%pchemsys(1)%ptr,molinout,vector,0,iserror)
!%--------------------------------------------------------------
20 continue 
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (cmob1,1,1,.false.)
call check_pointer_ (cmob2,1,1,.false.)
call check_pointer_ (vector,1,.false.)
if (iserror) goto 10
!%--------------------------------------------------------------
return
 
10 continue 
print *,'******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: compute_in_out_mol_'
print *, msg
print *,'******************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_kin_mol_cheproo &
   (this, &
    molkin, &
    ncomp, &
    dtime, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the total mol per component transferred due 
! kinetic reactions during time increment (dtime) 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)             :: this     ! Type CHEPROO variable. 

real*8, pointer, dimension(:)            :: molkin    

integer, intent(out)                     :: ncomp    ! Number of components. 

real*8, intent(in)                       :: dtime    ! Time increment. 

logical, intent(out)                     :: iserror  ! iserror=true, then there was an error. 
 
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
integer                              :: &
 inch
real*8, pointer                      :: &
 molkinith(:) => null()
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------
if (this%numchemsys==0) then
 msg='Error, not defined chemical system'
 goto 10
end if
!%--------------------------------------------------------------
call get_chem_info_ (this%pchemsys(1)%ptr,iserror,numbase=ncomp)
if (iserror) goto 10
!%--------------------------------------------------------------
!% Allocate and initialize mol kin vector 
!%--------------------------------------------------------------
call check_pointer_ (molkin,ncomp,.true.)
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Loop to nodal chemistry objects 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
do inch=1,this%numnch
 
 call get_chem_info_ (this%pnodalchemk1(inch)%ptr,iserror,molkin=molkinith)
 if (iserror) goto 20 
 molkin = molkin + molkinith * dtime 

end do
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Finishing loop to nodal chemistry objects 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
20 continue
!%--------------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------------- 
call check_pointer_ (molkinith,1,.false.)
if (iserror) goto 10
!%--------------------------------------------------------------
return
 
10 continue 
print *,'******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: compute_kin_mol_'
print *, msg
print *,'******************************'
iserror=.true.
return
 
end subroutine


!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_aq_sps_pos_cheproo &
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
 
type (t_cheproo), intent(in)                    :: this     ! Type CHEPROO variable. 

integer, intent(out)                             :: isps1 ! Global index of the first aqueous species. 

integer, intent(out)                             :: isps2 ! Global index of the last aqueous species. 

logical, intent(out)                             :: iserror ! iserror=true, there was an error. 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: By the moment we have only one chemical system so this subroutine
!   gives the
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

call get_aq_sps_pos_(this%pchemsys(1)%ptr,isps1,isps2,iserror)

!%--------------------------------------------------------------
return
 
10 continue 
print *,'******************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_aq_sps_pos_'
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
subroutine write_signature_cheproo &
  (this, &
   ioutput)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write the signature of CHEPROO object. 
!
!   $Arguments:
!

type(t_cheproo), intent(in)    :: this      ! Type CHEPROO variable 

integer, intent(in)            :: ioutput   ! Output unit 
 
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
  
!%------------------------------------------------------------------------
write(ioutput,*) '**************************************************'
write(ioutput,*) '**************************************************'
write(ioutput,*) '*                                                *'
write(ioutput,*) '*    @@@@@ @   @ @@@@@ @@@@@ @@@@@ @@@@@ @@@@@   *'
write(ioutput,*) '*    @     @   @ @     @   @ @   @ @   @ @   @   *'
write(ioutput,*) '*    @     @@@@@ @@@@@ @@@@@ @@@@@ @   @ @   @   *'  
write(ioutput,*) '*    @     @   @ @     @     @@    @   @ @   @   *'   
write(ioutput,*) '*    @@@@@ @   @ @@@@@ @     @ @   @@@@@ @@@@@   *'
write(ioutput,*) '*                                                *'
write(ioutput,*) '*        CHEmical PRocesses Object-Oriented      *'
write(ioutput,*) '*                                                *'
write(ioutput,*) '*                    Written by                  *'
write(ioutput,*) '*                                                *'
write(ioutput,*) '*             Sergio Andrés Bea Jofre            *'
write(ioutput,*) '*                                                *'
write(ioutput,*) '*              Scientific Assessors :            *'
write(ioutput,*) '*                                                *'
write(ioutput,*) '*                  Jesus Carrera                 *'
write(ioutput,*) '*                  Carlos Ayora                  *'
write(ioutput,*) '*                 Francesc Batlle                *'
write(ioutput,*) '*                  M.W. Saaltink                 *'
write(ioutput,*) '*                                                *'
write(ioutput,*) '*                                                *'
write(ioutput,*) '*                                                *'
write(ioutput,*) '*                     CSIC-UPC                   *'
write(ioutput,*) '*                 Barcelona-Spain                *'
write(ioutput,*) '*                       2008                     *'
write(ioutput,*) '*    Reference:                                  *' 
write(ioutput,*) '*              S.A. Bea, J. Carrera, C. Ayora,   *' 
write(ioutput,*) '*              F. Batlle and Saaltink, M.        *' 
write(ioutput,*) '*              A Fortran90 Object-Oriented module*' 
write(ioutput,*) '*              to solve chemical processes in the*'
write(ioutput,*) '*              Earth Sciences models. Computers &*'
write(ioutput,*) '*              Geosciences 35(6), 1098-1112      *' 
write(ioutput,*) '**************************************************'
write(ioutput,*) '**************************************************'
write(ioutput,*) '* CHEPROO is an object-oriented tool specialized *'
write(ioutput,*) '* in chemical processes. It presents features    *'
write(ioutput,*) '* such as:                                       *'
write(ioutput,*) '* - Extensibility.                               *'
write(ioutput,*) '* - Transportability.                            *'
write(ioutput,*) '* - Easy to link to other Earth Sciences         *'  
write(ioutput,*) '*   simulators.                                  *'
write(ioutput,*) '**************************************************'
write(ioutput,*) '**************************************************'
!&------------------------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%****************Public subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_eq_min_from_SeRe_cheproo &
   (this, &
    Sere, &
    nnch, &
	nsp, &
    dtime, &
    iscompzchanged, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:   Compute the concentrations of equilibrium minerals 
!                   from the SeRe
!
!   $Arguments:
!
 
type(t_cheproo), intent(inout):: this

real*8,dimension(:),pointer   :: Sere

real*8, intent(in)            :: dtime

integer, intent(in)           :: nsp

integer, intent(in)           :: nnch

logical, intent(out)          :: iserror ! If .true. there was error 

logical, intent(out)          :: iscompzchanged ! If .true. the components zone was changed 
 
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
 jthnch, &
 ithnch, &
 !ic, &
 !ncon, &
 !iband,&
 istart, &
 iend, &
 nraqu
real*8, pointer               :: &
 psetre(:) => null ()
logical                       :: &
 ischanged 
type(t_nodalchemistry), pointer        :: &
 pnchk => null (), &
 pnchk1 => null ()
 character(len=100)                   :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Check the number of nodal chemistry objects 
!%------------------------------------------------------------
if(nnch/=this%numnch) then
 msg='Error in number of nodal chemistry objects'
 goto 10 
end if
!%------------------------------------------------------------
!% Check SeRe size
!%------------------------------------------------------------
!we get the start and the end of the aqueous species in the concentration vector 
call get_aq_sps_pos_(this,istart,iend, iserror)
Nraqu=iend-istart+1

if(nnch*Nraqu/=size(SeRe)) then
 msg='Error in number of aqueous species'
 goto 10 
end if


allocate(psetre(nsp))
psetre=0

iscompzchanged=.false.
!%------------------------------------------------------------
! Set each nodal chemistry from Set re 
!%------------------------------------------------------------
do ithnch=1,nnch
 pnchk1 => this%pnodalchemk1(ithnch)%ptr 
 pnchk => this%pnodalchemk(ithnch)%ptr 
 psetre(istart:iend) = SeRe( Nraqu*(ithnch-1)+istart : Nraqu*(ithnch-1)+iend )
 call set_from_setre_ (pnchk1,pnchk,ischanged,psetre,nsp,dtime,iserror)
 if (iserror) goto 10
!%------------------------------------------------------------
! If change de components zone, then update derivates
!%------------------------------------------------------------
 if (ischanged) then
    call update_derivatives_ (pnchk1,iserror)
    if (iserror) goto 10
    iscompzchanged=.true. 
 end if
!%-----------------------------------------------------------
 
end do


!%------------------------------------------------------------
!% Nullify local pointers 
!%------------------------------------------------------------
pnchk => null ()
pnchk1 => null ()
deallocate(psetre)

if (iserror) goto 10
!%-------------------------------------------------------------
return
 
10 continue 
print *,'***************************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: compute_eq_min_from_SeRe_'
print *, msg
print *,'***************************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_compdef_cheproo &
   (this, &
    isbackwards, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy the component definition from one of the nodal chemistry
!                   group to the other
! isbackwards=true => Copy the component definition of k+1 in nodal chemistry in k.
! isbackwards=false => Copy the component definition of k in nodal chemistry in k+1. 
!
!   It's used when we want to evaluate the old concentrations in new component definiton
!   $Arguments:
!
 
type (t_cheproo), intent(inout)       :: this

logical, intent(in)                   :: isbackwards 
logical, intent(out)                  :: isError 
 
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
 inch, &
 hashcompz 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------


if (this%numnch>0) then
!%------------------------------------------------------------
! isbackward=true  then compDef Nodalchem(k)= compDef Nodalchem(k+1)
!%------------------------------------------------------------
 if (isbackwards) then 
  do inch=1,this%numnch
    !call get_comp_zone_index_(this%pnodalchemk1(inch)%ptr,ithcompz,iserror)
    call   get_chem_info_(this%pnodalchemk1(inch)%ptr,iserror,hashcompz=hashcompz)
     if (iserror) goto 10
    call set_(this%pnodalchemk(inch)%ptr,iserror,hashcompz=hashcompz)
     if (iserror) goto 10
  enddo
 else
!%------------------------------------------------------------
! isbackward=false  then compDef nodalchem(k+1)= compDef nodalchem(k)
!%------------------------------------------------------------
  do inch=1,this%numnch
    !call get_comp_zone_index_(this%pnodalchemk(inch)%ptr,ithcompz,iserror)
    call   get_chem_info_(this%pnodalchemk(inch)%ptr,iserror,hashcompz=hashcompz)
     if (iserror) goto 10
    call set_(this%pnodalchemk1(inch)%ptr,iserror,hashcompz=hashcompz)
     if (iserror) goto 10
  enddo
 end if
 
end if
!%------------------------------------------------------------
return
10 continue 
print *,'***************************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Error in copy_compdef_cheproo'
print *,'***************************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_cpri_water_cheproo &
   (this, &
    cpri, &
    npri, &
    water, &
    compDeffNch, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_cheproo), target       :: this

real*8, pointer, dimension(:)              :: cpri 

integer, intent(out)                       :: npri

logical, intent(out)                       :: iserror 

integer, intent(in)                        :: water

integer, intent(in)                        :: compDeffNch !nodale chemistry were the components are defined
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
!type(t_nodalchemistry), pointer    :: &
! pnch => null ()
 integer                           :: &
 hashcompz
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%----------------------------------------------------------
if (water<=0.or.water>this%numw) then
 msg='Error in water index, index='
 call add_ (msg,water) 
 goto 10
end if
!We get the component definiton from the nodal chemistry
!%----------------------------------------------------------
call get_chem_info_(this%pnodalchemk1(compDeffNch)%ptr,iserror,hashcompz=hashcompz)
call set_(this%pwater(water)%ptr,iserror,hashcompz=hashcompz)
!%----------------------------------------------------------
call get_chem_info_ (this%pwater(water)%ptr,iserror,cpri=cpri,npri=npri)
if (iserror) goto 10
!%----------------------------------------------------------
!now we put the water again with definition equal to 0
call set_(this%pwater(water)%ptr,iserror,hashcompz=0)
!pnch => null ()
!%----------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_cpri_water_'
print *,msg
print *,'**********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_umin_cheproo &
   (this, &
    umin, &
    npri, &
    ith, &
	jth, &
	isbackward, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return umin[npri]=Uith*cminjth, where ith and jth are referring
! to nodal chemistry index for U and cmin respectively. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this       ! Type CHEPROO variable

integer, intent(in)                   :: npri       ! Number of components

real*8, intent(inout), dimension(npri):: umin       ! Mineral part of the component 

integer, intent(in)                   :: ith        ! Nodal chemistry index for U

integer, intent(in)                   :: jth        ! Nodal chemistry index for cmin

logical, intent(in)                   :: isbackward ! if true, then use the nodalchemistry list in k
                                                    ! if false, then use the nodalchemistry list in k+1

logical, intent(out)                  :: iserror    ! iserror=true, then there was an error 
 
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
 npri1, &
 hashcompz
real*8, pointer                       :: &
 cmin(:) => null (), &
 umin1(:) => null ()
character(len=100)                    :: &
 msg 
type(t_nodalchemistry), pointer       :: &
 pnchith => null (), &
 pnchjth => null ()
type(t_chemicalsystem), pointer       :: &
 pchemsysith => null()
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check the ith and jth nodal chemistry indices
!%-------------------------------------------------------------
if (ith<=0.or.ith>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,ith) 
 goto 10
end if
if (jth<=0.or.jth>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,jth) 
 goto 10
end if
!%------------------------------------------------------------
!% Assign the corresponding nodal chemistries pointers
!%------------------------------------------------------------
if (isbackward) then
 pnchith => this%pnodalchemk(ith)%ptr
 pnchjth => this%pnodalchemk(jth)%ptr
else
 pnchith => this%pnodalchemk1(ith)%ptr
 pnchjth => this%pnodalchemk1(jth)%ptr
end if 
!%------------------------------------------------------------
!% Get components hash index of ith nodal chemistry object
!%------------------------------------------------------------
call get_chem_info_ (pnchith,iserror,hashcompz=hashcompz,npri=npri1,chemsys=pchemsysith)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Check the number of components of the ith components zone
!%------------------------------------------------------------
if (npri/=npri1) then
 msg='Error in number of components'
 goto 20 
end if
!%------------------------------------------------------------
!% Get concentrations vector of jth nodal chemistry object 
!%------------------------------------------------------------
call get_chem_info_ (pnchjth,iserror,cmin=cmin)
if (iserror) goto 20 

call make_lin_trf_ (pchemsysith,umin1,cmin,hashcompz,iserror)
if (iserror) goto 20 
!%-------------------------------------------------------------
!% Copy 
!%-------------------------------------------------------------
umin=umin1
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!% Deallocate local pointers
!%-------------------------------------------------------------
call check_pointer_ (cmin,1,.false.)
call check_pointer_ (umin1,1,.false.)
pnchith => null()
pnchjth => null ()
pchemsysith => null()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_umin_'
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
subroutine get_umineq_cheproo &
   (this, &
    umineq, &
    npri, &
    ith, &
	jth, &
	isbackward, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return umineq[npri]=Uith*cmineqjth, where ith and jth are referring
! to nodal chemistry index for U and cmin respectively. 
!   This service normally returns zero, because in the component definition all minerals
!   in equilibrium have been eliminated. But when the component definition change 
!   (for example a mineral disappear) it values may be different from zero. 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)          :: this       ! Type CHEPROO variable

integer, intent(in)                   :: npri       ! Number of components

real*8, intent(inout), dimension(npri):: umineq       ! mineral part of the component 

integer, intent(in)                   :: ith        ! Nodal chemistry index for U

integer, intent(in)                   :: jth        ! Nodal chemistry index for cmin

logical, intent(in)                   :: isbackward ! if true, then use the nodalchemistry list in k
                                                    ! if false, then use the nodalchemistry list in k+1

logical, intent(out)                  :: iserror    ! iserror=true, then there was an error 
 
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
 npri1, &
 hashcompz
real*8, pointer                       :: &
 cmineq(:) => null (), &
 umin1(:) => null ()
character(len=100)                    :: &
 msg 
type(t_nodalchemistry), pointer       :: &
 pnchith => null (), &
 pnchjth => null ()
type(t_chemicalsystem), pointer       :: &
 pchemsysith => null()
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check the ith and jth nodal chemistry indices
!%-------------------------------------------------------------
if (ith<=0.or.ith>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,ith) 
 goto 10
end if
if (jth<=0.or.jth>this%numnch) then
 msg='Error in nodal chemistry index='
 call add_ (msg,jth) 
 goto 10
end if
!%------------------------------------------------------------
!% Assign the corresponding nodal chemistries pointers
!%------------------------------------------------------------
if (isbackward) then
 pnchith => this%pnodalchemk(ith)%ptr
 pnchjth => this%pnodalchemk(jth)%ptr
else
 pnchith => this%pnodalchemk1(ith)%ptr
 pnchjth => this%pnodalchemk1(jth)%ptr
end if 
!%------------------------------------------------------------
!% Get components hash index of ith nodal chemistry object
!%------------------------------------------------------------
call get_chem_info_ (pnchith,iserror,hashcompz=hashcompz,npri=npri1,chemsys=pchemsysith)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Check the number of components of the ith components zone
!%------------------------------------------------------------
if (npri/=npri1) then
 msg='Error in number of components'
 goto 20 
end if
!%------------------------------------------------------------
!% Get concentrations vector of jth nodal chemistry object 
!%------------------------------------------------------------
call get_chem_info_ (pnchjth,iserror,cmineq=cmineq)
if (iserror) goto 20 

call make_lin_trf_ (pchemsysith,umin1,cmineq,hashcompz,iserror)
if (iserror) goto 20 
!%-------------------------------------------------------------
!% Copy 
!%-------------------------------------------------------------
umineq=umin1
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!% Deallocate local pointers
!%-------------------------------------------------------------
call check_pointer_ (cmineq,1,.false.)
call check_pointer_ (umin1,1,.false.)
pnchith => null()
pnchjth => null ()
pchemsysith => null()
if (iserror) goto 10 
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'CHEPROO:'
print *,'Name:', this%name
print *,'Service: get_umin_'
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
subroutine get_name_comp_ith_cheproo &
   (this, &
    namebase, &
    node, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of component for the it nodal chemistry
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                   :: this

integer, intent(in)                            :: node

character(len=100), dimension(:),pointer       :: namebase

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
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
!%--------------------------------------------------------------
call get_chem_info_(this%pnodalchemk1(node)%ptr,iserror,namebase=namebase)
if (iserror) goto 10 
!%--------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: get_name_comp_ith_'
print *,'********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_aqdensity_ith_cheproo &
    (this       , &
    iserror     , &
    ithnode     , &
    isbackwards , &
    density     , &
    ddensitydc)
    

type(t_cheproo),intent(in)              :: this         ! Type nodal chemistry variable. 

logical,intent(out)                     :: IsError      ! iserror=true, then there was an error

integer,intent(in)                      :: ithnode      ! Node were the density will be calculated

logical,intent(in)                      :: isbackwards  ! If true evaluate in K, if false evaluate in K1 nodal chemistry

real*8,intent(out)                      :: density      ! Density     

real*8,pointer,dimension(:),optional    :: ddensitydc   ! Density derivative

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 integer::phasecode
 type(t_pnodalchemistry)::pnodal
!-------------------------------------------------------------------------
!
!   $code
!

if (isbackwards) then
    pnodal%ptr=> this%pnodalchemk(ithnode)%ptr
else
    pnodal%ptr=> this%pnodalchemk1(ithnode)%ptr
endif

phasecode=ip_liquid_phase

if (present(ddensitydc)) then
    call compute_density_ith_(pnodal%ptr,phasecode,density,iserror,ddensitydc)
    if (isError) goto 10
else
    call compute_density_ith_(pnodal%ptr,phasecode,density,iserror)
    if (isError) goto 10
endif


!%----------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Cheproo:'
print *,'Service: compute_aqdensity_ith_'
print *,'***********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine add_f_df_dsa_ith_cheproo &
   (this, &
    jacobian, &
    residual, &
    nunk, &
    nbandjac, &
	nsp, &
    typesto, &
	nbox, &
    theta, &
    dtime, &
	posunk, &
	imethod, &
	iconnbox, &
	alpha, &
	phi, &
	tolunk, &
	tolres, &
	mxiter, &
    mxiterminchange, &
	isconvergence, &
	ithnch, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Add f and df in the jacobian and residual 
!
!   $Arguments:
!
 
type (t_cheproo), intent(in)                           :: this

integer, intent(in)                                    :: nunk

integer, intent(in)                                    :: nbandjac

integer, intent(in)                                    :: nbox

integer, intent(in)                                    :: nsp

integer, intent(in)                                    :: typesto

integer, intent(in)                                    :: imethod 

real*8, intent(in)                                     :: theta                     ! Temporal weight 

real*8, intent(in)                                     :: dtime                     ! Time increment
 
real*8, pointer, dimension(:,:)                        :: jacobian        ! Jacobian matrix

real*8, pointer, dimension(:)                          :: residual        ! Residual 

integer, intent(in), dimension(2,this%numnch)          :: posunk

integer, intent(in)                                    :: ithnch 

real*8, intent(in), dimension(nbox)                    :: alpha           ! Alpha for boxes 

real*8, intent(in), dimension(nbox)                    :: phi             ! Porosity for boxes 

integer, intent(in), dimension(nbox)                   :: iconnbox        ! Correspondence between boxes and nodes

real*8, intent(in)                                     :: tolunk

real*8, intent(in)                                     :: tolres 

integer, intent(in)                                    :: mxiter

integer, intent(in)                                    :: mxiterminchange 

logical, intent(out)                                   :: isconvergence   ! .true. if there was convergence 

logical, intent(out)                                   :: iserror         ! .true. if there is error 
 
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
type(t_vector)               :: &
 f
type(t_array)                :: &
 df
real*8                       :: &
 volii
real*8, pointer              :: &
 f1(:) => null (), &
 df1(:,:) => null ()
type(t_nodalchemistry), pointer :: &
 pmobk => null (), &
 pmobk1 => null (), &
 pboxk => null (), &
 pboxk1 => null ()
integer                         :: &
 ibox, &
 ncomp, &
 ipos1, &
 ipos2, &
 ndim  
!%------------------------------------------------------------
iserror=.false. 
msg=''
!%------------------------------------------------------------
!% Loop for nodal chemistry objects 
!%------------------------------------------------------------
pmobk1 => this%pnodalchemk1(ithnch)%ptr
pmobk  => this%pnodalchemk(ithnch)%ptr
ipos1=posunk(1,ithnch)
ipos2=posunk(2,ithnch)
!%------------------------------------------------------------
!% 
!%------------------------------------------------------------
call get_chem_info_ (pmobk1,iserror,npri=ncomp,volnch=volii)
if (iserror) goto 20 
!%------------------------------------------------------------
!% Allocate and initialice local variables 
!%------------------------------------------------------------
call check_pointer_ (f%vector,ncomp,.true.)
call check_pointer_ (df%array,ncomp,ncomp,.true.)
df%nrow=ncomp
df%ncol=ncomp 
!%------------------------------------------------------------
!% Loop for boxes 
!%------------------------------------------------------------   
do ibox=1,nbox
      
	  pboxk1 => this%prmrmtboxk1(iconnbox(ibox))%ptr
      pboxk  => this%prmrmtboxk(iconnbox(ibox))%ptr

      call compute_f_df_rmrmt_ &
          (pboxk1, &
           f1, &
           df1, &
           ndim, &
           pboxk, &
           pmobk1, &
           pmobk, &
           theta, &
           alpha(ibox), &
           phi(ibox), &
           dtime, &
           isconvergence, &
           tolunk, &
           tolres, &
           mxiter, &
           mxiterminchange, &
		   imethod, &
           iserror)
           

	   if (iserror.or..not.isconvergence) goto 20 
       
       f%vector = f%vector + f1
       df%array = df%array + df1
       

end do  
!%-------------------------------------------------------
!% Add f in the residual 
!%-------------------------------------------------------   
call add_ (residual(ipos1:ipos2),volii,f%vector,ncomp)
!%-------------------------------------------------------------
!% Add df in the jacobian matrix 
!%-------------------------------------------------------------
call add_jacobian_dsa_ (jacobian,nunk,ipos1,ipos1,nbandjac,df,volii,msg,iserror)
if (iserror) goto 10
!%-------------------------------------------------------
20 continue 
!%-------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------
call check_pointer_ (f1,1,.false.)
call check_pointer_ (f%vector,1,.false.)
call check_pointer_ (df1,1,1,.false.)
call check_pointer_ (df%array,1,1,.false.)
pmobk => null ()
pmobk1 => null ()
pboxk => null ()
pboxk1 => null ()
!%-------------------------------------------------------
if (iserror) goto 10 
!%-------------------------------------------------------
return
 
10 continue 
print *,'****************************'
print *,'CHEPROO:'
print *,'Name:',this%name
print *,'Service: add_f_df_dsa_ith_'
print *, msg
print *,'****************************'
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
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_cheproo
