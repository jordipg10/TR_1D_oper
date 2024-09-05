module m_nodalchemistry
!-------------------------------------------------------------------------
!
!>   $Description: Nodal chemistry class constitutes the basic operational 
!> entity of the geochemical process and contains the compositional information. 
!> Nodal chemistry objects may be used to simulate the chemistry of parcels of water.
!> These parcels may represent samples, batch experiments, portions of a contaminated 
!> soil or aquifer, etc. 
!> Therefore, this class contains all direct and derived state variables generally used 
!> in standard geochemical calculations (e.g. concentrations, activity coefficients, ionic strength, 
!> temperature, mass of water, kinetic reaction rates, surface area of minerals, derivatives 
!> of state variables, sorption sites).
!> In addition, each nodal chemistry object is associated to one chemical system object, to perform 
!> any chemical operation (see Figure \ref{fig:chemsys}a).
!> Methods offered by nodal chemistry class can be classified in four group of methods: A) internal 
!> maintenance (e.g. create, destroy, set, read/write), B) access to compositional information encapsulated it 
!> (e.g. access to concentrations or total concentrations), C) reaction path, D) setting compositional information. 
!> Four reaction path methods (group of methods C) can be applied upon a nodal chemistry object: 1) titration, 
!> 2) addition of mineral species, 3) evaporation in closed and open systems (i.e. the precipitated mineral phases 
!> are removed and not allowed to redissolve), 4) heating/cooling solutions.
!> The compositional information in a nodal chemistry object (group of methods D) can be defined from 1) 
!> concentrations of primary species, 2) total concentrations in all present phases, 3) general characteristics 
!> of the solution (i.e. given pH, equilibrium with some mineral or gas phase, total concentrations, etc.). 
!> 
!>   $Use: m_chemicalsystem.
!> m_general_tools_cheproo.
!> m_constants.
!> flib_xpath.
!> flib_sax.
!> 
!
!>   $Author: Sergio Andrés Bea Jofré 
!
!>   $License: UPC-CSIC 
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_chemicalsystem
use m_general_tools_cheproo 
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
private   
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
public:: &
create_ &                          !> Create the nodal chemistry object. 
,destroy_ &                        !> Destroy the nodal chemistry object. 
,read_backup_ &                    !
,read_xml_ &                       !> Read and set the nodal chemistry object from xml file. 
,read_txt_ &                       !> Read and set the nodal chemistry object from ascII file. 
,set_ &                            !> Set general attributtes in the nodal chemistry object. 
,set_from_cpri_ &                  !> Set the concentrations (and dependent variables) from concentration of primary species. 
,set_from_u_ &                     !> Set the concentrations (and dependent variables) from total concentrations in all phases. 
,set_from_setre_ &                 !> Set the concentration of equilibrium mineral species from the equilibrium term of the transport equations.
,set_from_uaq_ &                   !> Set the concentrations (and dependent variables) from total concentrattions.
,set_from_solution_type_ &         !> Set the concentrations (and dependent variables) from type of solution (e.g. according pH, total concentrations, etc.). 
,set_from_c_ &                     !> Set from concentrations vector
,add_ &                            !> Add species or surfaces (for adsorption) into nodal chemistry object.  
,check_compzone_ &                 !> Check components zone associated to nodal chemistry
,check_and_switch_base_ &          !> Change and switch chemical base. 
,get_dusktrk_ &                    !> Compute U*Skt*drk/dc1 
,get_dumob_ &                      !> Return the U*dcmob/dc1
,get_duads_ &                      !> Return the U*dcads/dc1
,get_iumob_ &                      !> Return concentration of the ith mobile components
,get_chem_info_ &                  !> Return general chemical information 
,get_from_setre_ &                 !> Get x as: Transpose(Se)*x=b
,compute_density_ith_ &            !> Compute density (and derivatives if asked) for any phase
,compute_f_df_rmrmt_ &             !> Compute f and df for reactive multi-rate mass transfer 
,equilibrate_ &                    !> Equilibrate the nodal chemistry solution wit any mineral phase. 
,reaction_path_ &                  !> Reaction path services 
,reset_ &                          !> Reset nodal chemistry. 
,write_backup_ &
,write_ &                          !> Write in ascii the attributtes encapsulated in a nodal chemistry object.
,write_sps_ &                      !> Write in ascii file according list of name of species.  
,write_tot_ &                      !> Write components in txt file 
,write_volfracmin_ &               !> Write in ascII the volumetric fraction of mineral species 
,write_sat_ &                      !> Write in ascII the saturation of the solution with respect to mineral species. 
,write_min_area_ &                 !> Write in ascII the reactive surface of minerals 
,write_chemical_parameters_ &      !> Write general chemical parameters (e.g. pH, ionic strength, water activity)
,update_ &                         !> Update porosity and yourself 
,update_derivatives_ &             !> Update derivates 
,update_kinetic_ &                 !> Update kinetic terms
,update_and_check_ &               !> Update the nodal chemistry object and check convergence:
>                                   !> 1) Update the nodal chemistry object according increment of concentration of primary 
								   !>    species (delta cpri), update the nodalchemistry list in k+1 and check convergence.
								   !> 2) Update the nodal list in k+1 according total concentrations, compute raction term for sia
								   !>    and check convergence.  
,update_r_sia_ &                   !> Update reaction term for SIA from transport solution
,assignment(=) &                   !> Copy the nodal chemistry object in other nodal chemistry object
,operator(+) &                     !> Mix two nodal chemistries and store in other nodal chemistry object. 
,operator(*) &                     !> Multiply one nodal chemistry object by a real value. 
,set_iswcompbal_ &                 !> Set in the chemical system object if the mass of water must be 
>                                   !> considered in the speciation services
,update_min_area_ &                !> Update the reactive surface of kinetic minerals  
,update_txoh_                      !> Update txoh
!%------------------------------------------------------------------------
!%-------------------------------------------------------------------
private:: &                
read_xml_loc_ &
,reaction_path_ev_dil_water1_nch &
,reaction_path_ev_dil_water2_nch &
,reaction_path_ev_dil_water3_nch &
,begin_element_handler &
,compute_f_df_rmrmt_1_nch & 
,compute_f_df_rmrmt_2_nch &
,pcdata &
,read_backup_nch &
,write_backup_nch
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!% Type pointer to nodal chemistry object 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
type, public:: t_pnodalchemistry

type(t_nodalchemistry), pointer:: ptr  

end type
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!% Type nodal chemistry 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
type, public::t_nodalchemistry
> 
private                         ::
>  
character(len=100)              :: name          !> Name of the nodal chemistry object 

type (t_chemicalsystem), pointer:: pchemsys      !> Pointer to chemical system object 

logical                         :: islockchemsys !> If true the chemical system object 
>                                                 !> was allocated for the nodal chemistry object 

real*8, pointer, dimension(:)   :: c             !> Molality vector [numsp]

real*8, pointer, dimension(:)   :: facc          !> Factor for multiply the c vector [numsp]

real*8, pointer, dimension(:)   :: g             !> Activity coefficients vector [numsp]

real*8, pointer, dimension(:)   :: sktrk         !> Changes due kinetic reactions [numsp]

real*8, pointer, dimension(:)   :: setre         !> Changes due equilibrium reactions [numsp]
> 
real*8, pointer, dimension(:,:) :: txoh          !> Total sites for adsorption [mol]

real*8, pointer, dimension(:,:) :: txoh0         !> Total sites for adsorption [mol]

real*8, pointer, dimension(:,:) :: capint        !> Internal capacitance for adsorption models [C mol-1]

real*8, pointer, dimension(:,:) :: capext        !> External capacitance for adsorption models [C mol-1]

real*8, pointer, dimension(:,:) :: spsurfarea    !> Specific surface area [m2 kg-1] ???? of the water?

integer, pointer, dimension(:)  :: idtxohmin     !> Index of txoh and mineral species [numsurf]

real*8, pointer, dimension(:)   :: alpha         !> Area of mineral species in m2 [numsp] 

real*8, pointer, dimension(:,:) :: alpha0        !> Initial concentration and area of mineral species [2,numsp]
> 
real*8                          :: temp          !> Temperature (in celcius) 

real*8                          :: ionstr        !> Ionic strength [M]

real*8                          :: volnch        !> total volume of the nodal chemistry object [m3]

real*8                          :: volgas        !> Gas volume associated to nodal chemistry object [m3]

real*8                          :: omgwfree      !> Mass of free water [kg]

real*8                          :: pgas          !> Gas preassure [MPa]

real*8                          :: pliq          !> Liquid preassure [MPa]

real*8                          :: liqdens       !> Liquid density [kg/m3]

logical                         :: ispgasconst   !> If true then the gas preassure is imposed 
> 
integer                         :: hashcompz     !> Hash index of the chemical system (identify the components definition zone) 

integer                         :: nchemiter     !> Number of iterations during chemical speciation 

logical                         :: isderstored   !> If true, the derivatives are stored
> 
real*8, pointer, dimension(:,:) :: dc            !> First derivatives of the concentrations 

real*8, pointer, dimension(:,:) :: d2c           !> Second derivatives of the concentrations

real*8, pointer, dimension(:,:) :: dsktrk        !> Derivatives of the kinetic changes

end type t_nodalchemistry
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface add_
> 
module procedure add_txoh_nch
module procedure add_species_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
> 
module procedure create_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
> 
module procedure destroy_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface check_compzone_
> 
module procedure check_compzone_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface check_and_switch_base_
> 
module procedure check_and_switch_base_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_r_sia_
> 
module procedure update_r_sia_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_derivatives_
> 
module procedure update_derivatives_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_kinetic_
> 
module procedure update_kinetic_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_comp_zone_index_
> 
module procedure get_comp_zone_index_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_f_df_rmrmt_
> 
module procedure compute_f_df_rmrmt_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_iumob_
> 
module procedure get_iumob_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_dusktrk_
> 
module procedure get_dusktrk_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_dumob_
> 
module procedure get_dumob_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_duads_
> 
module procedure get_duads_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_chem_info_
> 
module procedure get_chem_info_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface get_from_setre_
> 
module procedure get_from_setre_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_density_ith_
> 
module procedure compute_density_ith_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface equilibrate_
> 
module procedure equilibrate_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface reaction_path_
> 
module procedure reaction_path_nch
module procedure reaction_path_kinetic_nch 
module procedure reaction_path_temperature_nch
module procedure reaction_path_ev_dil_water1_nch
module procedure reaction_path_ev_dil_water2_nch
module procedure reaction_path_ev_dil_water3_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface operator (+)
> 
module procedure mix_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface operator (*)
> 
module procedure scalarbynodalchem_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_from_cpri_
> 
module procedure set_from_cpri_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_from_setre_
> 
module procedure set_from_setre_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_from_u_
> 
module procedure set_from_u_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_from_uaq_
> 
module procedure set_from_uaq_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_from_solution_type_
> 
module procedure set_from_solution_type_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_from_c_
> 
module procedure set_from_c_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface read_xml_
> 
module procedure read_xml_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface read_txt_
> 
module procedure read_txt_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
> 
module procedure set_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
> 
module procedure write_nch_info_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_sat_
> 
module procedure write_sat_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_min_area_
> 
module procedure write_min_area1_nch
module procedure write_min_area2_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_sps_
> 
module procedure write_sps_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_tot_
> 
module procedure write_tot_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_volfracmin_
> 
module procedure write_volfracmin_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_chemical_parameters_
> 
module procedure write_chemical_parameters_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_
> 
module procedure update_yourself_nch
module procedure update_porosity_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_and_check_
> 
module procedure update_from_delta_cpri_and_check_convergence_nch
module procedure update_r_and_check_convergence_sia_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface reset_
> 
module procedure reset_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_backup_
> 
module procedure write_backup_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface read_backup_
> 
module procedure read_backup_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_iswcompbal_
> 
module procedure set_iswcompbal_nch 
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment (=)
> 
module procedure copy_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_min_area_ 
> 
module procedure update_min_area_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_txoh_ 
> 
module procedure update_txoh_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%----------------Private Services------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
interface read_xml_loc_
> 
module procedure read_xml_loc_nch
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface begin_element_handler
> 
module procedure begin_element_handler
> 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface pcdata
> 
module procedure pcdata
> 
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
subroutine create_nch &
>  (this)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Create the nodal chemistry object. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout) :: this !> Type nodal chemistry variable. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
> 
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
!% 1) Nullify the chemical system 
!%------------------------------------------------------------
this%pchemsys => null ()
!%------------------------------------------------------------
!% 2) Nullify another pointers 
!%------------------------------------------------------------
this%c => null ()
this%facc => null ()
this%g => null ()
this%txoh => null ()
this%txoh0 => null ()
this%capint => null ()
this%capext => null ()
this%spsurfarea => null ()
this%idtxohmin => null ()
this%alpha => null ()
this%alpha0 => null ()
this%dc => null ()
this%sktrk => null ()
this%setre => null ()
this%dsktrk => null ()
this%d2c => null ()
!%------------------------------------------------------------
!% 3) Initializing variables
!%------------------------------------------------------------
this%name=' '
this%temp = 0.0d0
this%ionstr = 0.0d0
this%hashcompz = 0
this%isderstored=.false.
this%omgwfree=1.0d0
this%volnch=1.0d0
this%volgas=0.0d0
this%pgas=0.0d0
this%pliq=0.0d0
this%liqdens=1.0d3
this%islockchemsys=.false. 
this%nchemiter=0 
this%ispgasconst=.false. 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_nch &
>  (this)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Destroy the nodal chemistry object. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout) :: this !> Type nodal chemistry variable. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
real*8, parameter            :: &
> r0=0.0d0 
integer, parameter           :: &
> i0=0   
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
!% Destroy the chemical system 
!%------------------------------------------------------------
if (this%islockchemsys) then
> call destroy_ (this%pchemsys)
> deallocate (this%pchemsys)
end if    
this%pchemsys => null ()
this%islockchemsys=.false.
!%-----------------------------------------------------------
!% Deallocate pointers 
!%-----------------------------------------------------------
call check_pointer_ (this%c,1,.false.)
call check_pointer_ (this%facc,1,.false.)
call check_pointer_ (this%g,1,.false.)
call check_pointer_ (this%alpha,1,.false.)
call check_pointer_ (this%alpha0,1,1,.false.)
call check_pointer_ (this%sktrk,1,.false.)
call check_pointer_ (this%setre,1,.false.)
call check_pointer_ (this%txoh,1,1,.false.)
call check_pointer_ (this%txoh0,1,1,.false.)
call check_pointer_ (this%capint,1,1,.false.)
call check_pointer_ (this%capext,1,1,.false.)
call check_pointer_ (this%spsurfarea,1,1,.false.)
call check_pointer_ (this%idtxohmin,1,.false.)
!%-----------------------------------------------------------
!% If the derivatives are stored, then deallocate pointers 
!%-----------------------------------------------------------
if (this%isderstored) then
> call check_pointer_ (this%dc,1,1,.false.)
> call check_pointer_ (this%d2c,1,1,.false.)
> call check_pointer_ (this%dsktrk,1,1,.false.)
end if
!%-----------------------------------------------------------
this%isderstored=.false.
!%-----------------------------------------------------------
this%temp=r0
this%hashcompz=i0
this%name=' '
this%omgwfree=r0
this%volnch=r0
this%pgas=r0 
this%pliq=r0
this%liqdens=r0 
this%volgas=r0
this%ionstr=r0
this%nchemiter=i0 
this%ispgasconst=.false. 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine reset_nch &
>  (this)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Reset nodal chemistry object to zero. 
!
!>   $Arguments:
!
type(t_nodalchemistry), intent(inout):: this !> Type nodal chemistry variable. 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
> 
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
!% Check if the chemical system is associated 
!%------------------------------------------------------------
if (associated(this%pchemsys)) then
> this%c=0.0d0
> this%facc=1.0d0
> this%g=1.0d0
> this%txoh=0.0d0
> this%txoh0=0.0d0
> this%capint=0.0d0
> this%capext=0.0d0
> this%spsurfarea=0.0d0
> this%idtxohmin=0
!%------------------------------------------------------------
!% Only reset alpha, because alpha0 is constant 
!%------------------------------------------------------------
> this%alpha=0.0d0
> this%alpha0=0.0d0
> this%sktrk=0.0d0
> this%setre=0.0d0
> if (this%isderstored) then
>   this%dc=0.0d0
>   this%dsktrk=0.0d0
>   this%d2c=0.0d0
> else 
>   this%dc => null ()
>   this%dsktrk => null ()
>   this%d2c => null ()
> end if 
else 
> this%c => null ()
> this%facc => null ()
> this%g => null ()
> this%txoh => null ()
> this%txoh0 => null ()
> this%capint => null ()
> this%capext => null ()
> this%spsurfarea => null ()
> this%idtxohmin => null ()
> this%alpha => null ()
> this%alpha0 => null ()
> this%dc => null ()
> this%sktrk => null ()
> this%setre => null ()
> this%dsktrk => null ()
> this%d2c => null ()
end if 
!%------------------------------------------------------------
!% 3) Zeroing variables
!%------------------------------------------------------------
this%name=' '
this%temp = 0.0d0
this%ionstr = 0.0d0
this%hashcompz = 0
this%isderstored = .false.
this%omgwfree = 1.0d0
this%volnch = 1.0d0
this%volgas = 0.0d0
this%pgas = 0.0d0
this%pliq = 0.0d0
this%liqdens = 1.0d3
this%nchemiter = 0 
this%ispgasconst=.false. 
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_compzone_nch &
>  (nchk1, &
>   iserror, &
>   nchk)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Check components zone associated to nodal chemistry
!>   If optional variable (nchk) is present, then
!>   nodalchem in k is updated according new components zone
!>   evaluated for nodalchem in k+1.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)            :: nchk1   !> Type nodal chemistry variable. 

logical, intent(out)                             :: iserror !> iserror=true, then there was an error.

type(t_nodalchemistry), intent(inout), optional  :: nchk    !> Type nodal chemistry variable.  
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!------------------------------------------------------------------------- 
integer                                          :: &
> hashcompzk
character(len=100)                               :: &
> msg
integer                                          :: &
> ndim
logical                                          :: &
> havenchk
!-------------------------------------------------------------------------
!
!>   $code
!

msg=''
iserror=.false.
!%------------------------------------------------------------
!%Check optional arguments 
!%------------------------------------------------------------
havenchk=present(nchk)
!%------------------------------------------------------------
!% Check if the chemical system is associated 
!%------------------------------------------------------------
if (.not.associated(nchk1%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
!% 
!%------------------------------------------------------------
hashcompzk=nchk1%hashcompz
!%------------------------------------------------------------
ndim=size(nchk1%c)
call get_hashcompz_ (nchk1%pchemsys,nchk1%hashcompz,nchk1%c,nchk1%g,ndim,iserror)
if (iserror) then
> msg='Error when calling get_hashcompz_'
> goto 10
end if
!%-------------------------------------------------------------
!% 
!%-------------------------------------------------------------
if (hashcompzk/=nchk1%hashcompz) then
> if (havenchk) then
>   if (.not.associated(nchk%pchemsys)) then
>      msg='Error, not associated chemical system'
>      goto 10
>   end if
>   call set_ (nchk,iserror,hashcompz=nchk1%hashcompz) 
>   if (iserror) then
>     msg='Error when calling set_'
>    goto 10
>   end if
> end if
end if
!%------------------------------------------------------------
return
> 
10 continue 
print *,'************************'
print *,'Nodal Chemistry:'
print *,'Name:',nchk1%name
print *,'Service: check_compzone_'
print *,msg
print *,'************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_nch &
>  (this, &
>   namefile, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Read and set the nodal chemistry object from xml file. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)      :: this      !> Type nodal chemistry variable. 

character(len=*), intent(in)               :: namefile  !> Name and path of xml file. 

logical, intent(out)                       :: iserror   !> iserror=true, then there was an error.
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer:: &
>  iostat
type(xml_t)::fxml
character(len=100)          :: &
> name
character(len=100)          :: &
> msg
type(dictionary_t)          :: &
> attributes 
logical                     :: &
> isconvergence
!-------------------------------------------------------------------------
!
!>   $code
!
msg=''
iserror=.false.
!%---------------------------------------------------------------
!> Open the xml unit 
!%---------------------------------------------------------------
call open_xmlfile(namefile, fxml, iostat)
!%---------------------------------------------------------------
if (iostat/=0) then
> msg='Error when open file'
> call add_ (msg,namefile)
> goto 10
end if
!%----------------------------------------------------------------
!%----------------------------------------------------------------
!%----------------------------------------------------------------
print *,'=======> Reading Nodal Chemistry object'
!%----------------------------------------------------------------
call xml_parse(fxml, &
>         begin_element_handler = begin_element_handler, &
>         pcdata_chunk_handler  = pcdata, &
>         verbose = .false.)
!%--------------------------------------------------------
!> End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%---------------------------------------------------------------
call read_xml_loc_ (name,attributes,this,iserror=iserror,isconvergence=isconvergence)
if (iserror) goto 10 
if (.not.isconvergence) then
> msg='Convergence problems when specia water type' 
> goto 10
end if 
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
print *,'=======> Reading Nodal Chemistry object finished'
!%---------------------------------------------------------------
return
> 
10 continue 
print *,'******************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: read_xml_'
print *, msg
print *,'******************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_txt_nch &
>   (this, &
>    namefile, &
>    namethdb, &
>    namekindb, &
>    filebase, &
>    ioptxhl, &
>    iserror, &
	itypechemsys, &
>    iswriteinfo)

implicit none
!-------------------------------------------------------------------------
!
!>   $Description: !> Read and set the nodal chemistry object from ascII file . 
!
!>   $Arguments:
!
> 
type (t_nodalchemistry), intent(inout) :: this         !> Type nodal chemistry variable. 

character(len=*), intent(in)           :: namefile     !> Name of file (similar a root_che.inp of the retraso).

character(len=*), intent(in)           :: namethdb     !> Name of thermodynamic data base (at the moment master25.dat).

character(len=*), intent(in)           :: namekindb    !> Name of kinetic data base. 

character(len=*), intent(in)           :: filebase     !> Path of files. 

integer, intent(in)                    :: ioptxhl      !> If ioptxhl=2 Mass fraction will be computed (in this case, molecular weight are written fo the thermodynamic data base). 

logical, intent(out)                   :: iserror      !> iserror=true, then there was an error.

integer, intent(in), optional          :: itypechemsys !> Specializaton of the chemical system 
>                                                       !> itypechemsys=1 (Classic Chemical System)
													   !> itypechemsys=2 (Retraso Chemical System)

logical, intent(in), optional          :: iswriteinfo  !> If true all speciation calculations will be stored in a output unit. 
>  
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)             :: &
> msg, &
> name, &
> inamecomp, &
> iconstr, &
> label, &
> inamesps 
integer                        :: &
> lbase, &
> ncomp, &
> iicon, &
> ilabel, &
> ilabel1, &
> numsp, &
> numsites, &
> numsurf, &
> itypechemsys1, &
> isps, &
> nsurf, &
> i
logical                        :: &
> isconvergence, &
> isbe, &
> haveiswriteinfo, &
> haveitypechemsys, &
> isequilibrate   
real*8                         :: &
> temp, &
> icguess, &
> ictot, &
> cputime, &
> value1, &
> value2, &
> value3, &
> value4 
integer, parameter             :: &
> ndim=100 
real*8, pointer                :: &
> cguess(:) => null (), &
> ctot(:) => null (), &
> vector1(:) => null (), &
> vector2(:) => null (), &
> vector3(:) => null (), &
> vector4(:) => null ()
integer, pointer               :: &
> icon(:) => null ()
character(len=100), pointer    :: &
> namecomp(:) => null (), &
> unit(:) => null (), &
> constraint(:) => null ()
character(len=100), parameter  :: &
> label1='DEFINITION OF THE GEOCHEMICAL SYSTEM', &
> label2='INITIAL AND BOUNDARY WATER TYPES', &
> label3='INITIAL MINERAL ZONES', &
> label4='INITIAL SURFACE ADSORPTION ZONES'
real*8, parameter              :: &
> r1=1.0d0
integer, parameter             :: &
> iunit=1   
!-------------------------------------------------------------------------
!
!>   $code
!

!%----------------------------------------------------------
iserror=.false.
msg='' 
!%-----------------------------------------------------------
!% Check optional arguments 
!%-----------------------------------------------------------
haveiswriteinfo=present(iswriteinfo)
haveitypechemsys=present(itypechemsys)
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%----------------------------------------------------------- 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
print *,'=======> Reading Nodal Chemistry object'
!%-----------------------------------------------------------
!% Read the chemical system object 
!%-----------------------------------------------------------
name=namefile
call lastletter_ (lbase,filebase)
name=filebase(1:lbase)//name
open(unit=iunit,file=name,status='old',err=40)
!%-----------------------------------------------------------
!% Find the Title label 
!%-----------------------------------------------------------
call find_label_ (label1,iunit,isbe)
if (.not.isbe) then
>  msg='Error, not found the label DEFINITION OF THE GEOCHEMICAL SYSTEM'
>  goto 10
end if
close(unit=iunit)
!%-------------------------------------------------------------
!% Create, read and set the chemical system 
!%-------------------------------------------------------------
if (.not.associated(this%pchemsys)) then 
>   this%islockchemsys=.true.
!%-------------------------------------------------------------
!% Determine the chemical system specialization 
!%-------------------------------------------------------------
>   if (haveitypechemsys) then
>       itypechemsys1=itypechemsys
>   else
>       itypechemsys1=1
>   end if
!%-------------------------------------------------------------
!% Allocate, create and read the chemical system 
!%-------------------------------------------------------------
>   allocate(this%pchemsys)
>   call create_ (this%pchemsys)
>   call read_txt_ (this%pchemsys,namefile,namethdb,namekindb, &
>                   filebase,itypechemsys1,ioptxhl,iserror)
>   if (iserror) then
>      msg='Error when calling read_txt_'
>      goto 10
>   end if
!%-------------------------------------------------------------
!% Set if the Newton-Raphson information must be written 
!%-------------------------------------------------------------
>   if(haveiswriteinfo) then
>      call set_iswriteinfo_ (this%pchemsys,iswriteinfo,iserror)
>      if (iserror) goto 10 
>   end if 
!%-------------------------------------------------------------
!% Create, read and set the nodal chemistry object 
!%-------------------------------------------------------------
>   call get_chem_info_ (this%pchemsys,iserror,numsp=numsp,numsites=numsites, &
>                     numsurf=numsurf)
>   if (iserror) goto 10                     
!%-------------------------------------------------------------
!% Allocate attributtes in nodal chemistry object 
!%-------------------------------------------------------------
>   call check_pointer_ (this%c,numsp,.true.)
>   call check_pointer_ (this%facc,numsp,.true.)
>   call check_pointer_ (this%g,numsp,.true.)
>   call check_pointer_ (this%alpha,numsp,.true.)
>   call check_pointer_ (this%alpha0,2,numsp,.true.)
>   call check_pointer_ (this%sktrk,numsp,.true.)
>   call check_pointer_ (this%setre,numsp,.true.)
>   call check_pointer_ (this%txoh,numsites,numsurf,.true.)
>   call check_pointer_ (this%txoh0,numsites,numsurf,.true.)
>   call check_pointer_ (this%capint,numsites,numsurf,.true.)
>   call check_pointer_ (this%capext,numsites,numsurf,.true.)
>   call check_pointer_ (this%spsurfarea,numsites,numsurf,.true.)
>   call check_pointer_ (this%idtxohmin,numsurf,.true.)
>   this%g=r1
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
open(unit=iunit,file=name,status='old')
!%-----------------------------------------------------------
!% Find label corresponding to waters definition 
!%-----------------------------------------------------------
call find_label_ (label2,iunit,isbe)
if (.not.isbe) then
>  msg='Error, not found the label INITIAL AND BOUNDARY WATER TYPES'
>  goto 10
end if 
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (namecomp,ndim,.true.)
call check_pointer_ (icon,ndim,.true.)
call check_pointer_ (unit,ndim,.true.)
call check_pointer_ (cguess,ndim,.true.)
call check_pointer_ (ctot,ndim,.true.)
call check_pointer_ (constraint,ndim,.true.)
!%-----------------------------------------------------------
!% Initialice local variables 
!%-----------------------------------------------------------
unit='m'
ncomp=0
this%name='Water'
this%isderstored=.true. 
!%-----------------------------------------------------------
!% Read number of waters 
!%-----------------------------------------------------------
read (iunit,*) 
!%-----------------------------------------------------------
!% Read water definition 
!%-----------------------------------------------------------
read (iunit,*,err=30) ilabel,temp
read (iunit,*)
do 
>     read (iunit,*) inamecomp,iicon,icguess,ictot,iconstr
>     if (inamecomp=='*'.or.ncomp>=ndim) exit 
>     ncomp=ncomp+1
>     namecomp(ncomp)=inamecomp 
>     icon(ncomp)=iicon
>     cguess(ncomp)=icguess
>     ctot(ncomp)=ictot
>     constraint(ncomp)=iconstr 
end do 
!%-----------------------------------------------------------
!% Set the nodal chemistry according water definition 
!%-----------------------------------------------------------
print *,'=======> Starting speciation'
call set_from_solution_type_ &
>  (this, &
>   namecomp(1:ncomp), &
>   ctot(1:ncomp), &
>   cguess(1:ncomp), &
>   unit(1:ncomp), &
>   icon(1:ncomp), &
>   constraint(1:ncomp), &
>   ncomp, &
>   temp, &
>   isconvergence, &
>   cputime, &
>   iserror)
if (iserror) goto 20 
if (.not.isconvergence) then 
>   msg='Convergence problems when calling set_from_solution_'
>   goto 50 
end if 
print *,'===> CPU time consumed:',cputime,' [s]'
print *,'===> Number of iterations:',this%nchemiter
print *,'=======> Finishing speciation'
!%-----------------------------------------------------------
!% Find label corresponding to waters definition 
!%-----------------------------------------------------------
call find_label_ (label3,iunit,isbe)
if (isbe) then
!%-----------------------------------------------------------
!% Read the initial mineral zones 
!%-----------------------------------------------------------
read (iunit,*) ilabel
if (ilabel>0) then
> write(6,*) '==========> Reading mineral species' 
> isequilibrate=.false.
> isps=0 
> read (iunit,*)
> read (iunit,*)
> do 
>   isps=isps+1
>   read (iunit,*,err=60) inamesps,value1,value2
>   if (inamesps=='*') then
>    isps=isps-1
>    exit
>   end if 
>   write(6,*) inamesps
>   call add_ (this,inamesps,value1,'m3/m3rock',isequilibrate,iserror,alpha=value2,isconvergence=isconvergence)
>   if (iserror) goto 20
>   if (.not.isconvergence) then 
>     msg='Convergence problems during speciation when call add_ service for species'
>     call add_ (msg,inamesps)
>     goto 10 
>   end if 
> end do
> write(6,*) '==========> Reading mineral species finished'
end if
end if
!%-----------------------------------------------------------
!% Find label corresponding to surface zones
!%-----------------------------------------------------------
call find_label_ (label4,iunit,isbe)
if (isbe) then 
!%-----------------------------------------------------------
!% Get the total numbers of surfaces 
!%-----------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,numsurf=nsurf)
if (iserror) goto 20 
!%-----------------------------------------------------------
!% Read the initial output 
!%-----------------------------------------------------------
read (iunit,*) ilabel 
if (ilabel>0) then
>  call check_pointer_ (vector1,1,.true.)
>  call check_pointer_ (vector2,1,.true.)
>  call check_pointer_ (vector3,1,.true.)
>  call check_pointer_ (vector4,1,.true.)
>  write(6,*) '==========> Reading surface information'
>  do i=1,nsurf 
>     isequilibrate=.false. 
>     read (iunit,*,err=70) ilabel,ilabel1
>     if (ilabel1/=0) isequilibrate=.true. 
>     name='surface'
>     call add_ (name,i) 
>     read (iunit,*) inamesps,vector1(1)
>     read (iunit,*) inamesps,vector2(1)
>     read (iunit,*) inamesps,vector3(1)
>     read (iunit,*) inamesps,vector4(1)
>     call add_ (this,name,vector4,vector1,vector2,vector3,1,isequilibrate,iserror, &
>                isconvergence=isconvergence)   
>     if (iserror) goto 20
>     if (.not.isconvergence) then 
>       msg='Convergence problems during speciation when call add_ service for surface'
>       call add_ (msg,i)
>       goto 10 
>     end if      
>  end do
>  write(6,*) '==========> Reading surface zones finished' 
end if
end if
!%-----------------------------------------------------------
!%-----------------------------------------------------------
print *,'=======> Finishing reading Nodal Chemistry object'
!%-----------------------------------------------------------
!% Close unit 
!%-----------------------------------------------------------
close (unit=iunit)
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers
!%-----------------------------------------------------------
call check_pointer_ (namecomp,1,.false.)
call check_pointer_ (icon,1,.false.)
call check_pointer_ (unit,1,.false.)
call check_pointer_ (cguess,1,.false.)
call check_pointer_ (ctot,1,.false.)
call check_pointer_ (constraint,1,.false.)
call check_pointer_ (vector1,1,.false.)
call check_pointer_ (vector2,1,.false.)
call check_pointer_ (vector3,1,.false.)
call check_pointer_ (vector4,1,.false.)
if (iserror) goto 10 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
return
> 
10 continue 
print *,'*******************'
print *,'Nodal Chemistry:   '
print *,'Name:',this%name
print *,'Service: read_txt_ '
print *, msg
print *,'*******************'
iserror=.true.
return
30  msg='Error reading water' 
goto 10 
40  msg='Error when open file:'
call add_ (msg,namefile) 
goto 10
50 continue 
print *,'*******************'
print *,'Nodal Chemistry:   '
print *,'Name:',this%name
print *,'Service: read_txt_ '
print *, msg
print *,'*******************'
return 
60  msg='Error reading infomation of the mineral species' 
goto 10 
return 
70  msg='Error reading surface information for surface'
call add_ (msg,i) 
goto 10 
return 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_from_solution_type_nch &
>  (this, &
>   namecomp, &
>   ctot, &
>   cguess, &
>   unit, &
>   icon, &
>   constraint, &
>   ncomp, &
>   temp, &
>   isconvergence, &
>   cputime, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description:  Set the concentrations (and dependent variables) from 
!> the type of solution (e.g. according pH, total concentrations, etc.). 
!> The input number of componentes 
!
!>   $Arguments:
!
> 
type (t_nodalchemistry), intent(inout)         :: this            !> Type nodal chemistry variable. 

integer, intent (in)                           :: ncomp           !> Number of components 

integer, intent (in), dimension(ncomp)         :: icon            !> Flag indicating the type of constraint controlling the solute content which is given under ctot:
>                                                                  !> 0= the concentration of the species is assumed equal to ctot. 
														          !> 1= the concentration of the species is constrained by the total concentration of solute ctot, except for water.
														          !> 2= the concentration of the species is calculated through charge balance. 
														          !> 3= the activity of the species is assumed initially equal to ctot. 
														          !> 4= the concentration of the species is calculated from the equilibrium of the solution with any mineral specified in the constraint vector. 
														          !> 5= the activity of the species is calculated from the equilibrium with a partial pressure of the gas specified in the constraint vector. 

real*8, intent (in), dimension(ncomp)          :: cguess          !> Initial concentrations [ncomp]

real*8, intent (in), dimension(ncomp)          :: ctot            !> Total concentrations vector [ncomp]

character(len=*),intent (in), dimension(ncomp) :: namecomp        !> Name of components [ncomp]

character(len=*),intent (in), dimension(ncomp) :: unit            !> Unit of the total and initial concentrations. 

character(len=*),intent (in), dimension(ncomp) :: constraint      !> Name of mineral or gas used to constraint the concentration of the species. 

real*8, intent(in)                             :: temp            !> Temperature (in celcius).

logical, intent(out)                           :: iserror         !> iserror=true, then there was an error.

real*8, intent(out)                            :: cputime         !> CPU spent during speciation calculations

logical, intent(out)                           :: isconvergence   !> isconvergence=true, there was convergence in chemical speciation. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
> 
integer                     :: &
> i, &
> nmin, &
> ntxoh, &
> nsurf, &
> nsp
real*8                      :: &
> conc
logical                     :: &
> isrepeated, &
> isnonconvergence
character(len=100), pointer :: &
> namemin(:) => null ()
real*8, pointer             :: &
> si(:) => null (), &
> ctotloc(:) => null (), &
> cguessloc(:) => null ()
character(len=100)          :: &
> msg, &
> namerepeated 
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Check if the chemical system is associated
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
!% Check if the component name is repeated
!%------------------------------------------------------------
if (ncomp>0) then
> call find_repeated_(namecomp,isrepeated,ncomp,namerepeated)
> if (isrepeated) then
>  msg='Error, the component is repeated:'
>  call add_ (msg,namerepeated)
>  goto 10
> end if
end if
!%-----------------------------------------------------------
!% Set the temperature
!%-----------------------------------------------------------
this%temp=temp
call set_ (this,iserror,temp=this%temp)
if (iserror) goto 10 
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (ctotloc,ncomp,.true.)
call check_pointer_ (cguessloc,ncomp,.true.)
!%-----------------------------------------------------------
!% Change units 
!%-----------------------------------------------------------
do i=1,ncomp
> conc=cguess(i)
> call change_chem_unit_(this%pchemsys,conc,namecomp(i),unit(i),'m',iserror)
> cguessloc(i)=conc
!%----------------------
> conc=ctot(i)
> call change_chem_unit_(this%pchemsys,conc,namecomp(i),unit(i),'m',iserror)
> ctotloc(i)=conc
end do
!%-----------------------------------------------------------
!%-----------------------------------------------------------
if (associated(this%txoh)) then
> ntxoh=size(this%txoh,1)
> nsurf=size(this%txoh,2)
else
> ntxoh=0
> nsurf=0
end if
!%-----------------------------------------------------------
!% Specia according type of solution
!%-----------------------------------------------------------
nsp=size(this%c)
!%-----------------------------------------------------------
this%alpha=this%alpha/this%omgwfree
!%-----------------------------------------------------------
if (this%isderstored) then  !> If the derivate are stored 
> call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
>    this%alpha, &
>    nsp, &
>    nmin, &
>    si, &
>    namemin, &
>    ncomp, &
>    namecomp, &
>    icon, &
>    ctotloc, &
>    cguessloc, &
>    constraint, &
>    this%ionstr, &
>    this%txoh, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    ntxoh, &
>    nsurf, &
>    isconvergence, &
>    this%hashcompz, &
>    iserror, &
>    dc=this%dc, &
>    dsktrk=this%dsktrk, &
>    sktrk=this%sktrk, &
	cputime=cputime, &
>    nchemiter=this%nchemiter)
else
>  call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
>    this%alpha, &
>    nsp, &
>    nmin, &
>    si, &
>    namemin, &
>    ncomp, &
>    namecomp, &
>    icon, &
>    ctotloc, &
>    cguessloc, &
>    constraint, &
>    this%ionstr, &
>    this%txoh, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    ntxoh, &
>    nsurf, &
>    isconvergence, &
>    this%hashcompz, &
>    iserror, &
>    sktrk=this%sktrk, &
	cputime=cputime, &
>    nchemiter=this%nchemiter)
end if
!%-----------------------------------------------------------
this%alpha=this%alpha*this%omgwfree
!%-----------------------------------------------------------
if (iserror) then
> msg='Error when calling specia_'
> goto 20
end if
20 continue
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------- 
call check_pointer_ (si,1,.false.)
call check_pointer_ (namemin,1,.false.)
call check_pointer_ (ctotloc,1,.false.)
call check_pointer_ (cguessloc,1,.false.)
if (iserror) goto 10
!%---------------------------------------------------
return
> 
10 continue 
print *,'********************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: set_from_solution_type_'
print *, msg
print*, '********************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_nch &
>  (this, &
>   iserror, &
>   chemsys, &
>   isderstored, &
>   name, &
>   volnch, &
>   omgwfree, &
>   temp, &
>   hashcompz, &
>   volgas, &
>   pgas,&
>   pliq, &
>   liqdens)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: !> Set general attributtes in the nodal chemistry object. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)               :: this        !> Type nodal chemistry variable. 

logical, intent(out)                                :: iserror     !> If true there was a error. 

type(t_chemicalsystem), intent(in),optional,target  :: chemsys     !> Type chemical system variable 

logical, intent(in), optional                       :: isderstored !> isderstored=true, derivatives with respect to primary species are stored. 

integer, intent(in), optional                       :: hashcompz   !> Index of components definition.

real*8, intent(in), optional                        :: temp        !> Temperature (in celcius). 

real*8, intent(in), optional                        :: volnch      !> Volume associated to nodal chemistry object [m3]

real*8, intent(in), optional                        :: omgwfree    !> Mass of free water [kg]

real*8, intent(in), optional                        :: volgas      !> Gas volume associated to nodal chemistry object.

real*8, intent(in), optional                        :: pgas        !> Gas pressure imposed in the nodal chemistry objct [Mpa]

real*8, intent(in), optional                        :: pliq        !> Liquid pressure in the nodal chemistry objct [Mpa]

real*8, intent(in), optional                        :: liqdens     !> Liquid density [kg/m3]

character(len=*),intent(in), optional               :: name        !> Name of the nodal chemistry object. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)                           :: &
> msg
integer, pointer                             :: &
> ivector1(:) => null (), &
> ivector2(:) => null ()
integer                                      :: &
> numsp, &
> numsp1, &
> numsp2, &
> numsites, &
> numsurf
logical                                      :: &
> havechemsys, &
> haveisderstored, &
> havename, &
> havevolnch, &
> haveomgwfree, &
> havetemp, &
> havehashcompz, &
> havevolgas, &
> havepgas, &
> havepliq, &
> haveliqdens, &
> isconvergence 
!-------------------------------------------------------------------------
!
!>   $code
!
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check optional arguments
!%-------------------------------------------------------------
havechemsys=present(chemsys)
haveisderstored=present(isderstored)
havename=present(name)
havevolnch=present(volnch)
haveomgwfree=present(omgwfree)
havetemp=present(temp)
havehashcompz=present(hashcompz)
havevolgas=present(volgas)
havepgas=present(pgas)
havepliq=present(pliq)
haveliqdens=present(liqdens)
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!% Set the chemical system 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (havechemsys.and..not.associated(this%pchemsys,chemsys)) then
!%-------------------------------------------------------------
!% Pointer to new chemical system
!% If the chemical system was allocated previously for the nodal 
!% chemistry object, in first place destroy and allocate the 
!% chemical system. 
!%-------------------------------------------------------------
>   if (this%islockchemsys) then
>      call destroy_ (this%pchemsys)
>      deallocate (this%pchemsys)
>      this%pchemsys => null ()  
>      this%islockchemsys=.false.  
>   end if 
!%-------------------------------------------------------------
!% Pointer to new chemical system object 
!%-------------------------------------------------------------
>   this%pchemsys => chemsys
!%-------------------------------------------------------------
!% If the chemical system is different, then:
!%-------------------------------------------------------------
>    call get_chem_info_ (this%pchemsys,iserror,numsp=numsp, &
	                     numsites=numsites,numsurf=numsurf)
>    if (iserror) goto 10                       
!%-------------------------------------------------------------
!% Allocate nodal chemistry attributtes 
!%-------------------------------------------------------------
>    call check_pointer_ (this%c,numsp,.true.)
	call check_pointer_ (this%facc,numsp,.true.)
>    call check_pointer_ (this%g,numsp,.true.)
>    call check_pointer_ (this%alpha,numsp,.true.)
>    call check_pointer_ (this%alpha0,2,numsp,.true.)
>    call check_pointer_ (this%sktrk,numsp,.true.)
	call check_pointer_ (this%setre,numsp,.true.)
>    call check_pointer_ (this%txoh,numsites,numsurf,.true.)
>    call check_pointer_ (this%capint,numsites,numsurf,.true.)
>    call check_pointer_ (this%capext,numsites,numsurf,.true.)
>    call check_pointer_ (this%spsurfarea,numsites,numsurf,.true.)
	call check_pointer_ (this%txoh0,numsites,numsurf,.true.)
	call check_pointer_ (this%idtxohmin,numsurf,.true.)
!%-------------------------------------------------------------
!% Initialice activity coefficients vector to 1 
!%-------------------------------------------------------------
>    this%g=1.0d0 
	this%facc=1.0d0 
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%Set if the derivatives are stored 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (haveisderstored) then
>  this%isderstored=isderstored
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%Set the name to nodal chemistry object 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (havename) then
>  this%name=name
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%Set the liquid density 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (haveliqdens) then
>  this%liqdens=liqdens
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%Set the mass of free water in nodal chemistry object 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (haveomgwfree.and.associated(this%pchemsys)) then
>  this%omgwfree=omgwfree
>  call get_chem_info_ (this%pchemsys,iserror,nminsp=numsp1,ngassp=numsp2,&
	                   idminsp=ivector1,idgassp=ivector2)
>  if (iserror) goto 20 
>  this%facc=this%omgwfree 
>  if (numsp1>0) then
>     this%facc(ivector1)=this%volnch
>  end if 
>  if (numsp2>0) then
>     this%facc(ivector2)=this%volgas/(rgas*(this%temp+273.15d0))
>  end if 
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!%Set the volume associated to nodal chemistry object [m3rock]
!% In addition update the total values of reactive surface of
!% kinetic minerals, gas volume and total sites fo adsorption 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (havevolnch.and.associated(this%pchemsys)) then
> this%volgas=this%volgas*(volnch/this%volnch)
> this%alpha=this%alpha*(volnch/this%volnch)
> this%alpha0=this%alpha0*(volnch/this%volnch)
> this%txoh=this%txoh*(volnch/this%volnch)
> this%txoh0=this%txoh0*(volnch/this%volnch)
> this%volnch=volnch
> call get_chem_info_ (this%pchemsys,iserror,nminsp=numsp1,ngassp=numsp2,&
	                   idminsp=ivector1,idgassp=ivector2)
> if (iserror) goto 20 
> if (numsp1>0) then
>     this%facc(ivector1)=this%volnch
> end if 
> if (numsp2>0) then
>     this%facc(ivector2)=this%volgas/(rgas*(this%temp+273.15d0))
> end if 
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!> Set volgas [m3gas/m3rock]
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (havevolgas) then
> this%ispgasconst=.false.
> this%volgas=volgas*this%volnch
> call get_chem_info_ (this%pchemsys,iserror,ngassp=numsp1,idgassp=ivector1)
> if (iserror) goto 20 
> if (numsp1>0) then
>     this%facc(ivector1)=this%volgas/(rgas*(this%temp+273.15d0))
> end if 
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!> Set Pgas, then not impose the gas preassure 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (havepgas) then
> this%ispgasconst=.true. 
> this%pgas=pgas 
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!> Set Pliq
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (havepliq) then
>  this%pliq=pliq
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!% Set the temperature
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (havetemp.and.associated(this%pchemsys)) then
>    if (this%temp/=temp) then
	      this%temp=temp 
>          call get_chem_info_ (this%pchemsys,iserror,ngassp=numsp1, &
		                       idgassp=ivector1)
>          if (iserror) goto 20 
>          if (numsp1>0) then
>             this%facc(ivector1)=this%volgas/(rgas*(this%temp+273.15d0))
>          end if 
!%-------------------------------------------------------------
!% Update parameters in the chemical system object  
!%-------------------------------------------------------------
>          call update_ (this%pchemsys,this%temp,iserror)
>          if (iserror) then
>            msg='Error when calling update_'
>          goto 20
>          end if
	end if
end if
!%-------------------------------------------------------------
!%-------------------------------------------------------------
!% Update derivates is the components zone index is different 
!% to preview 
!%-------------------------------------------------------------
!%-------------------------------------------------------------
if (havehashcompz.and.hashcompz/=this%hashcompz) then
>    this%hashcompz=hashcompz
>    call update_derivatives_ (this,iserror)
>    if (iserror) then
>      msg='Error when calling update_derivatives_'
>      goto 20
>    end if
end if
!%-------------------------------------------------------------
20 continue
!%-------------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (ivector1,1,.false.)
call check_pointer_ (ivector2,1,.false.)
if (iserror) goto 10 
!%-------------------------------------------------------------
return
10 continue 
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: set_'
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
subroutine update_yourself_nch &
>  (this, &
>   iserror, &
>   isconvergence)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Update yourself the nodal chemistry.
!> Compute the speciation according total analytical concentrations. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout):: this          !> Type nodal chemistry variable. 

logical, intent(out)                 :: iserror       !> iserror=true, there was a error. 

logical, intent(out)                 :: isconvergence !> isconvergence=true, there was convergence during speciation calculations. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
real*8, pointer                      :: &
> u(:) => null ()
integer, pointer                     :: &
> idpri(:) => null ()
integer                              :: &
> npri, &
> numsp, &
> ntxoh, &
> nsurf
real*8                               :: &
> factoromgw
character(len=100)                   :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=' '
!%------------------------------------------------------------
!% Check if associated the chemical system 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>  msg='Error, not associated chemical system'
>  goto 10
end if
!%------------------------------------------------------------
!% Compute the total analytic concentrations 
!%------------------------------------------------------------
call make_lin_trf_(this%pchemsys,u,this%c,0,iserror)
if (iserror) then
>  msg='Error when calling make_lin_trf_'
>  goto 20
end if
!%-----------------------------------------------------------
call get_chem_info_(this%pchemsys,iserror,idbase=idpri)
if (iserror) then
>  msg='Error when calling get_chem_info_'
>  goto 20
end if
!%------------------------------------------------------------
!% Determine dimmensions 
!%------------------------------------------------------------
numsp=size(this%c)
npri=size(idpri)
ntxoh=size(this%txoh,1)
nsurf=size(this%txoh,2)
!%------------------------------------------------------------
factoromgw=1.0d0
this%alpha=this%alpha/this%omgwfree 
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
call update_ (this%pchemsys,this%temp,iserror)
if (iserror) goto 20 
!%-----------------------------------------------------------
if (this%isderstored) then !> If the derivate are stored 
> 
> call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	.true., &
	this%ispgasconst, &
>    this%ionstr, &
>    this%alpha, &
>    u, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    numsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    0.0d0, &
>    this%hashcompz, &
>    isconvergence, &
>    factoromgw, &
	this%volgas, &
	this%pgas, &
	msg, &
>    iserror, &
>    cguess=this%c(idpri), &
>    dc=this%dc, &
>    sktrk=this%sktrk, &
>    dsktrk=this%dsktrk, &
>    nchemiter=this%nchemiter)
> 
else
> 
> call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
>    .true., &
	this%ispgasconst, &
>    this%ionstr, &
>    this%alpha, &
>    u, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    numsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    0.0d0, &
>    this%hashcompz, &
>    isconvergence, &
>    factoromgw, &
	this%volgas, &
	this%pgas, &
	msg, &
>    iserror, &
>    cguess=this%c(idpri), &
>    sktrk=this%sktrk, &
>    nchemiter=this%nchemiter)
> 
end if
!%------------------------------------------------------------
this%alpha=this%alpha*this%omgwfree
!%------------------------------------------------------------
!% Update the mass of free water 
!%------------------------------------------------------------
this%omgwfree = this%omgwfree * factoromgw
!%------------------------------------------------------------
!% Update reactive surface of kinetic minerals 
!%------------------------------------------------------------
if (isconvergence) then   
>   call update_min_area_ (this,iserror)
end if 
!%------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!%Deallocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (u,1,.false.)
call check_pointer_ (idpri,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
> 
return
> 
10 continue 
print *,'*************************'
print *,'Nodal Chemistry:'
print *,'Name:', this%name
print *,'Service: update_yourself_'
print *, msg
print *,'*************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_derivatives_nch &
>  (this, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Update derivates (only if derivates are stored in the 
!> nodal chemistry)
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout) :: this    !> Type nodal chemistry variable. 

logical, intent(out)                  :: iserror !> iserror=true, then there was an error 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                              :: &
> npri, &
> nsp, &
> ntxoh, &
> nsurf
character(len=100)                   :: &
> msg
real*8                               :: &
> faccap 
type(t_nodalchemistry), pointer      :: &
> pnch => null ()
logical                              :: &
> isanomalous, &
> isupmxitergam
integer, pointer                     :: &
> idpri(:) => null () 
!-------------------------------------------------------------------------
!
!>   $code
!
iserror=.false.
msg=' '
!%------------------------------------------------------------
!% Check if the chemical system is associated 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>  msg='Error, not associated chemical system'
>  goto 10
end if
!%------------------------------------------------------------
!> Update derivates (Only if derivates are strored)
!%------------------------------------------------------------
if (this%isderstored) then
> 
> call get_chem_info_ &
>   (this%pchemsys, &
>    iserror, &
>    numbase=npri, &
>    numsp=nsp, &
>    numsites=ntxoh, &
>    numsurf=nsurf, &
>    idbase=idpri, &
>    hashcompz=this%hashcompz)
!%------------------------------------------------------------
!% Allocate and create a local nodal chemistry object 
!%------------------------------------------------------------ 
> allocate (pnch)
> call create_ (pnch)
!%------------------------------------------------------------
!% Copy in local nodal chemistry object 
!%------------------------------------------------------------
> pnch=this
> pnch%alpha=pnch%alpha/pnch%omgwfree
> pnch%alpha0=pnch%alpha0/pnch%omgwfree
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
call update_ (pnch%pchemsys,pnch%temp,iserror)
if (iserror) goto 20
!%-----------------------------------------------------------
!%  
!%-----------------------------------------------------------
call get_chem_info_ (pnch,iserror,faccap=faccap)
if (iserror) goto 20
!%-----------------------------------------------------------
!% Specia according cpri 
!%-----------------------------------------------------------
> call specia_ &
>   (pnch%pchemsys, &
>    pnch%temp, &
>    pnch%c, &
>    pnch%g, &
	pnch%c, &
>    pnch%alpha, &
>    nsp, &
>    pnch%c(idpri), &
>    npri, &
>    pnch%txoh/pnch%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    ntxoh, &
>    nsurf, &
>    0.0d0, &
>    pnch%volgas, &
>    pnch%ionstr, &
>    pnch%hashcompz, &
	faccap, &
>    isanomalous, &
>    isupmxitergam, &
>    .true., &
>    iserror, &
>    dc=this%dc, &
>    dsktrk=this%dsktrk)
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate and destroy the local nodal chemistry object 
!%------------------------------------------------------------
> call destroy_ (pnch)
> deallocate (pnch)
!%-------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------- 
> call check_pointer_ (idpri,1,.false.)
!%-------------------------------------------------------
> if (iserror) then
>  msg='Error when calling specia_'
>  goto 10
> end if
!%-------------------------------------------------------
!% Check if there was anomalous concentrations 
!%-------------------------------------------------------
> if (isanomalous) then
>  msg='Convergence problems, anomalous concentrations'
>  goto 10
> end if
!%-------------------------------------------------------
!% Check if there was convergence in gamma iterations
!%-------------------------------------------------------
> if (isupmxitergam) then
>  print *,'Convergence problems, number of gamma iterations exeeded'
> end if
!%-------------------------------------------------------
end if
!%-----------------------------------------------------------
> 
return
> 
10 continue 
print *,'***************************'
print *,'Nodal Chemistry:'
print *,'Name:', this%name
print *,'Service: update_derivatives_'
print *, msg
print *,'***************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_kinetic_nch &
>  (this, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Update kinetic term and derivatives with respect to primary species. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout) :: this    !> Type nodal chemistry variable. 

logical, intent(out)                  :: iserror !> iserror=true, there was an error. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                              :: &
> npri, &
> nsp, &
> ntxoh, &
> nsurf
character(len=100)                   :: &
> msg
type(t_nodalchemistry), pointer      :: &
> pnch => null ()
logical                              :: &
> isanomalous, &
> isupmxitergam
real*8, pointer                      :: &
> cold(:) => null ()
real*8                               :: &
> faccap 
integer, pointer                     :: &
> idpri(:) => null () 
real*8, parameter                    :: &
> r0=0.0d0
!-------------------------------------------------------------------------
!
!>   $code
!
iserror=.false.
msg=' '
!%------------------------------------------------------------
!% Check if the chemical system is associated 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>  msg='Error, not associated chemical system'
>  goto 10
end if
!%------------------------------------------------------------
!% Update derivates (Only if derivates are strored)
!%------------------------------------------------------------
call get_chem_info_ &
>   (this%pchemsys, &
>    iserror, &
>    numbase=npri, &
>    numsp=nsp, &
>    numsites=ntxoh, &
>    numsurf=nsurf, &
>    idbase=idpri, &
>    hashcompz=this%hashcompz)
> 
allocate (pnch)
call create_ (pnch)
pnch=this
call check_pointer_ (cold,nsp,.true.)
cold=pnch%c
pnch%alpha=pnch%alpha/pnch%omgwfree
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
call update_ (pnch%pchemsys,pnch%temp,iserror)
if (iserror) goto 20
!%-----------------------------------------------------------
!%  
!%-----------------------------------------------------------
call get_chem_info_ (pnch,iserror,faccap=faccap)
if (iserror) goto 20
!%-----------------------------------------------------------
if (this%isderstored) then
>   call specia_ &
>   (pnch%pchemsys, &
>    pnch%temp, &
>    pnch%c, &
>    pnch%g, &
	cold, &
>    pnch%alpha, &
>    nsp, &
>    pnch%c(idpri), &
>    npri, &
>    pnch%txoh/pnch%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    ntxoh, &
>    nsurf, &
>    r0, &
	pnch%volgas, &
>    pnch%ionstr, &
>    pnch%hashcompz, &
	faccap, &
>    isanomalous, &
>    isupmxitergam, &
>    .true., &
>    iserror, &
>    sktrk=this%sktrk, &
>    dsktrk=this%dsktrk)
else
>   call specia_ &
>   (pnch%pchemsys, &
>    pnch%temp, &
>    pnch%c, &
>    pnch%g, &
	cold, &
>    pnch%alpha, &
>    nsp, &
>    pnch%c(idpri), &
>    npri, &
>    pnch%txoh/pnch%omgwfree, &
>    pnch%capint, &
>    pnch%capext, &
>    pnch%spsurfarea, &
>    ntxoh, &
>    nsurf, &
>    r0, &
	pnch%volgas, &
>    pnch%ionstr, &
>    pnch%hashcompz, &
	faccap, &
>    isanomalous, &
>    isupmxitergam, &
>    .true., &
>    iserror, &
>    sktrk=this%sktrk)
end if 
!%-------------------------------------------------------
20 continue 
call destroy_ (pnch)
deallocate (pnch)
!%-------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------
call check_pointer_ (cold,1,.false.)
call check_pointer_ (idpri,1,.false.)
!%-------------------------------------------------------
if (iserror) then
>  msg='Error when calling specia_'
>  goto 10
end if
!%-------------------------------------------------------
if (isanomalous) then
>  msg='Convergence problems, anomalous concentrations'
>  goto 10
end if
!%-------------------------------------------------------
if (isupmxitergam) then
>  msg='Convergence problems, number of gamma iterations exeeded'
>  goto 10
end if
!%-----------------------------------------------------------
> 
return
> 
10 continue 
print *,'***************************'
print *,'Nodal Chemistry:'
print *,'Name:', this%name
print *,'Service: update_derivatives_'
print *, msg
print *,'***************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_porosity_nch &
>   (nchk1, &
>    porosity, &
>    nchk, &
>    iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Update the porosity. 
!> 
!> por (k+1)=por(k)+delta_por
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)   :: nchk1     !> Type nodal chemistry variable (in k+1). 

type(t_nodalchemistry), intent(in)   :: nchk      !> Type nodal chemistry variable (in k). 

real*8, intent(inout)                :: porosity  !> Porosity 

logical, intent(out)                 :: iserror   !> If true there was a error.  
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                              :: ndim

character(len=100)                   :: msg 
!-------------------------------------------------------------------------
!
!>   $code
!
> 
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check if the chemical system are associated 
!%------------------------------------------------------------
if (.not.associated(nchk1%pchemsys) &
>                  .or. &
>    .not.associated(nchk%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
ndim=size(nchk1%c)
call update_ &
>   (nchk1%pchemsys, &
>    porosity, &
>    nchk1%c, &
>    nchk%c, &
>    ndim, &
>    nchk1%omgwfree, &
>    nchk%omgwfree, &
>    iserror)
if (iserror) goto 10
!%------------------------------------------------------------
return
> 
10 continue 
print *, '*****************************'
print *, 'Nodal Chemistry:'
print *, 'Name:',nchk1%name
print *, 'Service: update_ (porosity)'
print *,  msg
print *, '*****************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine check_and_switch_base_nch &
>   (this, &
>    iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Check and change the chemical base defined in the chemical 
!> system 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout) :: this      !> Type nodal chemistry variable. 

logical, intent(out)                  :: iserror   !> iserror=true, then there was an error.  
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                              :: &
> nsp, &
> naqprisp, &
> nadsprisp, &
> isps1, &
> isps2
character(len=100)                   :: &
> msg 
character(len=100), pointer          :: &
> nameaqprisp(:) => null (), &
> nameadsprisp(:) => null ()
!-------------------------------------------------------------------------
!
!>   $code
!
> 
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check if the chemical system is associated 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
!% Get the new name of the primary and sorption primary 
!% species 
!%------------------------------------------------------------
nsp=size(this%c)
call get_new_chemical_base_ &
>   (this%pchemsys, &
>    nameaqprisp, &
	nameadsprisp, &
>    this%c, &
	nsp, &
>    iserror)
!%------------------------------------------------------------
!% Change the chemical base in the chemical system  
!%------------------------------------------------------------
naqprisp=size(nameaqprisp)
nadsprisp=0
call switch_base_ &
>   (this%pchemsys, &
>    naqprisp, &
	nadsprisp, &
>    nameaqprisp, &
	nameadsprisp, &
>    iserror)
if (iserror) goto 20 
!%------------------------------------------------------------
20 continue
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (nameaqprisp,1,.false.)
call check_pointer_ (nameadsprisp,1,.false.)
if (iserror) goto 10 
!%------------------------------------------------------------
return
> 
10 continue 
print *, '*****************************************'
print *, 'Nodal Chemistry:'
print *, 'Name:',this%name
print *, 'Service: check_and_change_chemical_base_'
print *,  msg
print *, '*****************************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_comp_zone_index_nch &
> (this, &
>  ithcompz, &
>  iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Return the index of the component defintion corresponding to nodal chemistry object. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)  :: this     !> Type nodal chemistry variable. 

integer, intent(out)                :: ithcompz !> Index of components zone definition. 

logical, intent(out)                :: iserror  !> iserror=true, there was an error. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond: 
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)                  :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Check if the cheical system is associated 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
ithcompz=this%hashcompz
!%------------------------------------------------------------
return
10 continue 
print *, '*****************************'
print *, 'Nodal Chemistry:'
print *, 'Name:',this%name
print *, 'Service: get_comp_zone_index_'
print *,  msg
print *, '*****************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_from_u_nch &
>  (this, &
>   u, &
>   dtime, &
>   isconvergence, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Set the concentrations (and dependent variables) according total concentrations in all phases. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout) :: this          !> Type nodal chemistry variable. 

real*8, intent(in), dimension(:)      :: u             !> Total concentrations in all phases [ncomp]. 

real*8, intent(in)                    :: dtime         !> Time increment [s]

logical, intent(out)                  :: iserror       !> iserror=true, there was a error. 

logical, intent(out)                  :: isconvergence !> iscovergence=true, there was convergence during speciation calculations. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer,pointer       :: &
> idpri(:) => null ()
integer               :: &
> npri, &
> ndim, &
> numsp, &
> ntxoh, &
> nsurf
logical               :: &
> isconv
real*8                :: &
> factoromgw
character(len=100)    :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
msg=' '
iserror=.false.
!%---------------------------------------------------------
!% Check if the chemical system is associated 
!%---------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>  msg='Error, not associated chemical system'
>  goto 10
end if
!%---------------------------------------------------------
!% Initialice variables
!%---------------------------------------------------------
factoromgw=1.0d0 
!%---------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,idbase=idpri)
if (iserror) then
> msg='Error when calling get_chem_info_'
> goto 20
end if
!%----------------------------------------------------------
!% Get the number of primary species
!%----------------------------------------------------------
npri=size(idpri)
!%----------------------------------------------------------
ndim=size(u)
ntxoh=size(this%txoh,1)
nsurf=size(this%txoh,2)
numsp=size(this%c)
!%------------------------------------------------------------
!% Check the number of components 
!%----------------------------------------------------------
ndim=size(u)
if (ndim/=npri) then
> msg='Error in number of components'
> goto 20
end if
!%----------------------------------------------------------
!% Get the remainder dimensions 
!%----------------------------------------------------------
ntxoh=size(this%txoh,1)
nsurf=size(this%txoh,2)
numsp=size(this%c)
!%----------------------------------------------------------
this%alpha=this%alpha/this%omgwfree
!%----------------------------------------------------------
!% Update the chemical system with the temperature
!%----------------------------------------------------------
call update_ (this%pchemsys,this%temp,iserror)
if (iserror) goto 20
!%----------------------------------------------------------
!% Specia from analytic concentrations 
!%----------------------------------------------------------
if (this%isderstored) then !> If the derivatives are stored 
> call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	.true., &
	this%ispgasconst, &
>    this%ionstr, &
>    this%alpha, &
>    u, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    numsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    dtime, &
>    this%hashcompz, &
>    isconv, &
>    factoromgw, &
	this%volgas, &
	this%pgas, &
	msg, &
>    iserror, &
>    cguess=this%c(idpri), &
>    dc=this%dc, &
>    sktrk=this%sktrk, &
>    dsktrk=this%dsktrk, &
>    nchemiter=this%nchemiter)
else
> call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	.true., &
	this%ispgasconst, &
>    this%ionstr, &
>    this%alpha, &
>    u, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    numsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    dtime, &
>    this%hashcompz, &
>    isconv, &
>    factoromgw, &
	this%volgas, &
	this%pgas, &
	msg, &
>    iserror, &
>    cguess=this%c(idpri), &
>    sktrk=this%sktrk, &
>    nchemiter=this%nchemiter)
end if
!%-----------------------------------------------------------
!% Correct the mass of free water 
!%-----------------------------------------------------------
this%alpha=this%alpha*this%omgwfree
!%-----------------------------------------------------------
!% Update the mass of free water 
!%-----------------------------------------------------------
this%omgwfree = this%omgwfree * factoromgw 
!%-----------------------------------------------------------
!% If there was convergence 
!%-----------------------------------------------------------
if (isconv) then  
>   isconvergence=.true. 
!%-----------------------------------------------------------
!% Update the reactive surface of kinetic minerals 
!%-----------------------------------------------------------
>   call update_min_area_ (this,iserror)
>   if (iserror) goto 20
else
>   isconvergence=.false.
end if
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (idpri,1,.false.)
!%------------------------------------------------------------
if (iserror) goto 10
!%------------------------------------------------------------
return
> 
> 
10 continue 
print *,'********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service:set_from_u_'
print *, msg
print *,'********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_iswcompbal_nch &
>   (this, &
>    iswcompbal, &
>    iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Set if the water component mass balance must be evaluated dring chemical speciation. 
!
!>   $Arguments:
!
> 
type (t_nodalchemistry), intent(inout)      :: this       !> Type nodal chemistry variable. 

logical, intent(in)                         :: iswcompbal !> iswcompbal=true, the water component is evaluated. 

logical, intent(out)                        :: iserror    !> iserror=true, there was an error. 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)              :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
iserror=.false.
msg='' 
!%---------------------------------------------------------
!% Check if the chemical system is associated. 
!%---------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>  msg='Error, not associated chemical system'
>  goto 10
end if
!%--------------------------------------------------------
call set_iswcompbal_ (this%pchemsys,iswcompbal,iserror)
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
return

10 continue 
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service:set_iswcompbal_'
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
subroutine set_from_c_nch &
>  (this, &
>   c, &
>   nsp, &
>   isspecia, &
>   iserror, &
>   dtime)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Set from c vector. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)  :: this     !> Type nodal chemistry variable. 

integer, intent(in)                    :: nsp      !> Number od species. 

real*8,intent(in), dimension(nsp)      :: c        !> Concentration vector. 

logical, intent(in)                    :: isspecia !> isspecia=true, the chemical speciation from total concentration in all phases is performed. 

logical, intent(out)                   :: iserror  !> iserror=true, there was an error. 

real*8, intent(in), optional           :: dtime    !> Time increment [s]
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer,pointer                        :: &
> idpri(:) => null ()
real*8, pointer                        :: &
> u(:) => null ()
integer                                :: &
> npri, &
> ndim, &
> nsploc, &
> ntxoh, &
> nsurf, &
> i, &
> isp
real*8                                 :: &
> dtimeloc, &
> factoromgw
logical                                :: &
> havedtime, &
> isconvergence
character(len=100)                     :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
msg=' '
iserror=.false.
!%---------------------------------------------------------
!% Check if the chemical system is associated
!%---------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>  msg='Error, not associated chemical system'
>  goto 10
end if
!%---------------------------------------------------------
!% Check optional arguments 
!%---------------------------------------------------------
havedtime=present(dtime)
!%---------------------------------------------------------
nsploc=size(this%c)
if (nsploc/=nsp) then
> msg='Error in number of species'
> goto 10
end if
!%---------------------------------------------------------
!> Copy c
!%---------------------------------------------------------
this%c=c
!%----------------------------------------------------------
!% Initialice local variables 
!%----------------------------------------------------------
factoromgw=1.0d0
!%----------------------------------------------------------
!> Make the chemical speciation (Optional)
!%----------------------------------------------------------
if (isspecia) then
>  if (havedtime) then
>   dtimeloc=dtime
>  else
>   dtimeloc=0.0d0
>  end if
!%----------------------------------------------------------
>  call get_chem_info_ (this%pchemsys,iserror,idbase=idpri,numbase=npri)
!%----------------------------------------------------------
!% Compute t=U*c
!%----------------------------------------------------------
>  call make_lin_trf_ (this%pchemsys,u,this%c,0,iserror)
!%----------------------------------------------------------
>  ntxoh=size(this%txoh,1)
>  nsurf=size(this%txoh,2)
!%------------------------------------------------------------
> this%alpha=this%alpha/this%omgwfree
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
> call update_ (this%pchemsys,this%temp,iserror)
> if (iserror) goto 20
!%-----------------------------------------------------------
!%------------------------------------------------------------
!% If the derivatives are stored
!%------------------------------------------------------------
>  if (this%isderstored) then
>    call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	.true., &
	this%ispgasconst, &
>    this%ionstr, &
>    this%alpha, &
>    u, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    nsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    dtimeloc, &
>    this%hashcompz, &
>    isconvergence, &
>    factoromgw, &
	this%volgas, &
	this%pgas, &
	msg, &
>    iserror, &
>    cguess=this%c(idpri), &
>    dc=this%dc, &
>    sktrk=this%sktrk, &
>    dsktrk=this%dsktrk, &
>    nchemiter=this%nchemiter)
>  else
>    call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	.true., &
	this%ispgasconst, &
>    this%ionstr, &
>    this%alpha, &
>    u, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    nsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    dtimeloc, &
>    this%hashcompz, &
>    isconvergence, &
>    factoromgw, &
	this%volgas, &
	this%pgas, &
	msg, &
>    iserror, &
>    cguess=this%c(idpri), &
>    sktrk=this%sktrk, &
>    nchemiter=this%nchemiter)
>   end if
>   this%alpha0=this%alpha0*this%omgwfree
>   this%alpha=this%alpha*this%omgwfree
!%-----------------------------------------------------------
!% Update the mass of free water 
!%-----------------------------------------------------------   
>   this%omgwfree = this%omgwfree * factoromgw
!%-----------------------------------------------------------
> if (.not.isconvergence) then
>  msg='Convergence problems'
>  goto 20
> end if   
!%-----------------------------------------------------------
!% Update reactive surface of kinetic minerals 
!%-----------------------------------------------------------
call update_min_area_ (this,iserror)
if (iserror) goto 20
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (idpri,1,.false.)
call check_pointer_ (u,1,.false.)
!%------------------------------------------------------------
> if (iserror) goto 10
end if
!%------------------------------------------------------------
return
> 
> 
10 continue 
print *,'********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service:set_from_c_'
print *, msg
print *,'********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_from_cpri_nch &
>  (this, &
>   nchk, &
>   cpri, &
>   dtime, &
>   isreset, &
>   isanomalous, &
>   isupmxitergam, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Set the concentrations (and dependent variables) according concentration of primary species. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)      :: this          !> Type nodal chemistry variable (in k+1). 

type(t_nodalchemistry), intent(in)         :: nchk          !> Type nodal chemistry variable (in k). 

real*8, intent(in), dimension(:)           :: cpri          !> Concentration of primary species. 

logical, intent(in)                        :: isreset       !> isreset=true, the previous concentration vector are zeroing. 

real*8, intent(in)                         :: dtime         !> Time increment [s]

logical, intent(out)                       :: iserror       !> iserror=true, there was an error.

logical, intent(out)                       :: isanomalous   !> isanomalous=true, there was anomalous concentrations during speciation calculation.

logical, intent(out)                       :: isupmxitergam !> isupmxitergam=true, the maximum number of gamma iterations was exceeded. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                    :: &
> nsp, &
> npri, &
> ntxoh, &
> nsurf
real*8                                     :: &
> faccap   !> Capillary correction for water activity 
real*8, pointer                            :: &
> alpha(:) => null ()
character(len=100)                         :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
!% Check if the chemical system is associated 
!%----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%----------------------------------------------------------
!% Set the  number of iterations to one 
!%----------------------------------------------------------
this%nchemiter=1 
!%-----------------------------------------------------------
nsp=size(this%c)
ntxoh=size(this%txoh,1)
nsurf=size(this%txoh,2)
npri=size(cpri)
!%-----------------------------------------------------------
call check_pointer_ (alpha,nsp,.true.)
alpha=this%alpha/this%omgwfree
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
call update_ (this%pchemsys,this%temp,iserror)
if (iserror) goto 20
!%-----------------------------------------------------------
!% Compute the capillary correction for water activity 
!%-----------------------------------------------------------
call get_chem_info_ (this,iserror,faccap=faccap)
if (iserror) goto 20
!%-----------------------------------------------------------
!%  Call specia_from_cpri_ in the chemicla system object
!%-----------------------------------------------------------
if (this%isderstored) then 
> call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	nchk%c*nchk%omgwfree/this%omgwfree, &  
>    alpha, &
>    nsp, &
>    cpri, &
>    npri, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    ntxoh, &
>    nsurf, &
>    dtime, &
	this%volgas, &
>    this%ionstr, &
>    this%hashcompz, &
	faccap, &    !> Capillary correction for water activity 
>    isanomalous, &
>    isupmxitergam, &
>    isreset, &
>    iserror, &
>    sktrk=this%sktrk, &
>    dc=this%dc, &
>    dsktrk=this%dsktrk)
else
>  call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	nchk%c*nchk%omgwfree/this%omgwfree, &
>    alpha, &
>    nsp, &
>    cpri, &
>    npri, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    ntxoh, &
>    nsurf, &
>    dtime, &
	this%volgas, &
>    this%ionstr, &
>    this%hashcompz, &
	faccap, &    !> Capillary correction for water activity 
>    isanomalous, &
>    isupmxitergam, &
>    isreset, &
>    iserror, &
>    sktrk=this%sktrk)
> end if 
20 continue 
!%----------------------------------------------------------
!% Deallocate local pointer 
!%----------------------------------------------------------
call check_pointer_ (alpha,1,.false.)
!%----------------------------------------------------------
if (iserror) then
> msg='Error when calling specia_'
> goto 10
end if
!%------------------------------------------------------------
return
> 
10 continue 
print *,'****************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service:set_from_cpri_'
print *, msg
print *,'****************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_from_delta_cpri_and_check_convergence_nch &
>  (this, &
>   nchk, &
>   isconvergence, &
>   isanomalous, &
>   isupmxitergam, &
>   mxerrorcpri, &
>   mxerrorres, &
>   namxerrorcpri, &
>   namxerrorres, &
>   deltacpri, &
>   residual, &
>   npri, &
>   tolcpri, &
>   tolres, &
>   factor, &
>   dtime, &
>   isreset, &
>   zero, &
>   mxconc, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Check the convergence computing the relative error and the residual error. 
!> After update the nodal chemistry object from cpri.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)       :: this            !> Type nodal chemistry variable (in k+1). 

type(t_nodalchemistry), intent(in)          :: nchk            !> Type nodal chemistry variable (in k). 

integer, intent(in)                         :: npri            !> Number of primary species

real*8, intent(in)                          :: tolcpri         !> Tolerance in unknowns 

real*8, intent(in)                          :: tolres          !> Tolerance in residual 

real*8, intent(in)                          :: dtime           !> Time increment [s]

real*8, intent(in)                          :: factor          !> Factor to scale the solution. 

real*8, intent(in)                          :: zero            !> Zero for concentrations 

real*8, intent(in), dimension(npri)         :: deltacpri       !> Newton Raphson solution 

real*8, intent(in), dimension(npri)         :: residual        !> Residual of Newton Raphson 

logical, intent(in)                         :: isreset         !> isreset=true, zeroing the concentration vector

real*8, intent(in)                          :: mxconc          !> Maximum concentration 

character(len=*), intent(out)               :: namxerrorcpri   !> Name of the components where the error is maximum in unknowns.

character(len=*), intent(out)               :: namxerrorres    !> Name of the components where the error is maximum in residual.

logical, intent(out)                        :: iserror         !> iserror=true, there was an error 

logical, intent(out)                        :: isconvergence   !> isconvergence=true, there was convergence 

logical, intent(out)                        :: isanomalous     !> isanomalous=true, there was anomalous concentrations. 

logical, intent(out)                        :: isupmxitergam   !> isupmxitergam=true, maximum number of gamma iteration 

real*8, intent(out)                         :: mxerrorcpri     !> Maximum relative error in unknowns. 

real*8, intent(out)                         :: mxerrorres      !> Maximum absolute error in residual 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
real*8, pointer              :: &
> error(:) => null (), &
> cpri(:) => null (), &
> delta(:) => null (), &
> res(:) => null ()
integer                      :: &
> npriloc, &
> wcindex, &
> i
character(len=100), pointer  :: &
> namepri(:) => null ()
integer, pointer             :: &
> icomp(:) => null ()
character(len=100)           :: &
> msg 
real*8, parameter            :: &
> r0=0.0d0, &
> r1=1.0d0
!-------------------------------------------------------------------------
!
!>   $code
!

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
!% Initializing variables 
!%----------------------------------------------------------
isconvergence=.false.
isanomalous=.false.
namxerrorcpri=''
namxerrorres=''
!%----------------------------------------------------------
!% Check if the chemical system is associated 
!%----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%----------------------------------------------------------
!% Get chemical information 
!%----------------------------------------------------------
call get_chem_info_ (this,iserror,npri=npriloc,cpri=cpri)
if (iserror) goto 20
!%----------------------------------------------------------
call get_chem_info_(this%pchemsys,iserror,namebase=namepri, &
>                    wcindex=wcindex,hashcompz=this%hashcompz)
if (iserror) goto 20
!%----------------------------------------------------------
!% Check the number of components 
!%----------------------------------------------------------
if (npriloc/=npri) then
> msg='Error in number of components'
> goto 20
end if 
!%----------------------------------------------------------
!% Allocate local pointers
!%----------------------------------------------------------
call check_pointer_ (delta,npri,.true.)
call check_pointer_ (error,npri,.true.)
call check_pointer_ (icomp,1,.true.)
call check_pointer_ (res,npri,.true.)
!%----------------------------------------------------------
!% Copy increment in unknowns and residual 
!%----------------------------------------------------------
delta=deltacpri
res=residual
!%----------------------------------------------------------
!% Compute relative error and scale according to factor 
!%----------------------------------------------------------
where (delta>r0)
>  error=factor*delta/((factor-r1)*cpri)
elsewhere
>  error=factor*delta/((r1/factor-r1)*cpri)
end where
!%----------------------------------------------------------
!% Not consider the water mol balance
!% Here must be carefully
!% Sometimes, the water primary species is different
!% than aqueous components 
!%----------------------------------------------------------
if (wcindex>0) error(wcindex)=r0
!%----------------------------------------------------------
!% Compute maximum error in components
!%----------------------------------------------------------
mxerrorcpri=maxval(dabs(error))
!%----------------------------------------------------------
if (mxerrorcpri>factor) then
> delta = delta * factor/mxerrorcpri
end if
!%----------------------------------------------------------
!% Update cpri with delta increment
!%----------------------------------------------------------
cpri = cpri + delta
!%----------------------------------------------------------
!% Check NAN numbers and zero concentrations 
!%----------------------------------------------------------
do i=1,npri
>  if (isnan(cpri(i))) then
>     isanomalous=.true.
>     goto 20      
>  end if 
>  if (cpri(i)<zero) then
>   cpri(i)=zero
>   delta(i)=r0
>   res(i)=r0
>  end if 
end do 
!%----------------------------------------------------------
!% Check negative concentrations and maximum concentrations
!% in cpri
!%----------------------------------------------------------
if (wcindex>0) then
>  cpri(wcindex)=r1
end if 
if (minval(cpri)<r0.or.maxval(cpri)>mxconc) then
> isanomalous=.true.
> goto 20 
end if
!%----------------------------------------------------------
!% Set the water concentration to 55.508 [mol/kgw] 
!%----------------------------------------------------------
if (wcindex>0) then
>  cpri(wcindex)=r1/kgwmol
end if 
!%----------------------------------------------------------
!% Specia from cpri
!%----------------------------------------------------------
call set_from_cpri_(this,nchk,cpri,dtime,isreset, &
>                    isanomalous,isupmxitergam,iserror)
if (iserror.or.isanomalous.or.isupmxitergam) goto 20
!%----------------------------------------------------------
!% Compute maximum error in components (only when 
!% deltacpri>zero
!%----------------------------------------------------------
where (dabs(cpri)>zero)
> error=delta/cpri
elsewhere
> error=delta
end where
!%----------------------------------------------------------
!% Not consider the water mol balance
!%----------------------------------------------------------
if (wcindex>0) error(wcindex)=r0
!%----------------------------------------------------------
mxerrorcpri=maxval(dabs(error))
!%----------------------------------------------------------
!% Storage the name of the components where the relative 
!% error is maximum
!%----------------------------------------------------------
icomp=maxloc(dabs(error))
namxerrorcpri=namepri(icomp(1))
!%----------------------------------------------------------
!% Change unit of residual [mol/s] => [mol/kgw]
!%----------------------------------------------------------
error=res*dtime/this%omgwfree
!%----------------------------------------------------------
!% Not consider the water mol balance ??????
!% Sometime, the water primary species not equal to 
!% water component 
!%----------------------------------------------------------
if (wcindex>0) then
>  error(wcindex)=r0
end if
!%----------------------------------------------------------
!> Compute maximum error in residual
!%----------------------------------------------------------
mxerrorres=maxval(dabs(error))
!%----------------------------------------------------------
!% Storage the name of the components where the error in 
!% residual is maximum
!%----------------------------------------------------------
icomp=maxloc(dabs(error))
namxerrorres=namepri(icomp(1))
!%----------------------------------------------------------
!% Determine convergence according relative error and residual
!%----------------------------------------------------------
isconvergence = (mxerrorcpri<=tolcpri) .and. (mxerrorres<=tolres)
!%----------------------------------------------------------
20 continue 
!%----------------------------------------------------------
!% Deallocate local pointers
!%----------------------------------------------------------
call check_pointer_ (error,1,.false.)
call check_pointer_ (delta,1,.false.)
call check_pointer_ (cpri,1,.false.)
call check_pointer_ (icomp,1,.false.)
call check_pointer_ (namepri,1,.false.)
call check_pointer_ (res,1,.false.)
if (iserror) goto 10
!%----------------------------------------------------------
return
> 
10 continue 
print *,'****************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service:update_and_check_'
print *, msg
print *,'****************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_r_and_check_convergence_sia_nch &
>   (this, &
>    nchk, &
>    rsia, &
>    utra, &
>    npri, &
>    mxerrsia, &
>    namemxerrsia, &
>    dtime, &
>    theta, &
>    tolsia, &
	factor, &
>    isconvsia, &
>    isconvchem, &
>    iter, &
>    isopensystem, &
>    zero, &
>    iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Compute the reaction term for SIA (rsia) and check 
!> convergence betwen transport and chemical step (if iteration > 1).
!
!>   $Arguments:
!
> 
type (t_nodalchemistry), intent(inout) :: this          !> Type nodal chemistry variable (in k+1). 

type (t_nodalchemistry), intent(in)    :: nchk          !> Type nodal chemistry variable (in k). 

integer, intent(in)                    :: npri          !> Number of components

integer, intent(in)                    :: iter          !> Iteration 

real*8, intent(in), dimension(npri)    :: utra          !> Transport step solution 

real*8, intent(out), dimension(npri)   :: rsia          !> Reaction term of SIA in k+1

real*8, intent(out)                    :: mxerrsia      !> Maximum error 

real*8, intent(in)                     :: dtime         !> Time increment 

real*8, intent(in)                     :: theta         !> Temporal weight 

real*8, intent(in)                     :: tolsia        !> Tolerance for SIA

real*8, intent(in)                     :: factor        !> Temporal weight 

real*8, intent(in)                     :: zero          !> Zero for SIA

logical, intent(out)                   :: iserror       !> If true there was a error.

logical, intent(out)                   :: isconvsia     !> If .true. there was convergence between transport and chemical step. 

logical, intent(out)                   :: isconvchem    !> If .true. there was convergence in cheistry step 

logical, intent(in)                    :: isopensystem  !> If .true. then system is considered open 

character(len=*), intent(out)          :: namemxerrsia  !> Name of components where error was maximum 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)                     :: &
> msg
integer                                :: &
> i, &
> j, &
> wcindex, &
> nminsp, &
> nsp, &
> npri1 
real*8                                 :: &
> mxvalue, &
> minvalue, &
> tol
real*8, pointer                        :: &
> ck1(:) => null (), &
> error(:) => null (), &
> uaq(:) => null (), &
> delta(:) => null ()
character(len=100), pointer            :: &
> namepri(:) => null ()
integer, pointer                       :: &
> idminsp(:) => null (), &
> icomp(:) => null ()
type (t_nodalchemistry), pointer       :: &
> nch => null () 
!-------------------------------------------------------------------------
!
!>   $code
!
iserror=.false.
msg=''
!%--------------------------------------------------------------
isconvsia=.false.
mxerrsia=0.0d0
namemxerrsia=''
!%--------------------------------------------------------------
!% Check if the chemical system is associated 
!%--------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%--------------------------------------------------------------
!% Allocate and create local nodal chemistry object 
!%--------------------------------------------------------------
allocate (nch)
call create_ (nch)
nch=nchk
!%--------------------------------------------------------------
!% Get information from chemical system 
!%--------------------------------------------------------------
call get_chem_info_(this%pchemsys,iserror,wcindex=wcindex,namebase=namepri,numbase=npri1)
if (iserror) goto 20 
if (npri1/=npri) then
> msg='Error in number of components'
> goto 20 
end if 
!%--------------------------------------------------------------
!% Allocate local pointers  
!%--------------------------------------------------------------
call check_pointer_ (delta,npri,.true.)
call check_pointer_ (error,npri,.true.)
call check_pointer_ (icomp,1,.true.)
!%--------------------------------------------------------------
!> 1) If the system is open zeroing mineral concentrations
!%--------------------------------------------------------------
if (isopensystem) then
> call get_chem_info_(this%pchemsys,iserror,idminsp=idminsp, &
>                     nminsp=nminsp,numsp=nsp)
> if (nminsp>0) then
>  call check_pointer_ (ck1,nsp,.true.)
>  ck1(idminsp)=this%c(idminsp)
>  this%c(idminsp)=0.0d0
>  nch%c(idminsp)=0.0d0
> end if
end if
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
!% Get uaq 
!%----------------------------------------------------------
call get_chem_info_ (this,iserror,uaq=uaq)
if (iserror) goto 20 
!%----------------------------------------------------------
!% Compute delta 
!%----------------------------------------------------------
delta=utra-uaq
!%----------------------------------------------------------
!% Compute error 
!%----------------------------------------------------------
where (delta>0.0d0)
>  error=factor*delta/((factor-1.0d0)*uaq)
elsewhere
>  error=factor*delta/((1.0d0/factor-1.0d0)*uaq)
end where
!%----------------------------------------------------------
!% Not consider the water 
!%----------------------------------------------------------
if (wcindex>0) error(wcindex)=0.0d0
!%----------------------------------------------------------
!% Compute maximum error in components
!%----------------------------------------------------------
mxerrsia=maxval(dabs(error))
!%----------------------------------------------------------
if (mxerrsia>factor) then
> delta = delta * factor/mxerrsia
end if
!%-------------------------------------------------------------
!% Update the solution 
!%-------------------------------------------------------------
rsia = uaq + delta
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
!%----------------------------------------------------------
!%-------------------------------------------------------------
!> 2) Update r sia
!%-------------------------------------------------------------
call update_r_sia_(this,nch,rsia,npri,dtime,theta,isconvchem,iserror)
if (iserror) then
>  msg='Error when calling update_r_sia_'
>  goto 10
end if
if (.not.isconvchem) goto 20 
!%----------------------------------------------------------
!% Get uaq 
!%----------------------------------------------------------
call get_chem_info_ (this,iserror,uaq=uaq)
if (iserror) goto 20 
!%--------------------------------------------------------------
!% Compute error 
!%--------------------------------------------------------------
where (uaq>0.0d0)
> error=dabs(delta/uaq)
end where 
!%--------------------------------------------------------------
!% Not consider water component 
!%--------------------------------------------------------------
if (wcindex>0) error(wcindex)=0.0d0
!%--------------------------------------------------------------
mxerrsia=maxval(error)
icomp=maxloc(error)
namemxerrsia=namepri(icomp(1))
!%--------------------------------------------------------------
!% Check convergence 
!%-------------------------------------------------------------- 
isconvsia=(mxerrsia<=tolsia)
!%--------------------------------------------------------------
20 continue 
!%--------------------------------------------------------------
!% Deallocate local pointers and types
!%--------------------------------------------------------------
call destroy_ (nch)
deallocate (nch)
call check_pointer_ (uaq,1,.false.)
call check_pointer_ (delta,1,.false.)
call check_pointer_ (error,1,.false.)
call check_pointer_ (namepri,1,.false.)
call check_pointer_ (icomp,1,.false.)
if (isopensystem.and.nminsp>0) then
> this%c=this%c+ck1
> call check_pointer_ (ck1,1,.false.)
> call check_pointer_ (idminsp,1,.false.)
end if
if (iserror) goto 10 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
return
> 
10 continue 
print *,'**************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: update_and_check_'
print *, msg
print *,'**************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_chem_info_nch &
>  (this, &
>   iserror, &
>   c, &
>   g, &
>   cmob, &
>   caq, &
>   cads, &
>   cmin, &
>   sktrk, &
>   setre, &
>   dc, &
>   dcmob, &
>   dcads, &
>   dsktrk, &
>   nummobph, &
>   numsp, &
>   omgwcryst, &
>   omgwfree, &
>   act, &
>   actwater, &
>   ionstr, &
>   cpri, &
>   npri, &
>   hashcompz, &
>   name, &
>   u, &
>   uaq, &
>   umob, &
>   uads, &
>   usktrk, &
>   molaq, &
>   molgas, &
>   molads, &
>   molprec, &
>   molkin, &
>   xsalt, &
>   xwater, &
>   dxwaterdc, &
>   chemsys, &
>   alk, &
>   simin, &
>   volnch, &
>   volgas, &
>   cmineq,&
>   namebase,&
>   pgas, &
>   faccap, &
>   temp, &
>   ph)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Return general chemical information about nodal chemistry
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)       :: this       !> Type nodal chemistry variable. 

type(t_chemicalsystem), optional, pointer:: chemsys    !> Type chemical system pointer. 

real*8, intent(out), optional            :: omgwfree   !> mass of free water

real*8, intent(out), optional            :: omgwcryst  !> mass fixed inside hydrated minerals

real*8, intent(out), optional            :: actwater   !> water activity

real*8, intent(out), optional            :: ionstr     !> Ionic strength 

real*8, intent(out), optional            :: xsalt      !> mass fraction of salt per mass of liquid

real*8, intent(out), optional            :: ph         !> ph

real*8, intent(out), optional            :: xwater     !> mass fraction of water per mass of liquid

real*8, intent(out), optional            :: faccap     !> Capillary correction for water activity 

real*8, intent(out), optional            :: temp       !> Temperature 

real*8, optional, pointer, dimension(:)  :: dxwaterdc  !> derivative of mass fraction of water per mass of liquid

real*8, intent(out), optional            :: alk        !> alkalinity

real*8, intent(out), optional            :: volnch     !> Volume associated to nodal chemistry object 

real*8, intent(out), optional            :: volgas     !> Volume of gas associated to nodal chemistry object 

real*8, intent(out), optional            :: pgas       !> Gas pressure imposed in the nodal chemistry object

real*8, optional, pointer, dimension(:)  :: c          !> concentrations

real*8, optional, pointer, dimension(:)  :: g          !> activity coefficients vector

real*8, optional, pointer, dimension(:,:):: cmob       !> mobile concentrations

real*8, optional, pointer, dimension(:)  :: caq        !> Concentrations of aqueous species 

real*8, optional, pointer, dimension(:)  :: cads       !> concentrations of adsorbed species

real*8, optional, pointer, dimension(:)  :: cmin       !> concentrations of mineral species

real*8, optional, pointer, dimension(:)  :: cmineq     !> concentrations of non kinetic mineral species

real*8, optional, pointer, dimension(:)  :: sktrk      !> change in concentrations due kinetic reactions

real*8, optional, pointer, dimension(:)  :: setre      !> change in concentrations due equilibrium reactions

real*8, optional, pointer, dimension(:,:):: dc         !> derivate of concentrations

real*8, optional, pointer, dimension(:,:):: dcmob      !> derivate of mobile concentrations

real*8, optional, pointer, dimension(:,:):: dcads      !> derivate of adsorbed species

real*8, optional, pointer, dimension(:,:):: dsktrk

real*8, optional, pointer, dimension(:)  :: act        !> activity of sepecies

real*8, optional, pointer, dimension(:)  :: cpri       !> primary concentration species

real*8, optional, pointer, dimension(:)  :: u          !> total concentrations

real*8, optional, pointer, dimension(:)  :: uaq        !> Total concentrations in aqueous phase

real*8, optional, pointer, dimension(:,:):: umob       !> Total concentrations in mobile phases 

real*8, optional, pointer, dimension(:)  :: uads       !> adsorbed components

real*8, optional, pointer, dimension(:)  :: usktrk     !> Vector usktrk[npri]

real*8, optional, pointer, dimension(:)  :: molaq      !> mol in the solution 

real*8, optional, pointer, dimension(:)  :: molgas     !> mol in the gas phases 

real*8, optional, pointer, dimension(:)  :: molads     !> mol adsorbed 

real*8, optional, pointer, dimension(:)  :: molprec    !> mol precipitated 

real*8, optional, pointer, dimension(:)  :: molkin     !> mol tranferred due kinetic reactions 

real*8, pointer, dimension(:), optional  :: simin      !> Saturation indices of mineral species

integer, intent(out), optional           :: numsp      !> number of species 

integer, intent(out), optional           :: nummobph   !> number of mobile phases 

integer, intent(out), optional           :: npri       !> number of components

integer, intent(out), optional           :: hashcompz  !> Hash index of the components zone 

logical, intent(out)                     :: iserror    !> If true there was a error.

character(len=100), intent(out), optional:: name       !> Name of the nodal chemistry object. 

character(len=100), dimension(:),pointer,optional:: namebase !name of the component definition

> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
> 
!-------------------------------------------------------------------------
!
!>   $code
!
> 
character(len=100)                       :: &
> msg
logical                                  :: &
> havec, &
> haveg, &
> havecmob, &
> havecaq, &
> havecads, &
> havecmin, &
> havecmineq, &
> havedcmob, &
> havedcads, &
> havedc, &
> havesktrk, &
> havesetre, &
> havedsktrk, &
> havenumsp, &
> havenummobph, &
> haveomgwfree, &
> haveomgwcryst, &
> haveact, &
> haveactwater, &
> haveionstr, &
> havenpri, &
> havecpri, &
> havehashcompz, &
> havename, &
> haveu, &
> haveumob, &
> haveuaq, &
> haveuads, &
> haveusktrk, &
> havemolaq, &
> havemolgas, &
> havemolads, &
> havemolprec, &
> havemolkin, &
> havexsalt, &
> havexwater, &
> havedxwaterdc, &
> havechemsys, &
> havesimin, &
> havealk, &
> havevolnch, &
> havevolgas, &
> havepgas, &
> isanomalous, &
> isupmxitergam, &
> havenamebase, &
> havefaccap, &
> havetemp, &
> haveph 
integer                     :: &
> ndim1, &
> ndim2, &
> ndim3, &
> ndim4, &
> ndim5, &
> ispw, &
> icw,  &
> ispsh, &
> nnpri
integer, pointer           :: &
> idloc(:) => null ()
real*8, pointer            :: &
> dcloc(:,:) => null (), &
> dsktrkloc(:,:) => null (), &
> vector1(:) => null (), &
> vector2(:) => null (), &
> array1(:,:) => null ()
character(len=100), pointer         :: &
> nameloc(:) => null ()
type(t_nodalchemistry), pointer :: &
> pnch => null ()
real*8::mass
!%-----------------------------------------------------------
>  msg=''
iserror=.false.
!%-----------------------------------------------------------
!% Check if the chemical system is associated
!%-----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>  msg='Error, not associated chemical system'
>  goto 10
end if
!%----------------------------------------------------------
ndim1=0
!%----------------------------------------------------------
!% Check optional arguments
!%----------------------------------------------------------
havec=present(c)
haveg=present(g)
havedc=present(dc)
havesktrk=present(sktrk)
havesetre=present(setre)
havecmob=present(cmob)
havecaq=present(caq)
havecads=present(cads)
havecmin=present(cmin)
havecmineq=present(cmineq)
havedcmob=present(dcmob)
havedcads=present(dcads)
havedsktrk=present(dsktrk)
havenumsp=present(numsp)
havenummobph=present(nummobph)
haveomgwfree=present(omgwfree)
haveomgwcryst=present(omgwcryst)
haveact=present(act)
haveactwater=present(actwater)
haveionstr=present(ionstr)
havecpri=present(cpri)
havenpri=present(npri)
havehashcompz=present(hashcompz)
havename=present(name)
haveu=present(u)
haveuaq=present(uaq)
haveumob=present(umob)
haveuads=present(uads)
haveusktrk=present(usktrk)
havemolaq=present(molaq)
havemolgas=present(molgas)
havemolads=present(molads)
havemolprec=present(molprec)
havemolkin=present(molkin)
havexsalt=present(xsalt)
havexwater=present(xwater)
havedxwaterdc=present(dxwaterdc)
havechemsys=present(chemsys)
havealk=present(alk)
havesimin=present(simin)
havevolnch=present(volnch)
havevolgas=present(volgas)
havenamebase=present(namebase)
havepgas=present(pgas)
havefaccap=present(faccap)
havetemp=present(temp) 
haveph=present(ph)
!%---------------------------------------------------------
!% Return dc, dsktrk or dcmob 
!%---------------------------------------------------------
if ((havedc.or.havedsktrk.or.havedcmob) &
>                     .and. &
>          .not.this%isderstored) then
> allocate(pnch)
> call create_ (pnch)
> pnch=this 
> ndim1=size(pnch%c)
> ndim2=size(pnch%txoh,1)
> ndim4=size(pnch%txoh,2)
> ndim5=size(pnch%alpha0,1)
> call get_chem_info_ (pnch%pchemsys,iserror,idbase=idloc,hashcompz=pnch%hashcompz)
> if (iserror) then
>   msg='Error when calling get_chem_info_'
>   goto 20
> end if
> ndim3=size(idloc)
> call check_pointer_ (vector1,ndim1,.true.)
> call check_pointer_ (vector2,ndim1,.true.)
> vector1=pnch%c
> vector2=pnch%alpha / pnch%omgwfree
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
> call update_ (pnch%pchemsys,pnch%temp,iserror)
> if (iserror) goto 20
!%-----------------------------------------------------------
> call specia_ &
>   (pnch%pchemsys, &
>    pnch%temp, &
>    pnch%c, &
>    pnch%g, &
	vector1, &
>    vector2, &
>    ndim1, &
>    pnch%c(idloc), &
>    ndim3, &
>    pnch%txoh/pnch%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    ndim2, &
>    ndim4, &
>    0.0d0, &
>    pnch%volgas, &
>    pnch%ionstr, &
>    pnch%hashcompz, &
	1.0d0, &
>    isanomalous, &
>    isupmxitergam, &
>    .false., &
>    iserror, &
>    dc=dcloc, &
>    dsktrk=dsktrkloc)
> call destroy_ (pnch)
> deallocate (pnch)  
> if (iserror.or.isanomalous.or.isupmxitergam) then
>   iserror=.true.
>   msg='Error when calling specia_'
>   goto 20
> end if
> 
> 
end if
!%---------------------------------------------------------
!% Return the concentrations vector
!%---------------------------------------------------------
if (havec) then
> ndim1=size(this%c)
> call check_pointer_ (c,ndim1,.true.)
> c=this%c
end if
!%---------------------------------------------------------
!% Return the ph
!%---------------------------------------------------------
if (haveph) then
> ph=7.0d0
> call get_sp_index_ (this%pchemsys,'h+',ispsh)
> if (ispsh>0) then
>    ph=-dlog10(this%c(ispsh)*this%g(ispsh)) 
> end if
end if
!%---------------------------------------------------------
!% Return pointer to chemical system object
!%---------------------------------------------------------
if (havechemsys) then
> chemsys => this%pchemsys
end if
!%---------------------------------------------------------
!% Return the temperature 
!%---------------------------------------------------------
if (havetemp) then
> temp = this%temp 
end if
!%---------------------------------------------------------
!% Return the activity coefficients vector
!%---------------------------------------------------------
if (haveg) then
> ndim1=size(this%g)
> call check_pointer_ (g,ndim1,.true.)
> g=this%g
end if
!%---------------------------------------------------------
!% Return the mass fraction of salt in liquid
!%---------------------------------------------------------
if (havexsalt) then
>  ndim1=size(this%c)
>  call compute_mass_salt_(this%pchemsys,xsalt,this%c,ndim1,iserror)
>  if (iserror) goto 20
>  xsalt=xsalt/(1.0d0+xsalt)
end if
!%---------------------------------------------------------
!% Return the water mass fraction in liquid
!%---------------------------------------------------------
if (havexwater) then
>  ndim1=size(this%c)
>  call compute_mass_salt_(this%pchemsys,xwater,this%c,ndim1,iserror)
>  if (iserror) goto 20
>  xwater=1.0d0/(1.0d0+xwater)
end if
!%---------------------------------------------------------
!% Return the ionic strength 
!%---------------------------------------------------------
if (haveionstr) then
>  ionstr=this%ionstr 
end if
!%---------------------------------------------------------
!% Return the derivative of water mass fraction in liquid
!%---------------------------------------------------------
if (havedxwaterdc) then
>  !1st we get the solute mass
>  ndim1=size(this%c)
>  call compute_mass_salt_(this%pchemsys,mass,this%c,ndim1,iserror)
>  
>  !2dly we get the derivative of solute mass
>   call compute_dmassSalt_dc_(this%pchemsys,this%dc,dxwaterdc,iserror)

>  !Mass fracction of water is w=1/(1+mass)
>  !So the derivative is dwdc= -(1+mass)**(-2)*(dmass/dc)
>  dxwaterdc=-1d0*(1d0+mass)**(-2d0)*dxwaterdc
end if
!%---------------------------------------------------------
!% Return the volume of the nodal chemistry object 
!%---------------------------------------------------------
if (havevolnch) then
>  volnch = this%volnch 
end if
!%---------------------------------------------------------
!% Return the volume of gas associated in the nodal 
!% chemistry object 
!%---------------------------------------------------------
if (havevolgas) then
>  volgas = this%volgas
end if
!%---------------------------------------------------------
!% Return the saturation indices of minerals [numsp]!!!!!
!%---------------------------------------------------------
if (havesimin) then
>  ndim1=size(this%c)
>  call compute_min_sat_index_ &
>   (this%pchemsys, &
>    simin, &
>    idloc, &
>    nameloc, &
>    ndim2, &
>    this%c, &
>    this%g, &
>    ndim1, &
>    iserror)
>  if (iserror) goto 20
end if
!%---------------------------------------------------------
!% Return dc[numsp,npri]
!%---------------------------------------------------------
if (havedc.and.this%isderstored) then
>  ndim1=size(this%dc,1)
>  ndim2=size(this%dc,2)
>  call check_pointer_ (dc,ndim1,ndim2,.true.)
>  dc=this%dc
else if(havedc.and..not.this%isderstored) then
>  ndim1=size(dcloc,1)
>  ndim2=size(dcloc,2)
>  call check_pointer_ (dc,ndim1,ndim2,.true.)
>  dc=dcloc
end if
!%---------------------------------------------------------
!% Return sktrk[numsp] array 
!%---------------------------------------------------------
if (havesktrk) then
> ndim1=size(this%sktrk)
> call check_pointer_ (sktrk,ndim1,.true.)
> sktrk=this%sktrk
end if
!%---------------------------------------------------------
!% Return setre[numsp] array
!%---------------------------------------------------------
if (havesetre) then
> ndim1=size(this%setre)
> call check_pointer_ (setre,ndim1,.true.)
> setre=this%setre
end if
!%---------------------------------------------------------
!% Return dsktrk[numsp,npri] array
!%---------------------------------------------------------
if (havedsktrk.and.this%isderstored) then
>  ndim1=size(this%dsktrk,1)
>  ndim2=size(this%dsktrk,2)
>  call check_pointer_ (dsktrk,ndim1,ndim2,.true.)
>  dsktrk=this%dsktrk
else if (havedsktrk.and..not.this%isderstored) then
>  ndim1=size(dsktrkloc,1)
>  ndim2=size(dsktrkloc,2)
>  call check_pointer_ (dsktrk,ndim1,ndim2,.true.)
>  dsktrk=dsktrkloc
end if
!%---------------------------------------------------------
!% Return dcmob[numsp,npri]
!%---------------------------------------------------------
if (havedcmob) then
> ndim4=size(this%c)
> if (this%isderstored) then
>  ndim5=size(this%dc,2)
>  call compute_dcmob_ &
>   (this%pchemsys, &
>    dcmob, &
>    ndim1, &
>    ndim2, &
>    ndim3, &
>    this%dc, &
>    ndim4, &
>    ndim5, &
>    iserror)
> else
>  ndim5=size(dcloc,2)
>  call compute_dcmob_ &
>   (this%pchemsys, &
>    dcmob, &
>    ndim1, &
>    ndim2, &
>    ndim3, &
>    dcloc, &
>    ndim4, &
>    ndim5, &
>    iserror)
> end if
> 
> if (iserror) then
>  msg='Error when calling compute_dcmob_'
>  goto 20
> end if
> 
end if
!%---------------------------------------------------------
!% Return dcads[numsp,npri]
!%---------------------------------------------------------
if (havedcads) then
> ndim4=size(this%c)
> if (this%isderstored) then
>  ndim5=size(this%dc,2)
>  call compute_dcads_ &
>   (this%pchemsys, &
>    dcads, &
>    ndim1, &
>    ndim2, &
>    this%dc, &
>    ndim4, &
>    ndim5, &
>    iserror)
>   if (iserror) then
>     msg='Error when calling compute_dcads_'
>     goto 20
>   end if
> else
>  ndim5=size(dcloc,2)
>  call compute_dcads_ &
>   (this%pchemsys, &
>    dcads, &
>    ndim1, &
>    ndim2, &
>    dcloc, &
>    ndim4, &
>    ndim5, &
>    iserror)
>   if (iserror) then
>     msg='Error when calling compute_dcads_'
>     goto 20
>   end if
> end if
> 
> 
> 
end if
!%---------------------------------------------------------
!% Return usktrk[npri]
!%---------------------------------------------------------
if (haveusktrk) then
> ndim2=size(this%c)
> call compute_usktrk_ &
>   (this%pchemsys, &
>    usktrk, &
>    ndim1, &
>    this%sktrk, &
>    ndim2, &
>    this%hashcompz, &
>    iserror)
> if (iserror) goto 20
end if
!%---------------------------------------------------------
!% Return cmob[numsp]
!%---------------------------------------------------------
if (havecmob) then
> ndim2=size(this%c)
> call  select_cmob_ &
>   (this%pchemsys, &
>    cmob, &
>    ndim1, &
>    ndim2, &
>    this%c, &
>    iserror)
> 
> if (iserror) then
>  msg='Error when calling select_cmob_'
>  goto 20
> end if
> 
end if
!%---------------------------------------------------------
!% Return caq[numsp]
!%---------------------------------------------------------
if (havecaq) then
> ndim1=size(this%c)
> call  select_caq_ &
>   (this%pchemsys, &
>    caq, &
>    ndim1, &
>    this%c, &
>    iserror)
> 
> if (iserror) then
>  msg='Error when calling select_caq_'
>  goto 20
> end if
> 
end if
!%---------------------------------------------------------
!% Return uaq[ncomp]
!%---------------------------------------------------------
if (haveuaq) then
> ndim1=size(this%c)
> call  select_caq_ (this%pchemsys,vector1,ndim1,this%c,iserror)
> if (iserror) then
>  msg='Error when calling select_caq_'
>  goto 20
> end if
> if (havehashcompz) then
>   ndim2=hashcompz
> else
>   ndim2=0
> end if 
> call make_lin_trf_(this%pchemsys,uaq,vector1,ndim2,iserror)
> if (iserror) goto 20
end if
!%---------------------------------------------------------
!% Return cmin[numsp]
!%---------------------------------------------------------
if (havecmin) then
> ndim2=size(this%c)
> call  select_cmin_ &
>   (this%pchemsys, &
>    cmin, &
>    ndim2, &
>    this%c, &
>    iserror)
> 
> if (iserror) then
>  msg='Error when calling select_cmin_'
>  goto 20
> end if
> 
end if
!%---------------------------------------------------------
!% Return cmineq[numsp]
!%---------------------------------------------------------
if (havecmineq) then
> ndim2=size(this%c)
> call  select_cmineq_ &
>   (this%pchemsys, &
>    cmineq, &
>    ndim2, &
>    this%c, &
>    iserror)
> 
> if (iserror) then
>  msg='Error when calling select_cmineq_'
>  goto 20
> end if
> 
end if
!%----------------------------------------------------------
!% Return cads
!%----------------------------------------------------------
if (havecads) then
> ndim2=size(this%c)
> call  select_cads_ &
>   (this%pchemsys, &
>    cads, &
>    ndim2, &
>    this%c, &
>    iserror)
> 
> if (iserror) then
>  msg='Error when calling select_cads_'
>  goto 20
> end if
> 
end if
!%----------------------------------------------------------
!% Return the mass of free water
!%----------------------------------------------------------
if (haveomgwfree) then
> omgwfree=this%omgwfree
end if
!%----------------------------------------------------------
!% Return the mass of water inside of hydrated minerals
!%----------------------------------------------------------
if (haveomgwcryst) then
> ndim1=size(this%c)
> call compute_omgwcryst_ &
>   (this%pchemsys, &
>    omgwcryst, &
>    this%omgwfree, &
>    this%c, &
>    ndim1, &
>    msg, &
>    iserror)
> if (iserror) goto 20
end if
!%----------------------------------------------------------
!% Return the activity vector 
!%----------------------------------------------------------
if (haveact) then
> call get_chem_info_ (this%pchemsys,iserror,wspindex=ispw)
> if (iserror) goto 20 
> ndim1=size(this%c)
> call check_pointer_ (act,ndim1,.true.)
> act=this%c*this%g
> if (ispw>0) act(ispw)=this%g(ispw)
end if
!%----------------------------------------------------------
!> Return the water activity 
!%----------------------------------------------------------
if (haveactwater) then
> actwater=0.0d0
> call get_chem_info_ (this%pchemsys,iserror,wspindex=ispw)
> if (iserror) goto 20 
> if (ispw>0) actwater=this%g(ispw)
end if
!%----------------------------------------------------------
!% Return the concentrations of primary species 
!% (only aqueous)
!%----------------------------------------------------------
if (havecpri) then
> call get_chem_info_(this%pchemsys,iserror,idbase=idloc, &
>                     hashcompz=this%hashcompz)
> if (iserror) goto 20                      
> ndim1=size(idloc)
> call check_pointer_ (cpri,ndim1,.true.)
> cpri=this%c(idloc)
end if
!%----------------------------------------------------------
!% Return the hash of the components zone
!%----------------------------------------------------------
if (havehashcompz) then
>  hashcompz=this%hashcompz
end if
!%----------------------------------------------------------
!% Return the name of nodal chemistry object 
!%----------------------------------------------------------
if (havename) then
>  name=this%name
end if
!%----------------------------------------------------------
!% Return the total concentrations in all phases 
!% In this case we consider the components definition of the 
!% nodal chemistry object. 
!%----------------------------------------------------------
if (haveu) then
> call make_lin_trf_(this%pchemsys,u,this%c,this%hashcompz,iserror)
> if (iserror) goto 20
end if
!%---------------------------------------------------------
!% Return the umob[npri,nmob]
!%---------------------------------------------------------
if (haveumob) then
> ndim1=size(this%c)
> call compute_umob_(this%pchemsys,umob,ndim2,ndim3,ndim1,this%c,this%hashcompz,iserror)
> if (iserror) goto 20
end if
!%---------------------------------------------------------
!% Return uads[npri]
!%---------------------------------------------------------
if (haveuads) then
> ndim1=size(this%c)
> call compute_uads_(this%pchemsys,uads,ndim2,this%c,ndim1,this%hashcompz,iserror)
> if (iserror) goto 20
end if
!%----------------------------------------------------------
!% Return the mol of components in the aqueous phase 
!% (in solution)
!%----------------------------------------------------------
if (havemolaq) then
> call get_chem_info_(this%pchemsys,iserror,numsp=ndim1)
> if (iserror) goto 20
> call select_caq_(this%pchemsys,vector1,ndim1,this%c,iserror)
> if (iserror) goto 20
> call make_lin_trf_(this%pchemsys,molaq,vector1,0,iserror,isueq=.false.)
> if (iserror) goto 20
> molaq = molaq*this%omgwfree 
end if
!%----------------------------------------------------------
!% Return the mol of components in the gas phases 
!%----------------------------------------------------------
if (havemolgas) then
> call get_chem_info_(this%pchemsys,iserror,numsp=ndim1)
> if (iserror) goto 20
> call select_cgas_(this%pchemsys,vector1,ndim1,this%c,iserror)
> if (iserror) goto 20
> call make_lin_trf_(this%pchemsys,molgas,vector1,0,iserror,isueq=.false.)
> if (iserror) goto 20
> molgas = molgas*this%omgwfree 
end if
!%----------------------------------------------------------
!% Return the adsorbed mol of components
!%----------------------------------------------------------
if (havemolads) then
> call get_chem_info_(this%pchemsys,iserror,numbase=ndim1,numsp=ndim2)
> if (iserror) goto 20
> call select_cads_(this%pchemsys,vector1,ndim2,this%c,iserror)
> if (iserror) goto 20
> call make_lin_trf_(this%pchemsys,molads,vector1,0,iserror,isueq=.false.)
> if (iserror) goto 20
> molads=molads*this%omgwfree
end if
!%----------------------------------------------------------
!> Return the kinetic changes of components per unit time 
!%----------------------------------------------------------
if (havemolkin) then
> call make_lin_trf_(this%pchemsys,molkin,this%sktrk,0,iserror,isueq=.true.)
> if (iserror) goto 20
> molkin=molkin*this%omgwfree
end if
!%----------------------------------------------------------
!% Return the precipitated mol of components
!%----------------------------------------------------------
if (havemolprec) then
> call get_chem_info_(this%pchemsys,iserror,numbase=ndim1,numsp=ndim2)
> if (iserror) goto 20
> call select_cmin_ (this%pchemsys,vector1,ndim2,this%c,iserror)
> if (iserror) goto 20
> call make_lin_trf_(this%pchemsys,molprec,vector1,0,iserror,isueq=.false.)
> if (iserror) goto 20
> molprec=molprec*this%omgwfree
end if
!%---------------------------------------------------------
!% Return the number of species
!%---------------------------------------------------------
if (havenumsp) then
> numsp=size(this%c)
end if
!%---------------------------------------------------------
!% Return the gas pressure imposed in the nodal chemistry 
!% object 
!%---------------------------------------------------------
if (havepgas) then
> pgas=this%pgas
end if
!%---------------------------------------------------------
!% Return the number of components  
!%---------------------------------------------------------
if (havenpri) then
>  call get_chem_info_(this%pchemsys,iserror,numbase=npri,hashcompz=this%hashcompz)
>  if (iserror) goto 20
end if
!%---------------------------------------------------------
!% Return the number of mobile phases 
!%---------------------------------------------------------
if (havenummobph) then
>  call get_chem_info_ (this%pchemsys,iserror,nummobph=nummobph)
>  if (iserror) goto 20 
end if
!%---------------------------------------------------------
!% Compute alkalinity of the nodal chemistry
!%---------------------------------------------------------
if (havealk) then
>  ndim1=size(this%c)
>  call compute_alkalinity_(this%pchemsys,alk,this%c,ndim1,iserror)
>  if (iserror) goto 20
end if
!%---------------------------------------------------------
!% Compute capillary correction for water activity 
!%---------------------------------------------------------
if (havefaccap) then
>  faccap=(this%pliq-this%pgas)*1.0d6*kgwmol/(rgas*(273.15d0+this%temp)*this%liqdens)
>  faccap=dexp(faccap)
end if
!%---------------------------------------------------------
!% Get the name of each component
!%---------------------------------------------------------
if (havenamebase) then
>       
>    call get_chem_info_(this%pchemsys,iserror,numbase=nnpri,hashcompz=this%hashcompz)
>    if (iserror) goto 20
>    allocate(namebase(nnpri))
>    if (havenamebase) then
>        call get_chem_info_(this%pchemsys,iserror,namebase=namebase,hashcompz=this%hashcompz)
>    else 
>        call get_chem_info_(this%pchemsys,iserror,namebase=namebase)
>    endif
>    
endif
!%----------------------------------------------------------
20 continue 
!%----------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------
call check_pointer_ (dcloc,1,1,.false.)
call check_pointer_ (dsktrkloc,1,1,.false.)
call check_pointer_ (vector1,1,.false.)
call check_pointer_ (vector2,1,.false.)
call check_pointer_ (array1,1,1,.false.)
call check_pointer_ (idloc,1,.false.)
call check_pointer_ (nameloc,1,.false.)
if (iserror) goto 10
!%----------------------------------------------------------
return
> 
10 continue 
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: get_chem_info_'
print *, msg
print *,'***********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_from_setre_nch &
>  (this, &
>   setre, &
>   nreact, &
>   nsp, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Set the concentration of equilibrium mineral species from the equilibrium term of the transport equation (Set re)
!> where Se is the stoichiometric matrix corresponding to equilibrium reactions, and re [mol/kgw/s] is the vector of equilibrium reaction rates. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout):: this           !> Type nodal chemistry variable (in k+1). 

integer, intent(in)                  :: nsp            !> Number of species

real*8, intent(inout), dimension(nsp):: setre          !> Set re [nsp]

logical, intent(out)                 :: iserror        !> iserror=true, there was a error. 

Integer, intent(out)                 :: nreact
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)                 :: &
> msg
integer                            :: &
> ndim 
!-------------------------------------------------------------------------
!
!>   $code
!
!%-----------------------------------------------------------
msg=''
iserror=.false.
!%-----------------------------------------------------------
!% Check if the chemical system is associated 
!%-----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,numsp=ndim)
if (iserror) goto 10
if (ndim/=nsp) then
> msg='Error in number of species'
> goto 10
end if
!%-----------------------------------------------------------

call compute_from_setre_ &
>   (this%pchemsys, &
>    this%c, &
>    setre, &
	nreact, &
>    nsp, &
>    iserror)

!%-----------------------------------------------------------
if (iserror) goto 10
!%-----------------------------------------------------------
!% Storage the vector setre 
!%-----------------------------------------------------------
return
10 continue 
print *,'**************************'
print *,'Nodal Chemistry:'
print *,'Namë:',this%name
print *,'Service:  get_from_setre_'
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
subroutine update_r_sia_nch &
>   (nchk1, &
>    nchk, &
>    utra, &
>    ncomp, &
>    dtime, &
>    theta, &
>    isconvergence, &
>    iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Update reaction term for SIA from transport solution. 
!> Compute reaction term for SIA
!%   according total aqueous components.
!%
!%
!%   ----------------
!%   unew(numaqprisp)= Total components in the aqueous phase.
!%
!%
!%   dt= Time interval.
!%   ith= Index of nodal chemistry object.
!%
!%
!%   Description:
!%   ------------
!%
!%   The service make the speciation according total analytic
!%   components for the ith nodal chemistry as following:
!%
!%   Tnew_ith = Told_ith + (Unew_ith - Uold_ith)
!%
!%   where
!%
!%   Tnew_ith(numaqprisp)= Total analytic in k+1 for ith
!%                         nodal chemistry
!%   Told_ith(numaqprisp)= Total analytic in k for ith
!%                         nodal chemistry
!%   Unew_ith(numaqprisp)= Total in the aqueous phase to k+1 for ith
!%                         nodal chemistry
!%   Uold_ith(numaqprisp)= Total analytic in k for ith
!%                         nodal chemistry
!%
!%   After speciation, is computed agin Unew*_ith and USktRk
!%   (if there are kinetic reaction in teh system), and R_ith
!%   is evaluated according
!%
!%   R_ith = - (Unew_ith - Unew*_ith)
!%             ---------------------- + USKtRk
!%                    dt
!%
!%   Warning: The new speciation computed is not stored in c!!!
!%
!%   Make by:
!%   Sergio Andrés Bea Jofré
!%   Polytechnical University of Catalonia
!%   2006
!
!>   $Arguments:
!
> 
type (t_nodalchemistry), intent(inout)  :: nchk1         !> Type nodal chemistry variable in k+1.

type (t_nodalchemistry), intent(in)     :: nchk          !> Type nodal chemistry variable in k. 

integer, intent(in)                     :: ncomp         !> Number of components 

real*8, intent(inout), dimension(ncomp) :: utra          !> Transport solution (in) rsia (out)

real*8, intent(in)                      :: dtime         !> Time increment [s].

real*8, intent(in)                      :: theta         !> Temporal weight.

logical, intent(out)                    :: iserror       !> iserror=true then there was an error.

logical, intent(out)                    :: isconvergence !> isconvergence=true, then there was convergence in the chemistry speciation.
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                               :: &
> i, &
> nmobph, &
> npri, &
> icw
real*8                                :: &
> factor  
real*8, pointer                       :: &
> rsia(:) => null (), & 
> usktrk(:) => null (), &
> u(:) => null (), &
> sktrk(:) => null (), &
> cmob(:,:) => null (), &
> c(:) => null (), &
> t(:) => null ()
character(len=100)                    :: &
> msg 
type(t_chemicalsystem), pointer       :: & 
> pchemsys => null ()
!-------------------------------------------------------------------------
!
!>   $code
!
iserror=.false.
msg=''
!%---------------------------------------------------------------
!% Check if the chemical system is associated in nodal chemistry 
!% objects in k and k+1 
!%---------------------------------------------------------------
if (.not.associated(nchk1%pchemsys).or..not.associated(nchk%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%---------------------------------------------------------------
!% Get chemical information 
!%---------------------------------------------------------------
call get_chem_info_ (nchk1%pchemsys,iserror,numbase=npri,wcindex=icw)
if (iserror) then
> msg='Error when calling get_chem_info_'
> goto 20
end if
!%--------------------------------------------------------------
!% Check number of components
!%--------------------------------------------------------------
if (npri/=ncomp) then
> msg='Error in number of primary species'
> goto 20
end if
!%--------------------------------------------------------------
!% Check the time increment 
!%--------------------------------------------------------------
if (dtime==0.0d0) then
> msg='Error, the time increment is 0.0d0'
> goto 20
end if
!%---------------------------------------------------------------
!% Allocate pointer 
!%---------------------------------------------------------------
call check_pointer_ (rsia,ncomp,.true.)
!%--------------------------------------------------------------
!% Get nodal chemistry information in k
!%--------------------------------------------------------------
pchemsys => nchk%pchemsys 
!%--------------------------------------------------------------
call get_chem_info_ (nchk,iserror,cmob=cmob,c=c,sktrk=sktrk)
if (iserror) goto 20
call make_lin_trf_ (pchemsys,usktrk,sktrk,0,iserror)
if (iserror) goto 20
call make_lin_trf_ (pchemsys,t,c,0,iserror)
if (iserror) goto 20
call make_lin_trf_ (pchemsys,u,cmob(:,1),0,iserror)
if (iserror) goto 20
!%----------------------------------------------------------------
!% Store the k term in rsia
!% rsia = (t - u)/dtime + (1.0d0 - theta) * usktrk
!%----------------------------------------------------------------
rsia = (t - u)/dtime + (1.0d0 - theta) * usktrk
!%-----------------------------------------------------------------
!> Correct the total mass of components 
!%-----------------------------------------------------------------
if (icw>0) utra(icw)=u(icw)
!%-----------------------------------------------------------------
!% Set total concentrations in nodal chemistry object in k+1
!%-----------------------------------------------------------------
call set_from_uaq_(nchk1,utra,ncomp,dtime,isconvergence,iserror)
if (iserror) goto 20
!%----------------------------------------------------------------
!% If there was covergence problems during speciation calculations, 
!% then return 
!%----------------------------------------------------------------
if (.not.isconvergence) goto 20
!%----------------------------------------------------------------
!% Get nodal chemical information
!%----------------------------------------------------------------
call get_chem_info_(nchk1,iserror,c=c,cmob=cmob,sktrk=sktrk) 
if (iserror) goto 20
!%----------------------------------------------------------------
pchemsys => nchk1%pchemsys
!%----------------------------------------------------------------
!% Compute in k+1
!% umob=U*cmob
!%----------------------------------------------------------------
call make_lin_trf_ (pchemsys,u,cmob(:,1),0,iserror)
if (iserror) goto 20
!%----------------------------------------------------------------
!% Compute k+1
!% t=U*c
!%----------------------------------------------------------------
call make_lin_trf_(pchemsys,t,c,0,iserror)
if (iserror) goto 20
!%----------------------------------------------------------------
!% Compute in k+1 
!% usktrk=U*sktrk
!%----------------------------------------------------------------
call make_lin_trf_(pchemsys,usktrk,sktrk,0,iserror)
if (iserror) goto 20
!%----------------------------------------------------------------
!% rsia = rsia - (t - u)/dtime + theta * usktrk
!%----------------------------------------------------------------
rsia = rsia - (t - u)/dtime + theta * usktrk
!%----------------------------------------------------------------
!% Copy rsia in utra
!%----------------------------------------------------------------
utra=rsia
!%----------------------------------------------------------------
20 continue      
!%---------------------------------------------------------------
!% Deallocate and nullify local pointers 
!%---------------------------------------------------------------
call check_pointer_ (rsia,1,.false.)
call check_pointer_ (u,1,.false.)
call check_pointer_ (usktrk,1,.false.)
call check_pointer_ (t,1,.false.)
call check_pointer_ (sktrk,1,.false.)
call check_pointer_ (cmob,1,1,.false.)
call check_pointer_ (c,1,.false.)
pchemsys => null ()
!%-----------------------------------------------------------
if (iserror) goto 10
!%-----------------------------------------------------------
return
> 
10 continue 
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',nchk1%name
print *,'Service: compute_r_sia_'
print *, msg
print *,'***********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_from_uaq_nch &
>   (this, &
>    unew, &
>    ncomp, &
>    dtime, &
>    isconvergence, &
>    iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Set the concentrations (and dependent variables) from total 
!> concentrattions (unew).
!
!>   The total concentrations (tnew) in all phases are computed as 
!
!>   tnew = told + (tnew - told)
!
!>   where
!
!>   tnew(ncomp)= Total concentrations in all phases in k+1.
!
!>   told(ncomp)= Total concentrations in all phases in k.
!
!>   unew(ncomp)= Total concentrations in k+1.

!>   uold(ncomp)= TTotal concentrations in k.
!
!>   $Arguments:
!
> 
type (t_nodalchemistry), intent(inout) :: this          !> Type nodal chemistry variable. 

integer, intent(in)                    :: ncomp         !> Number of components.  

real*8, intent(in), dimension(ncomp)   :: unew          !> Total concentration vector [ncomp]. 

real*8, intent(in)                     :: dtime         !> Time increment [s]

logical, intent(out)                   :: iserror       !> iserror=true, there was an error.

logical, intent(out)                   :: isconvergence !> isconvergence=true, there was convergence in speciation calculations. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
> nmobph, &
> i, &
> ncomploc
real*8, pointer                        :: &
> told(:) => null(), &
> tnew(:) => null(), &
> uold(:) => null(), &
> c(:) => null(), &
> cmob(:,:) => null()
character(len=100)                     :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
iserror=.false.
msg=' '
!%-------------------------------------------------------------
!% Check if the chemical system is associated 
!%-------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>   msg='Error, not associated chemical system'
>   goto 10
end if
!%-------------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,numbase=ncomploc)
if (iserror) then
> msg='Error when calling get_chem_info_'
> goto 10
end if
!%-------------------------------------------------------------
!% Check the number of components
!%-------------------------------------------------------------
if (ncomploc/=ncomp) then
>  msg='Error in number of components'
>  goto 10
end if
!%-------------------------------------------------------------
!% Allocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (tnew,ncomploc,.true.)
!%-------------------------------------------------------------
call get_chem_info_ (this,iserror,c=c,cmob=cmob)
if (iserror) goto 20
!%-------------------------------------------------------------
!% Compute told=U*cold
!%-------------------------------------------------------------
call make_lin_trf_ (this%pchemsys,told,c,0,iserror)
if (iserror) goto 20
!%-------------------------------------------------------------
!% Compute uold=U*cmobold
!%-------------------------------------------------------------
call make_lin_trf_ (this%pchemsys,uold,cmob(:,1),0,iserror)
if (iserror) goto 20
!%-------------------------------------------------------------
!% Compute tnew
!% tnew = told + unew - uold
!%-------------------------------------------------------------
tnew = told + unew - uold
!%-------------------------------------------------------------
!% Specia according total concentrations in all phases
!%-------------------------------------------------------------
call set_from_u_ (this,tnew,dtime,isconvergence,iserror)
if (iserror) then
>   msg='Error when calling set_from_u_'
>   goto 20
end if
!%-------------------------------------------------------------
20  continue 
!%-------------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (uold,1,.false.)
call check_pointer_ (told,1,.false.)
call check_pointer_ (tnew,1,.false.)
call check_pointer_ (c,1,.false.)
call check_pointer_ (cmob,1,1,.false.)
!%-------------------------------------------------------------
if (iserror) goto 10
!%-------------------------------------------------------------
return
> 
10 continue 
print *,'*********************'
print *,'Nodal Chemistry:'
print *,'Name:', this%name
print *,'Service: set_from_uaq_'
print *, msg
print *,'*********************'
iserror=.true.
return
> 
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_iumob_nch &
>  (this, &
>   iumob, &
>   naqcol, &
>   ngascol, &
>   nnonaqcol, &
>   ithcomp, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Return concentration of the ith mobile components
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in) :: this       !> Type nodal chemistry variable. 

real*8, pointer, dimension(:)      :: iumob 

integer, intent(out)               :: naqcol

integer, intent(out)               :: ngascol

integer, intent(out)               :: nnonaqcol

integer, intent(in)                :: ithcomp

logical, intent(out)               :: iserror    !> If true there was a error.
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                :: &
> numsp
character(len=100)     :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
iserror=.false.
msg='' 
!%------------------------------------------------------------
!% Check if the chemical system is associated
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
numsp=size(this%c)
call compute_iumob_ &
>   (this%pchemsys, &
>    iumob, &
>    naqcol, &
>    ngascol, &
>    nnonaqcol, &
>    this%c, &
>    numsp, &
>    ithcomp, &
>    this%hashcompz, &
>    iserror)
!%------------------------------------------------------------
> 
return
> 
10 continue 
print *,'*********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: get_iumob_'
print *, msg
print *,'*********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_name_nch &
>  (this, &
>   name)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Return the name of nodal chemistry object. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in) :: this !> Type nodal chemistry variable. 

character(len=*), intent(out)      :: name 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
> 
!-------------------------------------------------------------------------
!
!>   $code
!
name=this%name
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_dusktrk_nch &
>  (this, &
>   dusktrk, &
>   ncomp, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Compute U*Skt*drk/dc1 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in) :: this     !> Type nodal chemistry variable. 

integer, intent(out)               :: ncomp    !> Number of components. 

real*8, pointer, dimension(:,:)    :: dusktrk  !> Derivatives. 

logical, intent(out)               :: iserror  !> iserror=true, then there was an error. 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)           :: &
> msg
integer                :: &
> ndim 
!-------------------------------------------------------------------------
!
!>   $code
!
!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
!% Check if the chemical system is associated 
!%----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%----------------------------------------------------------
ndim=size(this%c)
!%----------------------------------------------------------
if (this%isderstored) then
> 
> ncomp=size(this%dsktrk,2)
> call compute_dusktrk_ &
>   (this%pchemsys, &
>    dusktrk, &
>    ncomp, &
>    this%dsktrk, &
>    ndim, &
>    this%hashcompz, &
>    iserror)
> 
else
> 
> call compute_dusktrk_ &
>   (this%pchemsys, &
>    dusktrk, &
>    ncomp, &
>    this%c, &
>    this%alpha, &
>    ndim, &
>    this%hashcompz, &
>    iserror)
> 
end if
!%---------------------------------------------------------
return
10 continue 
print *,'*********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: get_dusktrk_'
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
subroutine get_dumob_nch &
>  (this, &
>   dumob, &
>   nrow, &
>   ncol, &
>   nmobph, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Compute the derivate of the concentration with respect to 
!> primary species of ith nodal chemistry, evaluated in the components matrix 
!> corresponding to jth nodal chemistry.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in) :: this    !> Type nodal chemistry variable. 

real*8, pointer, dimension(:,:)    :: dumob   !> Derivatives of the mobile concentrations with respect to primary species

integer, intent(out)               :: nrow    !> Number components
> 
integer, intent(out)               :: ncol    !> Number of primary species 

integer, intent(out)               :: nmobph  !> Number of mobile phases 

logical, intent(out)               :: iserror !> iserror=true, then there was an error
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)     :: &
> msg
integer                :: &
> numsp, &
> naqpri 
!-------------------------------------------------------------------------
!
!>   $code
!
iserror=.false.
msg=''
!%----------------------------------------------------------
!% Check if the chemical system is associated 
!%----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%----------------------------------------------------------
numsp=size(this%c)
!%----------------------------------------------------------
!% Call the corresponding service in the chemical system
!% object 
!%----------------------------------------------------------
if (this%isderstored) then
> 
> naqpri=size(this%dc,2)
> call compute_dumob_ &
>   (this%pchemsys, &
>    dumob, &
>    nmobph, &
>    nrow, &
>    ncol, &
>    this%dc, &
>    numsp, &
>    naqpri, &
>    this%hashcompz, &
>    iserror)
> 
else
> 
> call compute_dumob_ &
>   (this%pchemsys, &
>    dumob, &
>    nrow, &
>    ncol, &
>    nmobph, &
>    this%c, &
>    numsp, &
>    this%hashcompz, &
>    iserror)
> 
end if
> 
if (iserror) then
> msg='Error when calling compute_dumob_'
> goto 10
end if
!%-----------------------------------------------------------
return
> 
10 continue 
print *,'*********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: get_dumob_'
print *, msg
print *,'*********************'
iserror=.true.
return
> 
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_duads_nch &
>  (this, &
>   duads, &
>   nrow, &
>   ncol, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Return the U*dcads/dc1
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in) :: this     !> Type nodal chemistry variable. 

real*8, pointer, dimension(:,:)    :: duads    !> Derivatives of the components according 

integer, intent(out)               :: nrow     !> Number of components

integer, intent(out)               :: ncol     !> Number of primary species

logical, intent(out)               :: iserror  !> iserror=true, then there was an error
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)     :: &
> msg
integer                :: &
> numsp 
!-------------------------------------------------------------------------
!
!>   $code
!
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Check if the chemical system is associated 
!%-----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------------
numsp=size(this%c)
call compute_duads_ &
>   (this%pchemsys, &
>    duads, &
>    nrow, &
>    ncol, &
>    numsp, &
>    this%c, &
>    this%hashcompz, &
>    iserror)
!%-----------------------------------------------------------
if (iserror) then
> msg='Error when call compute_duads_'
> goto 10
end if 
return
10 continue 
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: get_duads_'
print *, msg
print *,'***********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
function mix_nch &
>  (nch1, &
>   nch2)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Mix two nodal chemistries and store in other nodal chemistry object. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)   :: nch1     !> Type nodal chemistry variable. 

type(t_nodalchemistry), intent(in)   :: nch2     !> Type nodal chemistry variable. 

type(t_nodalchemistry)               :: mix_nch  !> Type nodal chemistry variable. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond: The type nodal chemistry object mix_nch must be created. 
!
!>   $Post-cond:
!
!>   $License: UPC-CSIC
!
!-------------------------------------------------------------------------
logical                           :: &
> iserror
character(len=100)                :: &
> msg
real*8, pointer                   :: &
> t(:) => null ()
logical                           :: &
> isconvergence 
!-------------------------------------------------------------------------
!
!>   $code
!
msg=''
!%----------------------------------------------------------
!% Check if the chemical systems are associated 
!%----------------------------------------------------------
!if (.not.associated(nch1%pchemsys,nch2%pchemsys)) then
!> msg='Error, diferent assigned chemical systems in nch1 and nch2'
!> goto 10
!end if
!%----------------------------------------------------------
!% Create nodal chemistry object 
!%----------------------------------------------------------
call create_ (mix_nch)
!%----------------------------------------------------------
!% Copy the nch1 and mix_nch
!%----------------------------------------------------------
mix_nch = nch1
!%----------------------------------------------------------
!% Add the mass of water 
!%----------------------------------------------------------
mix_nch%omgwfree = nch1%omgwfree + nch2%omgwfree
!%----------------------------------------------------------
mix_nch%volnch = nch1%volnch + nch2%volnch
!%----------------------------------------------------------
!% Add volume of gas
!%----------------------------------------------------------
mix_nch%volgas = nch1%volgas + nch2%volgas
!%----------------------------------------------------------
mix_nch%alpha0 = nch1%alpha0 + nch2%alpha0
!%----------------------------------------------------------
!% Add m2 of mineral 
!%----------------------------------------------------------
mix_nch%alpha = nch1%alpha + nch2%alpha
!%----------------------------------------------------------
!% Add CEC
!%----------------------------------------------------------
if (associated(mix_nch%txoh)) then
>  mix_nch%txoh = nch1%txoh + nch2%txoh
end if 
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
mix_nch%c = (nch1%c * nch1%omgwfree + nch2%c * nch2%omgwfree) / mix_nch%omgwfree
!%-----------------------------------------------------------
mix_nch%temp = (nch1%temp * nch1%omgwfree + nch2%temp * nch2%omgwfree) / mix_nch%omgwfree
!%-----------------------------------------------------------
mix_nch%pgas = (nch1%pgas * nch1%volgas/nch1%temp + nch2%pgas * nch2%volgas/nch2%temp)*mix_nch%temp/mix_nch%volgas
!%-----------------------------------------------------------
!There is no sence in mixing liquid presurre
mix_nch%pliq = 0
!%-----------------------------------------------------------

call make_lin_trf_ (mix_nch%pchemsys,t,mix_nch%c,0,iserror)
if (iserror) then 
>  msg='Error when call make_lin_trf_'
>  goto 20
end if 
!%-----------------------------------------------------------
!% Specia according total concentrations
!%-----------------------------------------------------------
call set_from_u_ (mix_nch,t,0.0d0,isconvergence,iserror)
if (iserror) goto 20 
!%-----------------------------------------------------------
if (.not.isconvergence) then
>  iserror=.true.
>  msg='Convergence problems when call set_from_u_'
>  goto 20
end if 
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------
call check_pointer_ (t,1,.false.)
if (iserror) goto 10 
!%----------------------------------------------------------
return
> 
10 continue 
print *,'*************************'
print *,'Nodal Chemistry:'
print *,'Name:',nch1%name
print *,'Service: (+)'
print *, msg
print *,'*************************'
iserror=.true.
stop 
> 
end function
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine reaction_path_nch &
>  (this, &
>   namesp, &
>   mol, &
>   amount, &
>   nstep, &
>   iserror, &
>   ioutput, &
>   namespout)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Similar to reaction path of phreeqc.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)                 :: this     !> Type nodal chemistry variable. 

integer, intent(in)                                   :: nstep    !> Number of steps 

character(len=*), intent(in)                          :: namesp   !> Name of the species

logical, intent(out)                                  :: iserror  !> If true there was error

real*8, intent(in)                                    :: mol      !> Mol of  added species 

real*8, intent(in)                                    :: amount

integer, intent(in), optional                         :: ioutput

character(len=*), intent(in), dimension(:), optional  :: namespout !> List of species to write 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                           :: &
> i, &
> npri, &
> nsp, &
> nsurf, &
> ntxoh, &
> istep, &
> nminsp, &
> nspout, &
> nw, &
> nwout, &
> isp, &
> isp1 
real*8                            :: &
> delta, &
> icputime1, &
> icputime2, &
> dtcputime, &
> sum, &
> factoromgw, &
> omgwfree0 
real*8, pointer                   :: &
> moladded(:) => null (), &
> cout(:,:) => null (), &
> actout(:,:) => null (), &
> t(:) => null (), &
> c(:) => null ()
integer, pointer                  :: &
> idminsp(:) => null (), &
> idspout(:) => null (), &
> idpri(:) => null ()
character(len=100), pointer       :: &
> namesploc(:) => null ()
logical, pointer                  :: &
> isokstep(:) => null ()
logical                           :: &
> haveioutput, &
> havenamespout, &
> isconvergence, &
> writehead
character(len=100)                :: &
> msg, &
> head
type(t_nodalchemistry), pointer   :: & 
> nch => null ()
!-------------------------------------------------------------------------
!
!>   $code
!
msg=' '
iserror=.false.
!%-----------------------------------------------------------------
!% Check if the chemical system object is associated to nodal 
!% chemistry object 
!%-----------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>   msg='Error, not associated chemical system'
>   goto 10
end if
!%-----------------------------------------------------------------
!% Check optional arguments 
!%-----------------------------------------------------------------
haveioutput=present(ioutput)
havenamespout=present(namespout)
!%-----------------------------------------------------------------
!% Check if the species was defined in the chemical system 
!%-----------------------------------------------------------------
call get_sp_index_(this%pchemsys,namesp,isp)
if (isp==0) then
> msg='Error, chemical species not defined in the chemical system:'
> call add_ (msg,namesp)
> goto 10 
end if  
!%-----------------------------------------------------------------
!% Allocate and create local nodal chemistry 
!%-----------------------------------------------------------------
allocate (nch)
call create_ (nch)
nch=this
omgwfree0=this%omgwfree
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror, &
>                     numsp=nsp,numsites=ntxoh,numsurf=nsurf, &
>                     idminsp=idminsp,nminsp=nminsp, &
>                     namesp=namesploc,wspindex=nw, &
>                     idbase=idpri,numbase=npri)
!%---------------------------------------------------------
!% Allocate and copy concentrations vector c 
!%---------------------------------------------------------
call check_pointer_ (c,nsp,.true.)
c=this%c 
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
nspout=0 
if (haveioutput.and.havenamespout) then
> nspout=size(namespout)
> nwout=0
> if (nspout>0) then
>  call check_pointer_ (idspout,nspout,.true.)
>  call check_pointer_ (cout,nspout,nstep,.true.)
>  call check_pointer_ (actout,nspout,nstep,.true.)
>  do i=1,nspout
>   call get_sp_index_(this%pchemsys,namespout(i),isp1)
>   idspout(i)=isp1
>   if (isp1==nw) nwout=i
>  end do
> end if
end if
!%-------------------------------------------------------
if (haveioutput) then
> 
> 
> 
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) 'Beginning of batch-reaction calculations'
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,1)  amount, 'moles in',nstep,'steps.'
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) 'Reactant         Relative moles'
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,2)  namesp,mol
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) 'Initial Speciation'
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) 'Name:',this%name
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     call write_ &
>      (this%pchemsys, &
>       this%c, &
>       this%g, &
>       this%temp, &
>       this%ionstr, &
>       nsp, &
>       ioutput, &
>       this%omgwfree, &
>       msg, &
>       iserror)
>     if (iserror) goto 20
> 
end if
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
> call update_ (nch%pchemsys,nch%temp,iserror)
> if (iserror) goto 20
!%-----------------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (moladded,nstep,.true.)
call check_pointer_ (isokstep,nstep,.true.)
!%-----------------------------------------------------------------
delta=amount*mol/real(nstep)
!%-----------------------------------------------------------------
sum=0.0d0
!%-----------------------------------------------------------------
call cpu_time (icputime1)
!%-----------------------------------------------------------------
do istep=1,nstep
>      sum=sum+delta
>      moladded(istep)=sum
>      c(isp)=c(isp)+delta 
>      nch%alpha=nch%alpha/nch%omgwfree
>      call make_lin_trf_(nch%pchemsys,t,c,0,iserror)
>      if (iserror) goto 20
!%-----------------------------------------------------------
	  if (this%isderstored) then
>        call specia_ &
>       (nch%pchemsys, &
	    nch%temp, &
>        nch%c, &
>        nch%g, &
	    .false., &
		nch%ispgasconst, &
>        nch%ionstr, &
>        nch%alpha, &
>        t, &
>        nch%txoh/nch%omgwfree, &
>        nch%capint, &
>        nch%capext, &
>        nch%spsurfarea, &
>        nsp, &
>        npri, &
>        ntxoh, &
>        nsurf, &
>        0.0d0, &
>        nch%hashcompz, &
>        isconvergence, &
>        factoromgw, &
		nch%volgas, &
	    nch%pgas, &
		msg, &
>        iserror, &
>        cguess=nch%c(idpri), &
>        sktrk=nch%sktrk, &
>        dc=nch%dc, &
>        dsktrk=nch%dsktrk, &
>        nchemiter=nch%nchemiter)
>      else 
>        call specia_ &
>       (nch%pchemsys, &
	    nch%temp, &
>        nch%c, &
>        nch%g, &
	    .false., &
		nch%ispgasconst, &
>        nch%ionstr, &
>        nch%alpha, &
>        t, &
>        nch%txoh/nch%omgwfree, &
>        nch%capint, &
>        nch%capext, &
>        nch%spsurfarea, &
>        nsp, &
>        npri, &
>        ntxoh, &
>        nsurf, &
>        0.0d0, &
>        nch%hashcompz, &
>        isconvergence, &
>        factoromgw, &
		nch%volgas, &
	    nch%pgas, &
		msg, &
>        iserror, &
>        cguess=nch%c(idpri), &
>        sktrk=nch%sktrk, &
>        nchemiter=nch%nchemiter)
>      end if 
>      nch%alpha=nch%alpha*nch%omgwfree
>      !%-------------------------------------------------------
>      !% Correct the mass of free water 
>      !%-------------------------------------------------------
>      nch%omgwfree = omgwfree0 * factoromgw
>      if (isconvergence) then
	   isokstep(istep)=.true.
	   call update_min_area_ (nch,iserror)
	   if (iserror) goto 20 
	   this=nch
	  else 
	   nch=this
	  end if 
	  if (iserror) goto 20 
!%-------------------------------------------------------
!> Write output 
!%------------------------------------------------------- 
>    if (haveioutput) then
>      if (havenamespout.and.nspout>0) then
>        cout(:,istep)=this%c(idspout)
>        actout(:,istep)=this%c(idspout)*this%g(idspout)
>        if (nwout>0) actout(nwout,istep)=this%g(nw)
>      end if
>      write(ioutput,*) 'Reaction step', istep
>      write(ioutput,3) '----------------------------------------'// &
>                       '----------------------------------------'// &
>                       '-------------'
>      if (isokstep(istep)) then  
>        call write_ (this,ioutput,iserror)
>      else
>        write(ioutput,*) 'FAILED STEP!!!!!!, step:',istep
>        write(ioutput,3) '----------------------------------------'// &
>                         '----------------------------------------'// &
>                         '-------------'
>      end if 
> 
>      if (iserror) goto 20
>    end if
!%-------------------------------------------------------
!%------------------------------------------------------- 
> 
end do
!%-------------------------------------------------------
call cpu_time (icputime2)
dtcputime=icputime2-icputime1
!%-------------------------------------------------------
!> Write the cpu time
!%-------------------------------------------------------
if (haveioutput) then
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,4) 'Total cpu time used:',dtcputime,' seconds'
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
end if
!%--------------------------------------------------------------------
!> Write evolution of concentration species vs. vs. mol added
!%--------------------------------------------------------------------
if (haveioutput.and.nspout>0) then
> head='Evolution of concentration species vs. vs. mol added'
> do i=1,nstep
>   if (isokstep(i)) then  
>     if (i==1) then
>        writehead=.true.
>     else
>        writehead=.false.
>     end if
	 call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       cout(1:nspout,i), &
>       nspout, &
>       moladded(i), &
>       writehead, &
>       writehead)
>    end if    
> end do
!%--------------------------------------------------------------------
!> Write evolution of activity species vs. vs. mol added
!%-------------------------------------------------------------------- 
> head='Evolution of activity species vs. vs. mol added'
> do i=1,nstep
>   if (isokstep(i)) then  
>     if (i==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
>     call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       actout(1:nspout,i), &
>       nspout, &
>       moladded(i), &
>       writehead, &
>       writehead)
>    end if 
> end do
end if
!%-------------------------------------------------------
!%-------------------------------------------------------
!%-------------------------------------------------------
20 continue 
!%-----------------------------------------------------------------
!% Destroy and deallocate local nodal chemistry 
!%-----------------------------------------------------------------
call destroy_ (nch)
deallocate (nch)
!%-------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------
call check_pointer_ (idspout,1,.false.)
call check_pointer_ (cout,1,1,.false.)
call check_pointer_ (actout,1,1,.false.)
call check_pointer_ (moladded,1,.false.)
call check_pointer_ (namesploc,1,.false.)
call check_pointer_ (idminsp,1,.false.)
call check_pointer_ (idpri,1,.false.)
call check_pointer_ (t,1,.false.)
call check_pointer_ (c,1,.false.)
call check_pointer_ (isokstep,1,.false.)
if (iserror) goto 10
!%-------------------------------------------------------
return
1 format (f10.2,a9,i5,a7)
2 format (a10,f10.2)
3 format (a95)
4 format (a20,e10.4,a8)
> 
10 continue 
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: reaction_path_'
print *, msg
print *,'***********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine reaction_path_kinetic_nch &
>  (this, &
>   time, &
>   nstep, &
>   theta, &
>   isconvergence, &
>   iserror, &
>   ioutput, &
>   namespout)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Kinetic reaction path 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)                  :: this            !> Type nodal chemistry variable. 

integer, intent(in)                                    :: nstep           !> Number of steps 

logical, intent(out)                                   :: iserror         !> If true there was error 

logical, intent(out)                                   :: isconvergence   !> isconvergence=true, there was convergence in the speciation 

real*8, intent(in)                                     :: time            !> Total time 

integer, intent(in), optional                          :: ioutput         !> Unit for output

real*8, intent(in)                                     :: theta 

character(len=*), intent(in), dimension(:), optional   :: namespout       !> Name of species for output
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                           :: &
> i, &
> npri, &
> nsp, &
> nsurf, &
> ntxoh, &
> istep, &
> nspout, &
> nw, &
> nwout, &
> nh, &
> nhout, &
> isp, &
> nminsp 
real*8                            :: &
> icputime1, &
> icputime2, &
> dtcputime, &
> sum, &
> dtime, &
> factoromgw
real*8, pointer                   :: &
> simin(:,:) => null (), &
> u(:) => null (), &
> ukin(:) => null (), &
> d(:) => null (), &
> cout(:,:) => null (), &
> actout(:,:) => null (), &
> timeloc(:) => null (), &
> molminsp(:,:) => null (), &
> areaminsp(:,:) => null (), &
> molaq(:) => null (), &
> molads(:) => null (), &
> molprec(:) => null (), &
> molinit(:) => null (), &
> molfin(:) => null (), &
> errorcomp(:) => null ()
integer, pointer                  :: &
> idpri(:) => null (), &
> idspout(:) => null (), &
> idminsp(:) => null ()
character(len=100), pointer       :: &
> namesploc(:) => null ()
logical                           :: &
> haveioutput, &
> havenamespout, &
> writehead
character(len=100)                :: &
> msg, &
> head 
!-------------------------------------------------------------------------
!
!>   $code
!
msg=' '
iserror=.false.
!%-------------------------------------------------------
!% Check if the chemical system is associated
!%-------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>   msg='Error, not associated chemical system'
>   goto 10
end if
!%-------------------------------------------------------
!% Check optional arguments
!%-------------------------------------------------------
haveioutput=present(ioutput)
havenamespout=present(namespout)
!%-------------------------------------------------------
factoromgw=1.0d0
nspout=0  
!%-------------------------------------------------------
call get_chem_info_ (this%pchemsys, &
>                     iserror, &
>                     numsp=nsp, &
					 numsites=ntxoh, &
					 numsurf=nsurf, &
>                     namesp=namesploc, &
					 wspindex=nw, &
					 idbase=idpri,&
					 numbase=npri, &
					 idminsp=idminsp, &
					 nminsp=nminsp)
if (iserror) goto 20
!%---------------------------------------------------------
!% Set in the chemical ssystem if the water components mass
!% balance must be considered 
!%---------------------------------------------------------
call set_iswcompbal_ (this%pchemsys,.true.,iserror)					 
if (iserror) goto 20
!%---------------------------------------------------------
if (haveioutput.and.havenamespout) then
> nwout=0
> nhout=0
> nspout=size(namespout)
> call get_sp_index_(this%pchemsys,'h+',nh)
> if (nspout>0) then
>  call check_pointer_ (idspout,nspout,.true.)
>  call check_pointer_ (cout,nspout,nstep,.true.)
>  call check_pointer_ (actout,nspout,nstep,.true.)
>  do i=1,nspout
>   call get_sp_index_(this%pchemsys,namespout(i),isp)
>   if (isp/=0) then 
>    idspout(i)=isp
>    if (isp==nw) nwout=i
>    if (isp==nh) nhout=i
>   else
>    msg='Error, not defined in the chemical system the species:'
>    call add_ (msg,namespout(i))
>    iserror=.true.
>    goto 20 
>   end if 
>  end do
> end if
end if
!%-------------------------------------------------------
if (haveioutput) then
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) 'Beginning of kinetic batch-reaction calculations'
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
	 write(ioutput,2) 'Total time:',time,'[s]' 
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) 'Initial Speciation'
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) 'Name:',this%name
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     call write_ (this,ioutput,iserror)
>     if (iserror) goto 20
> 
end if
!%-----------------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (timeloc,nstep,.true.)
call check_pointer_ (molinit,npri,.true.)
call check_pointer_ (molfin,npri,.true.)
call check_pointer_ (errorcomp,npri,.true.)
if (haveioutput.and.nminsp>0) then 
> call check_pointer_ (molminsp,nminsp,nstep,.true.)
> call check_pointer_ (areaminsp,nminsp,nstep,.true.)
end if
!%-----------------------------------------------------------------
!% Compute the initial total of mol
!%-----------------------------------------------------------------
if (haveioutput) then
> 
> call get_chem_info_ (this,iserror,molaq=molaq,molads=molads,molprec=molprec)
> 
> molinit = molaq + molads + molprec
> 
end if 
!%-----------------------------------------------------------------
!% Compute time increment 
!%-----------------------------------------------------------------
dtime=dabs(time)/real(nstep)
!%-----------------------------------------------------------------
sum=0.0d0
!%-----------------------------------------------------------------
call cpu_time (icputime1)
!%-----------------------------------------------------------------
d = this%omgwfree * ( u + (1.0d0-theta) * ukin * dtime)
!%-----------------------------------------------------------------
do istep=1,nstep
>      sum=sum+dtime
	  timeloc(istep) = sum 
>      this%alpha=this%alpha/this%omgwfree
>      call make_lin_trf_(this%pchemsys,u,this%c,0,iserror)
>      if (iserror) goto 20
>      call make_lin_trf_(this%pchemsys,ukin,this%sktrk,0,iserror)
>      if (iserror) goto 20
!%----------------------------------
!% Compute total in k 
!%----------------------------------
>      u = theta * ukin * dtime + d / this%omgwfree
	  call set_from_u_ (this,u,dtime,isconvergence,iserror)
	  if (.not.isconvergence.or.iserror) goto 20 
	  call update_min_area_ (this,iserror)
>      if (iserror) goto 20 
!%-------------------------------------------------------
!%-------------------------------------------------------
!%-------------------------------------------------------


>      d = this%omgwfree * ( u + (1.0d0-theta) * ukin * dtime)
	  
	  call make_lin_trf_(this%pchemsys,u,this%c,0,iserror)
>      if (iserror) goto 20
>      call make_lin_trf_(this%pchemsys,ukin,this%sktrk,0,iserror)
>      if (iserror) goto 20
!%-------------------------------------------------------
!> Write output (optional)
!%------------------------------------------------------- 
>      if (haveioutput) then
>        if (nminsp>0) then
>          molminsp(:,istep)=this%c(idminsp)*this%omgwfree
>          areaminsp(:,istep)=this%alpha(idminsp)
>        end if 
>        if (havenamespout.and.nspout>0) then
>          cout(:,istep)=this%c(idspout)
>          actout(:,istep)=this%c(idspout)*this%g(idspout)
>          if (nwout>0) actout(nwout,istep)=this%g(nw)
		  if (nhout>0) actout(nhout,istep)=-dlog10(this%c(nh)*this%g(nh))
>        end if
>        write(ioutput,3) '****************************************'// &
>                       '****************************************'// &
>                       '*************'
	    write(ioutput,*) 'Reaction step', istep
	    write(ioutput,2) 'Time increment',dtime,'[s]' 
>        write(ioutput,3) '----------------------------------------'// &
>                       '----------------------------------------'// &
>                       '-------------'
>        call write_ (this,ioutput,iserror)
>        if (iserror) goto 20
>      end if
!%-------------------------------------------------------
!%------------------------------------------------------- 
> 
end do
!%-------------------------------------------------------
call cpu_time (icputime2)
dtcputime=icputime2-icputime1
!%-------------------------------------------------------
!> Write the cpu time
!%-------------------------------------------------------
if (haveioutput) then
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,4) 'Total cpu time used:',dtcputime,' seconds'
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
end if 
!%--------------------------------------------------------------------
!> Write relative error in components 
!%--------------------------------------------------------------------
if (haveioutput) then
> 
> call get_chem_info_ (this,iserror,molaq=molaq,molads=molads,molprec=molprec)
> 
> head='Relative error in components'
> 
> molfin = molaq + molads + molprec
> 
> errorcomp = dabs((molinit - molfin) / molfin)
> 
> writehead=.true.
> 
> call write_ &
>      (ioutput, &
>       head, &
>       namesploc(idpri), &
>       errorcomp, &
>       npri, &
>       'error', &
>       writehead, &
>       writehead)
> 
> 
> 
end if
!%--------------------------------------------------------------------
!> Write evolution of concentration species vs. time 
!%--------------------------------------------------------------------
if (haveioutput.and.nminsp>0) then
> head='Evolution of mol of minerals vs. time'
> do istep=1,nstep
>     if (istep==1) then
>        writehead=.true.
>     else
>        writehead=.false.
>     end if
	 call write_ &
>      (ioutput, &
>       head, &
>       namesploc(idminsp), &
>       molminsp(:,istep), &
>       nminsp, &
>       timeloc(istep), &
>       writehead, &
>       writehead)
>  end do
>  
>  head='Evolution of area [m2] of minerals vs. time'
>  do istep=1,nstep
>     if (istep==1) then
>        writehead=.true.
>     else
>        writehead=.false.
>     end if
	 call write_ &
>      (ioutput, &
>       head, &
>       namesploc(idminsp), &
>       areaminsp(:,istep), &
>       nminsp, &
>       timeloc(istep), &
>       writehead, &
>       writehead)
>   end do
end if 

if (haveioutput.and.nspout>0) then
> head='Evolution of concentration species vs. time'
> do istep=1,nstep
>     if (istep==1) then
>        writehead=.true.
>     else
>        writehead=.false.
>     end if
	 call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       cout(1:nspout,istep), &
>       nspout, &
>       timeloc(istep), &
>       writehead, &
>       writehead)
> end do
!%--------------------------------------------------------------------
!> Write evolution of activity species vs. time 
!%-------------------------------------------------------------------- 
> head='Evolution of activity species vs. time'
> do istep=1,nstep
>     if (istep==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
>     call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       actout(1:nspout,istep), &
>       nspout, &
>       timeloc(istep), &
>       writehead, &
>       writehead)
> end do
end if

20 continue 
!%-------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------
call check_pointer_ (idspout,1,.false.)
call check_pointer_ (cout,1,1,.false.)
call check_pointer_ (actout,1,1,.false.)
call check_pointer_ (u,1,.false.)
call check_pointer_ (ukin,1,.false.)
call check_pointer_ (d,1,.false.)
call check_pointer_ (namesploc,1,.false.)
call check_pointer_ (idpri,1,.false.)
call check_pointer_ (timeloc,1,.false.)
call check_pointer_ (idminsp,1,.false.)
call check_pointer_ (molminsp,1,1,.false.)
call check_pointer_ (areaminsp,1,1,.false.)
call check_pointer_ (molaq,1,.false.)
call check_pointer_ (molads,1,.false.)
call check_pointer_ (molprec,1,.false.)
call check_pointer_ (molinit,1,.false.)
call check_pointer_ (molfin,1,.false.)
call check_pointer_ (errorcomp,1,.false.)
if (iserror) goto 10
!%-------------------------------------------------------
return
1 format (f10.2,a9,i5,a7)
2 format (a15,f10.2,a5)
3 format (a95)
4 format (a20,e10.4,a8)
> 
10 continue 
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: reaction_path_'
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
subroutine compute_f_df_rmrmt_nch &
>  (this, &
>   f, &
>   df, &
>   ndim, &
>   boxk, &
>   mobnodek1, &
>   mobnodek, &
>   theta, &
>   alpha, &
>   phi, &
>   dtime, &
>   isconvergence, &
>   tolunk, &
>   tolres, &
>   mxiter, &
>   mxiterminchange, &
>   method, &
>   iserror, &
>   fsp)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Compute f and df according to method 1, i.e., assuming exponential time behavior in the evolution of immobile concentrations. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)                  :: this            !> Immobile box at k+1. Type nodal chemistry variable. 

type(t_nodalchemistry), intent(inout)                  :: boxk            !> Immobile box at k. Type nodal chemistry variable.

type(t_nodalchemistry), intent(in)                     :: mobnodek1       !> mobile node at k+1. Type nodal chemistry variable (box) 

type(t_nodalchemistry), intent(in)                     :: mobnodek        !> mobile node at k. Type nodal chemistry variable (box) 

real*8, pointer, dimension(:)                          :: f               !> F [ndim]

real*8, pointer, dimension(:,:)                        :: df              !> dF [ndim,ndim]

integer, intent(out)                                   :: ndim            !> dimensions of F

real*8, intent(in)                                     :: dtime           !> Time increment (s)

real*8, intent(in)                                     :: theta           !> Temporal weight

real*8, intent(in)                                     :: alpha           !> Alpha for box 

real*8, intent(in)                                     :: phi             !> Porosity for box 

real*8, intent(in)                                     :: tolunk          !> Tolerance in the unknowns for Newton-Raphson

real*8, intent(in)                                     :: tolres          !> Tolerance in the residual for Newton-Raphson

integer, intent(in)                                    :: mxiter          !> Maximum number of iterations for Newton-Raphson

integer, intent(in)                                    :: mxiterminchange !> Maximum number of iterations for changes in mineral presence

logical, intent(out)                                   :: isconvergence   !> isconvergence=true, then there was convergence

logical, intent(out)                                   :: iserror         !> iserror=true, then there was an error

integer, intent(in)                                    :: method          !> Numerical method to solve immobile region. If method = 1 then an exponential , if method = 2 Newton-Raphson is employed

real*8, pointer, dimension(:), optional                :: fsp               

!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)              :: msg
integer, parameter              :: &
> method_1=1, &
> method_2=2
!%-------------------------------------------------------------------------
!% Initialice variables 
!%-------------------------------------------------------------------------
msg=''
iserror=.false. 
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
select case(method)

case(method_1)

>    call compute_f_df_rmrmt_1_nch(this,f,df,ndim,boxk,mobnodek1,mobnodek,theta,alpha,phi,dtime,isconvergence,tolunk,tolres,mxiter,mxiterminchange,iserror,fsp=fsp)

case(method_2)

>    call compute_f_df_rmrmt_2_nch(this,f,df,ndim,boxk,mobnodek1,mobnodek,theta,alpha,phi,dtime,isconvergence,tolunk,tolres,mxiter,mxiterminchange,iserror,fsp=fsp)

case default

>    iserror = .true.
>    msg = 'Error, reactive multi-rate mass transfer method not implemented'

end select

!%-----------------------------------------------------------------
20 continue  
if (iserror) goto 10 
!%-------------------------------------------------------
return
> 
10 continue 
print *,'****************************'
print *,'Nodal Chemistry:'
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
subroutine write_backup_nch &
>  (this, &
>   iobackup, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Write state variables in the backup file.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)      :: this     !> Type nodal chemistry variable. 

integer, intent(in)                     :: iobackup !> Output unit 

logical, intent(out)                    :: iserror  !> iserror=true, then there was an error 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                 :: &
> nsp, &
> npri, &
> nsurf, &
> nsite, &
> nmin, &
> i, &
> j 
integer, pointer                        :: &
> idpri(:) => null (), &
> idmin(:) => null () 
character(len=100)                      :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!

!%-----------------------------------------------------------
msg=''
iserror=.false.
!%-----------------------------------------------------------
!% Check if the chemical system is associated
!%-----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%--------------------------------------------------------------
!% Request chemical information to chemical system 
!%--------------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,numsp=nsp,numsites=nsite, &
>                     numsurf=nsurf,idminsp=idmin,nminsp=nmin, &
>                     idbase=idpri,numbase=npri)
!%--------------------------------------------------------------
!% 
!%--------------------------------------------------------------
write (iobackup,'(/,7e22.14)',ADVANCE='no') this%volnch, this%volgas, &
>       this%omgwfree, this%pgas, this%pliq, this%liqdens, this%temp
!%--------------------------------------------------------------
!% 
!%--------------------------------------------------------------
write (iobackup,'(i5)',ADVANCE='no') this%hashcompz
!%--------------------------------------------------------------
!% Write concentrations of primary species 
!%--------------------------------------------------------------
write (iobackup,'(<npri>e22.14)',ADVANCE='no') (this%c(idpri(i)),i=1,npri)
!%--------------------------------------------------------------
!% Write concentrations of mineral species 
!%--------------------------------------------------------------
write (iobackup,'(<nmin>e22.14)',ADVANCE='no') (this%c(idmin(i)),i=1,nmin)
!%--------------------------------------------------------------
!%  
!%--------------------------------------------------------------
write (iobackup,'(<nsp>e22.14)',ADVANCE='no') (this%alpha(i),i=1,nsp)
write (iobackup,'(<nsp>e22.14)',ADVANCE='no') (this%alpha0(1,i),i=1,nsp)
write (iobackup,'(<nsp>e22.14)',ADVANCE='no') (this%alpha0(2,i),i=1,nsp)
do j=1,nsurf
> write (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%txoh(i,j),i=1,nsite)
end do 
!%---------------------------------------------------------------
!% 
!%---------------------------------------------------------------
do j=1,nsurf
> write (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%txoh0(i,j),i=1,nsite)
end do
do j=1,nsurf
> write (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%capint(i,j),i=1,nsite)
end do
do j=1,nsurf
> write (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%capext(i,j),i=1,nsite)
end do
do j=1,nsurf
> write (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%spsurfarea(i,j),i=1,nsite)
end do
!%---------------------------------------------------------------
20 continue 
call check_pointer_ (idpri,1,.false.)
call check_pointer_ (idmin,1,.false.)
if (iserror) goto 10 
!%---------------------------------------------------------------
!%---------------------------------------------------------------
return

> 
10 continue 
print *,'*******************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: write_'
print *, msg
print *,'*******************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_backup_nch &
>  (this, &
>   iobackup, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Read state variables in the backup file.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)   :: this     !> Type nodal chemistry variable. 

integer, intent(in)                     :: iobackup !> Output unit 

logical, intent(out)                    :: iserror  !> iserror=true, then there was an error 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                 :: &
> nsp, &
> npri, &
> nsurf, &
> nsite, &
> nmin, &
> i, &
> j 
real*8                                  :: &
> faccap, &
> dummy 
integer, pointer                        :: &
> idpri(:) => null (), &
> idmin(:) => null () 
real*8, pointer                         :: &
> alpha(:) => null (), &
> c(:) => null () 
character(len=100)                      :: &
> msg 
logical                                 :: &
> isanomalous, &
> isupmxitergam  
!-------------------------------------------------------------------------
!
!>   $code
!
!%-----------------------------------------------------------
msg=''
iserror=.false.
!%-----------------------------------------------------------
!% Check if the chemical system is associated
!%-----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%--------------------------------------------------------------
!% Request information to chemical system 
!%--------------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,numsp=nsp,numsites=nsite, &
>                     numsurf=nsurf,idminsp=idmin,nminsp=nmin, &
>                     idbase=idpri,numbase=npri)
if (iserror) then
> msg='Error when call get_chem_info_ service'
> goto 20                    
end if 
!%--------------------------------------------------------------
!% 
!%--------------------------------------------------------------
read (iobackup,'(/,7e22.14)',ADVANCE='no') this%volnch, this%volgas, &
>       this%omgwfree, this%pgas, this%pliq, this%liqdens, this%temp 
!%--------------------------------------------------------------       
!% 
!%--------------------------------------------------------------       
read (iobackup,'(i5)',ADVANCE='no') this%hashcompz                        
!%--------------------------------------------------------------
!% Write concentrations of primary species 
!%--------------------------------------------------------------
read (iobackup,'(<npri>e22.14)',ADVANCE='no') (this%c(idpri(i)),i=1,npri)
!%--------------------------------------------------------------
!% Write concentrations of mineral species 
!%--------------------------------------------------------------
read (iobackup,'(<nmin>e22.14)',ADVANCE='no') (this%c(idmin(i)),i=1,nmin)
!%--------------------------------------------------------------
!%  
!%--------------------------------------------------------------
read (iobackup,'(<nsp>e22.14)',ADVANCE='no') (this%alpha(i),i=1,nsp)
read (iobackup,'(<nsp>e22.14)',ADVANCE='no') (this%alpha0(1,i),i=1,nsp)
read (iobackup,'(<nsp>e22.14)',ADVANCE='no') (this%alpha0(2,i),i=1,nsp)
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
do j=1,nsurf
> read (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%txoh(i,j),i=1,nsite)
end do 
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
do j=1,nsurf
> read (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%txoh0(i,j),i=1,nsite)
end do 
do j=1,nsurf
> read (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%capint(i,j),i=1,nsite)
end do 
do j=1,nsurf
> read (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%capext(i,j),i=1,nsite)
end do 
do j=1,nsurf
> read (iobackup,'(<nsite>e22.14)',ADVANCE='no') (this%spsurfarea(i,j),i=1,nsite)
end do 
!%---------------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------------
call check_pointer_ (alpha,nsp,.true.)
call check_pointer_ (c,nsp,.true.)
alpha=this%alpha/this%omgwfree
c=this%c
!%--------------------------------------------------------------
!% Request information to chemical system 
!%--------------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,idbase=idpri,numbase=npri, &
>                     hashcompz=this%hashcompz)
if (iserror) then
> msg='Error when call get_chem_info_ service'
> goto 20                    
end if    
!%-----------------------------------------------------------
!% Compute the capilar factor 
!%-----------------------------------------------------------
call get_chem_info_ (this,iserror,faccap=faccap)
if (iserror) goto 20
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
call update_ (this%pchemsys,this%temp,iserror)
if (iserror) goto 20
!%-----------------------------------------------------------
!%  Call specia_from_cpri_ in the chemicla system object
!%-----------------------------------------------------------
if (this%isderstored) then 
> call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	c, &  
>    alpha, &
>    nsp, &
>    c(idpri), &
>    npri, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    nsite, &
>    nsurf, &
>    0.0d0, &
	this%volgas, &
>    this%ionstr, &
>    this%hashcompz, &
	faccap, &    !> Capillary correction for water activity 
>    isanomalous, &
>    isupmxitergam, &
>    .false., &
>    iserror, &
>    sktrk=this%sktrk, &
>    dc=this%dc, &
>    dsktrk=this%dsktrk)
else
>  call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	c, &
>    alpha, &
>    nsp, &
>    c(idpri), &
>    npri, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    nsite, &
>    nsurf, &
>    0.0d0, &
	this%volgas, &
>    this%ionstr, &
>    this%hashcompz, &
	faccap, &    !> Capillary correction for water activity 
>    isanomalous, &
>    isupmxitergam, &
>    .false., &
>    iserror, &
>    sktrk=this%sktrk)
> end if 

if (iserror) goto 20  
20 continue 
!%---------------------------------------------------------------
!% Deallocate local pointers 
!%---------------------------------------------------------------
call check_pointer_ (idpri,1,.false.)
call check_pointer_ (idmin,1,.false.)
call check_pointer_ (alpha,1,.false.)
call check_pointer_ (c,1,.false.)
if (iserror) goto 10 
!%---------------------------------------------------------------
!%---------------------------------------------------------------
return

> 
10 continue 
print *,'**********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: read_backup_'
print *, msg
print *,'**********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_f_df_rmrmt_1_nch &
>  (this, &
>   f, &
>   df, &
>   ndim, &
>   boxk, &
>   mobnodek1, &
>   mobnodek, &
>   theta, &
>   alpha, &
>   phi, &
>   dtime, &
>   isconvergence, &
>   tolunk, &
>   tolres, &
>   mxiter, &
>   mxiterminchange, &
>   iserror, &
>   fsp)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Compute f and df according to method 1, i.e., assuming exponential time behavior in the evolution of immobile concentrations. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)                  :: this            !> Immobile box at k+1. Type nodal chemistry variable. 

type(t_nodalchemistry), intent(inout)                  :: boxk            !> Immobile box at k. Type nodal chemistry variable.

type(t_nodalchemistry), intent(in)                     :: mobnodek1       !> mobile node at k+1. Type nodal chemistry variable (box) 

type(t_nodalchemistry), intent(in)                     :: mobnodek        !> mobile node at k. Type nodal chemistry variable (box) 

real*8, pointer, dimension(:)                          :: f               !> F [ndim]

real*8, pointer, dimension(:,:)                        :: df              !> dF [ndim,ndim]

integer, intent(out)                                   :: ndim            !> dimensions of F

real*8, intent(in)                                     :: dtime           !> Time increment (s)

real*8, intent(in)                                     :: theta           !> Temporal weight

real*8, intent(in)                                     :: alpha           !> Alpha for box 

real*8, intent(in)                                     :: phi             !> Porosity for box 

real*8, intent(in)                                     :: tolunk          !> Tolerance in the unknowns for Newton-Raphson

real*8, intent(in)                                     :: tolres          !> Tolerance in the residual for Newton-Raphson

integer, intent(in)                                    :: mxiter          !> Maximum number of iterations for Newton-Raphson

integer, intent(in)                                    :: mxiterminchange !> Maximum number of iterations for changes in mineral presence

logical, intent(out)                                   :: isconvergence   !> isconvergence=true, then there was convergence

logical, intent(out)                                   :: iserror         !> iserror=true, then there was an error

real*8, pointer, dimension(:), optional                :: fsp             !> F [nsp]
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
type(t_nodalchemistry)          :: boxktheta   !> Immobile box at k+theta. Type nodal chemistry variable.

character(len=100)              :: msg
real*8, pointer, dimension(:,:) :: dc => null ()
real*8, pointer, dimension(:,:) :: dcmob => null ()
real*8, pointer, dimension(:,:) :: du => null ()
real*8, pointer, dimension(:,:) :: dutotdc1 => null ()
real*8, pointer, dimension(:,:) :: dsktrk => null ()
real*8, pointer, dimension(:)   :: sktrk => null ()
real*8, pointer, dimension(:)   :: c => null ()
real*8, pointer, dimension(:)   :: caq => null ()
real*8, pointer, dimension(:)   :: u => null ()
real*8, pointer, dimension(:)   :: c1k => null ()
real*8, pointer, dimension(:)   :: delc1 => null ()
real*8, pointer, dimension(:)   :: wsp => null ()
real*8, pointer, dimension(:,:) :: Abox => null ()
real*8, pointer, dimension(:,:) :: expdtthetaA => null ()
real*8, pointer, dimension(:,:) :: expdtA => null ()
real*8, pointer, dimension(:,:) :: Ipri => null ()
real*8, pointer, dimension(:)   :: bbox => null ()
real*8, pointer, dimension(:)   :: avbox => null ()
real*8, pointer, dimension(:)   :: c1ktheta => null ()
real*8, pointer, dimension(:)   :: c1k1 => null ()
integer, pointer, dimension(:)  :: indx => null ()
real*8, pointer, dimension(:)   :: setre => null ()
real*8, pointer, dimension(:,:) :: dc1imdc1m => null ()
integer, pointer, dimension(:)  :: ipiv => null ()
integer                         :: iter
integer                         :: iterminchange
integer                         :: npri
integer                         :: nprimob
integer                         :: nsp
integer                         :: i,j
integer                         :: info
integer, parameter              :: ideg = 6
integer                         :: lwsp
integer                         :: iexph
integer                         :: ns
integer                         :: iflag
logical                         :: isanomalous
logical                         :: isupmxitergam
logical                         :: iscompzchanged
logical                         :: havefsp
> 
!%-------------------------------------------------------------------------
!% Initialice variables 
!%-------------------------------------------------------------------------
msg=''
iserror=.false.
isconvergence=.false. 
!%-------------------------------------------------------------------------
!% Check optional arguments 
!%-------------------------------------------------------------------------
havefsp=present(fsp)
!%-------------------------------------------------------------------------
!% Check if the chemical system is associated in the nodal chemistry 
!% objects  
!%-----------------------------------------------------------------
if (.not.associated(this%pchemsys).or..not.associated(boxk%pchemsys)) then
>   msg='Error, not associated chemical system'
>   goto 10
end if
!%-----------------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,numsp=nsp)
!%-----------------------------------------------------------------
!%
!%-----------------------------------------------------------------
if (havefsp) then 
> call check_pointer_ (fsp,nsp,.true.)
end if 
!%-----------------------------------------------------------------
!% Create boxtheta 
!%-----------------------------------------------------------------
call create_ (boxktheta)


iterminchange = 0
iscompzchanged=.true.

do while (iscompzchanged.and.iterminchange.le.mxiterminchange)

>  isconvergence=.false.
>  
>  iterminchange = iterminchange + 1

!%-----------------------------------------------------------------
!% Dimension of vector and matrices
!%-----------------------------------------------------------------
>  call get_chem_info_ (this,iserror,npri=npri)
>  call check_pointer_ (indx,npri,.true.)
>  call check_pointer_ (ipiv,npri,.true.)
>  call check_pointer_ (delc1,npri,.true.)
>  call check_pointer_ (dutotdc1,npri,npri,.true.)
>  call check_pointer_ (Abox,npri,npri,.true.)
>  call check_pointer_ (Ipri,npri,npri,.true.)
>  call check_pointer_ (expdtthetaA,npri,npri,.true.)
>  call check_pointer_ (expdtA,npri,npri,.true.)
>  call check_pointer_ (bbox,npri,.true.)
>  call check_pointer_ (avbox,npri,.true.)
>  call check_pointer_ (c1ktheta,npri,.true.)
>  call check_pointer_ (c1k1,npri,.true.)

>  
>  lwsp = 4*npri*npri+ideg+1
>  
>  call check_pointer_ (wsp,lwsp,.true.)

!%-----------------------------------------------------------------
!% Iterative processes to calculate c1ktheta
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Initialization for Picard
!%-----------------------------------------------------------------
>  iter=0
!%-----------------------------------------------------------------  
!% We get c1,imm at time tk
!%-----------------------------------------------------------------
>  call get_chem_info_ (boxk,iserror,npri=npri,cpri=c1k)
!%-----------------------------------------------------------------
!% 
!%-----------------------------------------------------------------  
>  boxktheta=boxk
>  

>  call get_chem_info_ (boxktheta,iserror,cpri=c1ktheta)
>  if (iserror) goto 20
!%-----------------------------------------------------------------
!% Start Newton-Raphson iterative process 
!%-----------------------------------------------------------------
>  do while (.not.isconvergence.and.iter.le.mxiter)
>    
	iter = iter + 1

!%------------------------------------------------------------------------------------------
!%  We first calculate concentration c1imk+theta at iteration i
!%  Then we set intermediate nodal chemistry for immobile box at k+theta, boxktheta
!%------------------------------------------------------------------------------------------
>     
>     delc1 = c1ktheta
>    
!%------------------------------------------------------------------------------------------
!%  Now we calculate matrix A and vectors a and b that define the ODE for primary c1,imm, at k+theta
!%------------------------------------------------------------------------------------------
!
>     Abox  = 0.d0
>     bbox  = 0.d0
>     avbox = 0.d0
>     
>     call get_chem_info_ (boxktheta,iserror,c=c,caq=caq,sktrk=sktrk,dc=dc,dcmob=dcmob,dsktrk=dsktrk)
>     if (iserror) goto 20

>     call make_lin_trf_  (boxktheta%pchemsys,dutotdc1,dc,this%hashcompz,iserror)
>     if (iserror) goto 20
>     
>     call make_lin_trf_  (boxktheta%pchemsys,u,caq,this%hashcompz,iserror)
>     if (iserror) goto 20
>     bbox = bbox - alpha*u
>     
>     call make_lin_trf_  (boxktheta%pchemsys,du,dcmob(1:nsp,1:npri),this%hashcompz,iserror)
>     if (iserror) goto 20
>     Abox = Abox + alpha*du
>     
>     call make_lin_trf_  (boxktheta%pchemsys,u,sktrk,this%hashcompz,iserror)
>     if (iserror) goto 20
>     
>     bbox = bbox + u
>     
>     call make_lin_trf_  (boxktheta%pchemsys,du,dsktrk,this%hashcompz,iserror)
>     if (iserror) goto 20
>     
>     Abox = Abox - du
>     
>     !%-------------------------------------------
	 !% Dependence of mobile node at k+1 (mobnodek1)
>     !%-------------------------------------------
>     call get_chem_info_ (mobnodek1,iserror,caq=caq)
>     if (iserror) goto 20
>     
>     call make_lin_trf_  (mobnodek1%pchemsys,u,caq,this%hashcompz,iserror)
>     if (iserror) goto 20
>     avbox = avbox + u/dtime
>     !%-------------------------------------------
	 !% Dependence of mobile node at k (mobnodek)
	 !%-------------------------------------------
>     call get_chem_info_ (mobnodek,iserror,caq=caq)
>     if (iserror) goto 20
>     
>     call make_lin_trf_  (mobnodek%pchemsys,u,caq,this%hashcompz,iserror)
>     if (iserror) goto 20
>     bbox = bbox + alpha*u
>     avbox = avbox - u/dtime
>     avbox = alpha*avbox

>     !%--------------------------------------
>     !% Final algebra to get A, a and b
>     !%--------------------------------------
>     !% Calculates the LU decomposition of dutotdc1
	 !%--------------------------------------
>     call f07adf(npri,npri,dutotdc1,npri,indx,info)
>     if (info/=0) then 
	    iserror=.true.
		goto 20
	 end if
>     call f07aef('n',npri,1,dutotdc1,npri,indx,avbox,npri,info)
>     if (info/=0) then 
	    iserror=.true.
		goto 20
	 end if
>     bbox = bbox + matmul(Abox,c1ktheta)

>     call f07aef('n',npri,1,dutotdc1,npri,indx,bbox,npri,info)
>     if (info/=0) then 
	    iserror=.true.
		goto 20
	 end if
>     call f07aef('n',npri,npri,dutotdc1,npri,indx,Abox,npri,info)
	 if (info/=0) then 
	    iserror=.true.
		goto 20
	 end if

>     Abox = -Abox
>     
>     call dgpadm(ideg,npri,dtime*theta,Abox,npri,wsp,lwsp,ipiv,iexph,ns,iflag)
>     
>     do i=1,npri
>        do j=1,npri
>            expdtthetaA(i,j) = wsp(iexph+i-1+npri*j-npri)
>        end do
>     end do
>     
>     Ipri = 0.0d0
>     do i=1,npri
>       Ipri(i,i) = 1.d0
>     end do
>     
>     call f07adf(npri,npri,Abox,npri,indx,info)
	 if (info/=0) then 
	    iserror=.true.
		goto 20
	 end if
>     call f07aef('n',npri,npri,Abox,npri,indx,Ipri,npri,info)
>     if (info/=0) then 
	    iserror=.true.
		goto 20
	 end if
>     !%---------------------------------------
	 !Now the inverse of A is stored in Abox
	 !%---------------------------------------
>     Abox = Ipri
>     
>     Ipri = 0.d0
>     do i=1,npri
>       Ipri(i,i) = 1.d0
>     end do
>     
>     
!%-----------------------------------------------------------------
!% Update c1,imm, at k+theta and check convergence
!%-----------------------------------------------------------------
>          
>     
>     c1ktheta = matmul(expdtthetaA,c1k+matmul(Abox,matmul(Abox,avbox)+bbox))
>     c1ktheta = c1ktheta - matmul(Abox,matmul(Abox,avbox)+theta*dtime*avbox+bbox)
>   
>     delc1 = c1ktheta - delc1
>     
>     call set_from_cpri_ &
>    (boxktheta, &
>     boxk, &
>     c1ktheta, &
>     dtime*theta, &
>     .false., &
>     isanomalous, &
>     isupmxitergam, &
>     iserror)
>     if (iserror) goto 20
>     if (isanomalous.or.isupmxitergam) then
>       goto 20
>     endif

>     isconvergence = (maxval(dabs(delc1/c1ktheta)).lt.tolunk ) 

>  end do
>  
>  if (.not.isconvergence) then
>     goto 20
>  end if

!%-----------------------------------------------------------------
!% Calculate concentration of primary species at k+1
!%-----------------------------------------------------------------
!
>  call dgpadm(ideg,npri,dtime,Abox,npri,wsp,lwsp,ipiv,iexph,ns,iflag)

>  do i=1,npri
>     do j=1,npri
>         expdtA(i,j) = wsp(iexph+i-1+npri*j-npri)
>     end do
>  end do
>     
>  c1k1 = (c1ktheta - (1.d0 - theta)*c1k)/theta
>     
>  call set_from_cpri_ &
>    (this, &
>     boxk, &
>     c1k1, &
>     dtime, &
>     .false., &
>     isanomalous, &
>     isupmxitergam, &
>     iserror)
>     if (iserror) goto 20
>     if (isanomalous.or.isupmxitergam) then
>       goto 20
>     endif
!%-----------------------------------------------------------------
!% Calculate concentration of minerals at k+1
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Dimension of Setre
!%-----------------------------------------------------------------
>  call check_pointer_ (setre,nsp,.true.)
>  !setre = 0.0
!%-----------------------------------------------------------------
!% Dependence of immobile box at k+1 (this)
!%-----------------------------------------------------------------
>  call get_chem_info_ (this,iserror,c=c)
>  if (iserror) goto 20
>  setre = setre + c/dtime 
!%-----------------------------------------------------------------
!% Dependence of immobile box at k (boxk)
!%-----------------------------------------------------------------
>  call get_chem_info_ (boxk,iserror,c=c)
>  if (iserror) goto 20

>  
>  setre = setre - c/dtime
!%-----------------------------------------------------------------  
!% Dependence of immobile box at k+theta (boxktheta)
!%-----------------------------------------------------------------
>  call get_chem_info_ (boxktheta,iserror,c=c,caq=caq,sktrk=sktrk)
>  if (iserror) goto 20
> 
>  setre = setre + alpha*caq - sktrk
!%-----------------------------------------------------------------
!% Dependence of mobile node at k+1 (mobnodek1)
!%-----------------------------------------------------------------
>  call get_chem_info_ (mobnodek1,iserror,caq=caq)
>  if (iserror) goto 20
>  
>  setre = setre - alpha*theta*caq

!%-----------------------------------------------------------------
!% Dependence of mobile node at k (mobnodek)
!%-----------------------------------------------------------------
>  call get_chem_info_ (mobnodek,iserror,caq=caq)
>  if (iserror) goto 20
>  
>  setre = setre - alpha*(1.d0-theta)*caq
!%-----------------------------------------------------------------
!% The unit of setre [mol/s]
!%-----------------------------------------------------------------
>  setre = phi*setre 
!%-----------------------------------------------------------------
>  call set_from_setre_(this,boxk,iscompzchanged,setre,nsp,dtime,iserror)
>  if (iscompzchanged) then
>    call update_derivatives_ (this,iserror)
	if (iserror) goto 20 
	call set_ (boxk,iserror,hashcompz=this%hashcompz)
	if (iserror) goto 20
>  end if
>  
end do
!%-----------------------------------------------------------------
!% Maximum number of mxiterminchange has been reached
!%-----------------------------------------------------------------
if(iscompzchanged)then
>  isconvergence=.false.
>  goto 20
end if
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Calculate immobile-mobile flux, F
!%-----------------------------------------------------------------
!% We resize if necessary due to changes in minerals
!%-----------------------------------------------------------------
call get_chem_info_ (mobnodek1,iserror,npri=ndim)
call check_pointer_ (f,ndim,.true.)
call check_pointer_ (df,ndim,ndim,.true.)
!%-----------------------------------------------------------------
!% Dependence of immobile box at k+theta (boxktheta)
!%-----------------------------------------------------------------
call get_chem_info_ (boxktheta,iserror,caq=caq)
if (iserror) goto 20
call make_lin_trf_  (boxktheta%pchemsys,u,caq,mobnodek1%hashcompz,iserror)
if (iserror) goto 20 
!%-----------------------------------------------------------------
f = f - u
if (havefsp) then
> fsp = fsp - caq
end if
!%-----------------------------------------------------------------  
!% Dependence of mobile node at k+1 (mobnodek1)
!%-----------------------------------------------------------------
call get_chem_info_ (mobnodek1,iserror,caq=caq)
if (iserror) goto 20
call make_lin_trf_  (mobnodek1%pchemsys,u,caq,mobnodek1%hashcompz,iserror)
if (iserror) goto 20
!%-----------------------------------------------------------------
f = f + theta*u
if (havefsp) then
> fsp = fsp + theta*caq
end if
!%-----------------------------------------------------------------
!% Dependence of mobile node at k (mobnodek)
!%-----------------------------------------------------------------
call get_chem_info_ (mobnodek,iserror,caq=caq)
if (iserror) goto 20
call make_lin_trf_  (mobnodek%pchemsys,u,caq,mobnodek1%hashcompz,iserror)
if (iserror) goto 20
f = f + (1.d0-theta)*u
f = phi*alpha*f
if (havefsp) then
> fsp = fsp + (1.0d0-theta)*caq
> fsp = phi*alpha*fsp
end if
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Calculate derivative of immobile-mobile flux, dF
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Dependence of mobile node at k+1 (mobnodek1)
!%-----------------------------------------------------------------
call get_chem_info_ (mobnodek1,iserror,npri=nprimob,numsp=nsp,caq=caq,dcmob=dcmob)
if (iserror) goto 20
call make_lin_trf_  (mobnodek1%pchemsys,du,dcmob(1:nsp,1:nprimob),this%hashcompz,iserror)
if (iserror) goto 20
>    
call check_pointer_ (dc1imdc1m,npri,nprimob,.true.)

call f07aef('n',npri,npri,dutotdc1,npri,indx,du,npri,info)
if (info/=0) then 
>  iserror=.true.
>  goto 20
end if

dc1imdc1m = (alpha/dtime)*matmul(matmul(expdtthetaA-Ipri,Abox)-theta*dtime*Ipri,matmul(Abox,du))

call make_lin_trf_  (mobnodek1%pchemsys,du,dcmob(1:nsp,1:nprimob),mobnodek1%hashcompz,iserror)
if (iserror) goto 20

df = theta*du


call get_chem_info_ (boxktheta,iserror,npri=npri,numsp=nsp,caq=caq,dcmob=dcmob)
if (iserror) goto 20

call make_lin_trf_  (boxktheta%pchemsys,du,dcmob(1:nsp,1:npri),mobnodek1%hashcompz,iserror)
if (iserror) goto 20


df = df - matmul(du,dc1imdc1m)
df = phi*alpha*df
!%-----------------------------------------------------------------
20 continue  
!%-----------------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dcmob,1,1,.false.)
call check_pointer_ (du,1,1,.false.)
call check_pointer_ (dutotdc1,1,1,.false.)
call check_pointer_ (dsktrk,1,1,.false.)
call check_pointer_ (sktrk ,1,.false.)
call check_pointer_ (c ,1,.false.)
call check_pointer_ (caq ,1,.false.)
call check_pointer_ (u ,1,.false.)
call check_pointer_ (c1k ,1,.false.)
call check_pointer_ (c1ktheta,1,.false.)
call check_pointer_ (c1k1,1,.false.)
call check_pointer_ (delc1 ,1,.false.)
call check_pointer_ (setre ,1,.false.)
call check_pointer_ (indx ,1,.false.)
call check_pointer_ (wsp ,1,.false.)
call check_pointer_ (expdtA,1,1,.false.)
call check_pointer_ (expdtthetaA,1,1,.false.)
call check_pointer_ (Abox,1,1,.false.)
call check_pointer_ (Ipri,1,1,.false.)  
call check_pointer_ (bbox,1,.false.)
call check_pointer_ (avbox,1,.false.)

call destroy_ (boxktheta)


if (iserror) goto 10 
!%-------------------------------------------------------
return
> 
10 continue 
print *,'****************************'
print *,'Nodal Chemistry:'
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
subroutine compute_f_df_rmrmt_2_nch &
>  (this, &
>   f, &
>   df, &
>   ndim, &
>   boxk, &
>   mobnodek1, &
>   mobnodek, &
>   theta, &
>   alpha, &
>   phi, &
>   dtime, &
>   isconvergence, &
>   tolunk, &
>   tolres, &
>   mxiter, &
>   mxiterminchange, &
>   iserror, &
>   fsp)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Compute f and df (See equations....). 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)                  :: this            !> Immobile box at k+1. Type nodal chemistry variable. 

type(t_nodalchemistry), intent(inout)                  :: boxk            !> Immobile box at k. Type nodal chemistry variable.

type(t_nodalchemistry), intent(in)                     :: mobnodek1       !> mobile node at k+1. Type nodal chemistry variable (box) 

type(t_nodalchemistry), intent(in)                     :: mobnodek        !> mobile node at k. Type nodal chemistry variable (box) 

real*8, pointer, dimension(:)                          :: f               !> F [ndim]

real*8, pointer, dimension(:,:)                        :: df              !> dF [ndim,ndim]

integer, intent(out)                                   :: ndim            !> dimensions of F 

real*8, intent(in)                                     :: dtime           !> Time increment (s)

real*8, intent(in)                                     :: theta           !> Temporal weight

real*8, intent(in)                                     :: alpha           !> Alpha for box 

real*8, intent(in)                                     :: phi             !> Porosity for box 

real*8, intent(in)                                     :: tolunk          !> Tolerance in the unknowns for Newton-Raphson

real*8, intent(in)                                     :: tolres          !> Tolerance in the residual for Newton-Raphson

integer, intent(in)                                    :: mxiter          !> Maximum number of iterations for Newton-Raphson

integer, intent(in)                                    :: mxiterminchange !> Maximum number of iterations for changes in mineral presence

logical, intent(out)                                   :: isconvergence   !> isconvergence=true, then there was convergence

logical, intent(out)                                   :: iserror         !> iserror=true, then there was an error

real*8, pointer, dimension(:), optional                :: fsp
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)              :: msg
real*8, pointer, dimension(:,:) :: dc => null ()
real*8, pointer, dimension(:,:) :: dcmob => null ()
real*8, pointer, dimension(:,:) :: du => null ()
real*8, pointer, dimension(:,:) :: dsktrk => null ()
real*8, pointer, dimension(:)   :: sktrk => null ()
real*8, pointer, dimension(:)   :: c => null ()
real*8, pointer, dimension(:)   :: caq => null ()
real*8, pointer, dimension(:)   :: u => null ()
real*8, pointer, dimension(:)   :: c1 => null ()
real*8, pointer, dimension(:)   :: delc1 => null ()
real*8, pointer, dimension(:,:) :: jacob => null ()
real*8, pointer, dimension(:)   :: res => null ()
integer, pointer, dimension(:)  :: indx => null ()
real*8, pointer, dimension(:)   :: setre => null ()
real*8, pointer, dimension(:,:) :: dc1imdc1m => null ()
integer                         :: iter
integer                         :: iterminchange
integer                         :: npri
integer                         :: nprimob
integer                         :: nsp
integer                         :: i
integer                         :: info 
logical                         :: isanomalous
logical                         :: isupmxitergam
logical                         :: iscompzchanged
logical                         :: havefsp
> 
!%-------------------------------------------------------------------------
!% Initialice variables 
!%-------------------------------------------------------------------------
msg=''
iserror=.false.
isconvergence=.false.
!%-------------------------------------------------------------------------
!%  
!%-------------------------------------------------------------------------
havefsp=present(fsp)
!%-------------------------------------------------------------------------
!% Check if the chemical system is associated in the nodal chemistry 
!% objects  
!%-----------------------------------------------------------------
if (.not.associated(this%pchemsys).or..not.associated(boxk%pchemsys)) then
>   msg='Error, not associated chemical system'
>   goto 10
end if

iterminchange = 0
iscompzchanged=.true.



do while (iscompzchanged.and.iterminchange.le.mxiterminchange)
>  
>  isconvergence=.false. 
>  
>  iterminchange = iterminchange + 1

!%-----------------------------------------------------------------
!% Dimension of residual and jacobian
!%-----------------------------------------------------------------
>  call get_chem_info_ (this,iserror,npri=npri,numsp=nsp)
>  call check_pointer_ (res,npri,.true.)
>  call check_pointer_ (jacob,npri,npri,.true.)
>  call check_pointer_ (indx,npri,.true.)
>  call check_pointer_ (delc1,npri,.true.)
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Newton-Raphson iterative processes
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Initialization for NR
!%-----------------------------------------------------------------
>  iter=0
>  call get_chem_info_ (boxk,iserror,cpri=c1)
!%-----------------------------------------------------------------
!% Initialice with box k+theta 
!%-----------------------------------------------------------------  
>  !cprovi this=boxk
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
>  do while (.not.isconvergence.and.iter.le.mxiter)
>    iter = iter + 1
!%-----------------------------------------------------------------
!% Calculate residual and jacobian matrix
!%-----------------------------------------------------------------
>    res=0.0d0
>    jacob=0.0d0
!%-----------------------------------------------------------------
!% Dependence of immobile box at k+1 (this)
!%-----------------------------------------------------------------
>     call get_chem_info_ (this,iserror,c=c,caq=caq,sktrk=sktrk,dc=dc,dcmob=dcmob,dsktrk=dsktrk)
>     if (iserror) goto 20
> 
>     call make_lin_trf_  (this%pchemsys,u,c,this%hashcompz,iserror)
>     if (iserror) goto 20
>     res=res-u/dtime

>     call make_lin_trf_  (this%pchemsys,du,dc,this%hashcompz,iserror)
>     if (iserror) goto 20
>     jacob=jacob+du/dtime

>     call make_lin_trf_  (this%pchemsys,u,caq,this%hashcompz,iserror)
>     if (iserror) goto 20
>     res=res-alpha*theta*u

>     call make_lin_trf_  (this%pchemsys,du,dcmob(1:nsp,1:npri),this%hashcompz,iserror)
>     if (iserror) goto 20
>     jacob=jacob+alpha*theta*du
>   
>     call make_lin_trf_  (this%pchemsys,u,sktrk,this%hashcompz,iserror)
>     if (iserror) goto 20    
>     res=res+theta*u

>     call make_lin_trf_  (this%pchemsys,du,dsktrk,this%hashcompz,iserror)
>     if (iserror) goto 20
>     jacob=jacob-theta*du
!%-----------------------------------------------------------------   
!% Dependence of immobile box at k (boxk)
!%-----------------------------------------------------------------
>     call get_chem_info_ (boxk,iserror,c=c,caq=caq,sktrk=sktrk)
>     if (iserror) goto 20
>     call make_lin_trf_  (boxk%pchemsys,u,c,this%hashcompz,iserror)
>     if (iserror) goto 20
>     res=res+u/dtime

>     call make_lin_trf_  (boxk%pchemsys,u,caq,this%hashcompz,iserror)
>     if (iserror) goto 20
>     res=res-alpha*(1.0-theta)*u
>  
>     call make_lin_trf_  (boxk%pchemsys,u,sktrk,this%hashcompz,iserror)
>     if (iserror) goto 20
>     res=res+(1.0-theta)*u
!%-----------------------------------------------------------------
!% Dependence of mobile node at k+1 (mobnodek1)
!%-----------------------------------------------------------------
>     call get_chem_info_ (mobnodek1,iserror,caq=caq)
>     if (iserror) goto 20
>     call make_lin_trf_  (mobnodek1%pchemsys,u,caq,this%hashcompz,iserror)
>     if (iserror) goto 20
>     res=res+alpha*theta*u
!%-----------------------------------------------------------------
!% Dependence of mobile node at k (mobnodek)
!%-----------------------------------------------------------------
>     call get_chem_info_ (mobnodek,iserror,caq=caq)
>     if (iserror) goto 20
>     call make_lin_trf_  (mobnodek%pchemsys,u,caq,this%hashcompz,iserror)
>     if (iserror) goto 20
>     res=res+alpha*(1.0-theta)*u
!%-----------------------------------------------------------------
!% Solve linear system, update and convergence
!%-----------------------------------------------------------------
>     delc1 = res
>     call f07adf(npri,npri,jacob,npri,indx,info)
	 if (info/=0) then 
	    iserror=.true.
		goto 20
	 end if
>     call f07aef('n',npri,1,jacob,npri,indx,delc1,npri,info)
	 if (info/=0) then 
	    iserror=.true.
		goto 20
	 end if
>     c1 = c1 + delc1
>     call set_from_cpri_ &
>    (this, &
>     boxk, &
>     c1, &
>     dtime, &
>     .false., &
>     isanomalous, &
>     isupmxitergam, &
>     iserror)
>     if (iserror) goto 20
>     if (isanomalous.or.isupmxitergam) then
>       goto 20
>     endif
>     isconvergence = (maxval(dabs(delc1/c1)).lt.tolunk .and. maxval(dabs(res)).lt.tolres)

>  end do
>  
>  if (.not.isconvergence) then
>     goto 20
>  end if

!%-----------------------------------------------------------------
!% Calculate concentration of minerals
!%-----------------------------------------------------------------
!
!%-----------------------------------------------------------------
!% Dimension of Setre
!%-----------------------------------------------------------------
>  call check_pointer_ (setre,nsp,.true.)
!%-----------------------------------------------------------------
!% Dependence of immobile box at k+1 (this)
!%-----------------------------------------------------------------
>  call get_chem_info_ (this,iserror,c=c,caq=caq,sktrk=sktrk)
>  if (iserror) goto 20
> 
>  setre = setre + c/dtime + alpha*theta*caq - theta*sktrk
!%-----------------------------------------------------------------
!% Dependence of immobile box at k (boxk)
!%-----------------------------------------------------------------
>  call get_chem_info_ (boxk,iserror,c=c,caq=caq,sktrk=sktrk)
>  if (iserror) goto 20

>  setre = setre - c/dtime + alpha*(1.0-theta)*caq - (1.0-theta)*sktrk
!%-----------------------------------------------------------------
!% Dependence of mobile node at k+1 (mobnodek1)
!%-----------------------------------------------------------------
>  call get_chem_info_ (mobnodek1,iserror,caq=caq)
>  if (iserror) goto 20
>  
>  setre = setre - alpha*theta*caq
!%-----------------------------------------------------------------
!% Dependence of mobile node at k (mobnodek)
!%-----------------------------------------------------------------
>  call get_chem_info_ (mobnodek,iserror,caq=caq)
>  if (iserror) goto 20
>  
>  setre = setre - alpha*(1.0-theta)*caq
!%-----------------------------------------------------------------
!% The unit of setre [mol/s]
!%-----------------------------------------------------------------
>  setre = phi*setre 
!%-----------------------------------------------------------------
>  call set_from_setre_(this,boxk,iscompzchanged,setre,nsp,dtime,iserror)

>  if (iscompzchanged) then
>    call update_derivatives_ (this,iserror)
	if (iserror) goto 20 
	call set_ (boxk,iserror,hashcompz=this%hashcompz)
	if (iserror) goto 20
>  end if
>  
end do
!%-----------------------------------------------------------------
!% Maximum number of mxiterminchange has been reached
!%-----------------------------------------------------------------

if(iscompzchanged)then
>  isconvergence=.false.
>  goto 20
end if

!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Calculate immobile-mobile flux, F
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% We resize if necessary due to changes in minerals
!%-----------------------------------------------------------------
call get_chem_info_ (mobnodek1,iserror,npri=ndim)
if (havefsp) then
> 
> call check_pointer_ (fsp,nsp,.true.)
> 
end if
!%-----------------------------------------------------------------
call check_pointer_ (f,ndim,.true.)
call check_pointer_ (df,ndim,ndim,.true.)
!%-----------------------------------------------------------------
!% Dependence of immobile box at k+1 (this)
!%-----------------------------------------------------------------
call get_chem_info_ (this,iserror,caq=caq)
if (iserror) goto 20
call make_lin_trf_  (this%pchemsys,u,caq,mobnodek1%hashcompz,iserror)
if (iserror) goto 20
f = f - theta*u

if (havefsp) then

> fsp = fsp - theta*caq
> 
end if


!%-----------------------------------------------------------------

!%-----------------------------------------------------------------      
!% Dependence of immobile box at k (boxk)
!%-----------------------------------------------------------------
call get_chem_info_ (boxk,iserror,caq=caq)
if (iserror) goto 20
call make_lin_trf_  (boxk%pchemsys,u,caq,mobnodek1%hashcompz,iserror)
if (iserror) goto 20
f = f - (1.0-theta)*u

if (havefsp) then

> fsp = fsp - (1.0-theta)*caq
> 
end if     
!%-----------------------------------------------------------------

!%-----------------------------------------------------------------  
!% Dependence of mobile node at k+1 (mobnodek1)
!%-----------------------------------------------------------------
call get_chem_info_ (mobnodek1,iserror,caq=caq)
if (iserror) goto 20
call make_lin_trf_  (mobnodek1%pchemsys,u,caq,mobnodek1%hashcompz,iserror)
f = f + theta*u
if (iserror) goto 20

if (havefsp) then 

> fsp = fsp + theta*caq
> 
end if
!%-----------------------------------------------------------------

!%-----------------------------------------------------------------
!% Dependence of mobile node at k (mobnodek)
!%-----------------------------------------------------------------
call get_chem_info_ (mobnodek,iserror,caq=caq)
if (iserror) goto 20
call make_lin_trf_  (mobnodek%pchemsys,u,caq,mobnodek1%hashcompz,iserror)
if (iserror) goto 20
f = f + (1.0-theta)*u
f = phi*alpha*f

if (havefsp) then 

> fsp = fsp + (1.0d0-theta)*caq
> fsp = phi*alpha*fsp
> 
end if
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Calculate derivative of immobile-mobile flux, dF
!%-----------------------------------------------------------------
!% Dependence of mobile node at k+1 (mobnodek1)
!%-----------------------------------------------------------------
call get_chem_info_ (mobnodek1,iserror,npri=nprimob,numsp=nsp,caq=caq,dcmob=dcmob)
if (iserror) goto 20
call make_lin_trf_  (mobnodek1%pchemsys,du,dcmob(1:nsp,1:nprimob),this%hashcompz,iserror)
if (iserror) goto 20
du=-alpha*theta*du
>    
call check_pointer_ (dc1imdc1m,npri,nprimob,.true.)
call f07aef('n',npri,nprimob,jacob,npri,indx,du,npri,info)
if (info/=0) then 
>  iserror=.true.
>  goto 20
end if
dc1imdc1m=-du

call make_lin_trf_  (mobnodek1%pchemsys,du,dcmob(1:nsp,1:nprimob),mobnodek1%hashcompz,iserror)
if (iserror) goto 20

df = phi*alpha*theta*du

call get_chem_info_ (this,iserror,npri=npri,numsp=nsp,caq=caq,dcmob=dcmob)
if (iserror) goto 20
call make_lin_trf_  (this%pchemsys,du,dcmob(1:nsp,1:npri),mobnodek1%hashcompz,iserror)
if (iserror) goto 20

df = df -  phi*alpha*theta*matmul(du,dc1imdc1m)


!%-----------------------------------------------------------------
20 continue  
!%-----------------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dcmob,1,1,.false.)
call check_pointer_ (du,1,1,.false.)
call check_pointer_ (dsktrk,1,1,.false.)
call check_pointer_ (sktrk ,1,.false.)
call check_pointer_ (c ,1,.false.)
call check_pointer_ (caq ,1,.false.)
call check_pointer_ (u ,1,.false.)
call check_pointer_ (c1 ,1,.false.)
call check_pointer_ (delc1 ,1,.false.)
call check_pointer_ (jacob ,1,1,.false.)
call check_pointer_ (res ,1,.false.)
call check_pointer_ (setre ,1,.false.)
call check_pointer_ (indx ,1,.false.)
if (iserror) goto 10 
!%-------------------------------------------------------
return
> 
10 continue 
print *,'****************************'
print *,'Nodal Chemistry:'
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
subroutine equilibrate_nch &
>  (this, &
>   nameph, &
>   iserror, &
>   ioutput, &
>   isconvergence, &
>   siph)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Equlibrate the nodal chemistry with one mineral phase to any saturation state. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)      :: this             !> Type nodal chemistry variable. 

character(len=*), intent(in)               :: nameph           !> Name of the phase. 

logical, intent(out)                       :: iserror          !> iserror=true, then there was an error. 

integer, intent(in), optional              :: ioutput          !> Output unit. 

logical, intent(out), optional             :: isconvergence    !> isconvergence=true, then there was convergence in speciation calculations. 

real*8, intent(in), optional               :: siph             !> Saturation of the phase 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                    :: &
> isp, &
> i, &
> j, &
> npri, &
> nsp, &
> nsurf, &
> ntxoh
real*8, pointer                            :: &
> simin(:) => null ()
real*8                                     :: &
> siphloc, &
> icputime1, &
> icputime2, &
> dtcputime
logical                                    :: &
> haveioutput, &
> haveisconvergence, &
> havesiph, &
> isconv
character(len=100)                         :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
msg=' '
iserror=.false.
!%-----------------------------------------------------------------
!% Check optional arguments
!%-----------------------------------------------------------------
haveioutput=present(ioutput)
haveisconvergence=present(isconvergence)
havesiph=present(siph)
!%-----------------------------------------------------------------
!% Check if the chemical system is associated. 
!%-----------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
>   msg='Error, not associated chemical system'
>   goto 10
end if
!%-----------------------------------------------------------------
if (havesiph) then
> siphloc=siph
else
> siphloc=0.0d0
end if
!%-----------------------------------------------------------------
if (haveioutput) then
> 
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) '      Beginning equilibrate             '
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,2) 'Equilibrate with phase:',nameph
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) 'Initial Speciation'
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     write(ioutput,*) 'Saturation indice state:',siphloc
>     write(ioutput,3) '----------------------------------------'// &
>                      '----------------------------------------'// &
>                      '-------------'
>     call write_ (this,ioutput,iserror)
>      if (iserror) then
>       msg='Error when calling write_'
>       goto 10
>      end if
end if
!%-----------------------------------------------------------------
nsp=size(this%c)
ntxoh=size(this%txoh,1)
nsurf=size(this%txoh,2)
call check_pointer_ (simin,nsp,.true.)
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
> call update_ (this%pchemsys,this%temp,iserror)
> if (iserror) goto 20
!%-----------------------------------------------------------
call cpu_time (icputime1)
!%-----------------------------------------------------------------
if (this%isderstored) then 
> call specia_ &
>     (this%pchemsys, &
	  this%temp, &
>      this%c, &
>      this%g, &
>      simin, &
>      this%ionstr, &
>      nameph, &
>      siphloc, &
>      this%txoh, &
>      this%capint, &
>      this%capext, &
>      this%spsurfarea, &
>      nsp, &
>      nsurf, &
>      ntxoh, &
>      this%omgwfree, &
>      isconv, &
>      iserror, &
>      nchemiter=this%nchemiter, & 
>      sktrk=this%sktrk, &
	  dc=this%dc, &
	  dsktrk=this%dsktrk)
else
> call specia_ &
>     (this%pchemsys, &
	  this%temp, &
>      this%c, &
>      this%g, &
>      simin, &
>      this%ionstr, &
>      nameph, &
>      siphloc, &
>      this%txoh, &
>      this%capint, &
>      this%capext, &
>      this%spsurfarea, &
>      nsp, &
>      nsurf, &
>      ntxoh, &
>      this%omgwfree, &
>      isconv, &
>      iserror, &
>      nchemiter=this%nchemiter, &
>      sktrk=this%sktrk)
end if
if (.not.isconv) goto 20
call cpu_time (icputime2)
dtcputime=icputime2-icputime1
!%---------------------------------------------------------------
!> Write speciation information
!%---------------------------------------------------------------
if (haveioutput) then
>  call write_ &
>     (this%pchemsys, &
>      this%c, &
>      this%g, &
>      this%temp, &
>      this%ionstr, &
>      nsp, &
>      ioutput, &
>      this%omgwfree, &
>      msg, &
>      iserror, &
>      simin)
>  if (iserror) goto 20
!%---------------------------------------------------------------
!> Write the cpu time
!%---------------------------------------------------------------

>   write(ioutput,3) '----------------------------------------'// &
>                    '----------------------------------------'// &
>                    '-------------'
>   write(ioutput,4) 'Total cpu time used:',dtcputime,' [s]'
>   write(ioutput,3) '----------------------------------------'// &
>                    '----------------------------------------'// &
>                    '-------------'
> 
end if
!%-----------------------------------------------------------------
20 continue 
!%-----------------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------------
call check_pointer_ (simin,1,.false.)
if (iserror) goto 10
if (haveisconvergence.and..not.isconv) then
> isconvergence=.false.
> return
else if (haveisconvergence.and.isconv) then
> isconvergence=.true.
> return 
end if
if (.not.isconv) stop
!%-----------------------------------------------------------------
return

10 continue
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: equilibrate_'
print *, msg
print *,'***********************'
iserror=.true.
return
1 format (f10.2,a9,i5,a7)
2 format (a24,a15)
3 format (a95)
4 format (a20,e10.4,a8) 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine reaction_path_temperature_nch &
>  (this, &
>   temp1, &
>   temp2, &
>   nstep, &
>   iserror, &
>   issave, &
>   ioutput)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description:Define temperature during batch-reaction-steps
!
!>   temp1=   Temperature of the first reaction step (in celsius)
!>   temp2=   Temperature of the last reaction step (in celsius)
!>   nstep=   Number of steps
!>   issave=     Bolean. If is .true. then set the last temperature
!>            in the nodal chemistry object
!>   ioutput= Unit for to write output file (optional)
!
!
!
!>  Warning: Is necessary to define the chemical system as
!>           "non isoterm" for to use this public service.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)    :: this      !> Type nodal chemistry variable. 

real*8, intent(in)                       :: temp1     !> Temperature 1

real*8, intent(in)                       :: temp2     !> Temperature 2

integer, intent(in)                      :: nstep     !> Number of steps

logical, intent(in)                      :: issave    !> issave=true, then 

logical, intent(out)                     :: iserror   !> iserror=true, then there was an error

integer, intent(in), optional            :: ioutput   !> output unit 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                  :: &
> i
logical                                  :: &
> haveioutput
real*8                                   :: &
> tempi, &
> delta
type(t_nodalchemistry), pointer          :: &
> pnch => null ()
character(len=100)                       :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
> 

> 

!%-----------------------------------------------------
msg=' '
iserror=.false.
!%-----------------------------------------------------
!% Check optional arguments 
!%-----------------------------------------------------
haveioutput=present(ioutput)
!%-----------------------------------------------------
!% Check if the chemical system is associated. 
!%-----------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------
if (nstep==0) return
!%-----------------------------------------------------
!% Allocate and create local nodal chemistry and copy
!%-----------------------------------------------------
allocate (pnch)
call create_ (pnch)
pnch=this
!%-----------------------------------------------------
if (haveioutput) then
>     write(ioutput,*) '----------------------------------------'
>     write(ioutput,*) 'Beginning of batch-reaction calculations'
>     write(ioutput,*) '----------------------------------------'
>     write(ioutput,*) 'Initial Speciation'
>     write(ioutput,*) '------------------'
>     call write_ (pnch,ioutput,iserror)
>     if (iserror) then
>       msg='Error when calling write_'
>       goto 20
>     end if
end if
!%--------------------------------------------------
delta=(temp2-temp1)/real(nstep)
tempi=temp1
!%--------------------------------------------------
do i=1,nstep
> 
>   call set_ (pnch,iserror,temp=tempi)

>   if (haveioutput) then
>      write(ioutput,*) 'Reaction temperature step', i
>      write(ioutput,*) '-------------------------'
>      call write_ (pnch,ioutput,iserror)
>      if (iserror) then
>        msg='Error when calling write_'
>        goto 20
>      end if
>   end if

> 
>   tempi=tempi+delta
> 
end do
!%--------------------------------------------------
!% If issave=true, then copy the nodal chemistry 
!% object 
!%--------------------------------------------------
if (issave) this=pnch 
!%--------------------------------------------------
20 continue 
!%--------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------
call destroy_ (pnch)
deallocate (pnch)
pnch => null ()
if (iserror) goto 10
!%------------------------------------------------------------
return
> 
> 
10 continue 
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: reaction_path_'
print *, msg
print *,'***********************'
iserror=.true.
return
> 
1 format (f10.2,a9,i5,a7)
2 format (a10,f10.2)
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine add_txoh_nch &
>  (this, &
>   namesurf, &
>   txoh, &
>   capint, &
>   capext, &
>   spsurfarea, &
>   nsite, &
>   isequilibrate, &
>   iserror, &
>   isconvergence, &
>   dtime, &
>   namemin)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Add txoh in the nodal chemistry 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)   :: this           !> Type nodal chemistry variable. 

integer, intent(in)                     :: nsite          !> Number of sites 

character(len=*), intent(in)            :: namesurf       !> Number of txoh 

real*8, intent(in), dimension(nsite)    :: txoh           !> Total of sites [mol/kgw]

real*8, intent(in), dimension(nsite)    :: capint           !> Capacitance 1 

real*8, intent(in), dimension(nsite)    :: capext           !> Capacitance 2

real*8, intent(in), dimension(nsite)    :: spsurfarea     !> Specific surface area [m2 kgw -1]

logical, intent(out)                    :: iserror        !> iserror=true, then there was an error

logical, intent(in)                     :: isequilibrate  !> isequilibrate=true 

logical, intent(out), optional          :: isconvergence  !> If .true. there was convergence in chemistry 

real*8, intent(in), optional            :: dtime          !> Time increment 

character(len=*), intent(in), optional  :: namemin        !> Name of the mineral species that is related with CEC 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond: the total of sites (txoh) is expressed in [mol/m3 rock]
!
!>   $Post-cond:
!
!>   $License: CSIC-UPC
!
!-------------------------------------------------------------------------
integer                                 :: &
> i, &
> isps, &
> isite1, &
> isite2, &
> isurf, &
> nsiteloc, &
> nsp, &
> ntxoh, &
> nsurf, &
> npri, &
> nminsp
real*8                                 :: &
> dtimeloc, &
> factoromgw, &
> faccap 
integer, pointer                       :: &
> idpri(:) => null ()
real*8, pointer                        :: &
> u(:) => null (), &
> cold(:) => null ()
character(len=100), pointer            :: &
> nameminsp(:) => null () 
character(len=100)                     :: &
> msg
logical                                :: &
> isanomalous, &
> isupmxitergam, &
> isconv, &
> haveisconvergence, &
> havedtime, &
> havenamemin 
real*8, parameter                      :: &
> r0=0.0d0, &
> r1=1.0d0  
!-------------------------------------------------------------------------
!
!>   $code
!
> 
!%-----------------------------------------------------------------
iserror=.false.
msg=' '
!%-----------------------------------------------------------------
!% Check if the chemical system is associated 
!%-----------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------------------
!% Check optional arguments
!%-----------------------------------------------------------------
haveisconvergence=present(isconvergence)
havedtime=present(dtime)
havenamemin=present(namemin)
!%-----------------------------------------------------------------
!% Initialice variables
!%-----------------------------------------------------------------
isanomalous=.false.
isupmxitergam=.false.
isconv=.true.
!%-----------------------------------------------------------------
!% Check the name of the surface 
!%-----------------------------------------------------------------
call get_surf_index_ (this%pchemsys,isurf,nsiteloc,namesurf)
if (isurf==0) then
> msg='Error, not defined in the chemical system the surface:'
> call add_ (msg,namesurf)
> goto 10
end if
!%-----------------------------------------------------------------
!% Check the number of sites 
!%-----------------------------------------------------------------
if (nsite/=nsiteloc) then
> msg='Error in number of sites in surface:'
> call add_ (msg,namesurf)
> goto 10
end if
!%-----------------------------------------------------------------
!% If dtime is present as argument starge dt in dtloc
!%-----------------------------------------------------------------
if (havedtime) then
> dtimeloc=dtime
else
> dtimeloc=r0
end if
!%-----------------------------------------------------------------
!% Storage txoh in mol
!% We suppose that txoh is in mol m-3 of nodal chemistry 
!%-----------------------------------------------------------------
this%txoh(1:nsite,isurf)=txoh*this%volnch
this%txoh0(1:nsite,isurf)=this%txoh(1:nsite,isurf)
!%-----------------------------------------------------------------
!% 
!%-----------------------------------------------------------------
this%capint(1:nsite,isurf)=capint
this%capext(1:nsite,isurf)=capext
this%spsurfarea(1:nsite,isurf)=spsurfarea 
!%-----------------------------------------------------------------
!% Find the index of mineral species that depends of the surface 
!%-----------------------------------------------------------------
if (havenamemin) then
> call get_chem_info_ (this%pchemsys,iserror,nameminsp=nameminsp,nminsp=nminsp)
> do i=1,nminsp
>     if (namemin==nameminsp(i)) then
>       call get_sp_index_ (this%pchemsys,namemin,isps)
>       this%idtxohmin(isurf)=isps
	   exit 
	 end if 
> end do 
> call check_pointer_ (nameminsp,1,.false.) 
end if 
!%-----------------------------------------------------------------
!% Initialice 
!%-----------------------------------------------------------------
factoromgw=r1
!%-----------------------------------------------------------------
!% 
!%-----------------------------------------------------------------
this%alpha=this%alpha/this%omgwfree
!%-----------------------------------------------------------------
if (isequilibrate) then !> If not isequilibrate the solution specia from
> 
> call get_chem_info_ (this%pchemsys,iserror,numbase=npri, &
>                     idbase=idpri,numsp=nsp,numsites=ntxoh, &
>                     numsurf=nsurf)
> 
> if (iserror) then 
>   msg='Error when call get_chem_info_'
>   goto 20
> end if 
> 
> if (this%isderstored) then
> 
>  call make_lin_trf_(this%pchemsys,u,this%c,0,iserror)
>  
>  if (iserror) goto 20
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
> call update_ (this%pchemsys,this%temp,iserror)
> if (iserror) goto 20
!%-----------------------------------------------------------
>  call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	.true., &
	this%ispgasconst, &
>    this%ionstr, &
>    this%alpha, &
>    u, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    nsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    dtimeloc, &
>    this%hashcompz, &
>    isconv, &
>    factoromgw, &
	this%volgas, &
	this%pgas, &
	msg, &
>    iserror, &
>    cguess=this%c(idpri), &
>    sktrk=this%sktrk, &
>    dc=this%dc, &
>    dsktrk=this%dsktrk, &
>    nchemiter=this%nchemiter)
> 
> else
> 
>  call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	.true., &
	this%ispgasconst, &
>    this%ionstr, &
>    this%alpha, &
>    u, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    nsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    dtimeloc, &
>    this%hashcompz, &
>    isconv, &
>    factoromgw, &
	this%volgas, &
	this%pgas, &
	msg, &
>    iserror, &
>    cguess=this%c(idpri), &
>    sktrk=this%sktrk, &
>    nchemiter=this%nchemiter)
> 
> end if
> 
> this%alpha=this%alpha*this%omgwfree
> this%omgwfree = this%omgwfree * factoromgw 
> !%-------------------------------------------------------
> !% Update the reactive surface of kinetic minerals 
> !% (only if there was convergence in the speciation)
> !% Update the TXOH according when it is associated
> !% to mineral phase 
> !%-------------------------------------------------------
>   if (isconv) then  
>     call update_min_area_ (this,iserror)
>     if (iserror) goto 20
	 call update_txoh_ (this,iserror)
>     if (iserror) goto 20 
>   end if 
else    !> If not isequilibrate the solution specia straigthforward
> 
> call get_chem_info_ (this%pchemsys,iserror,numbase=npri, &
>                     idbase=idpri,numsp=nsp,numsites=ntxoh, &
>                     numsurf=nsurf,hashcompz=this%hashcompz)
> 
> call check_pointer_ (cold,nsp,.true.)
> cold=this%c
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
> call update_ (this%pchemsys,this%temp,iserror)
> if (iserror) goto 20
!%-----------------------------------------------------------
!%  
!%-----------------------------------------------------------
call get_chem_info_ (this,iserror,faccap=faccap)
if (iserror) goto 20
!%-----------------------------------------------------------
> if (this%isderstored) then
>  
>  call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	cold, &
>    this%alpha, &
>    nsp, &
>    this%c(idpri), &
>    npri, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    ntxoh, &
>    nsurf, &
>    dtimeloc, &
	this%volgas, &
>    this%ionstr, &
>    this%hashcompz, &
	faccap, &
>    isanomalous, &
>    isupmxitergam, &
>    .false., &
>    iserror, &
>    sktrk=this%sktrk, &
>    dc=this%dc, &
>    dsktrk=this%dsktrk)
>    
> else
> 
>  call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
	cold, &
>    this%alpha, &
>    nsp, &
>    this%c(idpri), &
>    npri, &
>    this%txoh/this%omgwfree, &
>    this%capint, &
>    this%capext, &
>    this%spsurfarea, &
>    ntxoh, &
>    nsurf, &
>    dtimeloc, &
	this%volgas, &
>    this%ionstr, &
>    this%hashcompz, &
	faccap, &
>    isanomalous, &
>    isupmxitergam, &
>    .false., &
>    iserror, &
>    sktrk=this%sktrk)
> 
> end if
> 
> this%alpha=this%alpha*this%omgwfree
> call check_pointer_ (cold,nsp,.false.)
end if
!%--------------------------------------------------
!% Deallocate local pointers
!%--------------------------------------------------
20 continue 
call check_pointer_ (idpri,1,.false.)
call check_pointer_ (u,1,.false.)
if (iserror) goto 10
!%--------------------------------------------------
if (haveisconvergence.and..not.isconv) then
>  isconvergence=.false.
>  return
end if
!%--------------------------------------------------
if (haveisconvergence.and.isanomalous) then
>  isconvergence=.false.
>  return
end if
!%--------------------------------------------------
if (haveisconvergence.and.isupmxitergam) then
>  isconvergence=.false.
>  return
end if
!%--------------------------------------------------
if (haveisconvergence) then
>  isconvergence=.true.
>  return
end if
!%--------------------------------------------------
if (isanomalous.or.isupmxitergam.or..not.isconv) stop
!%-----------------------------------------------------------------
return
> 
10 continue 
print *,'************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: add_txoh_'
print *, msg
print *,'************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine add_species_nch &
>  (this, &
>   namesp, &
>   mass, &
>   unit, &
>   isspecia, &
>   iserror, &
>   alpha, &
>   isconvergence, &
>   dtime)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Add species to nodal chemistry object and compute the 
!> chemical speciation (optional, equilibrate=true). 
!> Also other parameters related with the species can be added (e.g.
!> reactive surface, etc)
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)          :: this           !> Type nodal chemistry variable. 

character(len=*), intent(in)                   :: namesp         !> Name of species. 

character(len=*), intent(in)                   :: unit           !> Unit of mass. 

real*8, intent(in)                             :: mass           !> Mass of the species. 

logical, intent(in)                            :: isspecia       !> isspecia=true, then the nodal chemistry is speciated. 

logical, intent(out)                           :: iserror        !> iserror=true, then there was an error. 

logical, intent(out), optional                 :: isconvergence  !> isconvergence=true, there was convergence in speciation calculations. 

real*8, intent(in), optional                   :: alpha          !> Area [m2/m3 rock]

real*8, intent(in), optional                   :: dtime          !> Time increment [s]
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                        :: &
> isps, &
> npri, &
> ntxoh, &
> nsp, &
> nsurf, &
> hashcompz
logical                                        :: &
> isconv, &
> havealpha, &
> havedtime, &
> haveisconvergence, &
> isueq, &
> isequilibrate, &
> isunitm3, &
> isanomalous, &
> isupmxitergam 
real*8                                         :: &
> dtimeloc, &
> massloc, &
> factoromgw
real*8, pointer                                :: &
> t(:) => null (), &
> told(:) => null (), &
> c(:) => null (), &
> factor(:) => null ()
integer, pointer                               :: &
> idpri(:) => null ()
character(len=100)                             :: &
> msg
real*8, parameter                              :: &
> r0=0.0d0, &
> r1=1.0d0 
!-------------------------------------------------------------------------
!
!>   $code
!
!%-----------------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------------
!%Check optional arguments
!%-----------------------------------------------------------------
havealpha=present(alpha)
havedtime=present(dtime)
haveisconvergence=present(isconvergence)
!%-----------------------------------------------------------------
!% Check if the chemical system is associated 
!%-----------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------------------
!% Check if the species was defined in the chemical system 
!%-----------------------------------------------------------------
call get_sp_index_ (this%pchemsys,namesp,isps)
if (isps==0) then
> msg='Error, species not defined in the chemical system:'
> call add_ (msg,namesp)
> goto 10
end if
!%-----------------------------------------------------------------
!% Initialice variables 
!%-----------------------------------------------------------------
factoromgw=r1
isconv=.true. 
!%-----------------------------------------------------------------
!% Check if the mass of species is expressed for m3 rock 
!%-----------------------------------------------------------------
call find(unit,'/m3rock',isunitm3)
!%-----------------------------------------------------------------
!> Her determine if time increment (dtime) = 0.
!> If dtime > 0 then considers operations with equilibrium matrix components.
!> If dtime = 0 then considers operations with global matrix of components. 
!%-----------------------------------------------------------------
dtimeloc=r0
isueq=.false. 
if (havedtime.and.dtime>r0) then
> dtimeloc=dtime
> isueq=.true. 
end if
!%-----------------------------------------------------------------
!% Change the chemical units
!%-----------------------------------------------------------------
massloc=mass
call change_chem_unit_(this%pchemsys,massloc,namesp,unit,'m',iserror)
!%-----------------------------------------------------------------
!% Storage the mass of species in massloc and multiply of the 
!% volume associated to nodal chemistry
!%-----------------------------------------------------------------
if (isunitm3) then 
> massloc=massloc*this%volnch
> this%alpha0(1,isps)=massloc
> massloc=massloc/this%omgwfree 
else
> this%alpha0(1,isps)=massloc*this%omgwfree 
end if 
!%------------------------------------------------------------------
!% Add aditional information as area of mineral (expressed in 
!% [m2/m3 rock]
!%------------------------------------------------------------------
if (havealpha) then
> this%alpha(isps)=alpha*this%volnch
> this%alpha0(2,isps)=this%alpha(isps)
end if
!%-----------------------------------------------------------------
!% Get information to chemical system 
!%-----------------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,numbase=npri,idbase=idpri, &
>                      numsp=nsp,numsites=ntxoh,numsurf=nsurf) 
if (iserror) goto 20 
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
call update_ (this%pchemsys,this%temp,iserror)
if (iserror) goto 20                      
!%-----------------------------------------------------------------
!% If the nodal chemistry is equilibrated with mineral phases, then
!% compute the chemical speciation 
!%-----------------------------------------------------------------
if (isspecia) then
>     isequilibrate=.true. 
>     call check_pointer_ (factor,npri,.true.)
>     call check_pointer_ (c,nsp,.true.)
!%-----------------------------------------------------------------
!% Compute told=U*cold 
!%-----------------------------------------------------------------
>     call make_lin_trf_ (this%pchemsys,told,this%c,0,iserror,isueq=isueq)
>     if (iserror) goto 20 
>     c=this%c 
>     c(isps)=c(isps)+massloc
!%-----------------------------------------------------------------
!% Compute t=U*c 
!%----------------------------------------------------------------- 
>     call make_lin_trf_ (this%pchemsys,t,c,0,iserror,isueq=isueq)
>     call check_pointer_ (c,1,.false.)
>     if (iserror) goto 20 
>     factor=t/told 
>     this%c(idpri)=this%c(idpri)*factor 
>     call check_pointer_ (factor,1,.false.) 
!%-----------------------------------------------------------------
!% The water mass balance is considered
!%-----------------------------------------------------------------
>     call set_iswcompbal_ (this%pchemsys,.true.,iserror)
>     if (iserror) goto 20   
!%-----------------------------------------------------------------
!% Specia according total concentrations 
!%-----------------------------------------------------------------
>     this%alpha = this%alpha / this%omgwfree
!%-----------------------------------------------------------     
	 if (this%isderstored) then
>      
	   call specia_ &
>        (this%pchemsys, &
		 this%temp, &
>         this%c, &
>         this%g, &
>         isequilibrate, &
	     this%ispgasconst, &
>         this%ionstr, &
>         this%alpha, &
>         t, &
>         this%txoh/this%omgwfree, &
>         this%capint, &
>         this%capext, &
>         this%spsurfarea, &
>         nsp, &
>         npri, &
>         ntxoh, &
>         nsurf, &
>         dtimeloc, &
>         this%hashcompz, &
>         isconv, &
>         factoromgw, &
	     this%volgas, &
>         this%pgas, &
>         msg, &
>         iserror, &
>         cguess=this%c(idpri), &
>         sktrk=this%sktrk, &
>         dc=this%dc, &
>         dsktrk=this%dsktrk, &
>         nchemiter=this%nchemiter)

>      else
>    
	   call specia_ &
>        (this%pchemsys, &
		 this%temp, &
>         this%c, &
>         this%g, &
	     isequilibrate, &
	     this%ispgasconst, &
>         this%ionstr, &
>         this%alpha, &
>         t, &
>         this%txoh/this%omgwfree, &
>         this%capint, &
>         this%capext, &
>         this%spsurfarea, &
>         nsp, &
>         npri, &
>         ntxoh, &
>         nsurf, &
>         dtimeloc, &
>         this%hashcompz, &
>         isconv, &
>         factoromgw, &
	     this%volgas, &
	     this%pgas, &
	     msg, &
>         iserror, &
>         cguess=this%c(idpri), &
>         sktrk=this%sktrk, &
>         nchemiter=this%nchemiter)

>      end if
>       
	  this%alpha = this%alpha * this%omgwfree
>      this%omgwfree = this%omgwfree * factoromgw
>      
else !> Add the mass loc in the vector c and compute derivatives 
>      
	  
	  this%c(isps)=this%c(isps)+massloc  
>      nsp=size(this%c) 
>      hashcompz=this%hashcompz
>      call get_hashcompz_ (this%pchemsys,this%hashcompz,this%c, &
>                          this%g,nsp,iserror) 
>      if (hashcompz/=this%hashcompz) then
>         call update_derivatives_ (this,iserror)
>      end if 

end if 
!%-----------------------------------------------------------------
20 continue
!%-----------------------------------------------------------------
!%Deallocate local pointers 
!%----------------------------------------------------------------- 
call check_pointer_ (idpri,1,.false.)
call check_pointer_ (t,1,.false.)
call check_pointer_ (told,1,.false.)
if (iserror) goto 10
if(haveisconvergence.and.isconv) then
> isconvergence=.true.
> return
else if (haveisconvergence.and..not.isconv) then
> isconvergence=.false.
> return
end if
if (.not.isconv) stop
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
> 
return
> 
10 continue
print *,'**********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: add_species_'
print *, msg
print *,'**********************'
iserror=.true.
return
> 
> 
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
function scalarbynodalchem_nch &
>  (scalar, &
>   this)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Multiply the nodal chemistry attributes by a scalar.  
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)    :: this                  !> Type nodal chemistry variable. 

type(t_nodalchemistry)                :: scalarbynodalchem_nch !> Type nodal chemistry variable. 

real*8, intent(in)                    :: scalar                !> Scalar. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
logical                  :: &
iserror 
character(len=100)       :: &
msg
integer                  :: &
isps1, &
isps2
!-------------------------------------------------------------------------
!
!>   $code
!
!
!%------------------------------------------------------------------------
msg='' 
!%------------------------------------------------------------------------
!% Create nodal chemistry object 
!%------------------------------------------------------------------------
call create_ (scalarbynodalchem_nch)
!%------------------------------------------------------------------------
!% Copy the nodal chemistry object 
!%------------------------------------------------------------------------
scalarbynodalchem_nch = this
!%------------------------------------------------------------------------
!% Multiply the mass of water 
!%------------------------------------------------------------------------
scalarbynodalchem_nch%omgwfree = scalar*scalarbynodalchem_nch%omgwfree
!%------------------------------------------------------------------------
!% Multiply the volume
!%------------------------------------------------------------------------
scalarbynodalchem_nch%volnch = scalar*scalarbynodalchem_nch%volnch
!%------------------------------------------------------------------------
!% Multiply the volume of gas
!%------------------------------------------------------------------------
scalarbynodalchem_nch%volgas = scalar*scalarbynodalchem_nch%volgas
!%------------------------------------------------------------------------
!% Multiply the TXOH
!%------------------------------------------------------------------------
if (associated(scalarbynodalchem_nch%txoh)) then
>   scalarbynodalchem_nch%txoh = scalar*scalarbynodalchem_nch%txoh
end if 
!%------------------------------------------------------------------------
return
> 
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
subroutine write_nch_info_nch &
>  (this, &
>   ioutput, &
>   iserror, &
>   iswchsys)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Write the information about nodal chemistry object.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)   :: this     !> Type nodal chemistry variable. 

integer, intent(in)                     :: ioutput  !> Output unit 

logical, intent(out)                    :: iserror  !> iserror=true, then there was an error 

logical, intent(in), optional           :: iswchsys !> If .true to write the chemical system information 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                 :: &
> nsp
character(len=100)                      :: &
> msg
logical                                 :: &
> haveiswchsys 
!-------------------------------------------------------------------------
!
!>   $code
!
> 
!%-----------------------------------------------------------
msg=''
iserror=.false.
!%-----------------------------------------------------------
!% Check optional arguments
!%-----------------------------------------------------------
haveiswchsys=present(iswchsys)
!%-----------------------------------------------------------
!% Check if the chemical system is associated
!%-----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%--------------------------------------------------------------
!% Update the temperature in the chemical system 
!%--------------------------------------------------------------
call update_ (this%pchemsys,this%temp,iserror)
if (iserror) goto 10 
!%--------------------------------------------------------------
!> Write nodal chemistry information
!%--------------------------------------------------------------
write(ioutput,*) '-------------------------------------------------'
write(ioutput,*) '-        Nodal Chemistry information            -'
write(ioutput,*) '-------------------------------------------------'
write(ioutput,1) 'Name= ',this%name
write(ioutput,*) '-------------------------------------------------'
write(ioutput,3) 'volume:',this%volnch,'[m3]'
write(ioutput,*) '-------------------------------------------------'
!%-------------------------------------------------------------
!% Write the speciation information 
!%---------------------------------------------------------------
nsp=size(this%c)
call write_ &
>   (this%pchemsys, &
>    this%c, &
>    this%g, &
>    this%temp, &
>    this%ionstr, &
>    nsp, &
>    ioutput, &
>    this%omgwfree, &
>    msg, &
>    iserror, &
>    sktrk=this%sktrk)
!%-------------------------------------------------------------
if (iserror) goto 10
write(ioutput,4) 'number of iter. in chem. speciation:',this%nchemiter 
write(ioutput,2) '----------------------------------------------'// &
>                 '-----------------------------------------------'
!%-------------------------------------------------------------
!%--------------------------------------------------------------
!> Write the chemical system associated to nodal chemistry (optional)
!%--------------------------------------------------------------
if (haveiswchsys.and.iswchsys) then
> write(ioutput,2) '----------------------------------------------'// &
>                  '-----------------------------------------------'
> call write_ (this%pchemsys,ioutput,iserror)
> if (iserror) goto 10   
> write(ioutput,2) '----------------------------------------------'// &
>                  '-----------------------------------------------'
end if
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
return
1 format(a15,a100)
2 format(a95)
3 format(a8,e10.4,a5)
4 format(a36,i5)
> 
10 continue 
print *,'*******************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: write_'
print *, msg
print *,'*******************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_sps_nch &
>  (this, &
>   ioutput, &
>   typeoutput, &
>   namesp, &
>   nsp, &
>   iserror, &
>   iswritename, &
>   realitem, &
>   integeritem)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Write in txt file concetration of species according list (namesp).
!> For the case of H+ species, writes the pH. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)              :: this        !> Type nodal chemistry variable. 

integer, intent(in)                             :: nsp         !> Number of species 

character(len=*), intent(in), dimension(nsp)    :: namesp      !> List of name od species

integer, intent(in)                             :: ioutput     !> Output unit 

character(len=*), intent(in)                    :: typeoutput  !> Type of output

logical, intent(out)                            :: iserror     !> iserror=true, there was an error. 

logical, intent(in), optional                   :: iswritename !> If .true. then write the name species

real*8, intent(in), optional                    :: realitem    !> Real item (optional)   

integer, intent(in), optional                   :: integeritem !> Integer item (optional)
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License: 
!
!-------------------------------------------------------------------------
> 
integer                            :: &
> isps, &
> i, &
> ispsh, &
> ispsw, &
> item 
real*8                             :: &
> value
real*8, pointer                    :: &
> cloc(:) => null ()
character(len=100)                 :: &
> msg, &
> name 
character(len=100), pointer        :: & 
> namesploc(:) => null ()
logical                            :: &
> haveiswritename, &
> haverealitem, &
> haveintegeritem, &
> iswrite
!-------------------------------------------------------------------------
!
!>   $code
!
> 

!%------------------------------------------------------------
if (nsp==0) return
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
haveiswritename=present(iswritename)
haverealitem=present(realitem)
haveintegeritem=present(integeritem)
!%------------------------------------------------------------
if (haveiswritename) then
> iswrite=iswritename
else
> if (item<=1) then
>  iswrite=.true.
> else
>  iswrite=.false. 
> end if
end if
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check if the chemical system is associated
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (cloc,nsp,.true.)
call check_pointer_ (namesploc,nsp,.true.)
!%-----------------------------------------------------------
call get_sp_index_ (this%pchemsys,'h+',ispsh)
call get_sp_index_ (this%pchemsys,'h2o',ispsw)
!%-----------------------------------------------------------
do i=1,nsp
> 
>   call get_sp_index_ (this%pchemsys,namesp(i),isps)
> 
>   if (isps>0) then

>     select case(typeoutput)
>   
>     case ('activity')
>   
>           if (isps==ispsh) then
>              value=-dlog10(this%c(isps)*this%g(isps))
>   	          name='pH'
>           else if (isps==ispsw) then
>              value=this%g(isps)
		      name=namesp(i)
>           else
	          value=this%c(isps)*this%g(isps)
	          name=namesp(i)
	       end if
	  
>     case ('concentration')
>        
		   value=this%c(isps)
>           name=namesp(i)

>     case ('activity coefficient')
>        
		   value=this%g(isps)
>           name=namesp(i) 
>   
>     case ('mol/m3')
>        
		   value=this%c(isps)*this%omgwfree/this%volnch
>           name=namesp(i)

>     case ('kgr')
>        
		   value=this%c(i)*this%omgwfree
		   name=namesp(i)
		   call change_chem_unit_(this%pchemsys,value,name,'mol','kgr',iserror)
		   if (iserror) then 
		    iserror=.false. 
		    value=0.0d0 
		   end if 

>     case ('mol')
>        
		   value=this%c(isps)*this%omgwfree
>           name=namesp(i)		 

>     case ('kinetic changes')
>        
		   value=this%sktrk(isps)*this%omgwfree
>           name=namesp(i) 

>     case ('equilibrium changes')
>        
		   value=this%setre(isps)*this%omgwfree
>           name=namesp(i) 
>  
>     case default

>           msg='Error, not recognized the output type:'
		   call add_ (msg,typeoutput)
		   iserror=.true. 
		   goto 20

>     end select 
> 
>     cloc(i)=value
>     namesploc(i)=name

>  end if 
> 
end do
!%-----------------------------------------------------------
!% If real item is present then write with the real label 
!%-----------------------------------------------------------
if (haverealitem) then  
>  call write_ (ioutput,'',namesploc,cloc,nsp,realitem,iswrite,.false.)
>  goto 20 
end if
!%-----------------------------------------------------------
!% 
!%-----------------------------------------------------------
if (haveintegeritem) then
>  item=integeritem 
else
>  item=0 
end if
call write_ (ioutput,'',namesploc,cloc,nsp,item,iswrite,.false.)
!%-----------------------------------------------------------
20 continue
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (cloc,1,.false.)
call check_pointer_ (namesploc,1,.false.)
!%-----------------------------------------------------------
if (iserror) goto 10 
!%-----------------------------------------------------------
return
> 
> 
> 
10 continue 
print *,'*********************'
print *,'Nodal Chemistry:'
print *,'Namë:',this%name
print *,'Service:  write_'
print *, msg
print *,'********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_chemical_parameters_nch &
>  (this, &
>   ioutput, &
>   nameparam, &
>   nparam, &
>   iserror, &
>   iswritename, &
>   realitem, &
>   integeritem)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Write general chemical parameters specified for the external user
!> (e.g. pH, ionic strength, water activity)
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)              :: this        !> Type nodal chemistry variable. 

integer, intent(in)                             :: ioutput     !> Output unit 

integer, intent(in)                             :: nparam      !> Number of parameters  

character(len=*), dimension(nparam), intent(in) :: nameparam   !> Name of parameters 

logical, intent(out)                            :: iserror     !> iserror=true, there was an error. 

logical, intent(in), optional                   :: iswritename !> If .true. then write the name species

real*8, intent(in), optional                    :: realitem    !> Real item (optional) 

integer, intent(in), optional                   :: integeritem !> Integer item (optional)
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License: 
!
!-------------------------------------------------------------------------
> 
integer                            :: &
> isps, &
> i, &
> item 
real*8, pointer                    :: &
> value(:) => null ()
character(len=100)                 :: &
> msg
logical                            :: &
> haveiswritename, &
> haverealitem, &
> haveintegeritem, &
> iswrite
!-------------------------------------------------------------------------
!
!>   $code
!
> 

!%------------------------------------------------------------
if (nparam==0) return
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
haveiswritename=present(iswritename)
haverealitem=present(realitem)
haveintegeritem=present(integeritem)
!%------------------------------------------------------------
if (haveiswritename) then
> iswrite=iswritename
else
> if (item<=1) then
>  iswrite=.true.
> else
>  iswrite=.false. 
> end if
end if
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check if the chemical system is associated
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (value,nparam,.true.)
!%-----------------------------------------------------------
do i=1,nparam

>   select case(nameparam(i))
>   
>   case ('ionic strength','IONIC STRENGTH')

	    value(i)=this%ionstr

>   case ('pH','PH')

>        call get_sp_index_ (this%pchemsys,'h+',isps)
	    if (isps>0) then
	      value(i)=this%c(isps)*this%g(isps)
		  value(i)=-dlog10(value(i))
		end if

>   case ('temperature','TEMPERATURE')

>        value(i)=this%temp

>   case ('water activity','WATER ACTIVITY','aw','AW')

>        call get_sp_index_ (this%pchemsys,'h2o',isps)
	    if (isps>0) then
	      value(i)=this%g(isps)
		end if

>   case ('mass of water','mass water','MASS OF WATER','MASS WATER')

>        value(i)=this%omgwfree

>   case ('xsalt','XSALT')

>        call get_chem_info_ (this,iserror,xsalt=value(i))
		if (iserror) then
		  iserror=.false. 
		  value(i)=0.0d0
		end if  

>   case ('volume','VOLUME')

>        value(i)=this%volnch

>   case ('liquid density','LIQUID DENSITY')

>        value(i)=this%liqdens

>   case ('liquid pressure','LIQUID PRESSURE')

>        value(i)=this%pliq

>   case ('gas pressure','GAS PRESSURE')

>        value(i)=this%pgas
		
>   case ('capilar factor','CAPILAR FACTOR')

>        call get_chem_info_ (this,iserror,faccap=value(i))
		if (iserror) goto 20  
		
>   case default

>        msg='Error, not recognized the output type:'
		call add_ (msg,nameparam(i))
		iserror=.true. 
		goto 20

>   end select 
> 
end do
!%-----------------------------------------------------------
!% Write with real label 
!%-----------------------------------------------------------
if (haverealitem) then  
>  call write_ (ioutput,'',nameparam,value,nparam,realitem,iswrite,.false.)
>  goto 20 
end if
!%-----------------------------------------------------------
!% Wrtie with integer label 
!%-----------------------------------------------------------
if (haveintegeritem) then
>  item=integeritem 
else
>  item=0 
end if
call write_ (ioutput,'',nameparam,value,nparam,item,iswrite,.false.)
!%-----------------------------------------------------------
20 continue
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (value,1,.false.)
!%-----------------------------------------------------------
if (iserror) goto 10 
!%-----------------------------------------------------------
return
> 
10 continue 
print *,'*********************'
print *,'Nodal Chemistry:'
print *,'Namë:',this%name
print *,'Service:  write_'
print *, msg
print *,'********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_sat_nch &
>  (this, &
>   ioutput, &
>   item, &
>   iswritehead, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Write in ascII the saturation of the solution with respect to mineral species. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout) :: this         !> Type nodal chemistry variable. 

integer, intent(in)                   :: ioutput      !> Output unit.

integer, intent(in)                   :: item         !> Item to write when start the line 

logical, intent(in)                   :: iswritehead  !> iswritename=true, then the name of the mineral species are written. 

logical, intent(out)                  :: iserror      !> iserror=true, there was an error. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                            :: &
> nmin, &
> nsp
real*8, pointer                    :: &
> simin(:) => null ()
integer, pointer                   :: &
> idmin(:) => null ()
character(len=100), pointer        :: &
> namemin(:) => null ()
character(len=100)                 :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check if the chemical system is associated. 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
nsp=size(this%c)
!%-----------------------------------------------------------
!% Update the temperature in the chemical system 
!%-----------------------------------------------------------
call update_ (this%pchemsys,this%temp,iserror)
if (iserror) goto 20
!%-----------------------------------------------------------
!% Compute saturation indices 
!%-----------------------------------------------------------
call compute_min_sat_index_ &
>   (this%pchemsys, &
>    simin, &
>    idmin, &
>    namemin, &
>    nmin, &
>    this%c, &
>    this%g, &
>    nsp, &
>    iserror)
if (iserror) goto 20
!%------------------------------------------------------------
if (nmin==0) goto 20
!%------------------------------------------------------------
!% Compute Si = log( IAP/K)
!%------------------------------------------------------------
simin=dlog10(simin)
!%------------------------------------------------------------
!% Write saturation indices
!%------------------------------------------------------------
call write_ (ioutput,'',namemin,simin,nmin,item,iswritehead,.false.)
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (idmin,1,.false.)
call check_pointer_ (simin,1,.false.)
call check_pointer_ (namemin,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
> 
> 
10 continue 
print *,'************************'
print *,'Nodal Chemistry:'
print *,'Namë:',this%name
print *,'Service:  write_sat_'
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
subroutine write_min_area1_nch &
>  (this, &
>   ioutput, &
>   head, &
>   item, &
>   iswritehead, &
>   ism3nch, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Write in ascII the reactive surface of all mineral species 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in) :: this         !> Type nodal chemistry variable. 

integer, intent(in)                :: ioutput      !> Output unit.

integer, intent(in)                :: item         !> Item to write when start the line 

logical, intent(in)                :: iswritehead  !> iswritename=true, then the name of the mineral species are written. 

character(len=*), intent(in)       :: head         !> Head 

logical, intent(in)                :: ism3nch      !> If .true. the reactive surface are in m3 of nodal chemistry

logical, intent(out)               :: iserror !> iserror=true, there was an error. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                            :: &
> nmin
real*8                             :: &
> vol
integer, pointer                   :: &
> idmin(:) => null ()
character(len=100), pointer        :: &
> namesp(:) => null ()
character(len=100)                 :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check if the chemical system is associated. 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------------
!% Get chemical information 
!%-----------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,nminsp=nmin, &
>                     idminsp=idmin,namesp=namesp)
if (iserror) goto 20
if (nmin==0) goto 20
!%------------------------------------------------------------
!% 
!%------------------------------------------------------------
if (ism3nch) then
> vol=this%volnch
else
> vol=1.0d0
end if 
!%------------------------------------------------------------
!% Write surface reactive area 
!%------------------------------------------------------------
call write_ &
>  (ioutput, &
>   head, &
>   namesp(idmin), &
>   this%alpha(idmin)/vol, &
>   nmin, &
>   item, &
>   iswritehead, &
>   .false.)
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (idmin,1,.false.)
call check_pointer_ (namesp,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
> 
> 
10 continue 
print *,'************************'
print *,'Nodal Chemistry:'
print *,'Namë:',this%name
print *,'Service:  write_sat_'
print *, msg
print *,'************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_min_area2_nch &
>  (this, &
>   ioutput, &
>   head, &
>   item, &
>   iswritehead, &
>   ism3nch, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Write in ascII the reactive surface of all mineral species 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in) :: this        !> Type nodal chemistry variable. 

integer, intent(in)                :: ioutput     !> Output unit.

real*8, intent(in)                 :: item        !> Item to write when start the line 

logical, intent(in)                :: iswritehead !> iswritename=true, then the name of the mineral species are written. 

character(len=*), intent(in)       :: head        !> Head 

logical, intent(in)                :: ism3nch     !> If .true. the reactive surface are in m3 of nodal chemistry

logical, intent(out)               :: iserror     !> iserror=true, there was an error. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                            :: &
> nmin
real*8                             :: &
> vol
integer, pointer                   :: &
> idmin(:) => null ()
character(len=100), pointer        :: &
> namesp(:) => null ()
character(len=100)                 :: &
> msg 
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check if the chemical system is associated. 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------------
!% Get chemical information 
!%-----------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,nminsp=nmin, &
>                     idminsp=idmin,namesp=namesp)
if (iserror) goto 20
if (nmin==0) goto 20
!%------------------------------------------------------------
!% 
!%------------------------------------------------------------
if (ism3nch) then
> vol=this%volnch
else
> vol=1.0d0
end if 
!%------------------------------------------------------------
!% Write surface reactive area 
!%------------------------------------------------------------
call write_ &
>  (ioutput, &
>   head, &
>   namesp(idmin), &
>   this%alpha(idmin)/vol, &
>   nmin, &
>   item, &
>   iswritehead, &
>   .false.)
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (idmin,1,.false.)
call check_pointer_ (namesp,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
> 
> 
10 continue 
print *,'************************'
print *,'Nodal Chemistry:'
print *,'Namë:',this%name
print *,'Service:  write_sat_'
print *, msg
print *,'************************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_volfracmin_nch &
>  (this, &
>   ioutput, &
>   iserror, &
>   iswritename, &
>   realitem, &
>   integeritem)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Write in ascii the %volmin
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)              :: this        !> Type nodal chemistry variable. 

integer, intent(in)                             :: ioutput     !> Output unit 

logical, intent(out)                            :: iserror     !> iserror=true, there was an error. 

logical, intent(in), optional                   :: iswritename !> If .true. then write the name species

real*8, intent(in), optional                    :: realitem    !> Real item  

integer, intent(in), optional                   :: integeritem !> Integer item 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License: 
!
!-------------------------------------------------------------------------
integer                            :: &
> i, &
> item, &
> nmin 
real*8                             :: &
> vol, &
> voltot 
real*8, pointer                    :: &
> volmin(:) => null ()
character(len=100)                 :: &
> msg, &
> name 
integer, pointer                  :: &
> idmin(:) => null ()
character(len=100), pointer        :: & 
> namesp(:) => null ()
logical                            :: &
> haveiswritename, &
> haverealitem, &
> haveintegeritem, &
> iswrite
!-------------------------------------------------------------------------
!
!>   $code
!

!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
haveiswritename=present(iswritename)
haverealitem=present(realitem)
haveintegeritem=present(integeritem)
!%------------------------------------------------------------
if (haveiswritename) then
> iswrite=iswritename
else
> if (item<=1) then
>  iswrite=.true.
> else
>  iswrite=.false. 
> end if
end if
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check if the chemical system is associated
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
!% Get information of the chemical system 
!%------------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,idminsp=idmin,nminsp=nmin,namesp=namesp)
if (iserror) goto 20
!%-----------------------------------------------------------
!% If nmin=0 then return 
!%-----------------------------------------------------------
if (nmin==0) goto 20 
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (volmin,nmin,.true.)
voltot=0.0d0 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
do i=1,nmin

> vol=this%c(idmin(i))*this%omgwfree
> name=namesp(idmin(i))
> call change_chem_unit_ (this%pchemsys,vol,name,'m','m3',iserror)

> if (iserror) goto 20 

> volmin(i)=vol

> voltot=voltot+vol


end do
!%-----------------------------------------------------------
!% 
!%-----------------------------------------------------------
volmin=(volmin/this%volnch)
voltot=voltot/this%volnch
!%-----------------------------------------------------------
!% Write total of minerals  
!%-----------------------------------------------------------
!%write (ioutput,*) '----------------------------------------------------------'
!%write (ioutput,1) 'Total fractional volume of minerals:',voltot,'[m3/m3 rock]'
!%write (ioutput,*) '----------------------------------------------------------'
!%-----------------------------------------------------------
if (haverealitem) then  
>  call write_ (ioutput,'',namesp(idmin),volmin,nmin,realitem,iswrite,.false.)
end if
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
if (haveintegeritem) then
>  item=integeritem 
else
>  item=0 
end if
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
call write_ (ioutput,'',namesp(idmin),volmin,nmin,item,iswrite,.false.)
!%-----------------------------------------------------------
20 continue
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (volmin,1,.false.)
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (idmin,1,.false.)
!%-----------------------------------------------------------
if (iserror) goto 10 
!%-----------------------------------------------------------
return
> 
1 format(a37,f10.4,a12)
> 
10 continue 
print *,'*********************'
print *,'Nodal Chemistry:'
print *,'Namë:',this%name
print *,'Service:  write_'
print *, msg
print *,'********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_tot_nch &
>  (this, &
>   ioutput, &
>   iswritehead, &
>   iserror, &
>   integeritem, &
>   realitem)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Write in ascII file the total concentrations 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in) :: this         !> Type nodal chemistry variable. 

integer, intent(in)                :: ioutput      !> Output unit.

logical, intent(in)                :: iswritehead  !> iswritename=true, then the name of the mineral species are written. 

logical, intent(out)               :: iserror      !> iserror=true, there was an error. 

integer, intent(in), optional      :: integeritem  !> Item to write when start the line 

real*8, intent(in), optional       :: realitem     !> Item to write when start the line 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                            :: &
> ncomp, &
> nsp, &
> nmob
real*8, pointer                    :: &
> u(:) => null (), &
> cmob(:,:) => null ()
character(len=100), pointer        :: &
> namecomp(:) => null ()
character(len=100)                 :: &
> msg 
logical   :: &
> haveintegeritem, &
> haverealitem
!-------------------------------------------------------------------------
!
!>   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Check optional arguments
!%------------------------------------------------------------
haveintegeritem=present(integeritem)
haverealitem=present(realitem)
!%------------------------------------------------------------
!% Check if the chemical system is associated. 
!%------------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
!% Get information to chemical system 
!%------------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,namebase=namecomp,numbase=ncomp,numsp=nsp)
if (iserror) goto 20 
!%------------------------------------------------------------
call select_cmob_ (this%pchemsys,cmob,nmob,nsp,this%c,iserror)
!%---------------------------------------------------------
!% Compute u
!%---------------------------------------------------------
call make_lin_trf_ (this%pchemsys,u,cmob(1:nsp,1),0,iserror)
!%------------------------------------------------------------
!% Write components
!%------------------------------------------------------------
if (haveintegeritem) then  
>  call write_ &
>   (ioutput, &
>    '', &
>    namecomp, &
>    u, &
>    ncomp, &
>    integeritem, &
>    iswritehead, &
>   .false.)
end if 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
!%-----------------------------------------------------------
if (haverealitem) then  
>  call write_ &
>   (ioutput, &
>    '', &
>    namecomp, &
>    u, &
>    ncomp, &
>    realitem, &
>    iswritehead, &
>   .false.)
end if 
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (u,1,.false.)
call check_pointer_ (cmob,1,1,.false.)
call check_pointer_ (namecomp,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
> 
> 
10 continue 
print *,'************************'
print *,'Nodal Chemistry:'
print *,'Namë:',this%name
print *,'Service:  write_sat_'
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
subroutine set_from_setre_nch &
>  (this, &
>   nchk, &
>   iscompzchanged, &
>   setre, &
>   nsp, &
>   dtime, &
>   iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Set the concentration of equilibrium mineral species from the equilibrium term of the transport equation (Set re)
!> where Se is the stoichiometric matrix corresponding to equilibrium reactions, and re [mol/kgw/s] is the vector of equilibrium reaction rates. 
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout):: this           !> Type nodal chemistry variable (in k+1). 

type(t_nodalchemistry), intent(in)   :: nchk           !> Type nodal chemistry variable (in k). 

integer, intent(in)                  :: nsp            !> Number of species

real*8, intent(in)                   :: dtime          !> Time increment [s]

real*8, intent(inout), dimension(nsp):: setre          !> Set re [nsp]

logical, intent(out)                 :: iscompzchanged !> iscompzchanged=true, then the components definition was changed. 

logical, intent(out)                 :: iserror        !> iserror=true, there was a error. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)                 :: &
> msg
integer                            :: &
> ndim 
!-------------------------------------------------------------------------
!
!>   $code
!
!%-----------------------------------------------------------
msg=''
iserror=.false.
!%-----------------------------------------------------------
!% Check if the chemical system is associated 
!%-----------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%-----------------------------------------------------------
!% Check then number of species 
!%-----------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,numsp=ndim)
if (iserror) goto 10
if (ndim/=nsp) then
> msg='Error in number of species'
> goto 10
end if
!%-----------------------------------------------------------
!% Update the chemical system with the temperature
!%-----------------------------------------------------------
> call update_ (this%pchemsys,this%temp,iserror)
> if (iserror) goto 10
!%-----------------------------------------------------------
!% Specia equilibrium minerals according Set re the nodal 
!% chemistry in k+1
!%-----------------------------------------------------------
call specia_ &
>   (this%pchemsys, &
>    this%temp, &
>    this%c, &
>    this%g, &
>    nchk%c, & 
>    iscompzchanged, &
>    setre, &
>    nsp, &
>    dtime, &
>    this%hashcompz, &
>    this%omgwfree, &
>    nchk%omgwfree, &
>    iserror)
if (iserror) goto 10  
!%-----------------------------------------------------------
!% Storage the vector setre 
!%-----------------------------------------------------------
this%setre=setre
!%-----------------------------------------------------------
return
10 continue 
print *,'**************************'
print *,'Nodal Chemistry:'
print *,'Namë:',this%name
print *,'Service:  set_from_setre_'
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
subroutine copy_nch &
>  (targetobj, &
>   sourceobj)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Copy the nodal chemistry object in other nodal chemistry 
!> object
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(in)         :: sourceobj  !> Type nodal chemistry variable. 

type(t_nodalchemistry), target, intent(out):: targetobj  !> Type nodal chemistry variable. 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                    :: &
> ndim1, &
> ndim2 
!-------------------------------------------------------------------------
!
!>   $code
!
!%-----------------------------------------------------------
!% Copy concentration vector 
!%-----------------------------------------------------------
if (associated(sourceobj%c)) then
> ndim1=size(sourceobj%c)
> call check_pointer_ (targetobj%c,ndim1,.true.)
> targetobj%c = sourceobj%c
end if
!%-----------------------------------------------------------
!% Copy the factor c vector 
!%-----------------------------------------------------------
if (associated(sourceobj%facc)) then
> ndim1=size(sourceobj%facc)
> call check_pointer_ (targetobj%facc,ndim1,.true.)
> targetobj%facc = sourceobj%facc
end if
!%-----------------------------------------------------------
!% Copy activity coefficients vector
!%-----------------------------------------------------------
if (associated(sourceobj%g)) then
> ndim1=size(sourceobj%g)
> call check_pointer_ (targetobj%g,ndim1,.true.)
> targetobj%g = sourceobj%g
end if
!%-----------------------------------------------------------
!% Copy total site for sorption capacitances and specific 
!% surface area
!%-----------------------------------------------------------
if (associated(sourceobj%txoh)) then
> ndim1=size(sourceobj%txoh,1)
> ndim2=size(sourceobj%txoh,2)
> call check_pointer_ (targetobj%txoh,ndim1,ndim2,.true.)
> call check_pointer_ (targetobj%txoh0,ndim1,ndim2,.true.)
> call check_pointer_ (targetobj%idtxohmin,ndim2,.true.)
> call check_pointer_ (targetobj%capint,ndim1,ndim2,.true.)
> call check_pointer_ (targetobj%capext,ndim1,ndim2,.true.)
> call check_pointer_ (targetobj%spsurfarea,ndim1,ndim2,.true.)
> targetobj%txoh = sourceobj%txoh
> targetobj%txoh0 = sourceobj%txoh0
> targetobj%idtxohmin = sourceobj%idtxohmin
> targetobj%capint = sourceobj%capint
> targetobj%capext = sourceobj%capext
> targetobj%spsurfarea = sourceobj%spsurfarea
end if
!%-----------------------------------------------------------
!% Copy alpha 
!%-----------------------------------------------------------
if (associated(sourceobj%alpha)) then
> ndim1=size(sourceobj%alpha,1)
> call check_pointer_ (targetobj%alpha,ndim1,.true.)
> targetobj%alpha = sourceobj%alpha
end if
!%-----------------------------------------------------------
!% Copy alpha0 
!%-----------------------------------------------------------
if (associated(sourceobj%alpha0)) then
> ndim1=size(sourceobj%alpha0,1)
> ndim2=size(sourceobj%alpha0,2)
> call check_pointer_ (targetobj%alpha0,ndim1,ndim2,.true.)
> targetobj%alpha0 = sourceobj%alpha0
end if
!%------------------------------------------------------------
!% Copy sktrk vector 
!%------------------------------------------------------------
if (associated(sourceobj%sktrk)) then
> ndim1=size(sourceobj%sktrk)
> call check_pointer_ (targetobj%sktrk,ndim1,.true.)
> targetobj%sktrk = sourceobj%sktrk
end if
!%------------------------------------------------------------
!% Copy setre vector 
!%------------------------------------------------------------
if (associated(sourceobj%setre)) then
> ndim1=size(sourceobj%setre)
> call check_pointer_ (targetobj%setre,ndim1,.true.)
> targetobj%setre = sourceobj%setre
end if
!%-----------------------------------------------------------
!% If derivatives are sotored, then copy derivatives
!%-----------------------------------------------------------
targetobj%isderstored=sourceobj%isderstored
if (sourceobj%isderstored) then
> if (associated(sourceobj%dc)) then
>   ndim1=size(sourceobj%dc,1)
>   ndim2=size(sourceobj%dc,2)
>   call check_pointer_ (targetobj%dc,ndim1,ndim2,.true.)
>   targetobj%dc = sourceobj%dc
> end if
> 
> if (associated(sourceobj%dsktrk)) then
>   ndim1=size(sourceobj%dsktrk,1)
>   ndim2=size(sourceobj%dsktrk,2)
>   call check_pointer_ (targetobj%dsktrk,ndim1,ndim2,.true.)
>   targetobj%dsktrk = sourceobj%dsktrk
> end if
> 
> if (associated(sourceobj%d2c)) then
>   ndim1=size(sourceobj%d2c,1)
>   ndim2=size(sourceobj%d2c,2)
>   call check_pointer_ (targetobj%d2c,ndim1,ndim2,.true.)
>   targetobj%d2c = sourceobj%d2c
> end if
> 
end if
!%-----------------------------------------------------------
!% Copy remainder attributtes 
!%-----------------------------------------------------------
targetobj%name=sourceobj%name
targetobj%hashcompz=sourceobj%hashcompz
targetobj%temp=sourceobj%temp
targetobj%ionstr=sourceobj%ionstr
targetobj%omgwfree=sourceobj%omgwfree
targetobj%volnch=sourceobj%volnch
targetobj%volgas=sourceobj%volgas
targetobj%pgas=sourceobj%pgas
targetobj%pliq=sourceobj%pliq
targetobj%liqdens=sourceobj%liqdens
targetobj%nchemiter=sourceobj%nchemiter 
targetobj%ispgasconst=sourceobj%ispgasconst 
!%---------------------------------------------------------
!% Copy or assign the chemical system pointer 
!%---------------------------------------------------------
targetobj%islockchemsys=sourceobj%islockchemsys
if (targetobj%islockchemsys) then 
>  allocate (targetobj%pchemsys)
>  call create_ (targetobj%pchemsys)
>  targetobj%pchemsys=sourceobj%pchemsys
else
>  targetobj%pchemsys => sourceobj%pchemsys 
end if 
!%---------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine begin_element_handler (name,attributes)
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
> 
> call read_xml_loc_ (name,attributes)
> 
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine pcdata (name)
character(len=*)               :: &
> name
type(dictionary_t)             :: &
> attributes
> 
> 
call read_xml_loc_ (name,attributes)
> 
> 
> 
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc_nch  &
>  (name, &
>   attributes, &
>   this, &
>   iserror, &
>   isconvergence)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description:
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), optional:: this         !> Type nodal chemistry variable. 

character(len=*), intent(in)    :: name

type(dictionary_t), intent(in)  :: attributes

logical, intent(out), optional  :: iserror

logical, intent(out), optional  :: isconvergence
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer        :: &
> status, &
> n
character(len=100)             :: &
> id
character(len=300), save       :: &
> namefilechemsys, &
> namench
logical        :: &
> havethis
real*8         :: &
> reallocal(1)
integer        :: &
> integerlocal(1)
integer, save  :: &
> i, &
> isp, &
> numpri
real*8, save            :: &
> temp
real*8                  :: &
> cputime 
integer, pointer, save  :: &
> icon(:) => null ()
real*8, pointer, save   :: &
> cguess(:) => null (), &
> ctot(:) => null ()
character(len=100), pointer, save :: &
> unit(:) => null (), &
> namesp(:) => null (), &
> constraint(:) => null ()
logical, save                 :: &
> bechemsys, &
> benamefilechemsys, &
> isderstored
type(t_chemicalsystem), pointer:: &
> chemsys
character(len=100)            :: &
> msg
integer, parameter            :: &
> mxdim=100 
!-------------------------------------------------------------------------
!
!>   $code
!
msg=''
!%---------------------------------------------------------------
havethis=present(this)
call lowercase (name)
!%---------------------------------------------------------------
!%---------------------------------------------------------------
select case (havethis)
> 
case (.true.)
iserror=.false.
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
print *,'=======> Reading Nodal Chemistry object'
!%---------------------------------------------------------------
!> If in the file there is namefile of new chemical system then 
!> create new chemical system
!%---------------------------------------------------------------
>  if (bechemsys) then
>    allocate (chemsys)
>    call create_ (chemsys)
>    call read_xml_ (chemsys, namefilechemsys, iserror)
>    if (iserror) goto 20
>  end if
!%-------------------------------------------------
!> Set nodal chemistry object
!%-------------------------------------------------
> call set_ &
>  (this, &
>   iserror, &
>   chemsys=chemsys, &
>   isderstored=isderstored, &
>   name=namench)
!%-------------------------------------------------
!> Specia fron definition of solution
!%-------------------------------------------------
>  print *,'=======> Starting speciation'
>  call set_from_solution_type_ &
>  (this, &
>   namesp(1:numpri), &
>   ctot(1:numpri), &
>   cguess(1:numpri), &
>   unit(1:numpri), &
>   icon(1:numpri), &
>   constraint(1:numpri), &
>   numpri, &
>   temp, &
>   isconvergence, & 
>   cputime, &
>   iserror)
>  if (iserror) then
>     msg='Error when calling set_from_solution_type_'
>     goto 20
>  end if
>  if (.not.isconvergence) then
>     msg='Convergence problems when calling specia_from_solution_type_'
>     goto 20
>  end if
>  print *,'===> CPU time consumed:',cputime,' [s]'
>  print *,'===> Number of iterations:',this%nchemiter
>  print *,'=======> Finishing speciation'
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
print *,'=======> Finishing reading Nodal Chemistry object'
!%---------------------------------------------------------------
20  continue 
!%---------------------------------------------------------------
!% Deallocate local pointers 
!%---------------------------------------------------------------
>    call check_pointer_ (icon,1,.false.)
>    call check_pointer_ (cguess,1,.false.)
>    call check_pointer_ (ctot,1,.false.)
>    call check_pointer_ (unit,1,.false.)
>    call check_pointer_ (namesp,1,.false.)
>    call check_pointer_ (constraint,1,.false.)
>    if (iserror) goto 10
!%---------------------------------------------------------------
!%---------------------------------------------------------------
case default
> 
>  select case (name)
!%----------------
>  case ('nodalchemistry')
>    call check_pointer_ (icon,mxdim,.true.)
>    call check_pointer_ (cguess,mxdim,.true.)
>    call check_pointer_ (ctot,mxdim,.true.)
>    call check_pointer_ (unit,mxdim,.true.)
>    call check_pointer_ (namesp,mxdim,.true.)
>    call check_pointer_ (constraint,mxdim,.true.)
>    namefilechemsys='' 
	bechemsys=.false.
	isderstored=.false.
>    benamefilechemsys=.false.
>    isp=0
>    numpri=0
>    id=' '
>    call get_value (attributes,"name", id, status)
>    namench=id
>    id=' '
>    call get_value (attributes,"temp", id, status)
>    n=0
>    call build_data_array (id,reallocal,n)
>    temp = reallocal(1)
	id=' '
>    call get_value (attributes,"derivatestored", id, status)
>    if (id=='yes') then
	 isderstored=.true.
	end if
!%-----------------
> case ('chemicalsystem')
>    bechemsys=.true.
> case ('file')
>    if (bechemsys) then
>      benamefilechemsys=.true.
>      id=''
	  call get_value (attributes,"name", id, status)
>      namefilechemsys=id	   
>    end if
!%----------------
>  case ('species')
>    isp=isp+1
	numpri=isp
>    id=' '
>    call get_value (attributes,"name", id, status)
>    namesp(isp)=id
>    id=' '
>    call get_value (attributes,"cguess", id, status)
>    n=0
>    call build_data_array (id,reallocal,n)
>    cguess(isp) = reallocal(1)
>    id=' '
>    call get_value (attributes,"unit", id, status)
>    unit(isp)=id
>    if(status.lt.0.or.id.eq.' ') unit(isp)='m'
>    id=' '
>    call get_value (attributes,"icon", id, status)
>    if (status/=0) then
>     icon(isp)=0
>    else
>     select case (id)
>     case('tot')
>      icon(isp)=1
>     case('chgbal')
>      icon(isp)=2
>     case('act')
>      icon(isp)=3
>     case('eqgas')
>      icon(isp)=5
>     case('eqmin')
>      icon(isp)=4
>     case default
>      goto 10
>     end select
>    end if
>    if (icon(isp)/=0) then
>     id=' '
>     call get_value (attributes,"ctot", id, status)
>     n=0
>     call build_data_array (id,reallocal,n)
>     ctot(isp) = reallocal(1)
>     id=' '
>     call get_value (attributes,"constraint", id, status)
>     constraint(isp)=id
>    end if
> end select 

end select
!%----------------------------------------------------------
!%----------------------------------------------------------
return
> 
10 continue 
print *,'******************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: read_xml_'
print *, msg
print *,'******************'
iserror=.true.
return
> 
> 
> 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine reaction_path_ev_dil_water1_nch &
>  (this, &
>   ntime, &
>   nstep, &
>   iserror, &
>   ioutput, &
>   temprank)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Make the evaporation of the water not impossing the mineral equilibrium
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)            :: this !> Type nodal chemistry variable. 

integer, intent(in)                              :: nstep

real*8, intent(in)                               :: ntime

logical, intent(out)                             :: iserror

integer, intent(in), optional                    :: ioutput

real*8, intent(in), optional                     :: temprank(2) 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                          :: &
> i, &
> j, &
> nw, &
> ipos1, &
> ipos2, &
> nminsp, &
> istep, &
> nsp, &
> ntxoh, &
> nsurf
real*8, pointer                                 :: &
> simin(:,:) => null (), &
> progress(:) => null ()
character(len=100), pointer                     :: &
> namesp(:) => null ()
integer, pointer                                :: &
> idminsp(:) => null ()
real*8                                          :: &
> dtmh2o, &
> mh2o, &
> sum, &
> icputime1, &
> icputime2, &
> dtcputime
type(t_nodalchemistry), pointer                 :: &
> nch
logical                                         :: &
> haveioutput, &
> havetemprank, &
> isconvergence, &
> isevaporation, &
> writehead
character(len=100)                              :: &
> msg, &
> head
integer, parameter                              :: &
> mxstep=1000
real*8, parameter                               :: &
> expo=1.05d0 
!-------------------------------------------------------------------------
!
!>   $code
!
> 

> 

!%---------------------------------------------------------
msg=''
iserror=.false.
!%---------------------------------------------------------
!% Check if the chemical system is associated 
!%---------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%---------------------------------------------------------
isevaporation=.false.
!%---------------------------------------------------------
haveioutput=present(ioutput)
havetemprank=present(temprank)
!%---------------------------------------------------------
!% If ntime<0, then evaporation processes
!%---------------------------------------------------------
if (ntime<0.0d0) isevaporation=.true.
!%---------------------------------------------------------
allocate (nch)
call create_ (nch)
nch=this
!%---------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,wspindex=nw, &
>                     numsp=nsp,numsites=ntxoh,numsurf=nsurf, &
>                     idminsp=idminsp,nminsp=nminsp, &
>                     namesp=namesp)
if (iserror) goto 20
!%---------------------------------------------------------
!% Get the global position of the aqueous species in the 
!% concentration vector 
!%---------------------------------------------------------
call get_aq_sps_pos_ (this%pchemsys,ipos1,ipos2,iserror)
if (iserror) goto 20
!%---------------------------------------------------------
call check_pointer_ (simin,nsp,nstep,.true.)
call check_pointer_ (progress,nstep,.true.)
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
if (haveioutput) then
> 
> 
>     write(ioutput,*) '----------------------------------------'
>     write(ioutput,*) 'Beginning of batch-reaction calculations'
>     write(ioutput,*) '----------------------------------------'
>     if (isevaporation) then
>       write(ioutput,*) 'Solution evaporated:',-ntime,'times'
>     else
>       write(ioutput,*) 'Solution diluted:',ntime,'times'
>     end if
>     write(ioutput,*) '----------------------------------------'
> 
>     call write_ &
>     (nch%pchemsys, &
>      nch%c, &
>      nch%g, &
>      nch%temp, &
>      nch%ionstr, &
>      nsp, &
>      ioutput, &
>      nch%omgwfree, &
>      msg, &
>      iserror)
> 
end if
!%--------------------------------------------------------------------
!% Determine initial cpu time
!%--------------------------------------------------------------------
if (haveioutput) then
> call cpu_time (icputime1)
end if
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
if (isevaporation) then
> mh2o = -(nch%omgwfree + nch%omgwfree/ntime)
else
> mh2o = nch%omgwfree*ntime - nch%omgwfree
end if
dtmh2o = mh2o/real(nstep)
!%---------------------------------------------------------
sum=0.0d0
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!% Start loop for steps
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
do i=1,nstep
> 
>  sum = sum + dtmh2o
> 
>  progress(i)=this%omgwfree/(this%omgwfree+sum)
> 
!%--------------------
>  if (havetemprank) then
>   nch%temp=temprank(2)-temprank(1)
>   nch%temp=-(nch%temp/ntime)*progress(i)+temprank(1)
>   call set_ (nch,iserror,temp=nch%temp)
>   if (iserror) goto 20
>  end if
!%--------------------
>  call specia_ &
>   (nch%pchemsys, &
>    nch%temp, &
>    nch%c, &
>    nch%g, &
>    nch%ionstr, &
>    simin(:,i), &
>    ' ', &
>    0.0d0, &
>    dtmh2o, &
>    nch%txoh, &
>    nch%capint, &
>    nch%capext, &
>    nch%spsurfarea, &
>    nsp, &
>    nsurf, &
>    ntxoh, &
>    nch%omgwfree, &
>    isconvergence, &
>    iserror)
> 
> 
> if (iserror) then
>  msg='Error when call specia_'
>  goto 20
> end if
> 
> if (.not.isconvergence) then
>  msg='Convergence problems'
>  goto 20
> end if
> 
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
!% Write results
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
>  if (haveioutput) then
> 
>      if (isevaporation) then
>        write(ioutput,*) '----------------------------------'
>        write(ioutput,*) 'Reaction Evaporation step', i
>        write(ioutput,*) '----------------------------------'
>      else
>        write(ioutput,*) '----------------------------------'
>        write(ioutput,*) 'Reaction Dilution step', i
>        write(ioutput,*) '----------------------------------'
>      end if
> 
>       call write_ &
>     (nch%pchemsys, &
>      nch%c, &
>      nch%g, &
>      nch%temp, &
>      nch%ionstr, &
>      nsp, &
>      ioutput, &
>      nch%omgwfree, &
>      msg, &
>      iserror, &
>      simin(:,i))
> 
>      if (iserror) goto 20
>      
	  if (iserror) then
>         msg='Error when calling write_'
>         goto 10 
	  end if
> 
> 
> 
> 
>  end if
> 
> 
> 
> 
> 
end do
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!% Finish loop fro steps 
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%--------------------------------------------------------------------
!% Determine final cputime and write cpu time consumed
!%--------------------------------------------------------------------
if (haveioutput) then
> call cpu_time (icputime2)
> dtcputime=icputime2-icputime1
> write(ioutput,*) '----------------------------------'
> write(ioutput,*) 'CPU time consumed:',dtcputime,'[s]'
> write(ioutput,*) '----------------------------------'
end if
!%--------------------------------------------------------------------
!% Copy the nodal chemistry object 
!%--------------------------------------------------------------------
this=nch
!%--------------------------------------------------------------------
!% Destroy and deallocate the local chemistry object 
!%--------------------------------------------------------------------
call destroy_ (nch)
deallocate (nch)
!%--------------------------------------------------------------------
if (haveioutput.and.nminsp>0) then
> 
> if (isevaporation) then
>   head='Evolution of saturation indices vs. evaporation progress'
> else
>   head='Evolution of saturation indices vs. dilution progress'
> end if
> 
> simin=dlog10(simin)
> 
> do i=1,nstep
> 
>     if (i==1) then
>        writehead=.true.
>     else
>      writehead=.false.
>       end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namesp(idminsp), &
>       simin(idminsp,i), &
>       nminsp, &
>       progress(i), &
>       writehead, &
>       writehead)
> 
> 
> end do
> 
end if
20 continue 
!%--------------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------------
call check_pointer_ (simin,1,1,.false.)
call check_pointer_ (idminsp,1,.false.)
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (progress,1,.false.)
if (iserror) goto 10
!%--------------------------------------------------------------------
return
> 
10 continue 
print *,'************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: reaction_'
print *, msg
print *,'************************'
iserror=.true.
return
> 
1 format(a15,f10.1,a1)
2 format(a4,f7.2,a1)
3 format(a15,e15.7)
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_min_area_nch &
>   (this, &
>    iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Update the reactive surface of mineral 
!
!>   $Arguments:
!
> 
> 
type (t_nodalchemistry), intent(inout)          :: this     !> Type nodal chemistry variable. 

logical, intent(out)                            :: iserror  !> iserror=true, then there was an error

!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)                              :: &
> msg
!-------------------------------------------------------------------------
!
!>   $code
!
!%---------------------------------------------------------
!% Check if the chemical system is associated 
!%---------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
this%alpha=this%alpha/this%omgwfree
this%alpha0=this%alpha0/this%omgwfree
!%------------------------------------------------------------
!% Call the service implemented in the chemical system 
!%------------------------------------------------------------
call update_min_area_ &
>   (this%pchemsys, &
>    this%alpha, &
>    this%alpha0(2,:), &
>    this%c, &
>    this%alpha0(1,:), &
>    msg, &
>    iserror)
!%------------------------------------------------------------    
this%alpha=this%alpha*this%omgwfree
this%alpha0=this%alpha0*this%omgwfree    
if (iserror) goto 10 
!%------------------------------------------------------------
return
10 continue 
print *,'***************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: update_min_area_'
print *, msg
print *,'***************************'
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_txoh_nch &
>   (this, &
>    iserror)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Update the txoh  
!
!>   $Arguments:
!
> 
> 
type (t_nodalchemistry), intent(inout)          :: this     !> Type nodal chemistry variable. 

logical, intent(out)                            :: iserror  !> iserror=true, then there was an error

!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
character(len=100)                              :: &
> msg
!-------------------------------------------------------------------------
!
!>   $code
!
!%---------------------------------------------------------
!% Check if the chemical system is associated 
!%---------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%------------------------------------------------------------
this%txoh=this%txoh/this%omgwfree
this%txoh0=this%txoh0/this%omgwfree
!%------------------------------------------------------------
!% Call the service implemented in the chemical system 
!%------------------------------------------------------------
call update_txoh_ &
>   (this%pchemsys, &
>    this%txoh, &
>    this%txoh0, &
>    this%idtxohmin, &
>    this%c, &
>    this%alpha0(1,:), &
>    msg, &
>    iserror)
!%------------------------------------------------------------    
this%txoh=this%txoh*this%omgwfree
this%txoh0=this%txoh0*this%omgwfree    
if (iserror) goto 10 
!%------------------------------------------------------------
return
10 continue 
print *,'***************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: update_txoh_'
print *, msg
print *,'***************************'
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine reaction_path_ev_dil_water2_nch &
>  (this, &
>   ntime, &
>   nstep, &
>   fraccryst, &
>   isequilibrate, &
>   iswcompbal, &
>   iserror, &
>   ioutput, &
>   namespout, &
>   temprank, &
>   isconvergence)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Make the evaporation of the water impossing the mineral equilibrium.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)                :: this           !> Type nodal chemistry variable. 

integer, intent(in)                                  :: nstep

real*8, intent(in)                                   :: fraccryst      !> Cristalline fractionation 

real*8, intent(in)                                   :: ntime

logical, intent(in)                                  :: isequilibrate 

logical, intent(in)                                  :: iswcompbal 

logical, intent(out)                                 :: iserror

character(len=*), intent(in), optional, dimension(:) :: namespout

real*8, intent(in), optional, dimension(2)           :: temprank

integer, intent(in), optional                        :: ioutput

logical, intent(out), optional                       :: isconvergence 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                          :: &
> i, &
> j, &
> nw, &
> icw, &
> ipos1, &
> ipos2, &
> nminsp, &
> istep, &
> nsp, &
> ntxoh, &
> nsurf, &
> npri, &
> iter, &
> nspout, &
> isp, &
> nwout, &
> ntempstep
real*8, pointer                   :: &
> progress(:) => null (), &
> molminstep(:,:) => null (), &
> tk(:) => null (), &
> tk1(:) => null (), &
> simin(:,:) => null (), &
> siloc(:) => null (), &
> molinit(:) => null (), &
> molfin(:) => null (), &
> molaq(:) => null (), &
> molads(:) => null (), &
> molprec(:) => null (), &
> errorcomp(:) => null (), &
> cout(:,:) => null (), &
> actout(:,:) => null (), &
> masswfree(:) => null (), &
> chemparam(:,:) => null ()
character(len=100), pointer       :: &
> namesp(:) => null (), &
> nameparam(:) => null ()
integer, pointer                  :: &
> idminsp(:) => null (), &
> idpri(:) => null (), &
> idspout(:) => null ()
real*8                            :: &
> dtmh2o, &
> mh2o, &
> sum, &
> error, &
> icputime1, &
> icputime2, &
> dtcputime, &
> factoromgw
type(t_nodalchemistry), pointer   :: &
> nch => null (), &
> nchold => null ()
logical                           :: &
> haveioutput, &
> havenamespout, &
> havetemprank, &
> havepgas, &
> haveisconvergence, &
> isconv, &
> evaporation, &
> writehead
integer, parameter                :: &
> mxiter=100
character(len=100)                :: &
> msg, &
> head
integer, parameter                :: &
> nparam=2
real*8, parameter                 :: &
> expo=1.05d0, &
> r0=0.0d0, &
> r1=1.0d0
!-------------------------------------------------------------------------
!
!>   $code
!
> 
!%---------------------------------------------------------
msg=''
iserror=.false.
!%---------------------------------------------------------
!% Check if the chemical system is associated 
!%---------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%---------------------------------------------------------
evaporation=.false.
!%---------------------------------------------------------
!% Check optional arguments
!%---------------------------------------------------------
haveioutput=present(ioutput)
havenamespout=present(namespout)
havetemprank=present(temprank)
haveisconvergence=present(isconvergence)
!%---------------------------------------------------------
if (ntime<r0) evaporation=.true.
!%---------------------------------------------------------
!% Allocate and create local nodal chemistry object 
!%---------------------------------------------------------
allocate (nch)
allocate (nchold)
call create_ (nch)
call create_ (nchold)
!%---------------------------------------------------------
!% Copy information in nodal chemistry object 
!%---------------------------------------------------------
nch=this
!%---------------------------------------------------------
!% Update the nodal chemistry with the temperature 
!%---------------------------------------------------------
call update_ (nch%pchemsys,nch%temp,iserror)
if (iserror) goto 20
!%---------------------------------------------------------
!% Set if the water component mass balance must be evaluated
!%---------------------------------------------------------
call set_iswcompbal_ (nch%pchemsys,iswcompbal,iserror)
if (iserror) goto 10
!%---------------------------------------------------------
!% Get chemical information 
!%---------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,wspindex=nw,wcindex=icw, &
>                     numsp=nsp,numsites=ntxoh,numsurf=nsurf, &
>                     idminsp=idminsp,nminsp=nminsp, &
>                     namesp=namesp,numbase=npri, &
>                     idbase=idpri)
!%---------------------------------------------------------
!% Check the crystalline fractionation 
!%---------------------------------------------------------
if (fraccryst<r0.or.fraccryst>r1) then
>  msg='Error, fractionation crystalline factor must be contained between 0-1'
>  goto 20
end if
!%---------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------
call check_pointer_ (progress,nstep,.true.)
call check_pointer_ (masswfree,nstep,.true.)
call check_pointer_ (chemparam,nparam,nstep,.true.)
call check_pointer_ (nameparam,nparam,.true.)
call check_pointer_ (tk,npri,.true.)
call check_pointer_ (tk1,npri,.true.)
nameparam(1)='water activity'
nameparam(2)='ionic strength'
!%---------------------------------------------------------
!% Find global index of species for output 
!%---------------------------------------------------------
nspout=0 
nwout=0 
if (haveioutput.and.havenamespout) then
> nspout=size(namespout)
> if (nspout>0) then
>  call check_pointer_ (idspout,nspout,.true.)
>  call check_pointer_ (cout,nspout,nstep,.true.)
>  call check_pointer_ (actout,nspout,nstep,.true.)
>  do i=1,nspout
>   call get_sp_index_(this%pchemsys,namespout(i),isp)
>   idspout(i)=isp
>   if (isp==nw) nwout=i
>  end do
> end if
end if
!%---------------------------------------------------------
!% Allocate local pointers  
!%---------------------------------------------------------
if (haveioutput) then
> call check_pointer_ (errorcomp,npri,.true.)
> call check_pointer_ (molinit,npri,.true.)
> call check_pointer_ (molfin,npri,.true.)
end if
!%---------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------
if (nminsp>0) then
> call check_pointer_ (molminstep,nminsp,nstep,.true.)
> call check_pointer_ (simin,nminsp,nstep,.true.)
end if
!%---------------------------------------------------------
!% Compute initial moles in the system
!%---------------------------------------------------------
call get_chem_info_(this,iserror,molaq=molaq,molads=molads,molprec=molprec)
molinit = molaq + molads + molprec
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
if (haveioutput) then
> 
> 
>     write(ioutput,*) '----------------------------------------'
>     write(ioutput,*) 'Beginning of batch-reaction calculations'
>     write(ioutput,*) '----------------------------------------'
>     if (evaporation) then
>        write(ioutput,1) 'Solution evaporated:',-ntime,'times'
>     else
>      write(ioutput,1) 'Solution diluted:',ntime,'times'
>     end if
>      write(ioutput,*) '----------------------------------------'
>      write(ioutput,1) 'Performed in ',nstep,'steps'
>      write(ioutput,*) '----------------------------------------'
>      write(ioutput,1)'Fractionation crystalline factor',fraccryst
>      write(ioutput,*) '----------------------------------------'
!%--------------------------------------------------------------------
!% Write initial nodal chemistry information 
!%--------------------------------------------------------------------  
>     call write_ (nch,ioutput,iserror)
>     if (iserror) goto 20
> 
end if
> 
!%--------------------------------------------------------------------
!% Compute the initial cpu time
!%--------------------------------------------------------------------
if (haveioutput) then
> call cpu_time (icputime1)
end if
!%---------------------------------------------------------
!% Compute increment/decrement of mass of water
!%---------------------------------------------------------
if (evaporation) then
> mh2o = -(nch%omgwfree + nch%omgwfree/ntime)
else
> mh2o = nch%omgwfree*ntime - nch%omgwfree
end if
dtmh2o = mh2o/real(nstep)
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------  
sum=r0
factoromgw = r1
!%--------------------------------------------------------------------
!%Remove mineral phases according crystalline fracctionation factor
!%--------------------------------------------------------------------
if (nminsp>0) then
>  molminstep(:,1) = fraccryst * this%c(idminsp) * this%omgwfree
>  this%c(idminsp) = (r1-fraccryst) * this%c(idminsp)
end if
!%---------------------------------------------------------
!% Compute the total analytic concentrations in k 
!%---------------------------------------------------------
call make_lin_trf_ (this%pchemsys,tk,this%c,0,iserror)
!%-----------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%Start loop for steps 
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
istep=0 
do i=1,nstep
nchold=nch
!%--------------------------------------------------------------------
!% Compute increment/decrement of mass of water 
!%--------------------------------------------------------------------
>  dtmh2o=mh2o*expo**(real(i-nstep)/real(i)) - sum
!%--------------------------------------------------------------------
!% Compute accumulative mass of water 
!%--------------------------------------------------------------------
>  sum = sum + dtmh2o 
!%--------------------------------------------------------------------
!% Add in nodal chemistry object 
!%-------------------------------------------------------------------- 
>  nch%omgwfree = this%omgwfree + sum 
!%--------------------------------------------------------------------
!% Compute evaporation/dilution progress 
!%--------------------------------------------------------------------  
>  progress(i) = this%omgwfree / nch%omgwfree
!%--------------------------------------------------------------------
!% Update the temperature 
!%--------------------------------------------------------------------
>  if (havetemprank) then
>   nch%temp=temprank(2)-temprank(1)
>   nch%temp=-(nch%temp/ntime)*progress(i)+temprank(1)
>   call set_ (nch,iserror,temp=nch%temp)
>   if (iserror) goto 20
>  end if
!%--------------------------------------------------------------------
!% Compute the total analytic concentrations 
!%--------------------------------------------------------------------
>   tk1 = progress(i) *  tk
!%--------------------------------------------------------------------
!% Correction for total mol of water component
!%--------------------------------------------------------------------   
>   if (icw>0) then
>    tk1(icw) = tk1(icw) + sum/(kgwmol*nch%omgwfree)
>   end if 
!%--------------------------------------------------------------------
!% 
!%--------------------------------------------------------------------   
>   nch%alpha=nch%alpha/nch%omgwfree
!%--------------------------------------------------------------------
!% Make the speciation according total analytic concentrations 
!%--------------------------------------------------------------------  
>   call specia_ &
>   (nch%pchemsys, &
>    nch%temp, &
>    nch%c, &
>    nch%g, &
	isequilibrate, &
	nch%ispgasconst, &
>    nch%ionstr, &
>    nch%alpha, &
>    tk1, &
>    nch%txoh/nch%omgwfree, &
>    nch%capint, &
>    nch%capext, &
>    nch%spsurfarea, &
>    nsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    0.0d0, &
>    nch%hashcompz, &
>    isconv, &
>    factoromgw, &
	nch%volgas, &
	nch%pgas, &
	msg, &
>    iserror, &
>    simin=siloc, &
>    nchemiter=nch%nchemiter, &
>    cguess=nch%c(idpri))
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------    
	if (iserror) goto 20
>    if (.not.isconv) then
>      nch=nchold
>      msg='Convergence problems'
>      if (haveioutput) then
>         write (ioutput,'(a17,i5)') 'Evaporation step:',i
>         write (ioutput,*) msg(1:20)//' when evaporation progress is',progress(i)
>      end if
>      exit
>    else
>      istep=istep+1
>    end if
!%---------------------------------------------------------------------
!% Change the alpha vector to 
!%---------------------------------------------------------------------    
	nch%alpha=nch%alpha*nch%omgwfree
!%---------------------------------------------------------------------
!% Correct the mass of free water in the nodal chemistry object 
!%---------------------------------------------------------------------
>    nch%omgwfree = nch%omgwfree * factoromgw  
!%---------------------------------------------------------------------
!% Update the reactive surface of kinetic minerals
!%---------------------------------------------------------------------
>    call update_min_area_ (nch,iserror)
>    if (iserror) goto 20 
!%---------------------------------------------------------------------
!% Correct the evaporation progress 
!%---------------------------------------------------------------------    
>    progress(i) = progress(i) / factoromgw    
!%---------------------------------------------------------------------
!% Storage the mass of free water in the evaoration/dilution step 
!%---------------------------------------------------------------------    
>    masswfree(i) = nch%omgwfree
!%---------------------------------------------------------------------
!% Storage the water activity 
!%---------------------------------------------------------------------    
>    if (nw>0) then 
>     chemparam(1,i) = nch%g(nw)
>    end if 
!%---------------------------------------------------------------------
!% Storage the ionic strength 
!%---------------------------------------------------------------------    
>    chemparam(2,i) = nch%ionstr
!%---------------------------------------------------------------------
!% Write results
!%---------------------------------------------------------------------
> if (haveioutput) then
> 
>    if (nminsp>0) then
>     simin(:,i)=siloc(idminsp)
>     molminstep(:,i) = molminstep(:,i) + nch%c(idminsp) * nch%omgwfree
>    end if
> 
> 
>    if (havenamespout.and.nspout>0) then
>     cout(:,i)=nch%c(idspout)
>     actout(:,i)=nch%c(idspout)*nch%g(idspout)
>     if (nwout>0) actout(nwout,i)=nch%g(nw)
>    end if
> 
> 
>    if (evaporation) then
>     write(ioutput,*) '----------------------------------'
>     write(ioutput,*) 'Reaction Evaporation step', i
>     write(ioutput,*) '----------------------------------'
>    else
>     write(ioutput,*) '----------------------------------'
>     write(ioutput,*) 'Reaction Dilution step', i
>     write(ioutput,*) '----------------------------------'
>    end if
> 
> 
>    write(ioutput,1) 'Progress:',progress(i),'times'
>    write(ioutput,*) '----------------------------------'
!%--------------------------------------------------------------------
!% Write nodal chemistry information 
!%-------------------------------------------------------------------- 
>    call write_(nch,ioutput,iserror)
>    if (iserror) then
>      msg='Error when calling write_'
>      if (haveioutput) then
>       write (ioutput,*) msg
>      end if
>      goto 20
>    end if
> end if
> 
> 
end do
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!% Determine final cputime and write cpu time consumed
!%--------------------------------------------------------------------
if (haveioutput) then
> call cpu_time (icputime2)
> dtcputime=icputime2-icputime1
> write(ioutput,*) '----------------------------------'
> write(ioutput,1) 'CPU time consumed:',dtcputime,'[s]'
> write(ioutput,*) '----------------------------------'
end if
!%--------------------------------------------------------------------
!% Write information about relative error
!%--------------------------------------------------------------------
if (haveioutput) then
> 
> call get_chem_info_ (nch,iserror,molaq=molaq,molads=molads,molprec=molprec)
> 
> head='Relative error in components'
> 
> molfin = molaq + molads + molprec
> 
> errorcomp = dabs((molinit - molfin) / molfin)
> 
> writehead=.true.
> 
> call write_ (ioutput,head,namesp(idpri),errorcomp,npri,'error',writehead,writehead)
> 
> 
> 
end if
!%--------------------------------------------------------------------
!% Update the last information in nodal chemistry object 
!%--------------------------------------------------------------------
this=nch
!%--------------------------------------------------------------------
!% Destroy and deallocate nodal chemistry object 
!%--------------------------------------------------------------------
call destroy_ (nch)
call destroy_ (nchold)
deallocate (nch)
deallocate (nchold)
!%--------------------------------------------------------------------
!% Write mol of minerals vs. evaporation/dilution progress
!%--------------------------------------------------------------------
if (haveioutput.and.nminsp>0) then
> 
> head='Mol of minerals vs. evaporated/diluted times'

> 
> do i=1,istep
> 
>     if (i==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namesp(idminsp), &
>       molminstep(:,i), &
>       nminsp, &
>       progress(i), &
>       writehead, &
>       writehead)
> 
> 
> end do
!%--------------------------------------------------------------------
!% Write m3 of minerals vs. evaporation/dilution progress
!%--------------------------------------------------------------------
> 
> head='Volume of minerals (m3) vs. evaporated/diluted times'
> 
> 
> do i=1,istep
> 
>     do j=1,nminsp
>        call change_chem_unit_ (this%pchemsys,molminstep(j,i),namesp(idminsp(j)),'mol','m3',iserror)
>     end do
> 
>     if (i==1) then
>        writehead=.true.
>     else
>        writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namesp(idminsp), &
>       molminstep(:,i), &
>       nminsp, &
>       progress(i), &
>       writehead, &
>       writehead)
> 
> 
> end do
!%--------------------------------------------------------------------
!% Write saturation indices of minerals vs. evaporation/dilution 
!% progress
!%--------------------------------------------------------------------
> 
> head='Evolution of saturation indices vs. evaporated/diluted times'
> 
> 
> simin(1:nminsp,1:istep)=dlog10(simin(1:nminsp,1:istep))
> 
> do i=1,istep 
> 
>     if (i==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namesp(idminsp), &
>       simin(:,i), &
>       nminsp, &
>       progress(i), &
>       writehead, &
>       writehead)
> 
> 
> end do
> 
> 
end if
!%--------------------------------------------------------------------
!% Write concentration of species vs. dilution/evaporation progress
!%--------------------------------------------------------------------
if (haveioutput.and.nspout>0) then

> head='Evolution of concentration species vs. evaporated/diluted times'

> do i=1,istep
> 
>     if (i==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       cout(1:nspout,i), &
>       nspout, &
>       progress(i), &
>       writehead, &
>       writehead)
> end do

> head='Evolution of mol of species vs. evaporated/diluted times'
>  

> do i=1,istep
> 
>     if (i==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
> 
>     cout(1:nspout,i) = cout(1:nspout,i) * masswfree(i)

>     call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       cout(1:nspout,i), &
>       nspout, &
>       progress(i), &
>       writehead, &
>       writehead)
> end do
> 

> head='Evolution of activity of species vs. evaporated/diluted times'

> do i=1,istep
> 
>     if (i==1) then
>       writehead=.true.
>     else
>       writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       actout(1:nspout,i), &
>       nspout, &
>       progress(i), &
>       writehead, &
>       writehead)
> end do
> 
> 
> head='Evolution of chemical parameters vs. evaporated/diluted times'
> 
> do i=1,istep
> 
>     if (i==1) then
>       writehead=.true.
>     else
>       writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       nameparam, &
>       chemparam(1:nparam,i), &
>       nparam, &
>       progress(i), &
>       writehead, &
>       writehead)
> end do
> 
> 
> 
end if
20 continue 
!%--------------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------------
call check_pointer_ (masswfree,1,.false.)
call check_pointer_ (idminsp,1,.false.)
call check_pointer_ (simin,1,1,.false.)
call check_pointer_ (siloc,1,.false.)
call check_pointer_ (molminstep,1,1,.false.)
call check_pointer_ (idspout,1,.false.)
call check_pointer_ (cout,1,1,.false.)
call check_pointer_ (actout,1,1,.false.)
call check_pointer_ (errorcomp,1,.false.)
call check_pointer_ (molinit,1,.false.)
call check_pointer_ (molfin,1,.false.)
call check_pointer_ (molaq,1,.false.)
call check_pointer_ (molprec,1,.false.)
call check_pointer_ (molads,1,.false.)
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (chemparam,1,1,.false.)
call check_pointer_ (nameparam,1,.false.)
call check_pointer_ (progress,1,.false.)
call check_pointer_ (tk,1,.false.)
call check_pointer_ (tk1,1,.false.)
call check_pointer_ (idpri,1,.false.)
if (iserror) goto 10
if (haveisconvergence.and..not.isconv) then
>  isconvergence=.false.
>  return
else if (haveisconvergence.and.isconv) then
>  isconvergence=.true.
>  return
end if
if (.not.isconv) stop
!%--------------------------------------------------------------------
return
> 
10 continue 
print *,'************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: reaction_path_'
print *, msg
print *,'************************'
iserror=.true.
return
> 
1 format(a21,f10.5,a6)
2 format(a34,f10.5,a20)
> 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine reaction_path_ev_dil_water3_nch &
>  (this, &
>   nstep, &
>   fraccryst, &
>   isequilibrate, &
>   iswcompbal, &
>   rhair, &
>   time, &
>   unittime, &
>   kev, &
>   iserror, &
>   ioutput, &
>   namespout, &
>   isconvergence)
> 
implicit none
!-------------------------------------------------------------------------
!
!>   $Description: Make the evaporation of the water impossing the mineral equilibrium.
!
!>   $Arguments:
!
> 
type(t_nodalchemistry), intent(inout)                :: this           !> Type nodal chemistry variable. 

integer, intent(in)                                  :: nstep

real*8, intent(in)                                   :: fraccryst      !> Cristalline fractionation 

logical, intent(in)                                  :: isequilibrate 

real*8, intent(in), dimension(nstep)                 :: rhair          !> Relative humidity of the air

real*8, intent(in), dimension(nstep)                 :: time           !> Relative humidity of the air 

real*8, intent(in)                                   :: kev            !> mass transfer coefficient, that depends of wid speed 

logical, intent(in)                                  :: iswcompbal 

character(len=*)                                     :: unittime

logical, intent(out)                                 :: iserror

character(len=*), intent(in), optional, dimension(:) :: namespout

integer, intent(in), optional                        :: ioutput

logical, intent(out), optional                       :: isconvergence 
> 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
integer                                          :: &
> i, &
> j, &
> nw, &
> icw, &
> ipos1, &
> ipos2, &
> nminsp, &
> istep, &
> nsp, &
> ntxoh, &
> nsurf, &
> npri, &
> iter, &
> nspout, &
> isp, &
> nwout, &
> ntempstep
real*8, pointer                   :: &
> progress(:) => null (), &
> molminstep(:,:) => null (), &
> tk(:) => null (), &
> tk1(:) => null (), &
> simin(:,:) => null (), &
> siloc(:) => null (), &
> molinit(:) => null (), &
> molfin(:) => null (), &
> molaq(:) => null (), &
> molads(:) => null (), &
> molprec(:) => null (), &
> errorcomp(:) => null (), &
> cout(:,:) => null (), &
> actout(:,:) => null (), &
> chemparam(:,:) => null (), &
> masswfree(:) => null ()
character(len=100), pointer       :: &
> namesp(:) => null (), &
> nameparam(:) => null ()
integer, pointer                  :: &
> idminsp(:) => null (), &
> idpri(:) => null (), &
> idspout(:) => null ()
real*8                            :: &
> mh2o, &
> sum, &
> error, &
> icputime1, &
> icputime2, &
> dtcputime, &
> factoromgw, &
> aw, &
> deltime, &
> factime
type(t_nodalchemistry), pointer   :: &
> nch => null ()
logical                           :: &
> haveioutput, &
> havenamespout, &
> havepgas, &
> haveisconvergence, &
> isconv, &
> writehead
integer, parameter                :: &
> mxiter=100
character(len=100)                :: &
> msg, &
> head
integer, parameter                :: &
> nparam=2
real*8, parameter                 :: &
> r0=0.0d0, &
> r1=1.0d0, &
> pw=0.032d0  !> [atm]
!-------------------------------------------------------------------------
!
!>   $code
!
> 
!%---------------------------------------------------------
msg=''
iserror=.false.
!%---------------------------------------------------------
!% Check if the chemical system is associated 
!%---------------------------------------------------------
if (.not.associated(this%pchemsys)) then
> msg='Error, not associated chemical system'
> goto 10
end if
!%---------------------------------------------------------
!% Check optional arguments
!%---------------------------------------------------------
haveioutput=present(ioutput)
havenamespout=present(namespout)
haveisconvergence=present(isconvergence)
!%---------------------------------------------------------
!% Check the uni time
!%---------------------------------------------------------
if (unittime=='seconds') then
>  factime=r1
else if (unittime=='days') then
>  factime=86400.0d0
else if (unittime=='years') then
>  factime=3.1536d7
else
>  msg='Error, unit time not recognized'
>  goto 10
end if
!%---------------------------------------------------------
!% Allocate and create local nodal chemistry object 
!%---------------------------------------------------------
allocate (nch)
call create_ (nch)
!%---------------------------------------------------------
!% Copy information in nodal chemistry object 
!%---------------------------------------------------------
nch=this
!%---------------------------------------------------------
!% Update the nodal chemistry with the temperature 
!%---------------------------------------------------------
call update_ (nch%pchemsys,nch%temp,iserror)
if (iserror) goto 20
!%---------------------------------------------------------
!% Set if the water component mass balance must be evaluated
!%---------------------------------------------------------
call set_iswcompbal_ (nch%pchemsys,iswcompbal,iserror)
if (iserror) goto 10
!%---------------------------------------------------------
!% Get chemical information 
!%---------------------------------------------------------
call get_chem_info_ (this%pchemsys,iserror,wspindex=nw,wcindex=icw, &
>                     numsp=nsp,numsites=ntxoh,numsurf=nsurf, &
>                     idminsp=idminsp,nminsp=nminsp, &
>                     namesp=namesp,numbase=npri, &
>                     idbase=idpri)
!%---------------------------------------------------------
!% Check the crystalline fractionation 
!%---------------------------------------------------------
if (fraccryst<r0.or.fraccryst>r1) then
>  msg='Error, fractionation crystalline factor must be contained between 0-1'
>  goto 20
end if
!%---------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------
call check_pointer_ (progress,nstep,.true.)
call check_pointer_ (masswfree,nstep,.true.)
call check_pointer_ (chemparam,nparam,nstep,.true.)
call check_pointer_ (nameparam,nparam,.true.)
call check_pointer_ (tk,npri,.true.)
call check_pointer_ (tk1,npri,.true.)
nameparam(1)='water activity'
nameparam(2)='ionic strength'
!%---------------------------------------------------------
!% Find global index of species for output 
!%---------------------------------------------------------
nspout=0 
nwout=0 
if (haveioutput.and.havenamespout) then
> nspout=size(namespout)
> if (nspout>0) then
>  call check_pointer_ (idspout,nspout,.true.)
>  call check_pointer_ (cout,nspout,nstep,.true.)
>  call check_pointer_ (actout,nspout,nstep,.true.)
>  do i=1,nspout
>   call get_sp_index_(this%pchemsys,namespout(i),isp)
>   idspout(i)=isp
>   if (isp==nw) nwout=i
>  end do
> end if
end if
!%---------------------------------------------------------
!% Allocate local pointers  
!%---------------------------------------------------------
if (haveioutput) then
> call check_pointer_ (errorcomp,npri,.true.)
> call check_pointer_ (molinit,npri,.true.)
> call check_pointer_ (molfin,npri,.true.)
end if
!%---------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------
if (nminsp>0) then
> call check_pointer_ (molminstep,nminsp,nstep,.true.)
> call check_pointer_ (simin,nminsp,nstep,.true.)
end if
!%---------------------------------------------------------
!% Compute initial moles in the system
!%---------------------------------------------------------
call get_chem_info_(this,iserror,molaq=molaq,molads=molads,molprec=molprec)
molinit = molaq + molads + molprec
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
if (haveioutput) then
> 
> 
>     write(ioutput,*) '----------------------------------------'
>     write(ioutput,*) 'Beginning of batch-reaction calculations'
>     write(ioutput,*) '----------------------------------------'
>     write(ioutput,*) '----------------------------------------'
>     write(ioutput,1) 'Performed in ',nstep,'steps'
>     write(ioutput,1) 'Total time ',time(nstep),unittime
>     write(ioutput,*) '----------------------------------------'
>     write(ioutput,1)'Fractionation crystalline factor',fraccryst
>     write(ioutput,*) '----------------------------------------'
!%--------------------------------------------------------------------
!% Write initial nodal chemistry information 
!%--------------------------------------------------------------------  
>     call write_ (nch,ioutput,iserror)
>     if (iserror) goto 20
> 
end if
> 
!%--------------------------------------------------------------------
!% Compute the initial cpu time
!%--------------------------------------------------------------------
if (haveioutput) then
> call cpu_time (icputime1)
end if
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------  
sum=r0
factoromgw = r1
!%---------------------------------------------------------
!% Compute the total analytic concentrations in k 
!%---------------------------------------------------------
call make_lin_trf_ (this%pchemsys,tk,this%c,0,iserror)
!%-----------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
!%Start loop for steps 
!%---------------------------------------------------------
!%---------------------------------------------------------
!%---------------------------------------------------------
istep=0 
do i=1,nstep
> if (i==1) then 
>  deltime=time(i)
> else
>  deltime=time(i)-time(i-1)
> end if
> call get_chem_info_(nch,iserror,actwater=aw)
!%--------------------------------------------------------------------
!% Compute increment/decrement of mass of water 
!%--------------------------------------------------------------------
> mh2o=-kev*pw*(aw-rhair(i))*deltime 
!%--------------------------------------------------------------------
!% Compute accumulative mass of water 
!%--------------------------------------------------------------------
> sum = sum + mh2o 
!%--------------------------------------------------------------------
!% Compute evaporation/dilution progress 
!%--------------------------------------------------------------------  
> progress(i) = nch%omgwfree / (nch%omgwfree+mh2o)  
!%--------------------------------------------------------------------
!% Add in nodal chemistry object 
!%-------------------------------------------------------------------- 
> nch%omgwfree = nch%omgwfree + mh2o
!%---------------------------------------------------------
!% Compute the total analytic concentrations in k 
!%---------------------------------------------------------
call make_lin_trf_ (this%pchemsys,tk,nch%c,0,iserror)
!%--------------------------------------------------------------------
!% Compute the total analytic concentrations 
!%--------------------------------------------------------------------
tk1 = progress(i) *  tk
!%--------------------------------------------------------------------
!% Correction for total mol of water component
!%--------------------------------------------------------------------   
if (icw>0) then
>    tk1(icw) = tk1(icw) + (r1-progress(i))/kgwmol
end if 
!%--------------------------------------------------------------------
!% 
!%--------------------------------------------------------------------   
nch%alpha=nch%alpha/nch%omgwfree
!%--------------------------------------------------------------------
!% Make the speciation according total analytic concentrations 
!%--------------------------------------------------------------------  
call specia_ &
>   (nch%pchemsys, &
>    nch%temp, &
>    nch%c, &
>    nch%g, &
	isequilibrate, &
	nch%ispgasconst, &
>    nch%ionstr, &
>    nch%alpha, &
>    tk1, &
>    nch%txoh/nch%omgwfree, &
>    nch%capint, &
>    nch%capext, &
>    nch%spsurfarea, &
>    nsp, &
>    npri, &
>    ntxoh, &
>    nsurf, &
>    deltime*factime, &
>    nch%hashcompz, &
>    isconv, &
>    factoromgw, &
	nch%volgas, &
	nch%pgas, &
	msg, &
>    iserror, &
>    simin=siloc, &
>    nchemiter=nch%nchemiter, &
>    cguess=nch%c(idpri))
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------    
	if (iserror) goto 20
>    if (.not.isconv) then
>      msg='Convergence problems'
>      if (haveioutput) write (ioutput,*) msg(1:20)//' in time step:',time(i)
>      exit
>    else
>      istep=istep+1
>    end if
!%---------------------------------------------------------------------
!% Change the alpha vector to 
!%---------------------------------------------------------------------    
	nch%alpha=nch%alpha*nch%omgwfree
!%---------------------------------------------------------------------
!% Correct the mass of free water in the nodal chemistry object 
!%---------------------------------------------------------------------
>    nch%omgwfree = nch%omgwfree * factoromgw  
!%---------------------------------------------------------------------
!% Update the reactive surface of kinetic minerals
!%---------------------------------------------------------------------
>    call update_min_area_ (nch,iserror)
>    if (iserror) goto 20 
!%---------------------------------------------------------------------
!% Correct the evaporation progress 
!%---------------------------------------------------------------------    
>    progress(i) = progress(i) / factoromgw    
!%---------------------------------------------------------------------
!% Storage the mass of free water in the evaoration/dilution step 
!%---------------------------------------------------------------------    
>    masswfree(i) = nch%omgwfree
!%---------------------------------------------------------------------
!% Storage the water activity 
!%---------------------------------------------------------------------    
>    if (nw>0) then 
>     chemparam(1,i) = nch%g(nw)
>    end if 
!%---------------------------------------------------------------------
!% Storage the ionic strength 
!%---------------------------------------------------------------------    
>    chemparam(2,i) = nch%ionstr
>    write(*,*) '----------------------------------'
>    write(*,*) 'Reaction Evaporation step', i
>    write(*,*) '----------------------------------'
>    write(*,1) 'Time:',time(i),unittime
>    write(*,1) 'Progress:',progress(i),'times'
>    if (mh2o<r0) then
>     write(*,*) 'Evaporation'
>    else
>     write(*,*) 'Condensation'
>    end if  
>    write(*,*) '----------------------------------'
!%---------------------------------------------------------------------
!% Write results
!%---------------------------------------------------------------------
> if (haveioutput) then
> 
>    if (nminsp>0) then
>     simin(:,i)=siloc(idminsp)
>     molminstep(:,i) = molminstep(:,i) + nch%c(idminsp) * nch%omgwfree
>    end if
> 
> 
>    if (havenamespout.and.nspout>0) then
>     cout(:,i)=nch%c(idspout)
>     actout(:,i)=nch%c(idspout)*nch%g(idspout)
>     if (nwout>0) actout(nwout,i)=nch%g(nw)
>    end if
>  
>    write(ioutput,*) '----------------------------------'
>    write(ioutput,*) 'Reaction Evaporation step', i
>    write(ioutput,*) '----------------------------------'
>  
>    write(ioutput,1) 'Time:',time(i),unittime
>    write(ioutput,1) 'Progress:',progress(i),'times'
>    write(ioutput,*) '----------------------------------'
!%--------------------------------------------------------------------
!% Write nodal chemistry information 
!%-------------------------------------------------------------------- 
>    call write_(nch,ioutput,iserror)
>    if (iserror) then
>      msg='Error when calling write_'
>      if (haveioutput) then
>       write (ioutput,*) msg
>      end if
>      goto 20
>    end if
> end if
> 
!%--------------------------------------------------------------------
!%Remove mineral phases according crystalline fracctionation factor
!%--------------------------------------------------------------------
if (nminsp>0) then
>  molminstep(:,i) = fraccryst * nch%c(idminsp) * nch%omgwfree
>  nch%c(idminsp) = (r1-fraccryst) * nch%c(idminsp)
end if
> 
> 
end do
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!% Determine final cputime and write cpu time consumed
!%--------------------------------------------------------------------
if (haveioutput) then
> call cpu_time (icputime2)
> dtcputime=icputime2-icputime1
> write(ioutput,*) '----------------------------------'
> write(ioutput,1) 'CPU time consumed:',dtcputime,'[s]'
> write(ioutput,*) '----------------------------------'
end if
!%--------------------------------------------------------------------
!% Write information about relative error
!%--------------------------------------------------------------------
if (haveioutput) then
> 
> call get_chem_info_ (nch,iserror,molaq=molaq,molads=molads,molprec=molprec)
> 
> head='Relative error in components'
> 
> molfin = molaq + molads + molprec
> 
> errorcomp = dabs((molinit - molfin) / molfin)
> 
> writehead=.true.
> 
> call write_ (ioutput,head,namesp(idpri),errorcomp,npri,'error',writehead,writehead)
> 
> 
> 
end if
!%--------------------------------------------------------------------
!% Update the last information in nodal chemistry object 
!%--------------------------------------------------------------------
this=nch
!%--------------------------------------------------------------------
!% Destroy and deallocate nodal chemistry object 
!%--------------------------------------------------------------------
call destroy_ (nch)
deallocate (nch)
!%--------------------------------------------------------------------
!% Write mol of minerals vs. evaporation/dilution progress
!%--------------------------------------------------------------------
if (haveioutput.and.nminsp>0) then
> 
> head='Mol of minerals vs. time'

> 
> do i=1,istep
> 
>     if (i==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namesp(idminsp), &
>       molminstep(:,i), &
>       nminsp, &
>       time(i), &
>       writehead, &
>       writehead)
> 
> 
> end do
!%--------------------------------------------------------------------
!% Write m3 of minerals vs. evaporation/dilution progress
!%--------------------------------------------------------------------
> 
> head='Volume of minerals (m3) vs. time'
> 
> 
> do i=1,istep
> 
>     do j=1,nminsp
>        call change_chem_unit_ (this%pchemsys,molminstep(j,i),namesp(idminsp(j)),'mol','m3',iserror)
>     end do
> 
>     if (i==1) then
>        writehead=.true.
>     else
>        writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namesp(idminsp), &
>       molminstep(:,i), &
>       nminsp, &
>       time(i), &
>       writehead, &
>       writehead)
> 
> 
> end do
!%--------------------------------------------------------------------
!% Write saturation indices of minerals vs. evaporation/dilution 
!% progress
!%--------------------------------------------------------------------
> 
> head='Evolution of saturation indices vs. time'
> 
> 
> simin(1:nminsp,1:istep)=dlog10(simin(1:nminsp,1:istep))
> 
> do i=1,istep 
> 
>     if (i==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namesp(idminsp), &
>       simin(:,i), &
>       nminsp, &
>       time(i), &
>       writehead, &
>       writehead)
> 
> 
> end do
> 
> 
end if
!%--------------------------------------------------------------------
!% Write concentration of species vs. dilution/evaporation progress
!%--------------------------------------------------------------------
if (haveioutput.and.nspout>0) then

> head='Evolution of concentration species vs. time'

> do i=1,istep
> 
>     if (i==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       cout(1:nspout,i), &
>       nspout, &
>       time(i), &
>       writehead, &
>       writehead)
> end do

> head='Evolution of mol of species vs. time'
> 
> 
> cout=cout*this%omgwfree 

> do i=1,istep
> 
>     if (i==1) then
>      writehead=.true.
>     else
>      writehead=.false.
>     end if
> 
>     cout(1:nspout,i) = cout(1:nspout,i) * masswfree(i)

>     call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       cout(1:nspout,i), &
>       nspout, &
>       time(i), &
>       writehead, &
>       writehead)
> end do
> 

> head='Evolution of activity of species vs. time'

> do i=1,istep
> 
>     if (i==1) then
>       writehead=.true.
>     else
>       writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       namespout, &
>       actout(1:nspout,i), &
>       nspout, &
>       time(i), &
>       writehead, &
>       writehead)
> end do
> 
> 
> head='Evolution of chemical parameters vs. time'
> 
> do i=1,istep
> 
>     if (i==1) then
>       writehead=.true.
>     else
>       writehead=.false.
>     end if
> 
>     call write_ &
>      (ioutput, &
>       head, &
>       nameparam, &
>       chemparam(1:nparam,i), &
>       nparam, &
>       time(i), &
>       writehead, &
>       writehead)
> end do
> 
> 
> 
end if
20 continue 
!%--------------------------------------------------------------------
!% Deallocate local pointers 
!%--------------------------------------------------------------------
call check_pointer_ (masswfree,1,.false.)
call check_pointer_ (idminsp,1,.false.)
call check_pointer_ (simin,1,1,.false.)
call check_pointer_ (siloc,1,.false.)
call check_pointer_ (molminstep,1,1,.false.)
call check_pointer_ (idspout,1,.false.)
call check_pointer_ (cout,1,1,.false.)
call check_pointer_ (actout,1,1,.false.)
call check_pointer_ (errorcomp,1,.false.)
call check_pointer_ (molinit,1,.false.)
call check_pointer_ (molfin,1,.false.)
call check_pointer_ (molaq,1,.false.)
call check_pointer_ (molprec,1,.false.)
call check_pointer_ (molads,1,.false.)
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (chemparam,1,1,.false.)
call check_pointer_ (nameparam,1,.false.)
call check_pointer_ (progress,1,.false.)
call check_pointer_ (tk,1,.false.)
call check_pointer_ (tk1,1,.false.)
call check_pointer_ (idpri,1,.false.)
if (iserror) goto 10
if (haveisconvergence.and..not.isconv) then
>  isconvergence=.false.
>  return
else if (haveisconvergence.and.isconv) then
>  isconvergence=.true.
>  return
end if
if (.not.isconv) stop
!%--------------------------------------------------------------------
return
> 
10 continue 
print *,'************************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: reaction_path_'
print *, msg
print *,'************************'
iserror=.true.
return
> 
1 format(a21,f10.5,a6)
2 format(a34,f10.5,a20)
> 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_density_ith_nch &
>    (this       , &
>    phasecode   , &
>    density     , &
>    iserror     , &
>    ddensitydc)
>    

type(t_nodalchemistry),intent(in)       :: this         !Type nodal chemistry variable. 

integer,intent(in)                      :: phasecode    !Phase for which the density will be calculated

real*8,intent(out)                      :: density      !Density     

logical,intent(out)                     :: IsError

real*8,pointer,dimension(:),optional    :: ddensitydc !Density derivative

!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
> 
!-------------------------------------------------------------------------
!
!>   $code
!

if (present(ddensitydc)) then
>    call compute_density_(this%pchemsys,phasecode,this%pliq,this%temp,this%c,density,iserror,this%dc,ddensitydc)
>    if (isError) goto 10
else
>    call compute_density_(this%pchemsys,phasecode,this%pliq,this%temp,this%c,density,iserror)
>    if (isError) goto 10
endif


!%----------------------------------------------------------
return
> 
10 continue 
print *,'***********************'
print *,'Nodal Chemistry:'
print *,'Name:',this%name
print *,'Service: compute_density_ith_'
print *,'***********************'
iserror=.true.
return
> 
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_nodalchemistry
