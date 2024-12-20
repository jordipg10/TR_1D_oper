!> Chemistry class (main class): contains all chemical information and solves reactive mixing
module chemistry_Lagr_m
    use matrices_m
    use aqueous_chemistry_m
    use CV_params_m
    use chem_out_options_m
    use spatial_discr_1D_m
    implicit none
    save
    type, public :: chemistry_c
        integer(kind=4) :: option !> option for reading chemical data (1: CHEPROO-based, 2: PHREEQC, 3: PFLOTRAN)
        integer(kind=4) :: act_coeffs_model !> model to compute activity coefficients
        type(chem_system_c) :: chem_syst !> chemical system object
        integer(kind=4) :: num_target_waters=0 !> number of target waters
        type(aqueous_chemistry_c), allocatable :: target_waters(:) !> target waters
        type(aqueous_chemistry_c), allocatable :: target_waters_init(:) !> target waters initial
        integer(kind=4) :: num_target_waters_dom=0 !> number of initial target waters
        integer(kind=4) :: num_ext_waters=0 !> number of external waters
        integer(kind=4), allocatable :: ext_waters_indices(:) !> external waters indices
        integer(kind=4), allocatable :: dom_tar_wat_indices(:) !> domain target waters indices
        integer(kind=4) :: num_target_solids=0 !> number of target solids (<= num_target_waters)
        integer(kind=4) :: num_target_solids_dom=0 !> number of target solids init (<= num_target_waters_dom)
        type(solid_chemistry_c), allocatable :: target_solids(:) !> target solids
        type(solid_chemistry_c), allocatable :: target_solids_init(:) !> initial target solids
        integer(kind=4) :: num_target_gases=0 !> number of target gases
        type(gas_chemistry_c), allocatable :: target_gases(:) !> target gases
        type(gas_chemistry_c), allocatable :: target_gases_init(:) !> target gases init
        integer(kind=4) :: num_reactive_zones=0 !> number of reactive zones (<=num_target_solids)
        type(reactive_zone_c), allocatable :: reactive_zones(:) !> reactive zones
        integer(kind=4) :: Jac_flag !> model to compute Jacobians (0: incremental coefficnets, 1: analytical)
        type(CV_params_s) :: CV_params !> parameters for convergence
        type(chem_out_options_t) :: chem_out_options !> output results options variable
    contains
    !> Set
        procedure, public :: set_chem_syst
        procedure, public :: set_num_target_waters_dom
        procedure, public :: set_num_target_solids
        procedure, public :: set_num_ext_waters
        procedure, public :: set_target_waters
        procedure, public :: set_target_solids
        procedure, public :: set_target_gases
        procedure, public :: set_option
        procedure, public :: set_reactive_zones
        procedure, public :: set_Jac_flag
    !> Allocate
        procedure, public :: allocate_target_waters
        procedure, public :: allocate_dom_tar_wat_indices
        procedure, public :: allocate_ext_waters_indices
        procedure, public :: allocate_target_solids
        procedure, public :: allocate_target_gases
        procedure, public :: allocate_reactive_zones
    !> Read
        procedure, public :: read_target_waters_init
        procedure, public :: read_init_min_zones_CHEPROO
        procedure, public :: read_chemistry
        procedure, public :: read_chemistry_CHEPROO
        !procedure, public :: read_chemistry_PHREEQC
        procedure, public :: read_init_bd_rech_wat_types_CHEPROO
        procedure, public :: read_init_cat_exch_zones_CHEPROO
        procedure, public :: read_gas_bd_zones_CHEPROO
        procedure, public :: read_init_gas_zones_CHEPROO
    !> Initialisation
        procedure, public :: initialise_chemistry
    !> Write
        procedure, public :: write_chemistry
    !> Solve
        procedure, public :: solve_reactive_mixing !> main solver
        procedure, public :: solve_reactive_mixing_ideal !> main solver
        procedure, public :: solve_reactive_mixing_bis !> main solver
        procedure, public :: solve_reactive_mixing_BCs_dep_t !> main solver
    !> Link
        procedure, public :: link_target_waters_target_solids
        procedure, public :: link_target_waters_target_gases
        procedure, public :: link_target_solids_reactive_zone
        procedure, public :: link_target_gases_reactive_zone
        procedure, public :: link_target_waters_reactive_zone
    !> Check
        procedure, public :: check_new_reactive_zones
    end type
    
   
    
    interface
        
        subroutine solve_reactive_mixing_bis(this,root,unit,mixing_ratios,mixing_waters_indices,F_mat,time_discr,int_method_chem_reacts)
            import chemistry_c
            import real_array_c
            import int_array_c
            import diag_matrix_c
            import time_discr_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            integer(kind=4), intent(in) :: unit
            class(real_array_c), intent(in) :: mixing_ratios
            class(int_array_c), intent(in) :: mixing_waters_indices
            !class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
            real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
            class(time_discr_c), intent(in) :: time_discr !> time discretisation object
            integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
        end subroutine
        
        subroutine solve_reactive_mixing(this,root,unit,mixing_ratios,mixing_waters_indices,F_mat,time_discr,int_method_chem_reacts)
            import chemistry_c
            import real_array_c
            import int_array_c
            import diag_matrix_c
            import time_discr_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            integer(kind=4), intent(in) :: unit
            class(real_array_c), intent(in) :: mixing_ratios
            class(int_array_c), intent(in) :: mixing_waters_indices
            !class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
            real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
            class(time_discr_c), intent(in) :: time_discr !> time discretisation object
            integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
        end subroutine
        
        subroutine solve_reactive_mixing_ideal(this,root,unit,mixing_ratios,mixing_waters_indices,F_mat,time_discr,int_method_chem_reacts)
            import chemistry_c
            import real_array_c
            import int_array_c
            import diag_matrix_c
            import time_discr_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            integer(kind=4), intent(in) :: unit
            class(real_array_c), intent(in) :: mixing_ratios
            class(int_array_c), intent(in) :: mixing_waters_indices
            !class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
            real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
            class(time_discr_c), intent(in) :: time_discr !> time discretisation object
            integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
        end subroutine
        
        
        subroutine solve_reactive_mixing_BCs_dep_t(this,root,unit,mixing_ratios,mixing_waters_indices,time_discr_tpt,int_method_chem_reacts,spatial_discr_tpt,D,q,phi,anal_sol)
            import chemistry_c
            import spatial_discr_c
            import time_discr_c
            import real_array_c
            import int_array_c
            implicit none
        !> Arguments
            class(chemistry_c) :: this !> chemistry object
            character(len=*), intent(in) :: root
            integer(kind=4), intent(in) :: unit
            class(real_array_c), intent(in) :: mixing_ratios !> mixing ratios matrix
            class(int_array_c), intent(in) :: mixing_waters_indices !> matrix that contains indices of target waters that mix with each target water
            class(time_discr_c), intent(in) :: time_discr_tpt !> time discretisation object (used to solve transport)
            integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
            class(spatial_discr_c), intent(in) :: spatial_discr_tpt !> spatial discretisation object (used to solve transport)
            real(kind=8), intent(in) :: D
            real(kind=8), intent(in) :: q
            real(kind=8), intent(in) :: phi
            real(kind=8), external :: anal_sol
        end subroutine
        
        subroutine read_chemistry_PHREEQC(this,path_inp,path_DB,filename)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: path_inp
            character(len=*), intent(in) :: path_DB
            character(len=*), intent(in) :: filename
        end subroutine
        
        subroutine read_chemistry_CHEPROO(this,root,path_DB,unit_chem_syst_file,unit_loc_chem_file,unit_target_waters_dom_file,unit_output_file)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            character(len=*), intent(in) :: path_DB
            integer(kind=4), intent(in) :: unit_chem_syst_file
            !character(len=*), intent(in) :: chem_syst_file
            integer(kind=4), intent(in) :: unit_loc_chem_file
            !character(len=*), intent(in) :: loc_chem_file
            integer(kind=4), intent(in) :: unit_target_waters_dom_file
            !character(len=*), intent(in) :: target_waters_dom_file
            integer(kind=4), intent(in) :: unit_output_file
            !character(len=*), intent(in) :: output_file
        end subroutine
        
        subroutine read_chemistry(this,root,path_DB,unit_chem_syst_file,unit_loc_chem_file,unit_target_waters_dom_file,unit_output_file)
            import chemistry_c
            implicit none    
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            character(len=*), intent(in) :: path_DB
            integer(kind=4), intent(in) :: unit_chem_syst_file
            !character(len=*), intent(in) :: chem_syst_file
            integer(kind=4), intent(in) :: unit_loc_chem_file
            !character(len=*), intent(in) :: loc_chem_file
            integer(kind=4), intent(in) :: unit_target_waters_dom_file
            !character(len=*), intent(in) :: target_waters_dom_file
            integer(kind=4), intent(in) :: unit_output_file
            !character(len=*), intent(in) :: output_file
        end subroutine
                
        subroutine read_reactive_zones_Lagr(this,unit)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit
            !character(len=*), intent(in) :: filename
            !integer(kind=4), intent(in) :: line
        end subroutine
        
        subroutine link_target_waters_dom_reactive_zone(this,i,tar_wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: tar_wat_indices(:) 
        end subroutine
        
        subroutine link_target_waters_reactive_zone(this,i,dom_indices,ext_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: dom_indices(:) 
            integer(kind=4), intent(out), allocatable :: ext_indices(:) 
        end subroutine
        
        subroutine link_target_waters_target_solid(this,i,tw_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: tw_indices(:) 
        end subroutine
        
        subroutine link_target_waters_target_solids(this,tar_sol_indices,tar_wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: tar_sol_indices(:) !> target solid indices
            integer(kind=4), intent(inout), allocatable :: tar_wat_indices(:) 
        end subroutine
        
        subroutine link_target_waters_target_gases(this,tar_gas_indices,tar_wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: tar_gas_indices(:) !> target gas indices
            integer(kind=4), intent(inout), allocatable :: tar_wat_indices(:) 
        end subroutine
        
        subroutine link_target_solids_reactive_zone(this,i,tar_sol_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: tar_sol_indices(:)
        end subroutine
        
        subroutine link_target_gases_reactive_zone(this,i,tar_gas_indices) 
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: tar_gas_indices(:)
        end subroutine
        
       

        
        subroutine initialise_target_waters_dom(this,initial_water_types)
            import chemistry_c
            import aqueous_chemistry_c
            import chem_system_c
            implicit none
            class(chemistry_c) :: this
            class(aqueous_chemistry_c), intent(in) :: initial_water_types(:)
        end subroutine
        
        subroutine read_target_waters_init(this,unit,water_types,init_sol_types,init_gas_types,niter,CV_flag)
            import chemistry_c
            import aqueous_chemistry_c
            import solid_chemistry_c
            import gas_chemistry_c
            import reactive_zone_c
            import aq_phase_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            !class(aqueous_chemistry_c), intent(in) :: init_water_types(:)
            !class(aqueous_chemistry_c), intent(in) :: bd_water_types(:)
            type(aqueous_chemistry_c), intent(in) :: water_types(:)
            type(solid_chemistry_c), intent(in) :: init_sol_types(:)
            type(gas_chemistry_c), intent(in) :: init_gas_types(:)
            !type(reactive_zone_c), intent(in) :: react_zones(:) !> all possible reactive zones
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
            !class(aq_phase_c), intent(out), optional :: aq_phase_new
        end subroutine
        
        subroutine read_target_waters_dom_bis(this,unit,init_water_types,bd_water_types,init_sol_types,niter,CV_flag)
            import chemistry_c
            import aqueous_chemistry_c
            import solid_chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            class(aqueous_chemistry_c), intent(in) :: init_water_types(:)
            class(aqueous_chemistry_c), intent(in) :: bd_water_types(:)
            class(solid_chemistry_c), intent(in) :: init_sol_types(:)
            !real(kind=8), intent(in) :: tolerance
            !real(kind=8), intent(in) :: rel_tolerance
            !real(kind=8), intent(in) :: control_factor
            !integer(kind=4), intent(in) :: niter_max
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
       
        
        subroutine initialise_target_solids(this,n)
            import chemistry_c
            import solid_chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: n
        end subroutine
        
        subroutine initialise_ext_waters(this)
            import chemistry_c
            import chem_system_c
            implicit none
            class(chemistry_c) :: this
        end subroutine
        
       
        
        subroutine initialise_chemistry(this,path_DB,root,unit_chem_syst_file,unit_loc_chem_file,unit_target_waters_dom_file,unit_output_file)
            import chemistry_c
            import real_array_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            character(len=*), intent(in) :: path_DB
            integer(kind=4), intent(in) :: unit_chem_syst_file
            !character(len=*), intent(in) :: chem_syst_file
            integer(kind=4), intent(in) :: unit_loc_chem_file
            !character(len=*), intent(in) :: loc_chem_file
            integer(kind=4), intent(in) :: unit_target_waters_dom_file
            !character(len=*), intent(in) :: target_waters_dom_file
            integer(kind=4), intent(in) :: unit_output_file
            !character(len=*), intent(in) :: output_file
        end subroutine
        
       
        
        subroutine write_chemistry(this,unit)
            import chemistry_c
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit
            !character(len=*), intent(in) :: file_out
        end subroutine
        
        subroutine check_new_reactive_zones(this,i,tolerance)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in)  :: i !> reactive zone index
            real(kind=8), intent(in) :: tolerance !> for concentration of non-flowing species
        end subroutine
        
       
        
        subroutine read_init_bd_rech_wat_types_CHEPROO(this,unit,ind_wat_type,num_aq_prim_array,num_cstr_array,init_cat_exch_zones,wat_types,gas_chem)
            import chemistry_c
            import aqueous_chemistry_c
            import solid_chemistry_c
            import gas_chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            !type(aqueous_chemistry_c), intent(out), allocatable :: init_wat_types(:)
            !type(aqueous_chemistry_c), intent(out), allocatable :: bd_wat_types(:)
            !type(aqueous_chemistry_c), intent(out), allocatable :: rech_wat_types(:)
            !real(kind=8), intent(in) :: tolerance
            !real(kind=8), intent(in) :: rel_tolerance
            !real(kind=8), intent(in) :: control_factor
            !integer(kind=4), intent(in) :: niter_max
            integer(kind=4), intent(out), allocatable :: ind_wat_type(:)
            integer(kind=4), intent(out), allocatable :: num_aq_prim_array(:)
            integer(kind=4), intent(out), allocatable :: num_cstr_array(:)
            type(solid_chemistry_c), intent(inout) :: init_cat_exch_zones(:)
            type(aqueous_chemistry_c), intent(out), allocatable :: wat_types(:)
            type(gas_chemistry_c), intent(in),optional :: gas_chem !> chapuza
            !logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise 
        end subroutine
        
        subroutine read_init_min_zones_CHEPROO(this,unit,init_min_zones,reactive_zones)
            import chemistry_c
            import solid_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            type(solid_chemistry_c), intent(out), allocatable :: init_min_zones(:)
            type(reactive_zone_c), intent(out), allocatable, optional :: reactive_zones(:)
        end subroutine
        
        subroutine read_init_min_zones_CHEPROO_bis(this,unit,init_min_zones,reactive_zones)
            import chemistry_c
            import solid_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            type(solid_chemistry_c), intent(out), allocatable :: init_min_zones(:)
            type(reactive_zone_c), intent(out), allocatable, optional :: reactive_zones(:)
        end subroutine
        
        subroutine read_init_cat_exch_zones_CHEPROO(this,unit,init_cat_exch_zones,reactive_zones)
            import chemistry_c
            import solid_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            type(solid_chemistry_c), intent(out), allocatable :: init_cat_exch_zones(:)
            type(reactive_zone_c), intent(inout), allocatable, optional :: reactive_zones(:)
        end subroutine
        
        subroutine read_gas_bd_zones_CHEPROO(this,unit,gas_bd_zones,reactive_zones)
            import chemistry_c
            import gas_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            type(gas_chemistry_c), intent(out), allocatable :: gas_bd_zones(:)
            type(reactive_zone_c), intent(out), allocatable, optional :: reactive_zones(:)
        end subroutine
        
        subroutine read_init_gas_zones_CHEPROO(this,unit,gas_zones,reactive_zones)
            import chemistry_c
            import gas_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            type(gas_chemistry_c), intent(out), allocatable :: gas_zones(:)
            type(reactive_zone_c), intent(out), allocatable, optional :: reactive_zones(:)
        end subroutine
    end interface
    
    contains
        subroutine set_option(this,option)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: option
            if (option<1 .or. option>3) then
                error stop "Chemistry option not implemented yet"
            else if (option>1) then
                error stop "Chemistry option not fully implemented yet"
            else
                this%option=option
            end if
        end subroutine
        
        subroutine set_Jac_flag(this,Jac_flag)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: Jac_flag
            if (Jac_flag<0 .or. Jac_flag>1) then
                error stop "Chemistry attribute 'Jac_flag' must be 0 or 1"
            else
                this%Jac_flag=Jac_flag
            end if
        end subroutine
        
        subroutine set_chem_syst(this,chem_syst_obj)
            implicit none
            class(chemistry_c) :: this
            type(chem_system_c), intent(in) :: chem_syst_obj
            this%chem_syst=chem_syst_obj
        end subroutine
        
       
        
        subroutine set_num_target_waters_dom(this,num_target_waters_dom)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_target_waters_dom
            if (present(num_target_waters_dom)) then
                this%num_target_waters_dom=num_target_waters_dom
            else
                this%num_target_waters_dom=this%num_target_waters-this%num_ext_waters
            end if
        end subroutine
        
        subroutine set_num_target_solids(this,num_target_solids)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: num_target_solids
            this%num_target_solids=num_target_solids
        end subroutine
        
        subroutine set_num_ext_waters(this,num_ext_waters)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_ext_waters
            if (present(num_ext_waters)) then
                this%num_ext_waters=num_ext_waters
            else
                this%num_ext_waters=this%num_target_waters_dom
            end if
        end subroutine
        
        !subroutine set_ext_waters(this,ext_waters)
        !    implicit none
        !    class(chemistry_c) :: this
        !    class(aqueous_chemistry_c), intent(in) :: ext_waters(:)
        !    this%ext_waters=ext_waters
        !end subroutine
        
        subroutine set_target_solids(this,target_solids)
            implicit none
            class(chemistry_c) :: this
            type(solid_chemistry_c), intent(in) :: target_solids(:)
            this%target_solids=target_solids
            this%num_target_solids=size(target_solids)
        end subroutine
        
        subroutine set_target_gases(this,target_gases)
            implicit none
            class(chemistry_c) :: this
            type(gas_chemistry_c), intent(in) :: target_gases(:)
            this%target_gases=target_gases
            this%num_target_gases=size(target_gases)
        end subroutine
        
       subroutine allocate_target_waters(this,num_tar_wat)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_tar_wat
            if (present(num_tar_wat)) then
                this%num_target_waters=num_tar_wat
            end if
            allocate(this%target_waters(this%num_target_waters),this%target_waters_init(this%num_target_waters))
        end subroutine
        
        subroutine allocate_dom_tar_wat_indices(this,num_tar_wat_dom)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_tar_wat_dom
            if (present(num_tar_wat_dom)) then
                this%num_target_waters_dom=num_tar_wat_dom
            end if
            allocate(this%dom_tar_wat_indices(this%num_target_waters_dom))
        end subroutine
        
        subroutine allocate_target_solids(this,n)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: n
            if (present(n)) then
                this%num_target_solids=n
            end if
            allocate(this%target_solids(this%num_target_solids),this%target_solids_init(this%num_target_solids))
        end subroutine
        
        
        subroutine allocate_target_gases(this,n)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: n
            if (present(n)) then
                this%num_target_gases=n
            end if
            allocate(this%target_gases(this%num_target_gases),this%target_gases_init(this%num_target_gases))
        end subroutine
        
        subroutine allocate_ext_waters_indices(this,num_ext_waters)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_ext_waters
            if (present(num_ext_waters)) then
                this%num_ext_waters=num_ext_waters
            end if
            allocate(this%ext_waters_indices(this%num_ext_waters))
        end subroutine
        
       
        
        subroutine set_target_waters(this,target_waters)
            implicit none
            class(chemistry_c) :: this
            type(aqueous_chemistry_c), intent(in) :: target_waters(:)
            this%target_waters=target_waters
            this%num_target_waters=size(target_waters)
        end subroutine
        
        !subroutine set_target_waters_dom(this,target_waters_dom)
        !    implicit none
        !    class(chemistry_c) :: this
        !    class(aqueous_chemistry_c), intent(in) :: target_waters_dom(:)
        !    this%target_waters()=target_waters_dom
        !    this%num_target_waters_dom=size(target_waters_dom)
        !end subroutine
        
       subroutine allocate_reactive_zones(this,num_reactive_zones)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_reactive_zones
            if (present(num_reactive_zones)) then
                this%num_reactive_zones=num_reactive_zones
            end if
            if (allocated(this%reactive_zones)) then
                deallocate(this%reactive_zones)
            end if
            allocate(this%reactive_zones(this%num_reactive_zones))
       end subroutine
       
       subroutine set_reactive_zones(this,reactive_zones)
            implicit none
            class(chemistry_c) :: this
            type(reactive_zone_c), intent(in), optional :: reactive_zones(:)
            integer(kind=4) :: i,j,l
            logical :: flag
            integer(kind=4), allocatable :: rz_indices(:)
            if (present(reactive_zones)) then
                this%reactive_zones=reactive_zones
                this%num_reactive_zones=size(reactive_zones)
            else if (this%num_target_solids>0) then
                i=1
                j=2
                this%num_reactive_zones=1
                allocate(rz_indices(this%num_target_solids))
                rz_indices=0
                do
                        call compare_react_zones(this%target_solids(i)%reactive_zone,this%target_solids(j)%reactive_zone,flag)
                        if (flag==.true.) then
                            if (i<this%num_target_solids-1) then
                                i=i+1
                                j=i+1
                            else
                                exit
                            end if
                        else if (j<this%num_target_solids) then
                            j=j+1
                        else if (i<this%num_target_solids-1) then
                            this%num_reactive_zones=this%num_reactive_zones+1
                            rz_indices(i)=1
                            i=i+1
                            j=i+1
                        else
                            exit
                        end if
                end do
                rz_indices(j)=1
                call this%allocate_reactive_zones()
                l=1
                do i=1,this%num_target_solids
                    if (rz_indices(i)==1) then
                        call this%reactive_zones(l)%assign_react_zone(this%target_solids(i)%reactive_zone)
                        l=l+1
                    end if
                end do
            else if (this%num_target_gases>0) then
                i=1
                j=2
                this%num_reactive_zones=1
                allocate(rz_indices(this%num_target_gases))
                rz_indices=0
                do
                        call compare_react_zones(this%target_gases(i)%reactive_zone,this%target_gases(j)%reactive_zone,flag)
                        if (flag==.true.) then
                            if (i<this%num_target_gases-1) then
                                i=i+1
                                j=i+1
                            else
                                exit
                            end if
                        else if (j<this%num_target_gases) then
                            j=j+1
                        else if (i<this%num_target_gases-1) then
                            this%num_reactive_zones=this%num_reactive_zones+1
                            rz_indices(i)=1
                            i=i+1
                            j=i+1
                        else
                            exit
                        end if
                end do
                rz_indices(j)=1
                call this%allocate_reactive_zones()
                l=1
                do i=1,this%num_target_gases
                    if (rz_indices(i)==1) then
                        this%reactive_zones(l)=this%target_gases(i)%reactive_zone
                        l=l+1
                    end if
                end do
            else
                call this%allocate_reactive_zones(0)
            end if
       end subroutine
end module