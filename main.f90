program main
    use RT_1D_m
    implicit none 
!> Objects
    type(RT_1D_transient_c) :: my_RT_trans !> 1D transient reactive transport class
    type(transport_1D_transient_c) :: my_tpt_trans !> 1D transient transport class
    type(chemistry_c) :: my_chem !> chemistry class
!> Variables
    integer(kind=4) :: problem
    integer(kind=4) :: option_chem !> chemical input option (1: CHEPROO-based)
    integer(kind=4) :: opc_Jac !> Jacobian model (0: incremental coefficient, 1: analytically)
    integer(kind=4) :: option_tpt !> transport option (0: compute lambdas, 1: read lambdas)
    integer(kind=4) :: unit_tpt !> transport file unit
    integer(kind=4) :: unit_out !> output file unit
    integer(kind=4) :: unit_res !> results file unit
    integer(kind=4) :: unit_discr !> temporal discretisation file unit
    integer(kind=4) :: int_method_chem !> integration method for chemical reactions (1: Euler explicit, 2: Euler fully implicit)
    integer(kind=4) :: unit_chem_syst !> chemical system file unit
    integer(kind=4) :: unit_loc_chem !> local chemistry file unit
    integer(kind=4) :: unit_tw !> target waters file unit
    character(len=256) :: path_inp !> path for reading files
    character(len=:), allocatable :: path_inp_trim !> path for reading files trimmed
    character(len=256) :: path_DB !> path for database
    character(len=:), allocatable :: path_DB_trim !> path for database trimmed
    character(len=256) :: root !> root of problem to solve
    character(len=:), allocatable :: root_trim !> root of problem to solve trimmed
!**************************************************************************************************
!> Name of path containing chemical and transport information
    path_DB = 'C:\Users\Jordi\source\repos\jordipg10\TR_1D_oper\BBDD\' !> must be written by the user
    path_DB_trim = trim(path_DB)
!> File units (arbitrary)
    unit_chem_syst=1
    unit_loc_chem=2
    unit_tw=3
    unit_discr=10
    unit_tpt=11
    unit_out=4
    unit_res=6
!> Choose problem
    problem=1
    if (problem==1) then
        root='C:\Users\Jordi\source\repos\jordipg10\TR_1D_oper\examples\gypsum_eq\gypsum_eq' !> name of path containing user input (must be written by the user)
    else if (problem==2) then
        root='C:\Users\Jordi\source\repos\jordipg10\TR_1D_oper\examples\cc_anh_eq\cc_anh' !> name of path containing user input (must be written by the user)
    else if (problem==3) then
        root='C:\Users\Jordi\source\repos\jordipg10\TR_1D_oper\examples\cc_anh_kin\cc_anh' !> name of path containing user input (must be written by the user)
    else
        error stop "Problem not implemented yet"
    end if
    root_trim=trim(root)
!> Initialise transport
    option_tpt=1
    if (option_tpt==0) then
    !> we read transport data, BCs and discretisations
        !> in the explicit case, we also compute stability parameters
            call my_tpt_trans%initialise_transport_1D_transient_RT(root_trim)
        !> we allocate transport arrays
            call my_tpt_trans%allocate_arrays_PDE_1D()
        !> we compute transport arrays, including mixing ratios, and we impose BCs
            call my_tpt_trans%compute_mixing_ratios_Delta_t_homog()
        !> we set transport attribute in reactive transport object
            call my_RT_trans%set_transport_trans(my_tpt_trans)
        !> we choose and set integration method for chemical reactions
            int_method_chem=2
            call my_RT_trans%set_int_method_chem_reacts(int_method_chem)
    else if (option_tpt==1) then
    !> we read transport data for WMA
        call my_tpt_trans%read_transport_data_WMA(unit_tpt,root_trim)
    !> we set transport attribute in reactive transport object
        call my_RT_trans%set_transport(my_tpt_trans)
    !> we read temporal discretisation
        call my_RT_trans%read_time_discretisation(unit_discr,root_trim)
    else
        error stop "This option is not implemented yet"
    end if
!> Initialise chemistry
    option_chem=1 
    call my_chem%set_option(option_chem)
    opc_Jac=1
    call my_chem%set_Jac_flag(opc_Jac)
    !> we read chemistry
    call my_chem%read_chemistry(root_trim,path_DB_trim,unit_chem_syst,unit_loc_chem,unit_tw,unit_out)
!> We call the main solver
    call my_chem%solve_reactive_mixing(root_trim,unit_res,my_RT_trans%transport%mixing_ratios,my_RT_trans%transport%mixing_waters_indices,my_RT_trans%transport%F_mat%diag,my_RT_trans%transport%time_discr,my_RT_trans%int_method_chem_reacts)
!> We set chemistry attribute in reactive transport object
    call my_RT_trans%set_chemistry(my_chem)
!> We write data and results
    call my_RT_trans%write_RT_1D(unit_res,root_trim)
end program