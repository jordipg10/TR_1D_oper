program main
    use RT_1D_m, only: RT_1D_transient_c, transport_1D_transient_c, chemistry_c
    implicit none 
!> Objects
    type(RT_1D_transient_c) :: my_RT_trans !> 1D transient reactive transport class
    type(transport_1D_transient_c) :: my_tpt_trans !> 1D transient transport class
    type(chemistry_c) :: my_chem !> chemistry class
!> Variables
    integer(kind=4) :: problem !> problem to solve
    integer(kind=4) :: cons_opt !> lumping option
    integer(kind=4) :: rk_opt !> rk average option
    integer(kind=4) :: option_chem !> chemical input option (1: CHEPROO-based)
    integer(kind=4) :: opc_Jac !> Jacobian model (0: incremental coefficient, 1: analytically)
    integer(kind=4) :: option_tpt !> transport option (0: compute lambdas, 1: read lambdas)
    integer(kind=4) :: int_method_chem !> integration method for chemical reactions (1: Euler explicit, 2: Euler fully implicit)
    character(len=256) :: path_py !> path for Python
    character(len=256) :: path_inp !> path for reading files
    character(len=:), allocatable :: path_inp_trim !> path for reading files trimmed
    character(len=256) :: path_DB !> path for database
    character(len=:), allocatable :: path_DB_trim !> path for database trimmed
    character(len=256) :: root !> root of problem to solve
    character(len=:), allocatable :: root_trim !> root of problem to solve trimmed
!****************************************************************************************************************************
!> Name of path containing chemical and transport information
    write(*,*) "Path to databases?"
    read(*,*) path_DB !> must be written by the user
    path_DB_trim = trim(path_DB)
!> Choose problem
    write(*,*) "Path and root of problem to be solved?"
    read(*,*) root !> must be written by the user
    root_trim=trim(root)
!> Initialise transport
    write(*,*) "Choose transport option (0: compute lambdas, 1: read lambdas):"
    read(*,*) option_tpt !> must be written by the user
    if (option_tpt.eq.0) then !> compute lambdas
    !> we read transport data, BCs and discretisations
        !> in the explicit case, we also compute stability parameters
            call my_tpt_trans%initialise_transport_1D_transient_RT(root_trim)
        !> we allocate transport arrays
            call my_tpt_trans%allocate_arrays_PDE_1D()
            call my_tpt_trans%allocate_conc()
        !> we compute transport arrays, including mixing ratios, and we impose BCs
            call my_tpt_trans%compute_mixing_ratios_Delta_t_homog() !> missing the case of heterogenous time steps
        !> we set transport attribute in reactive transport object
            call my_RT_trans%set_transport_trans(my_tpt_trans)
        !> we choose and set integration method for chemical reactions
            !! 1: Euler explicit, 2: Euler fully implicit, 3: Crank-Nicolson
            write(*,*) "Choose integration method for chemical reactions (1: Euler explicit, 2: Euler fully implicit, 3:Crank-Nicolson):"
            read(*,*) int_method_chem !> must be written by the user
            call my_RT_trans%set_int_method_chem_reacts(int_method_chem)
    else if (option_tpt.eq.1) then !> read lambdas
    !> we set transport attribute in reactive transport object
        call my_RT_trans%set_transport_trans(my_tpt_trans) !> esto es un create en realidad
    !> we read temporal discretisation
        call my_RT_trans%read_time_discretisation(root_trim)
    !> we read transport data for WMA
        call my_RT_trans%transport%read_transport_data_WMA(root_trim)
    else
        error stop "This option is not implemented yet"
    end if
!> we read chemistry
    call my_chem%read_chemistry(root_trim,path_DB_trim)
!> We call the main solver
    if (my_chem%act_coeffs_model==0 .and. my_chem%lump_flag .eqv. .true.) then !> ideal with lumping
        call my_chem%solve_reactive_mixing_ideal_lump(root_trim,my_RT_trans%transport%mixing_ratios_conc,&
             my_RT_trans%transport%mixing_waters_indices,my_RT_trans%transport%mixing_waters_indices_dom,&
             my_RT_trans%transport%time_discr,my_RT_trans%int_method_chem_reacts)
    else if (my_chem%act_coeffs_model==0 .and. my_chem%lump_flag .eqv. .false.) then !> ideal without lumping
        call my_chem%solve_reactive_mixing_ideal_cons(root_trim,my_RT_trans%transport%mixing_ratios_conc,&
            my_RT_trans%transport%mixing_ratios_Rk_init,my_RT_trans%transport%mixing_waters_indices,&
            my_RT_trans%transport%mixing_waters_indices_dom,&
            my_RT_trans%transport%time_discr,my_RT_trans%int_method_chem_reacts,my_RT_trans%transport%mixing_ratios_Rk)
    else if (my_chem%act_coeffs_model>0 .and. my_chem%lump_flag .eqv. .false.) then !> non-ideal with lumping
        call my_chem%solve_reactive_mixing_lump(root_trim,my_RT_trans%transport%mixing_ratios_conc,&
            my_RT_trans%transport%mixing_waters_indices,my_RT_trans%transport%mixing_waters_indices_dom,&
            my_RT_trans%transport%time_discr,&
            my_RT_trans%int_method_chem_reacts)
    else !> non-ideal without lumping
        call my_chem%solve_reactive_mixing_cons(root_trim,my_RT_trans%transport%mixing_ratios_conc,&
            my_RT_trans%transport%mixing_ratios_Rk_init,my_RT_trans%transport%mixing_waters_indices,&
            my_RT_trans%transport%mixing_waters_indices_dom,&
            my_RT_trans%transport%time_discr,my_RT_trans%int_method_chem_reacts,my_RT_trans%transport%mixing_ratios_Rk)
    end if
!> We set chemistry attribute in reactive transport object
    call my_RT_trans%set_chemistry(my_chem)
!> We write data and results
    call my_RT_trans%write_RT_1D(root_trim)
end program