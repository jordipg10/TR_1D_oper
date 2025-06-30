module flow_transient_m
    use PDE_transient_m, only: PDE_1D_transient_c, spatial_discr_c, spatial_discr_rad_c
    use flow_props_heterog_m, only: flow_props_heterog_c, flow_props_heterog_conf_c
    use stab_params_flow_m, only: stab_params_flow_c
    implicit none
    save
    type, public, extends(PDE_1D_transient_c) :: flow_transient_c !> 1D transient flow subclass
        real(kind=8), allocatable :: head(:)                    !> hydraulic head
        real(kind=8), allocatable :: head_init(:)               !> initial hydraulic head
        type(flow_props_heterog_conf_c) :: flow_props_heterog       !> properties
        type(stab_params_flow_c) :: stab_params_flow            !> stability parameters
    contains
        procedure, public :: set_head
        procedure, public :: set_head_init
        procedure, public :: compute_F_mat_PDE=>compute_F_mat_flow
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_flow_transient
        procedure, public :: set_stab_params_flow
        procedure, public :: initialise_PDE=>initialise_flow_transient
        procedure, public :: compute_source_term_PDE=>compute_source_term_flow
        procedure, public :: compute_rech_mat_PDE=>compute_rech_mat_flow
        !procedure, public :: initialise_flow_1D_transient_RT
        procedure, public :: write_PDE_1D=>write_flow_transient
        procedure, public :: set_flow_props_heterog
        !procedure, public :: mass_balance_error_ADE_trans_PMF_evap
        !procedure, public :: mass_balance_error_ADE_trans_Dirichlet_evap
        !procedure, public :: mass_balance_error_ADE_trans_PMF_recharge
        !procedure, public :: mass_balance_error_ADE_trans_Dirichlet_recharge
        !procedure, public :: mass_balance_error_ADE_trans_PMF_discharge
        !procedure, public :: mass_balance_error_ADE_trans_Dirichlet_discharge
        !procedure, public :: check_Delta_t
        !procedure, public :: read_flow_data_WMA
        !procedure, public :: write_flow_data_WMA
        !procedure, public :: fund_sol_flow_eqn_1D
        procedure, public :: solve_PDE_1D=>solve_flow_transient
        procedure, public :: solve_flow_EE_Delta_t_homog
    end type
    
    interface
        subroutine compute_F_mat_flow(this)
            import flow_transient_c
            implicit none
            class(flow_transient_c) :: this
        end subroutine
        
        subroutine compute_trans_mat_flow_transient(this)
            import flow_transient_c
            implicit none
            class(flow_transient_c) :: this
        end subroutine
        
        subroutine initialise_flow_transient(this,root)
            import flow_transient_c
            implicit none
            class(flow_transient_c) :: this
            character(len=*), intent(in) :: root
        end subroutine
        
        subroutine solve_flow_transient(this,Time_out,output)
        import flow_transient_c
        implicit none
        class(flow_transient_c) :: this !> flow object
        real(kind=8), intent(in) :: Time_out(:)
        real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine solve_flow_EE_Delta_t_homog(this,Time_out,output)
        import flow_transient_c
        implicit none
        class(flow_transient_c) :: this !> flow object
        real(kind=8), intent(in) :: Time_out(:)
        real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        ! subroutine initialise_flow_transient_RT(this,root)
        !     import flow_transient_c
        !     implicit none
        !     class(flow_transient_c) :: this
        !     !character(len=*), intent(in) :: path
        !     character(len=*), intent(in) :: root
        !     !character(len=*), intent(in) :: file_BCs
        !     !character(len=*), intent(in) :: file_spatial_discr
        !     !character(len=*), intent(in) :: file_time_discr
        !     !character(len=*), intent(in) :: file_flow_props
        ! end subroutine
        
        subroutine write_flow_transient(this,Time_out,output)
            import flow_transient_c
            implicit none
            class(flow_transient_c), intent(in) :: this !> flow object
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
        end subroutine
        
        function mass_balance_error_ADE_trans_Dirichlet_recharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import flow_transient_c
            implicit none
            class(flow_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_Dirichlet_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import flow_transient_c
            implicit none
            class(flow_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_Dirichlet_evap(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import flow_transient_c
            implicit none
            class(flow_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_PMF_recharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import flow_transient_c
            implicit none
            class(flow_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_PMF_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import flow_transient_c
            implicit none
            class(flow_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_PMF_evap(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import flow_transient_c
            implicit none
            class(flow_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        subroutine read_flow_data_WMA(this,root)!,mixing_ratios)!,f_vec)!,flow_props,BCs,mesh,time_discr)
            import flow_transient_c
            !import flow_props_heterog_c
            !import BCs_t
            !import mesh_1D_Euler_homog_c
            !import time_discr_homog_c
            class(flow_transient_c) :: this
            !character(len=*), intent(in) :: path
            !integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: root
            !real(kind=8), intent(out), allocatable :: mixing_ratios(:,:)
            !real(kind=8), intent(out), allocatable :: f_vec(:)
            !type(flow_props_heterog_c), intent(out) :: flow_props
            !type(BCs_t), intent(out) :: BCs
            !!real(kind=8), intent(in) :: Delta_x
            !type(mesh_1D_Euler_homog_c), intent(out) :: mesh !> homogeneous Euler mesh 1D
            !type(time_discr_homog_c), intent(out) :: time_discr !> homogeneous time discretisation
        end subroutine
        
        subroutine read_discretisation_WMA(this,unit,root)
            import flow_transient_c
            class(flow_transient_c) :: this
            !character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: root
        end subroutine
        
        subroutine write_flow_data_WMA(this,unit)
            import flow_transient_c
            implicit none
            class(flow_transient_c), intent(in) :: this                 !> 1D transient flow object                         
            integer(kind=4), intent(in) :: unit                                 !> unit of output file
        end subroutine
    end interface
    
    contains
      
    subroutine set_head(this,head)
    class(flow_transient_c) :: this
    real(kind=8), intent(in) :: head(:)
    this%head=head
    end subroutine
    
    subroutine set_head_init(this,head)
    class(flow_transient_c) :: this
    real(kind=8), intent(in) :: head(:)
    this%head_init=head
    end subroutine
    
        subroutine set_stab_params_flow(this,stab_params_flow)
            implicit none
            class(flow_transient_c) :: this
            class(stab_params_flow_c), intent(in) :: stab_params_flow
            this%stab_params_flow=stab_params_flow
        end subroutine
        
        subroutine set_flow_props_heterog(this,flow_props_heterog)
            implicit none
            class(flow_transient_c) :: this
            class(flow_props_heterog_conf_c), intent(in) :: flow_props_heterog
            this%flow_props_heterog=flow_props_heterog
        end subroutine
        
        !subroutine set_conc_r_flag_flow(this)
        !    implicit none
        !    class(flow_transient_c) :: this
        !    
        !    integer(kind=4) :: i
        !    
        !    allocate(this%conc_r_flag(this%spatial_discr%Num_targets))
        !    this%conc_r_flag=0
        !    do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        !        if (this%flow_props_heterog%source_term(i)>0) then
        !            this%conc_r_flag(i)=1
        !        end if
        !    end do
        !end subroutine
        
        !subroutine check_Delta_t(this)
        !    implicit none
        !    class(flow_transient_c) :: this
        !    select type (time=>this%time_discr)
        !    type is (time_discr_homog_c)
        !        if (time%Delta_t>this%stab_params_flow%Delta_t_crit) then
        !            !print *, "Critical time step: ", this%stab_params_flow%Delta_t_crit
        !            !error stop "You must reduce time step to have stability"
        !        end if
        !    end select
        !end subroutine

        !function fund_sol_flow_eqn_1D(this,c0,M,Delta_x,x,mu,t) result(conc)
        !    class(flow_transient_c), intent(in) :: this
        !    real(kind=8), intent(in) :: c0
        !    real(kind=8), intent(in) :: M
        !    real(kind=8), intent(in) :: Delta_x
        !    real(kind=8), intent(in) :: x
        !    real(kind=8), intent(in) :: mu
        !    real(kind=8), intent(in) :: t
        !    !real(kind=8), intent(in) :: phi
        !    !real(kind=8), intent(in) :: D
        !    !real(kind=8), intent(in) :: q
        !    real(kind=8) :: conc
        !    
        !    real(kind=8), parameter :: eps=1d-12
        !
        !        conc=exp(5d-1*(this%flow_props_heterog%flux(1)/this%flow_props_heterog%dispersion(1))*&
        !            (x-mu-t*this%flow_props_heterog%flux(1)/(2*this%flow_props_heterog%porosity(1))))
        !        conc=conc*fund_sol_diff_eqn_1D(this,M,Delta_x,x,mu,t)+c0
        !end function
        
        subroutine Robin_Neumann_homog_BCs(this)
         !> Imposes Robin BC inflow & Neumann homogeneous BC outflow in transition matrix & sink/source term
            implicit none
            class(flow_transient_c) :: this
            
            integer(kind=4) :: n
            real(kind=8) :: q_1,q_32,q_n,q_inf,q_out,D_1,D_n,c_inf
            
            n=this%spatial_discr%Num_targets
            
            !!select type (this)
            !!type is (transport_1D_transient_c)
            !    q_1=this%tpt_props_heterog%flux(1)
            !    q_32=this%tpt_props_heterog%flux(2)
            !    q_n=this%tpt_props_heterog%flux(n)
            !    q_inf=this%BCs%flux_inf
            !    q_out=this%BCs%flux_out
            !    D_1=this%tpt_props_heterog%dispersion(1)
            !    D_n=this%tpt_props_heterog%dispersion(n)
            !    c_inf=this%BCs%conc_inf
            !    select type (mesh=>this%spatial_discr)
            !    type is (mesh_1D_Euler_homog_c)
            !        if (this%tpt_props_heterog%homog_flag.eqv..true. .and. mesh%scheme<3 .and. mesh%targets_flag==0) then !> CCFD
            !            this%trans_mat%super(1)=-q_1/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
            !            this%trans_mat%sub(n-1)= D_n/(mesh%Delta_x**2) + q_n/(2d0*mesh%Delta_x)
            !            this%trans_mat%diag(1)=this%trans_mat%diag(1)+(2d0*D_1-q_inf*mesh%Delta_x)*(q_1*mesh%Delta_x+2*D_1)/(&
            !            2*q_inf*mesh%Delta_x**3+4*D_1*mesh%Delta_x**2) - 2*D_1/(mesh%Delta_x**2)
            !            this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
            !            !this%source_term_PDE(1)=this%source_term_PDE(1)+q_inf*c_inf*(q_1*mesh%Delta_x+2*D_1)/(&
            !            !    q_inf*mesh%Delta_x**2+2*D_1*mesh%Delta_x)
            !            this%bd_mat(1)=q_inf*(q_1*mesh%Delta_x+2*D_1)/(&
            !                q_inf*mesh%Delta_x**2+2*D_1*mesh%Delta_x)
            !            this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
            !        else if ((mesh%scheme.eq.2 .and. mesh%targets_flag.eq.0) .AND. (this%tpt_props_heterog%homog_flag.eqv..false.))&
            !            then !> proposed scheme
            !            this%trans_mat%super(1)=-q_32/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
            !            this%trans_mat%sub(n-1)= q_n/(2d0*mesh%Delta_x) + D_n/(mesh%Delta_x**2)
            !            !this%trans_mat%diag(1)=this%trans_mat%diag(1)-(q_inf**2+D_1*(3*q_inf*mesh%Delta_x+2*D_1)/(mesh%Delta_x**2))/(q_inf*mesh%Delta_x+2*D_1) + q_32/(2*mesh%Delta_x)
            !            this%trans_mat%diag(1)=this%trans_mat%diag(1) - q_inf/mesh%Delta_x + q_32/(2*mesh%Delta_x) - &
            !                D_1/(mesh%Delta_x**2)
            !            this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
            !            !this%source_term_PDE(1)=this%source_term_PDE(1)+q_inf*c_inf/mesh%Delta_x
            !            this%bd_mat(1)=q_inf/mesh%Delta_x
            !            this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
            !        else if (mesh%targets_flag.eq.1 .AND. this%tpt_props_heterog%homog_flag.eqv..true.) then !> edge centred
            !            this%trans_mat%super(1)=2d0*D_1/(mesh%Delta_x**2)
            !            this%trans_mat%sub(n-1)=2d0*D_n/(mesh%Delta_x**2)
            !            this%trans_mat%diag(1)=this%trans_mat%diag(1) - q_inf*(q_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x) - &
            !                2*D_1/(mesh%Delta_x**2)
            !            this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
            !            !this%source_term_PDE(1)=this%source_term_PDE(1) + q_inf*(q_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x)
            !            this%bd_mat(1)=q_inf*(q_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x)
            !            this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
            !        end if
            !    end select
            !end select
        end subroutine
        
        subroutine compute_source_term_flow(this)
        implicit none
        class(flow_transient_c) :: this

        integer(kind=4) :: i

        !> Initialise source term
        allocate(this%source_term_PDE(this%spatial_discr%Num_targets))
        this%source_term_PDE=0d0

        !> Compute source term
        do i=1,this%spatial_discr%Num_targets
            this%source_term_PDE(i)=this%flow_props_heterog%source_term(i) !> incomplete source term
        end do

        end subroutine
        
        subroutine compute_rech_mat_flow(this)
        implicit none
        class(flow_transient_c) :: this
        this%rech_mat%diag=this%flow_props_heterog%source_term
        end subroutine
        
        subroutine Robin_Neumann_homog_BCs_flow(this)
         !> Imposes Robin BC inflow & Neumann homogeneous BC outflow in transition matrix & sink/source term
            implicit none
            class(flow_transient_c) :: this
            
            integer(kind=4) :: n
            real(kind=8) :: q_1,q_32,q_n,q_inf,q_out,D_1,D_n,c_inf
            
            n=this%spatial_discr%Num_targets
            
            !select type (this)
            !type is (transport_1D_transient_c)
                !q_1=this%tpt_props_heterog%flux(1)
                !q_32=this%tpt_props_heterog%flux(2)
                !q_n=this%tpt_props_heterog%flux(n)
                !q_inf=this%BCs%flux_inf
                !q_out=this%BCs%flux_out
                !D_1=this%tpt_props_heterog%dispersion(1)
                !D_n=this%tpt_props_heterog%dispersion(n)
                !c_inf=this%BCs%conc_inf
                select type (mesh=>this%spatial_discr)
                type is (spatial_discr_rad_c)
                    
                !type is (mesh_1D_Euler_homog_c)
                !    if (this%tpt_props_heterog%homog_flag.eqv..true. .and. mesh%scheme<3 .and. mesh%targets_flag==0) then !> CCFD
                !        this%trans_mat%super(1)=-q_1/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
                !        this%trans_mat%sub(n-1)= D_n/(mesh%Delta_x**2) + q_n/(2d0*mesh%Delta_x)
                !        this%trans_mat%diag(1)=this%trans_mat%diag(1)+(2d0*D_1-q_inf*mesh%Delta_x)*(q_1*mesh%Delta_x+2*D_1)/(&
                !        2*q_inf*mesh%Delta_x**3+4*D_1*mesh%Delta_x**2) - 2*D_1/(mesh%Delta_x**2)
                !        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                !        !this%source_term_PDE(1)=this%source_term_PDE(1)+q_inf*c_inf*(q_1*mesh%Delta_x+2*D_1)/(&
                !        !    q_inf*mesh%Delta_x**2+2*D_1*mesh%Delta_x)
                !        this%bd_mat(1)=q_inf*(q_1*mesh%Delta_x+2*D_1)/(&
                !            q_inf*mesh%Delta_x**2+2*D_1*mesh%Delta_x)
                !        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                !    else if ((mesh%scheme.eq.2 .and. mesh%targets_flag.eq.0) .AND. (this%tpt_props_heterog%homog_flag.eqv..false.))&
                !        then !> proposed scheme
                !        this%trans_mat%super(1)=-q_32/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
                !        this%trans_mat%sub(n-1)= q_n/(2d0*mesh%Delta_x) + D_n/(mesh%Delta_x**2)
                !        !this%trans_mat%diag(1)=this%trans_mat%diag(1)-(q_inf**2+D_1*(3*q_inf*mesh%Delta_x+2*D_1)/(mesh%Delta_x**2))/(q_inf*mesh%Delta_x+2*D_1) + q_32/(2*mesh%Delta_x)
                !        this%trans_mat%diag(1)=this%trans_mat%diag(1) - q_inf/mesh%Delta_x + q_32/(2*mesh%Delta_x) - &
                !            D_1/(mesh%Delta_x**2)
                !        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                !        !this%source_term_PDE(1)=this%source_term_PDE(1)+q_inf*c_inf/mesh%Delta_x
                !        this%bd_mat(1)=q_inf/mesh%Delta_x
                !        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                !    else if (mesh%targets_flag.eq.1 .AND. this%tpt_props_heterog%homog_flag.eqv..true.) then !> edge centred
                !        this%trans_mat%super(1)=2d0*D_1/(mesh%Delta_x**2)
                !        this%trans_mat%sub(n-1)=2d0*D_n/(mesh%Delta_x**2)
                !        this%trans_mat%diag(1)=this%trans_mat%diag(1) - q_inf*(q_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x) - &
                !            2*D_1/(mesh%Delta_x**2)
                !        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                !        !this%source_term_PDE(1)=this%source_term_PDE(1) + q_inf*(q_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x)
                !        this%bd_mat(1)=q_inf*(q_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x)
                !        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                !    end if
                end select
            !end select
         end subroutine

end module