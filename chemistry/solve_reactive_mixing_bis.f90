!> Solves reactive mixing problem
!> Computes concentrations, activities, activity coefficients, reaction rates and volumetric fractions of minerals (if present)
subroutine solve_reactive_mixing_bis(this,root,unit,mixing_ratios,mixing_waters_indices,F_mat,time_discr_tpt,int_method_chem_reacts)
    use chemistry_Lagr_m
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: root
    integer(kind=4), intent(in) :: unit
    class(real_array_c), intent(in) :: mixing_ratios !> mixing ratios matrix
    class(int_array_c), intent(in) :: mixing_waters_indices !> matrix that contains indices of target waters that mix with each target water
    real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
    class(time_discr_c), intent(in) :: time_discr_tpt !> time discretisation object (used to solve transport)
    integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
!> Variables
    integer(kind=4) :: i !> counter target waters
    integer(kind=4) :: j !> counter target solids 
    integer(kind=4) :: l !> counter reactive zones
    integer(kind=4) :: k !> counter time steps
    integer(kind=4) :: kk !> counter time steps chem_out_options
    integer(kind=4) :: ii !> counter target waters chem_out_options
    integer(kind=4) :: num_tar_wat !> number of target waters
    integer(kind=4) :: num_tar_sol !> number of target solids
    integer(kind=4), allocatable :: tar_gas_indices(:) !> indices target gases in each reactive zone
    integer(kind=4), allocatable :: tar_sol_indices(:) !> indices target solids in each reactive zone
    integer(kind=4), allocatable :: tar_wat_indices(:) !> indices target waters in each reactive zone
    REAL(KIND=8) :: time !> time
    REAL(KIND=8) :: Delta_t !> time step
    REAL(KIND=8), allocatable :: c_tilde(:) !> concentrations after mixing
    REAL(KIND=8), allocatable :: conc_old(:,:) !> concentrations before mixing
    REAL(KIND=8), allocatable :: conc_nc(:) !> concentrations variable activity species
    !REAL(KIND=8), allocatable :: conc_comp(:) !> concentrations components
    type(aqueous_chemistry_c), allocatable :: target_waters_new(:) !> target waters time step k+1
    type(aqueous_chemistry_c), allocatable :: target_waters_old(:) !> target waters time step k
    type(aqueous_chemistry_c), allocatable :: target_waters_old_old(:) !> target waters time step k-1
    type(reactive_zone_c) :: react_zone
    type(aq_phase_c), target :: aq_phase
!> Procedure pointers
    !> primary concentrations
    !procedure(get_c1_aq), pointer :: p_prim=>NULL()
    !> reactive mixing subroutines
    procedure(mixing_iter_comp), pointer :: p_solver=>null()
    !> mobile species mixing
    !procedure(compute_c_tilde), pointer :: p_c_tilde=>null()
    !p_c_tilde=>compute_c_tilde !> by default
!> We initialise target waters
    target_waters_old=this%target_waters
    target_waters_old_old=target_waters_old
    target_waters_new=target_waters_old

    time=0d0
    kk=2 !> counter time steps chem_out_options
    ii=1 !> counter target waters chem_out_options
    this%chem_out_options%time_steps(this%chem_out_options%num_time_steps)=time_discr_tpt%Num_time !> chapuza
            
    open(unit,file=root//'.out',form="formatted",access="sequential",status="unknown")
!> Time loop
    do k=1,time_discr_tpt%Num_time
        Delta_t=time_discr_tpt%get_Delta_t(k)
        time=time+Delta_t
        if (k==this%chem_out_options%time_steps(kk)) then
            write(unit,"(2x,'t = ',*(ES15.5))") time
            write(unit,"(20x,*(A15))") (this%chem_syst%aq_phase%aq_species(this%chem_out_options%ind_aq_species(j))%name, j=1,this%chem_out_options%num_aq_species)
        end if
    !> Target waters loop
        do i=this%num_ext_waters+1,this%num_target_waters
        !> We select reactive mixing subroutine depending on the nature of the chemical system and the methods to compute Jacobians and integrate in time
            if (this%chem_syst%num_kin_reacts>0) then !> equilibrium and kinetic reactions
                if (int_method_chem_reacts==1) then !> Euler explicit
                    if (this%act_coeffs_model==0) then !> ideal
                        p_solver=>water_mixing_iter_EE_eq_kin_ideal
                    else
                        p_solver=>water_mixing_iter_EE_eq_kin
                    end if
                else if (int_method_chem_reacts==2 .and. this%Jac_flag==1) then !> Euler fully implicit, analytical Jacobian
                    if (this%act_coeffs_model==0) then !> ideal
                        p_solver=>water_mixing_iter_EfI_eq_kin_anal_ideal
                    else
                        p_solver=>water_mixing_iter_EfI_eq_kin_anal
                    end if
                else
                    error stop "Integration method for chemical reactions not implemented yet"
                end if
            else if (react_zone%cat_exch_zone%num_surf_compl>0) then !> variable activity species are aqueous and solid
                p_solver=>mixing_iter_comp_exch !> only equilibrium reactions
            else !> faltan los gases en equilibrio
                if (this%act_coeffs_model==0) then !> ideal
                    p_solver=>mixing_iter_comp_ideal !> only equilibrium reactions
                else
                    p_solver=>mixing_iter_comp !> only equilibrium reactions
                end if
            end if
            !allocate(conc_nc(this%target_waters(i)%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
            allocate(conc_old(this%target_waters(i)%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,mixing_ratios%cols(i-this%num_ext_waters)%dim)) !> chapuza
            conc_old(:,1)=target_waters_old(i)%get_conc_nc() !> chapuza
            do j=1,mixing_waters_indices%cols(i-this%num_ext_waters)%dim
                conc_old(:,1+j)=target_waters_old(mixing_waters_indices%cols(i-this%num_ext_waters)%col_1(j))%get_conc_nc() !> chapuza
            end do
            call target_waters_new(i)%solve_reactive_mixing_iter(target_waters_old_old(i)%get_c1(),mixing_ratios%cols(i-this%num_ext_waters)%col_1,conc_old,F_mat(i-this%num_ext_waters),Delta_t,p_solver)
        !> Deallocate
            deallocate(conc_old)
        end do
    !> We update target waters
        target_waters_old_old=target_waters_old
        target_waters_old=target_waters_new
    end do
!> We set the new target waters to the chemistry object
    this%target_waters=target_waters_new
    close(unit)
end subroutine