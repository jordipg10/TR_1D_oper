!> Solves reactive mixing problem
!> Computes concentrations, activities, activity coefficients, reaction rates and volumetric fractions of minerals (if present)
subroutine solve_reactive_mixing(this,mixing_ratios,mixing_waters_indices,F_mat,time_discr_tpt,int_method_chem_reacts)
    use chemistry_Lagr_m
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    class(matrix_real_c), intent(in) :: mixing_ratios !> mixing ratios matrix
    class(matrix_int_c), intent(in) :: mixing_waters_indices !> matrix that contains indices of target waters that mix with each target water
    class(diag_matrix_c), intent(in) :: F_mat !> storage matrix (diagonal)
    class(time_discr_c), intent(in) :: time_discr_tpt !> time discretisation object (used to solve transport)
    integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
!> Variables
    integer(kind=4) :: i !> counter target waters
    integer(kind=4) :: j !> counter target solids 
    integer(kind=4) :: l !> counter reactive zones
    integer(kind=4) :: k !> counter time steps
    integer(kind=4) :: num_tar_wat !> number of target waters
    integer(kind=4) :: num_tar_sol !> number of target solids
    integer(kind=4), allocatable :: tar_gas_indices(:) !> indices target gases in each reactive zone
    integer(kind=4), allocatable :: tar_sol_indices(:) !> indices target solids in each reactive zone
    integer(kind=4), allocatable :: tar_wat_indices(:) !> indices target waters in each reactive zone
    REAL(KIND=8), allocatable :: c_tilde(:) !> concentrations after mixing
    REAL(KIND=8), allocatable :: conc_old(:,:) !> concentrations before mixing
    REAL(KIND=8), allocatable :: conc_nc(:,:) !> concentrations variable activity species
    REAL(KIND=8), allocatable :: conc_comp(:,:) !> concentrations components
    type(aqueous_chemistry_c), allocatable :: target_waters_new(:) !> target waters time step k+1
    type(aqueous_chemistry_c), allocatable :: target_waters_old(:) !> target waters time step k
    type(aqueous_chemistry_c), allocatable :: target_waters_old_old(:) !> target waters time step k-1
!> Procedure pointers
    !> primary concentrations
    procedure(get_c1_aq), pointer :: p_prim=>NULL()
    !> reactive mixing subroutines
    procedure(transport_iter_comp_EE_aq_chem), pointer :: p_solver=>null()
    !> mobile species mixing
    procedure(compute_c_tilde_aq_chem), pointer :: p_c_tilde=>null()
    p_c_tilde=>compute_c_tilde_aq_chem !> by default
!> We initialise target waters
    target_waters_old=this%target_waters
    !print *, target_waters_old(1)%gas_chemistry%concentrations
    target_waters_old_old=target_waters_old
    target_waters_new=target_waters_old
    !print *, target_waters_new(1)%gas_chemistry%concentrations
!> We select reactive mixing subroutine depending on the nature of the chemical system and the methods to compute Jacobians and integrate in time
    if (this%num_reactive_zones>0) then !> there are reactive zones
    !> Loop in reactive zones
        do l=1,this%num_reactive_zones
            if (this%chem_syst%num_kin_reacts>0) then !> equilibrium and kinetic reactions
                if (this%reactive_zones(l)%cat_exch_zone%num_surf_compl==0) then
                    p_prim=>get_c1_aq !> chapuza
                else
                    p_prim=>get_c1_exch
                end if
                if (int_method_chem_reacts==1) then !> Euler explicit
                    p_solver=>water_mixing_iter_EE_eq_kin
                else if (int_method_chem_reacts==2 .and. this%Jac_flag==1) then !> Euler fully implicit, analytical Jacobian
                    p_solver=>water_mixing_iter_EfI_eq_kin_anal
                else
                    error stop "Integration method for chemical reactions not implemented yet"
                end if
            !else if (this%reactive_zones(l)%cat_exch_zone%num_surf_compl==0 .AND. this%reactive_zones(l)%gas_phase%num_species==0) then !> all variable activity species are aqueous
            else if (this%reactive_zones(l)%cat_exch_zone%num_surf_compl==0) then !> all variable activity species are aqueous
                p_solver=>transport_iter_comp_EE_aq_chem !> only equilibrium reactions
                p_prim=>get_c1_aq
            !else if (this%reactive_zones(l)%cat_exch_zone%num_surf_compl>0 .AND. this%reactive_zones(l)%gas_phase%num_species==0) then !> variable activity species are aqueous and solid
            else if (this%reactive_zones(l)%cat_exch_zone%num_surf_compl>0) then !> variable activity species are aqueous and solid
                p_solver=>transport_iter_comp_exch_EE_aq_chem !> only equilibrium reactions
                p_prim=>get_c1_exch
            else !> faltan los gases en equilibrio
                p_solver=>transport_iter_comp_EE_aq_chem !> only equilibrium reactions
            end if
        !> We link indices of target solids (tar_sol_indices) associated to this reactive zone
            !call this%link_target_solids_reactive_zone(l,tar_sol_indices)
        !> We link indices of target gases (tar_gas_indices) associated to this reactive zone
            !call this%link_target_gases_reactive_zone(l,tar_gas_indices)
        !> We link indices of target waters (tar_wat_indices) associated to tar_sol_indices
            !call this%link_target_waters_target_solids(tar_sol_indices,tar_wat_indices)
        !> We link indices of target waters (tar_wat_indices) associated to tar_gas_indices
            !call this%link_target_waters_target_gases(tar_gas_indices,tar_wat_indices)
        !> Number of target waters mixing in this reactive zone
            !num_tar_wat=size(tar_wat_indices)
        !> Number of target solids associated to this reactive zone (we assume = num_tar_wat)
            !num_tar_sol=size(tar_sol_indices) 
        !> Allocate
            allocate(conc_nc(this%target_waters_init(tar_wat_indices(1))%speciation_alg%num_var_act_species,num_tar_wat))
            allocate(conc_comp(this%target_waters_init(tar_wat_indices(1))%speciation_alg%num_prim_species,num_tar_wat))
        !> Time loop
            do k=1,time_discr_tpt%Num_time
            !> Target waters loop
                do i=1,num_tar_wat
                    allocate(conc_old(target_waters_new(tar_wat_indices(i))%speciation_alg%num_var_act_species,mixing_ratios%cols(tar_wat_indices(i))%dim)) !> chapuza
                    conc_old(:,1)=target_waters_old(tar_wat_indices(i))%get_conc_nc() !> chapuza
                    do j=1,mixing_waters_indices%cols(tar_wat_indices(i))%dim
                        conc_old(:,1+j)=target_waters_old(mixing_waters_indices%cols(tar_wat_indices(i))%col_1(j))%get_conc_nc() !> chapuza
                    end do
                !> We solve mixing caused by transport
                    c_tilde=p_c_tilde(target_waters_old(tar_wat_indices(I)),mixing_ratios%cols(tar_wat_indices(i))%col_1,conc_old)
                !> We solve reactive mixing iteration
                    call p_solver(target_waters_new(tar_wat_indices(i)),p_prim(target_waters_old_old(tar_wat_indices(i))),target_waters_old(tar_wat_indices(i))%get_c2nc(),c_tilde,conc_nc(:,tar_wat_indices(i)),conc_comp(:,tar_wat_indices(i)),F_mat%diag(tar_wat_indices(i)),time_discr_tpt%get_Delta_t(k))
                !> We compute equilibrium reaction rates from mass balance equation
                    call target_waters_new(tar_wat_indices(i))%compute_r_eq_aq_chem(c_tilde(target_waters_new(tar_wat_indices(i))%speciation_alg%num_prim_species+1:target_waters_new(tar_wat_indices(i))%speciation_alg%num_var_act_species),time_discr_tpt%get_Delta_t(k),F_mat%diag(tar_wat_indices(i)))
                !> Chapuza
                    if (associated(target_waters_new(tar_wat_indices(i))%solid_chemistry)) then
                    !> We compute mass volumetric fractions of minerals from mass balance equation
                        call target_waters_new(tar_wat_indices(i))%solid_chemistry%compute_mass_bal_mins(time_discr_tpt%get_Delta_t(k))
                    !> We compute concentrations of minerals
                        call target_waters_new(tar_wat_indices(i))%solid_chemistry%compute_conc_minerals_iter(time_discr_tpt%get_Delta_t(k))
                    end if
                    if (associated(target_waters_new(tar_wat_indices(i))%gas_chemistry)) then
                        !print *, target_waters_new(tar_wat_indices(i))%gas_chemistry%concentrations
                    !> We compute concentrations of gases
                        call target_waters_new(tar_wat_indices(i))%gas_chemistry%compute_conc_gases_iter(time_discr_tpt%get_Delta_t(k),F_mat%diag(tar_wat_indices(i)),target_waters_new(tar_wat_indices(i))%volume)
                    !> We compute total volume of gas
                        !call target_waters_new(tar_wat_indices(i))%gas_chemistry%compute_vol_gas()
                    !> We compute activity coefficients of gases    
                        call target_waters_new(tar_wat_indices(i))%gas_chemistry%compute_log_act_coeffs_gases()
                    end if
                !> Deallocate
                    deallocate(c_tilde,conc_old)
                end do
            !> We update target waters
                target_waters_old_old=target_waters_old
                target_waters_old=target_waters_new
            end do
        !> Deallocate indices
            deallocate(conc_nc,conc_comp)
        end do
    else !> there are NO reactive zones
        p_c_tilde=>compute_c_tilde_aq_chem
        if (this%chem_syst%num_eq_reacts>0) then 
            if (this%chem_syst%num_kin_reacts>0) then !> equilibrium and kinetic reactions
                if (int_method_chem_reacts==1) then !> Euler explicit
                    p_solver=>water_mixing_iter_EE_eq_kin
                else if (int_method_chem_reacts==2 .and. this%Jac_flag==1) then !> Euler fully implicit
                    p_solver=>water_mixing_iter_EfI_eq_kin_anal
                else
                    error stop "Integration method for chemical reactions not implemented yet"
                end if
            else
                p_solver=>transport_iter_comp_EE_aq_chem !> only equilibrium reactions
            end if
        else if (this%chem_syst%num_kin_reacts>0) then !> only kinetic reactions 
            if (int_method_chem_reacts==1) then !> Euler explicit
                p_solver=>water_mixing_iter_EE_kin 
            else if (int_method_chem_reacts==2 .and. this%Jac_flag==1) then !> Euler fully implicit
                p_solver=>water_mixing_iter_EfI_kin_anal
            else
                error stop "Integration method for chemical reactions not implemented yet"
            end if
        else !> no chemical reactions
            p_solver=>transport_iter_species_EE_aq_chem !> conservative transport
        end if
        allocate(tar_wat_indices(this%num_target_waters))
        do i=1,this%num_target_waters
            tar_wat_indices(i)=i
        end do
    !> Time loop
        do k=1,time_discr_tpt%Num_time
        !> Target waters loop
            do i=1,num_tar_wat
                allocate(conc_old(target_waters_new(tar_wat_indices(i))%speciation_alg%num_var_act_species,mixing_ratios%cols(tar_wat_indices(i))%dim)) !> chapuza
                conc_old(:,1)=target_waters_old(tar_wat_indices(i))%get_conc_nc() !> chapuza
                do j=1,mixing_waters_indices%cols(tar_wat_indices(i))%dim
                    conc_old(:,1+j)=target_waters_old(mixing_waters_indices%cols(tar_wat_indices(i))%col_1(j))%get_conc_nc() !> chapuza
                end do
            !> We solve mixing caused by transport
                c_tilde=p_c_tilde(target_waters_old(tar_wat_indices(I)),mixing_ratios%cols(tar_wat_indices(i))%col_1,conc_old)
            !> We solve reactive mixing iteration
                call p_solver(target_waters_new(tar_wat_indices(i)),p_prim(target_waters_old_old(tar_wat_indices(i))),target_waters_old(tar_wat_indices(i))%get_c2nc(),c_tilde,conc_nc(:,tar_wat_indices(i)),conc_comp(:,tar_wat_indices(i)),F_mat%diag(tar_wat_indices(i)),time_discr_tpt%get_Delta_t(k))
            end do
        !> We update target waters
            target_waters_old_old=target_waters_old
            target_waters_old=target_waters_new
        end do
    !> Deallocate indices
        deallocate(tar_wat_indices)
    end if
!> We set the new target waters to the chemistry object
    this%target_waters=target_waters_new
end subroutine