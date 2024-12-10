!> Computes concentration of secondary species with variable activity from concentration of primary species explicitly using mass action law
!!  We assume the activity coefficients are ideal
!! This subroutine is meant to be used only when all primary species are aqueous
subroutine compute_c2nc_from_c1_aq_ideal(this,c2nc)
    use aqueous_chemistry_m
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(out) :: c2nc(:) !> secondary variable activity concentrations (must be already allocated)
    
    integer(kind=4) :: n_p,n_nc2_aq,n_e
    real(kind=8), allocatable :: log_c2nc(:)
    
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_nc2_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species

    log_c2nc=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log10(this%concentrations(1:n_p)))+this%solid_chemistry%reactive_zone%speciation_alg%logK_star
    c2nc=10**log_c2nc
    call this%update_conc_sec_aq_var_act_species(c2nc(1:n_nc2_aq))
    if (associated(this%gas_chemistry)) then
        if (this%gas_chemistry%reactive_zone%gas_phase%num_var_act_species>0) then !> chapuza
            call this%gas_chemistry%update_conc_gases(c2nc(n_nc2_aq+1:n_e)*this%volume) !> we update moles of gases
            !call this%gas_chemistry%compute_vol_gas() !> we compute total volume of gas
        end if
    end if
 end subroutine 