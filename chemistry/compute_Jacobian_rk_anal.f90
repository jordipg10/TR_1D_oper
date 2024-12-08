!> Computes Jacobian of kinetic reaction rates with respect to aqueous concentrations
subroutine compute_Jacobian_rk_anal(this,drk_dc)
    use aqueous_chemistry_m
    implicit none
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(out) :: drk_dc(:,:) !> Jacobian (must be allocated)
    
    integer(kind=4) :: i,k,n,m,l,rk_ind,cat_ind,p,inh_ind,DOC_ind
    integer(kind=4), allocatable :: indices(:)
    real(kind=8), allocatable :: act_cat(:),drk_dc_loc(:)
    real(kind=8) :: saturation

    drk_dc=0d0 !> by default
!> Linear kinetic reactions
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_drk_dc_lin(drk_dc(i,:))
    end do
!> Mineral kinetic reactions
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_min_kin_reacts
        indices=this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%indices_aq_phase
        allocate(drk_dc_loc(SIZE(indices)))
        saturation=this%compute_saturation_min(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i))
        call this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%compute_drk_dc_mineral(this%concentrations(this%indices_aq_phase(indices)),this%activities(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%params%cat_indices),saturation,this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,drk_dc_loc)
        drk_dc(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i,indices)=drk_dc_loc !> chapuza
        deallocate(indices,drk_dc_loc)
    end do
!> Redox kinetic reactions
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        indices=this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase
        allocate(drk_dc_loc(SIZE(indices)))
        call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_drk_dc_Monod(this%concentrations(this%indices_aq_phase(indices)),this%rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+this%solid_chemistry%reactive_zone%chem_syst%num_min_kin_reacts+i),drk_dc_loc)
        drk_dc(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+this%solid_chemistry%reactive_zone%chem_syst%num_min_kin_reacts+i,indices)=drk_dc_loc !> chapuza
        deallocate(drk_dc_loc)
    end do
end subroutine