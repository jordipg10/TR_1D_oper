!> Computes kinetic reaction rates associated to aqueous chemistry object
subroutine compute_rk(this)
    use aqueous_chemistry_m
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
!> Variables
    integer(kind=4) :: i,n,niter,rk_ind,l,index
    integer(kind=4), allocatable :: indices(:)
    real(kind=8) :: saturation
!> We compute linear kinetic reaction rates
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        index=this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1)
        call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_rk_lin(this%concentrations(this%indices_aq_phase(index)),this%rk(i))
    end do
!> We compute mineral kinetic reaction rates
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_min_kin_reacts
        indices=this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%indices_aq_phase
        saturation=this%compute_saturation_min(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i))
        call this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%compute_rk_mineral(this%activities(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%params%cat_indices),saturation,this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,this%rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i))
        deallocate(indices)
    end do
!> We compute Monod reaction rates
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        indices=this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase
        call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_rk_Monod(this%concentrations(this%indices_aq_phase(indices)),this%rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+this%solid_chemistry%reactive_zone%chem_syst%num_min_kin_reacts+i))
        deallocate(indices)
    end do
end subroutine