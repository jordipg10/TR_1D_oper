!> Computes concentration of minerals after a given time step
subroutine compute_conc_minerals_iter(this,Delta_t) 
    use solid_chemistry_m
!> Arguments
    implicit none
    class(solid_chemistry_c) :: this !> solid chemistry object
    !real(kind=8), intent(in) :: r_eq(:) !> equilibrium reaction rates at previous time step
    real(kind=8), intent(in) :: Delta_t !> time step
!> Variables
    integer(kind=4) :: i !> counter minerals
    real(kind=8), parameter :: eps=1d-16
!> Process
    !print *, this%concentrations, this%r_eq
    do i=1,this%reactive_zone%num_minerals
        if (abs(this%vol_fracts(this%reactive_zone%chem_syst%num_min_kin_reacts+i))<eps) then
            continue
        else
            this%concentrations(this%reactive_zone%chem_syst%num_min_kin_reacts+i)=this%concentrations(this%reactive_zone%chem_syst%num_min_kin_reacts+i)+Delta_t*dot_product(this%reactive_zone%stoich_mat_sol(:,this%reactive_zone%cat_exch_zone%num_surf_compl+i),this%r_eq(1:this%reactive_zone%num_minerals))/this%vol_fracts(this%reactive_zone%chem_syst%num_min_kin_reacts+i)
        end if
    end do
    do i=1,this%reactive_zone%chem_syst%num_min_kin_reacts
        if (abs(this%vol_fracts(i))<eps) then
            continue
        else
            this%concentrations(i)=this%concentrations(i)+Delta_t*dot_product(this%reactive_zone%chem_syst%stoich_mat_sol(:,this%reactive_zone%chem_syst%cat_exch%num_surf_compl+i),this%rk(1:this%reactive_zone%chem_syst%num_min_kin_reacts))/this%vol_fracts(i)
        end if
    end do
end subroutine