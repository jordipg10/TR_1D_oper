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
!> Process
    do i=1,this%reactive_zone%num_minerals
        this%concentrations(i)=this%concentrations(i)+Delta_t*dot_product(this%reactive_zone%stoich_mat_sol(:,this%reactive_zone%cat_exch_zone%num_surf_compl+i),this%r_eq(1:this%reactive_zone%num_minerals))/this%vol_fracts(i)
    end do
end subroutine