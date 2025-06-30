!!> Computes concentration of solids after a given time step
!!!> The mass balance equation is: \f$ phi_j*d{c_s}/dt=S_s^T*r_{eq_j}
!subroutine compute_conc_solids_iter_EE(this,r_eq,Delta_t) 
!    use solid_chemistry_Lagr_m
!    use vectors_m
!    implicit none
!    class(solid_chemistry_c) :: this
!    real(kind=8), intent(in) :: r_eq(:) !> equilibrium reaction rates at previous time step
!    real(kind=8), intent(in) :: Delta_t !> time step
!
!    integer(kind=4) :: i
!    
!    do i=1,this%reactive_zone%num_minerals
!        this%concentrations(i)=this%concentrations(i)+Delta_t*dot_product(this%reactive_zone%stoich_mat_sol(:,i),r_eq)/this%vol_fracts(i) !> esto está mal
!    end do
!           
!end subroutine