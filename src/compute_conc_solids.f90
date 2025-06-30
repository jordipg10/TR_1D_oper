!> Computes concentration of solids at a given time
!!> phi_j*dc_s/dt=S_s^T*r_j
subroutine compute_conc_solids(this,r_vec,time) 
    use solid_chemistry_m, only: solid_chemistry_c
    use vectors_m
    implicit none
    class(solid_chemistry_c) :: this
    real(kind=8), intent(in) :: r_vec(:) !> reaction rates
    real(kind=8), intent(in) :: time

    integer(kind=4) :: i
    
    do i=1,this%reactive_zone%num_solids
        this%concentrations(i)=this%concentrations(i)+time*dot_product(this%reactive_zone%stoich_mat_sol(:,i),r_vec)/&
        this%vol_fracts(i)
    end do
           
end subroutine