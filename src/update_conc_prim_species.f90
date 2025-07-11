!> Updates concentration aqueous and solid primary species in iterative method
subroutine update_conc_prim_species(this,c1,Delta_c1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: c1(:)
    real(kind=8), intent(inout) :: Delta_c1(:) !> must be already allocated
    
    integer(kind=4) :: i,n_p_aq,n_p
    real(kind=8), allocatable :: c1_old(:)
    
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    
    if (n_p/=size(c1)) error stop "Dimension error in update_conc_prim_species"
        
    c1_old=c1 !> old primary concentrations
    do i=1,n_p
        if (c1_old(i)+Delta_c1(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*c1_old(i)) then
            c1(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*c1_old(i)
        else if (c1_old(i)+Delta_c1(i)>=c1_old(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            c1(i)=c1_old(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            c1(i)=c1_old(i)+Delta_c1(i)
        end if
        Delta_c1(i)=c1(i)-c1_old(i)
    end do
    call this%set_conc_prim_species(c1)
end subroutine