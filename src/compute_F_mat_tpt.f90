subroutine compute_F_mat_tpt(this) !> diagonal matrix
!> F_ii=phi_i
    use transport_transient_m, only: transport_1D_transient_c
    implicit none
    class(transport_1D_transient_c) :: this
    
    integer(kind=4) :: n
    
    n=this%spatial_discr%Num_targets
     
    if (this%spatial_discr%adapt_ref==1) then
        deallocate(this%F_mat%diag)
        call this%F_mat%allocate_array(n)
    end if
    
    if (this%tpt_props_heterog%homog_flag .eqv. .true.) then
        this%F_mat%diag=this%tpt_props_heterog%porosity(1)
    else
        this%F_mat%diag=this%tpt_props_heterog%porosity
    end if
    ! print *, this%tpt_props_heterog%porosity
    ! print *, this%F_mat%diag
end subroutine