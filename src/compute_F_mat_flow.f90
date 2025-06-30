subroutine compute_F_mat_flow(this)
use flow_transient_m, only: flow_transient_c, spatial_discr_rad_c
    implicit none
    class(flow_transient_c) :: this
    
    integer(kind=4) :: i,n
    real(kind=8) :: r_i

    !> Compute the F matrix for the flow transient object
    !> This subroutine computes the F matrix for the flow transient object based on its properties

    n=this%spatial_discr%Num_targets
    if (this%spatial_discr%adapt_ref==1) then
        deallocate(this%F_mat%diag)
        call this%F_mat%allocate_array(n)
    end if
    select type (mesh=>this%spatial_discr)
    type is (spatial_discr_rad_c)
        if (this%dimless .eqv. .true.) then
            do i=1,n
                r_i=mesh%r_min+sum(mesh%Delta_r(1:i))-5d-1*mesh%Delta_r(i)
                this%F_mat%diag(i)=mesh%Delta_r(i)*(r_i**(mesh%dim-1))
            end do
        else
            do i=1,n
                r_i=mesh%r_min+sum(mesh%Delta_r(1:i))-5d-1*mesh%Delta_r(i)
                this%F_mat%diag(i)=mesh%Delta_r(i)*(r_i**(mesh%dim-1))*this%flow_props_heterog%storativity(i)
            end do
        end if
    end select
    
end subroutine compute_F_mat_flow