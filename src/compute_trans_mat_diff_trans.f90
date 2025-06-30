subroutine compute_trans_mat_diff_trans(this)
!> T: transition matrix (tridiagonal, symmetric, negative stoich_mat_react_zonemi-definite)
!> rows sum = 0 if r=0
!> F*dc/dt=T*c+g
    use diffusion_transient_m
    use spatial_discr_rad_m
    implicit none
    
    class(diffusion_1D_transient_c) :: this
    
    real(kind=8) :: r_i,r_i_12
    integer(kind=4) :: i,n,opcion
    
    n=this%spatial_discr%Num_targets

    if (this%spatial_discr%adapt_ref.eq.1) then
        deallocate(this%trans_mat%sub,this%trans_mat%diag,this%trans_mat%super)
        call this%allocate_trans_mat()
    end if
    select type (mesh=>this%spatial_discr)
    type is (spatial_discr_rad_c)
        if (this%dimless.eqv..true.) then
            do i=1,n-1
                r_i_12=sum(mesh%Delta_r(1:i))
                this%trans_mat%sub(i)=(r_i_12**(mesh%dim-1))/(0.5*(mesh%Delta_r(i)+mesh%Delta_r(i+1)))
            end do
        end if
    end select
    this%trans_mat%super=this%trans_mat%sub
    this%trans_mat%diag=-this%diff%diff_props_heterog%source_term_flag*this%diff%diff_props_heterog%source_term
    this%trans_mat%diag(2:n-1)=this%trans_mat%diag(2:n-1)-this%trans_mat%sub(1:n-2)-this%trans_mat%super(2:n-1)
end subroutine 