subroutine compute_source_term_PDE(this)
!> $g=E*c_r$
    use PDE_m, only: PDE_1D_c
    use transport_m, only: transport_1D_c, diffusion_1D_c
    use transport_transient_m, only: transport_1D_transient_c, diffusion_1D_transient_c
    implicit none
    class(PDE_1D_c) :: this
    
    !allocate(this%source_term_PDE(this%spatial_discr%Num_targets))
    !this%source_term_PDE=0d0 !> $g=0$ chapuza
    select type (this)
    class is (diffusion_1D_c)
        this%source_term_PDE=this%ext_mat%diag*this%conc_ext
    class is (diffusion_1D_transient_c)
        this%source_term_PDE=this%ext_mat%diag*this%conc_ext
    end select
    this%source_term_PDE(1)=this%bd_mat(1)*this%BCs%conc_inf
    this%source_term_PDE(this%spatial_discr%Num_targets)=this%bd_mat(2)*this%BCs%conc_out
    ! select type (this)
    ! type is (transport_1D_c)
    !     this%source_term_PDE=this%conc_r_flag*this%tpt_props_heterog%source_term*this%conc_ext
    ! type is (transport_1D_transient_c)
    !     this%source_term_PDE=this%conc_r_flag*this%tpt_props_heterog%source_term*this%conc_ext
    ! type is (diffusion_1D_c)
    !     this%source_term_PDE=this%diff_props_heterog%source_term*this%conc_ext
    ! type is (diffusion_1D_transient_c)
    !     this%source_term_PDE=this%diff_props_heterog%source_term*this%conc_ext
    ! end select

end subroutine