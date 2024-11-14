!!> Computes component concentrations after iteration of WMA method for equilibrium chemical system
!subroutine transport_iter_EI_aq_chem(this,mixing_ratios,mixing_waters,rk_mat,porosity,Delta_t)
!    use aqueous_chemistry_m
!    
!    implicit none
!    class(aqueous_chemistry_c) :: this
!    real(kind=8), intent(in) :: mixing_ratios(:,:)
!    class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)
!    real(kind=8), intent(in), optional :: rk_mat(:,:) !> rk matrix
!    real(kind=8), intent(in), optional :: porosity
!    real(kind=8), intent(in), optional :: Delta_t !> time step
!    
!    integer(kind=4) :: niter
!    logical :: CV_flag
!    real(kind=8), allocatable :: conc_comp(:)
!    
!    conc_comp=this%compute_u_tilde_aq_chem(mixing_ratios(:,1),mixing_waters)
!    call this%compute_c_nc_from_u_aq_Newton(niter,CV_flag) !> rezaei
!    call this%check_conc_aq_var_act_species(conc_comp)
!    call this%check_act_aq_species()
!end subroutine