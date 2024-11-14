!!> Computes component concentrations after iteration of WMA method for equilibrium chemical system
!!> Computes variable activity concentrations from component concentrations
!!> We assume there are aqueous and solid variable activity species
!subroutine transport_iter_exch_EE_aq_chem(this,this_old,c_tilde,porosity,Delta_t)
!    use aqueous_chemistry_m
!    implicit none
!!> Arguments
!    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
!    class(aqueous_chemistry_c), intent(in) :: this_old !> aqueous chemistry object at previous time step (nombre muy malo)
!    real(kind=8), intent(in) :: c_tilde(:)
!    real(kind=8), intent(in), optional :: porosity
!    real(kind=8), intent(in), optional :: Delta_t !> time step
!!> Variables
!    integer(kind=4) :: niter !> number of iterations in Newton speciation
!    logical :: CV_flag !> convergence flag
!    real(kind=8) :: mu=0d0 !> Newton initialistaion parameter
!    real(kind=8), allocatable :: conc_new(:) !> concentration variable activity species next time step
!    real(kind=8), allocatable :: conc_old(:) !> concentration variable activity species current time step
!    real(kind=8), allocatable :: conc_old_old(:) !> concentration variable activity species previous time step
!    real(kind=8), allocatable :: conc_comp(:) !> concentration components next time step
!!> Pre-process
!    conc_old=this%get_conc_nc()
!    conc_old_old=this_old%get_conc_nc()
!!> Process  
!    !> We compute aqueous component concentrations after mixing
!        conc_comp=MATMUL(THIS%speciation_alg%comp_mat,c_tilde)
!        !conc_comp=this%get_conc_comp_exch()
!    !> Loop until speciation converges
!        do
!        !> We initialise variable activity concentrations for Newton speciation
!            call initialise_iterative_method(conc_old_old,conc_old,mu,conc_new)
!        !> We update variable activity concentrations
!            call this%update_c_nc(conc_new)
!        !> We compute variable activity concentrations from component concentrations
!            call this%compute_c_nc_from_u_Newton(conc_comp,niter,CV_flag)
!        !> We check convergence
!            if (CV_flag==.false.) then !> NO CV
!                if (mu<1d0) then
!                    mu=mu+0.25
!                else
!                    error stop "Newton speciation does not converge"
!                end if
!            else
!                exit
!            end if
!        end do
!!> Post-process
!    !> We check concentratuibns
!        call this%check_conc_aq_var_act_species(conc_comp)
!        call this%check_act_aq_species()
!end subroutine