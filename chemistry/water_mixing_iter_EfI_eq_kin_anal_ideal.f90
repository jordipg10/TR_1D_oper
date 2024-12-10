!> Computes aqueous species concentrations after iteration of WMA using Euler fully implicit in chemical reactions for kinetic system
!! We assume all primary species are aqueous
!! The Jacobians are computed analytically
subroutine water_mixing_iter_EfI_eq_kin_anal_ideal(this,c1_old,c_tilde,conc_nc,porosity,Delta_t)
    use aqueous_chemistry_m
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    !class(aqueous_chemistry_c), intent(in) :: this_old !> aqueous chemistry object at previous time step (nombre muy malo)
    real(kind=8), intent(in) :: c1_old(:)
    !real(kind=8), intent(in) :: c2nc_ig(:)
    real(kind=8), intent(in) :: c_tilde(:)
    real(kind=8), intent(out) :: conc_nc(:)
    !real(kind=8), intent(out) :: conc_comp(:)
    real(kind=8), intent(in), optional :: porosity !> at this target
    real(kind=8), intent(inout), optional :: Delta_t !> time step
!> Variables
    real(kind=8), allocatable :: c1(:) !> aqueous primary concentrations
    integer(kind=4) :: n_p !> number of primary species
    integer(kind=4) :: k_div !> counter of time step divisions
    integer(kind=4) :: k !> counter of time steps
    integer(kind=4) :: niter !> number of Newton iterations
    real(kind=8) :: mu !> Newton initialisation parameter
    logical :: CV_flag !> convergence flag
!> Pre-process
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    !> We assign aqueous primary species concentrations
        c1=this%concentrations(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
    !> Newton initialisation parameter     
        mu=0d0 
!> Process
    !> Initial guess aquoeus primary concentrations
        call initialise_iterative_method(c1_old,c1,mu,this%concentrations(1:n_p))
        k=0
        k_div=0
        do
            do !> loop until convergence is reached
            !> We apply Newton method to compute aqueous concentrations
                call this%Newton_EfI_rk_eq_kin_aq_anal_ideal(c_tilde,porosity,Delta_t,conc_nc,niter,CV_flag)
            !> We check convergence
                if (CV_flag==.false.) then !> no CV
                    if (mu<1d0) then
                        mu=mu+0.25 !> we increase Newton initialisation parameter
                        call initialise_iterative_method(c1_old,c1,mu,this%concentrations(1:n_p))
                    else
                        mu=0d0
                        Delta_t=Delta_t/2d0
                        k_div=k_div+1
                        !print *, k_div
                        !error stop "Iterative method does not CV"
                    end if
                else 
                    k=k+1
                    exit
                end if
            end do
            if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
                exit
            end if
        end do
    !> We compute component concentrations
        !conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,conc_nc)
!< Post-process
    deallocate(c1)
end subroutine