!> Performs Newton method in reactive mixing iteration using Euler fully implicit in chemical reactions
!! The chemical system has equilibrium and kinetic reactions

!> Computes concentration of variable activity species in target j at time step k+1
    
!> We assume all primary species are aqueous
    
!> The Jacobians are computed analytically

!! We DO NOT apply lumping in this subroutine
    
!! We use option 2 for the time integration of kinetic reactions
    
    
subroutine Newton_EI_rk_eq_kin_aq_anal_ideal_opt2(this,u_tilde,u_rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: u_tilde(:) !> component concentrations after mixing
    !real(kind=8), intent(in) :: porosity !> in target j
    !real(kind=8), intent(in) :: rk_tilde_old(:) !> old kinetic reaction rates after mixing (de momento no se usa)
    real(kind=8), intent(in) :: u_rk_tilde(:) !> reaction part of components concentrations after mixing 
    !real(kind=8), intent(in) :: rk_tilde_new(:) !> new kinetic reaction rates after mixing (de momento no se usa)
    real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction amount in this target
    real(kind=8), intent(in) :: Delta_t !> (k+1)-th time step
    real(kind=8), intent(in) :: theta !> reaction time weighting factor
    real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
!> Variables
    real(kind=8), allocatable :: rk(:) !> kinetic reaction rates
    real(kind=8), allocatable :: drk_dc(:,:) !> Jacobian of kinetic reactions
    real(kind=8), allocatable :: dfk_dc1(:,:) !> Jacobian Newton residual
    real(kind=8), allocatable :: c2nc(:) !> secondary variable activity concentrations next iteration
    real(kind=8), allocatable :: Delta_c1(:) !> primary concentrations difference
    real(kind=8), allocatable :: rk_old(:) !> old kinetic reaction rates
    real(kind=8), allocatable :: rk_avg(:) !> averaged kinetic reaction rates
    real(kind=8), allocatable :: fk(:) !> Newton residual
    integer(kind=4) :: n_p !> number of primary species
    integer(kind=4) :: n_nc !> number of variable activity species
    logical :: zero_flag !> flag for zero concentration
!> Pre-process
    niter=0 !> we initialise number of iterations
    CV_flag=.false. !> we initialise convergence flag
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_nc=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    allocate(c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions),dfk_dc1(n_p,n_p),Delta_c1(n_p),&
        drk_dc(this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts,n_nc),&
        rk(this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts),fk(n_p))
    drk_dc=0d0 !> we initialise Jacobian of kinetic reactions
    dfk_dc1=0d0 !> we initialise Jacobian of Newton residual
!> Process
    !> We compute component concentrations after mixing
        !u_tilde=this%compute_u_tilde(c_tilde)
    !> We compute reaction part of components after mixing (llamalo de otra forma)
        !u_rk_tilde=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,Rk_tilde)
    !> We get old kinetic reaction rates
        rk_old=this%get_rk()
    !> Newton loop
        do 
            niter=niter+1 !> we update number of iterations
            if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then
                print *, inf_norm_vec_real(fk)
                print *, "Too many iterations in subroutine Newton_EI_rk_eq_kin_aq_anal_ideal_opt2"
                exit
            end if
            conc_nc(1:n_p)=this%get_c1() !> we update primary species concentrations
        !> Compute concentration of secondary variable activity species from concentration of primary species using mass action law
            call this%compute_c2nc_from_c1_aq_ideal(c2nc)
            conc_nc(n_p+1:n_nc)=c2nc !> chapuza
        !> We compute kinetic reaction rates and its Jacobian analitically
            call this%compute_rk_Jac_rk_anal(rk,drk_dc)
        !> We compute the ponderated average of kinetic reaction rates
            rk_avg=theta*rk+(1d0-theta)*rk_old
        !> Newton residual
            fk=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,conc_nc)-u_tilde-u_rk_tilde-Delta_t*mix_ratio_Rk*&
                matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk_avg)
        !> Check convergence
            if (inf_norm_vec_real(fk)<this%solid_chemistry%reactive_zone%CV_params%abs_tol) then !> CV reached
                CV_flag=.true.
                call this%compute_Rk_mean(theta,Delta_t) !> we update mean reaction amounts
                exit
            else
                call this%compute_dfk_dc1_aq_EfI_ideal(conc_nc(n_p+1:n_nc),drk_dc,Delta_t,theta,mix_ratio_Rk,dfk_dc1) !> computes Jacobian of Newton resiudal
                call LU_lin_syst(dfk_dc1,-fk,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c1) !> solves linear system dfk_dc1*Delta_c1=-fk, where c1_new=c1_old+Delta_c1
                !call Gauss_Jordan(dfk_dc1,-fk,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c1) !> solves linear system dfk_dc1*Delta_c1=-fk, where c1_new=c1_old+Delta_c1
                if (inf_norm_vec_real(Delta_c1/conc_nc(1:n_p))<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2)&
                    then !> we check relative tolerance
                    print *, "Newton solution not accurate enough"
                    exit
                else
                    call this%update_conc_aq_prim_species(Delta_c1,zero_flag) !> updates c1 and Delta_c1
                    if (zero_flag) then
                        return
                    end if
                end if
            end if
        end do
!> Post-process
    deallocate(fk,dfk_dc1,Delta_c1,drk_dc,rk,c2nc)
end subroutine