module analytical_solutions_transport_m
    use analytical_solutions_diffusion_m
    use transport_transient_m
    use transport_m
    implicit none
    save
    contains
        function fund_sol_tpt_eqn_1D(c0,M,Delta_x,x,mu,t,phi,D,q) result(conc)
            implicit none
            real(kind=8), intent(in) :: c0
            real(kind=8), intent(in) :: M
            real(kind=8), intent(in) :: Delta_x
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: mu
            real(kind=8), intent(in) :: t
            real(kind=8), intent(in) :: phi
            real(kind=8), intent(in) :: D
            real(kind=8), intent(in) :: q
            real(kind=8) :: conc
            
            real(kind=8), parameter :: eps=1d-12

                conc=exp(5d-1*(q/D)*(x-mu-t*q/(2*phi)))
                conc=conc*fund_sol_diff_eqn_1D(M,Delta_x,x,mu,t,phi,D)+c0
        end function
        
        function anal_sol_tpt_1D_stat_flujo_lin(this,x) result(conc)
        !> Analytical solution 1D stationary transport equation with:
        !>   D=cst
        !>   q(x)=-x
        !>   c(0)=1
        !>   c(L)=0
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: x
            real(kind=8) :: conc
            
            real(kind=8) :: D
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            
            D=this%tpt_props_heterog%dispersion(1)

            if (this%dimensionless.eqv..true.) then
                conc=-erf(x/sqrt(2d0))/erf(1d0/sqrt(2d0)) + 1
            else
                conc=-erf(x/sqrt(2d0*D))/erf(this%spatial_discr%measure/sqrt(2d0*D)) + 1d0
            end if
        end function
        
   
        function der_anal_sol_tpt_1D_stat_flujo_lin(this,x) result(der_conc)
        !> Derivative of analytical solution 1D stationary transport equation with:
        !>   D=cst
        !>   q(x)=-x
        !>   c(0)=1
        !>   c(L)=0
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: x
            real(kind=8) :: der_conc
            
            real(kind=8) :: D,L 
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            
            L=this%spatial_discr%measure !> 1D
            D=this%tpt_props_heterog%dispersion(1)
            der_conc=(-sqrt(2d0)/(erf(this%spatial_discr%measure/sqrt(2d0*D))*sqrt(pi*D)))*exp(-(x**2)/(2d0*D)) !> derivative analytical solution
        end function

        
    
        function der_anal_sol_tpt_1D_stat_flujo_cuad(this,x) result(der_conc)
        !> Derivative of analytical solution 1D stationary transport equation with:
        !>   D=cst
        !>   q(x)=-x^2
        !>   c(0)=1
        !>   c(L)=0
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: x
            real(kind=8) :: der_conc
            
            real(kind=8) :: D,L
            real(kind=8), parameter :: pi=4d0*atan(1d0),incompl_gamma_term=0.327336564991358
            
            D=this%tpt_props_heterog%dispersion(1)
            L=this%spatial_discr%measure !> 1D

            der_conc=-exp(-(x**3)/(3d0*D))*(3d0**(2d0/3d0))/(gamma(1d0/3d0)-incompl_gamma_term)
        end function
end module