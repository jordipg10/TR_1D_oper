module analytical_solutions_diffusion_m
    use diffusion_m
    use diffusion_transient_m
    implicit none
    save
    contains
        function fund_sol_diff_eqn_1D(M,Delta_x,x,mu,t,phi,D) result(conc)
        !> Fundamental solution of diffusion equation in 1D
            real(kind=8), intent(in) :: M
            real(kind=8), intent(in) :: Delta_x
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: mu
            real(kind=8), intent(in) :: t
            real(kind=8), intent(in) :: phi
            real(kind=8), intent(in) :: D
            real(kind=8) :: conc
            
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            real(kind=8), parameter :: eps=1d-12
            
            !if (mod(n,2).eqv.0) then
            !    mu=(n-1)/2d0
            !else
            !    mu=floor(n/2d0)
            !end if
            
            if (abs(t)<eps) then
                if (abs(x-mu)<eps) then
                    conc=M/Delta_x
                else
                    conc=0d0
                end if
            else
                conc=(M/sqrt(4d0*pi*D*t))*exp(-(25d-2*(x-mu)**2)/(D*t))
            end if
        end function
        
       
        
end module