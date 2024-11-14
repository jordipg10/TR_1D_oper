!> Computes Monod kinetic reaction rate gradient
!!> drk_dc=rk*dlogT/dc
!!>        =rk*k/(c*(k+c)) if electron acceptor/donor
!!>        =-rk/(k+c) if inhibitor
subroutine compute_drk_dc_Monod(this,conc,rk,drk_dc)
    use redox_kin_reaction_m
    implicit none
    class(redox_kin_c), intent(in) :: this
    real(kind=8), intent(in) :: conc(:) !> concentrations of relevant species
    real(kind=8), intent(in) :: rk
    real(kind=8), intent(out) :: drk_dc(:) !> must be allocated
    
    
    integer(kind=4) :: i,j,k,n,m,l,DOC_ind
    type(species_c), allocatable :: stoich_mat_react_zonec_species(:)

    drk_dc=0d0
!> Inhibitors
    if (this%params%n_inh>0) then
        do i=1,this%params%n_inh
            drk_dc(i)=-rk/(this%params%k_inh(i)+conc(i))
        end do
    end if
!> Electron accceptor & donor
    do i=1,2
        drk_dc(this%params%n_inh+i)=rk*this%params%k_M(i)/(conc(this%params%n_inh+i)*(this%params%k_M(i)+conc(this%params%n_inh+i)))
    end do
end subroutine