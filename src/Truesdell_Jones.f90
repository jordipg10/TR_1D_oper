!> Truesdell-Jones
!! \f$(0.1<I<0.7)\f$
subroutine Truesdell_Jones_sub(this,ionic_act,log_act_coeff)
    use aq_species_m
    implicit none
    class(aq_species_c) :: this
    real(kind=8), intent(in) :: ionic_act !> I
    real(kind=8), intent(out) :: log_act_coeff !> log_10(gamma)
    
    !log_act_coeff=-this%params_act_coeff%A*sqrt(ionic_act)*(this%valence**2)/(1d0+this%params_act_coeff%B*this%params_act_coeff%a_TJ*sqrt(ionic_act))+this%params_act_coeff%b_TJ*ionic_act
end subroutine