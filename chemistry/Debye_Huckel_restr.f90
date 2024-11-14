!> Debye-Huckel restricted
!> Very diluted solutions (\f$ I<0.01 \f$)
subroutine Debye_Huckel_restr(this,ionic_act,log_act_coeff)
    use aq_species_m
    implicit none
    class(aq_species_c), intent(in) :: this
    real(kind=8), intent(in) :: ionic_act !> ionic activity
    real(kind=8), intent(out) :: log_act_coeff !> log_10(gamma)
    
    !log_act_coeff=-this%params_act_coeff%A*sqrt(ionic_act)*this%valence**2
end subroutine