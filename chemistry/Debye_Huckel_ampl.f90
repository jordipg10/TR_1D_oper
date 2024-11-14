!> Debye-Huckel extended
!> Diluted solutions $(0.01 \le I \le 0.1)$
subroutine Debye_Huckel_ampl(this,ionic_act,log_act_coeff)
    use aq_species_m
    implicit none
    class(aq_species_c) :: this
    real(kind=8), intent(in) :: ionic_act !> ionic activity
    real(kind=8), intent(out) :: log_act_coeff !> log_10(gamma)
    
    !log_act_coeff=-(this%params_act_coeff%A*sqrt(ionic_act)*this%valence**2)/(1d0+this%params_act_coeff%ion_size_param*this%params_act_coeff%B*sqrt(ionic_act))
end subroutine