!> Gaines-Thomas convention module
module Gaines_Thomas_m
    use exch_sites_conv_m
    implicit none
    save
    type, public, extends(exch_sites_conv_c) :: Gaines_Thomas_c
        
    contains
        procedure, public :: compute_log_act_coeff_ads_cat=>compute_log_act_coeff_Gaines_Thomas
    end type
    
    contains
        !> Computes logarithm activity coefficient of surface complex using Gaines-Thomas convention
        !!> Activity of surface complex = equivalent fraction
        subroutine compute_log_act_coeff_Gaines_Thomas(this,exchangeable_cation,CEC,log_act_coeff)
            implicit none
            class(Gaines_Thomas_c) :: this
            class(species_c), intent(in) :: exchangeable_cation !> exchangeable cation
            real(kind=8), intent(in) :: CEC !> cation exchange capacity
            real(kind=8), intent(out) :: log_act_coeff !> log_10 activity coefficient surface complex
        
            log_act_coeff=log10(exchangeable_cation%valence*1d0)-log10(CEC)
        end subroutine
        

    
    
      
end module