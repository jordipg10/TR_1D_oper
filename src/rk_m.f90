!> !> Kinetic reaction rate module
!module rk_m
!>    use reaction_rate_m
!>    use kin_params_m
!>    use chem_system_m
!>    implicit none
!>    save
!>    type, public, extends(reaction_rate_c) :: rk_c
!>    contains
!>        procedure, public :: compute_rk
!>        procedure, public :: compute_drk_dc
!>        !procedure(compute_rk), public, deferred :: compute_rk
!>        !procedure(compute_drk_dc), public, deferred :: compute_drk_dc
!>        !procedure(read_rk), public, deferred :: read_rk
!>        !procedure(compute_drk_dc2), public, deferred :: compute_drk_dc2
!>    end type
!>    
!>    interface 
!>        subroutine compute_rk(this,conc,kin_params,i,saturation)
!>            import rk_c
!>            import kin_params_c
!>            implicit none
!>            class(rk_c) :: this
!>            real(kind=8), intent(in) :: conc(:)
!>            class(kin_params_c), intent(in) :: kin_params
!>            integer(kind=4), intent(in) :: i
!>            real(kind=8), intent(in), optional :: saturation
!>        end subroutine
!>        
!>        subroutine compute_drk_dc(this,chem_syst,kin_params,i,drk_dc,conc)
!>            import rk_c
!>            import chem_system_c
!>            import kin_params_c
!>            implicit none
!>            class(rk_c) :: this
!>            class(chem_system_c), intent(in) :: chem_syst
!>            class(kin_params_c), intent(in) :: kin_params
!>            integer(kind=4), intent(in) :: i
!>            real(kind=8), intent(out) :: drk_dc(:)
!>            real(kind=8), intent(in), optional :: conc(:)
!>        end subroutine
!>        
!>        !subroutine read_rk(this,filename)
!>        !>    import rk_c
!>        !>    implicit none
!>        !>    class(rk_c) :: this
!>        !>    character(len=*), intent(in) :: filename
!>        !end subroutine
!>    end interface
!end module