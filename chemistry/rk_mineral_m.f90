!module rk_mineral_m
!    use rk_m
!    use kin_mineral_m
!    implicit none
!    save
!    type, public, extends(rk_c) :: rk_mineral_c
!        type(kin_mineral_params_c) :: params
!        !real(kind=8) :: sigma
!        !!integer(kind=4) :: zeta
!        !real(kind=8) :: act_E ! E_a
!        !integer(kind=4) :: N
!        !real(kind=8), allocatable :: k(:)
!        !real(kind=8), allocatable :: a(:)
!        !real(kind=8), allocatable :: p(:,:)
!        real(kind=8) :: saturation ! Omega
!        !real(kind=8), allocatable :: theta(:)
!        !real(kind=8), allocatable :: eta(:)
!        !real(kind=8) :: temp
!    contains
!        procedure, public :: read_rk=>read_mineral
!        !procedure, public :: compute_far_from_eq_fct
!        procedure, public :: compute_rk=>compute_mineral
!        procedure, public :: compute_drk_dc=>compute_drk_dc_mineral
!        procedure, public :: set_parameters
!    end type
!        
!    interface
!        subroutine compute_mineral(this,conc,i)
!            import rk_mineral_c
!            implicit none
!            class(rk_mineral_c) :: this
!            real(kind=8), intent(in) :: conc(:)
!            integer(kind=4), intent(in), optional :: i
!        end subroutine
!        
!        subroutine compute_drk_dc_mineral(this,chem_syst,drk_dc,conc,i)
!            import rk_mineral_c
!            import chem_system_c
!            implicit none
!            class(rk_mineral_c) :: this
!            class(chem_system_c), intent(in) :: chem_syst
!            real(kind=8), intent(out) :: drk_dc(:)
!            real(kind=8), intent(in), optional :: conc(:)
!            integer(kind=4), intent(in), optional :: i
!        end subroutine
!    end interface
!    
!    contains
!        subroutine read_mineral(this,filename)
!            implicit none
!            class(rk_mineral_c) :: this
!            character(len=*), intent(in) :: filename
!            open(unit=1,file=filename,status='old',action='read')
!            read(1,"(/,I10)") this%params%temp, this%params%N, this%params%sigma, this%params%act_E, this%params%k, this%params%a, this%params%theta, this%params%eta
!            !read(1,*) this%sigma
!            !read(1,*) this%zeta
!            !read(1,*) this%act_E
!            !read(1,*) this%saturation
!            !read(1,*) this%k
!            !read(1,*) this%a
!            !!read(1,*) this%p_k
!            !read(1,*) this%theta
!            !read(1,*) this%eta
!            !read(1,*) p_k
!            close(1)
!        end subroutine
!        
!        !function compute_far_from_eq_fct(this) result(far_from_eq)
!        !    implicit none
!        !    class(rk_mineral_c) :: this
!        !    real(kind=8), allocatable :: far_from_eq(:)
!        !    far_from_eq=this%saturation**this%theta-1d0
!        !end function
!        
!        subroutine set_parameters(this,kin_mineral)
!            implicit none
!            class(rk_mineral_c) :: this
!            class(kin_mineral_c), intent(in) :: kin_mineral
!            this%params%sigma=kin_mineral%params%sigma
!            !this%zeta=kin_mineral%zeta
!            this%params%act_E=kin_mineral%params%act_E
!            this%params%N=kin_mineral%params%N
!            this%params%k=kin_mineral%params%k
!            this%params%a=kin_mineral%params%a
!            this%params%p=kin_mineral%params%p
!            this%saturation=kin_mineral%mineral%saturation
!            this%params%theta=kin_mineral%params%theta
!            this%params%eta=kin_mineral%params%eta
!            this%params%temp=kin_mineral%params%temp
!        end subroutine
!end module