!!> Surface complexation module
!module surf_compl_m
!>    use phase_m
!>    use aq_species_m
!>    implicit none
!>    save
!>    type, public, extends(phase_c) :: surface_c
!>        integer(kind=4) :: num_surf_compl
!>        type(species_c), allocatable :: surf_compl(:) !> surface complexes
!>    contains
!>        procedure, public :: set_num_surf_compl
!>        procedure, public :: set_surf_compl
!>        procedure, public :: allocate_surf_compl
!>        procedure, public :: compute_log_act_coeffs_surf_compl
!>        procedure, public :: compute_log_Jacobian_act_coeffs_surf_compl
!>    end type
!>    
!>    interface
!>        subroutine compute_log_act_coeffs_surf_compl(this,log_act_coeffs)
!>            import surface_c
!>            implicit none
!>            class(surface_c) :: this
!>            real(kind=8), intent(out) :: log_act_coeffs(:) !> must be allocated
!>        end subroutine
!>        
!>        subroutine compute_log_Jacobian_act_coeffs_surf_compl(this,log_act_coeffs,conc,log_Jacobian_act_coeffs)
!>            import surface_c
!>            implicit none
!>            class(surface_c) :: this
!>            real(kind=8), intent(in) :: log_act_coeffs(:)
!>            real(kind=8), intent(in) :: conc(:) !> concentration of species in a given target
!>            real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !> must be allocated
!>        end subroutine
!>    end interface
!>    
!>    contains
!>    
!>        subroutine set_num_surf_compl(this,num_surf_compl)
!>            implicit none
!>            class(surface_c) :: this
!>            integer(kind=4), intent(in), optional :: num_surf_compl
!>            if (present(num_surf_compl)) then
!>                this%num_surf_compl=num_surf_compl
!>            else
!>                this%num_surf_compl=size(this%surf_compl)
!>            end if
!>        end subroutine
!>        
!>        subroutine set_surf_compl(this,surf_compl)
!>            implicit none
!>            class(surface_c) :: this
!>            class(species_c), intent(in) :: surf_compl(:)
!>            
!>            if (allocated(this%surf_compl) .and. size(surf_compl)>this%num_surf_compl) then
!>                error stop "Number of surface complexes cannot be greater than number of species"
!>            else
!>                this%surf_compl=surf_compl
!>            end if
!>        end subroutine
!>        
!>        subroutine allocate_surf_compl(this,num_surf_compl)
!>            implicit none
!>            class(surface_c) :: this
!>            integer(kind=4), intent(in), optional :: num_surf_compl
!>            
!>            if (present(num_surf_compl) .and. num_surf_compl>=0) then
!>                allocate(this%surf_compl(num_surf_compl))
!>            else if (present(num_surf_compl)) then
!>                error stop "Number of surface complexes must be non-negative"
!>            else
!>                allocate(this%surf_compl(this%num_surf_compl))
!>            end if
!>        end subroutine
!end module