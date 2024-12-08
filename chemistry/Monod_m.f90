!> Redox kinetic reaction module
module redox_kin_reaction_m
    use kin_reaction_m
    use Monod_params_m
    use aq_phase_m
    implicit none
    save
    type, public, extends(kin_reaction_c) :: redox_kin_c
        type(Monod_params_c) :: params !> kinetic parameters Monod
    contains
    !> Set
        procedure, public :: set_Monod_params
    !> Compute
        procedure, public :: compute_drk_dc_Monod
        procedure, public :: compute_rk_Monod
    !> Rearrange
        procedure, public :: rearrange_indices_aq_phase_Monod
    !> Check
        !procedure, public :: check_indices_aq_phase_Monod
    !> Write
        procedure, public :: write_params=>write_Monod_params
        procedure, public :: write_reaction=>write_redox_reaction
    end type
    
    !type, public :: monod_type
    !>    PetscInt :: id
    !>    character(len=MAXWORDLENGTH) :: species_name
    !>    PetscReal :: half_saturation_constant
    !>    PetscReal :: threshold_concentration
    !>    type(monod_type), pointer :: next
    !>  end type monod_type
    
    interface
     
        subroutine compute_rk_Monod(this,conc,rk)
            import redox_kin_c
            implicit none
            class(redox_kin_c), intent(in) :: this
            real(kind=8), intent(in) :: conc(:)
            real(kind=8), intent(out) :: rk
        end subroutine

        
        subroutine compute_drk_dc_Monod(this,conc,rk,drk_dc)
            import redox_kin_c
            implicit none
            class(redox_kin_c), intent(in) :: this
            real(kind=8), intent(in) :: conc(:)
            real(kind=8), intent(in) :: rk
            real(kind=8), intent(out) :: drk_dc(:)
        end subroutine
        
      
    end interface
    
    contains

        
        subroutine set_Monod_params(this,Monod_params)
            implicit none
            class(redox_kin_c) :: this
            class(Monod_params_c), intent(in) :: Monod_params
            this%params=Monod_params
        end subroutine
        
        
        
        subroutine rearrange_indices_aq_phase_Monod(this,aq_phase_old,aq_phase_new)
            implicit none
            class(redox_kin_c) :: this
            class(aq_phase_c), intent(in) :: aq_phase_old
            class(aq_phase_c), intent(in) :: aq_phase_new
            
            integer(kind=4) :: i,ind_new
            type(aq_species_c) :: DOC
            logical :: flag
            
            do i=1,this%params%num_terms
                call aq_phase_new%is_species_in_aq_phase(aq_phase_old%aq_species(this%indices_aq_phase(i)),flag,ind_new)
                if (flag==.true.) then
                    this%indices_aq_phase(i)=ind_new
                else
                    error stop "Inhibitor/electron acceptor/electron donor is not in aqueous phase"
                end if
            end do
        end subroutine
        
        subroutine write_Monod_params(this,unit)
            import redox_kin_c
            implicit none
            class(redox_kin_c) :: this
            integer(kind=4), intent(in) :: unit !> file unit
            write(unit,*) this%params%k_inh
        end subroutine
        
        subroutine write_redox_reaction(this,unit)
            import redox_kin_c
            implicit none
            class(redox_kin_c) :: this
            integer(kind=4), intent(in) :: unit !> file unit
            call write_reaction_sup(this,unit)
            call this%write_params(unit)
        end subroutine
end module