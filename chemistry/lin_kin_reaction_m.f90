!> Linear kinetics module
!!> $A \to B$
!!> \f$r_k(c)=\lambda*c_\f$
module lin_kin_reaction_m
    use kin_reaction_m
    use aq_phase_m
    implicit none
    save
    type, public, extends(kin_reaction_c) :: lin_kin_reaction_c !> linear kinetic reaction subclass
        real(kind=8) :: lambda !> reaction rate
    contains
        procedure, public :: set_lambda
        procedure, public :: compute_rk_lin
        procedure, public :: compute_drk_dc_lin
        procedure, public :: append_lin_kin_reaction
        procedure, public :: set_index_aq_phase_lin
    end type
    
    interface
        subroutine compute_rk_lin(this,conc,rk)
            import lin_kin_reaction_c
            implicit none
            class(lin_kin_reaction_c), intent(in) :: this
            real(kind=8), intent(in) :: conc
            real(kind=8), intent(out) :: rk
        end subroutine
        
        subroutine compute_drk_dc_lin(this,drk_dc)
            import lin_kin_reaction_c
            implicit none
            class(lin_kin_reaction_c), intent(in) :: this
            real(kind=8), intent(out) :: drk_dc(:)
        end subroutine
        
        subroutine read_kin_lin(this)
            import lin_kin_reaction_c
            implicit none
            class(lin_kin_reaction_c) :: this
        end subroutine
        
        subroutine get_conc_kin_lin(this,species,conc,conc_kin,kin_ind)
            import lin_kin_reaction_c
            import species_c
            implicit none
            class(lin_kin_reaction_c), intent(in) :: this
            class(species_c), intent(in) :: species(:)
            real(kind=8), intent(in) :: conc(:) !> species concentrations
            real(kind=8), intent(out) :: conc_kin(:) !> concentration of species relevant to kinetic reaction rates
            integer(kind=4), intent(out), optional :: kin_ind(:) !> indices of species relevant to kinetic reaction rates
        end subroutine
        
    end interface
    
    contains
        subroutine set_lambda(this,lambda)
            import lin_kin_reaction_c
            implicit none
            class(lin_kin_reaction_c) :: this
            real(kind=8), intent(in) :: lambda
            this%lambda=lambda
        end subroutine
        
        subroutine append_lin_kin_reaction(this,kin_reactions)
            implicit none
            class(lin_kin_reaction_c), intent(in) :: this
            type(lin_kin_reaction_c), intent(inout), allocatable :: kin_reactions(:)
        
            type(lin_kin_reaction_c), allocatable :: aux(:)
        
            aux=kin_reactions
            deallocate(kin_reactions)
            allocate(kin_reactions(size(aux)+1))
            kin_reactions(1:size(aux))=aux
            kin_reactions(size(kin_reactions))=this
        end subroutine
        
        subroutine set_index_aq_phase_lin(this,aq_phase)
            implicit none
            class(lin_kin_reaction_c) :: this
            class(aq_phase_c), intent(in) :: aq_phase
            
            integer(kind=4) :: aq_species_ind
            logical :: flag
            
            allocate(THIS%indices_aq_phase(1))
            call aq_phase%is_species_in_aq_phase(this%species(1),flag,aq_species_ind)
            if (flag==.true.) then
                this%indices_aq_phase(1)=aq_species_ind
            else
                error stop "Linear species is not in aqueous phase"
            end if
        end subroutine
end module