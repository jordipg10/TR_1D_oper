!> Kinetic reaction module
module kin_reaction_m
    use reaction_m
    use kin_params_m
    implicit none
    save
    type, public, abstract, extends(reaction_c) :: kin_reaction_c !> Kinetic reaction abstract subclass
        integer(kind=4), allocatable :: indices_aq_phase(:) !> indices of species in aqueous phase object that control kinetic reaction rates
    contains
        procedure(write_params), public, deferred :: write_params
    end type
    
    type, public :: kin_reaction_poly_c !< clase ad hoc para crear vector punteros clase reaccion cinetica (cosas de Fortran)
        class(kin_reaction_c), pointer :: kin_reaction
    contains
        procedure, public :: set_kin_reaction
    end type
    
    abstract interface
        subroutine write_params(this,unit)
            import kin_reaction_c
            implicit none
            class(kin_reaction_c) :: this
            integer(kind=4), intent(in) :: unit !> file unit
        end subroutine
        
    end interface
    
    interface 
        
        
        
    end interface
    
    contains
        subroutine set_kin_reaction(this,kin_reaction)
            implicit none
            class(kin_reaction_poly_c) :: this
            class(kin_reaction_c), intent(in), target :: kin_reaction
            this%kin_reaction=>kin_reaction
        end subroutine
            
end module