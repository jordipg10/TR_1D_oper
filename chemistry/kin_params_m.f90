!> Kinetic parameters structure
module kin_params_m
    implicit none
    save
    type, public, abstract :: kin_params_c !> kinetic parameters abstract superclass (structure)
        real(kind=8) :: rate_cst !> reaction rate constant
    end type
end module