!> Solves reactive mixing problem. Computes state variables and reaction rates.
subroutine solve_RT_1D(this,root,unit)
    use RT_1D_m
    implicit none
    class(RT_1D_c) :: this
    character(len=*), intent(in) :: root
    integer(kind=4), intent(in) :: unit
    
    select type (this)
    class is (RT_1D_transient_c)
        call this%chemistry%solve_reactive_mixing(root,unit,this%transport%mixing_ratios,this%transport%mixing_waters_indices,this%transport%F_mat%diag,this%transport%time_discr,this%int_method_chem_reacts)
    end select
end subroutine