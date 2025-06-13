subroutine solve_write_PDE_1D(this,Time_out)
!> Calls subroutines that perform computations and writes results
    use transport_transient_m, only: PDE_1D_c
    implicit none
    !> Variables
    class(PDE_1D_c) :: this !> PDE object
    real(kind=8), intent(in) :: Time_out(:)
    
    integer(kind=4) :: Num_output
    real(kind=8), allocatable :: conc_out(:,:)
    
    Num_output=size(Time_out)
    allocate(conc_out(this%spatial_discr%Num_targets,Num_output))
    call this%solve_PDE_1D(Time_out,conc_out)
    call this%write_PDE_1D(Time_out,conc_out)
end subroutine
    
    


