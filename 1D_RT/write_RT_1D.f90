!> Writes data and results of 1D reactive transport model
subroutine write_RT_1D(this,unit,root,path_py)
    use RT_1D_m
    implicit none
    class(RT_1D_c), intent(in) :: this
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: root
    character(len=*), intent(in), optional :: path_py
    
    open(unit,file=root//'.out',status='unknown',form='formatted')
    call this%write_transport_data(unit)
    call this%chemistry%write_chemistry(unit)
    close(unit)
    if (present(path_py)) then
        call this%write_python(path_py)
    end if
end subroutine