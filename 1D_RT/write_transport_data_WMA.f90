!> Writes transport data from 1D reactive transport problem solved with the WMA
!> Output file is already opened
subroutine write_transport_data_WMA(this,unit)
    use transport_transient_m
    implicit none
    class(transport_1D_transient_c), intent(in) :: this                 !> 1D transient transport object                         
    integer(kind=4), intent(in) :: unit                                 !> unit of output file
    
    integer(kind=4) :: i,j,num_cells

    num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag !> number of cells
    !if (this%int_method_chem_reacts==1) then
    !    write(unit,"(2x,'Integration method chemical reactions:',10x,'Euler explicit',/)")
    !else if (this%int_method_chem_reacts==2) then
    !    write(unit,"(2x,'Integration method chemical reactions:',10x,'Euler fully implicit',/)")
    !else
    !    error stop "Integration method not implemented yet for RT"
    !end if
    write(unit,"(2x,'Number of targets:',I5/)") this%spatial_discr%Num_targets
    write(unit,"(/,2x,'Time step:'/)")
    write(unit,"(2x,ES15.5/)") this%time_discr%get_Delta_t()
    write(unit,"(2x,'Final time:'/)")
    write(unit,"(2x,ES15.5/)") this%time_discr%Final_time
    write(unit,"(/,2x,'Dimension + Mixing ratios (by rows) (including boundary and sink/source terms):'/)")
    do i=1,this%mixing_ratios%num_cols
        write(unit,"(2x,I5,*(ES15.5))") this%mixing_ratios%cols(i)%dim, (this%mixing_ratios%cols(i)%col_1(j), j=1,this%mixing_ratios%cols(i)%dim)
    end do
    write(unit,"(/,2x,'Mixing waters indices:'/)")
    do i=1,this%mixing_waters_indices%num_cols
        write(unit,"(2x,I5,*(I5))") this%mixing_waters_indices%cols(i)%col_1
    end do
end subroutine