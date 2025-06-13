!> Writes transport data from 1D reactive transport problem solved with the WMA
!> Output file is already opened
subroutine write_transport_data_WMA(this,unit)
    use transport_transient_m, only: transport_1D_transient_c
    implicit none
    class(transport_1D_transient_c), intent(in) :: this                 !> 1D transient transport object                         
    integer(kind=4), intent(in) :: unit                                 !> unit of output file
    
    integer(kind=4) :: i,j,num_cells

    num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag !> number of cells
    write(unit,"(2x,'Number of targets:',I5/)") this%spatial_discr%Num_targets
    write(unit,"(/,2x,'Time step:'/)")
    write(unit,"(2x,ES15.5/)") this%time_discr%get_Delta_t()
    write(unit,"(2x,'Final time:'/)")
    write(unit,"(2x,ES15.5/)") this%time_discr%Final_time
    write(unit,"(/,2x,'Dimension + Mixing ratios concentrations (by rows) (including boundary and sink/source terms):'/)")
    do i=1,this%mixing_ratios_conc%num_cols
        write(unit,"(2x,I5,*(ES15.5))") this%mixing_ratios_conc%cols(i)%dim, (this%mixing_ratios_conc%cols(i)%col_1(j), &
            j=1,this%mixing_ratios_conc%cols(i)%dim)
    end do
    write(unit,"(/,2x,'Initial mixing ratios reaction amounts (by rows):'/)")
    do i=1,this%mixing_ratios_Rk%num_cols
        write(unit,"(2x,*(ES15.5))") (this%mixing_ratios_Rk_init%cols(i)%col_1(j), &
            j=1,this%mixing_ratios_Rk%cols(i)%dim)
    end do
    write(unit,"(/,2x,'Final mixing ratios reaction amounts (by rows):'/)")
    do i=1,this%mixing_ratios_Rk%num_cols
        write(unit,"(2x,(ES15.5))") this%mixing_ratios_Rk%cols(i)%col_1(1)
        !write(unit,"(2x,I5,*(ES15.5))") this%mixing_ratios_Rk%cols(i)%dim, (this%mixing_ratios_Rk%cols(i)%col_1(j), &
        !    j=1,this%mixing_ratios_Rk%cols(i)%dim)
    end do
    write(unit,"(/,2x,'Mixing waters indices + no. upstream waters + no. downstream waters:'/)")
    do i=1,this%mixing_waters_indices%num_cols
        write(unit,"(2x,*(I5))") (this%mixing_waters_indices%cols(i)%col_1(j), &
            j=1,this%mixing_waters_indices%cols(i)%dim)
    end do
    write(unit,"(/,2x,'Reaction mixing waters indices + no. upstream waters + no. downstream waters:'/)")
    do i=1,this%mixing_waters_indices_dom%num_cols
        write(unit,"(2x,*(I5))") (this%mixing_waters_indices_dom%cols(i)%col_1(j), &
            j=1,this%mixing_waters_indices_dom%cols(i)%dim)
    end do
end subroutine write_transport_data_WMA