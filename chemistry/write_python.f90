!> Writes data and results necessary for Python unit testing
subroutine write_python(this,path)
    use RT_1D_m
    
    implicit none
    class(RT_1D_c), intent(in) :: this !> chemistry object
    character(len=*), intent(in) :: path !> path output files
    
    integer(kind=4) :: i,j,n,l,k
    integer(kind=4), allocatable :: tar_sol_indices(:),tar_wat_indices(:)
        
    select type(this)
    type is (RT_1D_transient_c)
    !> Writes initial aqueous concentrations
        open(999,file=trim(path)//'conc_aq_init.dat')
        do i=1,this%chemistry%chem_syst%aq_phase%num_species
            write(999,"(*(ES15.5))") (this%chemistry%target_waters_init(j)%concentrations(i), j=1,this%chemistry%num_target_waters)
        end do
        close(999)
    !> Writes initial aqueous activity coefficients
        open(9999,file=trim(path)//'gamma_aq_init.dat')
        do i=1,this%chemistry%chem_syst%aq_phase%num_species
            write(9999,"(*(ES15.5))") (10**(this%chemistry%target_waters_init(j)%log_act_coeffs(i)), j=1,this%chemistry%num_target_waters)
        end do
        close(9999)
    !> Writes mixing ratios
        open(998,file=trim(path)//'lambdas.dat')
        if (this%transport%time_discr%int_method==1) then
            do i=1,this%transport%mixing_ratios%num_cols
                write(998,"(*(F15.5))") (this%transport%mixing_ratios%cols(i)%col_1(j), j=1,this%transport%mixing_ratios%cols(i)%dim)
            end do
        else if (this%transport%time_discr%int_method==2) then
            do i=1,this%transport%mixing_ratios%num_cols
                write(998,"(*(F15.5))") (this%transport%mixing_ratios_mat(i,j), j=1,this%transport%mixing_ratios%num_cols)
            end do
        end if
        close(998)
    end select
end subroutine