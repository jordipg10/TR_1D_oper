!> This subroutine reads the time discretisation used for the WMA from a file ended in "_WMA_discr.dat"
subroutine read_time_discretisation(this,unit,root)
    use RT_1D_m
    class(RT_1D_c) :: this                                      !> transient 1D reactive transport object
    integer(kind=4), intent(in) :: unit                         !> unit of file
    character(len=*), intent(in) :: root                        !> root of file
    
    character(len=256) :: str,label
    character(len=:), allocatable :: str_trim
    type(time_discr_homog_c), target :: time_discr_homog
    type(time_discr_heterog_c), target :: time_discr_heterog
    class(time_discr_c), pointer :: time_discr=>NULL()
    type(BCs_t) :: BCs
    type(tpt_props_heterog_c) :: tpt_props
    integer(kind=4) :: i,j,BCs_label(2),comm_ind,targets_flag
    integer(kind=4), allocatable :: num_mix_waters(:),mix_wat_indices(:,:)
    real(kind=8), allocatable :: Delta_t(:)
    real(kind=8), parameter :: epsilon=1d-16
    logical :: flag_props
    type(diag_matrix_c) :: F_mat
    
    
    
    
    select type (this)
    type is (RT_1D_transient_c)
        open(unit,file=root//'_WMA_discr.dat',status='old',action='read')
        do 
            read(unit,*) label
            if (label=='end') then
                exit
            else if (label=='TIME') then
                read(unit,*) Num_time
                allocate(Delta_t(Num_time))
                read(unit,*) Delta_t
                if (inf_norm_vec_real(Delta_t-Delta_t(1))<EPSILON) then
                    call time_discr_homog%set_Delta_t_homog(Delta_t(1))
                    time_discr=>time_discr_homog
                else
                    call time_discr_heterog%set_Delta_t_heterog(Delta_t)
                    time_discr=>time_discr_heterog
                end if
                read(unit,*) int_method_chem_reacts
                call this%set_int_method_chem_reacts(int_method_chem_reacts)
                call time_discr%set_Num_time(Num_time)
                call time_discr%compute_Final_time()
                call this%transport%set_time_discr(time_discr)
            else
                continue
            end if
        end do
        close(unit)
    end select
end subroutine