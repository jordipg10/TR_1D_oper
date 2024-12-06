!> This subroutine reads the time discretisation used for the WMA from a file ended in "_WMA_discr.dat"
subroutine read_time_discretisation(this,unit,root)
    use RT_1D_m
    class(RT_1D_c) :: this                                      !> transient 1D reactive transport object
    integer(kind=4), intent(in) :: unit                         !> unit of file
    character(len=*), intent(in) :: root                        !> root of file
    
    integer(KIND=4) :: Num_time                                 !< number of time steps
    integer(KIND=4) :: int_method_chem_reacts                   !< integration method chemical reactions
    character(len=256) :: label
    type(time_discr_homog_c), target :: time_discr_homog
    type(time_discr_heterog_c), target :: time_discr_heterog
    class(time_discr_c), pointer :: time_discr=>NULL()
    real(kind=8), allocatable :: Delta_t(:)                     !< time step array
    
    select type (this)
    type is (RT_1D_transient_c)
    !> We open file with time discretisation data
        open(unit,file=root//'_WMA_discr.dat',status='old',action='read')
        do 
            read(unit,*) label
            if (label=='end') then
                exit
            else if (label=='TIME') then
            !> we read number of time steps
                read(unit,*) Num_time
                allocate(Delta_t(Num_time))
            !< we read time step or steps
                read(unit,*) Delta_t
                if (inf_norm_vec_real(Delta_t-Delta_t(1))<EPSILON) then !< if time step is uniform, then we use homogeneous subclass
                    call time_discr_homog%set_Delta_t_homog(Delta_t(1))
                    time_discr=>time_discr_homog
                else
                    call time_discr_heterog%set_Delta_t_heterog(Delta_t) !< if time step is not uniform, then we use heterogeneous subclass
                    time_discr=>time_discr_heterog
                end if
            !< we read time integration method for chemical reactions
                read(unit,*) int_method_chem_reacts
                call this%set_int_method_chem_reacts(int_method_chem_reacts)
            !< we set the number of time steps attribute
                call time_discr%set_Num_time(Num_time)
            !< we compute the final time
                call time_discr%compute_Final_time()
                call this%transport%set_time_discr(time_discr)
            else
                continue
            end if
        end do
        close(unit)
    end select
end subroutine