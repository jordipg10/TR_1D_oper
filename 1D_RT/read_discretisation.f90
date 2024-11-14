subroutine read_time_discretisation(this,unit,root)
    use RT_1D_m
    class(RT_1D_c) :: this
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: root
    
    character(len=256) :: str,label
    character(len=:), allocatable :: str_trim
    type(mesh_1D_Euler_homog_c) :: mesh
    type(time_discr_homog_c) :: time_discr
    type(BCs_t) :: BCs
    type(tpt_props_heterog_c) :: tpt_props
    integer(kind=4) :: i,j,BCs_label(2),comm_ind,targets_flag
    integer(kind=4), allocatable :: num_mix_waters(:),mix_wat_indices(:,:)
    real(kind=8) :: Delta_x,Delta_t,r,phi,D
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
                read(unit,*) Delta_t
                read(unit,*) int_method_chem_reacts
                call time_discr%set_Num_time(Num_time)
                call time_discr%set_Delta_t_homog(Delta_t)
                call time_discr%compute_Final_time()
                call this%set_int_method_chem_reacts(int_method_chem_reacts)
                call this%transport%set_time_discr(time_discr)
            else
                continue
            end if
        end do
        close(unit)
    end select
end subroutine