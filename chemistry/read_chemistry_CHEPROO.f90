!> Lectura quimica CHEPROO
subroutine read_chemistry_CHEPROO(this,root,path_DB,unit_chem_syst_file,unit_loc_chem_file,unit_target_waters_init_file,unit_output_file)
    use chemistry_Lagr_m
    implicit none
    class(chemistry_c) :: this
    character(len=*), intent(in) :: root
    character(len=*), intent(in) :: path_DB
    integer(kind=4), intent(in) :: unit_chem_syst_file
    !character(len=*), intent(in) :: chem_syst_file
    integer(kind=4), intent(in) :: unit_loc_chem_file
    !character(len=*), intent(in) :: loc_chem_file
    integer(kind=4), intent(in) :: unit_target_waters_init_file
    !character(len=*), intent(in) :: target_waters_init_file
    integer(kind=4), intent(in) :: unit_output_file
    !character(len=*), intent(in) :: output_file
    
    real(kind=8), allocatable :: conc_exch(:)
    integer(kind=4) :: unit,i,j,n_sp,num_minerals,k,n_paths,niter,niwtype, nbwtype, nrwtype,index,kin_react_type,nk,chem_syst_type,nischem
    integer(kind=4), allocatable :: n_tar_sol(:),n_tar_aq(:),ind_wat_type(:),num_aq_prim_array(:),num_cstr_array(:)
    real(kind=8) :: aux,conc,temp
    character(len=256) :: str,str1,str2,str3,str4,str5,Monod_name,label
    logical :: flag,CV_flag
    real(kind=8), parameter :: tolerance=1d-14, rel_tolerance=1d-9, control_factor=1d-1, epsilon=1d-16
    integer(kind=4), parameter :: niter_max=40
    
    type(water_type_c), allocatable :: init_wat_types(:),bd_wat_types(:),rech_wat_types(:)
    type(solid_type_c), allocatable :: init_min_zones(:),init_cat_exch_zones_bis(:),init_cat_exch_zones(:),bd_cat_exch_zones(:),init_sol_zones(:)
    type(gas_type_c), allocatable :: init_gas_zones(:)
    
    character(len=256), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:),solid_species_str(:),kin_react_names(:)
    type(species_c) :: species,water,constrain
    type(species_c), allocatable :: aq_species(:),cst_act_species(:),prim_species(:)
    type(mineral_c) :: mineral
    class(mineral_c), allocatable :: minerals(:)
    type(reactive_zone_c), allocatable :: reactive_zones(:),rz_mins(:),rz_surf(:)
    class(kin_params_c), pointer :: p_kin_params=>null()
    !type(Monod_params_c), allocatable, target :: Monod_params_array(:)
    !type(Monod_params_c), target :: Monod_params
    class(kin_reaction_c), pointer :: p_kin_react=>null()
    type(kin_reaction_ptr_c) :: kin_react_ptr
    class(kin_reaction_ptr_c), allocatable :: kin_reacts(:)
    !type(redox_kin_c), target :: Monod
    type(eq_reaction_c) :: eq_react
    type(eq_reaction_c), allocatable :: eq_reacts(:)
    !type(kin_lin_c), target :: linear
    !type(kin_mineral_c), target :: kin_mineral
    
    class(chem_system_c), pointer :: p_chem_syst=>null()
    
    type(aq_phase_c), target :: aq_phase_new
    
    !type(surface_c), allocatable :: init_cat_exch_zones(:)
    

    
!> Chemical system
    open(unit_chem_syst_file,file=root//'_sist_quim.dat',status='old',action='read')
    call this%chem_syst%read_chem_system_CHEPROO(path_DB,unit_chem_syst_file)
    close(unit_chem_syst_file)
!> Local chemistry
    open(unit_loc_chem_file,file=root//'_quim_loc.dat',status='old',action='read')
    do
        read(unit_loc_chem_file,*) label
        if (label=='end') then
            rewind(unit_loc_chem_file)
            exit
        !else if (label=='INITIAL AND BOUNDARY WATER TYPES') then
            !call this%read_init_bd_rech_wat_types_CHEPROO(unit_loc_chem_file,ind_wat_type,num_aq_prim_array,num_cstr_array)
        !else if (label=='INITIAL MINERAL ZONES') then
        !    call this%read_init_min_zones_CHEPROO(unit_loc_chem_file,init_min_zones,rz_mins)
        else if (label=='INITIAL SURFACE ADSORPTION ZONES') then
            call this%read_init_cat_exch_zones_CHEPROO(unit_loc_chem_file,init_cat_exch_zones_bis,reactive_zones)
        else if (label=='INITIAL GAS ZONES') then
            call this%read_init_gas_zones_CHEPROO(unit_loc_chem_file,init_gas_zones)
        else
            continue
        end if
    end do
    do
        read(unit_loc_chem_file,*) label
        if (label=='end') then
            exit
        else if (label=='INITIAL AND BOUNDARY WATER TYPES') then
            if (size(init_gas_zones)==1) then
                call this%read_init_bd_rech_wat_types_CHEPROO(unit_loc_chem_file,ind_wat_type,num_aq_prim_array,num_cstr_array,init_cat_exch_zones_bis,init_gas_zones(1)%gas_chem)
            else if (size(init_gas_zones)==0) then
                call this%read_init_bd_rech_wat_types_CHEPROO(unit_loc_chem_file,ind_wat_type,num_aq_prim_array,num_cstr_array,init_cat_exch_zones_bis)
            end if
        else if (label=='INITIAL MINERAL ZONES') then
            call this%read_init_min_zones_CHEPROO(unit_loc_chem_file,init_min_zones,reactive_zones)
        !else if (label=='INITIAL GAS ZONES') then
            !call this%read_init_gas_zones_CHEPROO(unit_loc_chem_file,init_gas_zones,reactive_zones)
        else
            continue
        end if
    end do
    close(unit_loc_chem_file)
    allocate(init_sol_zones(size(init_min_zones)+size(init_cat_exch_zones_bis)))
    do i=1,size(init_min_zones)
        init_sol_zones(i)=init_min_zones(i)
    end do
    do i=1,size(init_cat_exch_zones_bis)
        init_sol_zones(size(init_min_zones)+i)=init_cat_exch_zones_bis(i)
    end do
    do i=1,size(reactive_zones)
        call reactive_zones(i)%set_num_solids()
        call reactive_zones(i)%set_non_flowing_species()
        call reactive_zones(i)%set_eq_reactions()
        call reactive_zones(i)%set_stoich_mat_react_zone()
        call reactive_zones(i)%set_stoich_mat_sol_rz()
        !call reactive_zones(i)%chem_syst%set_stoich_mat_gas()
    end do
!> We verify pointers
    if (size(reactive_zones)>0) then
        do i=1,size(init_sol_zones)
            if (.not. associated(init_sol_zones(i)%solid_chem%reactive_zone%chem_syst)) then
                call init_sol_zones(i)%solid_chem%reactive_zone%set_chem_syst_react_zone(this%chem_syst)
                !print *, associated(init_sol_zones(i)%solid_chem%reactive_zone%chem_syst)
            end if
            call init_sol_zones(i)%solid_chem%set_reactive_zone(reactive_zones(i))
        end do
        do i=1,size(init_gas_zones)
            if (.not. associated(init_gas_zones(i)%gas_chem%reactive_zone%chem_syst)) then
                call init_gas_zones(i)%gas_chem%reactive_zone%set_chem_syst_react_zone(this%chem_syst)
                !print *, associated(init_sol_zones(i)%solid_chem%reactive_zone%chem_syst)
            end if
            call init_gas_zones(i)%gas_chem%set_reactive_zone(reactive_zones(i))
        end do
    end if
!> Target waters
    open(unit_target_waters_init_file,file=root//'_tar_wat.dat',status='old',action='read')
    call this%read_target_waters_init(unit_target_waters_init_file,this%wat_types,init_sol_zones,init_gas_zones,niter,CV_flag)
    close(unit_target_waters_init_file)
!> Output data
    call this%chem_out_options%read_chem_out_options(root,unit_output_file,this%target_waters)
    if (this%chem_syst%num_minerals>0 .or. this%chem_syst%gas_phase%num_species>0) then
        call this%set_reactive_zones()
    end if
end subroutine