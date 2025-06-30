!> Lectura quimica CHEPROO
subroutine read_chemistry_CHEPROO(this,root,path_DB,unit_chem_syst_file,unit_loc_chem_file,unit_target_waters_init_file,&
        unit_output_file)
    use chemistry_Lagr_m, only: chemistry_c, aqueous_chemistry_c, solid_chemistry_c, gas_chemistry_c
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
    
    integer(kind=4) :: i,ngrz,nmrz,nsrz
    integer(kind=4), allocatable :: ind_wat_type(:),num_aq_prim_array(:),num_cstr_array(:)
    character(len=256) :: label
    
    type(aqueous_chemistry_c), allocatable :: wat_types(:),bd_wat_types(:),rech_wat_types(:)
    type(solid_chemistry_c), allocatable :: init_min_zones(:),init_cat_exch_zones(:),init_sol_zones(:)
    type(gas_chemistry_c), allocatable :: init_gas_zones(:)
    
    !type(species_c) :: species,water,constrain
    !type(species_c), allocatable :: aq_species(:),cst_act_species(:),prim_species(:)
    !type(mineral_c) :: mineral
    !class(mineral_c), allocatable :: minerals(:)
    !type(reactive_zone_c), allocatable :: reactive_zones(:),rz_mins(:),rz_surf(:)
    !class(kin_params_c), pointer :: p_kin_params=>null()
    !type(Monod_params_c), allocatable, target :: Monod_params_array(:)
    !type(Monod_params_c), target :: Monod_params
    !class(kin_reaction_c), pointer :: p_kin_react=>null()
    !type(kin_reaction_poly_c) :: kin_react_ptr
    !class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
    !type(redox_kin_c), target :: Monod
    !type(eq_reaction_c) :: eq_react
    !type(eq_reaction_c), allocatable :: eq_reacts(:)
    logical :: flag_comp
    !type(kin_mineral_c), target :: kin_mineral
    
    !class(chem_system_c), pointer :: p_chem_syst=>null()
    
    !type(aq_phase_c), target :: aq_phase_new
    
    !type(surface_c), allocatable :: init_cat_exch_zones(:)
    
    flag_comp=.true.
    ngrz=0 !> number of gas reactive zones by default
    nmrz=0 !> number of mineral reactive zones by default
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
            call this%read_init_cat_exch_zones_CHEPROO(unit_loc_chem_file,init_cat_exch_zones)
            nsrz=size(init_cat_exch_zones)
            !print *, init_cat_exch_zones(1)%reactive_zone%num_non_flowing_species
        else if (label=='INITIAL AND BOUNDARY GAS ZONES') then
            call this%read_init_gas_zones_CHEPROO(unit_loc_chem_file,init_gas_zones,ngrz)
        else
            continue
        end if
    end do
!> Chapuza
    !print*, size(init_gas_zones)
    !print *, allocated(init_gas_zones)
    if (allocated(init_gas_zones)) then
        do i=1,size(init_gas_zones) 
            call init_gas_zones(i)%set_reactive_zone(this%reactive_zones(i))
        end do
    else if (size(init_gas_zones)>0) then !> Chapuza
        allocate(init_gas_zones(0))
        deallocate(init_gas_zones)
        !call init_gas_zones(1)%set_reactive_zone(this%reactive_zones(1))
    end if
    !print *, size(init_cat_exch_zones)
    !print *, allocated(init_cat_exch_zones)
    if (allocated(init_cat_exch_zones)) then
        !print *, init_cat_exch_zones(1)%reactive_zone%num_non_flowing_species
        do i=1,size(init_cat_exch_zones)
            call init_cat_exch_zones(i)%set_reactive_zone(this%reactive_zones(i))
        end do
    else if (size(init_cat_exch_zones)>0) then !> Chapuza
        allocate(init_cat_exch_zones(0))
        deallocate(init_cat_exch_zones)
        !call init_cat_exch_zones(1)%set_reactive_zone(this%reactive_zones(size(init_gas_zones)+1))
    end if
    do
        read(unit_loc_chem_file,*) label
        if (label=='end') then
            exit
        else if (label=='INITIAL AND BOUNDARY WATER TYPES') then
            if (size(init_gas_zones)==1) then !> Chapuza
                call this%read_init_bd_wat_types_CHEPROO(unit_loc_chem_file,&
                init_cat_exch_zones,init_gas_zones(1))
            else
                call this%read_init_bd_wat_types_CHEPROO(unit_loc_chem_file,&
                init_cat_exch_zones)
            end if
        else if (label=='INITIAL MINERAL ZONES') then
            if (size(init_cat_exch_zones)==1) then
                call this%read_init_min_zones_CHEPROO(unit_loc_chem_file,init_min_zones,nmrz,init_cat_exch_zones(1))
            else
                call this%read_init_min_zones_CHEPROO(unit_loc_chem_file,init_min_zones,nmrz)
            end if
            !> Chapuza
            if (allocated(init_gas_zones)) then
                do i=1,size(init_gas_zones)
                    call init_gas_zones(i)%set_reactive_zone(this%reactive_zones(i))
                end do
            end if
            !if (allocated(init_cat_exch_zones)) then
            !    do i=1,size(init_cat_exch_zones)
            !        call init_cat_exch_zones(i)%set_reactive_zone(this%reactive_zones(size(init_gas_zones)+i))
            !    end do
            !end if
        else
            continue
        end if
    end do
    close(unit_loc_chem_file)
    !do i=1,size(init_cat_exch_zones)
    !    call init_cat_exch_zones(i)%set_reactive_zone(this%reactive_zones(size(init_gas_zones)+1))
    !end do
!!> Chapuza
!    if (size(init_cat_exch_zones)>0) then
!        allocate(init_cat_exch_zones(this%num_init_wat_types))
!        if (size(init_cat_exch_zones)==1) then
!            do i=1,this%num_init_wat_types
!                init_cat_exch_zones(ind_wat_type(i))=init_cat_exch_zones(1)
!                call this%init_wat_types(ind_wat_type(i))%read_wat_type_CHEPROO(num_aq_prim_array(ind_wat_type(i)),num_cstr_array(ind_wat_type(i)),this%model,this%Jac_opt,unit,niter,CV_flag,init_cat_exch_zones(ind_wat_type(i))%solid_chem)
!                !init_cat_exch_zones(ind_wat_type(i))%solid_chem%concentrations(1)=conc_exch(num_aq_prim_array(ind_wat_type(i))+1)
!                !init_cat_exch_zones(ind_wat_type(i))%solid_chem%concentrations(2:init_cat_exch_zones(ind_wat_type(i))%solid_chem%reactive_zone%cat_exch_zone%num_surf_compl)=conc_exch(this%init_wat_types(ind_wat_type(i))%aq_chem%aq_phase%num_species+2:this%init_wat_types(ind_wat_type(i))%aq_chem%solid_chemistry%reactive_zone%speciation_alg%num_species)
!            end do
!            do i=1,this%num_bd_wat_types
!                init_cat_exch_zones(ind_wat_type(this%num_init_wat_types+i))=init_cat_exch_zones(1)
!                call this%bd_wat_types(ind_wat_type(this%num_init_wat_types+i))%read_wat_type_CHEPROO(num_aq_prim_array(ind_wat_type(this%num_init_wat_types+i)),num_cstr_array(ind_wat_type(this%num_init_wat_types+i)),this%model,this%Jac_opt,unit,niter,CV_flag,init_cat_exch_zones(ind_wat_type( i))%solid_chem)
!                !init_cat_exch_zones(ind_wat_type(this%num_init_wat_types+i))%solid_chem%concentrations(1)=conc_exch(num_aq_prim_array(ind_wat_type(i))+1)
!                !init_cat_exch_zones(ind_wat_type(i))%solid_chem%concentrations(2:init_cat_exch_zones(ind_wat_type(i))%solid_chem%reactive_zone%cat_exch_zone%num_surf_compl)=conc_exch(this%init_wat_types(ind_wat_type(i))%aq_chem%aq_phase%num_species+2:this%init_wat_types(ind_wat_type(i))%aq_chem%solid_chemistry%reactive_zone%speciation_alg%num_species)
!            end do
!        end if
!    else
!        do i=1,this%num_init_wat_types
!            call this%init_wat_types(ind_wat_type(i))%read_wat_type_CHEPROO(num_aq_prim_array(ind_wat_type(i)),num_cstr_array(ind_wat_type(i)),this%model,this%Jac_opt,unit,niter,CV_flag)
!        end do
!        do i=1,this%num_bd_wat_types
!            call this%bd_wat_types(ind_wat_type(this%num_init_wat_types+i))%read_wat_type_CHEPROO(num_aq_prim_array(ind_wat_type(this%num_init_wat_types+i)),num_cstr_array(ind_wat_type(this%num_init_wat_types+i)),this%model,this%Jac_opt,unit,niter,CV_flag)
!        end do
!    end if
!> Chapuza 
    !if (allocated(init_cat_exch_zones) .and. allocated(init_min_zones)) then
    !    allocate(init_sol_zones(size(init_min_zones)+size(init_cat_exch_zones)))
    !    do i=1,size(init_min_zones)
    !        init_sol_zones(size(init_cat_exch_zones)+i)=init_min_zones(i)
    !    end do
    !    do i=1,size(init_cat_exch_zones)
    !        init_sol_zones(i)=init_cat_exch_zones(i)
    !    end do
    if (allocated(init_cat_exch_zones) .and. (.not. allocated(init_min_zones))) then
       allocate(init_sol_zones(size(init_cat_exch_zones)))
       do i=1,size(init_cat_exch_zones)
           init_sol_zones(i)=init_cat_exch_zones(i)
       end do
       !nsrz=size(init_cat_exch_zones)
    else if (allocated(init_min_zones) .and. (.not. allocated(init_cat_exch_zones))) then
        allocate(init_sol_zones(size(init_min_zones)))
        do i=1,size(init_min_zones)
            init_sol_zones(i)=init_min_zones(i)
        end do
        nsrz=nmrz
    else if (allocated(init_cat_exch_zones) .and. allocated(init_min_zones)) then
        allocate(init_sol_zones(size(init_min_zones)))
        do i=1,size(init_min_zones)
            init_sol_zones(i)=init_min_zones(i)
        end do
        nsrz=nmrz
    else
        allocate(init_sol_zones(0)) !> Chapuza
        nsrz=0
    end if  
    
    
!!> Chapuza
!    do i=1,size(reactive_zones)
!        call reactive_zones(i)%set_num_solids()
!        call reactive_zones(i)%set_non_flowing_species()
!        call reactive_zones(i)%set_speciation_alg_dimensions(flag_comp)
!        call reactive_zones(i)%set_eq_reactions()
!        call reactive_zones(i)%set_stoich_mat_react_zone()
!        call reactive_zones(i)%set_stoich_mat_sol_rz()
!        call reactive_zones(i)%set_CV_params(this%CV_params) !> chapuza
!    end do
!> We verify pointers (chapuza)
    !if (size(reactive_zones)>0) then
    !    do i=1,size(init_sol_zones)
    !        if (.not. associated(init_sol_zones(i)%reactive_zone%chem_syst)) then
    !            call init_sol_zones(i)%reactive_zone%set_chem_syst_react_zone(this%chem_syst)
    !            !print *, associated(init_sol_zones(i)%solid_chem%reactive_zone%chem_syst)
    !        end if
    !        call init_sol_zones(i)%set_reactive_zone(reactive_zones(i))
    !    end do
    !    do i=1,size(init_gas_zones)
    !        if (.not. associated(init_gas_zones(i)%reactive_zone%chem_syst)) then
    !            call init_gas_zones(i)%reactive_zone%set_chem_syst_react_zone(this%chem_syst)
    !            !print *, associated(init_sol_zones(i)%solid_chem%reactive_zone%chem_syst)
    !        end if
    !        call init_gas_zones(i)%set_reactive_zone(reactive_zones(i))
    !    end do
    !end if
!> We set reactive zones attribute
    !call this%set_reactive_zones(reactive_zones)
!> Target waters
    call this%read_target_waters_init(root,init_sol_zones,init_gas_zones,nsrz,ngrz)
!> Output data
    call this%chem_out_options%read_chem_out_options(root,unit_output_file,this%target_waters)
end subroutine