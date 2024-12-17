subroutine read_init_min_zones_CHEPROO(this,unit,init_min_zones)
    use chemistry_Lagr_m
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: unit !> file
    type(solid_chemistry_c), intent(out), allocatable :: init_min_zones(:)
    !type(reactive_zone_c), intent(inout), allocatable, optional :: this%reactive_zones(:)
    
    integer(kind=4) :: num_surf_rz,num_gas_rz,num_rz_old,i,j,k,nmtype,nrwtype,icon,num_min_zones,num_mins_rz,num_mins_glob,num_mins_loc,num_mins_loc_eq,ind_rz,imtype,num_rz,min_ind,num_mins_var_eq,num_mins_cst_eq
    integer(kind=4), allocatable :: init_min_zones_indices(:),min_indices(:,:)
    character(len=256) :: str,constrain,label,min_name
    real(kind=8) :: guess,c_tot,temp,vol_frac,react_surf,conc
    logical :: min_flag,flag_comp
    type(mineral_c) :: mineral
    type(reactive_zone_c) :: aux_react_zone
    type(reactive_zone_c), allocatable :: react_zones(:),aux_react_zones(:)
    

    flag_comp=.true.

    read(unit,*) nmtype
    
    allocate(react_zones(nmtype)) !> chapuza
    
    allocate(min_indices(this%chem_syst%num_minerals,nmtype))
    
    !call react_zone%set_chem_syst_react_zone(this%chem_syst)
    !call react_zone%allocate_minerals_react_zone(this%chem_syst%num_minerals_eq)
    !!call react_zone%allocate_non_flowing_species()
    !call react_zone%set_num_solids()
    !call react_zone%set_minerals_react_zone()
    !call react_zone%set_non_flowing_species()
    
    if (nmtype>0) then
        num_rz=nmtype
        allocate(init_min_zones(nmtype))
        num_min_zones=0 !> counter mineral zones
        do
            read(unit,*) imtype, temp !> index of mineral zone, temperature (C)
            if (imtype<1 .or. imtype>nmtype) error stop
                read(unit,*) str
                if (str=='*') then
                    if (num_min_zones<nmtype) then
                        continue
                    else
                        exit
                    end if
                else if (index(str,'mineral')/=0) then
                    num_min_zones=num_min_zones+1
                    num_mins_loc=0 !> counter minerals in this zone
                    num_mins_loc_eq=0 !> counter minerals in eq in this zone
                    num_mins_cst_eq=0
                    num_mins_var_eq=0
                    do
                        read(unit,*) mineral%name
                        if (mineral%name=='*') then
                            call react_zones(imtype)%set_chem_syst_react_zone(this%chem_syst)
                            call react_zones(imtype)%set_CV_params(this%CV_params)
                            call react_zones(imtype)%allocate_minerals_react_zone(num_mins_loc_eq)
                            call react_zones(imtype)%allocate_non_flowing_species()
                            !call react_zones(imtype)%set_non_flowing_species()
                            call react_zones(imtype)%set_num_solids()
                            call react_zones(imtype)%set_num_mins_cst_act(num_mins_cst_eq)
                            call react_zones(imtype)%set_num_mins_var_act(num_mins_var_eq)
                            call react_zones(imtype)%set_speciation_alg_dimensions(flag_comp)
                            !> aqui habria que comprobar si se estan repitiendo zonas reactivas
                            call init_min_zones(imtype)%set_reactive_zone(react_zones(imtype))
                            call init_min_zones(imtype)%allocate_vol_fracts()
                            call init_min_zones(imtype)%allocate_react_surfaces()
                            call init_min_zones(imtype)%allocate_conc_solids()
                            call init_min_zones(imtype)%allocate_log_act_coeffs_solid_chem()
                            call init_min_zones(imtype)%allocate_activities()
                            call init_min_zones(imtype)%set_temp(temp+273.18) !> Kelvin
                            call init_min_zones(imtype)%allocate_var_act_species_indices(react_zones(imtype)%num_minerals_var_act)
                            call init_min_zones(imtype)%allocate_cst_act_species_indices(react_zones(imtype)%num_minerals_cst_act)
                            call init_min_zones(imtype)%set_indices_solids()
                            exit
                        else
                            call this%chem_syst%is_mineral_in_chem_syst(mineral,min_flag,min_ind)
                            if (min_flag==.true.) then
                                num_mins_loc=num_mins_loc+1
                                min_indices(num_mins_loc,num_min_zones)=min_ind
                                if (this%chem_syst%minerals(min_ind)%mineral%cst_act_flag==.true. .and. min_ind>this%chem_syst%num_min_kin_reacts) then
                                    num_mins_cst_eq=num_mins_cst_eq+1
                                    num_mins_loc_eq=num_mins_loc_eq+1
                                else if (this%chem_syst%minerals(min_ind)%mineral%cst_act_flag==.false. .and. min_ind>this%chem_syst%num_min_kin_reacts) then
                                    num_mins_var_eq=num_mins_var_eq+1
                                    num_mins_loc_eq=num_mins_loc_eq+1
                                else
                                    continue
                                end if
                            else
                                error stop
                            end if
                        end if
                    end do
                else
                    exit
                end if
            !end do
            if (num_min_zones==nmtype) exit
        end do
        rewind(unit)
    !> Second iteration
        num_min_zones=0 !> counter mineral zones
        do
            read(unit,*) label
            if (label=='INITIAL MINERAL ZONES') then
                read(unit,*) nmtype !> number of mineral zones
                do
                    read(unit,*) imtype, temp
                    if (imtype<1 .or. imtype>nmtype) error stop
                    read(unit,*) str
                    !print *, init_min_zones(imtype)%reactive_zone%num_minerals
                    if (index(str,'mineral')/=0) then
                        num_min_zones=num_min_zones+1
                        num_mins_loc=0 !> counter minerals in this zone
                        do
                            read(unit,*) mineral%name, vol_frac, react_surf !> react_surf: [m^2 min./m^3 rock]
                            if (mineral%name=='*') then
                                !call init_min_zones(imtype)%react_zones(imtype)%set_eq_reactions()
                                !call init_min_zones(imtype)%react_zones(imtype)%set_stoich_mat_react_zone()
                                exit
                            else
                                num_mins_loc=num_mins_loc+1
                                if (min_indices(num_mins_loc,num_min_zones)>this%chem_syst%num_min_kin_reacts) then
                                    init_min_zones(imtype)%reactive_zone%minerals(num_mins_loc)=this%chem_syst%minerals(min_indices(num_mins_loc,num_min_zones)) !> we set mineral
                                end if
                                init_min_zones(imtype)%vol_fracts(num_mins_loc)=vol_frac !> we set volumetric fraction
                                init_min_zones(imtype)%react_surfaces(num_mins_loc)=react_surf !> we set reactive surface
                                init_min_zones(imtype)%concentrations(num_mins_loc)=1d0 !> we assume minerals are pure phases
                                init_min_zones(imtype)%activities(num_mins_loc)=1d0 !> we assume minerals are pure phases
                            end if
                        end do
                        call init_min_zones(imtype)%reactive_zone%set_non_flowing_species()
                        call init_min_zones(imtype)%reactive_zone%set_eq_reactions()
                        call init_min_zones(imtype)%reactive_zone%set_stoich_mat_react_zone()
                        call init_min_zones(imtype)%reactive_zone%set_stoich_mat_sol_rz()
                        if (num_min_zones==nmtype) then
                            exit
                        end if
                    else
                        exit
                    end if
                end do
                exit
            else
                continue
            end if
        end do
        !if (present(this%reactive_zones)) then
            if (allocated(this%reactive_zones)) then
                num_rz_old=size(this%reactive_zones)
                allocate(aux_react_zones(num_rz_old))
                num_gas_rz=0
                do i=1,num_rz_old
                    call aux_react_zones(i)%assign_react_zone(this%reactive_zones(i))
                    if (this%reactive_zones(i)%gas_phase%num_gases_eq>0 .and. this%reactive_zones(i)%num_solids==0) then
                        num_gas_rz=num_gas_rz+1
                    end if
                end do
                num_surf_rz=(num_rz_old-num_gas_rz)/(1+num_gas_rz)
                num_rz=num_rz_old+nmtype*(1+num_gas_rz)
                !deallocate(this%reactive_zones)
                call this%allocate_reactive_zones(num_rz)
                do i=1,num_gas_rz+num_surf_rz
                    call this%reactive_zones(i)%assign_react_zone(aux_react_zones(i))
                end do
                do i=1,nmtype
                    call this%reactive_zones(i)%assign_react_zone(init_min_zones(i)%reactive_zone)
                    !call this%reactive_zones(num_gas_rz+num_surf_rz+i)%set_chem_syst_react_zone(this%chem_syst)
                    !call this%reactive_zones(num_gas_rz+num_surf_rz+i)%set_CV_params(this%CV_params)
                    !call this%reactive_zones(num_gas_rz+num_surf_rz+i)%allocate_minerals_react_zone(init_min_zones(i)%reactive_zone%num_minerals)
                    !this%reactive_zones(num_gas_rz+num_surf_rz+i)%minerals=init_min_zones(i)%reactive_zone%minerals
                    !call this%reactive_zones(num_gas_rz+num_surf_rz+i)%set_num_mins_cst_act(init_min_zones(i)%reactive_zone%num_minerals_cst_act)
                    !call this%reactive_zones(num_gas_rz+num_surf_rz+i)%set_num_mins_var_act(init_min_zones(i)%reactive_zone%num_minerals_var_act)
                end do
                do i=1,num_gas_rz
                    do j=1,num_surf_rz
                        call this%reactive_zones(num_gas_rz+nmtype+i*num_surf_rz+j)%assign_react_zone(aux_react_zones(num_gas_rz+num_surf_rz+(i-1)*num_surf_rz+j))
                    end do
                    do j=1,nmtype
                        !call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_chem_syst_react_zone(this%chem_syst)
                        call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%assign_react_zone(init_min_zones(i)%reactive_zone)
                        call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_gas_phase(this%reactive_zones(i)%gas_phase)
                        call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_non_flowing_species()
                        call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_speciation_alg_dimensions()
                        call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_eq_reactions()
                        call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_stoich_mat_react_zone()
                        call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_stoich_mat_sol_rz()
                        !call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%allocate_minerals_react_zone(init_min_zones(i)%reactive_zone%num_minerals)
                        !this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%minerals=init_min_zones(i)%reactive_zone%minerals
                        !call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_num_mins_cst_act(init_min_zones(i)%reactive_zone%num_minerals_cst_act)
                        !call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_num_mins_var_act(init_min_zones(i)%reactive_zone%num_minerals_var_act)
                    end do
                end do
            else
                !allocate(this%reactive_zones(nmtype))
                call this%allocate_reactive_zones(nmtype)
                do i=1,nmtype
                    call this%reactive_zones(i)%assign_react_zone(init_min_zones(i)%reactive_zone)
                end do
            end if
        !end if
    end if
 end subroutine