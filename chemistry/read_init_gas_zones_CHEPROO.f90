subroutine read_init_gas_zones_CHEPROO(this,unit,gas_zones)
    use chemistry_Lagr_m
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: unit !> file
    type(gas_chemistry_c), intent(out), allocatable :: gas_zones(:)
    !type(reactive_zone_c), intent(inout), allocatable, optional :: reactive_zones(:)
    
    integer(kind=4) :: num_rz,i,j,k,ngtype,igtype,nrwtype,gas_ind,num_gas_zones,num_gases_loc,num_gases_var,num_gases_cst,n_gas_eq,n_gas_kin,num_surf_rz
    integer(kind=4), allocatable :: ind_gases(:)
    character(len=256) :: str,constrain,label
    real(kind=8) :: guess,conc,temp,part_press,vol
    type(gas_c) :: gas
    type(gas_phase_c) :: gas_phase
    type(reactive_zone_c) :: react_zone
    type(reactive_zone_c), allocatable :: aux_react_zones(:)
    logical :: flag

    read(unit,*) ngtype !> number of gas zones
    
    allocate(gas_zones(ngtype))
    allocate(ind_gases(this%chem_syst%gas_phase%num_species)) !> chapuza

    num_gas_zones=0 !> counter gas zones
    
    do
        read(unit,*) igtype, temp, vol
        if (igtype<1 .or. igtype>ngtype) error stop
        read(unit,*) str
        if (str=='*') then
            if (num_gas_zones<ngtype) then
                continue
            else
                exit
            end if
        else if (index(str,'gas')/=0) then
            num_gas_zones=num_gas_zones+1
            num_gases_loc=0 !> counter gases in this zone
            num_gases_var=0
            num_gases_cst=0
            n_gas_eq=0
            n_gas_kin=0
            do
                read(unit,*) gas%name, part_press
                if (gas%name=='*') then
                    call react_zone%set_chem_syst_react_zone(this%chem_syst)
                    call react_zone%set_CV_params(this%CV_params)
                    call react_zone%gas_phase%allocate_gases(num_gases_loc)
                    call react_zone%gas_phase%set_num_var_act_species_phase(num_gases_var)
                    call react_zone%gas_phase%set_num_cst_act_species_phase(num_gases_cst)
                    call react_zone%gas_phase%set_num_gases_eq(n_gas_eq) !> 
                    call react_zone%gas_phase%set_num_gases_kin(n_gas_kin) !> 
                    call gas_zones(igtype)%set_reactive_zone(react_zone)
                    call gas_zones(igtype)%set_temp(temp+273.15) !> Kelvin
                    call gas_zones(igtype)%set_volume(vol)
                    call gas_zones(igtype)%allocate_partial_pressures()
                    call gas_zones(igtype)%allocate_conc_gases()
                    call gas_zones(igtype)%allocate_log_act_coeffs_gases()
                    call gas_zones(igtype)%allocate_var_act_species_indices(num_gases_var)
                    call gas_zones(igtype)%allocate_cst_act_species_indices(num_gases_cst)
                    exit
                else
                    call this%chem_syst%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                    if (flag==.true.) then
                        num_gases_loc=num_gases_loc+1
                        ind_gases(num_gases_loc)=gas_ind
                        if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag==.true. .and. gas_ind<=this%chem_syst%gas_phase%num_gases_eq) then
                            num_gases_cst=num_gases_cst+1
                            n_gas_eq=n_gas_eq+1
                        !> we modify equilibrium constant
                            this%chem_syst%eq_reacts(this%chem_syst%aq_phase%num_aq_complexes+this%chem_syst%num_minerals_eq+gas_ind)%eq_cst=this%chem_syst%eq_reacts(this%chem_syst%aq_phase%num_aq_complexes+this%chem_syst%num_minerals_eq+gas_ind)%eq_cst/part_press
                        else if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag==.false. .and. gas_ind<=this%chem_syst%gas_phase%num_gases_eq) then
                            num_gases_var=num_gases_var+1
                            n_gas_eq=n_gas_eq+1
                        else if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag==.true. .and. gas_ind>this%chem_syst%gas_phase%num_gases_eq) then
                            num_gases_cst=num_gases_cst+1
                            n_gas_kin=n_gas_kin+1
                        else
                            num_gases_var=num_gases_var+1
                            n_gas_kin=n_gas_kin+1
                        end if
                        !if (gas_ind<=this%chem_syst%gas_phase%num_gases_eq) then
                        !    n_gas_eq=n_gas_eq+1
                        !else
                        !    n_gas_kin=n_gas_kin+1
                        !end if
                    else
                        error stop
                    end if
                end if
            end do
            if (num_gas_zones==ngtype) exit
        else
            exit
        end if
    end do
    rewind(unit)
    num_gas_zones=0 !> counter gas zones
    !num_gases_glob=1 !> counter gases in chemical system
    do
        read(unit,*) label
        if (label=='INITIAL GAS ZONES') then
            read(unit,*) ngtype !> number of gas zones
            do i=1,ngtype
                read(unit,*) igtype, temp
                read(unit,*) str
                num_gas_zones=num_gas_zones+1
                !num_gases_loc=0 !> counter gases in this zone
                do j=1,gas_zones(igtype)%reactive_zone%gas_phase%num_species
                    read(unit,*) gas%name, part_press
                    call gas_zones(igtype)%reactive_zone%gas_phase%gases(j)%assign_species(this%chem_syst%gas_phase%gases(ind_gases(j))) !> we set gas
                    gas_zones(igtype)%activities(j)=part_press !> activities are partial pressures
                end do
                call gas_zones(igtype)%compute_conc_gases_ideal() !> we assume gases are ideal
                call gas_zones(igtype)%compute_pressure()
                call gas_zones(igtype)%compute_log_act_coeffs_gases()
                call gas_zones(igtype)%set_indices_gases()
                read(unit,*) str
            end do
        else if (label=='end') then
            backspace(unit)
            exit
        else
            continue
        end if
    end do
    !if (present(reactive_zones)) then
        if (allocated(this%reactive_zones)) then
            num_surf_rz=size(this%reactive_zones)
            allocate(aux_react_zones(num_surf_rz))
            do i=1,num_surf_rz
                call aux_react_zones(i)%assign_react_zone(this%reactive_zones(i))
            end do
            num_rz=num_surf_rz+ngtype*(1+num_surf_rz)
            !deallocate(this%reactive_zones)
            !allocate(this%reactive_zones(num_rz))
            call this%allocate_reactive_zones(num_rz)
            do i=1,num_surf_rz
                call this%reactive_zones(ngtype+i)%assign_react_zone(aux_react_zones(i))
            end do
            do i=1,ngtype
                call this%reactive_zones(i)%set_chem_syst_react_zone(gas_zones(i)%reactive_zone%chem_syst)
                call this%reactive_zones(i)%set_gas_phase(gas_zones(i)%reactive_zone%gas_phase)
                call gas_zones(i)%set_reactive_zone(this%reactive_zones(i))
            end do
            do i=1,ngtype
                do j=1,num_surf_rz
                    call this%reactive_zones(ngtype+i*num_surf_rz+j)%set_chem_syst_react_zone(gas_zones(I)%reactive_zone%chem_syst)
                    call this%reactive_zones(ngtype+i*num_surf_rz+j)%set_gas_phase(gas_zones(i)%reactive_zone%gas_phase)
                    call this%reactive_zones(ngtype+i*num_surf_rz+j)%set_cat_exch_zone(aux_react_zones(i)%cat_exch_zone)
                end do
            end do
        else
            !allocate(this%reactive_zones(ngtype))
            call this%allocate_reactive_zones(ngtype)
            do i=1,ngtype
                call this%reactive_zones(i)%assign_react_zone(gas_zones(i)%reactive_zone)
                call gas_zones(i)%set_reactive_zone(this%reactive_zones(i))
            end do
        end if
    !end if
end subroutine