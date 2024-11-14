!> Reads initial target waters and their associated target solids (if present)
!> We assume file has already been opened
subroutine read_target_waters_init(this,unit,water_types,init_sol_types,init_gas_types,niter,CV_flag)
    use chemistry_Lagr_m
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: unit !> file
    class(water_type_c), intent(in) :: water_types(:) !> water types
    !class(water_type_c), intent(in) :: bd_water_types(:) !> boundary water types
    class(solid_type_c), intent(in) :: init_sol_types(:) !> includes mineral and adsorption zones (in that order)
    class(gas_type_c), intent(in) :: init_gas_types(:) !> initial gas boundary zones (CHEPROO)
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    !class(aq_phase_c), intent(out), optional :: aq_phase_new
    
    integer(kind=4) :: counter_cols,ind,i,j,k,m,nwtype,num_tar_wat,tar_wat_ind,wtype,istype,nstype,nbwtype,bwtype,mix_wat_ind,ngzns,igzn,num_tar_wat_ext,aux_col
    real(kind=8), allocatable :: c_nc(:),u_init(:,:),c1_init(:,:),c2_init(:),c2_ig(:),gamma_2aq(:)
    character(len=256) :: label,str
    logical :: flag_comp,flag_surf,flag_aq_phase
    integer(kind=4), allocatable :: cols(:),aux_cols(:)
    !type(aq_phase_c), target :: aq_phase_new
    
    !aq_phase_copy=this%chem_syst%aq_phase
    
    nwtype=size(water_types)
    !nbwtype=size(bd_water_types)
    nstype=size(init_sol_types)
    ngzns=size(init_gas_types)
    
    flag_comp=.true.
    !if (num_surf_compl>0) then
    !>    flag_surf=.true.
    !else
    !>    flag_surf=.false.
    !end if
    !print *, this%chem_syst%aq_phase%wat_flag
    allocate(cols(2),aux_cols(2)) !> chapuza
    counter_cols=0 !> chapuza
    aux_cols=0 !> chapuza
    do
        read(unit,*) label
        if (label=='end') then
            exit
        else if (label=='TARGET WATERS') then
            read(unit,*) num_tar_wat
            call this%allocate_target_waters(num_tar_wat)
            read(unit,*) str
            if (str=='external waters') then
                read(unit,*) num_tar_wat_ext
                call this%set_num_ext_waters(num_tar_wat_ext)
                allocate(this%ext_waters_indices(num_tar_wat_ext))
                do i=1,this%num_ext_waters
                    read(unit,*) tar_wat_ind, wtype
                    if (tar_wat_ind<1 .or. tar_wat_ind>this%num_ext_waters) then
                        error stop
                    else if (wtype<1 .or. wtype>nwtype) then
                        error stop
                    else
                        this%target_waters(tar_wat_ind)=water_types(wtype)%aq_chem
                        call this%target_waters(tar_wat_ind)%set_chem_syst_aq_chem(this%chem_syst)
                        call this%target_waters(tar_wat_ind)%set_ind_aq_phase()
                        call this%target_waters(tar_wat_ind)%set_speciation_alg_dimensions(flag_comp)
                        call this%target_waters(tar_wat_ind)%compute_speciation_alg_arrays(flag_aq_phase,cols)
                        call this%target_waters(tar_wat_ind)%set_prim_species_indices()
                        call this%target_waters(tar_wat_ind)%set_sec_var_act_species_indices()
                        call this%target_waters(tar_wat_ind)%compute_U_SkT_prod()
                        call this%target_waters(tar_wat_ind)%allocate_reaction_rates_aq_chem()
                    end if
                end do
            end if
            read(unit,*) str
            if (str=='initial target waters') then
                call this%allocate_target_waters_init(this%num_target_waters-this%num_ext_waters)
                if (nstype>0) then
                    call this%allocate_target_solids(this%num_target_waters_init) !> we assume bijection with target waters
                end if
                if (ngzns>0) then
                    call this%allocate_target_gases(this%num_target_waters_init) !> we assume bijection with target waters
                end if
                !call this%allocate_ext_waters()
            !if (nstype>0 .and. nwtype>0) then
            !    read(unit,*) tar_wat_ind, wtype, istype
            !!> El if de abajo deberia estar en una subrutina
            !    if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters_init) then
            !        error stop
            !    else if (bwtype<1 .or. wtype>nwtype) then
            !        error stop
            !    else if (istype<1 .or. istype>nstype) then
            !        error stop
            !    else
            !        this%target_waters_init(tar_wat_ind)=bd_water_types(bwtype)%aq_chem
            !        call this%target_waters_init(tar_wat_ind)%set_chem_syst_aq_chem(this%chem_syst) !> chapuza
            !        call this%target_waters_init(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase) !> chapuza
            !        this%target_solids_init(tar_wat_ind)=init_sol_types(istype)%solid_chem
            !        !call this%target_waters_init(tar_wat_ind)%set_solid_chemistry(init_min_zones(imtype))
            !        !print *, init_sol_types(istype)%solid_chem%vol_fracts
            !        call this%target_solids_init(tar_wat_ind)%set_reactive_zone(init_sol_types(istype)%solid_chem%reactive_zone)
            !        !call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(1)) !> autentica chapuza ad hoc
            !        !call this%target_solids(tar_wat_ind)%allocate_concentrations()
            !        !print *, this%target_solids(tar_wat_ind)%reactive_zone%Se
            !        call this%target_waters_init(tar_wat_ind)%set_solid_chemistry(this%target_solids_init(tar_wat_ind))
            !        call this%target_waters_init(tar_wat_ind)%set_speciation_alg_dimensions(flag_comp)
            !        call this%target_waters_init(tar_wat_ind)%compute_speciation_alg_arrays()
            !        call this%target_waters_init(tar_wat_ind)%set_prim_species_indices()
            !        call this%target_waters_init(tar_wat_ind)%set_sec_var_act_species_indices()
            !        call this%target_waters_init(tar_wat_ind)%compute_U_SkT_prod()
            !        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_aq_from_c1_aq_expl() !> rezaei
            !        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_from_c1_aq_Picard(niter,CV_flag) !> rezaei
            !        !call this%target_waters_init(tar_wat_ind)%compute_activities_aq_var_act_species()
            !        !call this%target_waters_init(tar_wat_ind)%compute_pH()
            !        !call this%target_waters_init(tar_wat_ind)%compute_act_water()
            !        !call this%target_waters_init(tar_wat_ind)%compute_salinity()
            !        !call this%target_waters_init(tar_wat_ind)%compute_alkalinity()
            !        call this%target_waters_init(tar_wat_ind)%compute_molarities()
            !        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_from_u_Newton(tolerance,rel_tolerance,control_factor,niter_max,niter,CV_flag)
            !        if (init_sol_types(istype)%solid_chem%reactive_zone%cat_exch_zone%num_exch_cats>0) then
            !            call this%target_waters_init(tar_wat_ind)%compute_c2nc_from_c1_Picard(niter,CV_flag)
            !            c_nc=this%target_waters_init(tar_wat_ind)%get_c_nc_exch()
            !            call this%target_waters_init(tar_wat_ind)%compute_conc_comp(c_nc)
            !        else
            !            call this%target_waters_init(tar_wat_ind)%compute_conc_comp_aq()
            !        end if
            !        call this%target_waters_init(tar_wat_ind)%allocate_reaction_rates()
            !        
            !        !call this%ext_waters(tar_wat_ind)%set_aq_phase(this%target_waters_init(tar_wat_ind)%aq_phase)
            !        !call this%ext_waters(tar_wat_ind)%set_solid_chemistry(this%target_waters_init(tar_wat_ind)%solid_chemistry)
            !        !call this%ext_waters(tar_wat_ind)%set_speciation_alg(this%target_waters_init(tar_wat_ind)%speciation_alg)
            !        !call this%ext_waters(tar_wat_ind)%allocate_conc_aq_species()
            !        !call this%ext_waters(tar_wat_ind)%allocate_conc_comp_aq()
            !        !!call this%ext_waters(tar_wat_ind)%set_conc_aq_species()
            !        !if (init_sol_types(istype)%solid_chem%reactive_zone%cat_exch_zone%num_exch_cats>0) then
            !        !    !call this%ext_waters(tar_wat_ind)%compute_c2nc_from_c1_Picard(niter,CV_flag)
            !        !    c_nc=this%ext_waters(tar_wat_ind)%get_c_nc_exch()
            !        !    call this%ext_waters(tar_wat_ind)%compute_conc_comp(c_nc)
            !        !else
            !        !    call this%ext_waters(tar_wat_ind)%compute_conc_comp_aq()
            !        !end if
            !    end if
            !end if
            if (nstype>0 .or. ngzns>0) then
            !> El do de abajo deberia estar en una subrutina
                do i=1,this%num_target_waters_init
                    read(unit,*) tar_wat_ind, wtype, istype, igzn
                    if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
                        error stop
                    else if (wtype<1 .or. wtype>nwtype) then
                        error stop
                    else if (istype<0 .or. istype>nstype) then
                        error stop
                    else if (igzn<0 .or. igzn>ngzns) then
                        error stop
                    else
                        this%target_waters_init(tar_wat_ind-this%num_ext_waters)=water_types(wtype)%aq_chem
                        !if (aux_cols
                        !allocate(this%indices_aq_phase(this%chem_syst%aq_phase%num_species))
                        !> chapuza
                        if (counter_cols==0) then
                            call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_chem_syst_aq_chem(this%chem_syst)
                            call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_ind_aq_phase()
                        end if
                        !> chapuza
                        if (istype>0) then
                            this%target_solids(tar_wat_ind-this%num_ext_waters)=init_sol_types(istype)%solid_chem
                            call this%target_solids(tar_wat_ind-this%num_ext_waters)%set_reactive_zone(init_sol_types(istype)%solid_chem%reactive_zone)
                            call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_solid_chemistry(this%target_solids(tar_wat_ind-this%num_ext_waters))
                            !print *, this%target_waters_init(tar_wat_ind)%solid_chemistry%concentrations
                        end if
                        !> chapuza
                        if (igzn>0) then
                            call init_gas_types(igzn)%gas_chem%allocate_conc_gases() !> chapuza
                            call init_gas_types(igzn)%gas_chem%compute_vol_gas_act_coeffs() !> chapuza
                            call init_gas_types(igzn)%gas_chem%compute_conc_gases_ideal() !> chapuza
                            this%target_gases(tar_wat_ind-this%num_ext_waters)=init_gas_types(igzn)%gas_chem
                            call this%target_gases(tar_wat_ind-this%num_ext_waters)%set_reactive_zone(init_gas_types(igzn)%gas_chem%reactive_zone)
                            call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_gas_chemistry(this%target_gases(tar_wat_ind-this%num_ext_waters))
                            !print *, this%target_waters_init(tar_wat_ind)%gas_chemistry%concentrations
                        end if
                        call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_speciation_alg_dimensions(flag_comp)
                        call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%compute_speciation_alg_arrays(flag_aq_phase,cols)
                        !print *, this%target_waters_init(tar_wat_ind-this%num_ext_waters)%aq_phase%num_species
                        if (flag_aq_phase==.true.) then
                            counter_cols=counter_cols+1
                            aux_cols=cols
                            !this%target_waters_init(tar_wat_ind-this%num_ext_waters)%indices_aq_phase(cols(2))=cols(1)
                            !this%target_waters_init(tar_wat_ind-this%num_ext_waters)%indices_aq_phase(cols(1))=aux_col
                            do j=tar_wat_ind-this%num_ext_waters,this%num_target_waters_init
                                call this%target_waters_init(j)%set_chem_syst_aq_chem(this%chem_syst)
                                call this%target_waters_init(j)%set_ind_aq_phase()
                                !allocate(this%indices_aq_phase(this%chem_syst%aq_phase%num_species))
                                this%target_waters_init(j)%indices_aq_phase(cols(2))=cols(1)
                                this%target_waters_init(j)%indices_aq_phase(cols(1))=aux_cols(2)
                                !nullify(this%target_waters_init(j)%aq_phase)
                                !call this%target_waters_init(j)%set_aq_phase(aq_phase_new)
                            end do
                        end if
                        if (aux_cols(1)>0 .AND. aux_cols(2)>0) then
                            !> esto debería estar en una subrutina
                            !print *, allocated(this%target_waters_init(tar_wat_ind-this%num_ext_waters)%indices_aq_phase)
                            call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_ind_aq_phase()
                            this%target_waters_init(tar_wat_ind-this%num_ext_waters)%indices_aq_phase(aux_cols(2))=aux_cols(1)
                            this%target_waters_init(tar_wat_ind-this%num_ext_waters)%indices_aq_phase(aux_cols(1))=aux_cols(2)
                            this%target_waters_init(tar_wat_ind-this%num_ext_waters)%concentrations(aux_cols(1))=water_types(wtype)%aq_chem%concentrations(aux_cols(2))
                            this%target_waters_init(tar_wat_ind-this%num_ext_waters)%concentrations(aux_cols(2))=water_types(wtype)%aq_chem%concentrations(aux_cols(1))
                            this%target_waters_init(tar_wat_ind-this%num_ext_waters)%activities(aux_cols(1))=water_types(wtype)%aq_chem%activities(aux_cols(2))
                            this%target_waters_init(tar_wat_ind-this%num_ext_waters)%activities(aux_cols(2))=water_types(wtype)%aq_chem%activities(aux_cols(1))
                            this%target_waters_init(tar_wat_ind-this%num_ext_waters)%log_act_coeffs(aux_cols(1))=water_types(wtype)%aq_chem%log_act_coeffs(aux_cols(2))
                            this%target_waters_init(tar_wat_ind-this%num_ext_waters)%log_act_coeffs(aux_cols(2))=water_types(wtype)%aq_chem%log_act_coeffs(aux_cols(1))
                            !this%target_waters_init(tar_wat_ind-this%num_ext_waters)%log_Jacobian_act_coeffs(:,aux_cols(1))=water_types(wtype)%aq_chem%log_Jacobian_act_coeffs(:,aux_cols(2))
                            !this%target_waters_init(tar_wat_ind-this%num_ext_waters)%log_Jacobian_act_coeffs(:,aux_cols(2))=water_types(wtype)%aq_chem%log_Jacobian_act_coeffs(:,aux_cols(1))                        
                        end if
                        !print *, this%target_waters_init(tar_wat_ind-this%num_ext_waters)%aq_phase%wat_flag
                        call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_prim_species_indices()
                        call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_sec_var_act_species_indices()
                        call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%compute_U_SkT_prod()
                        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_aq_from_c1_aq_expl() !> rezaei
                        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_from_c1_aq_Picard(niter,CV_flag) !> rezaei
                        !call this%target_waters_init(tar_wat_ind)%compute_activities_aq_var_act_species()
                        !call this%target_waters_init(tar_wat_ind)%compute_pH()
                        !call this%target_waters_init(tar_wat_ind)%compute_act_water()
                        !call this%target_waters_init(tar_wat_ind)%compute_salinity()
                        !call this%target_waters_init(tar_wat_ind)%compute_alkalinity()
                        !call this%target_waters_init(tar_wat_ind)%compute_molarities()
                        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_from_u_Newton(tolerance,rel_tolerance,control_factor,niter_max,niter,CV_flag)
                        if (istype>0) then
                            if (init_sol_types(istype)%solid_chem%reactive_zone%cat_exch_zone%num_exch_cats>0) then
                                allocate(c2_ig(this%target_waters_init(tar_wat_ind-this%num_ext_waters)%speciation_alg%num_eq_reactions))
                                allocate(c2_init(this%target_waters_init(tar_wat_ind-this%num_ext_waters)%speciation_alg%num_eq_reactions))
                                c2_ig=1d-16
                                call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%compute_c2nc_from_c1_Picard(c2_ig,c2_init,niter,CV_flag)
                                c_nc=this%target_waters_init(tar_wat_ind-this%num_ext_waters)%get_conc_nc()
                                !call this%target_waters_init(tar_wat_ind)%compute_conc_comp(c_nc)
                                deallocate(c2_ig,c2_init)
                            else
                                !call this%target_waters_init(tar_wat_ind)%compute_conc_comp_aq()
                            end if
                        end if
                        call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%allocate_reaction_rates_aq_chem()
                        !call this%ext_waters(tar_wat_ind)%set_aq_phase(this%target_waters_init(tar_wat_ind)%aq_phase)
                        !call this%ext_waters(tar_wat_ind)%set_solid_chemistry(this%target_waters_init(tar_wat_ind)%solid_chemistry)
                        !call this%ext_waters(tar_wat_ind)%set_speciation_alg(this%target_waters_init(tar_wat_ind)%speciation_alg)
                        !call this%ext_waters(tar_wat_ind)%allocate_conc_aq_species()
                        !call this%ext_waters(tar_wat_ind)%allocate_conc_comp_aq()
                        !!call this%ext_waters(tar_wat_ind)%set_conc_aq_species()
                        !if (init_sol_types(istype)%solid_chem%reactive_zone%cat_exch_zone%num_exch_cats>0) then
                        !    !call this%ext_waters(tar_wat_ind)%compute_c2nc_from_c1_Picard(niter,CV_flag)
                        !    c_nc=this%ext_waters(tar_wat_ind)%get_c_nc_exch()
                        !    call this%ext_waters(tar_wat_ind)%compute_conc_comp(c_nc)
                        !else
                        !    call this%ext_waters(tar_wat_ind)%compute_conc_comp_aq()
                        !end if
                    end if
                end do
            end if
            !> El if de abajo deberia estar en una subrutina
            !if (nstype>0 .and. nbwtype>0) then
            !!    read(unit,*) tar_wat_ind, bwtype, istype
            !!    if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters_init) then
            !!        error stop
            !!    else if (bwtype<1 .or. bwtype>nbwtype) then
            !!        error stop
            !!    else if (istype<1 .or. istype>nstype) then
            !!        error stop
            !!    else
            !!        this%target_waters_init(tar_wat_ind)=bd_water_types(bwtype)%aq_chem
            !!        call this%target_waters_init(tar_wat_ind)%set_chem_syst_aq_chem(this%chem_syst) !> chapuza
            !!        call this%target_waters_init(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase) !> chapuza
            !!        !call this%target_waters_init(tar_wat_ind)%set_aq_phase(bd_water_types(bwtype)%aq_chem%aq_phase)
            !!        this%target_solids_init(tar_wat_ind)=init_sol_types(istype)%solid_chem
            !!        !call this%target_waters_init(tar_wat_ind)%set_solid_chemistry(init_min_zones(imtype))
            !!        !print *, init_sol_types(istype)%reactive_zone%Se
            !!        call this%target_solids_init(tar_wat_ind)%set_reactive_zone(init_sol_types(istype)%solid_chem%reactive_zone)
            !!        !call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(1)) !> autentica chapuza ad hoc
            !!        !call this%target_solids(tar_wat_ind)%allocate_concentrations()
            !!        !print *, this%target_solids(tar_wat_ind)%reactive_zone%Se
            !!        call this%target_waters_init(tar_wat_ind)%set_solid_chemistry(this%target_solids_init(tar_wat_ind))
            !!        call this%target_waters_init(tar_wat_ind)%set_speciation_alg_dimensions(flag_comp)
            !!        call this%target_waters_init(tar_wat_ind)%compute_speciation_alg_arrays()
            !!        call this%target_waters_init(tar_wat_ind)%set_prim_species_indices()
            !!        call this%target_waters_init(tar_wat_ind)%set_sec_var_act_species_indices()
            !!        call this%target_waters_init(tar_wat_ind)%compute_U_SkT_prod()
            !!        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_aq_from_c1_aq_expl() !> rezaei
            !!        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_from_c1_aq_Picard(niter,CV_flag) !> rezaei
            !!        !call this%target_waters_init(tar_wat_ind)%compute_activities_aq_var_act_species()
            !!        !call this%target_waters_init(tar_wat_ind)%compute_pH()
            !!        !call this%target_waters_init(tar_wat_ind)%compute_act_water()
            !!        !call this%target_waters_init(tar_wat_ind)%compute_salinity()
            !!        !call this%target_waters_init(tar_wat_ind)%compute_alkalinity()
            !!        !call this%target_waters_init(tar_wat_ind)%compute_molarities()
            !!        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_from_u_Newton(tolerance,rel_tolerance,control_factor,niter_max,niter,CV_flag)
            !!        call this%target_waters_init(tar_wat_ind)%compute_conc_comp_aq()
            !!        call this%target_waters_init(tar_wat_ind)%allocate_reaction_rates()
            !!        
            !!        call this%ext_waters(tar_wat_ind)%set_aq_phase(this%target_waters_init(tar_wat_ind)%aq_phase)
            !!        call this%ext_waters(tar_wat_ind)%set_solid_chemistry(this%target_waters_init(tar_wat_ind)%solid_chemistry)
            !!        call this%ext_waters(tar_wat_ind)%set_speciation_alg(this%target_waters_init(tar_wat_ind)%speciation_alg)
            !!        call this%ext_waters(tar_wat_ind)%allocate_conc_aq_species()
            !!        call this%ext_waters(tar_wat_ind)%allocate_conc_comp()
            !!        call this%ext_waters(tar_wat_ind)%set_conc_aq_species()
            !!        call this%ext_waters(tar_wat_ind)%compute_conc_comp_aq()
            !!    end if
            !    this%target_waters_init=this%target_waters_init
            !    this%target_solids=this%target_solids_init
            !    do i=1,this%num_target_waters_init
            !        call this%target_waters_init(i)%set_solid_chemistry(this%target_solids(i))
            !    end do
            !end if
            if (nstype==0 .and. ngzns==0) then
                read(unit,*) tar_wat_ind, wtype!, istype
            !> El if de abajo deberia estar en una subrutina
                if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters_init) then
                    error stop
                else if (wtype<1 .or. wtype>nwtype) then
                    error stop
                !else if (istype<1 .or. istype>nstype) then
                !    error stop
                else
                    this%target_waters_init(tar_wat_ind-this%num_ext_waters)=water_types(wtype)%aq_chem
                    call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_chem_syst_aq_chem(this%chem_syst) !> chapuza
                    call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_ind_aq_phase() !> chapuza
                    !this%target_solids_init(tar_wat_ind)=init_sol_types(istype)%solid_chem
                    !call this%target_waters_init(tar_wat_ind)%set_solid_chemistry(init_min_zones(imtype))
                    !print *, init_sol_types(istype)%solid_chem%vol_fracts
                    !call this%target_solids_init(tar_wat_ind)%set_reactive_zone(init_sol_types(istype)%solid_chem%reactive_zone)
                    !call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(1)) !> autentica chapuza ad hoc
                    !call this%target_solids(tar_wat_ind)%allocate_concentrations()
                    !print *, this%target_solids(tar_wat_ind)%reactive_zone%Se
                    !call this%target_waters_init(tar_wat_ind)%set_solid_chemistry(this%target_solids_init(tar_wat_ind))
                    call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_speciation_alg_dimensions(flag_comp)
                    call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%compute_speciation_alg_arrays(flag_aq_phase,cols)
                    call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_prim_species_indices()
                    call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%set_sec_var_act_species_indices()
                    call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%compute_U_SkT_prod()
                    !call this%target_waters_init(tar_wat_ind)%compute_c_nc_aq_from_c1_aq_expl() !> rezaei
                    !call this%target_waters_init(tar_wat_ind)%compute_c_nc_from_c1_aq_Picard(niter,CV_flag) !> rezaei
                    !call this%target_waters_init(tar_wat_ind)%compute_activities_aq_var_act_species()
                    !call this%target_waters_init(tar_wat_ind)%compute_pH()
                    !call this%target_waters_init(tar_wat_ind)%compute_act_water()
                    !call this%target_waters_init(tar_wat_ind)%compute_salinity()
                    !call this%target_waters_init(tar_wat_ind)%compute_alkalinity()
                    !call this%target_waters_init(tar_wat_ind)%compute_molarities()
                    !call this%target_waters_init(tar_wat_ind)%compute_c_nc_from_u_Newton(tolerance,rel_tolerance,control_factor,niter_max,niter,CV_flag)
                    !if (init_sol_types(istype)%solid_chem%reactive_zone%cat_exch_zone%num_exch_cats>0) then
                    !    call this%target_waters_init(tar_wat_ind)%compute_c2nc_from_c1_Picard(niter,CV_flag)
                    !    c_nc=this%target_waters_init(tar_wat_ind)%get_c_nc_exch()
                    !    call this%target_waters_init(tar_wat_ind)%compute_conc_comp(c_nc) 
                    !else
                    !    call this%target_waters_init(tar_wat_ind)%compute_conc_comp_aq()
                    !end if
                    call this%target_waters_init(tar_wat_ind-this%num_ext_waters)%allocate_reaction_rates_aq_chem()
                    !call this%ext_waters(tar_wat_ind)%set_aq_phase(this%target_waters_init(tar_wat_ind)%aq_phase)
                    !!call this%ext_waters(tar_wat_ind)%set_solid_chemistry(this%target_waters_init(tar_wat_ind)%solid_chemistry)
                    !call this%ext_waters(tar_wat_ind)%set_speciation_alg(this%target_waters_init(tar_wat_ind)%speciation_alg)
                    !call this%ext_waters(tar_wat_ind)%allocate_conc_aq_species()
                    !call this%ext_waters(tar_wat_ind)%allocate_conc_comp_aq()
                    !!call this%ext_waters(tar_wat_ind)%set_conc_aq_species()
                    !if (init_sol_types(istype)%solid_chem%reactive_zone%cat_exch_zone%num_exch_cats>0) then
                    !    !call this%ext_waters(tar_wat_ind)%compute_c2nc_from_c1_Picard(niter,CV_flag)
                    !    c_nc=this%ext_waters(tar_wat_ind)%get_c_nc_exch()
                    !    call this%ext_waters(tar_wat_ind)%compute_conc_comp(c_nc)
                    !else
                    !    call this%ext_waters(tar_wat_ind)%compute_conc_comp_aq()
                    !end if
                end if
                !do i=2,this%num_target_waters_init
                !    read(unit,*) tar_wat_ind, iwtype
                !    if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters_init) then
                !        error stop
                !    else if (iwtype<1 .or. iwtype>niwtype) then
                !        error stop
                !    else
                !        this%target_waters_init(tar_wat_ind)=init_water_types(iwtype)%aq_chem
                !        call this%target_waters_init(tar_wat_ind)%set_chem_syst_aq_chem(this%chem_syst)
                !        !if (this%target_waters_init(tar_wat_ind)%chem_syst%num_eq_reacts_homog>0) then
                !        !    call this%target_waters_init(tar_wat_ind)%set_speciation_alg_dimensions(flag_comp)
                !        !    call this%target_waters_init(tar_wat_ind)%compute_speciation_alg_arrays()
                !        !    !call this%target_waters_init(tar_wat_ind)%set_conc_sec_var_act_species() !> initial guess c_nc
                !        !    if (allocated(this%target_waters_init(tar_wat_ind)%conc_comp)) then
                !        !        !call this%target_waters_init(tar_wat_ind)%compute_c_nc_from_u_aq_Newton(niter,CV_flag)
                !        !    else
                !        !        call this%target_waters_init(tar_wat_ind)%compute_c2nc_from_c1_aq_Picard(niter,CV_flag)
                !        !        if (init_sol_types(istype)%solid_chem%reactive_zone%cat_exch_zone%num_exch_cats>0) then
                !        !            call this%target_waters_init(tar_wat_ind)%compute_c2nc_from_c1_Picard(niter,CV_flag)
                !        !            c_nc=this%target_waters_init(tar_wat_ind)%get_c_nc_exch()
                !        !            call this%target_waters_init(tar_wat_ind)%compute_conc_comp(c_nc)
                !        !        else
                !        !            call this%target_waters_init(tar_wat_ind)%compute_conc_comp_aq()
                !        !        end if
                !        !    end if
                !        !    call this%target_waters_init(tar_wat_ind)%compute_U_SkT_prod()
                !        !end if
                !        call this%target_waters_init(tar_wat_ind)%allocate_reaction_rates()
                !        !call this%ext_waters(tar_wat_ind)%set_aq_phase(this%target_waters_init(tar_wat_ind)%aq_phase)
                !        !call this%ext_waters(tar_wat_ind)%allocate_conc_aq_species()
                !        !call this%ext_waters(tar_wat_ind)%allocate_conc_comp_aq(this%target_waters_init(tar_wat_ind)%speciation_alg%num_aq_prim_species)
                !        !if (init_sol_types(istype)%solid_chem%reactive_zone%cat_exch_zone%num_exch_cats>0) then
                !        !    call this%ext_waters(tar_wat_ind)%compute_c2nc_from_c1_Picard(niter,CV_flag)
                !        !    c_nc=this%ext_waters(tar_wat_ind)%get_c_nc_exch()
                !        !    call this%ext_waters(tar_wat_ind)%compute_conc_comp(c_nc)
                !        !else
                !        !    call this%ext_waters(tar_wat_ind)%compute_conc_comp_aq()
                !        !end if
                !    end if
                !end do
            end if
            end if
        else
            continue
        end if
    end do
    this%target_waters(this%num_ext_waters+1:this%num_target_waters_init)=this%target_waters_init
end subroutine