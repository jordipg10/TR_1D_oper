!> Reads initial target waters and their associated target solids and/or gases
!> We assume file has already been opened
subroutine read_target_waters_init(this,unit,water_types,init_sol_types,init_gas_types,niter,CV_flag)
    use chemistry_Lagr_m
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: unit !> file
    type(aqueous_chemistry_c), intent(in) :: water_types(:) !> water types
    type(solid_chemistry_c), intent(in) :: init_sol_types(:) !> includes adsorption & mineral zones (in that order)
    type(gas_chemistry_c), intent(in) :: init_gas_types(:) !> initial gas zones
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    
    integer(kind=4) :: ind_ext,ind_dom,counter_cols,ind,i,j,k,m,nwtype,num_tar_wat,tar_wat_ind,wtype,istype,nstype,nbwtype,bwtype,mix_wat_ind,ngzns,igzn,num_tar_wat_ext,aux_col
    real(kind=8), allocatable :: c_nc(:),u_init(:,:),c1_init(:),c2_init(:),c2_ig(:),gamma_2aq(:)
    character(len=256) :: label,str
    logical :: flag_comp,flag_surf,flag_aq_phase,flag
    integer(kind=4), allocatable :: cols(:),aux_cols(:)
    type(solid_chemistry_c), target :: solid_chem !> default object
    type(reactive_zone_c), target :: react_zone   !> default object
    
    nwtype=size(water_types)
    nstype=size(init_sol_types)
    ngzns=size(init_gas_types)
    
    flag_comp=.true. !> by default
    
    allocate(cols(2),aux_cols(2)) !> chapuza
    counter_cols=0 !> chapuza
    aux_cols=0 !> chapuza
    
    ind_Ext=0 !> counter external waters
    ind_dom=0 !> counter domain waters
    do
        read(unit,*) label
        if (label=='end') then
            exit
        else if (label=='TARGET WATERS') then
            read(unit,*) num_tar_wat
            call this%allocate_target_waters(num_tar_wat)
            if (nstype>0) then
                call this%allocate_target_solids(this%num_target_waters) !> we assume bijection with target waters
            end if
            if (ngzns>0) then
                call this%allocate_target_gases(this%num_target_waters) !> we assume bijection with target waters
            end if
            read(unit,*) num_tar_wat_ext
            call this%allocate_ext_waters_indices(num_tar_wat_ext)
            call this%allocate_dom_tar_wat_indices(this%num_target_waters-this%num_ext_waters)
                do i=1,this%num_target_waters
                    read(unit,*) tar_wat_ind, wtype, istype, igzn, flag
                    if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
                        error stop
                    else if (wtype<1 .or. wtype>nwtype) then
                        error stop
                    else if (istype<0 .or. istype>nstype) then
                        error stop
                    else if (igzn<0 .or. igzn>ngzns) then
                        error stop
                    else if (flag==.true.) then
                        ind_ext=ind_ext+1
                        this%ext_waters_indices(ind_ext)=tar_wat_ind
                    else
                        ind_dom=ind_dom+1
                        this%dom_tar_wat_indices(ind_dom)=tar_wat_ind
                    end if
                    this%target_waters(tar_wat_ind)=water_types(wtype)
                    if (counter_cols==0) then
                        call this%target_waters(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase)
                        call this%target_waters(tar_wat_ind)%set_ind_aq_phase()
                    end if
                    if (istype>0) then
                        this%target_solids(tar_wat_ind)=init_sol_types(istype)
                        call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+istype))
                        call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind))
                        if (igzn>0) then
                            call init_gas_types(igzn)%allocate_conc_gases() !> chapuza
                            call init_gas_types(igzn)%compute_vol_gas_act_coeffs() !> chapuza
                            call init_gas_types(igzn)%compute_conc_gases_ideal() !> chapuza
                            this%target_gases(tar_wat_ind)=init_gas_types(igzn)
                            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+istype))
                            call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                        end if
                    else if (igzn>0) then
                        call init_gas_types(igzn)%allocate_conc_gases() !> chapuza
                        call init_gas_types(igzn)%compute_vol_gas_act_coeffs() !> chapuza
                        call init_gas_types(igzn)%compute_conc_gases_ideal() !> chapuza
                        this%target_gases(tar_wat_ind)=init_gas_types(igzn)
                        call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(igzn))
                        call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                        call solid_chem%set_reactive_zone(this%reactive_zones(igzn))
                        call this%target_waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
                    else
                        call solid_chem%set_reactive_zone(react_zone)
                        call this%target_waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
                    end if
                    call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
                    call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_aq_phase,cols)
                    if (this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats>0) then
                        c1_init=this%target_waters(tar_wat_ind)%get_c1()
                        allocate(c2_init(this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
                        if (THIS%act_coeffs_model==0) then
                            call this%target_waters(tar_wat_ind)%compute_c2nc_from_c1_ideal(c1_init,c2_init)
                        else
                            allocate(c2_ig(this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
                            c2_ig=1d-16
                            call this%target_waters(tar_wat_ind)%compute_c2nc_from_c1_Picard(c1_init,c2_ig,c2_init,niter,CV_flag)
                            deallocate(c2_ig)
                        end if
                        deallocate(c1_init,c2_init)
                    end if
                    if (flag_aq_phase==.true.) then
                        counter_cols=counter_cols+1
                        aux_cols=cols
                        do j=tar_wat_ind,this%num_target_waters
                            call this%target_waters(j)%set_aq_phase(this%chem_syst%aq_phase)
                            call this%target_waters(j)%set_ind_aq_phase()
                            this%target_waters(j)%indices_aq_phase(cols(2))=cols(1)
                            this%target_waters(j)%indices_aq_phase(cols(1))=aux_cols(2)
                        end do
                    end if
                    if (aux_cols(1)>0 .AND. aux_cols(2)>0) then
                        call this%target_waters(tar_wat_ind)%set_ind_aq_phase()
                        this%target_waters(tar_wat_ind)%indices_aq_phase(aux_cols(2))=aux_cols(1)
                        this%target_waters(tar_wat_ind)%indices_aq_phase(aux_cols(1))=aux_cols(2)
                        this%target_waters(tar_wat_ind)%concentrations(aux_cols(1))=water_types(wtype)%concentrations(aux_cols(2))
                        this%target_waters(tar_wat_ind)%concentrations(aux_cols(2))=water_types(wtype)%concentrations(aux_cols(1))
                        this%target_waters(tar_wat_ind)%activities(aux_cols(1))=water_types(wtype)%activities(aux_cols(2))
                        this%target_waters(tar_wat_ind)%activities(aux_cols(2))=water_types(wtype)%activities(aux_cols(1))
                        this%target_waters(tar_wat_ind)%log_act_coeffs(aux_cols(1))=water_types(wtype)%log_act_coeffs(aux_cols(2))
                        this%target_waters(tar_wat_ind)%log_act_coeffs(aux_cols(2))=water_types(wtype)%log_act_coeffs(aux_cols(1))
                    end if
                    call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod()
                    call this%target_waters(tar_wat_ind)%allocate_reaction_rates()
                end do
        else
            continue
        end if
    end do
    this%target_waters_init=this%target_waters
    this%target_solids_init=this%target_solids
    this%target_gases_init=this%target_gases
end subroutine