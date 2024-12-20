!> This subroutine reads a water type from CHEPROO data input
subroutine read_wat_type_CHEPROO(this,n_p_aq,num_cstr,model,Jac_flag,unit,niter,CV_flag,surf_chem)
    use aqueous_chemistry_m
    implicit none
    class(aqueous_chemistry_c) :: this
    integer(kind=4), intent(in) :: n_p_aq !> number of primary aqueous species
    integer(kind=4), intent(in) :: num_cstr !> number of constrains
    !integer(kind=4), intent(in) :: num_gas_cstr !> number of gas constrains
    integer(kind=4), intent(in) :: model !> activity coefficients model
    integer(kind=4), intent(in) :: Jac_flag !> 0: incremental coeffficinets, 1: analtical
    integer(kind=4), intent(in) :: unit !> file unit
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    type(solid_chemistry_c), intent(inout), optional :: surf_chem !> surface chemistry (if there are exchange reactions)
    !real(kind=8), intent(in), optional :: c1_surf !> chapuza
    !real(kind=8), intent(in), optional :: CEC !> chapuza
    !real(kind=8), intent(out), optional :: conc_exch(:) !> chapuza
    
    integer(kind=4) :: i,j,k,l,m,niwtype,nbwtype,nrwtype,gas_ind,min_ind,n_comp_aq,ind_cstr,n_gas_constr,n_aq_comp,icon,aq_sp_ind,ind_sp
    integer(kind=4), allocatable :: icons(:),indices_constrains(:),gas_indices(:),n_icons(:),prim_indices(:),cols(:)
    real(kind=8), allocatable :: ctots(:),c2_init(:),c2_ig(:),c1(:)
    character(len=256) :: prim_sp_name,label,aq_sp_name
    character(len=256), allocatable :: constrains(:)
    real(kind=8) :: temp,conc,ionic_act,guess,ctot
    logical :: flag_gas,flag_min,flag,flag_comp,flag_surf,flag_Se,flag_sp
    
    type(reactive_zone_c) :: react_zone
    type(gas_chemistry_c) :: gas_chem
    type(gas_c) :: gas
    type(mineral_c) :: mineral
    type(aq_species_c) :: aq_species
    type(species_c) :: constrain
    
    if (present(surf_chem)) then
        call this%set_solid_chemistry(surf_chem)
    end if
    
    call this%solid_chemistry%reactive_zone%allocate_non_flowing_species(num_cstr)
    
    !    read(unit,*) label
    !    if (label=='INITIAL AND BOUNDARY WATER TYPES') then
    !        i=0 !> counter initial water types
    !        read(unit,*) model
    !        read(unit,*) this%num_init_wat_types, this%num_bd_wat_types, this%num_rech_wat_types
        if (n_p_aq<0 .or. n_p_aq>this%aq_phase%num_species) then
            error stop
        else if (num_cstr<0 .or. num_cstr>n_p_aq) then
            error stop
        else
            allocate(prim_indices(n_p_aq),icons(n_p_aq),ctots(n_p_aq),indices_constrains(num_cstr))
            !if (this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_species>0) then
            !    call react_zone%set_gas_phase(this%solid_chemistry%reactive_zone%chem_syst%gas_phase)
            !    call gas_chem%set_reactive_zone(react_zone)
            !    call gas_chem%allocate_partial_pressures()
            !    call gas_chem%allocate_conc_gases()
            !    call gas_chem%allocate_log_act_coeffs_gases()
            !    call gas_chem%set_temp(this%temp)
            !    call gas_chem%set_volume(1d0) !> arbitrary
            !    call this%set_gas_chemistry(gas_chem)
            !end if
        end if
        !call this%allocate_conc_comp(n_p_aq)
        !call this%allocate_log_act_coeffs()
        !call this%set_log_act_coeffs() 
        !call this%allocate_activities_aq_species()
        read(unit,*) label
        k=1 !> counter primary species
        l=1 !> counter constrains
        m=1 !> counter gas constrains
        if (index(label,'guess')/=0) then !> 'icon, guess, ctot, constrain'
            i=i+1
            allocate(n_icons(4))
            n_icons=0 !> number of each icon option
            do 
                read(unit,*) aq_species%name, icon, guess, ctot, constrain%name
                if (aq_species%name=='*') then
                    exit
                else
                    call this%aq_phase%is_species_in_aq_phase(aq_species,flag,aq_sp_ind)
                    if (flag==.true.) then
                        prim_indices(k)=aq_sp_ind
                        !> Chapuza
                        if (icon==1) then
                            n_icons(1)=n_icons(1)+1
                            if (aq_species%name=='h2o') then
                                guess=1d0/18d-3
                                ctot=guess
                            end if
                        else if (icon==2) then
                            n_icons(2)=n_icons(2)+1
                            n_aq_comp=n_aq_comp+1
                        else if (icon==3) then
                            n_icons(3)=n_icons(3)+1
                            if (aq_species%name=='h+') then
                                call this%aq_phase%set_ind_proton(aq_sp_ind)
                                call this%set_pH(-log10(ctot))
                            end if
                        else if (icon==4) then
                            n_icons(4)=n_icons(4)+1
                            call this%solid_chemistry%reactive_zone%chem_syst%is_eq_reaction_in_chem_syst(constrain%name,flag,ind_cstr)
                            if (flag==.true.) then
                                indices_constrains(l)=ind_cstr
                                call this%solid_chemistry%reactive_zone%chem_syst%is_species_in_chem_syst(constrain,flag_sp,ind_sp)
                                if (flag_sp==.true.) then
                                    call THIS%solid_chemistry%reactive_zone%non_flowing_species(l)%assign_species(this%solid_chemistry%reactive_zone%chem_syst%species(ind_sp))
                                    l=l+1
                                    if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq_cst_act) then
                                        this%solid_chemistry%reactive_zone%num_minerals_cst_act=this%solid_chemistry%reactive_zone%num_minerals_cst_act+1
                                        this%solid_chemistry%reactive_zone%num_minerals=this%solid_chemistry%reactive_zone%num_minerals+1
                                    else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq+this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_gases_eq_cst_act) then
                                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act=this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1
                                    else if (ind_cstr<this%solid_chemistry%reactive_zone%chem_syst%num_eq_reacts-this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_gases_eq_var_act) then
                                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats=this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats+1
                                    else
                                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act=this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act+1
                                    end if
                                else
                                    error stop
                                end if
                                        
                                !call constrain%is_gas(flag_gas)
                                !if (flag_gas==.true.) then
                                !    call gas%set_name(constrain%name)
                                !    call this%gas_chemistry%reactive_zone%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                                !    this%gas_chemistry%activities(gas_ind)=ctot
                                !else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq) then
                                !    call THIS%solid_chemistry%reactive_zone%non_flowing_species(ind_cstr)%assign_species(this%solid_chemistry%reactive_zone%chem_syst%minerals(this%solid_chemistry%reactive_zone%chem_syst%num_min_kin_reacts+ind_cstr)%mineral)
                                !else
                                !    call THIS%solid_chemistry%reactive_zone%non_flowing_species(ind_cstr-this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq-this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_gases_eq_cst_act)%assign_species(this%solid_chemistry%reactive_zone%chem_syst%cat_exch%surf_compl(this%solid_chemistry%reactive_zone%chem_syst%num_min_kin_reacts+ind_cstr)%mineral)
                                !end if
                                !if (ind_cstr>this%solid_chemistry%reactive_zone%chem_syst%num_minerals_cst_act .AND. ind_cstr<this%solid_chemistry%reactive_zone%chem_syst%num_redox_eq_reacts) then
                                !else if (ind_cstr>this%solid_chemistry%reactive_zone%chem_syst%num_eq_reacts .AND. ind_cstr<this%solid_chemistry%reactive_zone%chem_syst%num_eq_reacts) then
                                !end if
                                !call this%solid_chemistry%reactive_zone%chem_syst%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                            else
                                error stop
                            end if
                        !else if (icon==5) then
                        !    n_icons(6)=n_icons(6)+1
                        else
                            error stop "icon option not implemented yet"
                        end if
                        this%concentrations(aq_sp_ind)=guess
                        icons(aq_sp_ind)=icon
                        ctots(aq_sp_ind)=ctot
                        !constrains(aq_sp_ind)=constrain
                        k=k+1 !> aqui hay que verificar dimension
                    else
                        error stop 
                    end if
                end if
            end do
            call this%solid_chemistry%reactive_zone%cat_exch_zone%compute_num_surf_compl()
            call this%solid_chemistry%reactive_zone%rearrange_non_flowing_species()
        end if
    do i=1,this%aq_phase%num_species
        call this%aq_phase%aq_species(i)%params_act_coeff%compute_csts(this%aq_phase%aq_species(i)%valence,this%params_aq_sol,model)
    end do
!> We set speciation algebra
    call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_comp(.false.)
    if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
        flag_surf=.true.
    else
        flag_surf=.false.
    end if
    call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_cat_exch(flag_surf)
    call this%solid_chemistry%reactive_zone%set_speciation_alg_dimensions()
    call this%solid_chemistry%reactive_zone%set_eq_reactions()
    call this%solid_chemistry%reactive_zone%set_stoich_mat_react_zone()
    call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,cols)
    
    call this%solid_chemistry%allocate_conc_solids()
    call this%solid_chemistry%allocate_activities()
    call this%solid_chemistry%allocate_equivalents()
    call this%solid_chemistry%allocate_log_act_coeffs_solid_chem()
            
    call this%aq_phase%set_ind_diss_solids()
!> Primary concentrations
    c1=this%get_c1()
!> Initial guess c2
    allocate(c2_init(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions),c2_ig(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
!> We compute initial concentrations (in molalities)
    !> aqui podrias usar polimorifsmo
    if (sum(n_icons(2:4))>0) then
        if (flag_surf==.true.) then
            if (Jac_flag==1) then
                call this%initialise_conc_anal_exch(icons,n_icons,indices_constrains,ctots,surf_chem,niter,CV_flag)
            else
                error stop
            end if
        else if (model==0) then
            call this%initialise_conc_anal_ideal(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
        else
            if (Jac_flag==0) then
                call this%initialise_conc_incr_coeff(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
            else if (Jac_flag==1) then
                call this%initialise_conc_anal(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
            else
                error stop
            end if
        end if
    else if (model==0) then
        !call this%compute_c2_from_c1_ideal(c1,c2_init)
    else
        c2_ig=1d-16 !> chapuza
        !call this%compute_c2_from_c1_Picard(c1,c2_ig,c2_init,niter,CV_flag)
    end if
end subroutine