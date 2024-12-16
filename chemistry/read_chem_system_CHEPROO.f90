!> Lectura sistema quimico basada en CHEPROO
!> Suponemos que el archivo ya ha sido abierto
subroutine read_chem_system_CHEPROO(this,path_DB,unit)
    use chem_system_m
    implicit none
    class(chem_system_c) :: this !> chemical system object
    character(len=*), intent(in) :: path_DB !> path to databases
    integer(kind=4), intent(in) :: unit !> file unit
    
    real(kind=8), allocatable :: Sk(:,:),logK(:),gamma_1(:),gamma_2(:)
    integer(kind=4) :: n_gas_eq_cst_act,n_gas_eq_var_act,n_gas_kin,exch_cat_valence,n_eq_homog,unit_master_25,unit_kinetics,unit_redox,i,j,num_sp,ind_var_act_sp,num_var_act_sp,num_cst_act_sp,k,num_aq_sp,num_sec_aq_sp,exch_cat_ind, n_min_kin, n_gas_eq,index,kin_react_type,n_r,num_mins,num_gases,num_surf_compl,num_exch_cats,n_eq,n_k,num_mins_eq_indices,n_redox,n_lin_kin,n_lin_eq,num_aq_compl,num_mins_eq,num_cst_act_gases,num_var_act_gases
    integer(kind=4) :: num_eq_cst_act_mins,num_var_act_mins,n_redox_eq,n_redox_kin
    integer(kind=4), allocatable :: n_tar(:),mins_eq_indices(:),gases_eq_indices(:),indices_lin_reacts(:,:)
    real(kind=8) :: aux,conc,temp,SI,lambda
    character(len=256) :: str,str1,str2,str3,str4,str5,Monod_name,file_kin_params,label
    character(len=:), allocatable :: str_block_trim,str_trim,valence_str,exch_cat_val,exch_cat_name,str_trim_1,str_trim_2
    logical :: flag,eq_label,exch_cat_flag,cst_act_label,flag_1,flag_2
    
    character(len=256), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:),solid_species_str(:),kin_react_names(:)
    type(species_c) :: species
    type(species_c), allocatable :: surf_compl(:)
    type(aq_species_c) :: exch_cat,aq_sp_1,aq_sp_2
    type(aq_species_c), allocatable :: aq_species(:),exch_cats(:),prim_species(:)
    type(mineral_c) :: mineral
    type(mineral_c), allocatable :: mins(:)
    type(surface_c) :: cat_exch_obj
    type(gas_c) :: gas
    type(gas_c), allocatable :: gases(:)
    class(kin_params_c), pointer :: p_kin_params=>null()
    class(kin_reaction_c), pointer :: p_kin_react=>null()
    type(kin_reaction_poly_c) :: kin_react_ptr
    class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
    type(eq_reaction_c) :: eq_react
    type(eq_reaction_c), allocatable :: eq_reacts(:)
    
    
    logical :: flag_comp,flag_surf
!> We initialise counters
    num_sp=0
    num_aq_sp=0
    num_aq_compl=0
    num_cst_act_sp=0
    num_cst_act_gases=0
    num_var_act_gases=0
    num_eq_cst_act_mins=0
    !num_var_act_mins=0
    num_var_act_sp=0
    num_surf_compl=0
    num_exch_cats=0
    num_gases=0
    num_mins=0
    num_mins_eq=0
    n_eq=0
    n_k=0
    n_min_kin=0
    n_gas_eq=0
    n_gas_eq_cst_act=0
    n_gas_eq_var_act=0
    n_gas_kin=0
    n_lin_kin=0
    n_lin_eq=0
    n_redox_eq=0
    n_redox_kin=0
!> First iteration of chemical system file
    do
        read(unit,*) label
        if (label=='end') then
            rewind(unit)
            exit
        else if (label=='PRIMARY AQUEOUS SPECIES') then !> suponemos ordenadas en primarias y secundarias
            do
                read(unit,*) str
                if (str=='*') then
                    exit
                else if (str=='h2o') then
                !> we assume water has constant activity (approx 1)
                    num_cst_act_sp=num_cst_act_sp+1
                    this%aq_phase%wat_flag=1
                    this%aq_phase%ind_wat=num_aq_sp+1
                else
                    num_var_act_sp=num_var_act_sp+1
                end if
                num_aq_sp=num_aq_sp+1
                num_sp=num_sp+1
            end do
        else if (label=='AQUEOUS COMPLEXES') then
            do
                read(unit,*) str
                if (str=='*') then
                    exit
                else if (str=='h2o') then
                !> we assume water has constant activity (approx 1)
                    num_cst_act_sp=num_cst_act_sp+1
                    this%aq_phase%wat_flag=1
                    this%aq_phase%ind_wat=num_aq_sp+1
                else
                    num_var_act_sp=num_var_act_sp+1
                end if
                num_aq_sp=num_aq_sp+1
                num_aq_compl=num_aq_compl+1
                num_sp=num_sp+1
                n_eq=n_eq+1
            end do            
        else if (label=='MINERALS') then
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                num_mins=num_mins+1
                num_sp=num_sp+1
                if (cst_act_label==.true. .AND. eq_label==.true.) then
                    n_eq=n_eq+1
                    num_mins_eq=num_mins_eq+1
                    num_eq_cst_act_mins=num_eq_cst_act_mins+1
                    num_cst_act_sp=num_cst_act_sp+1
                else if (cst_act_label==.false. .AND. eq_label==.true.) then
                    !num_var_act_mins=num_var_act_mins+1
                    num_var_act_sp=num_var_act_sp+1
                else if (cst_act_label==.true. .AND. eq_label==.false.) then
                    num_cst_act_sp=num_cst_act_sp+1
                    n_k=n_k+1
                    n_min_kin=n_min_kin+1
                else
                    !num_var_act_mins=num_var_act_mins+1
                    num_var_act_sp=num_var_act_sp+1
                    n_k=n_k+1
                    n_min_kin=n_min_kin+1
                end if
            end do 
        else if (label=='GASES') then
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                num_gases=num_gases+1
                num_sp=num_sp+1
                if (cst_act_label==.true. .AND. eq_label==.true.) then 
                    n_gas_eq=n_gas_eq+1
                    n_gas_eq_cst_act=n_gas_eq_cst_act+1
                    num_cst_act_gases=num_cst_act_gases+1
                    num_cst_act_sp=num_cst_act_sp+1
                else if (cst_act_label==.false. .AND. eq_label==.true.) then 
                    n_gas_eq=n_gas_eq+1
                    n_gas_eq_var_act=n_gas_eq_var_act+1
                    num_var_act_gases=num_var_act_gases+1
                    num_var_act_sp=num_var_act_sp+1
                else if (cst_act_label==.true. .AND. eq_label==.false.) then 
                    n_gas_kin=n_gas_kin+1
                    num_cst_act_gases=num_cst_act_gases+1
                    num_cst_act_sp=num_cst_act_sp+1
                else
                    n_gas_kin=n_gas_kin+1
                    num_var_act_gases=num_var_act_gases+1
                    num_var_act_sp=num_var_act_sp+1
                end if
            end do 
        else if (label=='SURFACE COMPLEXES') then
            do
                read(unit,*) str
                if (str=='*') exit
                str_trim=trim(str)
                num_surf_compl=num_surf_compl+1
                num_var_act_sp=num_var_act_sp+1
                num_sp=num_sp+1
                exch_cat_ind=index(str_trim,'-')
                if (exch_cat_ind>0 .and. exch_cat_ind<len(str_trim)) then
                    num_exch_cats=num_exch_cats+1
                    n_eq=n_eq+1
                end if
            end do
        else if (label=='REDOX REACTIONS') then
            do
                read(unit,*) str, eq_label
                if (str=='*') exit
                if (eq_label==.true.) then
                    n_redox_eq=n_redox_eq+1
                    n_eq=n_eq+1
                else
                    n_redox_kin=n_redox_kin+1
                    n_k=n_k+1
                end if
            end do
        else if (label=='LINEAR REACTIONS') then
            do
                read(unit,*) str1, str2, lambda, eq_label
                if (str1=='*') exit
                if (eq_label==.true.) then
                    n_lin_eq=n_lin_eq+1
                    n_eq=n_eq+1
                else
                    n_lin_kin=n_lin_kin+1
                    n_k=n_k+1
                end if
            end do
        else 
            continue
        end if
    end do
!> We set & allocate attributes
    !> Aqueous phase
    call this%aq_phase%allocate_aq_species(num_aq_sp)
    call this%aq_phase%allocate_ind_diss_solids()
    call this%aq_phase%set_num_aq_complexes(num_aq_compl)
    !> Gas phase
    call this%gas_phase%allocate_gases(num_gases)
    call this%gas_phase%set_num_gases_eq(n_gas_eq) !> 
    call this%gas_phase%set_num_gases_eq_cst_act(n_gas_eq_cst_act) !> 
    call this%gas_phase%set_num_gases_eq_var_act(n_gas_eq_cst_act) !> 
    call this%gas_phase%set_num_gases_kin(n_gas_kin) !> 
    call this%gas_phase%set_num_var_act_species_phase(num_var_act_gases)
    call this%gas_phase%set_num_cst_act_species_phase(num_cst_act_gases)
    !> Cation exchange
    call this%cat_exch%allocate_surf_compl(num_surf_compl)
    call this%cat_exch%allocate_exch_cat_indices(num_exch_cats)
    call this%cat_exch%set_num_exch_cats(num_exch_cats)
    !> Chemical system
    call this%allocate_cst_act_sp_indices(num_cst_act_sp)
    call this%allocate_var_act_sp_indices(num_var_act_sp)
    call this%allocate_species()
    call this%allocate_minerals(num_mins)
    call this%set_num_minerals_eq(num_mins_eq)
    call this%set_num_minerals_eq_cst_act(num_eq_cst_act_mins)
    call this%allocate_reacts(n_eq,n_k)
    call this%allocate_min_kin_reacts(n_min_kin)
    call this%set_num_redox_eq_reacts(n_redox_eq)
    call this%allocate_redox_kin_reacts(n_redox_kin)
    call this%allocate_lin_kin_reacts(n_lin_kin)
    call this%compute_num_reacts()
    call this%compute_num_solids()
    !> Speciatrion algebra    
    call this%speciation_alg%set_flag_comp(.false.)
    if (num_surf_compl>0) then
        flag_surf=.true.
    else
        flag_surf=.false.
    end if
    call this%speciation_alg%set_flag_cat_exch(flag_surf)
    call this%speciation_alg%set_dimensions(this%num_species,this%num_eq_reacts,this%num_cst_act_species,this%aq_phase%num_species,this%aq_phase%num_species-this%aq_phase%wat_flag,this%num_min_kin_reacts,num_gases-n_gas_eq)
!> Second iteration of chemical system file
    ind_var_act_sp=0 !< counter variable activity species
    i=0 !> counter aqueous species
    do
        read(unit,*) label
        if (label=='end') then
            exit
        else if (label=='PRIMARY AQUEOUS SPECIES' .OR. label=='AQUEOUS COMPLEXES') then !> suponemos ordenadas en primarias y secundarias
            do
                read(unit,*) str
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                call this%aq_phase%aq_species(i)%set_name(str_trim)
                if (str_trim=='h2o') then
                    call this%aq_phase%aq_species(i)%set_cst_act_flag(.true.)
                    !call this%species(this%num_var_act_species+1)%assign_species(this%aq_phase%aq_species(i))
                    this%cst_act_sp_indices(1)=i
                else
                    ind_var_act_sp=ind_var_act_sp+1
                    call this%aq_phase%aq_species(i)%set_cst_act_flag(.false.)
                    !call this%species(ind_var_act_sp)%assign_species(this%aq_phase%aq_species(i))
                    this%var_act_sp_indices(ind_var_act_sp)=i
                end if
                call this%species(i)%assign_species(this%aq_phase%aq_species(i))
            end do
        else if (label=='MINERALS') then
            !allocate(mins_eq_indices(0))
            i=0 !> counter minerals equilibrium
            j=0 !> counter minerals kinetic
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                str_trim=trim(str)
                if (eq_label==.true.) then
                    i=i+1
                    call this%minerals(this%num_min_kin_reacts+i)%set_phase_name(str_trim)
                    call this%minerals(this%num_min_kin_reacts+i)%mineral%set_name(str_trim)
                    call this%minerals(this%num_min_kin_reacts+i)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(this%num_min_kin_reacts+i)%mineral%set_valence(0)
                else
                    j=j+1
                    call this%minerals(j)%set_phase_name(str_trim)
                    call this%minerals(j)%mineral%set_name(str_trim)
                    call this%minerals(j)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(j)%mineral%set_valence(0)
                end if
            end do
        else if (label=='GASES') then
            !allocate(gases_eq_indices(0))
            i=0 !> counter gases eq
            j=0 !> counter gases kin
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                str_trim=trim(str)
                if (eq_label==.true.) then
                    i=i+1
                    call this%gas_phase%gases(i)%set_name(str_trim)
                    call this%gas_phase%gases(i)%set_cst_act_flag(cst_act_label)
                    call this%gas_phase%gases(i)%set_valence(0)
                else
                    j=j+1
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq+j)%set_name(str_trim)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq+j)%set_cst_act_flag(cst_act_label)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq+j)%set_valence(0)
                end if
            end do
        else if (label=='SURFACE COMPLEXES') then
            i=0 !> counter surface complexes
            j=0 !> counter exchangeable cations
            do
                read(unit,*) str
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                call this%cat_exch%surf_compl(i)%set_name(str_trim)
                call this%cat_exch%surf_compl(i)%set_cst_act_flag(.false.)
                exch_cat_ind=index(str_trim,'-')
                if (exch_cat_ind>0 .and. exch_cat_ind<len(str_trim)) then
                    j=j+1
                    if (exch_cat_ind>1) then
                        if (exch_cat_ind==2) then
                            exch_cat_valence=1
                            allocate(character(len=(len(str_trim)-exch_cat_ind+1)) :: exch_cat_name)
                            write(exch_cat_name(1:len(str_trim)-exch_cat_ind),"(A)") str_trim(exch_cat_ind+1:len(str_trim))
                            write(exch_cat_name(len(str_trim)-exch_cat_ind+1:len(str_trim)-exch_cat_ind+1),"(A1)") '+'
                        else
                            exch_cat_val=trim(str(2:exch_cat_ind-1))
                            read(exch_cat_val,*) exch_cat_valence
                            allocate(character(len=(len(str_trim)-exch_cat_ind+1+len(exch_cat_val))) :: exch_cat_name)
                            write(exch_cat_name(1:len(str_trim)-exch_cat_ind),"(A)") str_trim(exch_cat_ind+1:len(str_trim))
                            write(exch_cat_name(len(str_trim)-exch_cat_ind+1:len(str_trim)-exch_cat_ind+1),"(A1)") '+'
                            write(exch_cat_name(len(exch_cat_name)-len(exch_cat_val)+1:len(exch_cat_name)),"(A)") exch_cat_val
                        end if
                    end if
                    call exch_cat%set_name(exch_cat_name)
                    call exch_cat%set_valence(exch_cat_valence)
                    call this%aq_phase%is_species_in_aq_phase(exch_cat,exch_cat_flag,exch_cat_ind)
                    if (exch_cat_flag==.true.) then
                        this%cat_exch%exch_cat_indices(j)=exch_cat_ind
                    else
                        error stop "Exchangeable cation is not in the chemical system"
                    end if
                    deallocate(exch_cat_name)
                end if
                
            end do
        else if (label=='REDOX REACTIONS') then
            i=0 !> counter redox kinetic reactions
            j=0 !> counter redox equilibrium reactions
            do
                read(unit,*) str, eq_label
                if (str=='*') exit
                str_trim=trim(str)
                if (eq_label==.false.) then
                    i=i+1
                    call this%redox_kin_reacts(i)%set_react_name(str_trim)
                    call this%redox_kin_reacts(i)%set_react_type(4)
                    call this%kin_reacts(this%num_lin_kin_reacts+I)%set_kin_reaction(this%redox_kin_reacts(i))
                else
                    j=j+1
                    call this%eq_reacts(j)%set_react_name(str_trim)
                    call this%eq_reacts(j)%set_react_type(4)
                end if
            end do
        else if (label=='LINEAR REACTIONS') then
            i=0 !> counter linear kinetic reactions
            j=0 !> counter linear equilibrium reactions
            allocate(indices_lin_reacts(n_lin_kin+n_lin_eq,2)) !> first kinetic, then equilibrium
            do
                read(unit,*) str1, str2, lambda, eq_label
                if (str1=='*') exit
                str_trim_1=trim(str1)
                str_trim_2=trim(str2)
                call aq_sp_1%set_name(str_trim_1)
                call aq_sp_2%set_name(str_trim_2)
                if (eq_label==.false.) then
                    i=i+1
                    call this%aq_phase%is_species_in_aq_phase(aq_sp_1,flag_1,indices_lin_reacts(i,1))
                    call this%aq_phase%is_species_in_aq_phase(aq_sp_2,flag_2,indices_lin_reacts(i,2))
                    if (flag_1==.false. .or. flag_2==.false.) then
                        error stop
                    else
                        call this%lin_kin_reacts(i)%set_react_name(str_trim_1//'-->'//str_trim_2)
                        call this%lin_kin_reacts(i)%set_react_type(5)
                        call this%lin_kin_reacts(i)%allocate_reaction(2)
                        call this%lin_kin_reacts(i)%set_lambda(lambda)
                        call this%kin_reacts(i)%set_kin_reaction(this%lin_kin_reacts(i))
                    end if
                else
                    !j=j+1
                    !call this%eq_reacts(j)%set_react_name(str_trim_1//str_trim_2)
                    !call this%eq_reacts(j)%set_react_type(5)
                end if
            end do
        else 
            continue
        end if
    end do
!> We read databases
    !> Master25 (equilibrium reactions)
    unit_master_25=2
    call this%read_master25(path_DB,unit_master_25)
    !> Kinetics (minerals)
    unit_kinetics=3
    if (this%num_min_kin_reacts>0) then
        call this%read_kinetics_DB(path_DB,unit_kinetics)
    end if
    !> Monod reactions
    unit_redox=4
    if (this%num_redox_kin_reacts>0 .or. this%num_redox_eq_reacts>0) then
        call this%read_Monod_reacts(path_DB,unit_redox)
    end if
    !> Linear kinetcir eactions
    if (this%num_lin_kin_reacts>0) then
        do i=1,this%num_lin_kin_reacts
            call this%lin_kin_reacts(i)%set_all_species(this%aq_phase%aq_species(indices_lin_reacts(i,:)))
            call this%lin_kin_reacts(i)%set_stoichiometry([-1d0,1d0])
        end do
    end if
!> We rearrange species and equilibrium reactions
    call this%rearrange_species()
    call this%compute_z2()
    call this%rearrange_eq_reacts()
!> We set stoichiometric matrices
    call this%set_stoich_mat()
    call this%set_stoich_mat_gas()
    call this%set_stoich_mat_sol()
end subroutine