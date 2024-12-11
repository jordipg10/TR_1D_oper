!> Writes data and results of chemistry object
subroutine write_chemistry(this,unit)
    use chemistry_Lagr_m
    
    implicit none
    class(chemistry_c), intent(in) :: this !> chemistry object
    integer(kind=4), intent(in) :: unit !> file unit
    
    integer(kind=4) :: i,j,n,l,k,num_dom,num_ext
    integer(kind=4), allocatable :: dom_indices(:),ext_indices(:)
    

    write(unit,"(/,2x,'Aqueous species:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,A15)") this%chem_syst%aq_phase%aq_species(i)%name
    end do
    write(unit,"(2x,'Minerals (name + molar volume):'/)")
    do i=1,this%chem_syst%num_minerals
        write(unit,"(10x,A15,ES15.5)") this%chem_syst%minerals(i)%name, this%chem_syst%minerals(i)%mineral%mol_vol
    end do
    write(unit,"(2x,'Gases:'/)")
    do i=1,this%chem_syst%gas_phase%num_species
        write(unit,"(10x,A15)") this%chem_syst%gas_phase%gases(i)%name
    end do
    write(unit,"(2x,'Surface complexes:'/)")
    do i=1,this%chem_syst%cat_exch%num_surf_compl
        write(unit,"(10x,A15)") this%chem_syst%cat_exch%surf_compl(i)%name
    end do
    write(unit,"(2x,'Equilibrium reactions:'/)")
    do i=1,this%chem_syst%num_eq_reacts
        write(unit,"(10x,A30/)") this%chem_syst%eq_reacts(i)%name
    end do
    write(unit,"(2x,'Kinetic reactions:'/)")
    do i=1,this%chem_syst%num_kin_reacts
        call this%chem_syst%kin_reacts(i)%kin_reaction%write_reaction(unit)
    end do
    write(unit,"(2x,'Global stoichiometric matrix:',/)")
    do i=1,this%chem_syst%num_reacts
        write(unit,"(10x,*(F15.5))") (this%chem_syst%stoich_mat(i,j), j=1,this%chem_syst%num_species)
    end do
    write(unit,"(2x,'Equilibrium constants:',/)")
    do i=1,this%chem_syst%num_eq_reacts
        write(unit,"(10x,*(ES15.5))") this%chem_syst%eq_reacts(i)%eq_cst
    end do
    write(unit,"(/,2x,'Global component matrix:'/)")
    do i=1,this%chem_syst%speciation_alg%num_prim_species
        write(unit,"(10x,*(F15.5))") (this%chem_syst%speciation_alg%comp_mat(i,j), j=1,this%chem_syst%speciation_alg%num_var_act_species)
    end do
    write(unit,"(/,2x,'Concentration of external waters: (rows -> aqueous species, columns -> targets)'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%ext_waters_indices(j))%concentrations(i), j=1,this%num_ext_waters)
    end do
    write(unit,"(/,2x,'Initial concentration of aqueous species in the domain:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%concentrations(this%target_waters_init(j)%indices_aq_phase(i)), j=1,this%num_target_waters_dom)
    end do
    write(unit,"(/,2x,'Initial activity coefficients aqueous species domain:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (10**this%target_waters_init(this%dom_tar_wat_indices(j))%log_act_coeffs(this%target_waters_init(j)%indices_aq_phase(i)), j=1,this%num_target_waters_dom)
    end do
    if (this%chem_syst%aq_phase%ind_wat>0) then
        write(unit,"(/,2x,'Initial activity water domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%activities(this%chem_syst%aq_phase%ind_wat), j=1,this%num_target_waters_dom)
        write(unit,"(/,2x,'Initial salinity domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%salinity, j=1,this%num_target_waters_dom)
        write(unit,"(/,2x,'Activity water domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%activities(this%chem_syst%aq_phase%ind_wat), j=1,this%num_target_waters_dom)
        write(unit,"(/,2x,'Final salinity domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%salinity, j=this%num_ext_waters+1,this%num_target_waters_dom)
    end if
    if (this%chem_syst%aq_phase%ind_proton>0) then
        write(unit,"(/,2x,'Initial pH domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%pH, j=1,this%num_target_waters_dom)
        write(unit,"(/,2x,'Final pH domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%pH, j=1,this%num_target_waters_dom)
    end if
    write(unit,"(/,2x,'Aqueous equilibrium reaction rates:'/)")
    do i=1,this%chem_syst%aq_phase%num_aq_complexes
        write(unit,"(10x,*(ES15.5))") (this%target_waters(j)%r_eq(i), j=1,this%num_target_waters_dom)
    end do
    write(unit,"(/,2x,'Kinetic reaction rates:'/)")
    do i=1,this%chem_syst%num_kin_reacts
        write(unit,"(10x,*(ES15.5))") (this%target_waters(j)%rk(i), j=1,this%num_target_waters_dom)
    end do
    do l=1,this%num_reactive_zones
        call this%link_target_waters_reactive_zone(l,dom_indices,ext_indices)
        num_dom=size(dom_indices)
        num_ext=size(ext_indices)
        if (num_dom>0) then
            write(unit,"(/,2x,'Reactive zone',I5,':'/)") l
            write(unit,"(10x,'Non-flowing species:'/)")
            do i=1,this%reactive_zones(l)%num_non_flowing_species
                write(unit,"(20x,A15/)") this%reactive_zones(l)%non_flowing_species(i)%name
            end do
            write(unit,"(10x,'Number of equilibrium reactions:',I10/)") this%reactive_zones(l)%speciation_alg%num_eq_reactions
            write(unit,"(10x,'Se:',/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_eq_reactions
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%stoich_mat(i,j), j=1,this%reactive_zones(l)%speciation_alg%num_species)
            end do
            write(unit,"(10x,'Sk:',/)")
            do i=1,this%reactive_zones(l)%chem_syst%num_kin_reacts
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%chem_syst%Sk(i,j), j=1,this%reactive_zones(l)%chem_syst%num_species)
            end do
            write(unit,"(10x,'U*Sk^T:',/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_aq_prim_species
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%U_SkT_prod(i,j), j=1,this%reactive_zones(l)%chem_syst%num_kin_reacts)
            end do
            write(unit,"(/,10x,'Equilibrium constants:',/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_eq_reactions
                write(unit,"(10x,*(ES15.5))") this%reactive_zones(l)%eq_reactions(i)%eq_cst
            end do
            write(unit,"(/,10x,'Aqueous species:',/)")
            do i=1,this%chem_syst%aq_phase%num_species
                write(unit,"(20x,A15)") this%reactive_zones(l)%chem_syst%aq_phase%aq_species(this%target_waters(dom_indices(1))%indices_aq_phase(i))%name
            end do
            write(unit,"(/,10x,'Primary species:',/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_aq_prim_species
                write(unit,"(20x,A15)") this%reactive_zones(l)%chem_syst%aq_phase%aq_species(this%target_waters(dom_indices(1))%indices_aq_phase(i))%name
            end do
            write(unit,"(10x,'Component matrix:'/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_prim_species
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%speciation_alg%comp_mat(i,j), j=1,this%reactive_zones(l)%speciation_alg%num_var_act_species)
            end do           
            write(unit,"(/,10x,'Initial concentration of components:'/)")
                do i=1,this%reactive_zones(l)%speciation_alg%num_prim_species
                    write(unit,"(10x,*(ES15.5))") (dot_product(this%reactive_zones(l)%speciation_alg%comp_mat(i,:),this%target_waters_init(dom_indices(j))%get_conc_nc()), j=1,size(dom_indices))
                end do
            write(unit,"(/,10x,'Final concentration of components:'/)")
                do i=1,this%reactive_zones(l)%speciation_alg%num_prim_species
                    write(unit,"(10x,*(ES15.5))") (dot_product(this%reactive_zones(l)%speciation_alg%comp_mat(i,:),this%target_waters(dom_indices(j))%get_conc_nc()), j=1,size(dom_indices))
                end do
            write(unit,"(/,10x,'Final concentration of aqueous species:'/)")
            do i=1,this%chem_syst%aq_phase%num_species
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%concentrations(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Final activity coefficients aqueous species:'/)")
            do i=1,this%chem_syst%aq_phase%num_species
                write(unit,"(10x,*(ES15.5))") (10**this%target_waters(dom_indices(j))%log_act_coeffs(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Initial volumetric fractions of minerals:'/)")
            do i=1,this%reactive_zones(l)%num_minerals
                write(unit,"(10x,*(ES15.5))") (this%target_solids_init(dom_indices(j))%vol_fracts(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Volumetric fractions of minerals:'/)")
            do i=1,this%reactive_zones(l)%num_minerals
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%vol_fracts(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Initial concentration of minerals:'/)")
            do i=1,this%reactive_zones(l)%num_minerals
                write(unit,"(10x,*(ES15.5))") (this%target_solids_init(dom_indices(j))%concentrations(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Concentration of minerals:'/)")
            do i=1,this%reactive_zones(l)%num_minerals
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%concentrations(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Concentration of surface complexes:'/)")
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_surf_compl
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%concentrations(this%reactive_zones(l)%num_minerals+i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Initial concentration of gases:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_species
                write(unit,"(10x,*(ES15.5))") (this%target_gases_init(dom_indices(j))%concentrations(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Concentration of gases:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_species
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%gas_chemistry%concentrations(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Initial partial pressures of gases:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_species
                write(unit,"(10x,*(ES15.5))") (this%target_gases_init(dom_indices(j))%activities(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Partial pressures of gases:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_species
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%gas_chemistry%activities(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Initial volume of gas:'/)")
            if (this%reactive_zones(l)%gas_phase%num_species>0) then
                write(unit,"(10x,*(ES15.5))") (this%target_gases_init(dom_indices(j))%volume, j=1,size(dom_indices))
            end if
            write(unit,"(/,10x,'Volume of gas:'/)")
            if (this%reactive_zones(l)%gas_phase%num_species>0) then
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%gas_chemistry%volume, j=1,size(dom_indices))
            end if
            write(unit,"(/,10x,'Aqueous equilibrium reaction rates:'/)")
            do i=1,this%reactive_zones(l)%chem_syst%aq_phase%num_aq_complexes+this%reactive_zones(l)%chem_syst%num_redox_eq_reacts
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%r_eq(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Mineral equilibrium reaction rates:'/)")
            do i=1,this%reactive_zones(l)%num_minerals
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%r_eq(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Gas equilibrium reaction rates:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_species
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%gas_chemistry%r_eq(i), j=1,size(dom_indices))
            end do
            deallocate(dom_indices,ext_indices)
        end if
    end do
end subroutine 