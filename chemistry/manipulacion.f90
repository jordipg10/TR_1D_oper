!subroutine manipulation(this)
!    use chemistry_Lagr_m
!    implicit none
!    class(chemistry_c) :: this
!    
!    integer(kind=4) :: i,j,k,l,nwtype,icon,n_p_aq,gas_ind,min_ind,model
!    integer(kind=4), allocatable :: num_aq_prim_array(:),num_cstr_array(:)
!    character(len=256) :: prim_sp_name,constrain,label,name
!    real(kind=8) :: guess,c_tot,temp,conc
!    logical :: flag,flag_min,flag_comp,flag_surf
!    type(aq_phase_c) :: old_aq_phase
!
!!> We eliminate constant activity species from component matrix and we rearrange aqueous species and equilibrium reactions
!    call this%chem_syst%speciation_alg%set_flag_comp(.true.)
!    if (this%chem_syst%cat_exch%num_surf_compl>0) then
!        flag_surf=.true.
!    else
!        flag_surf=.false.
!    end if
!    call this%chem_syst%speciation_alg%set_flag_cat_exch(flag_surf)
!    call this%chem_syst%speciation_alg%compute_num_prim_species(this%chem_syst%num_min_kin_reacts)
!    call this%chem_syst%speciation_alg%compute_num_aq_var_act_species()
!    call this%chem_syst%rearrange_species()
!    old_aq_phase=this%chem_syst%aq_phase !> chapuza
!    call this%chem_syst%aq_phase%rearrange_aq_species()
!    call this%chem_syst%aq_phase%set_indices_aq_phase()
!    call this%chem_syst%rearrange_eq_reacts()
!    call this%chem_syst%set_stoich_mat()
!    call this%chem_syst%speciation_alg%compute_arrays(this%chem_syst%Se,this%chem_syst%get_eq_csts(),this%CV_params%zero)
!!> Chapuza
!    do i=1,this%num_wat_types
!        !> rearrange cocnentrations and activities
!        call this%wat_types(i)%aq_chem%rearrange_state_vars(old_aq_phase)
!    end do
!    !do i=1,this%num_init_wat_types
!    !    !> rearrange cocnentrations and activities
!    !    call this%init_wat_types(i)%aq_chem%rearrange_state_vars(old_aq_phase)
!    !end do
!    !do i=1,this%num_bd_wat_types
!    !    call this%bd_wat_types(i)%aq_chem%rearrange_state_vars(old_aq_phase)
!    !end do
!end subroutine