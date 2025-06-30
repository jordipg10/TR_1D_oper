!> This class contains the partial pressures (activities) of gases
module gas_chemistry_m
    use local_chemistry_m
    use reactive_zone_Lagr_m
    implicit none
    save
!**************************************************************************************************
    type, public, extends(local_chemistry_c) :: gas_chemistry_c !> gas chemistry subclass
        class(reactive_zone_c), pointer :: reactive_zone !> reactive zone containing gas phase
        integer(kind=4), allocatable :: ind_gases_eq_cst_act(:) !> indices of gases in equilibrium with constant activity in chemical system
    contains
    !> Set
        procedure, public :: set_reactive_zone
        procedure, public :: set_indices_gases
        procedure, public :: set_concentrations=>set_conc_gases
    !> Allocate
        procedure, public :: allocate_partial_pressures
        procedure, public :: allocate_conc_gases
        procedure, public :: allocate_log_act_coeffs_gases
        procedure, public :: allocate_ind_gases_eq_cst_act
    !> Compute
        procedure, public :: compute_log_act_coeffs_gases
        procedure, public :: compute_partial_pressures
        procedure, public :: compute_conc_gases_ideal
        procedure, public :: compute_conc_gases_iter
        procedure, public :: compute_vol_gas_conc
        procedure, public :: compute_vol_gas_act_coeffs
        procedure, public :: compute_pressure
    !> Update
        procedure, public :: update_conc_gases
    end type
    
    interface
        !subroutine read_gas_chem_init(this,filename,gas_zones,line,num_tar)
        !>    import gas_zone_c
        !>    import gas_chemistry_c
        !>    implicit none
        !>    class(gas_chemistry_c) :: this
        !>    character(len=*), intent(in) :: filename
        !>    class(gas_zone_c), intent(in) :: gas_zones(:)
        !>    integer(kind=4), intent(inout) :: line
        !>    integer(kind=4), intent(out) :: num_tar
        !end subroutine
        
        !subroutine compute_conc_gases_iter(this,r_vec,Delta_t)
        !    import gas_chemistry_c
        !    implicit none
        !    class(gas_chemistry_c) :: this
        !    real(kind=8), intent(in) :: r_vec(:) !> reaction rates
        !    !real(kind=8), intent(in) :: porosity
        !    real(kind=8), intent(in) :: Delta_t !> time step
        !end subroutine
        
        
    end interface
    
    
    
    contains
        
        !subroutine set_partial_pressures(this,partial_pressures)
        !>    implicit none
        !>    class(gas_chemistry_c) :: this
        !>    real(kind=8), intent(in) :: partial_pressures(:)
        !>    !if (this%gas_zone%num_species/=size(partial_pressures)) error stop "Dimension error in set_partial_pressures"
        !>    this%activities=partial_pressures
        !end subroutine
        
        subroutine allocate_partial_pressures(this) !< units are atm
            implicit none
            class(gas_chemistry_c) :: this
            allocate(this%activities(this%reactive_zone%gas_phase%num_species))
        end subroutine
        
        subroutine allocate_conc_gases(this) !> units are moles
            implicit none
            class(gas_chemistry_c) :: this
            if (allocated(this%concentrations)) then
                deallocate(this%concentrations)
            end if
            allocate(this%concentrations(this%reactive_zone%gas_phase%num_species))
        end subroutine
        
        subroutine allocate_log_act_coeffs_gases(this) !> 
            implicit none
            class(gas_chemistry_c) :: this
            allocate(this%log_act_coeffs(this%reactive_zone%gas_phase%num_species))
        end subroutine
        !
        !subroutine set_gas_zone(this,gas_zone)
        !>    implicit none
        !>    class(gas_chemistry_c) :: this
        !>    class(gas_zone_c), intent(in), target :: gas_zone
        !>    this%gas_zone=>gas_zone
        !end subroutine
        
        subroutine set_reactive_zone(this,reactive_zone)
            implicit none
            class(gas_chemistry_c) :: this
            class(reactive_zone_c), intent(in), target :: reactive_zone
            this%reactive_zone=>reactive_zone
        end subroutine
        
        subroutine set_indices_gases(this)
            implicit none
            class(gas_chemistry_c) :: this
            integer(kind=4) :: i,j,k 
            j=0
            k=0
            do i=1,this%reactive_zone%gas_phase%num_species
                if (this%reactive_zone%gas_phase%gases(i)%cst_act_flag==.false.) then
                    j=j+1
                    this%var_act_species_indices(j)=i
                else
                    k=k+1
                    this%cst_act_species_indices(k)=i
                end if
            end do
        end subroutine
        
       subroutine compute_log_act_coeffs_gases(this)
            implicit none
            class(gas_chemistry_c) :: this
            integer(kind=4) :: i
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            do i=1,this%reactive_zone%gas_phase%num_species
                this%log_act_coeffs(i)=log10(R*this%temp/this%volume)
            end do
       end subroutine
       
       subroutine compute_partial_pressures(this)
            implicit none
            class(gas_chemistry_c) :: this
            integer(kind=4) :: i
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            do i=1,this%reactive_zone%gas_phase%num_species
                this%activities(i)=this%concentrations(i)*R*this%temp/this%volume
            end do
       end subroutine
       
       subroutine compute_conc_gases_ideal(this) !> ideal gas equation
       !> Concentrations are expressed in moles
            implicit none
            class(gas_chemistry_c) :: this
            integer(kind=4) :: i
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            do i=1,this%reactive_zone%gas_phase%num_species
                this%concentrations(i)=this%activities(i)*this%volume/(this%temp*R)
            end do
       end subroutine
       
       subroutine compute_conc_gases_iter(this,Delta_t,porosity,wat_vol,r_aq) !> gas conservation equation
       !> Concentrations are expressed in moles
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: porosity !> porosity
            real(kind=8), intent(in) :: wat_vol !> water volume
            real(kind=8), intent(in), optional :: r_aq(:)
            integer(kind=4) :: i
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            
            if (present(r_aq)) then
                if (this%reactive_zone%gas_phase%num_gases_kin>0) then
                    do i=1,this%reactive_zone%gas_phase%num_cst_act_species
                        this%concentrations(this%cst_act_species_indices(i))=this%concentrations(this%cst_act_species_indices(i))+Delta_t*dot_product(this%reactive_zone%chem_syst%stoich_mat_gas(:,this%cst_act_species_indices(i)),r_aq)*wat_vol/porosity
                    end do
                else if (this%reactive_zone%gas_phase%num_gases_eq>0) then
                    do i=1,this%reactive_zone%gas_phase%num_cst_act_species
                        this%concentrations(this%cst_act_species_indices(i))=this%concentrations(this%cst_act_species_indices(i))+Delta_t*dot_product(this%reactive_zone%chem_syst%stoich_mat_gas(:,this%cst_act_species_indices(i)),[this%r_eq(i),r_aq])*wat_vol/porosity
                    end do
                end if
            end if
       end subroutine
       
       subroutine compute_vol_gas_conc(this)
       !> Computes total volume of gas from concentrations
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            this%volume=sum(this%concentrations)*R*this%temp/this%pressure !> we choose first gas species in gas phase arbitrarily
       end subroutine
       
       subroutine compute_vol_gas_act_coeffs(this)
       !> Computes total volume of gas from activity coefficients
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            this%volume=R*this%temp/(10**this%log_act_coeffs(1)) !> we choose first gas species in gas phase arbitrarily
       end subroutine
       
       subroutine compute_pressure(this)
       !> Computes total pressure of gas
            implicit none
            class(gas_chemistry_c) :: this
            this%pressure=sum(this%activities) !> we sum the partial pressures
       end subroutine
       
       subroutine update_conc_gases(this,conc_gases)
       !> Updates concentration of gases
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_gases(:)
            this%concentrations=conc_gases
       end subroutine
       
       subroutine set_conc_gases(this,conc)
       !> Sets concentration of gases
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:)
            if (size(conc)/=this%reactive_zone%gas_phase%num_species) error stop "Dimension error in set_conc_gases"
            this%concentrations=conc
       end subroutine
       
        !subroutine initialise_log_act_coeffs_gas_chem(this)
        !    implicit none
        !    class(gas_chemistry_c) :: this
        !    allocate(this%log_act_coeffs(this%reactive_zone%gas_phase%num_species))
        !    this%log_act_coeffs=0d0
        !end subroutine
        !
        !subroutine initialise_log_Jacobian_act_coeffs_gas_chem(this)
        !    implicit none
        !    class(gas_chemistry_c) :: this
        !    allocate(this%log_Jacobian_act_coeffs(this%reactive_zone%gas_phase%num_species,this%reactive_zone%gas_phase%num_species))
        !    this%log_Jacobian_act_coeffs=0d0
        !end subroutine
       
       subroutine allocate_ind_gases_eq_cst_act(this)
           class(gas_chemistry_c) :: this
           allocate(this%ind_gases_eq_cst_act(this%reactive_zone%gas_phase%num_gases_eq_cst_act))
       end subroutine
end module