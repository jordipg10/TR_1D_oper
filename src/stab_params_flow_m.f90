module stab_params_flow_m
    use stability_parameters_m, only: stab_params_c
    use properties_m, only: props_c
    use flow_props_heterog_m, only: flow_props_heterog_c
    implicit none
    save
    type, public, extends(stab_params_c) :: stab_params_flow_c !> 1D flow equation stability parameters subclass
        real(kind=8) :: beta !> flow stability parameter (beta=K*Delta_t/S_s*Delta_x^2)
    contains
        procedure, public :: compute_stab_params=>compute_stab_params_flow
    end type
    
    contains
        subroutine compute_stab_params_flow(this,props_obj,mesh_size,time_step)
            implicit none
            class(stab_params_flow_c) :: this
            class(props_c), intent(in) :: props_obj
            real(kind=8), intent(in) :: mesh_size
            real(kind=8), intent(in) :: time_step
                        
            select type (props=>props_obj)
            class is (flow_props_heterog_c)
                this%beta=maxval(props%hydr_cond)*time_step/(minval(props%spec_stor)*mesh_size**2)
            end select
            if (this%beta>5d-1) then
                print *, "Unstable transport", this%beta
                !error stop "Dispersion condition violated"
            end if
        end subroutine
        
        
end module