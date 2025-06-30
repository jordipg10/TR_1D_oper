!> This subroutine initialises a 1D transient transport object for a RT simulation
!! It reads:
!!      discretisation parameters
!!      transport properties
!!      BCs
!! It computes stability parameters (in the explicit case)
subroutine initialise_transport_1D_transient_RT(this,root)
    use BCs_m, only: BCs_t
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c, time_discr_c
    use transport_stab_params_m, only: stab_params_tpt_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c, spatial_discr_c
    use transport_transient_m, only: transport_1D_transient_c
    use transport_properties_heterog_m, only: tpt_props_heterog_c
    use char_params_m, only: char_params_c
    use char_params_tpt_m, only: char_params_tpt_c
    use vectors_m, only: inf_norm_vec_real
    implicit none

    class(transport_1D_transient_c) :: this
    character(len=*), intent(in) :: root
    
    type(tpt_props_heterog_c) :: my_props_tpt
    class(spatial_discr_c), pointer :: my_mesh=>null()
    type(mesh_1D_Euler_homog_c), target :: my_homog_mesh
    type(mesh_1D_Euler_heterog_c), target :: my_heterog_mesh
    class(time_discr_c), pointer :: my_time_discr=>null()
    type(time_discr_homog_c), target :: my_homog_time_discr
    type(time_discr_heterog_c), target :: my_heterog_time_discr
    type(BCs_t) :: my_BCs
    type(stab_params_tpt_c) :: my_stab_params_tpt
    type(char_params_tpt_c) :: my_char_params_tpt
    
    real(kind=8) :: q0,Delta_t,Delta_x,theta,measure,Final_time,x_1,x_2,x_3,x
    real(kind=8), allocatable :: c0(:),c_e(:),source_term_vec(:),porosity_vec(:),dispersion_vec(:),flux_vec(:),flux_coeffs(:)
    real(kind=8), allocatable :: d(:),e(:)
    integer(kind=4) :: parameters_flag,i,Num_cells,Num_time,info,n,adapt_ref_flag,scheme,int_method,r_flag,flux_ord,half_num_tar
    real(kind=8), parameter :: pi=4d0*atan(1d0), eps=1d-12, epsilon_x=1d-2, epsilon_t=1d-6
    character(len=200) :: filename
    logical :: evap,dimless
!****************************************************************************************************************************************************
!> Dimensionless form flag
    dimless=.false. !> esto habria que leerlo
    this%dimless=dimless
 !> Uniform mesh
    !my_mesh=>my_homog_mesh
    allocate(mesh_1D_Euler_homog_c :: my_mesh)
    call my_mesh%read_mesh(root//'_discr_esp.dat')
    call this%set_spatial_discr(my_mesh)
!> Boundary conditions
    call my_BCs%read_BCs(root//'_BCs.dat')
    if (my_BCs%BCs_label(1).eq.1 .and. my_BCs%BCs_label(2).eq.1 .and. this%spatial_discr%targets_flag.eq.0) then
        call my_BCs%read_Dirichlet_BCs(root//"_Dirichlet_BCs.dat")
        call my_BCs%read_flux_inf(root//"_flux_inflow.dat")
    else if (my_BCs%BCs_label(1).eq.3) then
        call my_BCs%read_Robin_BC_inflow(root//"_Robin_BC_inflow.dat")
    end if
    call this%set_BCs(my_BCs)
!> Uniform time discretisation
    !my_time_discr=>my_homog_time_discr
    allocate(time_discr_homog_c :: my_time_discr)
    call my_time_discr%read_time_discr(root//'_discr_temp.dat')
    call this%set_time_discr(my_time_discr)
!****************************************************************************************************************************************************
!> Transport properties
    call my_props_tpt%read_props(root//'_tpt_props.dat',this%spatial_discr)
    if (my_props_tpt%source_term_order.eq.0) then !> constant source term
        if (my_props_tpt%cst_flux_flag .eqv. .true.) then !> flux is constant
            call this%BCs%set_cst_flux_boundary(my_props_tpt%flux(1))
        else !> flux is linear
            call my_props_tpt%compute_flux_lin(this%BCs%flux_inf,this%spatial_discr,this%BCs%flux_out)
            !print *, my_props_tpt%flux
        end if
    else if (my_props_tpt%source_term_order>0) then !> flux is polynomic
    !> chapuza
        open(unit=2,file=root//"_flux_coeffs.dat",status='old',action='read')
        read(2,*) flux_ord
        allocate(flux_coeffs(flux_ord+1))
        read(2,*) flux_coeffs
        close(2)
        call my_props_tpt%set_source_term_order(flux_ord-1)
        call my_props_tpt%compute_flux_nonlin(flux_coeffs,this%spatial_discr,this%BCs%flux_out)
        call my_props_tpt%compute_source_term(this%spatial_discr,flux_coeffs)
    end if
    call my_props_tpt%set_source_term_flag(this%BCs)
    call this%set_tpt_props_heterog_obj(my_props_tpt)
!****************************************************************************************************************************************************
!> Stability parameters
    call my_stab_params_tpt%compute_stab_params(this%tpt_props_heterog,my_mesh%get_mesh_size(),my_time_discr%get_Delta_t())
    call this%set_stab_params_tpt(my_stab_params_tpt)
!**************************************************************************************************************************************************
    nullify(my_mesh,my_time_discr)
end subroutine