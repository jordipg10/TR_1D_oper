!> PDE module
!> $Tc_D+g=0$
!> $g=Ec_rech+Bc_BD$
module PDE_m
    use BCs_m, only: BCs_t
    use metodos_sist_lin_m
    use spatial_discr_1D_m, only: spatial_discr_c, mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    implicit none
    save
    type, public, abstract :: PDE_1D_c !> 1D PDE superclass
        class(spatial_discr_c), pointer :: spatial_discr !> spatial discretisation (polymorphic variable)
        type(BCs_t) :: BCs !> Boundary conditions
        type(tridiag_matrix_c) :: trans_mat !> Transition matrix (T) (tridiagonal)
        type(diag_matrix_c) :: rech_mat !> Recharge matrix (E) (diagonal)
        real(kind=8), allocatable :: bd_mat(:) !> Boundary matrix (B) (dim<=2)
        real(kind=8), allocatable :: source_term_PDE(:) !> g
        logical :: dimensionless !> Dimensionless PDE
        integer(kind=4) :: sol_method   !> 1: Numerical
                                        !> 2: Eigendecomposition
    contains
    !> set
        procedure, public :: set_spatial_discr
        procedure, public :: set_BCs
        procedure, public :: set_sol_method
    !> Allocate
        procedure, public :: allocate_trans_mat
        procedure, public :: allocate_rech_mat
        procedure, public :: allocate_bd_mat
        procedure, public :: allocate_source_term_PDE
        procedure, public :: allocate_arrays_PDE_1D=>allocate_arrays_PDE_1D_stat
    !> Update
        procedure, public :: update_trans_mat
    !> Initialise
        procedure(initialise_PDE), public, deferred :: initialise_PDE
    !> Computations
        procedure(compute_trans_mat_PDE), public, deferred :: compute_trans_mat_PDE
        procedure(write_PDE_1D), public, deferred :: write_PDE_1D
        procedure, public :: compute_source_term_PDE
        procedure, public :: compute_rech_mat_PDE
        procedure, public :: solve_PDE_1D
        procedure, public :: solve_write_PDE_1D
        procedure, public :: main_PDE
        procedure, public :: solve_PDE_1D_stat
    end type
!****************************************************************************************************************************************************
    abstract interface
        subroutine compute_trans_mat_PDE(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this
        end subroutine
        
        subroutine initialise_PDE(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this
        end subroutine
        
        
        subroutine write_PDE_1D(this,Time_out,output)
            import PDE_1D_c
            class(PDE_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
        end subroutine
    end interface
!****************************************************************************************************************************************************
    interface
        subroutine compute_source_term_PDE(this,k)
            import PDE_1D_c
            class(PDE_1D_c) :: this
            integer(kind=4), intent(in), optional :: k
        end subroutine
        
        subroutine compute_rech_mat_PDE(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this
            !integer(kind=4), intent(in), optional :: k
        end subroutine
        
        subroutine solve_PDE_1D(this,Time_out,output)
            import PDE_1D_c
            class(PDE_1D_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine solve_PDE_1D_stat(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this
        end subroutine
        
        subroutine solve_write_PDE_1D(this,Time_out)
            import PDE_1D_c
            class(PDE_1D_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
        end subroutine
        
        subroutine main_PDE(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this
        end subroutine
    end interface
!****************************************************************************************************************************************************
    contains
        subroutine set_spatial_discr(this,spatial_discr_obj)
            implicit none
            class(PDE_1D_c) :: this
            class(spatial_discr_c), intent(in), target :: spatial_discr_obj
            this%spatial_discr=>spatial_discr_obj
        end subroutine 
        
        subroutine set_BCs(this,BCs_obj)
            implicit none
            class(PDE_1D_c) :: this
            class(BCs_t), intent(in) :: BCs_obj
            this%BCs=BCs_obj
        end subroutine
        
        subroutine allocate_trans_mat(this)
            implicit none
            class(PDE_1D_c) :: this
            call this%trans_mat%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_rech_mat(this)
            implicit none
            class(PDE_1D_c) :: this
            call this%rech_mat%allocate_array(this%spatial_discr%Num_targets)
            this%rech_mat%diag=0d0 !> Initialize external matrix to zero
        end subroutine
        
        subroutine allocate_bd_mat(this)
            implicit none
            class(PDE_1D_c) :: this
            allocate(this%bd_mat(2))
            this%bd_mat=0d0 !> Initialize boundary matrix to zero
        end subroutine
        
        subroutine update_trans_mat(this,trans_mat)
            implicit none
            class(PDE_1D_c) :: this
            class(tridiag_matrix_c), intent(in)  :: trans_mat
            this%trans_mat=trans_mat
        end subroutine
        
        subroutine allocate_source_term_PDE(this)
            implicit none
            class(PDE_1D_c) :: this
            allocate(this%source_term_PDE(this%spatial_discr%Num_targets))
            this%source_term_PDE=0d0 !> Initialize source term to zero
        end subroutine
        
        subroutine set_sol_method(this,method)
            implicit none
            class(PDE_1D_c) :: this
            integer(kind=4), intent(in) :: method
            if (method<0 .or. method>2) error stop "Solution method not implemented"
            this%sol_method=method
        end subroutine
        
        subroutine allocate_arrays_PDE_1D_stat(this)
            implicit none
            class(PDE_1D_c) :: this
            call this%allocate_trans_mat()
            call this%allocate_rech_mat()
            call this%allocate_bd_mat()
            call this%allocate_source_term_PDE()
        end subroutine

        
end module 