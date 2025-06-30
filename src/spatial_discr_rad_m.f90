!> Radial spatial discretisation module
module spatial_discr_rad_m
    use spatial_discr_m
    implicit none
    save
    type, public, extends(spatial_discr_c) :: spatial_discr_rad_c
        integer(kind=4) :: dim !> dimension
        real(kind=8) :: r_max
        real(kind=8) :: r_min
        real(kind=8), allocatable :: Delta_r(:)
    contains
        procedure, public :: set_dim
        procedure, public :: set_r_max
        procedure, public :: set_r_min
        procedure, public :: set_Delta_r
        procedure, public :: read_mesh=>read_mesh_rad
        procedure, public :: get_mesh_size=>get_Delta_r
        procedure, public :: get_dim=>get_dim_rad
        procedure, public :: compute_r_max
        procedure, public :: compute_Delta_r
        procedure, public :: compute_measure=>compute_measure_rad
        procedure, public :: refine_mesh=>refine_mesh_rad
        procedure, public :: compute_dimless_mesh=>compute_dimless_mesh_rad
    end type
    
    interface
        subroutine refine_mesh_rad(this,conc,conc_ext,rel_tol)
            import spatial_discr_rad_c
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol !> relative tolerance
        end subroutine
    end interface
    
    contains
        subroutine set_dim(this,dim)
            implicit none
            class(spatial_discr_rad_c) :: this
            integer(kind=4), intent(in) :: dim
            this%dim=dim
        end subroutine
        
        subroutine set_r_max(this,r_max)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), intent(in) :: r_max
            this%r_max=r_max
        end subroutine
        
        subroutine set_r_min(this,r_min)
        implicit none
        class(spatial_discr_rad_c) :: this
        real(kind=8), intent(in) :: r_min
        this%r_min=r_min
        end subroutine
        
        subroutine set_Delta_r(this,Delta_r)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), intent(in) :: Delta_r(:)
            this%Delta_r=Delta_r
        end subroutine
        
        
        
        subroutine read_mesh_rad(this,filename)
            implicit none
            class(spatial_discr_rad_c) :: this
            character(len=*), intent(in) :: filename
            open(unit=1,file=filename,status='old',action='read')
            read(1,*) this%scheme
            read(1,*) this%targets_flag
            read(1,*) this%Num_targets
            read(1,*) this%r_min
            read(1,*) this%r_max
            read(1,*) this%adapt_ref
            this%Num_targets_defined=.true.
        end subroutine
        
        function get_Delta_r(this,i) result(Delta_r)
            implicit none
            class(spatial_discr_rad_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: Delta_r
            if (present(i)) then
                Delta_r=this%Delta_r(i)
            else
                error stop "Heterogeneous mesh"
            end if
        end function
        
        function get_dim_rad(this) result(dim)
            implicit none
            class(spatial_discr_rad_c) :: this
            integer(kind=4) :: dim
            dim=this%dim
        end function
        
        subroutine compute_r_max(this)
            implicit none
            class(spatial_discr_rad_c) :: this
            this%r_max=this%r_min+sum(this%Delta_r)
        end subroutine
        
        subroutine compute_Delta_r(this)
            implicit none
            class(spatial_discr_rad_c) :: this
            if (.not. allocated(this%Delta_r)) then
                allocate(this%Delta_r(this%Num_targets-this%targets_flag))
            end if
            this%Delta_r=(this%r_max-this%r_min)/(this%Num_targets-this%targets_flag)
        end subroutine
        
        subroutine compute_measure_rad(this)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            if (this%dim == 1) then
                this%measure=this%r_max-this%r_min !> Length in 1D
            else if (this%dim == 2) then
                this%measure=pi*(this%r_max**2-this%r_min**2) !> Area in 2D
            else if (this%dim == 3) then
                this%measure=(4d0/3d0)*pi*(this%r_max**3-this%r_min**3) !> Volume in 3D
            else
                error stop "Dimension not implemented yet"
            end if
        end subroutine
        
        subroutine compute_dimless_mesh_rad(this,char_measure)
        implicit none
        class(spatial_discr_rad_c) :: this
        real(kind=8), intent(in) :: char_measure !> Characteristic measure
        real(kind=8) :: r_max_dimless,r_min_dimless
        this%r_max=this%r_max/char_measure
        this%r_min=this%r_min/char_measure
        this%Delta_r=this%Delta_r/char_measure
        this%measure=this%measure/char_measure
        end subroutine
end module
