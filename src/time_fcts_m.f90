module time_fct_m
implicit none
save
type, public :: time_fct_c !> Time function class
    integer(kind=4) :: ntime = 0 !> Number of time steps
    !integer(kind=4), allocatable :: time_ind(:) !> Time step indices
    real(kind=8), allocatable :: time_step(:) !> Time step sizes
contains
procedure :: set_ntime !> Set number of time steps
procedure :: allocate_time_step !> allocate time step indices
procedure :: set_time_step !> Set time step sizes
!procedure :: get_time_ind => time_fcts_get_time_ind !> Get time step indices
procedure :: get_time_step !> Gets time step size
end type time_fct_c

type, public, extends(time_fct_c) :: time_fct_int_c !> integer time function subclass
    integer(kind=4), allocatable :: time_series(:) !> Time series
contains
procedure :: set_time_series=>set_time_series_int !> Set time series
end type time_fct_int_c

type, public, extends(time_fct_c) :: time_fct_real_c !> real time function subclass
    real(kind=8), allocatable :: time_series(:) !> Time series
contains
procedure :: set_time_series=>set_time_series_real !> Set time series
procedure :: read_time_series=>read_time_series_real !> read time series
procedure :: allocate_time_series=>allocate_time_series_real !> allocate time series
end type time_fct_real_c

    contains
    
    subroutine set_ntime(this,ntime)
    class(time_fct_c) :: this
    integer(kind=4), intent(in) :: ntime
    if (ntime < 0) then
        error stop "Error: ntime must be non-negative."
    end if
    this%ntime = ntime
    end subroutine 
    
    subroutine set_time_step(this,time_step)
    class(time_fct_c) :: this
    real(kind=8), intent(in) :: time_step(:)
    if (minval(time_step) < 0d0) then
        error stop "Error: attribute time_step must be non-negative."
    end if
    this%ntime = size(time_step)
    this%time_step = time_step
    end subroutine set_time_step
    
    function get_time_step(this,ind) result(time_step)
    class(time_fct_c), intent(in) :: this
    integer(kind=4), intent(in), optional :: ind !> Index of time step
    real(kind=8) :: time_step
    
    if (present(ind)) then
        if (ind < 1 .or. ind > this%ntime) then
            error stop "Error: Index out of bounds."
        end if
        time_step = this%time_step(ind)
    else
        time_step = this%time_step(1) ! Default to first time step
    end if
    end function get_time_step
    
    subroutine set_time_series_int(this,time_series)
    class(time_fct_int_c) :: this
    integer(kind=4), intent(in) :: time_series(:)
    this%ntime = size(time_series)
    this%time_series = time_series
    end subroutine set_time_series_int
    
    subroutine set_time_series_real(this,time_series)
    class(time_fct_real_c) :: this
    real(kind=8), intent(in) :: time_series(:)
    this%ntime = size(time_series)
    this%time_series = time_series
    end subroutine set_time_series_real
    
    subroutine read_time_series_real(this,filename)
    class(time_fct_real_c) :: this
    character(len=*), intent(in) :: filename
    
    integer(kind=4) :: i !> Loop index
    open(unit=1, file=filename, status='old', action='read') !> Open the file for reading
    read(1,*) this%ntime !> Read number of time steps
    call this%allocate_time_series() !> Allocate time series array
    do i = 1, this%ntime !> Read each time step
        read(1,*) this%time_series(i)
    end do
    close(1) !> Close the file
    end subroutine read_time_series_real
    
    subroutine allocate_time_step(this,ntime)
    class(time_fct_c) :: this
    integer(kind=4), intent(in), optional :: ntime
    if (present(ntime)) then
        if (ntime < 0) then
            error stop "Error: ntime must be non-negative."
        end if
        this%ntime = ntime
    end if
    allocate(this%time_step(this%ntime))
    end subroutine allocate_time_step
    
    subroutine allocate_time_series_real(this,ntime)
    class(time_fct_real_c) :: this
    integer(kind=4), intent(in), optional :: ntime
    if (present(ntime)) then
        if (ntime < 0) then
            error stop "Error: ntime must be non-negative."
        end if
        this%ntime = ntime
    end if
    allocate(this%time_series(this%ntime))
    end subroutine allocate_time_series_real
    
end module time_fct_m