module m_parentaqueousphase
!-------------------------------------------------------------------------
!
!   $Description: Parent aqueous phase 
!
!   $Use: 
! m_parentphase.
! m_general_tools_cheproo.
! m_constants. 
!
!   $Author: Sergio Andr�s Bea Jofr�
!
!   $License: CSIC-UPC
!
!-------------------------------------------------------------------------
use m_species
use m_parentphase
use m_general_tools_cheproo
use m_constants_cheproo
implicit none
!%--------------------------------------------------------
!%--------------------------------------------------------
private                          ::
!%--------------------------------------------------------
!% Public services 
!%--------------------------------------------------------
public:: &
create_ &                    ! Create parent aqueous phase object 
,destroy_ &                  ! Destroy parent aqueous phase object 
,set_ &                      ! Set attributes in the parent aqueous phase 
,set_pparent_ &              ! Set the parent phase 
,compute_ionstr_ &           ! Compute the ionic strength (eq. 4)
,compute_dionstr_ &          ! Compute the derivate of the ionic strength (eq. 21) 
,compute_sum_c_ &            ! Compute M (eq. 6)
,compute_sum_absz_c_ &       ! Compute Z (eq. 5)
,compute_sum_zcat_ccat_ &    ! Compute Zcat
,compute_sum_dc_ &           ! Compute dM (eq. 23)
,compute_mass_ &             ! Compute the total mass of solutes   
,compute_volid_ &            ! Compute the total volume of the aqueous phase 
,compute_dmass_dc_ &         ! Compute the derivative of mass of solutes
,compute_density_ &          ! Compute density (and derivatives if asked)
,update_ &                   ! Update parameters that depends of the temperature 
,write_ &                    ! Write the parent aqueous phase 
,assignment(=)               ! Copy an object in other object  
!%------------------------------------------------------------------------
!% Constant parameters 
!%------------------------------------------------------------------------
real*8, public, parameter:: &
dielcoeff1=0.256983619d+3, &
dielcoeff2=-0.847552262d+0, &
dielcoeff3=0.794750869d-3, &
dielcoeff4=0.347887496d-6, &
dielcoeff5=0.721027523d-9
!%------------------------------------------------------------------------
!%Phase density model codes
!%------------------------------------------------------------------------
integer,parameter:: &
ip_lin_p_t_s=1 ,&   !Lineal expression: pressure, temperature, salinity
ip_exp_p_t_s=2      !Exponential expression: pressure, temperature, salinity

!%------------------------------------------------------------------------
!% Type definition 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
type, public::t_parentaqueousphase
  
type (t_parentphase), pointer    :: pp    ! Parent phase pointer 
 
integer                          :: ithe  ! Index of electron species

integer                          :: ithw  ! Index of water species

integer                          :: ithh  ! Index of hydrogen species 

integer                          :: densCode ! Desity model code

real*8,dimension(:),pointer             :: paramValues  ! List of parameters values
character(len=ip_namelenght),dimension(:),pointer  :: paramNames   ! List of parameters names

 
end type t_parentaqueousphase
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
!% Interfaces 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
 
module procedure create_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
 
module procedure destroy_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
 
module procedure set_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_pparent_
 
module procedure set_pparent_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_ionstr_
 
module procedure compute_ionstr_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dionstr_
 
module procedure compute_dionstr_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_density_
 
module procedure compute_density_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_sum_c_
 
module procedure compute_sum_c_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_sum_absz_c_
 
module procedure compute_sum_absz_c_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_sum_zcat_ccat_
 
module procedure compute_sum_zcat_ccat_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_sum_dc_
 
module procedure compute_sum_dc_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_mass_
 
module procedure compute_mass_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_volid_
 
module procedure compute_volid_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dmass_dc_
 
module procedure compute_dmass_dc_paqph
 
end interface

!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_
 
module procedure update_temp_param_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
 
module procedure write_paqph
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment (=)
 
module procedure copy_paqph
 
end interface
!%------------------------------------------------------------
!%------------------------------------------------------------
contains
!%------------------------------------------------------------
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_paqph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create parent aqueous phase object
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(inout):: this   ! Parent aqueous phase type     
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
this%ithw=0
this%ithe=0
this%ithh=0
this%densCode=0

nullify(this%paramValues)
nullify(this%paramNames)
!%-----------------------------------------------------------
!% Allocate and create the parent phase 
!%-----------------------------------------------------------
allocate (this%pp)
call create_ (this%pp)

!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_paqph &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy parent aqueous phase object
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(inout)     :: this  ! Parent aqueous phase type     
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
this%pp => null ()
call check_pointer_ (this%paramNames,size(this%paramNames),.false.)
call check_pointer_ (this%paramValues,size(this%paramValues),.false.)

!%-----------------------------------------------------------
this%ithw=0
this%ithe=0
this%ithh=0
this%densCode=0
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_paqph &
   (this, &
    pp, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the parent phase 
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(inout):: this      ! Parent aqueous phase type     

type(t_parentphase), intent(in), target  :: pp        ! Parent phase type

logical, intent(out)                     :: iserror   ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer        :: &
 i 
character(len=100) &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
iserror=.false. 
msg=''
!%-----------------------------------------------------------
this%pp => pp
!%-----------------------------------------------------------
this%pp%numprop=1
call check_pointer_ (this%pp%prop,this%pp%numprop,.true.)
call check_pointer_ (this%pp%nameprop,this%pp%numprop,.true.)
call check_pointer_ (this%pp%modelprop,this%pp%numprop,.true.)
this%pp%nameprop(1)='dielectriccte'
this%pp%modelprop(1)='Helgeson and Kirkham(1974)'
this%pp%prop(1)=0.0d0
!%-----------------------------------------------------------
this%ithw=0
this%ithe=0
this%ithh=0
!%-----------------------------------------------------------
!% Check if some aqueous species is water, proton or electron, 
!% save the index
!%-----------------------------------------------------------
do i=1, this%pp%numsp
 
 select case (this%pp%pspecies(i)%ptr%name)
 
 case ('h2o', 'H2O')
 
  this%ithw = i
 
 case ('e-', 'E-')
 
  this%ithe = i
 
 case ('h+', 'H+')
 
  this%ithh = i
 
 end select
 
end do
!%------------------------------------------------------------
!% Update parameters that depend of the temperature 
!%------------------------------------------------------------
call update_ (this,this%pp%tempref,iserror) 
if (iserror) goto 10
!%------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: set_'
print *, msg
print *,'*******************************'
iserror=.true.
return 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_pparent_paqph &
   (this, &
    pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set pasrent phase object 
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(inout) :: this ! Parent aqueous phase type     

type(t_parentphase), intent(in), target   :: pp   ! Parent phase type     
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%-----------------------------------------------------------
this%pp => pp 
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_ionstr_paqph &
   (this, &
    ionstr, &
    c, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the ionic strength (I) eq. (4) 
!
! It is defined as 
!
! I = 0.5 * sum(ci*(zi)^2)
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in)         :: this      ! Parent Aqueous phase type

real*8, intent(out)                            :: ionstr    ! Ionic strength 

real*8, intent(in), dimension(this%pp%numsp)   :: c         ! Molality vector

logical, intent(out)                           :: iserror   ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8               :: &
 charge, &
 sum 
integer              :: &
 i
character(len=100)   :: &
 msg 

!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
sum=0.0d0
!%---------------------------------------------------------- 
do i=1,this%pp%numsp
 
    call get_prop_ (this%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
 
    if (iserror) goto 10
 
    if(i/=this%ithe.and.i/=this%ithw) sum = sum + charge * charge * c(i)
 
end do
!%----------------------------------------------------------
ionstr = 0.5d0 * sum
!%----------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_ionstr_'
print *, msg
print *,'*******************************'
iserror=.true.
return 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dionstr_paqph &
   (this, &
    dionstr, &
    dc, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the derivative of the ionic strength (eq. 21). 
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in) :: this      ! Parent Aqueous Phase type 

real*8, pointer, dimension(:)          :: dionstr   ! Derivative of the ionic strength 

real*8, intent(in), dimension(:,:)     :: dc        ! Derivatives of the molalities 

logical, intent(out)                   :: iserror   ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8               :: &
 charge
integer              :: &
 i, &
 ndim
character(len=100)   :: &
 msg
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
ndim=size(dc,2)
!%------------------------------------------------------------
!% Allocate and zeroing pointer 
!%------------------------------------------------------------
call check_pointer_ (dionstr,ndim,.true.)
!%------------------------------------------------------------
do i=1,this%pp%numsp
 
    call get_prop_ (this%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
    if (iserror) goto 10
    if(i/=this%ithe.and.i/=this%ithw) then
        dionstr = dionstr + charge * charge * dc(i,:)
    end if
 
end do
!%-----------------------------------------------------------
dionstr  = 0.5d0 * dionstr
!%-----------------------------------------------------------
return
 
10 continue 
print *,'*************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_dionstr_'
print *, msg
print *,'*************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sum_c_paqph &
   (this, &
    m, &
    c)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute M (eq. 6)
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in)          :: this     ! Parent aqueous phase type

real*8, intent(in), dimension(this%pp%numsp)    :: c        ! molality vector 

real*8, intent(out)                             :: m        ! M (eq. 6)
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
m=0.0d0
!%----------------------------------------------------------
do i=1,this%pp%numsp
 
 if(i/=this%ithw.and.i/=this%ithe) then
   m = m + c(i)
 end if
 
end do
!%-----------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_mass_paqph &
   (this, &
    mass, &
    mol, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the total mass of solutes  
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in) :: this     ! Parent aqueous phase type

integer, intent(in)                    :: nsp

real*8, intent(out)                    :: mass     ! Mass of salt [kgr] or [kgr/kgr water]

real*8, intent(in), dimension(nsp)     :: mol      ! Molality vector 

logical, intent(out)                   :: iserror  ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)   :: &
 msg
integer              :: &
 isps
real*8               :: &
 mw
character(len=20), parameter :: &
 nameprop='molweight' 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Check the number of species 
!%-----------------------------------------------------------
if (nsp/=this%pp%numsp) then
 msg='Error in number of species'
 goto 10 
end if
!%-----------------------------------------------------------
!% Zeroing 
!%-----------------------------------------------------------
mass=0.0d0
!%-----------------------------------------------------------
do isps=1,nsp
 if (isps/=this%ithw.and.isps/=this%ithe) then
   call get_prop_ (this%pp%pspecies(isps)%ptr,mw,nameprop,msg,iserror)
   if (iserror) goto 10
   mass=mass+mol(isps)*mw*kgrgr
 end if
end do
!%------------------------------------------------------------
return
10 continue 
print *,'*********************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_mass_'
print *, msg
print *,'*********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_volid_paqph &
   (this, &
    volid, &
    c, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the total volume of the solutes
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in) :: this     ! Parent aqueous phase type

integer, intent(in)                    :: nsp

real*8, intent(out)                    :: volid    ! Mass of salt [kgr] or [kgr/kgr water]

real*8, intent(in), dimension(nsp)     :: c        ! Molality vector 

logical, intent(out)                   :: iserror  ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)   :: &
 msg
integer              :: &
 i 
real*8               :: &
 molvol  
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Check the number of species 
!%-----------------------------------------------------------
if (nsp/=this%pp%numsp) then
 msg='Error in number of species'
 goto 10 
end if
!%-----------------------------------------------------------
!% 
!%-----------------------------------------------------------
if (this%ithw==0) then
  msg='Error, water species must be defined in the phase'  
  iserror=.true.
  goto 10  
else
  call get_prop_ (this%pp%pspecies(this%ithw)%ptr,molvol,'molvol',msg,iserror)
  if (iserror) goto 10 
  volid=1.0d3*molvol
end if
!%------------------------------------------------------------
!% Compute the contribution of solutes in the volume of the 
!% solution 
!%------------------------------------------------------------
do i=1,this%pp%numsp
  
  if (i/=this%ithe.and.i/=this%ithw ) then 
    call get_prop_ (this%pp%pspecies(i)%ptr,molvol,'molvol',msg,iserror)
    if (iserror) goto 10 
    volid=volid+c(i)*molvol 
  end if
  
end do
!%------------------------------------------------------------
return
10 continue 
print *,'*********************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_volid_'
print *, msg
print *,'*********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dmass_dc_paqph &
   (this, &
    dm, &
    dc)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Compute derivative of Mass of salt
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in)         :: this        ! Parent aqueous phase type

real*8, intent(in), dimension(:,:)             :: dc          ! Derivatives of the molalities.     

real*8, pointer, dimension(:)                  :: dm          ! Derivative of Mass

 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
 integer              :: &
 ndim, &
 j
character(len=100)   :: &
 msg
integer              :: &
 isps
real*8               :: &
 mw
logical             :: &
iserror 
character(len=20), parameter :: &
 nameprop='molweight' 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
ndim=size(dc,2)
!%-----------------------------------------------------------
!% Allocate and zeroing pointer
!%-----------------------------------------------------------
call check_pointer_ (dm,ndim,.true.)
!%-----------------------------------------------------------

do isps=1,this%pp%numsp
 if (isps/=this%ithw.and.isps/=this%ithe) then
   call get_prop_ (this%pp%pspecies(isps)%ptr,mw,nameprop,msg,iserror)
   if (iserror) goto 10
   mw=mw*kgrgr
   do j=1,ndim
     dm(j)=dm(j)+mw*dc(isps,j)
   enddo
 end if
end do
!%------------------------------------------------------------
return
!%------------------------------------------------------------
return
10 continue 
print *,'*********************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_dmass_dc_'
print *, msg
print *,'*********************'
iserror=.true.
return


end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sum_absz_c_paqph &
   (this, &
    z, &
    c, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute Z (eq. 5)
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in)         :: this     ! Parent aqueous phase type

real*8, intent(in), dimension(this%pp%numsp)   :: c        ! Molality vector

real*8, intent(out)                            :: z        ! Z (eq. 5)  

logical, intent(out)                           :: iserror  ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer              :: &
 i
real*8               :: &
 charge
character(len=100)   :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
z=0.0d0
!%-----------------------------------------------------------
do i=1,this%pp%numsp
 
 call get_prop_ (this%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
 if (iserror) goto 10
 if(i/=this%ithw.and.i/=this%ithe) then
   z = z + dabs(charge) * c(i)
 end if
 
end do
 
!%----------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_sum_absz_c_'
print *, msg
print *,'*******************************'
iserror=.true.
return 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sum_zcat_ccat_paqph &
   (this, &
    zcat, &
    c, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute Z (eq. 5)
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in)         :: this     ! Parent aqueous phase type

real*8, intent(in), dimension(this%pp%numsp)   :: c        ! Molality vector

real*8, intent(out)                            :: zcat     !   

logical, intent(out)                           :: iserror  ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer              :: &
 i
real*8               :: &
 charge
character(len=100)   :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
zcat=0.0d0
!%-----------------------------------------------------------
do i=1,this%pp%numsp
 
 call get_prop_ (this%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
 if (iserror) goto 10
 if(i/=this%ithw.and.i/=this%ithe.and.charge>0.0d0) then
   zcat = zcat + charge * c(i)
 end if
 
end do
 
!%----------------------------------------------------------
return
 
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%name
print *,'Service: compute_sum_zcat_ccat_'
print *, msg
print *,'*******************************'
iserror=.true.
return 
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_sum_dc_paqph &
   (this, &
    dm, &
    dc)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Compute derivative of M (eq. 23)
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in)         :: this        ! Parent aqueous phase type

real*8, intent(in), dimension(:,:)             :: dc          ! Derivatives of the molalities.     

real*8, pointer, dimension(:)                  :: dm          ! Derivative of Z (eq. 22).
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
 integer              :: &
 ndim, &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%-----------------------------------------------------------
ndim=size(dc,2)
!%-----------------------------------------------------------
!% Allocate and zeroing pointer
!%-----------------------------------------------------------
call check_pointer_ (dm,ndim,.true.)
!%-----------------------------------------------------------
do i=1,this%pp%numsp
 
 if(i/=this%ithw.and.i/=this%ithe) then
   dm = dm + dc(i,:)
 end if
 
end do
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_paqph &
   (this, &
    temp, &
	iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update in the parent aqueous phase the parameters that 
! depends of the temperature. 
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(inout)    :: this      ! Parent aqueous phase type

real*8, intent(in)                           :: temp      ! Temperature [C]

logical, intent(out)                         :: iserror   ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                       :: &
 tk
integer                      :: &
 i
character(len=20)            :: &
 nameprop='dielectriccte' 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
iserror=.false. 
!%------------------------------------------------------------
!% Change to Kelvin 
!%------------------------------------------------------------
tk=temp+273.15d0
!%------------------------------------------------------------
!% Update the dielectric constant (Helgeson and Kirkham(1974)
!%------------------------------------------------------------
do i=1,this%pp%numprop
 if (this%pp%nameprop(i)==nameprop) then
 this%pp%prop(1)= dielcoeff1 + (dielcoeff2 + (dielcoeff3 &
                   +(dielcoeff4-dielcoeff5*tk)*tk)*tk)*tk
 end if
end do
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_paqph &
   (this, &
    ioutput, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write parent aqueous phase 
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in) :: this     ! Parent aqueous phase type

integer, intent(in)                    :: ioutput

logical, intent(out)                   :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
iserror=.false.
!%------------------------------------------------------------
!% Write parent phase 
!%------------------------------------------------------------
call write_ (this%pp,ioutput,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_paqph &
  (targetobj, &
   sourceobj)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy an object in other object.  
!
!   $Arguments:
!
 
type(t_parentaqueousphase), intent(in) :: sourceobj   ! Parent aqueous phase type

type(t_parentaqueousphase), intent(out):: targetobj   ! Parent aqueous phase type 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!   $code
!
!%---------------------------------------------------------------
targetobj%ithe = sourceobj%ithe
targetobj%ithh = sourceobj%ithh
targetobj%ithw = sourceobj%ithw
targetobj%densCode=sourceobj%densCode
targetobj%pp = sourceobj%pp

if (associated(sourceobj%paramNames))then
    call check_pointer_ (targetobj%paramNames,size(sourceobj%paramNames),.true.)
    targetobj%paramNames=sourceobj%paramNames
endif
if (associated(sourceobj%paramValues))then
    call check_pointer_ (targetobj%paramValues,size(sourceobj%paramValues),.true.)
    targetobj%paramValues=sourceobj%paramValues
endif

!%---------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_density_paqph &
    (this        ,&
     pliq        ,&
     temp        ,&
     c           ,&
     density     ,&
     iserror     ,&
     dc          ,&
     ddensitydc)


type(t_parentaqueousphase)                      :: this         !Type parent aquous phase

real*8,intent(in)                               :: pliq         !liquid pressure

real*8,intent(in)                               :: temp         !temperature

real*8,dimension(this%pp%numsp), intent(in)     :: c            !concentration array

real*8,intent(out)                              :: density      !Density     

logical,intent(out)                             :: IsError      !iserror=true, then there was an error

real*8,dimension(:,:),optional                  :: dc           !concentration derivatives (needed to calculate density derivatives)

real*8,pointer,dimension(:),optional            :: ddensitydc   !Density derivative

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                               :: &
 msg =""


real*8,save ::&

refdens=-1   ,&
refpre=-1    ,&
alpha=-1     ,&
beta=-1      ,&
gamma=-1

integer:: i
real*8::salinity
character(len=ip_namelenght) :: findname  
!-------------------------------------------------------------------------
!
!   $code
!


!TEMPORAL:this must be added to xml read !!!!!
this%densCode=ip_lin_p_t_s
call check_pointer_ (this%paramNames,5,.true.)
call check_pointer_ (this%paramValues,5,.true.)

this%paramNames(1)='refDensity'
this%paramNames(2)='refPressure'
this%paramNames(3)='alpha'
this%paramNames(4)='beta'
this%paramNames(5)='gamma'



this%paramValues(1)= 1002.6
this%paramValues(2)= 0.1
this%paramValues(3)= -3.4d-4
this%paramValues(4)= 4.5d-4
this%paramValues(5)= 0.6923


select case (this%densCode)

        case(0) !no model at all
            msg="Density can not be computed because no model has been selected"
            goto 10
        case(ip_lin_p_t_s)


            !we search for the parameters for this law
            if (refdens.eq.-1) then
                findname="refDensity"
                call findInArray_(this%paramNames,findname,i)
                if (i.eq.0) then
                    msg ="'refDensity' not found in phase parameter list"
                    goto 10
                endif
                refdens=this%paramValues(i)
             endif
             if (refpre.eq.-1) then
                findname="refPressure"
                call findInArray_(this%paramNames,findname,i)
                if (i.eq.0) then
                    msg ="'refPressure' not found in phase parameter list"
                    goto 10
                endif
                refpre=this%paramValues(i)
             endif
             if (alpha.eq.-1) then
                findname="alpha"
                call findInArray_(this%paramNames,findname,i)
                if (i.eq.0) then
                    msg ="'alpha' not found in phase parameter list"
                    goto 10
                endif
                alpha=this%paramValues(i)
             endif             
             if (beta.eq.-1) then
                findname="beta"
                call findInArray_(this%paramNames,findname,i)
                if (i.eq.0) then
                    msg ="'beta' not found in phase parameter list"
                    goto 10
                endif
                beta=this%paramValues(i)
             endif    
             if (gamma.eq.-1) then
                findname="gamma"
                call findInArray_(this%paramNames,findname,i)
                if (i.eq.0) then
                    msg ="'gamma' not found in phase parameter list"
                    goto 10
                endif
                gamma=this%paramValues(i)
             endif              
        
            
            !we ask for the salinity
            call compute_mass_(this,salinity, c,size(c),iserror)

            if (iserror) goto 10
            
            !now we build the densit
            
            density=refdens*(1+beta*(pliq-refpre)+alpha*temp+gamma*salinity)
            
            !now calculate derivatives, if we have to
            if  (present(ddensitydc) ) then
                if (.not.present(dc)) then
                    msg="dc are needed for calcutate density derivatives"
                    isError=.true.
                    goto 10
                endif 
           !we ask for the derivative of the salinity     
           call compute_dmass_dc_(this,ddensitydc,dc)
           
           !ddensitydc=refdens*(gamma*dsalinity)
           
           ddensitydc=refdens*gamma*ddensitydc
                
            end if
            
    case default
        msg="Not recognized code for Density model"
        goto 10
        
end select

!%----------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Parent AqueousPhase:'
print *,'Service: compute_density_'
print *, msg
print *,'***********************'
iserror=.true.
return



end subroutine
!%***************************************************************
!%***************************************************************
!%***************************************************************
end module m_parentaqueousphase