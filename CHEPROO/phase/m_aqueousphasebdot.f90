module m_aqueousphasebdot
!-------------------------------------------------------------------------
!
!   $Description: Represent the aqueous phase with Truesdell Jones model for compute the activity coefficients
!
!   $Use: use m_parentaqueousphase
! use m_species
! use m_general_tools_cheproo
! use m_constants
!
!   $Author:
!
!   $License:
!
!-------------------------------------------------------------------------
use m_parentaqueousphase
use m_species
use m_general_tools_cheproo
use m_constants_cheproo
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
private   ::
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
public:: &
create_ &
,destroy_ &
,set_ &
,set_pparent_ &
,compute_act_coeff_ &
,compute_dact_coeff_ &
,update_ &
,assignment(=) &
,write_
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
real*8, private, parameter      :: &
alt1 = 4.910300000D-01, &
alt2 = 5.808571429D-04, &
alt3 = 5.527142857D-06, &
alt4 = -4.857142857D-09, &
alt5 = 0.000000000D+00, &
aht1 = 6.440000000D-01, &
aht2 = -3.436166667D-03, &
aht3 = 4.408833333D-05, &
aht4 = -1.691333333D-07, &
aht5 = 2.766666667D-10, &
blt1 = 3.247000000D-01, &
blt2 = 1.309285714D-04, &
blt3 = 5.502380952D-07, &
blt4 = -1.095238095D-09, &
blt5 = 0.000000000D+00, &
bht1 = 3.302000000D-01, &
bht2 = -1.650000000D-05, &
bht3 = 1.991666667D-06, &
bht4 = -7.400000000D-09, &
bht5 = 1.133333333D-11
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
type, public::t_aqueousphasebdot
 
private                                ::
 
type (t_parentaqueousphase), pointer   :: pp      ! Parent aqueous phase 

real*8                                 :: agamma  

real*8                                 :: bgamma
 
end type t_aqueousphasebdot
 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
 
module procedure create_aqphbdot
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
 
module procedure destroy_aqphbdot
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
 
module procedure set_aqphbdot
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_pparent_
 
module procedure set_pparent_aqphbdot
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_act_coeff_
 
module procedure compute_act_coeff_aqphbdot
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dact_coeff_
 
module procedure compute_dact_coeff_aqphbdot
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_
 
module procedure update_temp_param_aqphbdot
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
 
module procedure write_aqphbdot
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment (=)
 
module procedure copy_aqphbdot
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
Contains
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_aqphbdot &
(this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create aqueous phase object with Truesdell Jones thermodynamic model
!
!   $Arguments:
!
 
type(t_aqueousphasebdot), intent(inout)::this 
 
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
this%pp => null ()
this%agamma = 0.0d0
this%bgamma = 0.0d0
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_aqphbdot &
(this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy aqueous phase object with Truesdell Jones thermodynamic model
!
!   $Arguments:
!
 
type(t_aqueousphasebdot), intent(inout):: this
  
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical                              :: &
iserror 
!-------------------------------------------------------------------------
!
!   $code
!
 


!%------------------------------------------------------------
this%pp => null ()
this%agamma = 0.0d0
this%bgamma = 0.0d0
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_aqphbdot &
(this, &
pp, &
iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_aqueousphasebdot), intent(inout)            :: this

type(t_parentaqueousphase), intent(in), target     :: pp 

logical, intent(out)                               :: iserror 
 
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
 

iserror=.false. 

!%------------------------------------------------------------
this%pp => pp
!%------------------------------------------------------------
!% Update the specialization according reference temperature 
!%------------------------------------------------------------
call update_ (this,this%pp%pp%tempref,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_pparent_aqphbdot &
(this, &
pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_aqueousphasebdot), intent(inout):: this

type(t_parentaqueousphase), target     :: pp
 
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
this%pp => pp
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_aqphbdot &
(this, &
temp, &
iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_aqueousphasebdot), intent(inout)  :: this

real*8, intent(in)                       :: temp

logical, intent(out)                     :: iserror 
 
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
iserror=.false.
!%-----------------------------------------------------------
if (temp<=100.d0) then
 
this%agamma= alt1+(alt2+(alt3+(alt4+alt5*temp)*temp)*temp)*temp
 
this%bgamma= blt1+(blt2+(blt3+(blt4+blt5*temp)*temp)*temp)*temp
 
else
 
this%agamma= aht1+(aht2+(aht3+(aht4+aht5*temp)*temp)*temp)*temp
 
this%bgamma= bht1+(bht2+(bht3+(bht4+bht5*temp)*temp)*temp)*temp
 
 
end if
 
!%-----------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_act_coeff_aqphbdot &
(this, &
g, &
c, &
iserror, &
ionstr)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_aqueousphasebdot), intent(in)      :: this

real*8, pointer                           :: g(:)

real*8, intent(in)                        :: c(this%pp%pp%numsp)

real*8, intent(out), optional             :: ionstr 

logical, intent(out)                      :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                        :: &
str2, &
gam, &
a0, &
z, &
tots, &
bdot
integer                       :: &
i, &
j
logical                       :: &
haveionstr
character(len=100)   :: &
msg
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
call check_pointer_ (g,this%pp%pp%numsp,.true.)
g=1.0
!%----------------------------------------------------------
haveionstr = present (ionstr)
!%----------------------------------------------------------
call compute_ionstr_ &
(this%pp, &
str2, &
c, &
iserror)
if (iserror) goto 10 
!%----------------------------------------------------------
 
do i=1,this%pp%pp%numsp
gam=0.0D0
if(i/=this%pp%ithe.and.i/=this%pp%ithw) then
 
call get_prop_ (this%pp%pp%pspecies(i)%ptr,a0,'ionsize', &
msg,iserror)
if (iserror) goto 10
call get_prop_ (this%pp%pp%pspecies(i)%ptr,z,'charge', &
msg,iserror)
if (iserror) goto 10
call get_prop_ (this%pp%pp%pspecies(i)%ptr,bdot,'bdot', &
msg,iserror)
if (iserror) goto 10
gam=1.0D0+this%bgamma*a0*DSQRT(DABS(str2))
gam=-this%agamma*z*z*DSQRT(DABS(str2))/gam
gam=gam+bdot*str2
g(i)=10.0d0**gam
end if
end do
 
!%------Activity coefficient of water (Garrels and Christ, 1965)
if (this%pp%ithw.gt.0) then
 
call compute_sum_c_ &
(this%pp, &
tots, &
c)
 
g(this%pp%ithw) = 1.0d0 - kgwmol * tots
 
end if
!%----------------------------Activity coefficient for electron
if (this%pp%ithe>0) g(this%pp%ithe) = 1.0D0
!%------------------------------------------------------------
if (haveionstr) ionstr=str2
!%------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%pp%name
print *,'Service: compute_act_coeff_'
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
subroutine compute_dact_coeff_aqphbdot &
(this,          &
 dg,            &
 c,             &
 dc, &
 iserror, &
 dionstr)

implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the derivates of the activity coefficients of
!%   aqueous species according Extended Debye-H�ckel model.
!%   The derivate of the activity of the water species is
!%   computed according the derivate of Garrels and Christ
!%   (1965)
!
!   $Arguments:
!
 
type(t_aqueousphasebdot), intent(in)        :: this

real*8, pointer                             :: dg(:,:)

real*8, intent(in)                          :: c(this%pp%pp%numsp)

real*8, intent(in)                          :: dc(:,:)

real*8, pointer, optional                   :: dionstr(:) 

logical, intent(out)                        :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                        :: &
str2, &
relstr, &
gam, &
a0, &
z, &
dum, &
den, &
bdot
 
integer                       :: &
i, &
j, &
ndim
 
real*8, pointer               :: &
dstr2(:) => null (), &
g(:) => null (), &
sum(:) => null ()
logical                       :: &
havedionstr
character(len=100)            :: &
msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------------
msg=''
iserror=.false.
!%-----------------------------------------------------------------
!% Check option al arguments 
!%-----------------------------------------------------------------
havedionstr=present(dionstr)
!%----------------------------------------------------------------- 
ndim=size(dc,2)
call check_pointer_ (dg,this%pp%pp%numsp,ndim,.true.)
!%----------------------------------------------------------------- 
!% Compute gammasaqt is not present as argument
!%----------------------------------------------------------------- 
call compute_act_coeff_ (this,g,c,iserror)
if (iserror) goto 10 
!%----------------------------------------------------------------- 
!%-----------------------------------------compute ionic strength
!%----------------------------------------------------------------- 
call compute_ionstr_ (this%pp,str2,c,iserror)
if (iserror) goto 10 
!%----------------------------------------------------------------- 
!% Compute derivatives of the ionic strength
!%----------------------------------------------------------------- 
call compute_dionstr_(this%pp,dstr2,dc,iserror)
if (iserror) goto 10 
!%---------------------------------------------------------------
relstr=dsqrt(dabs(str2))
!%-----------------------------------------------------------------
do i=1,this%pp%pp%numsp
 
call get_prop_ (this%pp%pp%pspecies(i)%ptr,z,'charge', msg,iserror)
if (iserror) goto 10
call get_prop_ (this%pp%pp%pspecies(i)%ptr,a0,'ionsize', msg,iserror)
if (iserror) goto 10
call get_prop_ (this%pp%pp%pspecies(i)%ptr,bdot,'bdot', msg,iserror)
if (iserror) goto 10

if(i/=this%pp%ithw.and.i/=this%pp%ithe) then
 
den= 1.0d0+this%bgamma*a0*relstr
den=den*den*2.0d0*relstr
dum=this%agamma*z*z/den-bdot
dum=-2.30258509299405d0*g(i)*dum
dg(i,:)=dum*dstr2
 
else if (i==this%pp%ithw) then
 
 call compute_sum_dc_ (this%pp,sum,dc)
 
 dg(i,:) = - kgwmol * sum
 
end if
 
end do
!%-------------------------------------------------------------
!% Allocate pointer 
!%-------------------------------------------------------------
if (havedionstr) then
 call check_pointer_ (dionstr,ndim,.true.)
 dionstr=dstr2
end if
!%-------------------------------------------------------------
!% Deallocate local pointer 
!%-------------------------------------------------------------
call check_pointer_ (dstr2,1,.false.)
call check_pointer_ (g,1,.false.)
call check_pointer_ (sum,1,.false.)
!%-------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'Phase:'
print *,'Name:', this%pp%pp%name
print *,'Service: compute_dact_coeff_'
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
subroutine copy_aqphbdot &
(targetobj, &
 sourceobj)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy a object in other species object.  
!
!   $Arguments:
!
 
type(t_aqueousphasebdot), intent(in) :: sourceobj

type(t_aqueousphasebdot), intent(out):: targetobj
 
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
!-------------------------------------------------------------------------
!
!   $code
!
 


!%-----------------------------------------------------------
targetobj%agamma = sourceobj%agamma
targetobj%bgamma = sourceobj%bgamma
!%-----------------------------------------------------------
targetobj%pp => null ()
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_aqphbdot &
(this, &
ioutput, &
iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_aqueousphasebdot), intent(in)   :: this

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
character(len=100)      :: &
msg 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
write (ioutput,*) '------------------------------------'
write (ioutput,*) '-         Aqueous Phase            -'
write (ioutput,*) '-    Activity model: BDOT          -'
write (ioutput,*) '-            extended              -'
write (ioutput,*) '-Water activity computed according -'
write (ioutput,*) '-    Garrels and Christ (1965)     -'
!%-----------------------------------------------------------
!% Write the parent aqueous phase 
!%-----------------------------------------------------------
call write_ (this%pp,ioutput,iserror)
!%-----------------------------------------------------------
write (ioutput,*) '------------------------------------'
!%-----------------------------------------------------------
return
 
10 continue 
print *,'************************'
print *,'Phase:'
print *,'Name:',this%pp%pp%name
print *,'Service: write_'
print *, msg
print *,'************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module 