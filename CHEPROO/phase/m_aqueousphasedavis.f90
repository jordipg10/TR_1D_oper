module m_aqueousphasedavis
!-------------------------------------------------------------------------
!
!   $Description: Aqueous phase object for dilute solutions.
! Compute activity coefficients according Debye-H�ckel aproximation.
! The water activity is computed according Garrels and Christ (1925).  
!
!   $Use: 
! m_species,
! m_parentaqueousphase.
! m_general_tools_cheproo.
! m_constants. 
!
!   $Author: Sergio Andr�s Bea Jofr� 
!
!   $License: UPC-CSIC
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_species
use m_parentaqueousphase
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%------------------------------------------------------------------
!%------------------------------------------------------------------
private                        ::
!%------------------------------------------------------------------
!%------------------------------------------------------------------
public                         :: &
create_ &                  ! Create the object
,destroy_ &                ! Destroy the object 
,set_ &                    ! Set attributes in the object  
,set_pparent_ &            ! Set pointer to parent aqueous phase object
,compute_act_coeff_ &      ! Compute activity coefficients 
,compute_dact_coeff_ &     ! Compute derivate of the activity coefficients 
,update_ &                 ! Update parameters that are depending of the temperature 
,write_ &                  ! Write in ascii the attributes encapsulated in the object. 
,assignment(=)             ! Copy an object in other object 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
real*8, private, parameter     :: &
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
bht5 = 1.133333333D-11, &
DLT1 = 1.740000000D-02, &
DLT2 = 1.509047619D-03, &
DLT3 = -2.605904762D-05, &
DLT4 = 1.382857143D-07, &
DLT5 = 0.000000000D+00, &
DHT1 = 1.090000000D-01, &
DHT2 = -1.483333333D-03, &
DHT3 = 1.173333333D-05, &
DHT4 = -3.466666667D-08, &
DHT5 = 2.666666667D-11
!%------------------------------------------------------------------------
!% Type definition 
!%------------------------------------------------------------------------
type, public::t_aqueousphasedavis
 
private                                ::
 
type (t_parentaqueousphase), pointer   :: pp         ! Parent aqueous phase 
 
real*8                                 :: agamma     ! Debye-H�ckel parameter

real*8                                 :: bgamma     ! Debye-H�ckel parameter

real*8                                 :: bdot       ! Debye-H�ckel parameter
 
end type t_aqueousphasedavis
 
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface create_
 
module procedure create_aqphdavis
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface destroy_
 
module procedure destroy_aqphdavis
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_
 
module procedure set_aqphdavis
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface set_pparent_
 
module procedure set_pparent_aqphdavis
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_act_coeff_
 
module procedure compute_act_coeff_aqphdavis
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface compute_dact_coeff_
 
module procedure compute_dact_coeff_aqphdavis
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface update_
 
module procedure update_temp_param_aqphdavis
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface write_
 
module procedure write_aqphdavis
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
interface assignment (=)
 
module procedure copy_aqphdavis
 
end interface
!%------------------------------------------------------------------------
!%------------------------------------------------------------------------
Contains
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_aqphdavis &
(this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create the object
!
!   $Arguments:
!
 
type(t_aqueousphasedavis), intent(inout):: this  ! Type aqueous phase Debye-H�ckel
 
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
!%------------------------------------------------------------
this%agamma = 0.0d0
this%bgamma = 0.0d0
this%bdot = 0.0d0
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_aqphdavis &
(this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy the object 
!
!   $Arguments:
!
 
type(t_aqueousphasedavis), intent(inout)     :: this  ! Type aqueous phase Debye-H�ckel
 
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
!%------------------------------------------------------------
this%agamma = 0.0d0
this%bgamma = 0.0d0
this%bdot = 0.0d0
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_aqphdavis &
(this, &
pp, &
iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set attributes in the object
!
!   $Arguments:
!
 
type(t_aqueousphasedavis), intent(inout)         :: this    ! Type aqueous phase Debye-H�ckel

type(t_parentaqueousphase), intent(in), target   :: pp      ! Type parent aqueous phase variable 

logical, intent(out)                             :: iserror ! iserror=true, then there was an error
 
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
!% Update the parent phase according to reference temperature 
!%------------------------------------------------------------
call update_ (this,this%pp%pp%tempref,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_pparent_aqphdavis &
(this, &
pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set pointer to parent aqueous phase object
!
!   $Arguments:
!

type(t_aqueousphasedavis), intent(inout)            :: this  ! Type aqueous phase (Debye-H�ckel model) variable

type(t_parentaqueousphase), intent(in), target      :: pp    ! Type parent aqueous phase variable 
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
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_temp_param_aqphdavis &
(this, &
temp, &
iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update parameters that are depending of the temperature 
! in the object 
!
!   $Arguments:
!
 
 
type(t_aqueousphasedavis), intent(inout) :: this     ! Type aqueous phase (Debye-H�ckel model) variable

real*8, intent(in)                       :: temp     ! Temperature [C] 

logical, intent(out)                     :: iserror  ! iserror=true, then there was an error
 
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

!%-------------------------------------------------------------
iserror=.false.
!%-------------------------------------------------------------
if (temp<=100.d0) then
 
this%agamma= alt1+(alt2+(alt3+(alt4+alt5*temp)*temp)*temp)*temp
 
this%bgamma= blt1+(blt2+(blt3+(blt4+blt5*temp)*temp)*temp)*temp
 
this%bdot= dlt1+(dlt2+(dlt3+(dlt4+dlt5*temp)*temp)*temp)*temp
 
else
 
this%agamma= aht1+(aht2+(aht3+(aht4+aht5*temp)*temp)*temp)*temp
 
this%bgamma= bht1+(bht2+(bht3+(bht4+bht5*temp)*temp)*temp)*temp
 
this%bdot= dht1+(dht2+(dht3+(dht4+dht5*temp)*temp)*temp)*temp
 
end if
 
!%-------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_aqphdavis &
(this, &
ioutput, &
iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Write in ascii the attributes encapsulated in the object. 
!
!   $Arguments:
!
 
type (t_aqueousphasedavis), intent(in) :: this     ! Type aqueous phase (Debye-H�ckel model) variable

integer, intent(in)                    :: ioutput  ! Output unit

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
write (ioutput,*) '-    Activity model: Debye-H�ckel  -'
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
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_act_coeff_aqphdavis &
(this, &
g, &
c, &
iserror, &
ionstr)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute activity coefficients with Davis model
!
!   $Arguments:
!
 
type(t_aqueousphasedavis), intent(in)           :: this    ! Type aqueous phase (Debye-H�ckel model) variable

real*8, pointer, dimension(:)                   :: g       ! Activity coefficients vector 

real*8, intent(in), dimension(this%pp%pp%numsp) :: c       ! Concentrations vector 

real*8, intent(out), optional                   :: ionstr  ! Ionic strength 

logical, intent(out)                            :: iserror ! iserror=true, then there was an error
 
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
ionstrloc, &
gam, &
a0, &
charge, &
tots
integer                       :: &
i, &
j
logical                       :: &
haveionstr
character(len=100)            :: &
msg 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%----------------------------------------------------------
call check_pointer_ (g,this%pp%pp%numsp,.true.)
g=1.0
!%----------------------------------------------------------
!% Check optional arguments 
!%----------------------------------------------------------
haveionstr = present (ionstr)
!%----------------------------------------------------------
!% Compute ionic strength 
!%----------------------------------------------------------
call compute_ionstr_ (this%pp,ionstrloc,c,iserror)
if (iserror) goto 10 
!%----------------------------------------------------------
do i=1,this%pp%pp%numsp
 gam=0.0d0
 if(i/=this%pp%ithe.and.i/=this%pp%ithw) then
  call get_prop_ (this%pp%pp%pspecies(i)%ptr,a0,'ionsize',msg,iserror)
  if (iserror) goto 10
  call get_prop_ (this%pp%pp%pspecies(i)%ptr,charge,'charge',msg,iserror)
  if (iserror) goto 10
  gam=1.0d0+this%bgamma*a0*dsqrt(dabs(ionstrloc))
  gam=-this%agamma*charge*charge*dsqrt(dabs(ionstrloc))/gam
  gam=gam+this%bdot*ionstrloc
  g(i)=10.0d0**gam
 end if
end do
!%---------------------------------------------------------- 
!% Water activity is computed from expression reported by 
!% Garrels and Christ (1965)
!%----------------------------------------------------------
if (this%pp%ithw>0) then
 call compute_sum_c_ (this%pp,tots,c)
 g(this%pp%ithw) = 1.0d0 - kgwmol * tots
end if
!%----------------------------------------------------------
!% Activity coefficient for electron
!%----------------------------------------------------------
if (this%pp%ithe>0) g(this%pp%ithe) = 1.0d0
!%-----------------------------------------------------------
if (haveionstr) ionstr=ionstrloc
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
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dact_coeff_aqphdavis &
(this, & 
dg,  & 
c,   & 
dc,  & 
iserror, & 
dionstr, & 
dtemp, &
g)                 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the derivates of the activity coefficients 
! of aqueous species according Extended Debye-H�ckel model.
! Derivatives of the water activity is computed deriving the expression 
! reported by Garrels and Christ (1965).
!
!   $Arguments:
!
 
type(t_aqueousphasedavis), intent(in)                     :: this      ! Type aqueous phase (Debye-H�ckel model) variable

real*8, pointer, dimension(:,:)                           :: dg        ! Derivatives of the activity coefficients 

real*8, intent(in), dimension(this%pp%pp%numsp)           :: c         ! Concentration vector 

real*8, intent(in), optional, dimension(this%pp%pp%numsp) :: g         ! Activity coefficients vector 

real*8, intent(in), dimension(:,:)                        :: dc        ! Derivatives of the concentrations 

real*8, pointer, optional, dimension(:)                   :: dionstr   ! Derivatives of the ionic strength 

real*8, intent(in), optional, dimension(:)                :: dtemp     ! Derivatives of the temperature. 

logical, intent(out)                                      :: iserror   ! iserror=true, then there was an error
 
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
den
integer                       :: &
i, &
j, &
ndimder
real*8, pointer               :: &
dstr2(:) => null (), &
gloc(:) => null (), &
sumdc(:) => null ()
logical                       :: &
havedionstr, &
haveg, &
havedtemp 
character(len=100)   :: &
msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%---------------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------------
!% Check optional arguments 
!%---------------------------------------------------------------
havedionstr=present(dionstr)
haveg=present(g)
havedtemp=present(dtemp) 
!%---------------------------------------------------------------
ndimder=size(dc,2)
call check_pointer_ (dg,this%pp%pp%numsp,ndimder,.true.)
!%---------------------------------------------------------------
!% Compute gloc is not present as argument
!%---------------------------------------------------------------
if (haveg) then
 call check_pointer_ (gloc,this%pp%pp%numsp,.true.)
 gloc=g
else
 call compute_act_coeff_(this,gloc,c,iserror)
 if (iserror) goto 20   
end if
!%---------------------------------------------------------------
!% Compute ionic strength
!%---------------------------------------------------------------
call compute_ionstr_ (this%pp,str2,c,iserror)
if (iserror) goto 20
!%---------------------------------------------------------------
!% Compute derivatives of the ionic strength
!%---------------------------------------------------------------
call compute_dionstr_ (this%pp,dstr2,dc,iserror)
if (iserror) goto 20 
!%---------------------------------------------------------------
relstr=dsqrt(dabs(str2))
!%---------------------------------------------------------------
do i=1,this%pp%pp%numsp
!%---------
 call get_prop_ (this%pp%pp%pspecies(i)%ptr,z,'charge',msg,iserror)
 if (iserror) goto 20
!%---------
 call get_prop_ (this%pp%pp%pspecies(i)%ptr,a0,'ionsize',msg,iserror)
 if (iserror) goto 20
!%---------
if(i/=this%pp%ithw.and.i/=this%pp%ithe) then
 den= 1.0d0+this%bgamma*a0*relstr
 den=den*den*2.0d0*relstr
 dum=this%agamma*z*z/den-this%bdot
 dum=-2.30258509299405d0*gloc(i)*dum
 dg(i,:)=dum*dstr2
else if (i==this%pp%ithw) then
 call compute_sum_dc_ (this%pp,sumdc,dc)
 dg(i,:) = - kgwmol * sumdc
 call check_pointer_ (sumdc,1,.false.)
end if
!%---------
end do
!%-------------------------------------------------------------
if (havedionstr) then
 call check_pointer_ (dionstr,ndimder,.true.)
 dionstr=dstr2
end if
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (dstr2,1,.false.)
call check_pointer_ (gloc,1,.false.)
if (iserror) goto 10
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
!%***************Public subroutines***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_aqphdavis &
(copied, &
this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy an object in other object  
!
!   $Arguments:
!
 
type(t_aqueousphasedavis), intent(in) :: this    ! Type aqueous phase (Debye-H�ckel model) variable

type(t_aqueousphasedavis), intent(out):: copied  ! Type aqueous phase (Debye-H�ckel model) variable
 
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
copied%agamma = this%agamma
copied%bgamma = this%bgamma
copied%bdot   = this%bdot
!%-----------------------------------------------------------
return
end subroutine
!%***********************************************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
end module 