module m_chemicalsystemclassic
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Use:
!
!   $Author: Sergio Andrés Bea Jofré 
!
!   $License:
!
!-------------------------------------------------------------------------
use m_parentchemicalsystem
use m_reaction
use m_species
use m_phase
use m_surface
use flib_xpath
use flib_sax
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------
!%-------------------------------------------------------
private      ::
!%-------------------------------------------------------
!%-------------------------------------------------------
public       :: &
create_ &
,destroy_ &
,compute_dcmob_ &
,compute_dcads_ &
,compute_usktrk_ &
,compute_umob_ &
,compute_uads_ &
,compute_dumob_ &
,compute_duads_ &
,compute_dusktrk_ &
,compute_iumob_ &
,make_lin_trf_ &
,set_ &
,set_parent_ &
,write_
!%------------------------------------------------------
!%------------------------------------------------------
type, public::t_chemicalsystemclassic
 
private                                 ::
 
type(t_parentchemicalsystem), pointer   :: pp
 
 
end type t_chemicalsystemclassic
 
!%---------------------------------------Public interfaces
!%--------------------------------------------------------
interface create_
 
module procedure create_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface destroy_
 
module procedure destroy_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_usktrk_
 
module procedure compute_usktrk1_chemsys_classic
module procedure compute_usktrk2_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_umob_
 
module procedure compute_umob_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_uads_
 
module procedure compute_uads_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dcmob_
 
module procedure compute_dcmob1_chemsys_classic
module procedure compute_dcmob2_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dcads_
 
module procedure compute_dcads_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dumob_
 
module procedure compute_dumob1_chemsys_classic
module procedure compute_dumob2_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_duads_
 
module procedure compute_duads_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dusktrk_
 
module procedure compute_dusktrk1_chemsys_classic
module procedure compute_dusktrk2_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_iumob_
 
module procedure compute_iumob_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface make_lin_trf_
 
module procedure make_lin_trf_vector_chemsys_classic
module procedure make_lin_trf_array_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_ 
 
module procedure set_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_parent_ 
 
module procedure set_parent_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface write_
 
module procedure write_file_chemsys_classic
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
contains
 
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_chemsys_classic &
   (this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create the classic chemical system
!
!   $Arguments:
!
 
type (t_chemicalsystemclassic), intent(inout) :: this 
 
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
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_chemsys_classic &
   (this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy the classic chemical system
!
!   $Arguments:
!
 
type (t_chemicalsystemclassic), intent(inout)   :: this 
 
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
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_chemsys_classic &
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
 
type (t_chemicalsystemclassic), intent(inout)    :: this

type (T_parentchemicalsystem), intent(in), target:: pp 

logical, intent(out)                             :: iserror
 
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
!%-------------------------------------------Pointer to Parent
this%pp => pp
!%------------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_parent_chemsys_classic &
   (this, &
    pp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set parent 
!
!   $Arguments:
!
 
type (t_chemicalsystemclassic) , intent(inout)        :: this

type (t_parentchemicalsystem), intent(in), target     :: pp

logical, intent(out)                                  :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)              :: &
msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------Set the parent chemical system
this%pp => pp
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: set_parent_'
print *, msg
 print *,'*******************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%   Give the icomp aqueous components in the mobile phases
!%   (aqueous, gas, and non-aqueous phases
!%
!%************************************************************
subroutine compute_iumob_chemsys_classic &
   (this, &
    iumob, &
    naqcol, &
    ngascol, &
    nnonaqcol, &
    c, &
    numsp, &
    icomp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
 
type (t_chemicalsystemclassic):: &
 this
real*8, pointer        :: &
 iumob(:)
integer                :: &
 numsp, &
 icomp, &
 naqcol, &
 ngascol, &
 nnonaqcol
real*8                 :: &
 c(numsp)
logical                :: &
 iserror
!%------------------------------------------------------------
call compute_iumob_ &
   (this%pp, &
    iumob, &
    naqcol, &
    ngascol, &
    nnonaqcol, &
    c, &
    numsp, &
    icomp, &
    iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_usktrk1_chemsys_classic &
   (this, &
    usktrk, &
    ncomp, &
    c, &
    g, &
    alpha, &
    numsp, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemclassic), intent(in) :: this

real*8, pointer                            :: usktrk(:)

integer, intent(out)                       :: ncomp

integer, intent(in)                        :: numsp

real*8, intent(in)                         :: c(numsp)

real*8, intent(in)                         :: g(numsp)

real*8, intent(in)                         :: alpha(numsp)

logical, intent(out)                       :: iserror 
 
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
 
 

!%----------------------------------------------------------
call compute_usktrk_ &
   (this%pp, &
    usktrk, &
    ncomp, &
    c, &
    g, &
    alpha, &
    numsp, &
    iserror)
!%----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_usktrk2_chemsys_classic &
   (this, &
    usktrk, &
    ncomp, &
    sktrk, &
    numsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemclassic), intent(in)       :: this

integer, intent(out)                             :: ncomp

integer, intent(in)                              :: numsp

real*8, intent(in)                               :: sktrk(numsp)

real*8, pointer                                  :: usktrk(:)

logical, intent(out)                             :: iserror 
 
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
 
 

!%---------------------------------------------------------
call compute_usktrk_ &
   (this%pp, &
    usktrk, &
    ncomp, &
    sktrk, &
    numsp, &
    iserror)
!%----------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%
!%
!%************************************************************
subroutine compute_umob_chemsys_classic &
   (this, &
    umob, &
    naqpri, &
    nmobph, &
    numsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
 
type (t_chemicalsystemclassic)       :: &
 this
integer                               :: &
 nmobph, &
 naqpri, &
 numsp
real*8                                :: &
 c(numsp)
real*8, pointer                       :: &
 umob(:,:)
logical                               :: &
 iserror
!%------------------------------------------------------------
 call compute_umob_ &
   (this%pp, &
    umob, &
    naqpri, &
    nmobph, &
    numsp, &
    c, &
    iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%   Compute Uads*cads
!%
!%************************************************************
subroutine compute_uads_chemsys_classic &
   (this, &
    uads, &
    naqpri, &
    c, &
    numsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
 
type (t_chemicalsystemclassic)       :: &
 this
real*8                                :: &
 c(:)
real*8, pointer                       :: &
 uads(:)
integer                               :: &
 naqpri, &
 numsp
logical                               :: &
 iserror
!%--------------------------------------------------------------
call compute_uads_ &
   (this%pp, &
    uads, &
    naqpri, &
    c, &
    numsp, &
    iserror)
!%--------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine make_lin_trf_vector_chemsys_classic &
   (this, &
    vnew, &
    vold, &
    iserror, &
	isueq)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemclassic), intent(in)  :: this

real*8, pointer                             :: vnew(:)

real*8, intent(in)                          :: vold(:)

logical, intent(out)                        :: iserror 

logical, intent(in), optional               :: isueq
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)             :: &
 msg
integer                        :: &
 ndim1 
logical                        :: &
 haveisueq
real*8, pointer                :: &
 u(:,:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
ndim1=size(vold)
if (ndim1/=this%pp%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
haveisueq=present(isueq)
!%------------------------------------------------------------
if (haveisueq.and..not.isueq) then
 u => this%pp%u
else
 u => this%pp%ueq 
end if
!%------------------------------------------------------------
call check_pointer_ (vnew,this%pp%numaqprisp,.true.)
vnew=matmul(u,vold)
!%------------------------------------------------------------
u => null ()
!%------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: make_lin_trf_'
print *, msg
print *,'***********************'
iserror=.true.
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine make_lin_trf_array_chemsys_classic &
   (this, &
    anew, &
    aold, &
    iserror, &
	isueq)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemclassic), intent(in)  :: this

real*8, pointer                             :: anew(:,:)

real*8, intent(in)                          :: aold(:,:)

logical, intent(out)                        :: iserror 

logical, intent(in), optional               :: isueq
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                        :: &
 ndim1, &
 ndim2
character(len=100)             :: &
 msg 
real*8, pointer                :: &
 u(:,:) => null ()
logical                        :: &
 haveisueq
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
ndim1=size(aold,1)
ndim2=size(aold,2)
if (ndim2==0) then
 msg='Error, number of columns of input matrix = 0'
 goto 10
end if
!%----------------------------------------------------------
if (ndim1/=this%pp%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------
!% Check optional arguments 
!%------------------------------------------------------------
haveisueq=present(isueq)
!%------------------------------------------------------------
if (haveisueq.and..not.isueq) then
 u => this%pp%u
else
 u => this%pp%ueq 
end if
!%----------------------------------------------------------
call check_pointer_ (anew,this%pp%numaqprisp,ndim2,.true.)
anew=matmul(u,aold)
!%----------------------------------------------------------
u => null ()
!%----------------------------------------------------------
return
!%-----------------------------------------------------------
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: make_lin_trf_'
print *,msg
print *,'***********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dusktrk1_chemsys_classic &
   (this, &
    dusktrk, &
    ncomp, &
    c, &
    alpha, &
    numsp, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemclassic), intent(in) :: this

real*8, pointer                :: dusktrk(:,:)

integer, intent(in)            :: numsp

real*8, intent(in)             :: c(numsp)

real*8, intent(in)             :: alpha(numsp)

integer, intent(out)           :: ncomp

logical, intent(out)           :: iserror 
 
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
call compute_dusktrk_ &
   (this%pp, &
    dusktrk, &
    ncomp, &
    c, &
    alpha, &
    numsp, &
    iserror)
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%****************Public subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%
!%
!%************************************************************
subroutine compute_dusktrk2_chemsys_classic &
   (this, &
    dusktrk, &
    naqpri, &
    dsktrk, &
    numsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
 
type (t_chemicalsystemclassic)      :: &
 this
integer                              :: &
 numsp, &
 naqpri
real*8                               :: &
 dsktrk(numsp,naqpri)
real*8, pointer                      :: &
 dusktrk(:,:)
logical                              :: &
 iserror
!%----------------------------------------------------------------
call compute_dusktrk_ &
   (this%pp, &
    dusktrk, &
    naqpri, &
    dsktrk, &
    numsp, &
    iserror)
!%------------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_duads_chemsys_classic &
   (this, &
    duads, &
    n1, &
    n2, &
    numsp, &
    c, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute Uadsjth*dcads/dc1ith
!
!   $Arguments:
!
 
type (t_chemicalsystemclassic), intent(in):: this

integer, intent(in)                   :: numsp

real*8, intent(in)                    :: c(numsp)

real*8, pointer                       :: duads(:,:)

integer, intent(out)                  :: n1

integer, intent(out)                  :: n2

logical, intent(out)                  :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                               :: &
 i, &
 j, &
 numreact, &
 ipos1, &
 ipos2, &
 ireact, &
 isp, &
 icw, &
 ispw
integer, pointer                      :: &
 idreact(:) => null ()
real*8, pointer                       :: &
 dc(:,:) => null (), &
 dg(:,:) => null (), &
 g(:) => null (), &
 dionstr(:) => null (), &
 dcloc(:) => null (), &
 cloc(:) => null ()
real*8                                :: &
 ionstr, &
 cd, &
 gd
character(len=100)                    :: &
 msg
logical                               :: &
 isanomalous, &
 isupmxitergam 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
if (numsp.ne.this%pp%numsp) then
 msg='Error, diferent number of species defined in the'// &
     ' chemical system'
 goto 10
end if
!%------------------------------------------------------
n1=this%pp%numaqprisp
n2=this%pp%numaqprisp
!%------------------------------------------------------
call check_pointer_ (duads,n1,n1,.true.)
call check_pointer_ (dc,this%pp%numsp,this%pp%numaqprisp,.true.)
call check_pointer_ (dg,this%pp%numsp,this%pp%numaqprisp,.true.)
call check_pointer_ (g,this%pp%numsp,.true.)
call check_pointer_ (cloc,this%pp%numsp,.true.)
cloc=c
g=1.0d0
!%------------------------------------------------------
do i=1,this%pp%numaqprisp
 dc(this%pp%idaqprisp(i),i)=1.0d0
end do
!%-----------------------------
call compute_secondaries_ &
   (this%pp, &
    cloc, &
    g, &
    dc, &
    dg, &
    this%pp%aqphindex, &
    0, &
    ionstr, &
    dionstr, &
	1.0d0, &
	.true., &
    isanomalous, &
    isupmxitergam, &
    .true., &
    msg, &
    iserror)
 
      if (iserror.or.isupmxitergam) goto 20
 
!%----------------------------
do i=1,this%pp%numsurf
 
 
 call get_idreaction_ &
     (this%pp, &
      idreact, &
      numreact, &
      this%pp%aqphindex, &
      i, & 
      .true.)
 
 do j=1,numreact
 
  ireact = idreact(j)
  isp = this%pp%idreactsp(ireact)
  cd = c(isp)
  gd = g(isp)
 
  call compute_dx_ &
  (this%pp%preaction(ireact)%ptr, &
   dcloc, &
   c(this%pp%idaqprisp), &
   g(this%pp%idaqprisp), &
   gd, &
   dc(this%pp%idaqprisp,:), &
   dg(this%pp%idaqprisp,:), &
   dg(isp,:), &
   iserror, &
   cd)
 
   if(iserror) goto 20
   dc(isp,:) = dcloc
 
 end do
 
call get_iposspsurf_ (this%pp, i, ipos1, ipos2)
 
duads=duads+matmul(this%pp%u(:,ipos1:ipos2),dc(ipos1:ipos2,:))
 
end do
!%-----------------------------------------------------------
20 continue 
call check_pointer_ (idreact,1,.false.)
call check_pointer_ (dcloc,1,.false.)
call check_pointer_ (g,1,.false.)
call check_pointer_ (cloc,1,.false.)
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dg,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
!%-----------------------------------------------------------
return
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_duads_'
print *,msg
print *,'**************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%   Compute dcmob/dc1
!%
!%************************************************************
subroutine compute_dcmob1_chemsys_classic &
   (this, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    iserror, &
    dg, &
    g)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
 
 
type (t_chemicalsystemclassic)       :: &
 this
integer                               :: &
 numsp
real*8                                :: &
 c(numsp)
real*8, pointer                       :: &
 dcmob(:,:)
integer, intent(out)                  :: &
 nrow, &
 ncol, &
 nmobph
logical                               :: &
 iserror
real*8, pointer, optional             :: &
 dg(:,:), &
 g(:)
!%------------------------------------------------------------
call compute_dcmob_ &
   (this%pp, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    iserror, &
    dg, &
    g)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%   Compute dcmob/dc1
!%
!%************************************************************
subroutine compute_dcmob2_chemsys_classic &
   (this, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
 
 
type (t_chemicalsystemclassic)       :: &
 this
integer                               :: &
 numsp, &
 naqpri
real*8                                :: &
 dc(numsp,naqpri)
real*8, pointer                       :: &
 dcmob(:,:)
integer, intent(out)                  :: &
 nrow, &
 ncol, &
 nmobph
logical                               :: &
 iserror
 
character(len=300)                    :: &
 msg
!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
call compute_dcmob_ &
   (this%pp, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    iserror)
!%-----------------------------------------------------------
return
 
10 print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dcmob_'
print *,msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%   Compute U*dcmob/dc1
!%
!%************************************************************
subroutine compute_dumob1_chemsys_classic &
   (this, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
 
 
type (t_chemicalsystemclassic)       :: &
 this
integer                               :: &
 numsp, &
 nmobph, &
 nrow, &
 ncol
real*8                                :: &
 c(numsp)
real*8, pointer                       :: &
 dumob(:,:)
logical                               :: &
 iserror
 
character(len=300)                    :: &
 msg
!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
call compute_dumob_ &
   (this%pp, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    iserror)
!%------------------------------------------------------
return
10 print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dumob_'
print *,msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%   Compute dcmob/dc1
!%
!%************************************************************
subroutine compute_dcads_chemsys_classic &
   (this, &
    dcads, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
 
 
type (t_chemicalsystemclassic)       :: &
 this
integer                               :: &
 numsp, &
 naqpri
real*8                                :: &
 dc(numsp,naqpri)
real*8, pointer                       :: &
 dcads(:,:)
integer, intent(out)                  :: &
 nrow, &
 ncol
logical                               :: &
 iserror
 
character(len=100)                    :: &
 msg
!%------------------------------------------------------
call compute_dcads_ &
   (this%pp, &
    dcads, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    iserror)
!%-----------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%
!%   Compute U*dcmob/dc1
!%
!%************************************************************
subroutine compute_dumob2_chemsys_classic &
   (this, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
 
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
 
 
 
type (t_chemicalsystemclassic)         :: &
 this
integer                               :: &
 numsp, &
 naqpri
real*8                                :: &
 dc(numsp,naqpri)
real*8, pointer                       :: &
 dumob(:,:)
integer, intent(out)                  :: &
 nrow, &
 ncol, &
 nmobph
logical                               :: &
 iserror
 
character(len=300)                    :: &
 msg
!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
call compute_dumob_ &
   (this%pp, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    iserror)
!%-----------------------------------------------------------
return
 
10 print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dumob_'
print *,msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_file_chemsys_classic &
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
 
type (t_chemicalsystemclassic), intent(in)  :: this

integer, intent(in)                         :: ioutput

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
 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

!%------------------------------------------------------------
iserror=.false.
!%------------------------------------------------------------
call write_ (this%pp,ioutput,iserror)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_chemicalsystemclassic
