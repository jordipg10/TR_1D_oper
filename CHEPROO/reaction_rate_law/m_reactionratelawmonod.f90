module m_reactionratelaw_monod
!-------------------------------------------------------------------------
!
!   $Description: Represent one reaction rate law
!
!   $Use: use m_species
! use m_parentreactionratelaw
! use flib_xpath
! use flib_sax
! use m_general_tools_cheproo
! use m_constants
!
!   $Author: Sergio Andrés Bea Jofré 
!
!   $License:
!
!-------------------------------------------------------------------------
use m_species
use m_parentreactionratelaw
use flib_xpath
use flib_sax
use m_general_tools_cheproo
use m_constants_cheproo
!%------------------------------------------------------------------
!%------------------------------------------------------------------
private                      ::
!%------------------------------------------------------------------
!%------------------------------------------------------------------
public                      :: &
create_ &         ! Create reaction rate law monod object 
,read_xml_ &      ! Read object from xml file 
,destroy_ &       ! Destroy reaction rate law monod object
,set_ &
,compute_rk_ &    ! Compute reaction rate 
,compute_drk_ &   
,write_ &
,rewrite_ &
,assignment(=)
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
private :: &
read_xml_loc_ &
,begin_element_handler 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Type definition 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
type, public:: t_reactionratelaw_monod
 
private                                ::
 
type(t_parentreactionratelaw), pointer :: pp       ! Parent reaction rate law 
 
real*8                                 :: mxrate   ! Maximum effective utilization rate

real*8                                 :: minrate  ! Minimum rate

real*8                                 :: theta    ! exponent coefficient for omega

real*8                                 :: eta      ! exponent coefficient of equilibrium term

real*8, pointer, dimension(:,:)        :: ki       ! Reaction rates [numsp,numterm]
 
 
end type t_reactionratelaw_monod
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface create_
 
module procedure create_rrlmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface read_xml_
 
module procedure read_xml_rrlmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface destroy_
 
module procedure destroy_rrlmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_
 
module procedure set_rrlmonod
module procedure set_parent_rrlmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_rk_
 
module procedure compute_rk_rrlmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_drk_
 
module procedure compute_drk_rrlmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface write_
 
module procedure write_rrlmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface rewrite_
 
module procedure rewrite_rrlmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface assignment(=)
 
module procedure copy_rrlmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Privates subroutines
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface read_xml_loc_
 
module procedure read_xml_loc_rrllmonod
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface begin_element_handler
 
module procedure begin_element_handler
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
contains
!%-----------------------------------------------------------------
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_rrlmonod &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create reaction rate law object monod 
!
!   $Arguments:
!
 
type(t_reactionratelaw_monod), intent(inout):: this 
 
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
this%pp => null()
!%------------------------------------------------------------
this%mxrate=0.0d0
this%minrate=0.0d0
this%eta=0.0d0
this%theta=0.0d0
!%------------------------------------------------------------
!% Nullify pointers 
!%------------------------------------------------------------
this%ki => null ()
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_rrlmonod &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy reaction rate law object 
!
!   $Arguments:
!
 
type(t_reactionratelaw_monod), intent(inout):: this 
 
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
this%pp => null()
!%------------------------------------------------------------
this%mxrate=0.0d0
this%minrate=0.0d0
this%eta=0.0d0
this%theta=0.0d0
!%------------------------------------------------------------
!% Deallocate 
!%------------------------------------------------------------
call check_pointer_ (this%ki,1,1,.false.)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_rrlmonod &
  (this, &
   pp, &
   namefile, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read reaction rate law from xml file. 
!
!   $Arguments:
!
 
type(t_reactionratelaw_monod), intent(inout)    :: this

type(t_parentreactionratelaw), intent(in),target:: pp

character(len=*), intent(in)                    :: namefile

logical, intent(out)                            :: iserror 
 
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
 iostat
type(xml_t):: &
 fxml
type(dictionary_t)                    :: &
 attributes
character(len=100)                    :: &
 name, &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------
this%pp => pp
!%--------------------------------------------------------
!% Open xml file 
!%--------------------------------------------------------
call open_xmlfile(namefile, fxml, iostat)
if (iostat /= 0) goto 10
!%--------------------------------------------------------
call xml_parse(fxml,begin_element_handler=begin_element_handler)
!%--------------------------------------------------------
! End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%--------------------------------------------------------
call read_xml_loc_(name,attributes,this,msg,iserror)
if (iserror) goto 10
!%--------------------------------------------------------
return
10 continue 
print *,'**********************'
print *,'Reaction Rate Law:'
print *,'Name:',this%pp%name
print *,'Service: read_xml_'
print *, msg
print *,'**********************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_rrlmonod &
  (this, &
   pp, &
   mxrate, &
   minrate, &
   theta, &
   eta, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the reaction rate law type monod
!
!   $Arguments:
!
 
type(t_reactionratelaw_monod), intent(inout)      :: this      ! Type reaction rate law monod

type(t_parentreactionratelaw), intent(in), target :: pp        ! Pointer to parent reaction rate law

real*8, intent(in)                                :: mxrate    ! maximum rate

real*8, intent(in)                                :: minrate   ! minimun rate

real*8, intent(in)                                :: theta     ! equilibrium exponent for omega 

real*8, intent(in)                                :: eta       ! equilibrium exponent

logical, intent(out)                              :: iserror   ! iserror=true, then there was an error 
 
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
 i
logical                               :: &
 isbe
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
this%pp => pp
!%------------------------------------------------------------
this%mxrate=mxrate
this%minrate=minrate
this%eta=eta
this%theta=theta
!%------------------------------------------------------------
!% Allocate attrsp array
!%------------------------------------------------------------
if (this%pp%numsp>0) then
  call check_pointer_ (this%ki,this%pp%numsp,this%pp%numterm,.true.)
  this%ki=this%pp%attrsp
end if 
!%------------------------------------------------------------
!% Check term type 
!%------------------------------------------------------------
do i=1,this%pp%numterm
 isbe=.false.
 select case (this%pp%typeterm(i))
 case ('monod','nonmonod','inhibitor')
  isbe=.true.
 end select
 if (.not.isbe) then
   msg='Error, term type non recognized:'
   call add_ (msg,this%pp%typeterm(i))
   goto 10
 end if
end do
!%------------------------------------------------------------
!% Update according the reference temperature
!%------------------------------------------------------------
call update_ (this%pp,this%pp%tempref,iserror)
if (iserror) goto 10
!%------------------------------------------------------------
return
 
10 continue 
print *,'*********************'
print *,'Reaction Rate Law:'
print *,'Namë:',this%pp%name
print *,'Service:  set_'
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
subroutine set_parent_rrlmonod &
  (this, &
   pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set parent of reaction rate law
!
!   $Arguments:
!
 
type(t_reactionratelaw_monod), intent(inout)        :: this

type(t_parentreactionratelaw), intent(in), target   :: pp 
 
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
subroutine write_rrlmonod &
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
 
type(t_reactionratelaw_monod), intent(in), target :: this

integer, intent(in)                               :: ioutput 

logical, intent(out)                              :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                :: &
 i, &
 j 
real*8, pointer        :: &
 value => null ()
character(len=100), pointer :: &
 name => null ()
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.  
!%------------------------------------------------------------
write (ioutput,*) "--------------------------------------"
write (ioutput,*) "Chemical Reaction Rate Law Information"
write (ioutput,*) "        (monod type)                  "
write (ioutput,*) "--------------------------------------"
write (ioutput,4) this%pp%name
write (ioutput,*) "--------------------------------------"
write (ioutput,2) "Activation energy:",this%pp%ea
write (ioutput,5) "Maximum substrate rate:",this%mxrate
write (ioutput,5) "Minimum substrate rate:",this%minrate
write (ioutput,5) "Equil. exponent for omega ",this%theta
write (ioutput,5) "Equil. exponent",this%eta
write (ioutput,*) "--------------------------------------"
do i=1,this%pp%numterm
 select case (this%pp%typeterm(i))
 case ('monod')
  write (ioutput,1) "Monod Terms"
  do j=1,this%pp%numsp
    value => this%pp%attrsp(j,i) 
    if (value/=0.0d0) then
     name => this%pp%pspecies(j)%ptr%name
	 write (ioutput,3) name,"half sat.const:",value
    end if
  end do
 case ('nonmonod')
  write (ioutput,*) "--------------------------------------"
  write (ioutput,1) "No-Monod Terms"
  do j=1,this%pp%numsp
    value => this%pp%attrsp(j,i) 
    if (value/=0.0d0) then
     name => this%pp%pspecies(j)%ptr%name
	 write (ioutput,3) name,"exponent:",value
    end if
  end do
 case ('inhibitor')
  write (ioutput,*) "--------------------------------------"
  write (ioutput,1) "Inhibitor Terms"
  do j=1,this%pp%numsp
    value => this%pp%attrsp(j,i) 
    if (value/=0.0d0) then
	 name => this%pp%pspecies(j)%ptr%name
     write (ioutput,3) name,"inhib.const.:",value
    end if
  end do
 end select
 
 
end do
write (ioutput,*) "--------------------------------------"
!%------------------------------------------------------------
!% Nullify local pointers 
!%------------------------------------------------------------
name => null ()
value => null ()
!%------------------------------------------------------------
return
 
1 format (a15)
4 format (a40)
5 format (a25,e15.7)
2 format(a25,f6.2)
3 format(a10,a15,e15.7)
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_rk_rrlmonod &
  (this, &
   rk, &
   omega, &
   c, &
   nsp, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute rk according monod equation
!
!   $Arguments:
!
 
type(t_reactionratelaw_monod),intent(in) :: this

real*8, intent(out)                      :: rk

real*8, intent(in)                       :: omega

integer, intent(in)                      :: nsp

real*8, intent(in)                       :: c(nsp)

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
character(len=100)                       :: &
 msg
integer                                  :: &
 isp, &
 iterm, &
 iou 
real*8                                   :: &
 t, &
 ki, &
 fomega, &
 dfomega, &
 omtol1, &
 omtol2, &
 ompsi, &
 thresh
real*8, parameter                        :: &
 r0=0.0d0, &
 r1=1.0d0
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
!% Initialice local variables 
!%------------------------------------------------------------
omtol1=r0
omtol2=r0
ompsi=r0
fomega=r1
dfomega=r0
thresh=r1
!%------------------------------------------------------------
!% Initialice variables 
!%------------------------------------------------------------
rk=r0
!%------------------------------------------------------------
if (this%pp%numsp>0.and.nsp/=size(this%ki,1)) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------
!% Compute rk
!%------------------------------------------------------------
rk=this%pp%exp_ea_rt*this%mxrate
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
if (this%pp%numsp>0) then   
  do isp=1,nsp
 
     do iterm=1,this%pp%numterm
       ki=this%ki(isp,iterm)
       t=c(isp)
       if (ki/=r0) then
          select case (this%pp%typeterm(iterm))
          case ('monod')
             rk=rk*(t/(t+ki))
          case ('nonmonod')
             rk=rk*(t**ki)
          case ('inhibitor')
             rk=rk*(ki/(t+ki))
          end select
       end if
     end do
 
  end do
end if 
!%------------------------------------------------------------
!% Compute equilibrium term 
!%------------------------------------------------------------
if (this%eta/=r0) then
  call omegaterm_ (fomega,dfomega,omega,this%theta,this%eta,thresh,omtol1,omtol2,ompsi,iou)
!%------------------------------------------------------------
!% Add equilibrium term 
!%------------------------------------------------------------
  rk = rk * fomega 
end if
!%------------------------------------------------------------
!% the rate is restricted to minimun rate
!%------------------------------------------------------------
rk=rk-this%minrate
!%------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Reaction Rate Law:'
print *,'Name:',this%pp%name
print *,'Service: compute_rk_'
print *, msg
print *,'**********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_drk_rrlmonod &
  (this, &
   drk, &
   omega, &
   domega, &
   c, &
   dc, &
   nsp, &
   ndimder, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute derivate of reaction rate with respect to primary concentration
!
!   $Arguments:
!
 
type(t_reactionratelaw_monod), intent(in)  :: this

real*8, pointer                            :: drk(:)

real*8, intent(in)                         :: omega

integer, intent(in)                        :: nsp

integer, intent(in)                        :: ndimder

real*8, intent(in)                         :: c(nsp)

real*8, intent(in)                         :: domega(ndimder)

real*8, intent(in),target                  :: dc(nsp,ndimder)

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
 
character(len=100)                         :: &
 msg
real*8                                     :: &
 aux, &
 rk, &
 ki, &
 t, &
 omtol1, &
 omtol2, &
 ompsi, &
 fomega, &
 dfomega, &
 thresh 
real*8, pointer                            :: &
 dt(:) => null ()
integer                                    :: &
 isp, &
 iterm, &
 iou, &
 i
real*8, parameter                          :: &
 r0=0.0d0, &
 r1=1.0d0 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
if (this%pp%numsp>0.and.nsp/=size(this%ki,1)) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------
!% Initialice local variables 
!%------------------------------------------------------------
omtol1=r0
omtol2=r0
ompsi=r0
fomega=r1
dfomega=r0
thresh=r1
!%------------------------------------------------------------
! Compute rk
!%------------------------------------------------------------
call compute_rk_(this,rk,omega,c,nsp,iserror)
if (iserror) goto 10
!%------------------------------------------------------------
call check_pointer_ (drk,ndimder,.true.)
!%------------------------------------------------------------
! Compute drk
!%------------------------------------------------------------
if (this%pp%numsp>0) then    
   do isp=1,nsp
 
     do iterm=1,this%pp%numterm
        ki=this%ki(isp,iterm)
        t=c(isp)
        dt => dc(isp,:)
        if (ki/=r0) then
          select case (this%pp%typeterm(iterm))
          case ('monod')
             aux=ki/(t*ki+t*t)
             drk=drk+aux*dt
          case ('nonmonod')
             aux=ki/t
             drk=drk+aux*dt
          case ('inhibitor')
             aux=-(r1/(ki+t))
             drk=drk+aux*dt
          end select
         end if
     end do
 
   end do
end if 
!%------------------------------------------------------------
!% Add equilibrium term
!%------------------------------------------------------------
if (this%eta/=r0) then
  call omegaterm_ (fomega,dfomega,omega,this%theta,this%eta,thresh,omtol1,omtol2,ompsi,iou)
!%------------------------------------------------------------
!% Add equilibrium term 
!%------------------------------------------------------------
  if (fomega/=r0) then 
    drk = drk + (dfomega/fomega) * domega
  end if  
end if 
!%------------------------------------------------------------
!%------------------------------------------------------------
!%------------------------------------------------------------
drk=rk*drk
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
dt => null ()
!%------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Reaction Rate Law:'
print *,'Name:',this%pp%name
print *,'Service: compute_drk_'
print *, msg
print *,'**********************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_rrlmonod &
  (copied, &
   this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy the reaction rate law object
!
!   $Arguments:
!
 
type(t_reactionratelaw_monod), intent(in)  :: this

type(t_reactionratelaw_monod), intent(out) :: copied 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer             :: &
 ndim1, &
 ndim2
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
copied%mxrate=this%mxrate
copied%minrate=this%minrate
copied%eta=this%eta
copied%theta=this%theta
!%-----------------------------------------------------------
if (associated(this%ki)) then
 ndim1=size(this%ki,1)
 ndim2=size(this%ki,2)
 call check_pointer_ (copied%ki,ndim1,ndim2,.true.)
 copied%ki=this%ki
end if
!%-----------------------------------------------------------
return
end subroutine
!%***********************************************************
!%***************Public subroutine***************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
subroutine rewrite_rrlmonod &
  (this, &
   namesp, &
   nsp, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_reactionratelaw_monod), intent(inout), target   :: this

integer, intent(in)                                    :: nsp

character(len=*), intent(in), dimension(nsp)           :: namesp

logical, intent(out)                                   :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer               :: &
 i, &
 j, &
 isps
character(len=100)    :: &
 msg 
character(len=100), pointer :: &
 name => null ()
!-------------------------------------------------------------------------
!
!   $code
!
!%--------------------------------------------------------------
iserror=.false.
msg=''
!%--------------------------------------------------------------
!% Warning, if the number of species is 0, in the reaction rate 
!% law, then return 
!%--------------------------------------------------------------
if (this%pp%numsp==0) return 
!%--------------------------------------------------------------
!% Check if the number of species in the reaction rate law is >
!% than the total of species in the list 
!%--------------------------------------------------------------
if (nsp<this%pp%numsp) then
 msg='Error in number of species'
 goto 10
end if 
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!%---------------------------------------------------------------
call check_pointer_ (this%ki,nsp,this%pp%numterm,.true.)
!%---------------------------------------------------------------
isps=0
do i=1,nsp
 
  do j=1,this%pp%numsp
    name => this%pp%pspecies(j)%ptr%name
    if(name==namesp(i)) then
       isps=isps+1
       this%ki(i,:) = this%pp%attrsp(j,:)
	   exit  
    end if
    
  end do
  
end do
!%---------------------------------------------------------------
!% Count the number of founded species 
!%---------------------------------------------------------------
if (isps/=this%pp%numsp) then
 msg='Catalytic species defined in the reaction rate law not'// &
     ' defined in the species list'
 iserror=.true.
 goto 20
end if 
!%---------------------------------------------------------------
!% Nullify local pointers 
!%---------------------------------------------------------------
20 continue 
name => null () 
!%--------------------------------------------------------------- 
return
 
10 continue 
print *,'*******************'
print *,'Reaction:'
print *,'Name:', this%pp%name
print *,'Service: rewrite_'
print *, msg
print *,'*******************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine begin_element_handler (name, attributes)
 
 
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
 
 
 call read_xml_loc_ (name, attributes)
 
return
end subroutine
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc_rrllmonod &
 (name, &
  attributes, &
  this, &
  msg, &
  iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
character(len=*), intent(in)                           :: name

type(dictionary_t), intent(in)                         :: attributes

type (t_reactionratelaw_monod), intent(inout), optional:: this

logical, intent(out), optional                         :: iserror

character(len=*), optional                             :: msg 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                         :: &
 reallocal(1)
integer, save                  :: &
 iexp, &
 isp
integer, pointer, save         :: &
 numsp(:)
real*8, save                   :: &
 ea, &
 mxrate, &
 minrate, &
 eta, &
 theta 
real*8, pointer, save          :: &
 attrterm(:,:), &
 attrspterm(:,:)
logical                        :: &
 havethis, &
 haveiserror, &
 iserror1
logical, save                  :: &
 isareadep
integer                        :: &
 i, &
 n, &
 status
character(len=100)                :: &
 id
character(len=100), pointer,save  :: &
 namesp(:,:)
character(len=100), save          :: &
 namerrlaw
character(len=100), pointer, save :: &
 typeterm(:)
integer, parameter                :: &
 mxdim=20 
!-------------------------------------------------------------------------
!
!   $code
! 
havethis=present(this)
haveiserror=present(iserror)
 
select case (havethis)
 
case (.true.)
 iserror1=.false.
 if (haveiserror) iserror=.false. 
 call set_ &
  (this%pp, &
   namerrlaw, &
   numsp(1:iexp), &
   attrterm, &
   attrspterm(:,1:iexp), &
   namesp(1:mxdim,1:iexp), &
   typeterm, &
   iexp, &
   0, &
   mxdim, &
   ea, &
   isareadep, &
   iserror1)
   
 if (iserror1.and.haveiserror) then
   iserror=.true.
   return 
 else if(iserror1.and..not.haveiserror) then
   stop 
 end if
 
 call set_ (this,this%pp,mxrate,minrate,theta,eta,iserror1)
 
 if (iserror1.and.haveiserror) then
   iserror=.true.
   return 
 else if(iserror1.and..not.haveiserror) then
   stop 
 end if
 
 call check_pointer_ (numsp,1,.false.)
 call check_pointer_ (attrspterm,1,1,.false.)
 call check_pointer_ (namesp,1,1,.false.)
 call check_pointer_ (typeterm,1,.false.)
 
 
case default
 
 select case (name)
 
 case ('reactionratelaw')
 
     call check_pointer_ (numsp,mxdim,.true.)
     call check_pointer_ (attrspterm,mxdim,mxdim,.true.)
     call check_pointer_ (namesp,mxdim,mxdim,.true.)
     call check_pointer_ (typeterm,mxdim,.true.)
     iexp=0
	 mxrate=0.0d0
	 minrate=0.0d0
	 eta=0.0d0
	 theta=0.0d0 
	 isareadep=.false. 
     id=''
     call get_value (attributes,"name", id, status)
     namerrlaw=id
     id=''
     call get_value (attributes,"ea", id, status)
     n=0
	 reallocal=0.0d0
     call build_data_array (id,reallocal,n)
     ea=reallocal(1)
     id=''
     call get_value (attributes,"attr1", id, status)
     if (status==0) then 
	    n=0
		reallocal=0.0d0
        call build_data_array (id,reallocal,n)
        mxrate=reallocal(1)
	 end if
     id=''
     call get_value (attributes,"attr2", id, status)
     if (status==0) then 
	    n=0
        reallocal=0.0d0
		call build_data_array (id,reallocal,n)
        minrate=reallocal(1)
	 end if 
     id=''
     call get_value (attributes,"attr3", id, status)
     if (status==0) then  
	    n=0
        reallocal=0.0d0
		call build_data_array (id,reallocal,n)
        theta=reallocal(1)
	 end if
	 id=''
     call get_value (attributes,"attr4", id, status)
     if (status==0) then  
	    n=0
        reallocal=0.0d0
		call build_data_array (id,reallocal,n)
        eta=reallocal(1)
	 end if 
	 id=''
     call get_value (attributes,"areadep", id, status)
     select case (id)
     case ('yes','YES','y','Y')
       isareadep=.true.
     end select 
 
 case ('term')
 
     iexp=iexp+1
     isp=0
     id=''
     call get_value (attributes,"type", id, status)
     typeterm(iexp)=id
 
 case ('species')
 
     isp=isp+1
     numsp(iexp)=isp
     call get_value (attributes,"name", id, status)
     namesp(isp,iexp)=id
     call get_value (attributes,"attr1", id, status)
     n=0
     call build_data_array (id,reallocal,n)
     attrspterm(isp,iexp)=reallocal(1)
 
 case default
 
 end select
 
end select
!%-------------------------------------------------------------
 
return
 end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_reactionratelaw_monod
