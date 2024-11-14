module m_reactionratelaw_lasaga
!-------------------------------------------------------------------------
!
!   $Description: Represent one reaction rate law
!
!   $Use:
!
!   $Author: Sergio Andrés Bea Jofré 
!
!   $License:
!
!-------------------------------------------------------------------------
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_species
use m_parentreactionratelaw
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%------------------------------------------------------------------
!%------------------------------------------------------------------
private                     ::
!%------------------------------------------------------------------
!%------------------------------------------------------------------
public                      :: &
create_ &          ! Create the object 
,destroy_ &        ! Destroy the reaction rate law object. 
,read_xml_ &       ! Read the reaction rate law object from xml file. 
,set_ &            ! Set attributes in the reaction rate law object 
,compute_rk_ &     ! Compute rk 
,compute_drk_ &    ! Compute drk 
,write_ &          ! Write in ascii the different attributes encapsulated in the object.
,rewrite_ &        ! 
,assignment(=)     ! Copy the reaction rate law object in other reaction rate law object. 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
private:: &
begin_element_handler &
,read_xml_loc_
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Type definition 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
type, public:: t_reactionratelaw_lasaga
 
private                                     ::

type(t_parentreactionratelaw), pointer      :: pp         ! Parent reaction rate law
 
real*8                                      :: omegathr   ! Supersaturation log (q/k) necessary for precipitation 

real*8, pointer, dimension(:,:)             :: pcat       ! Stoichiometric matrix of catalytic species [numsp,numterm]                     
 
end type t_reactionratelaw_lasaga
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Public services 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface create_
 
module procedure create_rrllasaga
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface destroy_
 
module procedure destroy_rrllasaga
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface read_xml_
 
module procedure read_xml_rrllasaga
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_
 
module procedure set_rrllasaga
module procedure set_parent_rrllasaga
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_rk_
 
module procedure compute_rk_rrllasaga
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface compute_drk_
 
module procedure compute_drk_rrllasaga
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface write_
 
module procedure write_rrllasaga
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface rewrite_
 
module procedure rewrite_rrllasaga
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface assignment(=)
 
module procedure copy_rrllasaga
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Private Services
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface begin_element_handler
 
module procedure begin_element_handler
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface read_xml_loc_
 
module procedure read_xml_loc_rrllasaga
 
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
subroutine create_rrllasaga &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_reactionratelaw_lasaga):: this 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond: Create reaction rate object (Lasaga type)
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
this%omegathr=0.0d0
!%------------------------------------------------------------
this%pcat => null ()
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_rrllasaga &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy reaction rate law object 
!
!   $Arguments:
!
 
type(t_reactionratelaw_lasaga):: this     ! Type reaction rate law object 
 
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
this%omegathr = 0.0d0
!%------------------------------------------------------------
call check_pointer_ (this%pcat,1,1,.false.)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_rrllasaga &
  (this, &
   pp, &
   namefile, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Read from xml file the reaction rate law object. 
!
!   $Arguments:
!
 
type(t_reactionratelaw_lasaga), intent(inout)        :: this

type(t_parentreactionratelaw), intent(in), target    :: pp 

character(len=*), intent(in)                         :: namefile

logical, intent(out)                                 :: iserror    ! iserror=true, then there was an error 
 
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
call xml_parse(fxml,begin_element_handler = begin_element_handler)
!%--------------------------------------------------------
!% End and close xml file
!%--------------------------------------------------------
call endfile_xmlfile(fxml)
call close_xmlfile(fxml)
!%--------------------------------------------------------
call read_xml_loc_ (name,attributes,this,msg,iserror)
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
subroutine set_rrllasaga &
  (this, &
   pp, &
   omegathr, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the reaction rate law object type lasaga
!
!   $Arguments:
!
 
type(t_reactionratelaw_lasaga), intent(inout)        :: this       ! Type reaction rate law lasaga 

type(t_parentreactionratelaw), intent(in), target    :: pp         ! Parent reaction rate law 

real*8, intent(in)                                   :: omegathr

logical, intent(out)                                 :: iserror    ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer         :: &
 i 
character(len=100)       :: &
 msg  
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
iserror=.false.
msg='' 
!%------------------------------------------------------------
!% Pointer to parent reaction rate object 
!%------------------------------------------------------------
this%pp => pp
!%------------------------------------------------------------
!% Check if theta parameter < 0 
!%------------------------------------------------------------
if (this%pp%numsp>0) then
 call check_pointer_ (this%pcat,this%pp%numsp,this%pp%numterm,.true.)
 this%pcat=this%pp%attrsp
end if
!%------------------------------------------------------------
do i=1,this%pp%numterm 
 if (this%pp%attrterm(3,i)<0.0d0) then 
   msg='Error, kinetic order lower than 0.0d0 in term:'
   call add_ (msg,i) 
   goto 10 
 end if 
end do  
!%------------------------------------------------------------
this%omegathr=omegathr
!%------------------------------------------------------------
!% Update according reference temperature
!%------------------------------------------------------------
call update_ (this%pp,this%pp%tempref,iserror)
if (iserror) goto 10 
!%------------------------------------------------------------
return
10 continue 
print *,'**********************'
print *,'Reaction Rate Law:'
print *,'Name:',this%pp%name
print *,'Service: set_'
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
subroutine set_parent_rrllasaga &
  (this, &
   pp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set parent of reaction rate law
!
!   $Arguments:
!
 
type(t_reactionratelaw_lasaga), intent(inout)        :: this

type(t_parentreactionratelaw), intent(in), target    :: pp
 
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
subroutine write_rrllasaga &
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
 
type(t_reactionratelaw_lasaga), intent(in), target  :: this

integer, intent(in)                                 :: ioutput 

logical, intent(out)                                :: iserror 
 
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
 value1=> null (), &
 value2 => null (), &
 value3 => null (), &
 value4 => null ()
character(len=100), pointer :: &
 name => null ()
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
!%------------------------------------------------------------
name => this%pp%name
value1 => this%pp%ea
value2 => this%pp%tempref
value3 => this%pp%exp_ea_rt
value4 => this%omegathr
!%------------------------------------------------------------
write (ioutput,*) "--------------------------------------"
write (ioutput,*) "Chemical Reaction Rate Law Information"
write (ioutput,*) "         (Lasaga type)                "
write (ioutput,*) "--------------------------------------"
write (ioutput,1) name
write (ioutput,4) "Ea=",value1
write (ioutput,4) "Temp.Ref.=",value2
write (ioutput,4) "e(-Ea/RT)=",value3
write (ioutput,4) "Omegathr=",value4
!%------------------------------------------------------------
!% Check if the reaction rate law depends of the reactive 
!% surface
!%------------------------------------------------------------
if (.not.this%pp%isareadep) then
 write (ioutput,*) "(Reaction rate not depend of the area)"
end if 
do i=1,this%pp%numterm
!%------------------------------------------------------------
 value1 => this%pp%attrterm(1,i)
 value2 => this%pp%attrterm(2,i)
 value3 => this%pp%attrterm(3,i)
 write (ioutput,*) "-------------"
 write (ioutput,5) "Exp. Term=",i
 write (ioutput,*) "-------------"
 write (ioutput,7) "Rate Cte.=",value1
 write (ioutput,4) "Theta=",value2
 write (ioutput,4) "Eta=",value3
 do j=1,this%pp%numsp
   value1 => this%pp%attrsp(j,i)
   if (value1/=0.0d0) then
     name => this%pp%pspecies(j)%ptr%name 
     write (ioutput,6) "Catalyser=",value1,name 
   end if
 end do
end do
!%------------------------------------------------------------
!% Nullify local pointers 
!%------------------------------------------------------------
value1 => null ()
value2 => null ()
value3 => null ()
name => null ()
!%------------------------------------------------------------
return
 
1 format (a40)
4 format(a10,f10.4)
5 format(a10,i5)
6 format(a10,f6.2,a10)
7 format(a10,e10.3,a10)
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
subroutine compute_rk_rrllasaga &
  (this, &
   rk, &
   omega, &
   c, &
   g, &
   nsp, &
   cold, &
   area, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute reaction rate. 
!
!   $Arguments:
!
 
type(t_reactionratelaw_lasaga), intent(in) :: this

real*8, intent(out)                        :: rk

integer, intent(in)                        :: nsp

real*8, intent(in)                         :: omega     ! IAP/Ke
 
real*8, intent(in)                         :: area      ! Area of mineral 

real*8, intent(in)                         :: cold

real*8, intent(in), dimension(nsp)         :: c

real*8, intent(in), dimension(nsp)         :: g

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
integer                :: &
 j, &
 k, &
 iou
real*8                 :: &
 threshmin, &
 rate, &
 theta, &
 eta, &
 coeff, &
 omegathr, &
 omtol1, &
 omtol2, &
 ompsi 
real*8, pointer        :: &
 catterm(:) => null (), &
 fdeltag(:) => null (), &
 dfdeltag(:) => null ()
character(len=100)     :: &
 msg 
logical                :: &
 isbe 
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------------
!% If the mineral area is not considered then always compute the 
!% reaction rate 
!%-------------------------------------------------------------------
if (.not.this%pp%isareadep  &
       .or.  &
 (this%pp%isareadep.and.cold>0.0d0)) then
 isbe=.true. 
else
 isbe=.false.
end if 
!%-------------------------------------------------------------------
!% Zeroing variables 
!%-------------------------------------------------------------------
rk = 0.0d0
omtol1=0.d0
omtol2=0.0d0 
ompsi=0.0d0
!%--------------------------------------------------------------------
!% Check the number of species
!%--------------------------------------------------------------------
if(this%pp%numsp>0.and.nsp/=size(this%pcat,1)) then
 msg='Error in number of species'
 goto 10
end if
!%--------------------------------------------------------------------
!% Control of the lower limit of thresh(i)
!%--------------------------------------------------------------------
if (omtol2>(this%omegathr-1.0d0)) then
 omegathr=1.d0+omtol2
else
 omegathr=this%omegathr
end if
!%--------------------------------------------------------------------
!% Mineral no exists + subsaturation: rk=0 and return
!%--------------------------------------------------------------------
if(omega<omegathr.and..not.isbe)then
 return
else if(omega>1.0d0.and.isbe) then
threshmin = 1.0d0
!%--------------------------------------------------------------------
!% Mineral exists+subsat. and mineral does not exist+supers.
!%--------------------------------------------------------------------
else
 threshmin = omegathr
end if
!%--------------------------------------------------------------------
!% Allocate local pointers 
!%--------------------------------------------------------------------
call check_pointer_ (catterm,this%pp%numterm,.true.)
call check_pointer_ (fdeltag,this%pp%numterm,.true.)
call check_pointer_ (dfdeltag,this%pp%numterm,.true.)
!%--------------------------------------------------------------------
if(this%pp%numsp>0) then  ! If there are catalytic species 
  catterm=matmul(transpose(this%pcat),dlog(c*g))  
  catterm=dexp(catterm)
else ! If there aren't catalytic species 
  catterm=1.0d0 
end if 
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
!%--------------------------------------------------------------------
do k=1,this%pp%numterm
      rate=this%pp%attrterm(1,k)
      theta=this%pp%attrterm(2,k)
      eta=this%pp%attrterm(3,k)
!%--------------------------------------------------------------------
!% Warning, eta always must be >0 !!!!!!!!!!!!!!!!!?????????
!%--------------------------------------------------------------------
      if (eta>0.0d0) then
        call omegaterm_ (fdeltag(k),dfdeltag(k),omega,threshmin,theta,eta,omtol1,omtol2,ompsi,iou)
      else 
        fdeltag(k)=1.0d0
        dfdeltag(k)=0.0d0
      end if
!%-------------------------------------------------------------------
      rk=rk+rate*catterm(k)*fdeltag(k) 
!%-------------------------------------------------------------------
 
end do
!%--------------------------------------------------------------------
!% Multiply for the Arrenius equation 
!%-------------------------------------------------------------------
rk=-rk*this%pp%exp_ea_rt
!%--------------------------------------------------------------------
!% If the raction rate depends of the area [m2/kgw]
!% Now the reaction rate is rk[mol/kgw/s]
!%-------------------------------------------------------------------
if (this%pp%isareadep) rk=rk*area
!%-------------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------------------
call check_pointer_ (catterm,1,.false.)
call check_pointer_ (fdeltag,1,.false.)
call check_pointer_ (dfdeltag,1,.false.)
if (iserror) goto 10
!%-------------------------------------------------------------------
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
subroutine compute_drk_rrllasaga &
  (this, &
   drk, &
   omega, &
   domega, &
   c, &
   dc, &
   g, &
   dg, &
   nsp, &
   ndimder, &
   area, &
   isbe, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute drk 
!
!   $Arguments:
!
 
type(t_reactionratelaw_lasaga), intent(in):: this

integer, intent(in)                       :: nsp

integer, intent(in)                       :: ndimder

real*8, pointer, dimension(:)             :: drk

real*8, intent(in), dimension(ndimder)    :: domega

real*8, intent(in)                        :: omega

real*8, intent(in)                        :: area       ! mineral area in [m2]

logical, intent(in)                       :: isbe       ! If true then mineral is present 

logical, intent(out)                      :: iserror    ! iserror=true, then there was an error 

real*8, intent(in), dimension(nsp)        :: c

real*8, intent(in), dimension(nsp)        :: g

real*8, intent(in), dimension(nsp,ndimder):: dc

real*8, intent(in), dimension(nsp,ndimder):: dg
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8                 :: &
 threshmin, &
 rate, &
 eta, &
 theta, &
 omegathr, &
 omtol1, &
 omtol2, &
 ompsi
integer                :: &
 iou, &
 k, &
 j, &
 i
real*8, pointer        :: &
 fdeltag(:) => null (), &
 dfdeltag(:) => null (), &
 catterm(:) => null (), & 
 dcatterm(:,:) => null (), &
 dact(:,:) => null ()
character(len=100)     :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------------
iou=1
!%------------------------------------------------------------------
!% Allocate reaction rate derivatives pointer and zeroing 
!%------------------------------------------------------------------
call check_pointer_ (drk,ndimder,.true.)
!%------------------------------------------------------------------
if(this%pp%numsp>0.and.nsp/=size(this%pcat,1)) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------------
!% If the mineral is not present then return (only if the reaction 
!% rate depends of the area)
!% If the mineral area is not considered then always compute the 
!% reaction rate 
!%------------------------------------------------------------------
if (this%pp%isareadep.and..not.isbe) return 
!%------------------------------------------------------------------
!% Initialize variables
!%------------------------------------------------------------------
omtol1=0.d0
omtol2=0.0d0
ompsi=0.0d0 
!%------------------------------------------------------------------
!% Allocate local pointers 
!%------------------------------------------------------------------
call check_pointer_ (catterm,this%pp%numterm,.true.)
call check_pointer_ (dcatterm,this%pp%numterm,ndimder,.true.)
call check_pointer_ (fdeltag,this%pp%numterm,.true.)
call check_pointer_ (dfdeltag,this%pp%numterm,.true.)
!%------------------------------------------------------------------
!% Compute activities and derivatives of these activities 
!% (only if there are catalytic species, nsp>0) 
!%------------------------------------------------------------------
if (this%pp%numsp>0) then
  call check_pointer_ (dact,nsp,ndimder,.true.)
  do i=1,nsp
   dact(i,:) = dc(i,:)/c(i) + dg(i,:)/g(i)
  end do
  catterm=matmul(transpose(this%pcat),dlog(c*g))
  catterm=dexp(catterm)
  dcatterm=matmul(transpose(this%pcat),dact)
  do k=1,this%pp%numterm
   dcatterm(k,:)=dcatterm(k,:)*catterm(k) 
  end do
  call check_pointer_ (dact,1,1,.false.)
else
  catterm=1.0d0
  dcatterm=0.0d0 
end if
!%------------------------------------------------------------------
!% Control of the lower limit of thresh(i)
!%------------------------------------------------------------------
if (omtol2>(this%omegathr-1.0d0)) then
 omegathr=1.d0+omtol2
else
 omegathr=this%omegathr
end if
!%------------------------------------------------------------------
if(omega>1.0d0.and.isbe) then
 threshmin = 1.0d0
!%------------------------------------------------------------------
!% Mineral exists+subsat. and Mineral does not exist+supers.
!%------------------------------------------------------------------
else
 threshmin = omegathr
end if
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!% Compute derivatives 
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
!%-------------------------------------------------------------------
do k=1,this%pp%numterm
  rate=this%pp%attrterm(1,k)
  theta=this%pp%attrterm(2,k)
  eta=this%pp%attrterm(3,k)
  if (eta>0.0d0) then
     call omegaterm_ (fdeltag(k),dfdeltag(k),omega,threshmin,theta,eta,omtol1,omtol2,ompsi,iou)
  else if (eta==0.0d0) then
     fdeltag(k)=1.0d0
     dfdeltag(k)=0.0d0
  end if
  drk=drk+rate*(dcatterm(k,:)*fdeltag(k)+catterm(k)*dfdeltag(k)*domega)
end do
!%--------------------------------------------------------------------
!% Multiply for Arrenius equation 
!%-------------------------------------------------------------------
drk=-drk*this%pp%exp_ea_rt
!%--------------------------------------------------------------------
!% If the raction rate depends of the area 
!%-------------------------------------------------------------------
if (this%pp%isareadep) drk=drk*area 
!%-------------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------------------
call check_pointer_ (catterm,1,.false.)
call check_pointer_ (dcatterm,1,1,.false.)
call check_pointer_ (fdeltag,1,.false.)
call check_pointer_ (dfdeltag,1,.false.)
if (iserror) goto 10
!%-------------------------------------------------------------------
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
subroutine copy_rrllasaga &
  (copied, &
   this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy an reaction rate law (lasaga type) in other reaction
! rate law. 
!
!   $Arguments:
!
type(t_reactionratelaw_lasaga), intent(in)  :: this

type(t_reactionratelaw_lasaga), intent(out) :: copied 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer            :: &
 ndim1, &
 ndim2 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
copied%omegathr=this%omegathr
!%-----------------------------------------------------------
if (this%pp%numsp>0) then
 ndim1=size(this%pcat,1)
 ndim2=size(this%pcat,2)
 call check_pointer_ (copied%pcat,ndim1,ndim2,.true.)
 copied%pcat=this%pcat
end if 
!%-----------------------------------------------------------
return
end subroutine
!%***********************************************************
!%***************Public subroutine***************************
!%***********************************************************
!%***********************************************************
!%***********************************************************
subroutine rewrite_rrllasaga &
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
 
type(t_reactionratelaw_lasaga), intent(inout), target  :: this

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
!% 
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
!% Allocate stoichiometric matrix of catalityc species 
!%---------------------------------------------------------------
call check_pointer_ (this%pcat,nsp,this%pp%numterm,.true.)
!%---------------------------------------------------------------
isps=0
do i=1,nsp
 
  do j=1,this%pp%numsp
    name => this%pp%pspecies(j)%ptr%name
    if(name==namesp(i)) then
       isps=isps+1
       this%pcat(i,:) = this%pp%attrsp(j,:)
	   exit  
    end if
    
  end do
  
end do
!%---------------------------------------------------------------
!% Count the number of founded species 
!%---------------------------------------------------------------
if (isps/=this%pp%numsp) then
 msg='Catalytic species defined in the reaction rate law not defined in the species list'
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
subroutine begin_element_handler (name, attributes)
 
 
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes
!%------------------------------------------------------------
 call read_xml_loc_ (name, attributes)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutines**************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine read_xml_loc_rrllasaga &
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
 
character(len=*), intent(in)             :: name

type(dictionary_t), intent(in)           :: attributes

type (t_reactionratelaw_lasaga), optional:: this

logical, optional                        :: iserror

character(len=*), optional               :: msg 
 
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
 omegathr
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
character(len=100)             :: &
 id
character(len=100), pointer,save:: &
 namesp(:,:)
character(len=100), save       :: &
 namerrlaw
character(len=100), pointer    :: &
 typeterm(:)
integer, parameter             :: &
 mxdim=20 
!-------------------------------------------------------------------------
!
!   $code
!

 

!%---------------------------------------------------------------
havethis=present(this)
haveiserror=present(iserror)
 
select case (havethis)
 
case (.true.)
 iserror1=.false. 
 if (haveiserror) iserror=.false.
 call check_pointer_ (typeterm,iexp,.true.)
 typeterm='catalyst'
 call set_ &
  (this%pp, &
   namerrlaw, &
   numsp(1:iexp), &
   attrterm(1:3,1:iexp), &
   attrspterm(:,1:iexp), &
   namesp(1:mxdim,1:iexp), &
   typeterm, &
   iexp, &
   3, &
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
 
 
 call set_(this,this%pp,omegathr,iserror1)
 
 if (iserror1.and.haveiserror) then
   iserror=.true.
   return 
 else if(iserror1.and..not.haveiserror) then
   stop 
 end if
 
 call check_pointer_ (typeterm,1,.false.)
 call check_pointer_ (numsp,1,.false.) 
 call check_pointer_ (attrterm,1,1,.false.)
 call check_pointer_ (attrspterm,1,1,.false.)
 call check_pointer_ (namesp,1,1,.false.)
 
case default
 
 select case (name)
 
 case ('reactionratelaw')
 
  call check_pointer_ (numsp,mxdim,.true.)
  call check_pointer_ (attrterm,mxdim,mxdim,.true.)
  call check_pointer_ (attrspterm,mxdim,mxdim,.true.)
  call check_pointer_ (namesp,mxdim,mxdim,.true.)
  isareadep=.true. 
  iexp=0 
  omegathr=0.0d0 
  id='' 
  call get_value (attributes,"name", id, status)
  namerrlaw=id
  id=''
  call get_value (attributes,"ea", id, status)
  n=0
  call build_data_array (id,reallocal,n)
  ea=reallocal(1)
  id='' 
  call get_value (attributes,"attr1", id, status)
  if (status==0) then 
   n=0
   call build_data_array (id,reallocal,n)
   omegathr=reallocal(1)
  end if 
  id=''
  call get_value (attributes,"areadep", id, status)
  select case (id)
  case ('not','NOT','n','N')
   isareadep=.false.
  end select 
 case ('term')
 
  iexp=iexp+1
  isp=0
  n=0
  call build_data_array (id,reallocal,n)
  attrterm(1,iexp)=reallocal(1)
  id=''
  call get_value (attributes,"attr2", id, status)
  n=0
  call build_data_array (id,reallocal,n)
  attrterm(2,iexp)=reallocal(1)
  id=''
  call get_value (attributes,"attr3", id, status)
  n=0
  call build_data_array (id,reallocal,n)
  attrterm(3,iexp)=reallocal(1)
    
 case ('species')
 
  isp=isp+1
  numsp(iexp)=isp
  id=''
  call get_value (attributes,"name", id, status)
  namesp(isp,iexp)=id
  id='' 
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
end module m_reactionratelaw_lasaga
