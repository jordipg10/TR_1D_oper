module m_parentreactionratelaw
!-------------------------------------------------------------------------
!
!   $Description: Represent the parent reaction rate law
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
create_ &             ! Create reaction rate law object. 
,destroy_ &           ! Destroy reaction rate law object. 
,set_ &               ! Set attributtes in the reaction rate law object. 
,set_pspecies_ &      ! Set pointer to species in the reaction rate law object. 
,update_ &            ! Update parameters that depend of the temperature in the reaction rate object. 
,get_namesp_ &        ! Return the name of species in the reaction rate law
,assignment(=) &      ! Copy the reaction rate law object in other reaction rate law object.
,omegaterm_           ! 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!% Type definition 
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
type, public:: t_parentreactionratelaw
 
character(len=100)                          :: name             ! Name of reaction rate law
 
type (t_pspecies), pointer, dimension(:)    :: pspecies         
 
character(len=100), pointer, dimension(:)   :: typeterm         ! type of terms [nterm]
 
real*8, pointer, dimension(:,:)             :: attrsp           ! general attributes per species [nsp,nterm]

real*8, pointer, dimension(:,:)             :: attrterm         ! general attributes per term [nattrterm,nterm]
 
real*8                                      :: ea               ! Activation energy

real*8                                      :: exp_ea_rt        ! Arrenius equation

real*8                                      :: tempref          ! Temperature of reference

real*8                                      :: pressref         ! Pressure of reference
 
integer                                     :: numsp            ! Number of species

integer                                     :: numterm          ! Number of terms in the reaction rate law 

integer                                     :: nattrterm        ! Number of attributes per term
 
logical, pointer, dimension(:)              :: islocksps        ! islocksps=true, then the species objects were created for the reaction rate law object. 

logical                                     :: isareadep        ! 
 
end type t_parentreactionratelaw
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface create_
 
module procedure create_prrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface destroy_
 
module procedure destroy_prrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_
 
module procedure set_prrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface set_pspecies_
 
module procedure set_pspecies_prrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface update_
 
module procedure update_temp_param_prrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface get_namesp_
 
module procedure get_namesp_prrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface omegaterm_
 
module procedure omegaterm_prrlaw
 
end interface
!%-----------------------------------------------------------------
!%-----------------------------------------------------------------
interface assignment(=)
 
module procedure copy_prrlaw
 
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
subroutine create_prrlaw &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create the parent reaction rate law object
!
!   $Arguments:
!
 
type(t_parentreactionratelaw), intent(inout) :: this 
 
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
this%name = ' '
this%ea = 0.0d0
this%exp_ea_rt = 0.0d0
this%numsp=0
this%numterm=0
this%nattrterm=0
this%tempref=25.0d0
this%pressref=0.0d0
this%isareadep=.true. 
!%------------------------------------------------------------
this%typeterm => null ()
this%attrsp => null ()
this%attrterm => null ()
this%islocksps => null ()
!%------------------------------------------------------------
this%pspecies => null ()
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_prrlaw &
  (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy parent reaction rate law 
!
!   $Arguments:
!
 
type(t_parentreactionratelaw), intent(inout) :: this    ! Type reaction rate law variable 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                      :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
this%name = ' '
!%------------------------------------------------------------
if (this%numsp>0) then
 do i=1,this%numsp
  if (this%islocksps(i)) then
   call destroy_ (this%pspecies(i)%ptr)
   deallocate(this%pspecies(i)%ptr) 
   this%islocksps(i)=.false. 
  end if  
  this%pspecies(i)%ptr => null ()
 end do
 deallocate (this%pspecies)
end if
this%numsp=0
!%------------------------------------------------------------
this%ea = 0.0d0
this%exp_ea_rt = 0.0d0
this%numterm=0
this%nattrterm=0
this%tempref=0.0d0
this%pressref=0.0d0
this%isareadep=.false. 
!%------------------------------------------------------------
!% Deallocate pointers
!%------------------------------------------------------------
call check_pointer_ (this%typeterm,1,.false.)
call check_pointer_ (this%attrsp,1,1,.false.)
call check_pointer_ (this%attrterm,1,1,.false.)
call check_pointer_ (this%islocksps,1,.false.)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_prrlaw &
  (this, &
   name, &
   nspterm, &
   attrterm, &
   attrspterm, &
   namespterm, &
   typeterm, &
   nterm, &
   nattrterm, &
   mxsp, &
   ea, &
   isareadep, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set parent reaction rate law
!
!   $Arguments:
!
 
type(t_parentreactionratelaw), intent(inout)             :: this           ! Type parent reaction rate law object 

integer, intent(in)                                      :: nterm          ! Number of terms 

integer, intent(in)                                      :: nattrterm      ! Number of attributes by term 

integer, intent(in)                                      :: mxsp           ! Maximum number of species 

character(len=*), intent(in)                             :: name           ! Name of the reaction rate law object 

character(len=*), intent(in), dimension(mxsp,nterm)      :: namespterm     ! Name of species by term

character(len=*), intent(in), dimension(nterm)           :: typeterm       ! Type term 

real*8, intent(in)                                       :: ea             ! Activation energy

real*8, intent(in), dimension(mxsp,nterm)                :: attrspterm     ! Attibutes of species

real*8, intent(in), dimension(nattrterm,nterm)           :: attrterm       ! Attributes by term 

integer, intent(in), dimension(nterm)                    :: nspterm        ! Number of species by term [nterm]

logical, intent(in)                                      :: isareadep      ! 

logical, intent(out)                                     :: iserror        ! iserror=true, then there was an error
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer           :: &
 i, &
 j, &
 k, &
 nsp
logical           :: &
 isbe
real*8, pointer           :: &
 array (:,:) => null ()
character(len=100), pointer:: &
 namesp(:) => null ()
integer, parameter:: &
 ndim=100 
character(len=100)         :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
!% Allocate local pointers 
!%----------------------------------------------------------
call check_pointer_ (namesp,mxsp,.true.)
call check_pointer_ (array,ndim,ndim,.true.)
!%----------------------------------------------------------
this%numterm=nterm
this%nattrterm=nattrterm
this%ea=ea
this%name=name
this%isareadep=isareadep 
!%----------------------------------------------------------
if (this%numterm>0) then
 call check_pointer_ (this%typeterm,nterm,.true.)
 this%typeterm=typeterm
end if
!%----------------------------------------------------------
nsp=0
isbe=.false.
do i=1,this%numterm
  do j=1,nspterm(i)
   if(nsp==0) then
    nsp=1
    namesp(nsp) = namespterm (j,i)
    array (nsp,i) = attrspterm (j,i)	
   end if
   do k=1,nsp
    if(namesp(k) == namespterm (j,i)) then
      array (k,i) = attrspterm (j,i)
      isbe=.true.
      exit
    end if
   end do
   if (.not.isbe) then
     nsp=nsp+1
     namesp (nsp) = namespterm (j,i)
     array (nsp,i) = attrspterm (j,i)
   end if
   isbe=.false.
  end do
 end do
!%------------------------------------------------------------
this%numsp=nsp
if (this%numsp>0) then
 call check_pointer_ (this%islocksps,this%numsp,.true.)
 this%islocksps=.true.  
 allocate (this%pspecies(this%numsp))
 do i=1,this%numsp
  allocate (this%pspecies(i)%ptr)
  call create_ (this%pspecies(i)%ptr)
  call set_name_ (this%pspecies(i)%ptr,namesp(i))
 end do
 call check_pointer_ (this%attrsp,this%numsp,this%numterm,.true.)
 this%attrsp=array
end if
!%-------------------------------------------------------------
if (this%numterm>0) then
 call check_pointer_(this%attrterm,this%nattrterm,this%numterm,.true.)
 this%attrterm=attrterm
end if
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (array,1,1,.false.)
!%------------------------------------------------------------
return
10 continue       
print *,'**********************'
print *,'Reaction Rate Law:'
print *,'Name:',this%name
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
subroutine set_pspecies_prrlaw &
  (this, &
   species)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set the pointer to species object 
!
!   $Arguments:
!
 
 
type(t_parentreactionratelaw), intent(inout) :: this       ! Type parent reaction rate law object 

type(t_species), intent(in), target          :: species    ! Type species object 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                      :: &
 isp 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
do isp=1,this%numsp
 if (this%pspecies(isp)%ptr%name==species%name) then
  if (this%islocksps(isp)) then 
   call destroy_ (this%pspecies(isp)%ptr)
   deallocate (this%pspecies(isp)%ptr)
   this%islocksps(isp)=.false.  
  end if 
  this%pspecies(isp)%ptr => species
  exit
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
subroutine update_temp_param_prrlaw &
  (this, &
   temp, &
   iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Solve the arrenius equation
!
!   $Arguments:
!
 
type(t_parentreactionratelaw), intent(inout)  :: this      ! Type reaction rate law object 

real*8, intent(in)                            :: temp      ! Temperature in celcius.

logical, intent(out)                          :: iserror   ! iserror=true, then there was an error. 
 
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
 aux, &
 tempk 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.  
!%------------------------------------------------------------
!% Update the Arrenius equation according temperature 
!%------------------------------------------------------------
tempk=temp+273.15d0
aux=(1.0d0/tempk-1.0d0/298.15d0)/rgas
this%exp_ea_rt = dexp(-this%ea*aux)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_namesp_prrlaw &
  (this, &
   namesp, &
   nsp)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the name of species in the reaction rate law
!
!   $Arguments:
!
 
type(t_parentreactionratelaw), intent(in) :: this

character(len=*), pointer                 :: namesp(:)

integer, intent(out)                      :: nsp 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                      :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
call check_pointer_ (namesp,this%numsp,.true.)
!%------------------------------------------------------------
nsp=this%numsp
if (nsp>0) then
 do i=1,this%numsp
  namesp(i)=this%pspecies(i)%ptr%name
 end do
end if
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_prrlaw &
  (copied, &
   this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy parent reaction rate law object 
!
!   $Arguments:
!
 
type(t_parentreactionratelaw), intent(in) :: this

type(t_parentreactionratelaw), intent(out):: copied 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                                   :: &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
copied%name = this%name
!%------------------------------------------------------------
copied%ea=this%ea
copied%exp_ea_rt=this%exp_ea_rt
copied%tempref=this%tempref
copied%pressref=this%pressref 
copied%isareadep=this%isareadep 
!%------------------------------------------------------------
copied%numsp = this%numsp
copied%numterm=this%numterm
copied%nattrterm=this%nattrterm
!%------------------------------------------------------------
if (this%numsp>0) then
 call check_pointer_(copied%islocksps,copied%numsp,.true.)
 copied%islocksps=this%islocksps
 call check_pointer_ (copied%attrsp,copied%numsp,copied%numterm,.true.)
 copied%attrsp=this%attrsp
 allocate(copied%pspecies(copied%numsp)) 
 do i=1,this%numsp
  if (this%islocksps(i)) then
    allocate(copied%pspecies(i)%ptr)
	call create_ (copied%pspecies(i)%ptr)
	copied%pspecies(i)%ptr=this%pspecies(i)%ptr  
  else
    copied%pspecies(i)%ptr => this%pspecies(i)%ptr   
  end if 
 end do 
end if
!%------------------------------------------------------------
if (this%numterm>0) then
 call check_pointer_(copied%typeterm,copied%numterm,.true.)
 copied%typeterm=this%typeterm
 call check_pointer_(copied%attrterm,copied%nattrterm,copied%numterm,.true.)
 copied%attrterm=this%attrterm
end if
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine omegaterm_prrlaw &
   (OMEGTERM, DOMEGTERM,    OMEGA,    OMTHR,    THETA, &
         ETA,    OMTOL1,   OMTOL2,      PSI,      IOU)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:!%Computes the omega(Sat. index) dependent term (including sign) of the
!%reaction rate equation (after smoothing discontinuities). Also computes
!%its derivative with repect to omega.
!%
!%DESCRIPTION
!%
!%The omega dependent term of the reaction rate equation ends up as follow
!%
!%   if OMEGA < 1-OMTOL1
!%      OMEGTERM = - (1-OMEGA**THETA)**ETA
!%
!%   if 1-OMTOL1 < OMEGA < 0
!%      Cubic interpolation ensuring continuity and
!%      derivability at 1-OMTOL1 and 1.
!%
!%   if 0 < OMEGA < OMTHR-OMTOL2
!%      OMEGTERM = 0
!%
!%   if OMTHR-OMTOL2 < OMEGA < OMTHR+OMTOL2
!%      Cubic interpolation ensuring continuity and
!%      derivability at OMTHR-OMTOL2 and OMTHR+OMTOL2
!%
!%   if OMEGA < OMTHR+OMTOL2
!%      OMEGTERM = (1-OMEGA**THETA)**ETA
!%
!%DEFINITION OF VARIABLES
!%
!%
!
!%
!%     OMEGA........ Saturation index (must be positive)
!%     OMTHR ..... Threshold of omega for precipitation. Except for smooth
!%                   precipitation only occurs when OMEGA > OMTHR (must be
!%                   greter than 1)
!%     THETA ....... Exponent of OMEGA in the reaction rate equation
!%                   (must be positive)
!%     ETA ......... Exponent of (1 - OMEGA**theta) in the reaction rate
!%                   equation (must be positive)
!%     OMTOL1 ...... Width of the smmothing portion at the dissolution sid
!%                   (must be zero or positive less than 1)
!%     OMTOL2 ...... Width of the smmothing portion at the dissolution sid
!%                   (must be zero or positive less than OMTHR-1)
!%     PSI ......... Relaxation factor for precipitation (must be between
!%                   and 1; 0.5 recommended)
!%     IOU ......... Main ouitput unit number
!%
 
!%
!%     OMEGTERM .... Omega dependent term of the reaction rate equation
!%     DOMEGTERM ... Derivative of OMEGTERM with respect to OMEGA.
!%
 
!%
!%     O1 .......... omega at the beginning of the smoothing portion
!%     T1 .......... term at the beginning of the smoothing portion
!%     S1 .......... slope at the beginning of the smoothing portion
!%     XO .......... omega during interpolation
!%     A,B  ........ coefficients of cubic interpolation
!%     IERR ........ error counter
!%
!%
!%    Prepared by Carrera on jan,25,1997
!%    Modified by C. Ayora on April, 16, 1997 and sep,23, 1997
!%************************************************************
!
!   $Arguments:
!
 
REAL*8   :: &
  OMEGA,DOMEGTERM,OMEGTERM,OMTHR,THETA,ETA,OMTOL1,OMTOL2, &
  PSI
INTEGER   :: &
 IOU 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
REAL*8     :: &
  O1,T1,S1,XO,U2,U1,A,B,DUM,OTH
INTEGER    :: &
 IERR 
INTEGER, PARAMETER :: &
 THLOC=0.0D0
!-------------------------------------------------------------------------
!
!   $code
!
 
! EXTERNAL VARIABLES

! LOCAL VARIABLES

 
 
!%.................................................. Initial checks
IERR=0
if ( OMEGA .LT. 0D0 ) then
   WRITE(IOU,9001)
9001    FORMAT(' Error in OMEGTERM_REAC_RATE_:'/ &
          ' Saturation index (OMEGA) must be positive.')
   IERR=IERR+1
end if
 
!%     if ( OMTHR .LE. 1D0 ) then
if ( OMTHR .LT. 1D0 ) then
   WRITE(IOU,9002)
9002    FORMAT(' Error in OMEGTERM_REAC_RATE_:'/ &
          ' Threshold of Saturation index (OMTHR) must be ', &
          'greater or equal than 1.')
   IERR=IERR+1
   STOP
end if
 
if ( (THETA.LT.0D0) .OR. (ETA.LT.0D0) ) then
   WRITE(IOU,9003) THETA,ETA
9003    FORMAT(' Error in OMEGTERM_REAC_RATE_:'/ &
          ' THETA=',E10.3,' ETA=',E10.3,' they must be positive.')
   IERR=IERR+1
   STOP
end if
 
if ( (OMTOL1.LT.0D0) .OR. (OMTOL1.GE.1D0) ) then
   WRITE(IOU,9004) OMTOL1
9004    FORMAT(' Error in OMEGTERM_REAC_RATE:'/ &
          ' OMTOL1=',E10.3, &
          ' It must be zero or positive less than 1.')
   IERR=IERR+1
   STOP
end if
 
if ( OMTOL2.LT.0D0 ) then
   WRITE(IOU,9005) OMTOL2
9005    FORMAT(' Error in OMEGTERM_REAC_RATE:'/ &
          ' OMTOL2=',E10.3, &
          ' It must be zero or positive less than OMTHR-1.')
   IERR=IERR+1
   stop
end if
 
!%Modified C. Ayora
!%Condition for numerical stability during precipitation
if (OMEGA.GT.1.0D0 .AND. OMTOL2.LT.THLOC) OMTOL2=THLOC
!%End
 
!%...................................................... Computes omega te
 
if ( OMEGA.LE.1D0-OMTOL1 ) then          ! Dissolution, no interpola
 
!%C.Ayora: numerical stability (6/23/97)
   OTH=OMEGA**THETA
   dum= 1.0D0-OTH
   if(dum.GE.THLOC) then
     OMEGTERM = -(1D0-OTH)**ETA
     DOMEGTERM = ETA*THETA*(OTH/OMEGA)*((1D0-OTH)**(ETA-1D0))
   else
     OMEGTERM=0.D0
     DOMEGTERM=1.D0
   end if
!%End
 
else if ( OMEGA.GT.1D0-OMTOL1 .AND. OMEGA.LT.1D0) then !Interpola
 
   O1 = 1D0-OMTOL1
   OTH= O1**THETA
   T1 = -(1D0-OTH)**ETA
   S1 = ETA*THETA*(OTH/O1)*((1D0-OTH)**(ETA-1D0))
   O1 = -OMTOL1
   A = (S1*O1 - 2D0*T1)/O1/O1/O1
   B = (3D0*T1 - S1*O1)/O1/O1
   XO= OMEGA - 1D0
   OMEGTERM = (A*XO+B)*XO*XO
   DOMEGTERM = (3D0*A*XO+2D0*B)*XO
 
else if ( OMEGA.GE.1D0 .AND. OMEGA.LE.OMTHR-OMTOL2) then
                                                      ! zero reactio
   OMEGTERM = 0D0
   DOMEGTERM = 0.0D0
 
else if ( OMEGA.GT.OMTHR-OMTOL2 .AND. OMEGA.LT.OMTHR+OMTOL2) then
 
                                                          ! Interpol
   O1 = OMTHR+OMTOL2
   OTH=O1**THETA
   T1 = (OTH-1D0)**ETA
   U1 = OTH/O1
   U2 = T1/(OTH-1D0)
   S1 = ETA*THETA*U1*U2
!%rovi
   dum=1.0D0+PSI*T1
   S1=S1/dum-T1*PSI*S1/dum/dum
   T1=T1/dum
!%rovi
   O1 = 2*OMTOL2
   A = (S1*O1 - 2D0*T1)/O1/O1/O1
   B = (3D0*T1 - S1*O1)/O1/O1
   XO= OMEGA - (OMTHR-OMTOL2)
   OMEGTERM = (A*XO+B)*XO*XO
   DOMEGTERM = (3D0*A*XO+2D0*B)*XO
else                                     ! Precipitation, no inerpol
 
   OTH = OMEGA**THETA
   OMEGTERM = (OTH-1D0)**ETA
   DOMEGTERM = ETA*THETA*(OTH/OMEGA)*(OMEGTERM/(OTH-1D0))
!%Relax precipitation   !C. Ayora 4/15/97)
   dum = 1.0D0+PSI*OMEGTERM
   DOMEGTERM= DOMEGTERM/dum-OMEGTERM*PSI*DOMEGTERM/dum/dum
   OMEGTERM = OMEGTERM/dum
!%Relax precipitation: End
end if
 
return
end subroutine 
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_parentreactionratelaw
