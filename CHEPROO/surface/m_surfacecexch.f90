module m_surface_cexch
!-------------------------------------------------------------------------
!
!   $Description: Represent one interface with cation exchange model
!
!   $Use:
!
!   $Author: Sergio Andrés Bea Jofré
!
!   $License:
!
!-------------------------------------------------------------------------
use m_parentsurface
use m_species 
use flib_xpath
use flib_sax
use m_general_tools_cheproo
use m_constants_cheproo
!%------------------------------------------------------
!%------------------------------------------------------
private                       ::
!%------------------------------------------------------
!%------------------------------------------------------
public                        :: &
create_ &
,destroy_ &
,set_ &
,set_parent_ &
,init_ &
,update_sk_ &
,get_num_sk_ &
,get_num_xna_ &
,get_xna_index_ &
,get_name_xna_ &
,update_ &
,write_ &
,assignment(=) &
,change_to_mol_
!%----------------------------------------------------
!%----------------------------------------------------
integer, parameter             :: &
numsk=1
!%------------------------------------------------------------
!%------------------------------------------------------------
type, public::t_surface_cexch
 
private                             ::
 
type (t_parentsurface), pointer     :: pp        ! Poiter to parent surface
 
real*8, dimension(:,:), pointer     :: cec       ! Cation exchange capacity 
 
logical                             :: lockcec
 
end type t_surface_cexch
!%----------------------------------------------------
!%----------------------------------------------------
interface create_
 
module procedure create_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface destroy_
 
module procedure destroy_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface set_
 
module procedure set_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface set_parent_
 
module procedure set_parent_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface init_
 
module procedure init_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface update_sk_
 
module procedure update_sk_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_fun_
 
module procedure compute_fun_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface compute_dfun_
 
module procedure compute_dfun_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_num_sk_
 
module procedure get_num_sk_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_num_xna_
 
module procedure get_num_xna_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_xna_index_
 
module procedure get_xna_index_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface get_name_xna_
 
module procedure get_name_xna_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface update_
 
module procedure update_xna_in_cd_surfcexch
module procedure update_tempdepparam_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface write_
 
module procedure write_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface assignment(=)
 
module procedure copy_surfcexch
 
end interface
!%----------------------------------------------------
!%----------------------------------------------------
interface change_to_mol_
 
module procedure change_to_mol_surfcexch
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
contains
!%--------------------------------------------------------
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine create_surfcexch &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create surface object with cation exchange model
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(inout)    :: this 
 
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
this%cec => null ()
this%lockcec=.false.
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_surfcexch &
   (this)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy surface object with cation-exchange model
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(inout)    :: this 
 
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
this%lockcec=.false.
!%------------------------------------------------------------
call check_pointer_ (this%cec,1,1,.false.)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_fun_surfcexch &
   (this, &
    fun, &
    beta)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
type(t_surface_cexch), intent(in)   :: this

real*8, intent(out)                 :: fun(this%pp%numsite)

real*8, intent(in)                  :: beta(this%pp%numsp) 
 
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
 j
real*8, parameter  :: &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!

 

!%------------------------------------------------------------
fun=r0
!%------------------------------------------------------------
do i=1,this%pp%numsite
 do j=1,this%pp%numsp
  fun(i)=fun(i)+this%cec(i,j)*beta(j)
 end do
end do
!%------------------------------------------------------------
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dfun_surfcexch &
   (this, &
    dfun, &
    sk, &
    beta)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(in)   :: this

real*8, intent(in)                  :: sk(this%pp%numsite)

real*8, intent(in)                  :: beta(this%pp%numsp)

real*8, intent(out)                 :: dfun(this%pp%numsite,this%pp%numsite) 
 
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
 i, &
 j 
real*8, parameter  :: &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!

!%------------------------------------------------------------
dfun=r0
!%------------------------------------------------------------
do i=1,this%pp%numsite
 do j=1,this%pp%numsp
  dfun(i,i)=dfun(i,i)+this%cec(i,j)*this%cec(i,j)*beta(j)/sk(i)
 end do
end do
!%------------------------------------------------------------
return
 
10 continue 
print *,'surface:'
print *,'Service: compute_dtxoh_'
print *,'Error'
stop
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_num_sk_surfcexch &
   (this, &
    numskloc)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(in)   :: this

integer, intent(out)                :: numskloc 
 
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
 numskloc=numsk*this%pp%numsite
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_num_xna_surfcexch &
   (this, &
    numxna)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(in)   :: this

integer, intent(out)                :: numxna 
 
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
 numxna=this%pp%numsite
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_name_xna_surfcexch &
   (this, &
    namexna, &
    numxna)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
type(t_surface_cexch), intent(in)        :: this

integer, intent(out)                     :: numxna

character(len=*), pointer                :: namexna(:) 
 
 
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
 ipri 
!-------------------------------------------------------------------------
!
!   $code
!


!%------------------------------------------------------------
call check_pointer_ (namexna,this%pp%numsite,.true.)
do i=1,this%pp%numsite
 ipri=this%pp%idxoh(i)
 namexna(i)=this%pp%pspecies(ipri)%ptr%name
end do
numxna=this%pp%numsite
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_xna_index_surfcexch &
   (this, &
    idxna, &
    numxna)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(in)  :: this

integer, intent(out)               :: numxna

integer, pointer                   :: idxna(:) 
 
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
call check_pointer_ (idxna,this%pp%numsite,.true.)
idxna=this%pp%idxoh
numxna=this%pp%numsite
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_xna_in_cd_surfcexch &
   (this, &
    cd, &
    sk)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(in) :: this

real*8, intent(in)                :: sk(this%pp%numsite*numsk)

real*8, intent(out)               :: cd(this%pp%numsp) 
 
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
 i, &
 isk 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
isk=0
do i=1,this%pp%numsite
 cd(this%pp%idxoh(i))=sk(isk+1)
 isk=isk+numsk
end do
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_tempdepparam_surfcexch &
   (this, &
    temp, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Update properties that depends of the temperature
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(inout)  :: this

real*8, intent(in)                    :: temp

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
 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%------------------------------------------------------------
iserror=.false.
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_surfcexch &
   (this, &
    pp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set surface object 
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(inout)      :: this

type(t_parentsurface), intent(in), target :: pp 

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
integer                  :: &
 i, &
 j, &
 isp
character(len=100)       :: &
 nameprop, &
 msg
real*8                   :: &
 z 
!-------------------------------------------------------------------------
!
!   $code
!
!%----------------------------------------------------------
msg=''
iserror=.false.
!%----------------------------------------------------------
this%pp => pp
!%----------------------------------------------------------
if (this%pp%numsite>0.and.this%pp%numsp>0) then
 call check_pointer_ (this%cec,this%pp%numsite,this%pp%numsp,.true.)
!%----------------------------------------------------------
 isp=0
 nameprop='chargecexch'
 do i=1,this%pp%numsite
  do j=1,this%pp%numspsite(i)
   isp=isp+1
   call get_prop_ (this%pp%pspecies(isp)%ptr,z,nameprop,msg,iserror)
   if (iserror) goto 10
   this%cec(i,isp)=z
  end do
 end do
 this%lockcec=.true.
end if
!%----------------------------------------------------------
 
return
 
10 continue      
print *,'********************'
print *,'Surface:'
print *,'Name:', this%pp%name
print *,'Service: set_'
print *, msg
print *,'********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_parent_surfcexch &
   (this, &
    parent)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(inout)         :: this

type(t_parentsurface), intent(in), target    :: parent 
 
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
this%pp => parent
!%----------------------------------------------------------
 
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine init_surfcexch &
   (this, &
    sk, &
    cd, &
    txoh)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Init sk and cd with cec
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(in) :: this

real*8, intent(in)                :: txoh(this%pp%numsite)

real*8, intent(inout)             :: cd(this%pp%numsp)

real*8, pointer                   :: sk(:) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, parameter     :: &
r0=0.0d0, &
r1=1.0d0
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
call check_pointer_ (sk,this%pp%numsite*numsk,.true.)
sk=r1
!%------------------------------------------------------------
call update_ (this,cd,sk)
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine update_sk_surfcexch &
   (this, &
    sk, &
    beta, &
	dbeta, &
    cec, &
    isconvergence, &
    isupmxiter, &
    iter, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:  Compute the xna concentrations in each sites defined
!    in the surface
!
!   $Arguments:
!
 
 
type(t_surface_cexch), intent(in)        :: this

real*8, intent(out)                      :: sk(this%pp%numsite*numsk)

real*8, intent(in)                       :: dbeta(this%pp%numsp,this%pp%numsite*numsk)

real*8                                   :: beta(this%pp%numsp)

real*8, intent(in)                       :: cec(this%pp%numsite)

integer, intent(in)                      :: iter

logical, intent(out)                     :: isconvergence 
 
logical, intent(out)                     :: isupmxiter
 
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
real*8, pointer :: &
 jacobian(:,:) => null (), &
 residual(:) => null (), &
 error(:) => null ()
integer, pointer:: &
 indx(:) => null ()
integer         :: &
 numtotsk, &
 i, &
 j, &
 ipri 
real*8          :: &
 dd, &
 mxerrorunk, &
 mxerrorres, &
 coeff
character(len=100) :: &
 msg 
real*8, parameter  :: &
 r0=0.0d0, &
 r1=1.0d0
!-------------------------------------------------------------------------
!
!   $code
!
!%----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
isconvergence=.false.
isupmxiter=.false.
!%-----------------------------------------------------------
!% Check the number of iterations
!%-----------------------------------------------------------
if (iter>maxiterads) then
 isupmxiter=.true.
 return
end if
!%-----------------------------------------------------------
!% Compute the total number of unknowns
!%-----------------------------------------------------------
numtotsk=this%pp%numsite*numsk
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (jacobian,numtotsk,numtotsk,.true.)
call check_pointer_ (residual,numtotsk,.true.)
call check_pointer_ (indx,numtotsk,.true.)
call check_pointer_ (error,numtotsk,.true.)
!%----------------------------------------------------------
!% Build jacobian and residual
!%----------------------------------------------------------
residual=r1
do i=1,this%pp%numsite
 do j=1,this%pp%numsp
   ipri=this%pp%idxoh(i)
   coeff=this%cec(i,j)
   if (coeff/=r0) then
      residual(i)=residual(i)-beta(j)
      if (j==ipri) then
       jacobian(i,i)=jacobian(i,i)+r1
	  else
	   jacobian(i,i)=jacobian(i,i)+dbeta(j,i)
	  end if 
   end if
 end do
end do
!%-----------------------------------------------------------
!% Compute maximum error in residual
!%-----------------------------------------------------------
mxerrorres=maxval(dabs(residual))
!%-----------------------------------------------------------
!% Solve the lineal system
!%-----------------------------------------------------------
call ludcmp (jacobian,numtotsk,numtotsk,indx,dd,msg,iserror)
if (iserror) goto 20 
call lubksb (jacobian,numtotsk,numtotsk,indx,residual)
!%-----------------------------------------------------------
!% Compute relative error
!%-----------------------------------------------------------
error=residual/sk
error=dabs(error)
mxerrorunk=maxval(error)
!%-----------------------------------------------------------
!% Scale the solution 
!%-----------------------------------------------------------
if (mxerrorunk>facmaxads) then
  residual = residual * facmaxads/mxerrorunk
end if
!%-----------------------------------------------------------
!% Update the solution 
!%-----------------------------------------------------------
sk = sk + residual
!%-----------------------------------------------------------
!% Check convergence 
!%-----------------------------------------------------------
isconvergence=(mxerrorres<=tolresads.and.mxerrorunk<=tolunkads) 
!%----------------------------------------------------------- 
call update_ (this,beta,sk)
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
20 continue 
call check_pointer_ (jacobian,1,1,.false.)
call check_pointer_ (residual,1,.false.)
call check_pointer_ (indx,1,.false.)
call check_pointer_ (error,1,.false.)
if (iserror) goto 10 
!%-----------------------------------------------------------
return
 
10 continue  
print *,'***************************'
print *,'Surface:'
print *,'Name:',this%pp%name
print *,'Service: update_sk_'
print *, msg
print*, '***************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine write_surfcexch &
   (this, &
    ioutput, &
    msg, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_surface_cexch), intent(in) :: this

integer, intent(in)                :: ioutput    

logical, intent(out)               :: iserror

character(len=*), intent(out)      :: msg 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                            :: &
 i, &
 j 
real*8                             :: &
 zexch
!-------------------------------------------------------------------------
!
!   $code
!
 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
write (ioutput,*) '----------------------------------------------'
write (ioutput,*) '           Surface information                '
write (ioutput,*) '         Model: Cation Exchange               '
!%-----------------------------------------------------------
!% Write the parent surface 
!%-----------------------------------------------------------
call write_ (this%pp,ioutput,msg,iserror)
!%-----------------------------------------------------------
!% Write the CEC definition 
!%-----------------------------------------------------------
!%-----------------------------------------------------------
write (ioutput,*) '----------------------------------------------'
write (ioutput,*) '      CEC definition by sites       '
write(ioutput,1) (this%pp%pspecies(j)%ptr%name,j=1,this%pp%numsp)
do i=1,this%pp%numsite
  write(ioutput,2) (this%cec(i,j),j=1,this%pp%numsp)
end do
!%-----------------------------------------------------------
if (iserror) goto 10
!%-----------------------------------------------------------
write (ioutput,*) '----------------------------------------------'
!%-----------------------------------------------------------
return
10 continue  
print *,'***************************'
print *,'Surface:'
print *,'Name:',this%pp%name
print *,'Service: write_'
print *, msg
print*, '***************************'
iserror=.true.
return
1 format (7x,<this%pp%numsp>a10) 
2 format (<this%pp%numsp>f10.1) 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine change_to_mol_surfcexch &
   (this, &
    beta, &
    cec, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Change from equivalent fraction to molality
!
!   $Arguments:
!
 
 
type (t_surface_cexch), intent(in)                :: this    ! Surface cation-exchange object 

real*8, intent(inout), dimension(this%pp%numsp)   :: beta    ! Equivalent fraction/molality of exchange complexes

real*8, intent(in), dimension(this%pp%numsite)    :: cec     ! Cation-exchange capacity 

logical, intent(out)                              :: iserror ! If true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)        :: &
 msg
integer                 :: &
 i, &
 j
real*8                 :: &
 zexch 
real*8, parameter      :: &
 r0=0.0d0
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
do i=1,this%pp%numsite
 do j=1,this%pp%numsp
  zexch=this%cec(i,j)
  if (zexch/=r0) then
   beta(j)=beta(j)*cec(i)/zexch
  end if
 end do
end do
!%-----------------------------------------------------------
return
10 continue 
print *,'***************************'
print *,'Surface:'
print *,'Name:',this%pp%name
print *,'Service: change_to_mol_'
print *, msg
print*, '***************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_surfcexch &
   (targetobj, &
    sourceobj)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Copy surface with cation exchange model
!
!   $Arguments:
!
 
type(t_surface_cexch), intent(in)  :: sourceobj

type(t_surface_cexch), intent(out) :: targetobj 
 
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
 ndim1, &
 ndim2 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
if (sourceobj%lockcec) then
 targetobj%lockcec=sourceobj%lockcec
 ndim1=size(sourceobj%cec,1)
 ndim2=size(sourceobj%cec,2)
 call check_pointer_ (targetobj%cec,ndim1,ndim2,.true.)
 targetobj%cec=sourceobj%cec
end if
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
!%************************************************************
!%************************************************************
!%************************************************************
end module m_surface_cexch
