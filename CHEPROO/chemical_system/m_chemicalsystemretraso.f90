module m_chemicalsystemretraso
!-------------------------------------------------------------------------
!
!   $Description: Represent the chemical system according for reactive transport formulation of Saaltink et. al.(1998)
!
!   $Use: use m_parentchemicalsystem
! use m_reaction
! use m_species
! use m_phase
! use m_surface
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
!%-------------------------------------------------------------------------
!% Modules corresponding to CHEPROO project
!%-------------------------------------------------------------------------
use m_parentchemicalsystem
use m_reaction
use m_species
use m_phase
use m_surface
use m_general_tools_cheproo
use m_constants_cheproo
!%-------------------------------------------------------------------------
!% Modules corresponding to xml parser 
!%-------------------------------------------------------------------------
use flib_xpath
use flib_sax
!%------------------------------------------------------------
!%------------------------------------------------------------
private   ::
!%------------------------------------------------------------
!% Public interfaces 
!%------------------------------------------------------------
public:: &
create_ &                            ! Create a retraso chemical system object. 
,specia_from_cpri_ &       
,specia_from_solution_type_ &
,specia_from_u_ &
,specia_aqueous_phase_ &
,specia_eqmin_from_setre_ &
,switch_base_ &                       ! Switch the chemical base. 
,destroy_ &
,set_ &
,set_parent_ &
,compute_umob_ &
,compute_iumob_ &
,compute_uads_ &
,compute_dcmob_ &
,compute_dcads_ &
,compute_dumob_ &
,compute_duads_ &
,compute_usktrk_ &
,compute_dusktrk_ &
,make_lin_trf_ &
,write_ &
,get_hashcompz_ &
,get_chem_info_ &
,assignment(=)
!%--------------------------------------------------------
!% Private interfaces 
!%--------------------------------------------------------
private :: &
compute_dc1_dc11_ &
,compute_c1_from_c11_ &
,get_elim_index_ &
,index_ &
,add_elimination_matrix_ &
,compute_elimination_matrix_
!%------------------------------------------------------
!% In this private variable, to storage the elimination
!% matrices 
!%------------------------------------------------------
type, private::t_elimination
 
 real*8, pointer, dimension(:,:)              :: elim          ! Elimination matrices
 
 integer, pointer, dimension(:)               :: idrprisp      ! Global index of reduced primary species
 
 integer, pointer, dimension(:)               :: idnrprisp     ! Global index of non-reduced primary species
 
 integer, pointer, dimension(:)               :: ideqminreact  ! Global index of equilibrium minerals
 
 integer                                      :: nrprisp       ! Number of reduced primary species
 
 integer                                      :: nnrprisp      ! Number of non-reduced primary species
 
 integer                                      :: hashindex     ! Hash index of the elimination matrix
 
 integer                                      :: wcindex       ! Index of water component

end type
!%------------------------------------------------------
!%------------------------------------------------------
!% Type definition 
!%------------------------------------------------------
!%------------------------------------------------------
type, public::t_chemicalsystemretraso
 
private                                     ::
 
type(t_parentchemicalsystem), pointer       :: pp         ! Pointer to parent chemical system
 
type(t_elimination), pointer, dimension(:)  :: elim       ! Elimination matrix information
 
integer                                     :: numelim    ! Number of elimination matrixes
 
end type t_chemicalsystemretraso
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!% Public interfaces
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
interface create_
 
module procedure create_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface specia_from_cpri_
 
module procedure specia_from_cpri_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface specia_from_u_
 
module procedure specia_from_u_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface specia_eqmin_from_setre_
 
module procedure specia_eqmin_from_setre_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface specia_aqueous_phase_
 
module procedure specia_aqueous_phase_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface destroy_
 
module procedure destroy_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface switch_base_
 
module procedure switch_base_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_
 
module procedure set_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface set_parent_
 
module procedure set_parent_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface specia_from_solution_type_
 
module procedure specia_from_solution_type_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_umob_
 
module procedure compute_umob_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_iumob_
 
module procedure compute_iumob_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_uads_
 
module procedure compute_uads_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dcmob_
 
module procedure compute_dcmob1_chemsysretraso
module procedure compute_dcmob2_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dcads_
 
module procedure compute_dcads_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dumob_
 
module procedure compute_dumob1_chemsysretraso
module procedure compute_dumob2_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_duads_
 
module procedure compute_duads_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_usktrk_
 
module procedure compute_usktrk1_chemsysretraso
module procedure compute_usktrk2_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dusktrk_
 
module procedure compute_dusktrk1_chemsysretraso
module procedure compute_dusktrk2_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface make_lin_trf_
 
module procedure make_lin_trf_vector_chemsysretraso
module procedure make_lin_trf_array_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface write_
 
module procedure write_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_chem_info_
 
module procedure get_chem_info_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_hashcompz_
 
module procedure get_hashcompz_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface assignment (=)
 
module procedure copy_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
!% Private interfaces 
!%--------------------------------------------------------
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_dc1_dc11_
 
module procedure compute_dc1_dc11_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_c1_from_c11_
 
module procedure compute_c1_from_c11_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface get_elim_index_
 
module procedure get_elim_index_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface index_
 
module procedure index_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface add_elimination_matrix_
 
module procedure add_elimination_matrix_chemsysretraso
 
end interface
!%--------------------------------------------------------
!%--------------------------------------------------------
interface compute_elimination_matrix_
 
module procedure compute_elimination_matrix_chemsysretraso
 
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
subroutine create_chemsysretraso &
   (this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Create chemical system according Saaltink et al. (1998)
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(inout) :: this   ! Type retraso chemical system. 

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
this%elim => null ()
this%numelim=0
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine destroy_chemsysretraso &
   (this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Destroy chemical system RETRASO object
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(inout):: this 
 
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
 i 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
this%pp => null ()
!%------------------------------------------------------------
if (this%numelim>0) then
 do i=1,this%numelim
  if (this%elim(i)%nnrprisp>0) then
   call check_pointer_ (this%elim(i)%ideqminreact,1,.false.)
   call check_pointer_ (this%elim(i)%idnrprisp,1,.false.)
  end if
  call check_pointer_ (this%elim(i)%idrprisp,1,.false.)
  call check_pointer_ (this%elim(i)%elim,1,1,.false.)
 end do
 deallocate (this%elim)
 this%elim => null ()
end if
!%------------------------------------------------------------
this%numelim=0
!%------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_chemsysretraso &
   (this, &
    pp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Set chemical system type RETRASO
!%    ---------------
!%    pp => Parent chemical system pointer
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso) , intent(inout)     :: this

type (t_parentchemicalsystem), intent(in), target  :: pp

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
integer                              :: &
 i
real*8, pointer                      :: &
 colu(:,:) => null ()
character(len=100)                   :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------
!% Set the parent chemical system
!%---------------------------------------------------------
this%pp => pp
!%---------------------------------------------------------
!% Compute the first set of elimination matrix not 
!% corresponding to any set of equilibrium minerals
!%---------------------------------------------------------
this%numelim=1
allocate(this%elim(1))
this%elim(1)%elim => null ()
this%elim(1)%idrprisp => null ()
this%elim(1)%idnrprisp => null ()
this%elim(1)%ideqminreact => null ()
this%elim(1)%hashindex=0
this%elim(1)%nrprisp=this%pp%numaqprisp
this%elim(1)%nnrprisp=0
call check_pointer_(this%elim(1)%elim,this%pp%numaqprisp,this%pp%numaqprisp,.true.)
call check_pointer_(this%elim(1)%idrprisp,this%pp%numaqprisp,.true.)
!%---------------------------------------------------------
do i=1,this%pp%numaqprisp
 this%elim(1)%idrprisp(i)=i
end do
!%---------------------------------------------------------
!% Build the first elimination matrix (identity matrix)
!%---------------------------------------------------------
call compute_elimination_matrix_ &
   (this, &
    colu, &
    this%elim(1)%elim, &
    this%elim(1)%nnrprisp, &
    this%elim(1)%nrprisp, &
    msg, &
    iserror)
if (iserror) goto 10
!%------------------------------------------------------------
if (this%pp%wcindex>0) then
 this%elim(1)%wcindex=this%pp%wcindex
else
 this%elim(1)%wcindex=0
end if
!%------------------------------------------------------------
return
 
10 continue 
print *,'*****************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: set_'
print *, msg
 print *,'*****************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine set_parent_chemsysretraso &
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
 
type (t_chemicalsystemretraso) , intent(inout)     :: this

type (t_parentchemicalsystem), intent(in), target  :: pp

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
character(len=100)        :: &
msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
msg=''
!%-----------------------------------------------------------
!% Set the parent chemical system
!%-----------------------------------------------------------
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
!%****************Public subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine specia_from_solution_type_chemsysretraso &
   (this, &
    c, &
    g, &
    alpha, &
    numsp, &
    nummin, &
    simin, &
    namemin, &
    numcomp, &
    namecomp, &
    icon, &
    ctot, &
    cguess, &
    constraint, &
    ionstr, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    numtxoh, &
    numsurf, &
    isconvergence, &
    hashcompz, &
    iserror, &
    dc, &
	dg, & 
    sktrk, &
    dsktrk, &
	cputime, &
    nchemiter)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Specia according solution type. 
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(inout)    :: this

integer, intent(in)                :: numcomp

integer, intent(in)                :: numtxoh

integer, intent(in)                :: numsurf

integer, intent(in)                :: numsp

integer, intent(out)               :: nummin

integer, intent(out)               :: hashcompz

integer, intent(in)                :: icon(numcomp)

real*8, intent(out)                :: c(numsp)

real*8, intent(out)                :: g(numsp)

real*8, intent(in)                 :: alpha(numsp)

real*8, intent(in)                 :: cguess(numcomp)

real*8, intent(in)                 :: ctot(numcomp)

real*8, intent(in)                 :: txoh(numtxoh,numsurf)

real*8, intent(in)                 :: cap1(numtxoh,numsurf)

real*8, intent(in)                 :: cap2(numtxoh,numsurf)

real*8, intent(in)                 :: spsurfarea(numtxoh,numsurf)

real*8, pointer                    :: simin(:)

character(len=*), intent(in)       :: constraint(numcomp)

character(len=*), intent(in)       :: namecomp(numcomp)

character(len=100), pointer        :: namemin(:)

real*8, intent(out)                :: ionstr

logical, intent(out)               :: iserror

logical, intent(out)               :: isconvergence

real*8, pointer, optional          :: dc(:,:)

real*8, pointer, optional          :: dg(:,:)

real*8, pointer, optional          :: sktrk(:)

real*8, pointer, optional          :: dsktrk(:,:) 

integer, intent(out), optional     :: nchemiter

real*8, intent(out), optional      :: cputime               ! CPU time consumed during speciation 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)          :: &
 msg
logical                     :: &
 havedc, &
 havedg, &
 havedsktrk, &
 isbe
integer                     :: &
 ielim, &
 nrpri
real*8, pointer             :: &
 c1(:) => null (), &
 g1(:) => null (), &
 dc1(:,:) => null (), &
 array(:,:) => null (), &
 dgloc(:,:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%-----------------------------------------------------------
msg=''
iserror=.false.
!%-----------------------------------------------------------
!% Check optional arguments 
!%-----------------------------------------------------------
havedc=present(dc)
havedg=present(dg)
havedsktrk=present(dsktrk)
!%-----------------------------------------------------------
!% Specia from solution type 
!%-----------------------------------------------------------
call specia_from_solution_type_ &
   (this%pp, &
    c, &
    g, &
    alpha, &
    numsp, &
    nummin, &
    simin, &
    namemin, &
    numcomp, &
    namecomp, &
    icon, &
    ctot, &
    cguess, &
    constraint, &
    ionstr, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    numtxoh, &
    numsurf, &
    isconvergence, &
    iserror, &
    dc=dc, &
    dg=dgloc, &
    sktrk=sktrk, &
    dsktrk=dsktrk, &
	cputime=cputime, &
    nchemiter=nchemiter)

if (iserror.or..not.isconvergence) goto 20 
!%-----------------------------------------------------------
! Compute index of the components zone
!%-----------------------------------------------------------
call get_hashcompz_ (this,hashcompz,c,g,numsp,iserror)
call get_elim_index_ (this,ielim,hashcompz,isbe)
if (iserror) then
 msg='Error when calling get_hashcompz_'
 goto 20
end if
nrpri=this%elim(ielim)%nrprisp
!%-----------------------------------------------------------
! If there is derivates as arguments
!%-----------------------------------------------------------
if (havedc.or.havedg.or.havedsktrk) then
  call check_pointer_ (c1,this%pp%numaqprisp,.true.)
  call check_pointer_ (g1,this%pp%numaqprisp,.true.)
  call check_pointer_ (dc1,this%pp%numaqprisp,this%elim(ielim)%nrprisp,.true.)
  c1=c(this%pp%idaqprisp)
  g1=g(this%pp%idaqprisp)
  call compute_dc1_dc11_ &
    (this, &
     dc1, &
     c1, &
     g1, &
     dgloc(this%pp%idaqprisp,:), &
     this%elim(ielim)%nnrprisp, &
     this%elim(ielim)%nrprisp, &
     this%elim(ielim)%idrprisp, &
     this%elim(ielim)%idnrprisp, &
     this%elim(ielim)%ideqminreact, &
	 this%elim(ielim)%wcindex)
end if
!%-----------------------------------------------------------
!% Optional 
!%-----------------------------------------------------------
if (havedc) then
 call check_pointer_ (array,numsp,this%pp%numaqprisp,.true.)
 array=dc
 call check_pointer_ (dc,numsp,nrpri,.true.)
 dc=matmul(array,dc1)
end if
!%-----------------------------------------------------------
!% Optional 
!%-----------------------------------------------------------
if (havedg) then
 call check_pointer_ (dg,numsp,nrpri,.true.)
 dg=matmul(dgloc,dc1)
end if
!%-----------------------------------------------------------
!% Optional 
!%-----------------------------------------------------------
if (havedsktrk) then
 call check_pointer_ (array,numsp,this%pp%numaqprisp,.true.)
 array=dsktrk
 call check_pointer_ (dsktrk,numsp,nrpri,.true.)
 dsktrk=matmul(array,dc1)
end if
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
20 continue 
call check_pointer_ (array,1,1,.false.)
call check_pointer_ (c1,1,.false.)
call check_pointer_ (dc1,1,1,.false.)
call check_pointer_ (g1,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
if (iserror) goto 10 
!%-----------------------------------------------------------
return
 
10 continue 
print *,'********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: specia_'
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
subroutine specia_from_u_chemsysretraso &
   (this, &
    temp, &
    c, &
    g, &
	iseqmin, &
	iseqgas, &
    ionstr, &
    alpha, &
    t, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    nsp, &
    npri, &
    ntxoh, &
    nsurf, &
    dtime, &
    hashcompz, &
    isconvergence, &
    factoromgw, &
	volgas, &
	pgas, &
	msg, &
    iserror, &
    cguess, &
    dc, &
    sktrk, &
    dsktrk, &
    dg, &
    simin, &
    nchemiter, &
    thetakin)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Make the chemical speciation from analytic concentrations
!% 
!%    ----------------------
!%    this=                Chemical system variable
!%    c[nsp]=            Concentrations vector in mol/Kgw
!%    g[nsp]=            Activity coefficients vector
!%    ionstr=                Ionic Strength
!%    alpha[nsp]=         Updated area of minerals in m2min/kgw
!%    alpha0[nsp] =       Initial area of minerals in m2min/kgw
!%    t[npri]=             Analytic total concentrations in mol/kgw
!%    cguess[npri]=        Initial concentrations of primary species for
!%                         start iterative precess
!%    txoh[ntxoh,nsurf]=   Total of sites for adsorption mol/kgw
!%
!%    nsp=               Number of species
!%    npri=                Number of aqueous primary species
!%    ntxoh=               Number of total sites for adsorption
!%    nsurf=               Number of surfaces
!%    dt=                  Time increment in seconds
!%    hashcompz=          Components zone index
!%    isnonconvergence=      .true. if there is converegence problems
!%    iserror=             .true. if there is error
!%
 
!%    ------------------
!%    temp=       Temperature in celsius
!%    dc_ext=     Derivate of the concentrations with respect to
!%                primary species
!%    sktrk_ext=  Variation of the concentrations due kinetic
!%                reactions
!%    dsktrk_ext= Derivate of the variation of the concentrations
!%                for kinetic reactions with respect to
!%                primary species
!%    dg_ext=     Derivate of the activity coefficients with respect to
!%                primary species
!%    ioutput=    Unit for to write newton raphson information
!%    mh2o
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(inout) :: this

real*8, intent(in)                            :: temp             ! Temperature

real*8, intent(in)                            :: dtime

real*8, intent(out)                           :: ionstr

real*8, intent(out)                           :: factoromgw

real*8, intent(in)                            :: volgas           ! Volume og gas or gas preassure

real*8, intent(in)                            :: pgas             ! Gas preassure

logical, intent(in)                           :: iseqgas          ! If true, then gas constraint is the gas preassure

integer, intent(in)                           :: nsp              ! Number of species 

integer, intent(in)                           :: ntxoh

integer, intent(in)                           :: nsurf

integer, intent(in)                           :: npri

integer, intent(inout)                        :: hashcompz

real*8, intent(inout), dimension(nsp)         :: c 

real*8, intent(inout), dimension(nsp)         :: g

real*8, intent(in), dimension(nsp)            :: alpha

real*8, intent(in), dimension(npri)           :: t

real*8, intent(in), dimension(ntxoh,nsurf)    :: txoh

real*8, intent(in), dimension(ntxoh,nsurf)    :: cap1

real*8, intent(in), dimension(ntxoh,nsurf)    :: cap2

real*8, intent(in), dimension(ntxoh,nsurf)    :: spsurfarea

logical, intent(out)                          :: iserror

character(len=*), intent(out)                 :: msg 

logical, intent(out)                          :: isconvergence

logical, intent(in)                           :: iseqmin

real*8, pointer, optional, dimension(:,:)     :: dc

real*8, pointer, optional, dimension(:)       :: sktrk

real*8, pointer, optional, dimension(:,:)     :: dsktrk

real*8, pointer, optional, dimension(:,:)     :: dg

real*8, pointer, optional, dimension(:)       :: simin

integer, intent(out), optional                :: nchemiter

real*8, intent(in), optional, dimension(npri) :: cguess 

real*8, intent(in), optional                  :: thetakin 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical                     :: &
 havedc, &
 havedsktrk, &
 havedg, &
 be
real*8, pointer             :: &
 dgloc(:,:) => null (), &
 dc1(:,:) => null (), &
 array(:,:) => null (), &
 c1(:) => null (), &
 g1(:) => null ()
integer                     :: &
 ielim 
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
havedc=present(dc)
havedsktrk=present(dsktrk)
havedg=present(dg)
!%---------------------------------------------------------------
call specia_from_u_ &
   (this%pp, &
    temp, &
    c, &
    g, &
	iseqmin, &
	iseqgas, &
    ionstr, &
    alpha, &
    t, &
    txoh, &
    cap1, &
    cap2, &
    spsurfarea, &
    nsp, &
    npri, &
    ntxoh, &
    nsurf, &
    dtime, &
    isconvergence, &
    factoromgw, &
	volgas, &
	pgas, &
	msg, &
    iserror, &
    cguess=cguess, &
    dc=dc, &
    sktrk=sktrk, &
    dsktrk=dsktrk, &
    dg=dgloc, &
    simin=simin, &
    nchemiter=nchemiter, &
    thetakin=thetakin)
 
if (iserror.or..not.isconvergence) goto 20
!%---------------------------------------------------------------
!%---------------------------------------------------------------
call get_hashcompz_ (this,hashcompz,c,g,this%pp%numsp,iserror)
call get_elim_index_ (this,ielim,hashcompz,be)
if (iserror) then
      msg='Error when calling get_hashcompz_'
      goto 20
end if
!%---------------------------------------------------------------
! If hashcompz equal 0
!%---------------------------------------------------------------
if (hashcompz==0) then
 if (havedg) then
  call check_pointer_ (dg,size(dgloc,1),size(dgloc,2),.true.)
  dg=dgloc
 end if
!%---------------------------------------------------------------
! If hashcompz not equal 0
!%---------------------------------------------------------------
else
 
 if (havedc.or.havedsktrk.or.havedg) then
 
  call check_pointer_ (c1,this%pp%numaqprisp,.true.)
  call check_pointer_ (g1,this%pp%numaqprisp,.true.)
  call check_pointer_ (dc1,this%pp%numaqprisp,this%elim(ielim)%nrprisp,.true.)
  c1=c(this%pp%idaqprisp)
  g1=g(this%pp%idaqprisp)
  call compute_dc1_dc11_ &
    (this, &
     dc1, &
     c1, &
     g1, &
     dgloc(this%pp%idaqprisp,:), &
     this%elim(ielim)%nnrprisp, &
     this%elim(ielim)%nrprisp, &
     this%elim(ielim)%idrprisp, &
     this%elim(ielim)%idnrprisp, &
     this%elim(ielim)%ideqminreact, &
	 this%elim(ielim)%wcindex)
 
 
 end if
!%---------------------------------------------------------------
!%--------------------------------------------------------------- 
 call check_pointer_ (array,this%pp%numsp,this%pp%numaqprisp,.true.)
!%---------------------------------------------------------------
!%---------------------------------------------------------------
 if (havedc) then
  array=dc
  call check_pointer_ (dc,this%pp%numsp,this%elim(ielim)%nrprisp,.true.)
  dc=matmul(array,dc1)
 end if
!%---------------------------------------------------------------
!% Compute dg 
!%---------------------------------------------------------------
 if (havedg) then
  array=dgloc
  call check_pointer_ (dg,this%pp%numsp,this%elim(ielim)%nrprisp,.true.)
  dg=matmul(array,dc1)
 end if
!%---------------------------------------------------------------
!% Compute dsktrk 
!%---------------------------------------------------------------
 if (havedsktrk) then
  array=dsktrk
  call check_pointer_ (dsktrk,this%pp%numsp,this%elim(ielim)%nrprisp,.true.)
  dsktrk=matmul(array,dc1)
 end if
 
end if
!%---------------------------------------------------------------
20 continue 
!%---------------------------------------------------------------
!% Deallocate local pointers 
!%---------------------------------------------------------------
call check_pointer_ (dgloc,1,1,.false.)
call check_pointer_ (array,1,1,.false.)
call check_pointer_ (dc1,1,1,.false.)
call check_pointer_ (c1,1,.false.)
call check_pointer_ (g1,1,.false.)
if (iserror) goto 10
!%---------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: specia_from_u_'
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
subroutine specia_eqmin_from_setre_chemsysretraso &
   (this, &
    ck1,  & 
    g,    & 
    ck,   & 
    iscompzchanged, &
    setre, &
    nsp, &
    dtime, &
    hashcompz, &
    omgwfreek1, &
    omgwfreek, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(inout):: this

real*8, intent(in)                           :: dtime          ! Time increment 

integer, intent(in)                          :: nsp            ! Number of species 

integer, intent(inout)                       :: hashcompz

real*8, intent(inout), dimension(nsp)        :: ck1            ! Vector of concentrations in  

real*8, intent(in), dimension(nsp)           :: ck             ! concentration in k.       

real*8, intent(inout), dimension(nsp)        :: setre          ! Changes due equilibrium reactions [mol/s] 

real*8, intent(in), dimension(nsp)           :: g              ! Activity coefficients vector

real*8, intent(in)                           :: omgwfreek1     ! Mass of free water in k+1. 

real*8, intent(in)                           :: omgwfreek      ! Mass of free water in k. 

logical, intent(out)                         :: iserror        ! iserror=true, there was an error. 

logical, intent(out)                         :: iscompzchanged 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)          :: &
 msg
character(len=100), pointer :: &
 namesp(:) => null (), &
 namesp1(:) => null ()
integer                     :: &
 nreact, &
 ielim, &
 i, &
 j, &
 isps, &
 ireact, &
 ireact1, &
 isps1, &
 isps2, &
 hashcompz1
logical                     :: &
 isbe
real*8                      :: &
 value 
real*8, pointer             :: &
 b(:) => null (), &
 stqi(:) => null ()
integer, pointer            :: &
 idreact(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!
 

!%---------------------------------------------------------------
iserror=.false.
msg=''
!%---------------------------------------------------------------
iscompzchanged=.false.
!%---------------------------------------------------------------
!% Check the number of species 
!%---------------------------------------------------------------
if (nsp/=this%pp%numsp) then
 msg='Error in number of species' 
 goto 10 
end if 
!%---------------------------------------------------------------
!% Allocate local pointers 
!%---------------------------------------------------------------
call check_pointer_ (idreact,this%pp%numreact,.true.)
call check_pointer_ (b,nsp,.true.)
call check_pointer_ (namesp1,nsp,.true.)
call get_chem_info_ (this%pp,msg,iserror,namesp=namesp)
if (iserror) goto 20
!%---------------------------------------------------------------
!% Check if the components zone was defined previously 
!%---------------------------------------------------------------
call get_elim_index_ (this,ielim,hashcompz,isbe)
if (.not.isbe) then
 msg='Error, not defined components zone'
 call add_ (msg,ielim)
 goto 20
end if
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!% Choose the equilibrium reactions 
!%---------------------------------------------------------------
!%---------------------------------------------------------------
!% Determine and storage indices of equilibrium reactions 
!%---------------------------------------------------------------
nreact=0 
do i=1,this%pp%numreact
 isps=this%pp%idreactsp(i)
 if (.not.this%pp%iskinreact(i).and.ck(isps)>0.0d0) then
  nreact=nreact+1
  idreact(nreact)=i
 end if
end do
!%---------------------------------------------------------------
!% Add equilibrium reactions corresponding of mineral zone
!%---------------------------------------------------------------
isbe=.false.
do i=1,this%elim(ielim)%nnrprisp
 ireact=this%elim(ielim)%ideqminreact(i)
 do j=1,nreact
  if (idreact(j)==ireact) then
   isbe=.true.
   exit
  end if
 end do
 
 if (isbe) then
  isbe=.false.
 else
  nreact=nreact+1
  idreact(nreact)=ireact
 end if
end do
!%---------------------------------------------------------------
!% If there aren't equilibrium reactions return
!%---------------------------------------------------------------
if (nreact>0) then 
!%---------------------------------------------------------------
!% Select the changes in the aqueous species 
!%---------------------------------------------------------------
   call get_iposspsph_ (this%pp,this%pp%aqphindex,isps1,isps2)
   isps=isps2-isps1+1
   b(1:isps)=setre(isps1:isps2)
   namesp1(1:isps)=namesp(isps1:isps2)
!%---------------------------------------------------------------
!% Select the equations corresponding to sorption complexes
!%---------------------------------------------------------------
   do i=1,this%pp%numsurf 
      call get_iposspsurf_ (this%pp,i,isps1,isps2)
      do j=isps1,isps2
	   isps=isps+1
       b(isps)=setre(j)
       namesp1(isps)=namesp(j)
	  end do
   end do
!%---------------------------------------------------------------
!% Not consider changes in the water species !!!!!
!%---------------------------------------------------------------
 where (namesp1=='h2o')
   namesp1=''
   b=0.0d0
 end where
!%----------------------------------------------------------------
!% Compute re from Set re
!%----------------------------------------------------------------
 call compute_r_from_stqtr_ &
   (this%pp, &
    b, &
    namesp1, &
    isps, &
    idreact, &
    nreact, &
    msg, &
    iserror)
 if (iserror) goto 20
!%---------------------------------------------------------------
!% Update the equilibrium mineral concentrations
!% Warning some concentrations can be negatives 
!%---------------------------------------------------------------
 setre=0.0d0 
 do i=1,nreact
    ireact=idreact(i)
    stqi => this%pp%stq(ireact,:)
    do j=1,this%elim(ielim)%nnrprisp
       ireact1=this%elim(ielim)%ideqminreact(j)
       if (ireact1==ireact) then
          isps=this%pp%idreactsp(ireact)
          ck1(isps)=ck(isps)*omgwfreek+stqi(isps)*b(i)*dtime
          ck1(isps)=ck1(isps)/omgwfreek1
       end if
    end do
    setre=setre+stqi*b(i)
 end do
 setre=setre/omgwfreek1
end if 
!%---------------------------------------------------------------
!% Compute the component zone index
!%---------------------------------------------------------------
hashcompz1=hashcompz
call get_hashcompz_ (this,hashcompz,ck1,g,nsp,iserror)
if (iserror) goto 20
!%---------------------------------------------------------------
!% Check if the components zone was changed
!%---------------------------------------------------------------
iscompzchanged=(hashcompz/=hashcompz1)
!%---------------------------------------------------------------
!% Zeroing negative concentrations 
!% CUIDADO SACAR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! cprovi 
!%---------------------------------------------------------------
where(ck1<0.0d0)
 ck1=0.0d0
end where 
!%---------------------------------------------------------------
20 continue 
!%---------------------------------------------------------------
!% Deallocate local pointers 
!%---------------------------------------------------------------
call check_pointer_ (idreact,1,.false.)
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (namesp1,1,.false.)
call check_pointer_ (b,1,.false.)
stqi => null ()
if (iserror) goto 10
!%---------------------------------------------------------------
return
 
10 continue 
print *,'*********************************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: specia_eqmin_from_setre_'
print *, msg
print *,'*********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_umob_chemsysretraso &
   (this, &
    umob, &
    npri, &
    nmobph, &
    nsp, &
    c, &
    hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in) :: this

integer, intent(out)                       :: nmobph

integer, intent(out)                       :: npri

integer, intent(in)                        :: nsp

integer, intent(in)                        :: hashcompz

real*8, intent(in), dimension(nsp)         :: c 

real*8, pointer, dimension(:,:)            :: umob

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
integer                               :: &
 ielim
logical                               :: &
 isbe
character(len=300)                    :: &
 msg
real*8, pointer                       :: &
 mloc(:,:) 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
call get_elim_index_ (this, ielim,hashcompz,isbe)
if (.not.isbe) then
 msg='Error, not defined components zone'
 call add_ (msg, hashcompz)
 goto 10
end if
!%-------------------------------------------------------------
call compute_umob_ (this%pp,umob,npri,nmobph,nsp,c,iserror)
if (iserror) goto 10 
!%-------------------------------------------------------------
call check_pointer_ (mloc,this%pp%numaqprisp,nmobph,.true.)
 
mloc=umob
 
npri=this%elim(ielim)%nrprisp
 
call check_pointer_ (umob,npri,nmobph,.true.)
 
umob=matmul(this%elim(ielim)%elim,mloc)
!%-------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------- 
call check_pointer_ (mloc,1,1,.false.)
!%-------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: compute_umob_'
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
subroutine compute_uads_chemsysretraso &
   (this, &
    uads, &
    npri, &
    c, &
    nsp, &
    hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
 
type (t_chemicalsystemretraso), intent(in)       :: this

integer, intent(in)                              :: nsp

integer, intent(in)                              :: hashcompz

integer, intent(out)                             :: npri

real*8, intent(in), dimension(nsp)               :: c 

real*8, pointer, dimension(:)                    :: uads 

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
integer                               :: &
 ielim
character(len=100)                    :: &
 msg
real*8, pointer                       :: &
 vloc(:)
logical                               :: &
 isbe 
!-------------------------------------------------------------------------
!
!   $code
!
 

 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call get_elim_index_ (this, ielim, hashcompz,isbe)
if (.not.isbe) then
 msg='Error, not defined components zone'
 call add_ (msg, hashcompz)
 goto 10
end if
!%------------------------------------------------------------
call compute_uads_(this%pp,uads,npri,c,nsp,iserror)
 
if (iserror) return
 
npri=this%elim(ielim)%nrprisp
 
call check_pointer_ (vloc,this%pp%numaqprisp,.true.)
 
vloc=uads
 
call check_pointer_ (uads,npri,.true.)
 
uads=matmul(this%elim(ielim)%elim,vloc)
 
call check_pointer_ (vloc,1,.false.)
 
!%------------------------------------------------------------
return
 
10 continue 
print *,'**********************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: compute_umob_'
print *,msg
print *,'**********************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dcmob1_chemsysretraso &
   (this, &
    dcmob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    hashcompz, &
    iserror, &
    dg, &
    g)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute dcmob/dc1
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in)       :: this

integer, intent(in)                              :: numsp

integer, intent(in)                              :: hashcompz

real*8, intent(in), dimension(numsp)             :: c

real*8, pointer, dimension(:,:)                  :: dcmob

integer, intent(out)                             :: nrow

integer, intent(out)                             :: ncol

integer, intent(out)                             :: nmobph

logical, intent(out)                             :: iserror

real*8, pointer, dimension(:,:), optional        :: dg

real*8, pointer, dimension(:), optional          :: g 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 
 

 

!%------------------------------------------------------------
msg=''
iserror=.false.
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
 
10 continue 
print *,'**************************'
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
subroutine compute_dcmob2_chemsysretraso &
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
!   $Description:Compute dcmob/dc1
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso)       :: &
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
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
 
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 
 


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
 
10 continue 
print *,'**************************'
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
subroutine compute_dumob1_chemsysretraso &
   (this, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    c, &
    numsp, &
    hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute the derivate of the concentration with respect 
! to primary species of $ith$ nodal chemistry, evaluated in the components 
! matrix corresponding to $jth$ nodal chemistry.
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso)       :: &
 this
integer                               :: &
 numsp, &
 hashcompz
real*8                                :: &
 c(numsp)
real*8, pointer                       :: &
 dumob(:,:)
integer, intent(out)                  :: &
 nrow, &
 ncol, &
 nmobph
logical                               :: &
 iserror 
 
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
 ipos1, &
 ipos2, &
 ielim, &
 ipri1, &
 ipri2
real*8, pointer                       :: &
 dcmob(:,:) => null(), &
 dg(:,:) => null(), &
 g(:) => null(), &
 dc1(:,:) => null(), &
 mloc(:,:) => null()
logical                               :: &
 be
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 
 

 

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
    c, &
    numsp, &
    iserror, &
    dg=dg, &
    g=g)
 
if (iserror) then
 msg='Error when calling compute_dcmob_)'
 goto 20
end if
!%---------------------------------------------------------------
call get_elim_index_ (this, ielim,hashcompz, be)
if (.not.be) then
 msg='Error, not defined components zone'
 call add_ (msg,hashcompz)
 goto 20
end if
!%---------------------------------------------------------------
call check_pointer_ (mloc,nrow*nmobph,ncol,.true.)
mloc=dumob
!%---------------------------------------------------------------
nrow=this%elim(ielim)%nrprisp
ncol=this%elim(ielim)%nrprisp
call check_pointer_ (dumob,nrow*nmobph,ncol,.true.)
call check_pointer_ (mloc,this%pp%numaqprisp,nrow,.true.)
call check_pointer_ (dc1,this%pp%numaqprisp,nrow,.true.)
!%--------------------------------------Get aqueous species index
call compute_dc1_dc11_ &
    (this, &
     dc1, &
     c(this%pp%idaqprisp), &
     g(this%pp%idaqprisp), &
     dg(this%pp%idaqprisp,:), &
     this%elim(ielim)%nnrprisp, &
     this%elim(ielim)%nrprisp, &
     this%elim(ielim)%idrprisp, &
     this%elim(ielim)%idnrprisp, &
     this%elim(ielim)%ideqminreact, &
	 this%elim(ielim)%wcindex)
!%--------------------------------------------------------------
ipri1=1
ipri2=nrow
ipos1=1
ipos2=this%pp%numsp
mloc=matmul(this%pp%ueq,matmul(dcmob(ipos1:ipos2,:),dc1))
dumob(ipri1:ipri2,:)=matmul(this%elim(ielim)%elim,mloc)
!%--------------------------------------------------------------
!%----------------------------------------by gas phase
ipri1=nrow+1
ipos1=this%pp%numsp
ipos2=0
do i=1,this%pp%numgasph
 ipos1=ipos1+1
 ipos2=ipos2+this%pp%numsp
 ipri2=ipri2+nrow
 mloc=matmul(this%pp%ueq,matmul(dcmob(ipos1:ipos2,:),dc1))
 dumob(ipri1:ipri2,:)=matmul(this%elim(ielim)%elim,mloc)
 ipri1=ipri2+1
 ipos1=ipos2
end do
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (dg,1,1,.false.)
call check_pointer_ (g,1,.false.)
call check_pointer_ (dc1,1,1,.false.)
call check_pointer_ (dcmob,1,1,.false.)
call check_pointer_ (mloc,1,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dumob_'
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
subroutine compute_dumob2_chemsysretraso &
   (this, &
    dumob, &
    nmobph, &
    nrow, &
    ncol, &
    dc, &
    numsp, &
    naqpri, &
    hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute U*dcmob/dc1
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso)       :: &
 this
integer                               :: &
 numsp, &
 naqpri, &
 hashcompz
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
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                    :: &
 msg
integer                               :: &
 ielim, &
 ipos1, &
 ipos2, &
 ipri1, &
 ipri2, &
 i
logical                               :: &
 be
real*8, pointer                       :: &
 mloc(:,:) => null()
!-------------------------------------------------------------------------
!
!   $code
!
 
 
 

 

!%------------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------------
call get_elim_index_ (this, ielim,hashcompz, be)
if (.not.be) then
 msg='Error, not defined components zone'
 call add_ (msg,hashcompz)
 goto 10
end if
!%------------------------------------------------------------
if (numsp.ne.this%pp%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------------
if (naqpri.ne.this%elim(ielim)%nrprisp) then
 msg='Error in number of primary species'
 goto 10
end if
!%------------------------------------------------------------
nrow=this%elim(ielim)%nrprisp
ncol=this%elim(ielim)%nrprisp
nmobph=1+this%pp%numgasph
!%------------------------------------------------------
call check_pointer_ (dumob,nmobph*nrow,ncol,.true.)
call check_pointer_ (mloc,this%pp%numaqprisp,ncol,.true.)
!%--------------------------------------for aqueous phase
call get_iposspsph_ (this%pp, this%pp%aqphindex, ipos1, ipos2)
ipri1=1
ipri2=nrow
mloc=matmul(this%pp%ueq(:,ipos1:ipos2),dc(ipos1:ipos2,:))
dumob(ipri1:ipri2,:)=matmul(this%elim(ielim)%elim,mloc)
!%----------------------------------------by gas phase
ipri1=nrow+1
do i=1,this%pp%numgasph
 
 call get_iposspsph_ (this%pp, this%pp%idgasph(i), ipos1, ipos2)
 ipri2=ipri2+nrow
 mloc=matmul(this%pp%ueq(:,ipos1:ipos2),dc(ipos1:ipos2,:))
 dumob(ipri1:ipri2,:)=matmul(this%elim(ielim)%elim,mloc)
 ipri1=ipri2+1
 
end do
!%------------------------------------------------------------
20 continue 
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (mloc,1,1,.false.)
if (iserror) goto 10
!%------------------------------------------------------------
return
 
10 continue 
print *,'***********************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dumob_'
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
subroutine compute_duads_chemsysretraso &
   (this, &
    duads, &
    nrow, &
    ncol, &
    nsp, &
    c, &
    hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute U*dcads/dc11ith
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in):: this

integer, intent(in)                       :: nsp

real*8, intent(in), dimension(nsp)        :: c 

real*8, pointer, dimension(:,:)           :: duads 

integer, intent(out)                      :: nrow

integer, intent(out)                      :: ncol

integer, intent(in)                       :: hashcompz

logical, intent(out)                      :: iserror      ! iserror=true, then there was an error. 
 
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
 isps1, &
 isps2, &
 ielim, &
 ireact
integer, pointer                      :: &
 idreact(:)
real*8, pointer                       :: &
 dc(:,:) => null(), &
 dg(:,:) => null(), &
 g(:) => null(), &
 cloc(:) => null(), &
 dionstr(:) => null(), &
 dc1(:,:) => null(), &
 dcloc(:) => null(), &
 duads_dc1(:,:) => null()
real*8                                :: &
 ionstr
logical                               :: &
 isbe, &
 isanomalous, &
 isupmxitergam
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------
iserror=.false.
msg=''
!%------------------------------------------------------
!% Check the number of species 
!%------------------------------------------------------
if (nsp/=this%pp%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%------------------------------------------------------
!% Check the components zone 
!%------------------------------------------------------
call get_elim_index_ (this,ielim,hashcompz,isbe)
if (.not.isbe) then
 msg='Error, not defined elimination matrix'
 call add_ (msg, ielim)
 goto 10
end if
!%------------------------------------------------------
nrow=this%pp%numaqprisp
ncol=this%elim(ielim)%nrprisp
call check_pointer_ (duads,nrow,ncol,.true.)
!%------------------------------------------------------
!% If numsurf=0, then return 
!%------------------------------------------------------
if (this%pp%numsurf==0) return 
!%------------------------------------------------------
!% Allocate local pointers 
!%------------------------------------------------------
call check_pointer_ (duads_dc1,nrow,nrow,.true.)
call check_pointer_ (dc,this%pp%numsp,this%pp%numaqprisp,.true.)
call check_pointer_ (dg,this%pp%numsp,this%pp%numaqprisp,.true.)
call check_pointer_ (g,this%pp%numsp,.true.)
call check_pointer_ (cloc,this%pp%numsp,.true.)
call check_pointer_ (dc1,this%pp%numaqprisp,this%elim(ielim)%nrprisp, .true.)
g=1.0d0
cloc=c
!%--------------------------------------by aqueous phase
call add_diagonal_ (dc(this%pp%idaqprisp,1:this%pp%numaqprisp),1.0d0,this%pp%numaqprisp)
!%-----------------------------------------------
!% Compute secondaries 
!%-----------------------------------------------
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
!%----------------------------------------------------
!% Compute dc1/dc11 
!%---------------------------------------------------- 
call compute_dc1_dc11_ &
    (this, &
     dc1, &
     cloc(this%pp%idaqprisp), &
     g(this%pp%idaqprisp), &
     dg(this%pp%idaqprisp,:), &
     this%elim(ielim)%nnrprisp, &
     this%elim(ielim)%nrprisp, &
     this%elim(ielim)%idrprisp, &
     this%elim(ielim)%idnrprisp, &
     this%elim(ielim)%ideqminreact, &
	 this%elim(ielim)%wcindex)
!%----------------------------------------------------
!%----------------------------------------------------
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
 
  call compute_dx_ &
  (this%pp%preaction(ireact)%ptr, &
   dcloc, &
   c(this%pp%idaqprisp), &
   g(this%pp%idaqprisp), &
   g(this%pp%idreactsp(ireact)), &
   dc(this%pp%idaqprisp,:), &
   dg(this%pp%idaqprisp,:), &
   dg(this%pp%idreactsp(ireact),:), &
   iserror, &
   c(this%pp%idreactsp(ireact)))
 
   if (iserror) goto 20
 
   dc(this%pp%idreactsp(ireact),:) = dcloc
 
 end do
 
 call get_iposspsurf_ (this%pp,i,isps1,isps2)
 
 duads_dc1=matmul(this%pp%ueq(:,isps1:isps2),dc(isps1:isps2,:))
 
 duads=duads+matmul(duads_dc1,dc1)
 
end do
 
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (idreact,1,.false.)
call check_pointer_ (dcloc,1,.false.)
call check_pointer_ (duads_dc1,1,1,.false.)
call check_pointer_ (g,1,.false.)
call check_pointer_ (cloc,1,.false.)
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dg,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (dc1,1,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_duads_'
print *, msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_usktrk1_chemsysretraso &
   (this, &
    usktrk, &
    naqpri, &
    c, &
    g, &
    area, &
    numsp, &
    hashcompz, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute E*U*Skt*rk
!
!   $Arguments:
!
 
 
type (t_chemicalsystemretraso):: &
 this
real*8, pointer              :: &
 usktrk(:)
integer                      :: &
 hashcompz, &
 naqpri, &
 numsp
real*8                       :: &
 c(numsp), &
 g(numsp), &
 area(numsp)
logical                      :: &
 iserror 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer              :: &
 vloc(:) => null ()
integer                      :: &
 ielim
logical                      :: &
 be
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
call get_elim_index_ (this, ielim, hashcompz, be)
if (.not.be) then
 msg='Error, not defined components zone'
 call add_ (msg,hashcompz)
 goto 10
end if
!%------------------------------------------------------------
call compute_usktrk_ &
   (this%pp, &
    usktrk, &
    naqpri, &
    c, &
    g, &
    area, &
    numsp, &
    iserror)
 
if (iserror) return
!%------------------------------------------------------------
call check_pointer_ (vloc,naqpri,.true.)
vloc=usktrk
naqpri=this%elim(ielim)%nrprisp
call check_pointer_ (usktrk,naqpri,.true.)
usktrk=matmul(this%elim(ielim)%elim,vloc)
!%------------------------------------------------------------
call check_pointer_ (vloc,1,.false.)
!%------------------------------------------------------------
return
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_usktrk_'
print *, msg
print *,'**************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_usktrk2_chemsysretraso &
   (this, &
    usktrk, &
    ncomp, &
    sktrk, &
    nsp, &
    hashcompz, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute E*U*Skt*rk
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in) :: this

real*8, pointer                            :: usktrk(:)

integer, intent(in)                        :: hashcompz

integer, intent(out)                       :: ncomp
 
integer, intent(in)                        :: nsp

real*8, intent(in)                         :: sktrk(nsp)

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
real*8, pointer              :: &
 vloc(:) => null()
integer                      :: &
 ielim
logical                      :: &
 be
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
!%------------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
call get_elim_index_ (this, ielim, hashcompz, be)
if (.not.be) then
 msg='Error, not defined components zone'
 call add_ (msg,hashcompz)
 goto 10
end if
!%------------------------------------------------------------
call compute_usktrk_ &
   (this%pp, &
    usktrk, &
    ncomp, &
    sktrk, &
    nsp, &
    iserror)
 
if (iserror) return
!%------------------------------------------------------------
call check_pointer_ (vloc,ncomp,.true.)
vloc=usktrk
ncomp=this%elim(ielim)%nrprisp
call check_pointer_ (usktrk,ncomp,.true.)
usktrk=matmul(this%elim(ielim)%elim,vloc)
!%------------------------------------------------------------
call check_pointer_ (vloc,1,.false.)
!%------------------------------------------------------------
return
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_usktrk_'
print *, msg
print *,'**************************'
iserror=.true.
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dusktrk1_chemsysretraso &
   (this, &
    dusktrk, &
    ncomp, &
    c, &
    alpha, &
    nsp, &
    hashcompz, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute E*U*Skt*drk/dc11
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in) :: this

integer, intent(in)                        :: hashcompz

integer, intent(out)                       :: ncomp

integer, intent(in)                        :: nsp

real*8, pointer                            :: dusktrk(:,:)

real*8, intent(in)                         :: c(nsp)

real*8, intent(in)                         :: alpha(nsp)

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
real*8, pointer              :: &
 cloc(:) => null (), &
 g(:) => null (), &
 dc(:,:) => null (), &
 dg(:,:) => null (), &
 dionstr(:) => null (), &
 dc1(:,:) => null (), &
 mloc(:,:) => null ()
real*8                       :: &
 ionstr
integer                      :: &
 i, &
 ielim
logical                      :: &
 be, &
 isanomalous, &
 isupmxitergam
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
call get_elim_index_ (this, ielim, hashcompz, be)
if (.not.be) then
 msg='Error, not defined components zone'
 call add_ (msg,hashcompz)
 goto 10
end if
!%-------------------------------------------------------------
call compute_dusktrk_ &
   (this%pp, &
    dusktrk, &
    ncomp, &
    c, &
    alpha, &
    nsp, &
    iserror)
 
if (iserror) return
 
!%----------------------------------------------------------
call check_pointer_ (cloc,nsp,.true.)
call check_pointer_ (g,nsp,.true.)
call check_pointer_ (dc,nsp,this%pp%numaqprisp,.true.)
call check_pointer_ (dg,nsp,this%pp%numaqprisp,.true.)
call check_pointer_ (dc1,this%pp%numaqprisp,this%elim(ielim)%nrprisp,.true.)
g=1.0d0
cloc=c
!%----------------------------------------------------------
do i=1,this%pp%numaqprisp
 dc (this%pp%idaqprisp(i),i)=1.0d0
end do
!%----------------------------------------------------------
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
 
  if (iserror.or.isupmxitergam) goto 10
!%--------------------------------------------------------
call compute_dc1_dc11_ &
    (this, &
     dc1, &
     cloc(this%pp%idaqprisp), &
     g(this%pp%idaqprisp), &
     dg(this%pp%idaqprisp,1:this%pp%numaqprisp), &
     this%elim(ielim)%nnrprisp, &
     this%elim(ielim)%nrprisp, &
     this%elim(ielim)%idrprisp, &
     this%elim(ielim)%idnrprisp, &
     this%elim(ielim)%ideqminreact, &
	 this%elim(ielim)%wcindex)
 
call check_pointer_ &
     (mloc, &
      this%pp%numaqprisp, &
      this%pp%numaqprisp, &
      .true.)
 
mloc=dusktrk
ncomp=this%elim(ielim)%nrprisp
call check_pointer_ (dusktrk,ncomp,ncomp,.true.)
dusktrk=matmul(this%elim(ielim)%elim,matmul(mloc,dc1))
!%----------------------------------------------------------
20 continue 
call check_pointer_ (mloc,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dg,1,1,.false.)
call check_pointer_ (g,1,.false.)
call check_pointer_ (dc1,1,1,.false.)
call check_pointer_ (cloc,1,.false.)
if (iserror) goto 10
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dusktrk_'
print *, msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dusktrk2_chemsysretraso &
   (this, &
    dusktrk, &
    naqpri, &
    dsktrk, &
    numsp, &
    hashcompz, &
    iserror)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Compute U*Skt*drk/dc11
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso):: &
 this
integer                      :: &
 hashcompz, &
 naqpri, &
 numsp
real*8, pointer              :: &
 dusktrk(:,:)
real*8                       :: &
 dsktrk(numsp,naqpri)
logical                      :: &
 iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer              :: &
 mloc(:,:) => null ()
real*8                       :: &
 ionstr
integer                      :: &
 i, &
 ielim
logical                      :: &
 be, &
 isanomalous
character(len=100)           :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%----------------------------------------------------------
iserror=.false.
msg=''
!%----------------------------------------------------------
call get_elim_index_ (this, ielim, hashcompz, be)
if (.not.be) then
 msg='Error, not defined components zone'
 call add_ (msg,hashcompz)
 goto 10
end if
!%----------------------------------------------------------
if (naqpri.ne.this%elim(ielim)%nrprisp) then
 msg='Error, number of primary aqueous species different to defined in the chemical system'
 goto 10
end if
!%-----------------------------------------------------------
call check_pointer_ (mloc,this%pp%numaqprisp,naqpri,.true.)
call check_pointer_ (dusktrk,naqpri,naqpri,.true.)
mloc=matmul(this%pp%ueq,dsktrk)
dusktrk=matmul(this%elim(ielim)%elim,mloc)
!%----------------------------------------------------------
call check_pointer_ (mloc,1,1,.false.)
!%-----------------------------------------------------------
return
 
10 continue 
print *,'**************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: compute_dusktrk_'
print *, msg
print *,'**************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine copy_chemsysretraso &
   (copied, &
    this)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in) :: this

type (t_chemicalsystemretraso), intent(out):: copied 
 
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
 ndim1, &
 ndim2
!-------------------------------------------------------------------------
!
!   $code
!
!%-----------------------------------------------------------
copied%numelim=this%numelim
!%-----------------------------------------------------------
if (copied%numelim>0) then
 allocate (copied%elim(copied%numelim))
 do i=1,copied%numelim
   ndim1=size(this%elim(i)%elim,1)
   ndim2=size(this%elim(i)%elim,2)
   copied%elim(i)%nrprisp=this%elim(i)%nrprisp 
   copied%elim(i)%nnrprisp=this%elim(i)%nnrprisp 
   copied%elim(i)%hashindex=this%elim(i)%hashindex
   copied%elim(i)%wcindex=this%elim(i)%wcindex
   copied%elim(i)%elim => null ()
   call check_pointer_ (copied%elim(i)%elim,ndim1,ndim2,.true.)
   copied%elim(i)%idrprisp => null () 
   call check_pointer_ (copied%elim(i)%idrprisp,ndim1,.true.)
   ndim1=copied%elim(i)%nnrprisp 
   if (ndim1>0) then 
	copied%elim(i)%ideqminreact => null ()
	call check_pointer_ (copied%elim(i)%ideqminreact,ndim1,.true.)
    copied%elim(i)%idnrprisp => null ()
    call check_pointer_ (copied%elim(i)%idnrprisp,ndim1,.true.)
   else
    copied%elim(i)%ideqminreact => null ()
    copied%elim(i)%idnrprisp => null ()
   end if 
   copied%elim(i)%elim=this%elim(i)%elim
   copied%elim(i)%idrprisp=this%elim(i)%idrprisp
   if (copied%elim(i)%nnrprisp>0) then 
	copied%elim(i)%ideqminreact=this%elim(i)%ideqminreact
    copied%elim(i)%idnrprisp=this%elim(i)%idnrprisp
   end if 
 end do
end if
!%-----------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine make_lin_trf_vector_chemsysretraso &
   (this, &
    vnew, &
    vold, &
    hashcompz, &
    iserror, &
	isueq)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in)  :: this

real*8, pointer, dimension(:)               :: vnew

real*8, intent(in), dimension(:)            :: vold 

integer, intent(in)                         :: hashcompz

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
 ielim, &
 ndim1
character(len=100)             :: &
 msg
logical                        :: &
 isbe, &
 haveisueq
real*8, pointer                :: &
 vloc(:) => null (), &
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
call get_elim_index_ (this, ielim,hashcompz,isbe)
 
if (.not.isbe) then
 msg='Error, not defined components zone'
 call add_ (msg,ielim)
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
ndim1=this%elim(ielim)%nrprisp
call check_pointer_ (vnew,ndim1,.true.)
call check_pointer_ (vloc,this%pp%numaqprisp,.true.)
vloc=matmul(u,vold)
vnew=matmul(this%elim(ielim)%elim,vloc)
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (vloc,1,.false.)
!%------------------------------------------------------------
!% Nullify local pointers 
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
subroutine make_lin_trf_array_chemsysretraso &
   (this, &
    anew, &
    aold, &
    hashcompz, &
    iserror, &
	isueq)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in)  :: this

real*8, pointer, dimension(:,:)             :: anew 

real*8, intent(in), dimension(:,:)          :: aold 

integer, intent(in)                         :: hashcompz

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
 n1, &
 n2, &
 ielim
character(len=100)             :: &
 msg
logical                        :: &
 be, &
 haveisueq
real*8, pointer                :: &
 aloc(:,:) => null (), &
 u(:,:) => null ()
!-------------------------------------------------------------------------
!
!   $code
!
 
 


!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
n1=size(aold,1)
if (n1.ne.this%pp%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%-------------------------------------------------------------
n2=size(aold,2)
if (n2==0) then
 msg='Error, number of column of input matrix = 0'
 goto 10
end if
!%-------------------------------------------------------------
call get_elim_index_ (this, ielim, hashcompz, be)
if (.not.be) then
 msg='Error, not defined components zone'
 call add_ (msg,ielim)
 goto 10
end if
!%------------------------------------------------------------
haveisueq=present(isueq)
!%------------------------------------------------------------
if (haveisueq.and..not.isueq) then
 u => this%pp%u
else
 u => this%pp%ueq 
end if
!%-------------------------------------------------------------
n1=this%elim(ielim)%nrprisp
call check_pointer_ (anew,n1,n2,.true.)
call check_pointer_ (aloc,this%pp%numaqprisp,n2,.true.)
aloc=matmul(u,aold)
anew=matmul(this%elim(ielim)%elim,aloc)
!%------------------------------------------------------------
!%Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (aloc,1,1,.false.)
!%------------------------------------------------------------
u => null ()
!%------------------------------------------------------------
return
!%------------------------------------------------------------
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
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_elim_index_chemsysretraso &
   (this, &
    elimindex, &
    hashcompz, &
    isbe)
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in) :: this

integer, intent(in)                        :: hashcompz

integer, intent(out)                       :: elimindex

logical, intent(out)                       :: isbe 
 
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
 i, &
 hashcompzloc 
!-------------------------------------------------------------------------
!
!   $code
!
!%-------------------------------------------------------------
isbe=.false.
elimindex=0
do i=1,this%numelim
 hashcompzloc=this%elim(i)%hashindex
 if (hashcompz==hashcompzloc) then
  elimindex=i
  isbe=.true.
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
subroutine specia_from_cpri_chemsysretraso &
   (this, &
    temp, &
    c, &
    g, &
	cold, &
    alpha, &
    nsp, &
    cprir, &
    nprir, &
    txoh, &
    c1, &
    c2, &
    spsurfarea, &
    ntxoh, &
    nsurf, &
    dtime, &
	volgas, &
    ionstr, &
    hashcompz, &
	faccap, &
    isanomalous, &
    isupmxitergam, &
    isreset, &
    iserror, &
    dc, &
    dg, &
    sktrk, &
    dsktrk)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Specia from concentration of primary species. 
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in)  :: this             ! Type chemical system (retraso) variable 

real*8, intent(in)                          :: temp             ! Temperature [c]

integer, intent(in)                         :: hashcompz         ! Index of components zone 

integer, intent(in)                         :: nsp              ! Number of species

integer, intent(in)                         :: nprir            ! Number of reduced primary species 

integer, intent(in)                         :: ntxoh            ! Total number of sites

integer, intent(in)                         :: nsurf            ! Number of surfaces 

real*8, intent(inout), dimension(nsp)       :: c                ! Vector of concentration of species 

real*8, intent(inout), dimension(nsp)       :: g                ! Vector of activity coefficients 

real*8, intent(in), dimension(nsp)          :: cold             ! Concentration vector in k

real*8, intent(in), dimension(nsp)          :: alpha            ! Reactive surface of minerals 

real*8, intent(in), dimension(nprir)        :: cprir            ! Concentration of primary species 

real*8, intent(in), dimension(ntxoh,nsurf)  :: txoh             ! Total of sites for adsorption 

real*8, intent(in), dimension(ntxoh,nsurf)  :: c1               ! 

real*8, intent(in), dimension(ntxoh,nsurf)  :: c2               ! 

real*8, intent(in), dimension(ntxoh,nsurf)  :: spsurfarea       ! 

real*8, intent(in)                          :: dtime            ! Time increment 

real*8, intent(in)                          :: faccap           ! Capillary correction for water activity 

real*8, intent(in)                          :: volgas           ! Gas volume 

logical, intent(in)                         :: isreset          ! isreset=true, then zeroing the vector c 

real*8, intent(out)                         :: ionstr           ! Ionic strength 

logical, intent(out)                        :: iserror          ! iserror=true, then there was an error 

logical, intent(out)                        :: isanomalous      ! isanomalous=true, then there was an anomalous concentrations 

logical, intent(out)                        :: isupmxitergam    ! isupmxitergam=true, then gamma iterations was exeeded  

real*8, pointer, optional, dimension(:)     :: sktrk            ! Change in the species due kinetic reactins [numsp]

real*8, pointer, optional, dimension(:,:)   :: dc               ! Derivatives of the concentrations with respect to primary species

real*8, pointer, optional, dimension(:,:)   :: dg               ! Derivatives of the activity coefficients with respect to primary species

real*8, pointer, optional, dimension(:,:)   :: dsktrk           ! d sktrk / d c1 
 
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
 i, &
 j, &
 iter, &
 ielim, &
 ipri, &
 isp, &
 neqmin, &
 isps1, &
 isps2
real*8, pointer                :: &
 dcloc(:,:) => null (), &
 dgloc(:,:) => null (), &
 dionstr(:) => null (), &
 dsktrkloc(:,:) => null (), &
 cpri(:) => null (), &
 gpri(:) => null (), &
 dcpri(:,:) => null (), &
 cold1(:) => null (), &
 error(:) => null ()
real*8                        :: &
 mxerror 
logical                       :: &
 isbe, &
 isderivates, &
 isconvergence, &
 havedc, &
 havedg, &
 havedsktrk
character(len=100)            :: &
 msg                    ! Error message 
real*8, parameter             :: &
 tolrel=1.0d-5, &       ! Relative tolerance
 r0=0.0d0, &
 r1=1.0d0
!-------------------------------------------------------------------------
!
!   $code
!

!%--------------------------------------------------------------
iserror=.false.
msg=' ' 
!%--------------------------------------------------------------
!% Check optional arguments 
!%--------------------------------------------------------------
havedc=present(dc)
havedg=present(dg)
havedsktrk=present(dsktrk)
!%--------------------------------------------------------------
!% Finf the elimination matrix 
!%--------------------------------------------------------------
call get_elim_index_ (this,ielim,hashcompz,isbe)
if (.not.isbe) then
 msg='Error, not defined components zone'
 call add_ (msg,hashcompz)
 goto 10
end if
!%--------------------------------------------------------------
!% Check number of primary species  
!%--------------------------------------------------------------
if (nprir/=this%elim(ielim)%nrprisp) then
 msg='Error in number of components'
 goto 10
end if
!%--------------------------------------------------------------
!% Check negative concentrations in primary species
!%--------------------------------------------------------------
if (minval(cprir)<r0) then
 msg='Error, negative concentrations in primary species'
 isanomalous=.true.
 goto 20
end if
!%--------------------------------------------------------------
!% Allocate local pointers 
!%--------------------------------------------------------------
call check_pointer_ (cpri,this%pp%numaqprisp,.true.)
call check_pointer_ (gpri,this%pp%numaqprisp,.true.)
call check_pointer_ (error,this%pp%numsp,.true.)
call check_pointer_ (cold1,this%pp%numsp,.true.)
!%--------------------------------------------------------------
!% Get firts and last global index of aqueous species
!%--------------------------------------------------------------
call get_iposspsph_ (this%pp,this%pp%aqphindex,isps1,isps2)
!%--------------------------------------------------------------
!% If reset zeroing c=0
!%--------------------------------------------------------------
if (isreset) c=r0
!%--------------------------------------------------------------
cpri(this%elim(ielim)%idrprisp(:))=cprir
gpri=g(this%pp%idaqprisp)
c(this%pp%idaqprisp)=cpri
!%--------------------------------------------------------------
!% For derivatives 
!%--------------------------------------------------------------
if (havedc.or.havedsktrk.or.havedg) then
  isderivates=.true.
!%--------------------------------------------------------------
!% Allocate local pointers 
!%--------------------------------------------------------------
  call check_pointer_ (dcloc,this%pp%numsp,this%pp%numaqprisp,.true.)
  call check_pointer_ (dgloc,this%pp%numsp,this%pp%numaqprisp,.true.)
!%--------------------------------------------------------------
!% Put 1 in the diagonal 
!%--------------------------------------------------------------
  do i=1,this%pp%numaqprisp
    dcloc(this%pp%idaqprisp(i),i)=r1
  end do
 
  call check_pointer_(dcpri,this%pp%numaqprisp,this%elim(ielim)%nrprisp,.true.)

else
!%--------------------------------------------------------------
!% Then not compute derivates 
!%--------------------------------------------------------------
  isderivates=.false. 
end if
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!% Start picard for gamma iterations 
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
!%--------------------------------------------------------------
iter=0 
do 

 iter=iter+1
 
 cold1(isps1:isps2)=c(isps1:isps2)
 
 call compute_c1_from_c11_ &
   (this, &
    cpri, &
    gpri, &
    this%elim(ielim)%nnrprisp, &
    this%elim(ielim)%nrprisp, &
    this%elim(ielim)%idrprisp, &
    this%elim(ielim)%idnrprisp, &
    this%elim(ielim)%ideqminreact, &
	this%elim(ielim)%wcindex, &
	msg, &
	iserror)

 if (iserror) goto 20 
 
 c(this%pp%idaqprisp)=cpri 
!%---------------------------------------------------------------------
!% Compute aqueous complexes 
!%---------------------------------------------------------------------
 call compute_secondaries_ &
   (this%pp, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    this%pp%aqphindex, &
    0, &
    ionstr, &
    dionstr, &
	faccap, &
	.false., &
    isanomalous, &
    isupmxitergam, &
    isderivates, &
    msg, &
    iserror)
 
 if (iserror.or.isanomalous) goto 20
 
 gpri=g(this%pp%idaqprisp)
!%---------------------------------------------------------------------
!% Compute relative error in the aqueous species
!%---------------------------------------------------------------------
 error(isps1:isps2)=(c(isps1:isps2)-cold1(isps1:isps2))/c(isps1:isps2) 
 mxerror=maxval(dabs(error(isps1:isps2)))
!%---------------------------------------------------------------------
!% Check convergence
!%---------------------------------------------------------------------
 isconvergence=(mxerror<=tolrel)
!%---------------------------------------------------------------------
!% Check the maximum number of gamma iterations 
!%---------------------------------------------------------------------
 isupmxitergam=(iter>mxitergam)
 
 if (isconvergence) exit 
 
 if (isupmxitergam) goto 20 
  
 
end do   
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
!% Finishing Picard for gamma iterations 
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
!%---------------------------------------------------------------------
!% If computing derivates is .true. then compute dc1/dc11.
!% Depend of the hash index
!%---------------------------------------------------------------------
if (isderivates) then
   call compute_dc1_dc11_ &
     (this, &
      dcpri, &
      cpri, &
      gpri, &
      dgloc(this%pp%idaqprisp,:), &
      this%elim(ielim)%nnrprisp, &
      this%elim(ielim)%nrprisp, &
      this%elim(ielim)%idrprisp, &
      this%elim(ielim)%idnrprisp, &
      this%elim(ielim)%ideqminreact, &
	  this%elim(ielim)%wcindex)
end if
!%-------------------------------------------------------------------
!% For gas phases
!%-------------------------------------------------------------------
do i=1,this%pp%numgasph
 call compute_secondaries_ &
   (this%pp, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    this%pp%aqphindex, &
    this%pp%idgasph(i), &
    ionstr, &
    dionstr, &
	1.0d0, &
	.false., &
    isanomalous, &
    isupmxitergam, &
    isderivates, &
    msg, &
    iserror)
 
 if (iserror.or.isanomalous) goto 20
 
 if (isderivates) then    
	call get_iposspsph_ (this%pp,this%pp%idgasph(i),isps1,isps2)
    c(isps1:isps2)=c(isps1:isps2)*volgas/(rgas*(temp+273.15d0))
    dcloc(isps1:isps2,1:this%pp%numaqprisp)= &
    dcloc(isps1:isps2,1:this%pp%numaqprisp)*volgas/(rgas*(temp+273.15d0))
 end if 

end do
!%------------------------------------------------------------------
!% Solve equilibrium for surfaces
!%------------------------------------------------------------------
do i=1,this%pp%numsurf
 
 call compute_secondaries_ &
   (this%pp, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    txoh(:,i), &
    c1(:,i), &
    c2(:,i), &
    spsurfarea(:,i), &
    ntxoh, &
    ionstr, &
    i, &
    isderivates, &
    msg, &
    iserror)
 
    if (iserror) goto 20
 
 
end do
!%----------------------------------------------------------------
!% Compute kinetic
!%----------------------------------------------------------------
call compute_kinetic_ &
   (this%pp, &
    c, &
    g, &
    dcloc, &
    dgloc, &
    alpha, &
    cold, &
    dtime, &
    msg, &
    iserror, &
    sktrk=sktrk, &
    dsktrk=dsktrkloc)
 
if (iserror) goto 20
!%-----------------------------------------------------------------------
!% Check NAN numbers in concentration and activity coefficients vectors
!%-----------------------------------------------------------------------    
do i=1,this%pp%numsp 
 if (isnan(c(i)).or.isnan(g(i))) then
  isanomalous=.true.
  goto 20
 end if
end do
!%--------------------------------------------------------------
!% If have dc
!%--------------------------------------------------------------
if (havedc) then
 call check_pointer_ (dc,this%pp%numsp,nprir,.true.)
 dc=matmul(dcloc,dcpri)
end if
!%--------------------------------------------------------------
!% If have dg
!%--------------------------------------------------------------
if (havedg) then
 call check_pointer_ (dg,this%pp%numsp,nprir,.true.)
 dg=matmul(dgloc,dcpri)
end if
!%--------------------------------------------------------------
!% If have dsktrk
!%--------------------------------------------------------------
if (havedsktrk) then
 call check_pointer_ (dsktrk,this%pp%numsp,this%elim(ielim)%nrprisp,.true.)
 dsktrk=matmul(dsktrkloc,dcpri)
 do i=1,this%pp%numsp 
  do j=1,this%elim(ielim)%nrprisp
   if (isnan(dsktrk(i,j))) then
    isanomalous=.true.
    goto 20
   end if
 end do
end do
end if
!%---------------------------------------------------------------
20 continue 
!%---------------------------------------------------------------
!% Deallocate local pointers 
!%---------------------------------------------------------------
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (dcpri,1,1,.false.)
call check_pointer_ (dcloc,1,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
call check_pointer_ (cpri,1,.false.)
call check_pointer_ (gpri,1,.false.)
call check_pointer_ (cold1,1,.false.)
call check_pointer_ (error,1,.false.)
call check_pointer_ (dsktrkloc,1,1,.false.)
if (iserror) goto 10
!%---------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: specia_from_cpri_'
print *, msg
print *,'********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine specia_aqueous_phase_chemsysretraso &
   (this, &
    c, &
    g, &
    nsp, &
    cpri, &
    npri, &
    ionstr, &
    hashcompz, &
    isanomalous, &
    isupmxitergam, &
    isreset, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Specia the aqueous phase 
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in) :: this

integer, intent(in)                        :: hashcompz

integer, intent(in)                        :: nsp

integer, intent(in)                        :: npri

real*8, intent(in), dimension(npri)        :: cpri

logical, intent(in)                        :: isreset

real*8, intent(inout), dimension(nsp)      :: c

real*8, intent(inout), dimension(nsp)      :: g 

real*8, intent(out)                        :: ionstr

logical, intent(out)                       :: iserror

logical, intent(out)                       :: isanomalous

logical, intent(out)                       :: isupmxitergam

!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)          :: &
 msg 
logical                     :: &
 isbe 
integer                     :: &
 i, &
 j, &
 isps, &
 isps1, &
 isps2, &
 ithelim, &
 naqx, &
 nunk, &
 iter, &
 ireact  
logical                     :: &
 isconv    
real*8, pointer             :: &
 s1(:,:) => null () , &
 s2(:,:) => null () , &
 jacobian(:,:) => null (), &
 residual(:) => null (), &
 dc(:,:) => null (), &
 dg(:,:) => null (), &
 dgloc(:,:) => null (), &
 gloc(:) => null (), &
 logk(:) => null ()
integer, pointer            :: &
 idaqx(:) => null ()
type(t_phase), pointer      :: &
 phase => null ()
type(t_reaction), pointer   :: &
 reaction 
!-------------------------------------------------------------------------
!
!   $code
!
msg=''
iserror=.false.
!%-----------------------------------------------------------
isanomalous=.false. 
isupmxitergam=.false. 
!%-----------------------------------------------------------
! Check if the components zone was defined previously 
!%-----------------------------------------------------------
call get_elim_index_ (this,ithelim,hashcompz,isbe)
if (.not.isbe) then
 msg='Error, not defined components zone'
 call add_ (msg,ithelim)
 goto 10
end if 
!%-----------------------------------------------------------
!% Check the number of species 
!%-----------------------------------------------------------
if (nsp/=this%pp%numsp) then
 msg='Error in number of species'
 goto 10
end if
!%-----------------------------------------------------------
!% Check the number of components associated to components 
!% zone 
!%-----------------------------------------------------------
if (npri/=this%elim(ithelim)%nrprisp) then
 msg='Error in number of components'
 goto 10
end if
!%-----------------------------------------------------------
!% Pointer to aqueous phase 
!% Get the global poistion of aqueous species 
!%-----------------------------------------------------------
phase => this%pp%pphase(this%pp%aqphindex)%ptr
call get_iposspsph_ (this%pp,this%pp%aqphindex,isps1,isps2)
!%-----------------------------------------------------------
!% Return the global indices of aqueous equilibrium reactions 
!%-----------------------------------------------------------
call get_idreaction_(this%pp,idaqx,naqx,this%pp%aqphindex, &
                     0,.false.)
!%-----------------------------------------------------------
!% Determine the number of unknowns 
!%-----------------------------------------------------------
nunk=this%elim(ithelim)%nnrprisp+naqx
!%-----------------------------------------------------------
if (nunk==0) goto 20 
!%-----------------------------------------------------------
!% Zeroing the concentrations vector (optional) 
!%-----------------------------------------------------------
if (isreset) c = 0.0d0 
!%-----------------------------------------------------------
!% Allocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (jacobian,nunk,nunk,.true.)
call check_pointer_ (residual,nunk,.true.)
call check_pointer_ (dc,this%pp%numsp,nunk,.true.)
call check_pointer_ (dg,this%pp%numsp,nunk,.true.)
call check_pointer_ (s1,nunk,npri,.true.)
call check_pointer_ (s2,nunk,nunk,.true.)
call check_pointer_ (logk,nunk,.true.)
!%-----------------------------------------------------------
!% Fill c and S1
!%-----------------------------------------------------------
do i=1,this%elim(ithelim)%nrprisp
 isps=this%elim(ithelim)%idrprisp(i)
 isps=this%pp%idaqprisp(isps)
 c(isps)=cpri(i) 
 do j=1,this%elim(ithelim)%nnrprisp
  ireact=this%elim(ithelim)%ideqminreact(j)
  s1(j,i)=this%pp%stq(ireact,isps)
 end do 
end do 
!%-----------------------------------------------------------
!% Fill dc and S2
!%-----------------------------------------------------------
do i=1,this%elim(ithelim)%nnrprisp
 isps=this%elim(ithelim)%idnrprisp(i)
 isps=this%pp%idaqprisp(isps)
 dc(isps,i)=1.0d0
 do j=1,this%elim(ithelim)%nnrprisp
  ireact=this%elim(ithelim)%ideqminreact(j)
  reaction => this%pp%preaction(ireact)%ptr
  call get_lnk_ (reaction,logk(i))
  s2(i,j)=this%pp%stq(ireact,isps)
 end do 
end do 
do i=1,naqx
 isps=idaqx(i)
 isps=this%pp%idreactsp(isps)
 dc(isps,i)=1.0d0
 do j=1,this%elim(ithelim)%nnrprisp
  ireact=idaqx(j)
  s2(this%elim(ithelim)%nnrprisp+j,this%elim(ithelim)%nnrprisp+i)= &
  this%pp%stq(ireact,isps)
 end do 
end do 
!%-----------------------------------------------------------
!% Start the iterative processes (Newton Raphson)
!%-----------------------------------------------------------
iter=0 
isconv=.false. 
do 
 iter=iter+1 
!%-----------------------------------------------------------
!% Compute activity coefficiets 
!%-----------------------------------------------------------
 call compute_act_coeff_(phase,gloc,c(isps1:isps2),iserror,ionstr)
 g(isps1:isps2)=gloc 
 if (iserror) goto 20 
 call compute_dact_coeff_ &
   (phase, &
    dgloc, &
    c(isps1:isps2), &
    dc(isps1:isps2,:), &
    iserror, &
    g=g(isps1:isps2))
 dg(isps1:isps2,:)=dgloc 
 if (iserror) goto 20 

 if (isconv) exit 

end do 
!%-----------------------------------------------------------
20 continue 
!%-----------------------------------------------------------
!% Deallocate local pointers 
!%-----------------------------------------------------------
call check_pointer_ (s1,1,1,.false.)
call check_pointer_ (s2,1,1,.false.)
call check_pointer_ (jacobian,1,1,.false.)
call check_pointer_ (residual,1,.false.)
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dg,1,1,.false.)
call check_pointer_ (gloc,1,.false.)
call check_pointer_ (dgloc,1,1,.false.)
call check_pointer_ (idaqx,1,.false.)
call check_pointer_ (logk,1,.false.)
phase => null () 
reaction => null () 
if (iserror) goto 10 
!%-----------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:', this%pp%name
print *,'Service: specia_aqueous_phase_'
print *, msg
print *,'********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_chem_info_chemsysretraso &
   (this, &
    hashcompz, &
    msg, &
    iserror, &
    numbase, &
    idbase, &
    namebase, &
    wcindex)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return chemical information
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in)           :: this

character(len=*), intent(out)                        :: msg

integer, intent(in)                                  :: hashcompz

logical, intent(out)                                 :: iserror

integer, pointer, optional, dimension(:)             :: idbase 

integer, intent(out), optional                       :: numbase

integer, intent(out), optional                       :: wcindex

character(len=100), pointer, optional, dimension(:)  :: namebase
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
logical                                   :: &
 havenumbase, &
 haveidbase, &
 havenamebase, &
 havewcindex, &
 isbe
integer                                   :: &
 ielim, &
 i, &
 npri, &
 ipri
integer, pointer                          :: &
 ptr(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!

!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
!% Check optional arguments
!%-------------------------------------------------------------
havenumbase=present(numbase)
haveidbase=present(idbase)
havenamebase=present(namebase)
havewcindex=present(wcindex)
!%-------------------------------------------------------------
call get_elim_index_ (this,ielim,hashcompz,isbe)
if (.not.isbe) then
    msg='Error, not defined components zone:'
    call add_ (msg,hashcompz)
    goto 10
end if
!%-------------------------------------------------------------
!% Give name of chemical base
!%-------------------------------------------------------------
if (havenamebase) then
  npri=this%elim(ielim)%nrprisp
  call check_pointer_ (namebase,npri,.true.)
  ptr =>  this%elim(ielim)%idrprisp
  do i=1,npri
   ipri=this%pp%idaqprisp(ptr(i))
   namebase(i)=this%pp%pspecies(ipri)%ptr%name
  end do
  ptr => null ()
end if
!%-------------------------------------------------------------
!% Give the number of base
!%-------------------------------------------------------------
if (havenumbase) then
 numbase=this%elim(ielim)%nrprisp
end if
!%-------------------------------------------------------------
!% Give index (pointer) to chemical base
!%-------------------------------------------------------------
if (haveidbase) then
  npri=this%elim(ielim)%nrprisp
  call check_pointer_ (idbase,npri,.true.)
  ptr => this%elim(ielim)%idrprisp
  idbase = this%pp%idaqprisp(ptr)
  ptr => null ()
end if
!%-------------------------------------------------------------
!% Give water component index
!%-------------------------------------------------------------
if (havewcindex) then
 wcindex=this%elim(ielim)%wcindex
end if
!%-------------------------------------------------------------
return
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_iumob_chemsysretraso &
   (this, &
    iumob, &
    naqcol, &
    ngascol, &
    nnonaqcol, &
    c, &
    numsp, &
    ithcomp, &
	hashcompz, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the ith mobile component
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in) :: this

real*8, pointer                            :: iumob(:)

integer, intent(in)                        :: numsp

integer, intent(in)                        :: ithcomp

integer, intent(in)                        :: hashcompz

integer, intent(out)                       :: naqcol

integer, intent(out)                       :: ngascol

integer, intent(out)                       :: nnonaqcol

real*8, intent(in)                         :: c(numsp)

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
logical             :: &
 isbe
character(len=100)  :: &
 msg
integer             :: &
 ithelim 
real*8, pointer     :: &
 u2(:,:) => null ()
integer             :: &
 isp1, &
 isp2, &
 i 
!-------------------------------------------------------------------------
!
!   $code
!
iserror=.false.
msg=''
!%------------------------------------------------------------
call get_elim_index_ (this,ithelim,hashcompz,isbe) 
if (.not.isbe) then
 msg='Error, not defined components zone'
 goto 10
end if 
!%------------------------------------------------------------
!% Check the number of species 
!%------------------------------------------------------------
if (numsp/=this%pp%numsp) then
 msg='Error, different number defined in the chemical system'
 goto 10
end if
!%------------------------------------------------------------
if (ithcomp<=0.or.ithcomp>this%elim(ithelim)%nrprisp) then
 msg='Error in component index'
 goto 10
end if
!%-------------------------------------------------------------
!% in the aqueous phase
!%-------------------------------------------------------------
call get_iposspsph_ (this%pp,this%pp%aqphindex,isp1,isp2)
naqcol=1
ngascol=this%pp%numgasph
nnonaqcol=0
!%-------------------------------------------------------------
call check_pointer_ (u2,this%elim(ithelim)%nrprisp,this%pp%numsp,.true.)
u2=matmul(this%elim(ithelim)%elim,this%pp%ueq)
!%-------------------------------------------------------------
call check_pointer_ (iumob,naqcol+ngascol+nnonaqcol,.true.)
iumob(naqcol)=dot_product(u2(ithcomp,isp1:isp2),c(isp1:isp2))
!%-------------------------------------------------------------
!% in the gas phases
!%-------------------------------------------------------------
do i=1,this%pp%numgasph
 call get_iposspsph_ (this%pp,this%pp%idgasph(i),isp1,isp2)
 iumob(naqcol+i)=dot_product(u2(ithcomp,isp1:isp2),c(isp1:isp2))
end do
!%------------------------------------------------------------
!% Deallocate local pointers 
!%------------------------------------------------------------
call check_pointer_ (u2,1,1,.false.)
!%------------------------------------------------------------
return
10 continue 
print *,'*******************************'
print *,'Chemical System:'
print *,'Namw:', this%pp%name
print *,'Service: compute_iumob_'
print *,msg
print *,'*******************************'
iserror=.true.
return
return
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dcads_chemsysretraso &
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
!   $Description: Compute dcads/dc1
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in)     :: this

integer, intent(in)                            :: numsp

integer, intent(in)                            :: naqpri

real*8, intent(in)                             :: dc(numsp,naqpri)

real*8, pointer                                :: dcads(:,:)

integer, intent(out)                           :: nrow

integer, intent(out)                           :: ncol

logical, intent(out)                           :: iserror 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
character(len=100)                    :: &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

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
subroutine write_chemsysretraso &
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
 
type (t_chemicalsystemretraso), intent(in) :: this

integer, intent(in)                        :: ioutput

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
 i, &
 j, &
 ipri, &
 isp, &
 k
real*8, pointer                 :: &
 u(:,:) => null ()
character (len=100), pointer    :: &
 namesp(:) => null ()
character (len=100)             :: &
 namepri, &
 msg 
!-------------------------------------------------------------------------
!
!   $code
!

!%-------------------------------------------------------------
iserror=.false.
msg=''
!%-------------------------------------------------------------
write(ioutput,*)'--------------------------------------------'
write(ioutput,*)' Chemical System (type retraso) Information'
write(ioutput,*)'--------------------------------------------'
!%-------------------------------------------------------------
!% Write parent chemical system
!%-------------------------------------------------------------
call write_ (this%pp,ioutput,iserror)
if (iserror) goto 10 
!%-------------------------------------------------------------
!% Get the name of species 
!%-------------------------------------------------------------
call get_chem_info_ (this%pp,msg,iserror,namesp=namesp)
!%-------------------------------------------------------------
write(ioutput,*)'---------------------------'
write(ioutput,*)'    Components matrices    '
write(ioutput,*)'---------------------------'
write(ioutput,4)'ncomp=',this%numelim
write(ioutput,*)'---------------------------'
!%-------------------
do j=1,this%numelim
 call check_pointer_ (u,this%elim(j)%nrprisp,this%pp%numsp,.true.)
 u=matmul(this%elim(j)%elim,this%pp%ueq)
 write(ioutput,*)'---------------------------'
 write(ioutput,4)'Comp=',j
 write(ioutput,*)'---------------------------'
 write(ioutput,*)'       Mineral set         '
 write(ioutput,*)'---------------------------'
 do i=1,this%elim(j)%nnrprisp
   isp=this%elim(j)%ideqminreact(i)
   isp=this%pp%idreactsp(isp)
   write(ioutput,3) namesp(isp)
 end do
 write(ioutput,*)'---------------------------'
 write(ioutput,*)'    Reduced Primary species'
 write(ioutput,*)'---------------------------'
 
 do i=1,this%elim(j)%nrprisp
  ipri=this%elim(j)%idrprisp(i)
  isp=this%pp%idaqprisp(ipri)
  write(ioutput,3) namesp(isp)
 end do
!%-------------------
if (this%elim(j)%nnrprisp>0) then
 
  write(ioutput,*)'----------------------------'
  write(ioutput,*)' Non Reduced Primary species'
  write(ioutput,*)'----------------------------'
 
 do i=1,this%elim(j)%nnrprisp
  ipri=this%elim(j)%idnrprisp(i)
  isp=this%pp%idaqprisp(ipri)
  write(ioutput,3) namesp(isp)
 end do
 
end if
!%------------------
write(ioutput,*)'------------------------------------------------'
write(ioutput,*)'    Components matrix'
write(ioutput,*)'------------------------------------------------'
 
 write(ioutput,2) 'species','=',(namesp(k),k=1,this%pp%numsp)
 
 do i=1,this%elim(j)%nrprisp
  ipri=this%elim(j)%idrprisp(i)
  isp=this%pp%idaqprisp(ipri)
  namepri=this%pp%pspecies(isp)%ptr%name
  write(ioutput,1) namepri,'=',(u(i,k),k=1,this%pp%numsp)
 end do
 
write(ioutput,*)'------------------------------------------------'
write(ioutput,*)'------------------------------------------------'
 
end do
!%-------------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (namesp,1,.false.)
call check_pointer_ (u,1,1,.false.)
!%-------------------------------------------------------------
return
1 format (a5,a3,<this%pp%numsp>e10.3)
2 format (a10,a3,<this%pp%numsp>a10)
3 format (a10)
4 format (a10,i5)
10 continue 
print *,'****************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: write_'
print *, msg
print*, '****************************'
iserror=.true.
return 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
!% This subroutine calculates the kernel to eliminate minerals in such a
!% way that Elim*stqB=0, using Singular Value Decomposition.
!% Inputs:  INDP     = Index for mineral concentrations not being zero
!%          IZ       = Elimination zone indicator
!%          NAQX     = Number of secondary species
!%          this%numminact  = Number of mineral concentrations not being z
!%          NPRI     = Number of primary species
!%          stqB     = Stoiquimetric matrix (only part of minerals)
!% Outputs: Elim    = Kernel
!% Author: Maarten W. retraso
!% Date: 15 may 1996
!%       16 oct 1996 SVD used instead of DLSVRR of IMSL
!%************************************************************
subroutine compute_elimination_matrix_chemsysretraso &
   (this, &
    stqb, &
    rkern, &
    nminact, &
    nrprisp, &
    msg, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in) :: this

real*8, pointer                            :: rkern(:,:)

integer, intent(in)                        :: nminact

integer, intent(in)                        :: nrprisp

real*8, intent(in)                         :: stqb(:,:)

logical, intent(out)                       :: iserror

character(len=*), intent(out)              :: msg 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                         :: &
 i,j, npri,irank,ierr,ii
integer iord(this%pp%numaqprisp)
real*8  work(this%pp%numaqprisp,this%pp%numaqprisp), &
        vv(this%pp%numaqprisp,this%pp%numaqprisp), &
        uu(this%pp%numaqprisp,this%pp%numaqprisp), &
        sv(this%pp%numaqprisp), &
        aa(this%pp%numaqprisp,this%pp%numaqprisp) 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

 
!%----------------------------------------------------------------------
msg=''
iserror=.false.
!%----------------------------------------------------------------------
npri=this%pp%numaqprisp
call check_pointer_ (rkern,npri-nminact,npri,.true.)
!%----------------------------------------------------------------------
if (nminact.eq.0) then
 
!%------------------------------------------------------------KERN = I
  do i=1,nrprisp
    do j=1,npri
      if (i.eq.j) then
        rkern(i,j) = 1.0
      else
        rkern(i,j) = 0.0
      end if
    end do
  end do
 
else
  iord=0 
  vv=0.0d0
  aa=0.0d0
  sv=0.0d0
  uu=0.0d0
  work=0.0d0
!%------------------------------------------Singular Value Decomposition
  do I=1,NMINACT
    do J=1,NPRI
      AA(I,J) = stqb(j,i)
    end do
  end do
  IERR = 0
  call SVD(npri,NMINACT,NPRI,AA,SV,.FALSE.,UU,.TRUE.,VV,IERR,WORK)
 
!%-------------------------------------------Search zero singular values
  IRANK = 0
  do I=1,NPRI
    IF (DABS(SV(I)).LT.1.0E-5) THEN
      IRANK = IRANK + 1
      IORD(IRANK) = I
    end if
  end do
  if (IRANK.NE.(NPRI-NMINACT)) then
   msg='Error, irank not equal to npri-nminact'
   goto 10
  end if
 
!%--------------------Construct rkern from right hand singular values
  do I=1,NPRI-NMINACT
    do J=1,NPRI
      II = IORD(I)
      rkern(I,J) = VV(J,II)
    end do
  end do
 
end if
!%----------------------------------------------------------------------
return
10 iserror=.true.
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine switch_base_chemsysretraso &
   (this, &
    naqprisp, &
	nadsprisp, &
    nameaqprisp, &
	nameadsprisp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Change the chemical base. 
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(inout)          :: this           ! Type chemical system (retraso) variable

integer, intent(in)                                    :: naqprisp       ! Number of aqueous primary species 

integer, intent(in)                                    :: nadsprisp      ! Number of adsorption primary species 

character(len=*), intent(in), dimension(naqprisp)      :: nameaqprisp    ! New set of aqueous primary species  

character(len=*), intent(in), dimension(nadsprisp)     :: nameadsprisp   ! New set of adsorption primary species 

logical, intent(out)                                   :: iserror        ! iserror=true, then there was an error 
 
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
 
!%---------------------------------------------------------
msg=''
iserror=.false.
!%------------------------------------------------------------
!% Call the corresponding service in the parent chemical 
!% system 
!%------------------------------------------------------------
call switch_base_ &
   (this%pp, &
    naqprisp, &
	nadsprisp, &
    nameaqprisp, &
	nameadsprisp, &
    iserror)
!%------------------------------------------------------------
return
 
10 continue 
print *,'********************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: set_chemical_base_'
print *,'********************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Public subroutine****************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine get_hashcompz_chemsysretraso &
   (this, &
    hashcompz, &
    c, &
    g, &
    nsp, &
    iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Return the hash index corresponding to components zone.
! 
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(inout):: this        ! Type retraso chemical system 

integer, intent(out)                         :: hashcompz    ! Hash index of components zone 

integer, intent(in)                          :: nsp         ! Number of species

real*8, intent(in), dimension (nsp)          :: c           ! Vector of the concentrations

real*8, intent(in), dimension (nsp)          :: g           ! Vector of activity coefficients 

logical, intent(out)                         :: iserror     ! iserror=true, then there was an error 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
integer                         :: &
 numeqmin, &
 numeqminglob, &
 numminglob, &
 i, &
 j, &
 imin, &
 ireact, &
 isp, &
 nnrprisp, &
 nrprisp
logical, pointer                :: &
 isbesps(:)=> null (), &
 isbeneg(:)=> null ()
character(len=100), pointer     :: &
 nameminglob(:) => null ()
logical                         :: &
 isbe, &
 isanomalous, &
 isupmxitergam
real*8, pointer                 :: &
 elim(:,:) => null (), &
 stqm(:,:) => null (), &
 colu(:,:) => null (), &
 si(:) => null (), &
 stq(:) => null (), &
 dc(:,:) => null (), &
 dionstr(:) => null (), &
 gloc(:) => null ()
integer, pointer                :: &
 ideqminglob(:) => null (), &
 ideqmin(:) => null (), &
 indr(:) => null (), &
 indnr(:) => null (), &
 ideqminreact(:) => null ()
character(len=100)              :: &
 msg
real*8                          :: &
 ionstr 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%-------------------------------------------------------------
msg=''
iserror=.false.
!%-------------------------------------------------------------
!% Check the number of species 
!%-------------------------------------------------------------
if (nsp/=this%pp%numsp) then
  msg='Error in number of species'
  goto 10
end if
!%-------------------------------------------------------------
!% Allocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (si,nsp,.true.)
call check_pointer_ (gloc,nsp,.true.)
call check_pointer_ (isbesps,nsp,.true.)
call check_pointer_ (isbeneg,nsp,.true.)
call check_pointer_ (indr,this%pp%numaqprisp,.true.)
call check_pointer_ (indnr,this%pp%numaqprisp,.true.)
gloc=g
hashcompz=0
!%-------------------------------------------------------------
call get_chem_info_ (this%pp,msg,iserror,ideqminsp=ideqminglob,neqminsp=numeqminglob)
if (iserror) goto 20
!%-------------------------------------------------------------
!% If there are not equilibrium mineral, then return 
!%-------------------------------------------------------------
if (numeqminglob==0) goto 20
!%--------------------------------------------------------------
!% Allocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (ideqmin,numeqminglob,.true.)
!%--------------------------------------------------------------
!% Compute saturation indices 
!%--------------------------------------------------------------
si=c
do i=1,this%pp%numminph
  call compute_secondaries_ &
     (this%pp, &
      si, &
      gloc, &
      dc, &
      dc, &
      this%pp%aqphindex, &
      this%pp%idminph(i), &
      ionstr, &
      dionstr, &
	  1.0d0, &
	  .false., &
      isanomalous, &
      isupmxitergam, &
      .false., &
      msg, &
      iserror)
 
   if (iserror.or.isupmxitergam.or.isanomalous) goto 20
 
end do
!%-------------------------------------------------------------
!% Check negative concentrations 
!%-------------------------------------------------------------
isbeneg=.false. 
where(c>0.0d0)
  isbesps=.true.
elsewhere (c<0.0d0)
  isbesps=.false.
  isbeneg=.true.
elsewhere (c==0.0d0)
  isbesps=.false.
end where
!%-------------------------------------------------------------
!% Check saturation indices of minerals 
!%-------------------------------------------------------------
numeqmin=0
isbe=.false.
do i=1,numeqminglob
 imin=ideqminglob(i)
 if (si(imin)>=1.0d0+this%pp%deltasatmin.and..not.isbesps(imin).and..not.isbeneg(imin)) then
   isbe=.true.
 else if (isbesps(imin)) then
   isbe=.true.
 end if
 if (isbe) then
   numeqmin=numeqmin+1
   ideqmin(numeqmin)=imin
   isbe=.false.
 end if
end do
!%----------------------------------------------------------
! Compute hash and check if exist the elimination matrix
!%-------------------------------------------------------------
call compute_hash_ (hashcompz,ideqmin(1:numeqmin),numeqmin,1)
!%----------------------------------------------------------
!% 
!%----------------------------------------------------------
if (numeqmin>0) then
 call check_pointer_ (ideqminreact,numeqmin,.true.)
 do i=1,numeqmin
  do j=1,this%pp%numreact
   if (ideqmin(i)==this%pp%idreactsp(j)) then
    ideqminreact(i)=j
    exit
   end if
  end do
 end do
end if
!%----------------------------------------------------------
call get_elim_index_ (this,i,hashcompz,isbe)
if (isbe) goto 20
!%----------------------------------------------------------
!% If not be the elimination matrix
!%----------------------------------------------------------
if (numeqmin>0) then
 call check_pointer_ (stqm,numeqmin,this%pp%numaqprisp,.true.)
 call check_pointer_ (colu,this%pp%numaqprisp,numeqmin,.true.)
 do i=1,numeqmin
  ireact=ideqminreact(i)
  call get_stq_ (this%pp%preaction(ireact)%ptr,stq)
  stqm(i,:)=stq
  isp=this%pp%idreactsp(ireact)
  colu(:,i)= this%pp%ueq(:,isp)
 end do 
!%-------------------------------------------------------------
!% Compute reduced and non reduced indices 
!%-------------------------------------------------------------
 call index_(numeqmin,indnr,indr,nnrprisp,nrprisp,stqm,this%pp%numaqprisp,this%pp%wcindex,msg,iserror)
 if (iserror) goto 20
else
  nrprisp=this%pp%numaqprisp
  nnrprisp=0
  do i=1,this%pp%numaqprisp
   indr(i)=i
  end do
end if
!%-------------------------------------------------------------
!% Compute elimination matrices 
!%-------------------------------------------------------------
call compute_elimination_matrix_ &
   (this, &
    colu, &
    elim, &
    nnrprisp, &
    nrprisp, &
    msg, &
    iserror)
if (iserror) goto 20
!%-------------------------------------------------------------
!% Add the elimination matrix 
!%-------------------------------------------------------------
call add_elimination_matrix_ &
   (this, &
    elim, &
    nrprisp, &
    nnrprisp, &
    indr, &
    indnr, &
    ideqminreact, &
    hashcompz)
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
!% Deallocate local pointers 
!%-------------------------------------------------------------
call check_pointer_ (elim,1,1,.false.)
call check_pointer_ (colu,1,1,.false.)
call check_pointer_ (ideqmin,1,.false.)
call check_pointer_ (ideqminreact,1,.false.)
call check_pointer_ (indnr,1,.false.)
call check_pointer_ (indr,1,.false.)
call check_pointer_ (ideqminglob,1,.false.)
call check_pointer_ (stq,1,.false.)
call check_pointer_ (stqm,1,1,.false.)
call check_pointer_ (si,1,.false.)
call check_pointer_ (isbesps,1,.false.)
call check_pointer_ (isbeneg,1,.false.)
call check_pointer_ (nameminglob,1,.false.)
call check_pointer_ (dc,1,1,.false.)
call check_pointer_ (dionstr,1,.false.)
call check_pointer_ (gloc,1,.false.)
if (iserror) goto 10
!%-------------------------------------------------------------
return
 
10 continue 
print *,'****************************'
print *,'Chemical System:'
print *,'Name:',this%pp%name
print *,'Service: get_hashcompz_'
print *, msg
print*, '****************************'
iserror=.true.
return
 
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_c1_from_c11_chemsysretraso &
   (this, &
    c1, &
    g1, &
    nnrprisp, &
    nrprisp, &
    idrprisp, &
    idnrprisp, &
    ideqminreact, &
	wcindex, &
	msg, &
	iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: This subroutine calculates the non reduced primary concentrations
!% from the reduced primary concentrations by means of the equations for
!% mineral equilibrium.
!% Inputs:  CP      = Primary concentrations (reduced)
!%          EKM     = Equilibrium constant of mineral reactions
!%          GAMP    = Activity coefficient of primary species
!%          INDNR   = Index for non reduced primary concentrations
!%          INDP    = Index for mineral concentrations not being zero
!%          INDR    = Index for reduced primary concentrations
!%          IZ      = Elimination zone indicator
!%          NB      = Number of base
!%          NINDNR  = Number of non reduced primary concentrations
!%          NINDP   = Number of mineral concentrations not being zero
!%          NINDR   = Number of reduced primary concentrations
!%          NPRI    = Number of primary species
!%          stqm    = Stoichiomatric matrix for minerals
!% Outputs: CP      = Primary concentrations (non-reduced)
!% Author: Maarten W. retraso
!% Date: 16 may 1996
!% Reformed to Fortran 90 by Sergio Andrés Bea Jofré
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in)  :: this

integer, intent(in)                         :: nrprisp

integer, intent(in)                         :: nnrprisp

integer, intent(in)                         :: idrprisp(nrprisp)

integer, intent(in)                         :: idnrprisp(:)

integer, intent(in)                         :: ideqminreact(:)

real*8, intent(inout)                       :: c1(this%pp%numaqprisp)

real*8, intent(in)                          :: g1(this%pp%numaqprisp) 

integer, intent(in)                         :: wcindex 

logical, intent(out)                        :: iserror

character(len=*), intent(out)               :: msg 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer                 :: &
 a(:,:) => null ()
real*8, pointer                 :: &
 b(:) => null (), &
 stq(:) => null ()
real*8                          :: &
 logk, &
 dd
integer                         :: &
 i, &
 ireact
integer, pointer                :: &
 indx(:) => null ()
real*8                          :: &
 cw
!-------------------------------------------------------------------------
!
!   $code
!

!%-----------------------------------------------------------
iserror=.false.
msg=' ' 
!%-----------------------------------------------------------
if (nnrprisp==0) return
!%-----------------------------------------------------------
!% If the water is a secundary species 
!%-----------------------------------------------------------
if (wcindex>0) then
 cw=c1(this%pp%wcindex)
 c1(this%pp%wcindex)=1.0d0
end if
!%-----------------------------------------------------------
call check_pointer_ (a,nnrprisp,nnrprisp,.true.)
call check_pointer_ (b,nnrprisp,.true.)
call check_pointer_ (indx,nnrprisp,.true.)
!%-------------------------------------------------------------
!% Build the system equations
!%-------------------------------------------------------------
do i=1,nnrprisp
 
 ireact = ideqminreact(i)
 
 call get_lnk_ (this%pp%preaction(ireact)%ptr, logk)
 call get_stq_ (this%pp%preaction(ireact)%ptr, stq)
 
 b(i)=dot_product(stq(1:this%pp%numaqprisp),dlog(g1))
 b(i)=b(i)+dot_product(stq(idrprisp),dlog(c1(idrprisp)))
 b(i)=b(i)-logk
 a(i,:)=-stq(idnrprisp)
 
end do
!%-------------------------------------------------------------
!% Solve the linear system
!%-------------------------------------------------------------
call ludcmp (a,nnrprisp,nnrprisp,indx,dd,msg,iserror)
if (iserror) goto 20  
call lubksb (a,nnrprisp,nnrprisp,indx,b)
!%-------------------------------------------------------------
!% Convert logarithms
!%-------------------------------------------------------------
c1(idnrprisp) = exp (b)
!%-------------------------------------------------------------
20 continue 
!%-------------------------------------------------------------
if (wcindex>0) c1(this%pp%wcindex)=cw
!%-------------------------------------------------------------
call check_pointer_ (a,1,1,.false.)
call check_pointer_ (b,1,.false.)
call check_pointer_ (indx,1,.false.)
call check_pointer_ (stq,1,.false.)
!%-------------------------------------------------------------------
return
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dc1_dc11_chemsysretraso &
    (this, &
     dc1, &
     c1, &
     g1, &
     dg1, &
     nnrprisp, &
     nrprisp, &
     idrprisp, &
     idnrprisp, &
     ideqminreact, &
	 wcindex)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: This subroutine calculates the derivatives of the non-reduced primary
!% concentration with respect to the reduced primary concentration.
!% Inputs:  CP      = Concentration of primary species
!%          DGAMP   = Derivatives of GAMP with respect to CP
!%          GAMP    = Activity coefficients of the primary species
!%          INDNR   = Index for non reduced primary concentrations
!%          INDP    = Index for mineral concentrations not being zero
!%          INDR    = Index for reduced primary concentrations
!%          IZ      = Elimination zone indicator
!%          NB      = Number of base
!%          NINDNR  = Number of non reduced primary concentrations
!%          NINDP   = Number of mineral concentrations not being zero
!%          NINDR   = Number of reduced primary concentrations
!%          stqm    = Stoichiomatric matrix for minerals
!% Outputs: DCNRE   = Derivatives of red. wrt non-red prim. species
!% Author: Maarten W. retraso
!% Date: 16 may 1996
!%       31 may 1996 derivatives of gammas added
!% Reformes to Fortran 90 by Sergio Andrés Bea (2006)
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(in)  :: this

integer, intent(in)                         :: nrprisp

integer, intent(in)                         :: nnrprisp

real*8, intent(out)                         :: dc1(this%pp%numaqprisp,nrprisp)

real*8, intent(in)                          :: c1(this%pp%numaqprisp)

real*8, intent(in)                          :: g1(this%pp%numaqprisp)

real*8, intent(in)                          :: dg1(this%pp%numaqprisp,this%pp%numaqprisp)

integer, intent(in)                         :: idrprisp(nrprisp)

integer, intent(in)                         :: idnrprisp(:)

integer, intent(in)                         :: ideqminreact(:) 

integer, intent(in)                         :: wcindex 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer                :: &
 a(:,:) => null (), &
 b(:,:) => null (), &
 dc12(:,:) => null (), &
 stq(:) => null (), &
 c1loc(:) => null ()
integer                        :: &
 ipri, &
 jpri, &
 kpri, &
 i, &
 j, &
 k, &
 n, &
 m, &
 ireact, &
 mmin, &
 info
integer, pointer               :: &
 ipiv(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!


!%---------------------------------------------------------------------
dc1=0.0d0
do i=1,nrprisp
 dc1(idrprisp(i),i)=1.0d0
end do
if (nnrprisp==0) return
!%---------------------------------------------------------------------
call check_pointer_ (c1loc,this%pp%numaqprisp,.true.)
c1loc=c1
if (wcindex>0) c1loc(this%pp%wcindex)=1.0d0
!%---------------------------------------------------------------------
call check_pointer_ (a,nnrprisp,nnrprisp,.true.)
call check_pointer_ (b,nnrprisp,nrprisp,.true.)
call check_pointer_ (ipiv,nnrprisp,.true.)
call check_pointer_ (dc12,nnrprisp,nrprisp,.true.)
!%----------------------------------------------------------------------
 
 
do i=1,nnrprisp
 
  ireact = ideqminreact(i)
 
  call get_stq_ (this%pp%preaction(ireact)%ptr, stq)
!%----------------------------------------------------------------------
!% Matrix of the linear system
!%----------------------------------------------------------------------
  do j=1,nnrprisp
    jpri = idnrprisp(j)
    a(i,j) = stq(jpri)
    do kpri=1,this%pp%numaqprisp
      a(i,j) = a(i,j) + stq(kpri)*dg1(kpri,jpri)*c1loc(jpri)/g1(kpri)
    end do
  end do
!%---------------------------------------------------------------------- 
!% Right hand sides of the linear system
!%----------------------------------------------------------------------
  do j=1,nrprisp
    jpri = idrprisp(j)
    b(i,j) = - stq(jpri)
    do kpri=1,this%pp%numaqprisp
      b(i,j) = b(i,j) - stq(kpri)*dg1(kpri,jpri)*c1loc(jpri)/g1(kpri)
    end do
  end do
 
 
 
end do
!%----------------------------------------------------------------------
!% Solve the linear system
!%----------------------------------------------------------------------
n = nnrprisp
m = nrprisp
mmin = nnrprisp
 
call f07adf(n,n,a,mmin,ipiv,info)
call f07aef('n',n,m,a,mmin,ipiv,b,mmin,info)
!%---------------------------------------------------------------------- 
!% Convert logarithms
!%----------------------------------------------------------------------
do i=1,nnrprisp
  ipri = idnrprisp(i)
  do j=1,nrprisp
    jpri = idrprisp(j)
    dc12(i,j) = b(i,j)*c1loc(ipri)/c1loc(jpri)
  end do
end do
 
dc1(idnrprisp,:) = dc12
!%----------------------------------------------------------------------
!% Deallocate local pointers 
!%----------------------------------------------------------------------
call check_pointer_ (a,1,1,.false.)
call check_pointer_ (b,1,1,.false.)
call check_pointer_ (dc12,1,1,.false.)
call check_pointer_ (stq,1,.false.)
call check_pointer_ (ipiv,1,.false.)
call check_pointer_ (c1loc,1,.false.)
!%----------------------------------------------------------------------
return
10 format (a5,i5,a2,i5,a3,e10.4)
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine compute_dc1_dc11b_chemsysretraso &
    (this, &
     dc1, &
     c1, &
     g1, &
     dg1, &
     nnrprisp, &
     nrprisp, &
     idrprisp, &
     idnrprisp, &
     ideqminreact, &
	 wcindex)
implicit none
!-------------------------------------------------------------------------
!
!   $Description: This subroutine calculates the derivatives of the non-reduced primary
! concentration with respect to the reduced primary concentration
!
!% Programmed in Fortran90 by Sergio Andrés Bea Jofré (2008) 
!
!   $Arguments:
!
type (t_chemicalsystemretraso), intent(in)  :: this

integer, intent(in)                         :: nrprisp

integer, intent(in)                         :: nnrprisp

real*8, intent(out)                         :: dc1(this%pp%numaqprisp,nrprisp)

real*8, intent(in)                          :: c1(this%pp%numaqprisp)

real*8, intent(in)                          :: g1(this%pp%numaqprisp)

real*8, intent(in)                          :: dg1(this%pp%numaqprisp,this%pp%numaqprisp)

integer, intent(in)                         :: idrprisp(nrprisp)

integer, intent(in)                         :: idnrprisp(:)

integer, intent(in)                         :: ideqminreact(:) 

integer, intent(in)                         :: wcindex
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
real*8, pointer               :: &
 a(:,:) => null (), &
 b(:,:) => null (), &
 ism(:) => null (), &
 sm(:,:) => null (), &
 prod(:,:) => null (), &
 prod2(:,:) => null (), &
 c1loc(:) => null ()
integer                        :: &
 ipri, &
 i, &
 n, &
 m, &
 ireact, &
 info
integer, pointer              :: &
 ipiv(:) => null () 
!-------------------------------------------------------------------------
!
!   $code
!
!%---------------------------------------------------------------------
dc1=0.0d0 
do i=1,nrprisp
 dc1(idrprisp(i),i)=1.0d0
end do
if (nnrprisp==0) return
!%---------------------------------------------------------------------
call check_pointer_ (c1loc,this%pp%numaqprisp,.true.)
c1loc=c1
if (wcindex>0) c1loc(this%pp%wcindex)=1.0d0
!%-------------------------------------------------
call check_pointer_ (ipiv,nnrprisp,.true.)
call check_pointer_ (a,nnrprisp,nnrprisp,.true.)
call check_pointer_ (prod,nrprisp,nrprisp,.true.)
call check_pointer_ (prod2,this%pp%numaqprisp,nrprisp,.true.)
call check_pointer_ (sm,nnrprisp,this%pp%numaqprisp,.true.)
call check_pointer_ (b,nnrprisp,nrprisp,.true.)
!%----------------------------------------------------------------------
do i=1,nnrprisp
  ireact = ideqminreact(i)
  call get_stq_ (this%pp%preaction(ireact)%ptr, ism)
  sm(i,:)=ism
end do
!%----------------------------------------------------------------------
do i=1,nrprisp
 ipri=idrprisp(i)
 prod(i,i)=1.0d0/c1loc(ipri)
end do
!%----------------------------------------------------------------------
do i=1,this%pp%numaqprisp
 prod2(i,:)=dg1(i,idrprisp)/g1(i)
end do
!%----------------------------------------------------------------------
b=matmul(sm(:,idrprisp),prod)+matmul(sm,prod2)
!%----------------------------------------------------------------------
a(:,:)=-sm(:,idnrprisp)
!%----------------------------------------------------------------------
n = nnrprisp
m = nrprisp
call f07adf(n,n,a,n,ipiv,info)
call f07aef('n',n,m,a,n,ipiv,b,n,info)
!%----------------------------------------------------------------------
do i=1,nnrprisp
 dc1(idnrprisp(i),:)=c1loc(idnrprisp(i))*b(i,:)
end do
!%----------------------------------------------------------------------
20 continue 
call check_pointer_ (a,1,1,.false.)
call check_pointer_ (ipiv,1,.false.)
call check_pointer_ (prod,1,1,.false.)
call check_pointer_ (prod2,1,1,.false.)
call check_pointer_ (c1loc,1,.false.)
call check_pointer_ (sm,1,1,.false.)
call check_pointer_ (b,1,1,.false.)
call check_pointer_ (ism,1,.false.)
!%----------------------------------------------------------------------
return
10 format (a5,i5,a2,i5,a3,e10.4)
end subroutine
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine index_chemsysretraso &
 (NMINACT, &
  INDNR, &
  INDR, &
  NINDNR, &
  NINDR, &
  STQM, &
  NPRI, &
  wcindex, &
  msg, &
  iserror)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description:This subroutine calculates the indices of reduced, non-reduced
!% primary concentrations and mineral concentrations.
!% Inputs:  INDP     = Index for mineral concentrations not being zero
!%          IZ       = Elimination zone indicator
!%          NB       = Number of base
!%          NPRI     = Number of primary species
!%          NMINACT  = Number of mineral concentrations not being zero
!%          stqm     = Stoiquiometric matrix of minerals
!% Outputs: INDNR    = Index for non reduced primary concentrations
!%          INDR     = Index for reduced primary concentrations
!%          NINDNR   = number of non reduced primary concentrations
!%          NINDR    = number of reduced primary concentrations
!% Author: Maarten W. retraso
!% Date: 15 may 1996
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

INTEGER   NINDNR, &
          NINDR, &
          NPRI, &
          wcindex, &
          nminact, &
          ipri, &
          imin, &
          jmin, &
          ncount, &
          jpri, &
          ifail
 
INTEGER              :: &
          INDNR(NPRI), &
          INDR(NPRI)
 
REAL*8               :: &
          AA(NMINACT,NMINACT), &
          IORD1(NPRI), &
          IORD2(NMINACT), &
          NUMMIN(NPRI), &
          WKSPCE(NMINACT), &
          STQM(NMINACT,NPRI)
character(len=100)                      :: &
 msg
logical                                 :: &
 iserror
 
real*8                                  :: &
 det
real*8, parameter                       :: &
 zero=1.0d-5
!%----------------------------------------------------------------------
msg=''
iserror=.false.
!%----------------------------------------------------------------------
!% Calculate number of minerals per component
!%----------------------------------------------------------------------
do ipri=1,npri
  nummin(ipri) = 0
  do imin=1,nminact
    if (dabs(stqm(imin,ipri))>=1.0d-5) then
      nummin(ipri) = nummin(ipri) + 1
    end if
  end do
end do
!%----------------------------------------------------------------------
!% Calculate components in order of preference
!%----------------------------------------------------------------------
ncount = 0
do imin=1,nminact
  do jpri=1,npri
    if (nummin(jpri)==imin) then
      ncount = ncount + 1
      iord1(ncount) = jpri
    end if
  end do
end do
do jpri=1,npri
  if (nummin(jpri)==0) then
    ncount = ncount + 1
    iord1(ncount) = jpri
  end if
end do
!%----------------------------------------------------------------------
!% Initialize IORD2
!%----------------------------------------------------------------------
do imin=1,nminact
  iord2(imin) = imin
end do
!%----------------------------------------------------------------------
10 do imin=1,nminact
  ipri=iord1(iord2(imin))
!%----------------------------------------------------------------------
!%----------------------------------------------------------------------
!%provi-----------------------------------------------------------------
!%----------------------------------------------------------------------
!%----------------------------------------------------------------------
!% No deja que el agua sea una secundaria. Esto da muchos problemas
!% si se lo deja arbitrario 
!%----------------------------------------------------------------------
 if (ipri==wcindex) then
   
   det=zero/10.0d0

   goto 30

end if
!%----------------------------------------------------------------------
!%----------------------------------------------------------------------
!%----------------------------------------------------------------------
!%----------------------------------------------------------------------
!%----------------------------------------------------------------------
!%----------------------------------------------------------------------
end do
!%----------------------------------------------------------------------
!% Calculate determinant
!%----------------------------------------------------------------------
do imin=1,nminact
  do jmin=1,nminact
    aa(imin,jmin) = stqm(imin,iord1(iord2(jmin)))
  end do
end do
ifail = 1
call f03aaf(aa,nminact,nminact,det,wkspce,ifail)
if (ifail>1) then
  msg='Error calculating determinant (subrout. index)'
  iserror=.true.
  return
end if
!%---------------------------------------------------------------------
30      if (dabs(det)<zero) then
!%----------------------------------------------------------------------
!% If determinant is zero, calculate other IORD2 and go back
!%----------------------------------------------------------------------
  do imin=1,nminact-1
    if(iord2(imin)<iord2(imin+1)-1) then
      iord2(imin)=iord2(imin)+1
      do jmin=1,imin-1
        iord2(jmin)=jmin
      end do
      goto 10
    end if
  end do
  if (iord2(nminact)<npri-1) then
    do imin=1,nminact-1
      iord2(imin)=imin
    end do
    iord2(nminact) = iord2(nminact) + 1
    goto 10
  end if
 
else
 
!%----------------------------------------------------------------------
!% If determinant is not zero, calculate INDR and  INDNR and quit
!%----------------------------------------------------------------------
  nindnr = nminact
  nindr = npri - nminact
  do imin=1,nminact
    indnr(imin) = iord1(iord2(imin))
  end do
  ncount = 0
  do ipri=1,npri
    do imin=1,nminact
      if(indnr(imin)==ipri) goto 20
    end do
    ncount = ncount + 1
    indr(ncount) = ipri
20  continue
  end do
 
  return
 
end if
!%----------------------------------------------------------------------
!% Error: no linear independent set can be found
!%----------------------------------------------------------------------
msg='Error, non linear independent set of primary species can be found'
iserror=.true.
!%-----------------------------------------------------------------
return
END SUBROUTINE
!%************************************************************
!%***************Private subroutine***************************
!%************************************************************
!%************************************************************
!%************************************************************
subroutine add_elimination_matrix_chemsysretraso &
   (this, &
    elim, &
    nrprisp, &
    nnrprisp, &
    idrprisp, &
    idnrprisp, &
    ideqminreact, &
    hashindex)
 
implicit none
!-------------------------------------------------------------------------
!
!   $Description: Add elimination matrix
!
!   $Arguments:
!
 
type (t_chemicalsystemretraso), intent(inout) :: this

integer, intent(in)                           :: nrprisp

integer, intent(in)                           :: nnrprisp

integer, intent(in)                           :: hashindex

integer, intent(in)                           :: idrprisp(nrprisp)

integer, intent(in)                           :: idnrprisp(nnrprisp)

integer, intent(in)                           :: ideqminreact(nnrprisp)

real*8, intent(in)                            :: elim(nrprisp,this%pp%numaqprisp) 
 
!-------------------------------------------------------------------------
!
!   $Pre-cond:
!
!   $Post-cond:
!
!   $License:
!
!-------------------------------------------------------------------------
type(t_elimination), pointer    :: &
 elimination(:) => null ()
integer                         :: &
 i, &
 ndim1, &
 ndim2, &
 ndim3 
!-------------------------------------------------------------------------
!
!   $code
!
 
 

 

!%-----------------------------------------------------------
if (this%numelim>0) then
!%-----------------------------------------------------------
 allocate (elimination(this%numelim))
 
 do i=1,this%numelim
  elimination(i)%nrprisp=this%elim(i)%nrprisp
  elimination(i)%nnrprisp=this%elim(i)%nnrprisp
  elimination(i)%wcindex=this%elim(i)%wcindex
  elimination(i)%hashindex=this%elim(i)%hashindex
  ndim1=elimination(i)%nrprisp
  ndim2=this%pp%numaqprisp
  ndim3=elimination(i)%nnrprisp
  allocate(elimination(i)%elim(ndim1,ndim2))
  elimination(i)%elim=this%elim(i)%elim
  if (elimination(i)%nnrprisp>0) then
   allocate (elimination(i)%ideqminreact(ndim3))
   allocate (elimination(i)%idnrprisp(ndim3))
   elimination(i)%idnrprisp=this%elim(i)%idnrprisp
   elimination(i)%ideqminreact=this%elim(i)%ideqminreact
  end if
  allocate (elimination(i)%idrprisp(ndim1))
  elimination(i)%idrprisp=this%elim(i)%idrprisp
 end do
 
 do i=1,this%numelim
  if (this%elim(i)%nnrprisp>0) then
   deallocate (this%elim(i)%ideqminreact)
   deallocate (this%elim(i)%idnrprisp)
  end if
  deallocate (this%elim(i)%idrprisp)
  deallocate (this%elim(i)%elim)
 end do
 
 deallocate (this%elim)
 this%numelim=this%numelim+1
 allocate(this%elim(this%numelim))
 
 do i=1,this%numelim-1
  this%elim(i)%nrprisp=elimination(i)%nrprisp
  this%elim(i)%nnrprisp=elimination(i)%nnrprisp
  this%elim(i)%wcindex=elimination(i)%wcindex
  this%elim(i)%hashindex=elimination(i)%hashindex
  ndim1=this%elim(i)%nrprisp
  ndim2=this%pp%numaqprisp
  ndim3=this%elim(i)%nnrprisp
  allocate(this%elim(i)%elim(ndim1,ndim2))
  this%elim(i)%elim=elimination(i)%elim
 
  if (this%elim(i)%nnrprisp>0) then
   allocate (this%elim(i)%ideqminreact(ndim3))
   allocate (this%elim(i)%idnrprisp(ndim3))
   this%elim(i)%idnrprisp=elimination(i)%idnrprisp
   this%elim(i)%ideqminreact=elimination(i)%ideqminreact
  else
   this%elim(i)%idnrprisp => null ()
   this%elim(i)%ideqminreact => null ()
  end if
 
  allocate (this%elim(i)%idrprisp(ndim1))
  this%elim(i)%idrprisp=elimination(i)%idrprisp
 end do
 
 
 this%elim(this%numelim)%nrprisp=nrprisp
 this%elim(this%numelim)%nnrprisp=nnrprisp
 this%elim(this%numelim)%hashindex=hashindex
 
 ndim1=this%pp%numaqprisp
 
 allocate(this%elim(this%numelim)%elim(nrprisp,ndim1))
 this%elim(this%numelim)%elim=elim
 if (nnrprisp>0) then
   allocate (this%elim(this%numelim)%ideqminreact(nnrprisp))
   allocate (this%elim(this%numelim)%idnrprisp(nnrprisp))
   this%elim(this%numelim)%idnrprisp=idnrprisp
   this%elim(this%numelim)%ideqminreact=ideqminreact
 else
   this%elim(this%numelim)%idnrprisp => null ()
   this%elim(this%numelim)%ideqminreact => null ()
 end if
 
  allocate (this%elim(this%numelim)%idrprisp(nrprisp))
  this%elim(this%numelim)%idrprisp=idrprisp
 
 
 
 do i=1,this%numelim-1
  if (elimination(i)%nnrprisp>0) then
   deallocate (elimination(i)%ideqminreact)
   deallocate (elimination(i)%idnrprisp)
  end if
  deallocate (elimination(i)%idrprisp)
  deallocate (elimination(i)%elim)
 end do
 
 deallocate (elimination)
 
 
 this%elim(this%numelim)%wcindex=0
 if (this%pp%wcindex>0) then
  do i=1,this%elim(this%numelim)%nrprisp
   if (this%elim(this%numelim)%idrprisp(i)==this%pp%wcindex) then
     this%elim(this%numelim)%wcindex=i
     exit
   end if
  end do
 end if
 
 
!%----------------------------------------------------------- 
else
!%----------------------------------------------------------- 
 
 this%elim(1)%nrprisp=nrprisp
 this%elim(1)%nnrprisp=nnrprisp
 this%elim(1)%hashindex=hashindex
 allocate(this%elim(1)%elim &
  (this%elim(1)%nrprisp,this%pp%numaqprisp))
  this%elim(1)%elim=elim
 if (nnrprisp>0) then
   allocate (this%elim(1)%ideqminreact(nnrprisp))
   allocate (this%elim(1)%idnrprisp(nnrprisp))
   this%elim(1)%idnrprisp=idnrprisp
   this%elim(this%numelim)%ideqminreact=ideqminreact
 else
   this%elim(1)%idnrprisp => null ()
   this%elim(1)%ideqminreact => null ()
 end if
 
  allocate (this%elim(1)%idrprisp(nrprisp))
  this%elim(1)%idrprisp=idrprisp
 
 
  this%elim(1)%wcindex=0
  if (this%pp%wcindex>0) then
  do i=1,this%elim(1)%nrprisp
   if (this%elim(1)%idrprisp(i)==this%pp%wcindex) then
     this%elim(1)%wcindex=i
     exit
   end if
  end do
 end if
 
 
end if
 
 
 
 
!%-------------------------------------------------------------------
return
end subroutine
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
!%************************************************************
end module m_chemicalsystemretraso
