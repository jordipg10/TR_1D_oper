!subroutine read_txt_pchemsys &
!>   (this, &
!>    namefile, &
!	namethdb, &
!>    namekindb, &
!	filebase, &
!	ioptxhl, &
!>    iserror)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Read the chemical system from txt file. 
!!> Read the chemical input file performed by Retraso 
!!
!!>   $Arguments:
!!
!> 
!> 
!type (t_parentchemicalsystem), intent(inout):: this
!
!character(len=*), intent(in)                :: namefile     !> Name of chemical system definition file
!
!character(len=*), intent(in)                :: namethdb     !> Name of thermodynamic database
!
!character(len=*), intent(in)                :: namekindb    !> Name of kinetic database
!
!character(len=*), intent(in)                :: filebase
!
!integer, intent(in)                         :: ioptxhl      !> Option to compute slat mass fraction 
!
!logical, intent(out)                        :: iserror 
!
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!character(len=100)  :: &
!> msg
!character(len=100)   :: &
!> namtyp
!character(len=100), pointer:: &
!> napri(:) => null(), & 
!> naaqx(:) => null(), &
!> naaqt(:) => null(), &
!> nagas(:) => null(), &
!> naads(:) => null(), &
!> namin(:) => null(), &
!> label(:) => null(), &
!> labelm(:) => null(), &
!> naadsmod(:) => null()
!integer, pointer          :: &
!> idmeq(:) => null(), &
!> nadsmod(:) => null()
!real*8                    :: &
!> tempisoterm 
!integer                   :: &
!> iact, &
!> iconv, &
!> ngas, &
!> naqx, &
!> naqt, &
!> nads, &
!> npri, &
!> naqpri, &
!> nadspri, &
!> nmin, &
!> itemp, &
!> nmineq, &
!> nminkin, &
!> nmod, &
!> lbase
!integer, parameter        :: &
!> mxsp=50, &
!> mxlabel=5
!character(len=100):: &
!> name
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!!%--------------------------------------------------------
!iserror=.false.
!msg=''
!!%---------------------------------------------------------
!!%Open the unit
!!%---------------------------------------------------------
!name=namefile
!call lastletter_ (lbase,filebase)
!name=filebase(1:lbase)//name
!open(unit=1,file=name,status='old',err=30)
!!%---------------------------------------------------------
!!%Allocate local pointers 
!!%---------------------------------------------------------
!call check_pointer_ (napri,mxsp,.true.)
!call check_pointer_ (naaqx,mxsp,.true.)
!call check_pointer_ (naaqt,mxsp,.true.)
!call check_pointer_ (nagas,mxsp,.true.)
!call check_pointer_ (namin,mxsp,.true.)
!call check_pointer_ (label,mxlabel,.true.)
!call check_pointer_ (naadsmod,mxsp,.true.)
!call check_pointer_ (labelm,mxlabel,.true.)
!call check_pointer_ (naads,mxsp,.true.)
!call check_pointer_ (idmeq,mxsp,.true.)
!call check_pointer_ (nadsmod,mxsp,.true.)
!!%---------------------------------------------------------
!!%Read the aqueous species defined in the chemical system
!!%---------------------------------------------------------
!call read_aqsp_sys_ (itemp,iact,iconv,naqt,npri,tempisoterm, &
!>                     label,naaqt,mxsp,mxlabel,iserror)
!if (iserror) goto 20          
!!%---------------------------------------------------------
!!%Read mineral and gases defined in the chemical system
!!%---------------------------------------------------------
!call read_miga_sys_ &
!>  (ngas,nmin,nmineq,nminkin,naqt,naaqt,idmeq,labelm,nagas,namin, &
!>   mxlabel,mxsp,iserror)
!if (iserror) goto 20  
!!%---------------------------------------------------------
!!% Compute the number of primary aqueous species 
!!%---------------------------------------------------------
!naqx=naqt-npri
!napri(1:npri)=naaqt(1:npri)
!naaqx(1:naqx)=naaqt(npri+1:npri+naqx)
!!%---------------------------------------------------------
!!%Read the surface complexes defined in the chemical system
!!%---------------------------------------------------------
!call read_ads_sys_ (nads,nmod,npri,napri,naads,naadsmod, &
!>                    nadsmod,mxsp,mxlabel,iserror)
!if (iserror) goto 20   
!!%---------------------------------------------------------
!!% set the chemical system object according thermodynamic
!!% data base 
!!%---------------------------------------------------------
!nadspri=nmod 
!naqpri=npri-nadspri
!select case (namethdb)
!case ('master25.dat','MASTER25.DAT')
!>  call set_from_master25_ &
!>   (this, &
!>    tempisoterm, &
!>    iact, &
!>    iconv, &
!>    ioptxhl, &
!>    naqpri, &
!>    nadspri, & 
!>    naqx, &
!>    nmin, &
!>    ngas, &
!>    nads, &
!>    nmod, &
!>    nadsmod, &
!>    napri, &
!>    naaqx, &
!>    namin, &
!>    nagas, &
!>    naads, &
!>    naadsmod, &
!>    idmeq, &
!>    namekindb, &
!	filebase, &
!>    iserror)
!case ('phreeqc.dat','PHREEQC.DAT')
!>  call set_from_phreeqc_ &
!>   (this, &
!>    tempisoterm, &
!>    iact, &
!>    iconv, &
!>    ioptxhl, &
!>    naqpri, &
!>    nadspri, & 
!>    naqx, &
!>    nmin, &
!>    ngas, &
!>    nads, &
!>    nmod, &
!>    nadsmod, &
!>    napri, &
!>    naaqx, &
!>    namin, &
!>    nagas, &
!>    naads, &
!>    naadsmod, &
!>    idmeq, &
!>    namekindb, &
!	filebase, &
!>    iserror)    
!> case default
!>  msg='Error, not defined thermodynamic data base:'
!>  call add_ (msg,namethdb)
!>  iserror=.true.
!>  goto 20 
!> end select 
!!%---------------------------------------------------------
!!% Close unit
!!%---------------------------------------------------------
!close(unit=1)
!!%---------------------------------------------------------
!20 continue
!!%---------------------------------------------------------
!!% Deallocate local pointers 
!!%--------------------------------------------------------- 
!call check_pointer_ (napri,1,.false.)
!call check_pointer_ (naaqx,1,.false.)
!call check_pointer_ (naaqt,1,.false.)
!call check_pointer_ (nagas,1,.false.)
!call check_pointer_ (namin,1,.false.)
!call check_pointer_ (label,1,.false.)
!call check_pointer_ (naads,1,.false.)
!call check_pointer_ (labelm,1,.false.)
!call check_pointer_ (naadsmod,1,.false.)
!call check_pointer_ (idmeq,1,.false.)
!if (iserror) goto 10 
!return
!> 
!10 continue 
!print *,'******************'
!print *,'Chemical System:'
!print *,'Name:',this%name
!print *,'service: read_txt_'
!print *, msg
!print *,'******************'
!iserror=.true.
!return
!> 
!30 continue
!msg='Error when open file:'
!call add_ (msg,namefile)
!goto 10  
!
!
!end subroutine