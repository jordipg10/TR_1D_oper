!module m_species
!!-------------------------------------------------------------------------
!!
!!>   $Description: The species class describes the intrinsic properties of a geochemical species.
!!> The actual intrinsic properties depend on the nature of the species. Aqueous
!!> species are described by their electrical charge, ion size or molecular weight.
!!> On the other hand, mineral species are described by the molar volume. The
!!> only operations of the species class are devoted to access these attributes.
!!
!!>   $Use: use flib_xpath
!!> use flib_sax
!!> use m_general_tools_cheproo
!!> use m_constants
!!
!!>   $Author: Sergio Andr�s Bea Jofr� 
!!
!!>   $License: CSIC-UPC
!!
!!-------------------------------------------------------------------------
!!%-------------------------------------------------------------------------
!!% Modules corresponding to CHEPROO project
!!%-------------------------------------------------------------------------
!use m_general_tools_cheproo
!use m_constants_cheproo
!!%-------------------------------------------------------------------------
!!% Modules corresponding to xml 
!!%-------------------------------------------------------------------------
!use flib_xpath
!use flib_sax
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!private         ::
!!%------------------------------------------------------------------------
!!% Public services 
!!%------------------------------------------------------------------------
!public:: &
!create_ &        !> Create species object.
!,read_xml_ &     !> Read species object from xml file
!,destroy_ &      !> Destroy species object 
!,get_name_ &     !> Return the name of the species
!,get_prop_ &     !> Return the value of some species property from its name 
!,add_prop_ &     !> Add new species property 
!,set_prop_ &     !> Set the name and values of the species properties.
!,set_name_ &     !> Set the name of the species object. 
!,write_ &        !> Write on the ioutput unit all attributtes of the species object. 
!,assignment(=) & !> Copy a species object in other species object. 
!,verify_         !> Verify different attributtes encapsulated in the species object. 
!!%------------------------------------------------------------------------
!!% Private services 
!!%------------------------------------------------------------------------
!private                   :: &
!read_xml_loc_  &
!,begin_element_handler
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!% Type pointer to species object 
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!type, public:: t_pspecies
!
!type(t_species), pointer:: ptr
!
!end type
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!% Type variable definition 
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!type, public::t_species
!> 
!character(len=100)                        :: name     !> Name of the species 
!> 
!real*8, pointer, dimension(:)             :: prop     !> Values of the species properties [numprop]
!> 
!character(len=100), pointer, dimension(:) :: nameprop !> Name of species properties [numprop]
!
!character(len=100), pointer, dimension(:) :: unitprop !> Unit of the properties [numprop]
!> 
!integer                                   :: numprop  !> Number of species properties
!> 
!end type t_species
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface create_
!> 
!module procedure create_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface read_xml_
!> 
!module procedure read_xml_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface destroy_
!> 
!module procedure destroy_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface get_name_
!> 
!module procedure get_name_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface get_prop_
!> 
!module procedure get_prop_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface add_
!> 
!module procedure add_prop_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface set_prop_
!> 
!module procedure set_prop_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface set_name_
!> 
!module procedure set_name_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface write_
!> 
!module procedure write_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface assignment (=)
!> 
!module procedure copy_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface verify_
!> 
!module procedure verify1_sps
!module procedure verify2_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!% Private subroutines 
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface read_xml_loc_
!> 
!module procedure read_xml_loc_sps
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!interface begin_element_handler
!> 
!module procedure begin_element_handler
!> 
!end interface
!!%------------------------------------------------------------------------
!!%------------------------------------------------------------------------
!contains
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine create_sps &
!>  (this)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Create species object. 
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(inout) :: this !> Type species variable 
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!> 
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!> 
!
!!%------------------------------------------------------------
!this%prop => null ()
!this%nameprop => null ()
!this%unitprop => null ()
!!%------------------------------------------------------------
!this%name=' '
!this%numprop=0
!!%------------------------------------------------------------
!return
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine read_xml_sps &
!>  (this, &
!>   namefile, &
!>   iserror)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Read species object from xml file
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(inout)      :: this     !> Type species variable 
!
!character(len=*), intent(in)        :: namefile !> Path and name of .xml file 
!
!logical, intent(out)                :: iserror  !> iserror=true there as an error
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!integer:: &
!> iostat
!type(xml_t):: &
!> fxml
!character(len=100)  :: &
!> name
!type(dictionary_t)  :: &
!> attributes
!character(len=100)  :: &
!> msg 
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!
!iserror=.false.
!msg=''
!!%-----------------------------------------------------------
!!% Open xml file
!!%-----------------------------------------------------------
!call open_xmlfile(namefile, fxml, iostat)
!if (iostat /= 0) then
!> msg='Error when open file:'
!> call add_ (msg,namefile)
!> goto 10
!end if
!!%-----------------------------------------------------------
!call xml_parse(fxml,begin_element_handler = begin_element_handler)
!!%-----------------------------------------------------------
!!% End and close xml file
!!%-----------------------------------------------------------
!call endfile_xmlfile(fxml)
!call close_xmlfile(fxml)
!!%-----------------------------------------------------------
!!% 
!!%-----------------------------------------------------------
!call read_xml_loc_ (name, attributes,this,msg,iserror)
!if (iserror) goto 10
!!%-----------------------------------------------------------
!> 
!return
!> 
!10 continue 
!print *,'********************'
!print *,'Species:'
!print *,'Name:', this%name
!print *,'Service: read_xml_'
!print *, msg
!print *,'********************'
!iserror=.true.
!return
!> 
!> 
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine destroy_sps &
!>  (this)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Destroy species object
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(inout) :: this !> Type species variable 
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!> 
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!
!this%name= ' '
!this%numprop=0
!!%-----------------------------------------------------------
!!% Deallocate pointers 
!!%-----------------------------------------------------------
!call check_pointer_ (this%nameprop,1,.false.)
!call check_pointer_ (this%unitprop,1,.false.)
!call check_pointer_ (this%prop,1,.false.)
!!%-----------------------------------------------------------
!return
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine get_name_sps &
!>  (this, &
!>   name)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Return the name of the species. 
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(in)   :: this  !> Type species variable. 
!
!character(len=*), intent(out) :: name  !> Name of the species. 
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!> 
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!> 
!
!!%-------------------------------------------------------------
!name=this%name
!!%-------------------------------------------------------------
!return
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine get_prop_sps &
!>  (this, &
!>   value, &
!>   nameprop, &
!>   msg, &
!>   iserror)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Return the value of some species property from its name.
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(in)  :: this     !> Type species variable. 
!
!real*8, intent(out)          :: value    !> Value of the species property. 
!
!character(len=*), intent(in) :: nameprop !> Name of the species property. 
!
!logical, intent(out)         :: iserror  !> iserror=true, there was an error. 
!
!character(len=*), intent(out):: msg      !> Error message 
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!integer          :: &
!> i
!logical          :: &
!> isbe
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!msg=''
!iserror=.false.
!!%-----------------------------------------------------------
!!% Initialice variables
!!%-----------------------------------------------------------
!isbe=.false.
!value=0.0d0
!!%-----------------------------------------------------------
!!% Find the species properties 
!!%-----------------------------------------------------------
!do i=1,this%numprop
!> 
!>  if (nameprop==this%nameprop(i)) then
!>   value=this%prop(i)
!>   isbe=.true. 
!>   exit
!>  end if
!> 
!end do
!!%------------------------------------------------------------
!!% Check if the species properties was found
!!%------------------------------------------------------------ 
!if (.not.isbe) then
!>  msg='Error, species property not found:'
!>  call add_ (msg,nameprop)
!>  goto 10
!end if
!!%------------------------------------------------------------
!return
!> 
!10 continue 
!print *,'********************'
!print *,'Species:'
!print *,'Name:', this%name
!print *,'Service: get_prop_'
!print *, msg
!print *,'********************'
!iserror=.true.
!return
!> 
!> 
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine set_prop_sps &
!>  (this, &
!>   prop, &
!>   nameprop, &
!>   unitprop, &
!>   nprop)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Set the name and values of the species properties.
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(inout)                  :: this     !> Type species variable
!
!integer, intent(in)                             :: nprop    !> Number of properties
!
!real*8, intent(in), dimension(nprop)            :: prop     !> Values of the species properties 
!
!character(len=*), intent(in), dimension(nprop)  :: nameprop !> Name of the species properties
!
!character(len=*), intent(in), dimension(nprop)  :: unitprop !> Name of the species properties
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!this%numprop = nprop
!!%------------------------------------------------------------
!!% Allocate pointers 
!!%------------------------------------------------------------
!call check_pointer_ (this%prop,this%numprop,.true.)
!call check_pointer_ (this%nameprop,this%numprop,.true.)
!call check_pointer_ (this%unitprop,this%numprop,.true.)
!!%------------------------------------------------------------
!!% Set the name and value of species properties 
!!%------------------------------------------------------------
!this%prop=prop
!this%nameprop=nameprop
!this%unitprop=unitprop
!!%------------------------------------------------------------
!return
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine add_prop_sps &
!>  (this, &
!>   prop, &
!>   nameprop, &
!>   unitprop, &
!>   iserror)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Add species property.
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(inout)                  :: this     !> Type species variable
!
!real*8, intent(in)                              :: prop     !> Values of the species property 
!
!character(len=*), intent(in)                    :: nameprop !> Name of the species property
!
!character(len=*), intent(in)                    :: unitprop !> Name of the species property
!
!logical, intent(out)                            :: iserror  !> If iserror=.true., then there was an error
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!real*8, pointer             :: &
!> proploc(:) => null ()
!character(len=100), pointer :: & 
!> nameproploc(:) => null ()
!character(len=100), pointer :: & 
!> unitproploc(:) => null ()
!integer                     :: &
!> i 
!!-------------------------------------------------------------------------
!!
!!>   $code
!
!
!iserror=.false. 
!!%------------------------------------------------------------
!!% Check if the species properties was previously defined 
!!%------------------------------------------------------------
!do i=1,this%numprop
!> if (this%nameprop(i)==nameprop) then
!>    this%prop(i)=prop
!>    return  
!> end if
!end do
!!%------------------------------------------------------------
!!% Allocate pointers 
!!%------------------------------------------------------------
!call check_pointer_ (proploc,this%numprop,.true.)
!call check_pointer_ (nameproploc,this%numprop,.true.)
!call check_pointer_ (unitproploc,this%numprop,.true.)
!nameproploc=this%nameprop
!unitproploc=this%unitprop
!proploc=this%prop
!this%numprop=this%numprop+1
!call check_pointer_ (this%prop,this%numprop,.true.)
!call check_pointer_ (this%nameprop,this%numprop,.true.)
!call check_pointer_ (this%unitprop,this%numprop,.true.)
!this%nameprop(1:this%numprop-1)=nameproploc
!this%unitprop(1:this%numprop-1)=unitproploc
!this%prop(1:this%numprop-1)=proploc
!!%------------------------------------------------------------
!!% Add the new property 
!!%------------------------------------------------------------
!this%nameprop(this%numprop)=nameprop
!this%unitprop(this%numprop)=unitprop
!this%prop(this%numprop)=prop
!20 continue 
!call check_pointer_ (proploc,1,.false.)
!call check_pointer_ (nameproploc,1,.false.)
!call check_pointer_ (unitproploc,1,.false.)
!!%------------------------------------------------------------
!return
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine set_name_sps &
!>  (this, &
!>   name)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Set the name of the species object. 
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(inout)   :: this !> Type species object. 
!
!character(len=*), intent(in)     :: name !> Name of the species object 
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!> 
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!> 
!
!!%------------------------------------------------------------
!this%name = name
!!%------------------------------------------------------------
!return
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine write_sps &
!>  (this, &
!>   ioutput, &
!>   iserror)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Write on the ioutput unit all attributtes of the species object. 
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(in)  :: this    !> Type species variable.
!
!integer, intent(in)          :: ioutput !> Output unit
!
!logical, intent(out)         :: iserror !> iserror=true, there was an error. 
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!integer        :: &
!> i 
!character(len=100)         :: &
!> msg 
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!> 
!!%------------------------------------------------------------
!iserror=.false.
!msg=''
!!%------------------------------------------------------------
!write (ioutput,*) "------------------------------------------"
!write (ioutput,*) "               Species object             "
!write (ioutput,*) "------------------------------------------"
!write (ioutput,1) this%name
!write (ioutput,*) "------------------------------------------"
!do i=1,this%numprop
!>  write (ioutput,3) this%nameprop(i),"=",this%prop(i),this%unitprop(i)
!end do
!write (ioutput,*) "------------------------------------------"
!!%------------------------------------------------------------
!return
!1 format (a20)
!2 format (a8,i5)
!3 format(a20,a1,f10.2,1x,a20)
!10 continue 
!print *,'********************'
!print *,'Species:'
!print *,'Name:', this%name
!print *,'Service: write_'
!print *, msg
!print *,'********************'
!iserror=.true.
!return
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine copy_sps &
!>  (targetobj, &
!>   sourceobj)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description: Copy a spcecies object in other species object.  
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(in)  :: sourceobj   !> Type species variable
!
!type(t_species), intent(out) :: targetobj   !> Type species variable
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond: The species must be previously created.
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!integer        :: &
!> i 
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!> 
!!%-----------------------------------------------------------
!targetobj%name = sourceobj%name
!targetobj%numprop = sourceobj%numprop
!!%-----------------------------------------------------------
!if (targetobj%numprop>0) then
!> call check_pointer_ (targetobj%prop,targetobj%numprop,.true.)
!> call check_pointer_ (targetobj%nameprop,targetobj%numprop,.true.)
!> call check_pointer_ (targetobj%unitprop,targetobj%numprop,.true.)
!> targetobj%prop=sourceobj%prop
!> targetobj%nameprop=sourceobj%nameprop
!> targetobj%unitprop=sourceobj%unitprop
!end if
!!%-----------------------------------------------------------
!return
!end subroutine
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%***************Private subroutines**************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine begin_element_handler (name, attributes)
!> 
!character(len=*), intent(in)   :: name
!type(dictionary_t), intent(in) :: attributes
!integer                        :: status
!> 
!> 
!> call read_xml_loc_ (name, attributes)
!> 
!> 
!return
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine read_xml_loc_sps (name, attributes,this,msg,iserror)
!implicit none 
!type (t_species), optional :: this
!character(len=*), intent(in)   :: name
!type(dictionary_t), intent(in) :: attributes
!logical, optional              :: iserror
!character(len=*), optional     :: msg
!> 
!integer        :: status,n
!real*8         :: reallocal(1)
!logical        :: havethis
!character(len=100)      ::id
!integer, save           ::numprop,iprop
!real*8, save            ::value(20)
!character(len=100), save::namesp, phasetype, nameprop(20),unitprop(20),empty
!!%---------------------------------------------------------------
!havethis=present(this)
!!%---------------------------------------------------------------
!if(havethis) then
!!%------------------------------
!iserror=.false.
!msg=''
!!%------------------------------
!empty=''
!!%------------------------------
!> this%numprop=numprop
!> call check_pointer_ (this%prop,numprop,.true.)
!> call check_pointer_ (this%nameprop,numprop,.true.)
!> call check_pointer_ (this%unitprop,numprop,.true.)
!> this%prop=value(1:numprop)
!> this%nameprop=nameprop(1:numprop)
!> this%unitprop=unitprop(1:numprop)
!> if (namesp==empty) then
!>  msg='Error, not defined the name of species'
!>  goto 10
!> end if
!> this%name=namesp
!> 
!> 
!else
!> 
!> call lowercase (name)
!> 
!> select case(name)
!> 
!> 
!> case ('species')
!>   iprop=0
!>   id=' '
!>   call get_value (attributes,"name", id, status)
!>   namesp=id
!>   id=' '
!>   call get_value (attributes,"phasetype", id, status)
!>   phasetype=id
!> 
!> case ('property')
!> 
!>   iprop=iprop+1
!>   numprop=iprop
!>   id=' '
!>   call get_value (attributes,"name", id, status)
!>   nameprop(iprop)=id
!>   id=' '
!>   call get_value (attributes,"unit", id, status)
!>   unitprop(iprop)=id
!>   id=' '
!>   call get_value (attributes,"value", id, status)
!>   n=0
!>   call build_data_array (id,reallocal,n)
!>   value(iprop)=reallocal(1)
!> 
!> end select
!> 
!end if
!!%-------------------------------------------------------------
!return
!> 
!10 iserror=.true.
!return
!> 
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine verify1_sps &
!>  (this, &
!>   name, &
!>   propsp, &
!>   nameprop, &
!>   numprop, &
!>   nerror, &
!>   errormsg)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description:  Verify different attributtes encapsulated in the species object. 
!!
!!>   $Arguments:
!!
!> 
!type(t_species), intent(in)                    :: this
!
!integer, intent(out)                           :: nerror
!
!integer, optional                              :: numprop
!
!real*8, pointer, optional                      :: propsp(:)
!
!character(len=*), intent(in), optional         :: name
!
!character(len=100), pointer, optional          :: nameprop(:)
!
!character(len=*), pointer                      :: errormsg(:) 
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!logical :: &
!> havename, &
!> havepropsp, &
!> havenameprop, &
!> havenumprop
!character(len=300):: &
!> errormsgloc(50)
!integer           :: &
!> i 
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!> 
!
!> 
!if(associated(errormsg)) deallocate (errormsg)
!> 
!> havename=present(name)
!> havepropsp=present(propsp)
!> havenameprop=present(nameprop)
!> havenumprop=present(numprop)
!> 
!> 
!nerror=0
!errormsgloc=' '
!!%--------------------------------------------------
!if (havename) then
!>  if (this%name.ne.name) then
!>   nerror=nerror+1
!>   errormsgloc(nerror)='error en name'
!>  end if
!end if
!!%--------------------------------------------------
!if (havenumprop) then
!>  if (this%numprop.ne.numprop) then
!>   nerror=nerror+1
!>   errormsgloc(nerror)='error en numprop'
!>  end if
!end if
!!%--------------------------------------------------
!if (havenameprop) then
!select case (associated(nameprop))
!> 
!case (.true.)
!> if (.not.associated(this%nameprop)) then
!>  nerror=nerror+1
!>  errormsgloc(nerror)='error en nameprop'
!> else
!>  if(size(nameprop).eq.this%numprop) then
!>  do i=1,this%numprop
!>   if (nameprop(i).ne.this%nameprop(i)) then
!>    nerror=nerror+1
!>    errormsgloc(nerror)='error en nameprop'
!>   end if
!>  end do
!>  else
!>    nerror=nerror+1
!>    errormsgloc(nerror)='error en nameprop'
!>  end if
!> end if
!case default
!> if (associated(this%nameprop)) then
!>  nerror=nerror+1
!>  errormsgloc(nerror)='error en nameprop'
!> end if
!end select
!end if
!!%--------------------------------------------------
!if (havepropsp) then
!select case (associated(propsp))
!> 
!case (.true.)
!> if (.not.associated(this%prop)) then
!>  nerror=nerror+1
!>  errormsgloc(nerror)='error en propsp'
!> else
!>  if(size(propsp).eq.this%numprop) then
!>   do i=1,this%numprop
!>    if (propsp(i).ne.this%prop(i)) then
!>     nerror=nerror+1
!>     errormsgloc(nerror)='error en propsp'
!>    end if
!>   end do
!>  else
!>    nerror=nerror+1
!>    errormsgloc(nerror)='error en propsp'
!>  end if
!>  end if
!case default
!> if (associated(this%prop)) then
!>  nerror=nerror+1
!>  errormsgloc(nerror)='error en propsp'
!> end if
!end select
!end if
!!%--------------------------------------------------
!if (nerror.ne.0) then
!> allocate (errormsg(nerror))
!> errormsg=errormsgloc(1:nerror)
!end if
!> 
!return
!end subroutine
!!%************************************************************
!!%***************Public subroutine****************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!subroutine verify2_sps &
!>  (this, &
!>   species, &
!>   nerror, &
!>   errormsg)
!> 
!implicit none
!!-------------------------------------------------------------------------
!!
!!>   $Description:  Verify different attributtes encapsulated in the species object. 
!!
!!>   $Arguments:
!!
!> 
!type(t_species)                        :: this
!
!type(t_species), intent(in)            :: species
!
!integer                                :: nerror
!
!character(len=*), pointer              :: errormsg(:) 
!> 
!!-------------------------------------------------------------------------
!!
!!>   $Pre-cond:
!!
!!>   $Post-cond:
!!
!!>   $License:
!!
!!-------------------------------------------------------------------------
!> 
!!-------------------------------------------------------------------------
!!
!!>   $code
!!
!> 
!> 
!
!!%------------------------------------------------------------
!call verify_ &
!> (this, &
!>  name=species%name, &
!>  propsp=species%prop, &
!>  nameprop=species%nameprop, &
!>  numprop=species%numprop, &
!>  nerror=nerror, &
!>  errormsg=errormsg)
!!%------------------------------------------------------------
!return
!end subroutine
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!!%************************************************************
!end module m_species