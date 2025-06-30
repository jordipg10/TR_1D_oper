!subroutine get_conc_kin_Monod(this,species,conc,conc_kin,kin_ind)
!    use Monod_m
!    use aq_species_m
!    implicit none
!    class(redox_kin_c), intent(in) :: this
!    class(aq_species_c), intent(in) :: species(:) !> aqueous species (chapuza)
!    real(kind=8), intent(in) :: conc(:) ! aqueous species concentrations
!    real(kind=8), intent(out) :: conc_kin(:) ! concentrations relevant to compute kinetic reaction rate
!    integer(kind=4), intent(out), optional :: kin_ind(:) ! indices of conc_kin
!            
!    integer(kind=4) :: i,j,k,l,m,p,n,DOC_ind
!    real(kind=8), parameter :: epsilon=1d-6
!    type(species_c) :: DOC ! CH2O
!    logical :: DOC_flag
!            
!    
!    if (size(conc)/=size(species)) error stop "Dimension error in get_conc_kin_Monod"
!                        
!    i=1
!    j=1
!    k=1
!    if (this%params%n_inh>0) then
!        do
!            if (species(i)%name==this%params%inhibitors(j)%name) then
!                conc_kin(k)=conc(i)
!                if (present(kin_ind)) then
!                    kin_ind(k)=i
!                end if
!                if (j<this%params%n_inh) then
!                    j=j+1
!                    i=1
!                end if
!                k=k+1
!            end if
!            if (i<size(species)) then
!                i=i+1
!            else
!                exit
!            end if
!        end do
!    end if
!    i=1
!    j=1
!    k=1
!    l=1
!    m=1
!    do
!        if (species(i)%name==this%params%TEAs(l)%name) then
!            conc_kin(this%params%n_inh+m)=conc(i) 
!            if (present(kin_ind)) then
!                kin_ind(this%params%n_inh+m)=i
!            end if
!            if (l<this%params%n_M) then
!                l=l+1
!                i=1
!            end if
!            m=m+1
!        !else if (species(i)%name=='ch2o' .and. species(i)%name/=this%params%TEAs(2)%name) then ! autentica chapuza
!        !    conc_kin(this%params%n_t+1)=conc(i)
!        !    if (present(kin_ind)) then
!        !        kin_ind(this%params%n_t+1)=i
!        !    end if
!        else if (species(i)%name=='ch2o(aq)') then
!            conc_kin(this%params%n_t+1)=conc(i)
!            if (present(kin_ind)) then
!                kin_ind(this%params%n_t+1)=i
!            end if
!        else
!            continue
!        end if
!        if (i<size(species)) then
!            i=i+1
!        else
!            exit
!        end if
!    end do
!end subroutine