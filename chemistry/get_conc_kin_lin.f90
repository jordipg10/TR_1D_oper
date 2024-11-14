subroutine get_conc_kin_lin(this,conc,conc_kin,kin_ind)
    use lin_kin_reaction_m
    implicit none
    class(lin_kin_reaction_c), intent(in) :: this
    real(kind=8), intent(in) :: conc(:) ! aqueous species concentrations
    real(kind=8), intent(out) :: conc_kin(:) ! concentrations relevant to compute kinetic reaction rate
    integer(kind=4), intent(out), optional :: kin_ind(:) ! indices of conc_kin
            
    integer(kind=4) :: i,j,k,l,m,p,n,DOC_ind
    real(kind=8), parameter :: epsilon=1d-6
    type(species_c) :: DOC ! CH2O
    logical :: DOC_flag
    
    !if (size(conc)/=this%num_species) error stop "Dimension error in get_conc_kin_lin"
    
    conc_kin=0d0
    kin_ind=0
    
    conc_kin(1)=conc(1)
    kin_ind(1)=1
end subroutine