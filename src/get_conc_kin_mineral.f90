subroutine get_conc_kin_mineral(this,species,conc,conc_kin,kin_ind)
    use kin_mineral_m
    implicit none
    class(kin_mineral_c), intent(in) :: this
    class(species_c), intent(in) :: species(:)
    real(kind=8), intent(in) :: conc(:) ! aqueous species concentrations
    real(kind=8), intent(out) :: conc_kin(:) ! concentrations relevant to compute kinetic reaction rate
    integer(kind=4), intent(out), optional :: kin_ind(:) ! indices of conc_kin
            
    integer(kind=4) :: j,l,m,p,n,DOC_ind
    real(kind=8), parameter :: epsilon=1d-6
    type(species_c) :: DOC ! CH2O
    logical :: DOC_flag
            
    l=1
    m=1
    n=1
    if (size(conc)/=size(species)) error stop "Dimension error in get_conc_kin_mineral"
                    
   error stop "Subroutine 'get_conc_kin_mineral' not implemented yet"
end subroutine