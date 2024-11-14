!subroutine set_rk(this,rk_type,n)
!! rk_type=1 for linear rk
!!         2 for affine rk
!!         3 for Monod
!!         4 for mineral
!    !use rk_aff_m
!    use RT_1D_m
!    use Monod_m
!    use rk_mineral_m
!    implicit none
!    class(chemistry_c) :: this
!    integer(kind=4), intent(in) :: rk_type
!    integer(kind=4), intent(in) :: n
!    
!    real(kind=8), allocatable :: lambda(:),conc_eq(:,:)
!    class(rk_lin_c), allocatable :: rk_lin(:,:)
!    class(rk_aff_c), allocatable :: rk_aff(:,:)
!    class(rk_redox_kin_c), allocatable :: rk_Monod(:,:)
!    integer(kind=4) :: i,j,nk
!    
!    nk=this%chem_syst%get_num_kin_reactions()
!    
!    if (rk_type==0) then
!        continue
!    else if (rk_type<3) then
!        allocate(rk_lin(nk,n))
!        allocate(lambda(nk))
!        lambda=1d3
!        do i=1,nk
!            do j=1,n
!                call rk_lin(i,j)%set_lambda(lambda(i))
!                !if (j==3) then
!                !    call rk_lin(i,j)%set_lambda(lambda(i))
!                !else
!                !    call rk_lin(i,j)%set_lambda(0d0)
!                !end if
!            end do
!        end do 
!        if (rk_type==2) then
!            allocate(rk_aff(nk,n))
!            allocate(conc_eq(nk,n))
!            conc_eq=0d0
!            do i=1,nk
!                do j=1,n
!                    rk_aff(i,j)%lambda=rk_lin(i,j)%lambda
!                    call rk_aff(i,j)%set_conc_eq(conc_eq(i,j))
!                end do
!            end do
!            this%rk=rk_aff
!        else
!            this%rk=rk_lin
!        end if
!    else if (rk_type==3) then
!        allocate(rk_Monod(nk,n))
!        do i=1,nk
!            do j=1,n
!                call rk_Monod(i,j)%read_rk('monod_O2.dat')
!            end do
!        end do
!        this%rk=rk_Monod
!    end if
!    
!end subroutine