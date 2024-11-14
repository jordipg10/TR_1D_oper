subroutine Thomas_induction(A,b,i,x_i,tol,x_k)    
    !> A: tridiagonal matrix
    !> b: independent term
    !> x: solution of linear system

    use matrices_m
    implicit none
    class(tridiag_matrix_c), intent(in) :: A
    real(kind=8), intent(in) :: b(:)
    integer(kind=4), intent(in) :: i !> row index
    real(kind=8), intent(in) :: x_i
    real(kind=8), intent(in) :: tol !> tolerance
    real(kind=8), intent(out) :: x_k(:) !> must be allocated
    
    integer(kind=4) :: k
    real(kind=8), allocatable :: c_star(:),d_star(:), matrix(:,:)
    
    !allocate(c_star(this%dim-1),d_star(this%dim))
    !c_star(1)=A%super(1)/A%diag(1)
    !d_star(1)=b(1)/A%diag(1)
    !do i=2,this%dim-1
    !>    c_star(i)=A%super(i)/(A%diag(i)-A%sub(i-1)*c_star(i-1))
    !>    d_star(i)=(b(i)-A%sub(i-1)*d_star(i-1))/(A%diag(i)-A%sub(i-1)*c_star(i-1))
    !end do
    !d_star(this%dim)=(b(this%dim)-A%sub(i-1)*d_star(this%dim-1))/(A%diag(i)-A%sub(i-1)*c_star(i-1))
    !x(this%dim)=d_star(this%dim)
    !x(this%dim-1)=d_star(this%dim-1)-c_star(this%dim-1)*x(this%dim)
    !do k=2,this%dim-1
    !>    sum=0d0
    !>    do j=0,k-1
    !>        prod=1d0
    !>        do l=1,k-j
    !>            prod=prod*c_star(i-j-l)
    !>        end do
    !>        sum=sum+prod*d_star(i-j)
    !>    x(this%dim-k)=d_star(this%dim-k)+sum()
    !end do
    !deallocate(c_star,d_star)
end subroutine