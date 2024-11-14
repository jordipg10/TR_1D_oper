!subroutine compute_conc_solid_species(this,k) !> F*dc/dt=S_c^T*r
!>    use RT_1D_m
!>    implicit none
!>    class(RT_1D_c) :: this
!>    integer(kind=4), intent(in), optional :: k
!>            
!>    integer(kind=4) :: i,j
!>    real(kind=8), allocatable :: rk_mat(:,:)
!>            
!>    select type (this)
!>    type is (RT_1D_transient_c)
!>        select type (transport=>this%transport)
!>        class is (diffusion_transient_c)
!>            select type (chem_syst=>this%chem_syst)
!>            type is (chem_system_eq_c)
!>                !allocate(rk_mat(chem_syst%num_reactions,transport%spatial_discr%Num_targets))
!>                !rk_mat=1d0
!>                !rk_mat(:,1)=1d0
!>                if (.not. present(k)) then
!>                    do i=1,chem_syst%num_solid_species
!>                        !allocate(this%solid_zone%minerals_eq(i)%conc_solid(transport%spatial_discr%Num_targets))
!>                        allocate(this%solid_zone%solid_species(i)%conc_solid(transport%spatial_discr%Num_targets))
!>                        do j=1,transport%spatial_discr%Num_targets
!>                            !this%solid_zone%minerals_eq(i)%conc_solid(j)=this%solid_zone_init%minerals_eq(i)%conc_solid(j)+transport%time_discr%Final_time*dot_product(chem_syst%Sc(:,i),rk_mat(:,j))/transport%F_mat%diag(j)
!>                            this%solid_zone%solid_species(i)%conc_solid(j)=this%solid_zone_init%solid_species(i)%conc_solid(j)+transport%time_discr%Final_time*dot_product(chem_syst%Sc(:,i),this%r_eq(:,j)%reaction_rate)/transport%F_mat%diag(j)
!>                        end do
!>                        !print *, this%solid_zone%solid_species(i)%conc_solid
!>                    end do
!>                end if
!>            end select
!>        end select
!>    end select
!end subroutine