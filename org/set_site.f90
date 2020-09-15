subroutine set_site(fem)
  use module_fem
  implicit none
  type(fem_set) :: fem
  real(8) :: sx, sz, xi, zeta
  real(8), allocatable :: n(:,:)
  integer :: dof, style
  integer :: ind, i, is, ie

  dof = fem%dof

  write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++"

  do is = 1, fem%nsite
     sx = fem%site(is)%x
     sz = fem%site(is)%z

     do ie = 1, fem%nelem
        call check_inside(sx,sz,fem%element(ie),xi,zeta,ind)
        if(ind == 1) then

           style = fem%element(ie)%style
           fem%site(is)%style = style
           allocate(fem%site(is)%node(style))

           do i = 1, style
              fem%site(is)%node(i)%p => fem%element(ie)%node(i)%p              
           end do

           allocate(n(dof,style*dof))
           call mk_n(dof,style,xi,zeta,n)

           allocate(fem%site(is)%coeff(dof,style*dof))
           fem%site(is)%coeff = n 

           allocate(fem%site(is)%disp(dof),fem%site(is)%dispm(dof))
           allocate(fem%site(is)%vel(dof))

           fem%site(is)%disp = 0.d0
           fem%site(is)%dispm = 0.d0
           fem%site(is)%vel = 0.d0

           write(*,*) " Site element No.", ie

           deallocate(n)
           exit
        end if
     end do
  end do


end subroutine set_site
