subroutine update_time(fem,dt,slip)
  use module_fem
  implicit none
  type(fem_set) :: fem
  real(8) :: dt, slip

  call update_time_core(fem,dt,slip)
  call update_site(fem,dt)

end subroutine update_time

!----------------------------------------------------------------!
subroutine update_time_core(fem,dt,slip)
  use module_fem
  implicit none
  type(fem_set) :: fem
  type(node_set), pointer :: node
  real(8) :: u, um
  real(8) :: dt, slip
  integer :: dof, style
  integer :: is, ie, j

  dof = fem%dof

  node => fem%nodelist%head
  do while(associated(node%next_node))
     node => node%next_node
     node%node_param%force(1:dof) = 0.d0
  end do

  do is = 1, fem%nsrc
     call mk_source(fem%src(is),dof,slip)
  end do

  do ie = 1, fem%nelem
     style = fem%element(ie)%style

     call mk_ku(fem%element(ie),dof,style)
     call mk_cv(fem%element(ie),dof,style)
  end do

  node => fem%nodelist%head
  do while(associated(node%next_node))
     node => node%next_node

     do j = 1, dof
        if(node%node_param%freedom(j) == 0) then
           node%node_param%u(j) = 0.d0
           node%node_param%um(j) = 0.d0
           node%node_param%v(j) = 0.d0

        else
           um = node%node_param%um(j)
           u = node%node_param%u(j)
           
           node%node_param%u(j) = (node%node_param%mass(j)*(2.d0*u - um) &
                + 0.5d0*dt*node%node_param%c(j)*um - dt*dt*node%node_param%force(j)) &
                * node%node_param%inv_mc(j)

           node%node_param%v(j) = (node%node_param%u(j) - um) / (2.d0*dt)
           node%node_param%um(j) = u
        end if
     end do          
  end do

end subroutine update_time_core

!----------------------------------------------------------------!
subroutine mk_source(src,dof,slip)
  use module_fem
  implicit none
  type(src_set) :: src
  type(node_set), pointer :: node
  real(8) :: slip
  integer :: dof
  integer :: in, i, j

  in = 1
  do i = 1, src%style
     node => src%node(i)%p
     do j = 1, dof
        node%node_param%force(j) = node%node_param%force(j) + src%coeff(in)*slip
        in = in + 1
     end do
!     write(*,*) node%node_param%force(1:3)
  end do

end subroutine mk_source

!----------------------------------------------------------------!
subroutine mk_ku(element,dof,style)
  use module_fem
  implicit none
  type(element_set) :: element
  type(node_set), pointer :: node
  integer :: dof, style
  real(8) :: u(dof*style), ku(dof*style)
  integer :: i, j, in

  in = 1
  do i = 1, style
     node => element%node(i)%p
     do j = 1, dof
        u(in) = node%node_param%u(j)
        in = in + 1
     end do
  end do

  ku = matmul(element%k,u)

  in = 1
  do i = 1, style
     node => element%node(i)%p
     do j = 1, dof
        node%node_param%force(j) = node%node_param%force(j) + ku(in)
        in = in + 1
     end do
  end do

end subroutine mk_ku

!----------------------------------------------------------------!
subroutine mk_cv(element,dof,style)
  use module_fem
  implicit none
  type(element_set) :: element
  type(node_set), pointer :: node
  integer :: dof, style
  real(8) :: v(dof*style), cv(dof*style)
  integer :: i, j, in

  in = 1
  do i = 1, style
     node => element%node(i)%p
     do j = 1, dof
        v(in) = node%node_param%v(j)
        in = in + 1
     end do
  end do

  do i = 1, style*dof
     cv(i) = 0.d0
     do j = 1, style*dof
        if(i == j) cycle
        cv(i) = cv(i) + element%c(i,j)*v(j)
     end do
  end do

  in = 1
  do i = 1, style
     node => element%node(i)%p
     do j = 1, dof
        node%node_param%force(j) = node%node_param%force(j) + cv(in)
        in = in + 1
     end do
  end do

end subroutine mk_cv

!----------------------------------------------------------------!
subroutine update_site(fem,dt)
  use module_fem
  implicit none
  type(fem_set) :: fem
  type(node_set), pointer :: node
  real(8) :: dt
  real(8), allocatable :: u(:), v(:)
  integer :: dof, style
  integer :: is, in, i, j

  dof = fem%dof

  do is = 1, fem%nsite
     style = fem%site(is)%style
     allocate(u(dof*style))
     allocate(v(dof*style))

     in = 1
     do i = 1, style
        node => fem%site(is)%node(i)%p
        do j = 1, dof
           u(in) = node%node_param%u(j)
           v(in) = node%node_param%v(j)
           in = in + 1
        end do
     end do

     fem%site(is)%dispm = fem%site(is)%disp
     fem%site(is)%disp = matmul(fem%site(is)%coeff,u)
     fem%site(is)%vel = matmul(fem%site(is)%coeff,v)

     deallocate(u,v)
  end do

end subroutine update_site
