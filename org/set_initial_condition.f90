subroutine set_initial_condition(fem)
  use module_fem
  implicit none
  type(fem_set) :: fem
  type(node_set), pointer :: node
  integer :: dof

  node => fem%nodelist%head
  do while(associated(node%next_node))
     node => node%next_node

     dof = node%node_param%dof

     allocate(node%node_param%u(dof))
     allocate(node%node_param%um(dof))
     node%node_param%u = 0.d0
     node%node_param%um = 0.d0

     allocate(node%node_param%v(dof))
     node%node_param%v = 0.d0

     allocate(node%node_param%mass(dof))
     node%node_param%mass = 0.d0

     allocate(node%node_param%c(dof))
     node%node_param%c = 0.d0

     allocate(node%node_param%inv_mc(dof))
     node%node_param%inv_mc = 0.d0

     allocate(node%node_param%force(dof))
     node%node_param%force = 0.d0
  end do

end subroutine set_initial_condition
