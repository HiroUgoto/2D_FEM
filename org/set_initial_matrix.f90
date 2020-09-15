subroutine set_initial_matrix(fem)
  use module_fem
  implicit none
  type(fem_set) :: fem
  type(node_set), pointer :: node
  integer :: ie, in, id, i

  do ie = 1, fem%nelem
     call mk_local_stiffness_matrix(fem%element(ie))
     call mk_local_mass_matrix(fem%element(ie))
     call mk_local_damping_matrix(fem%element(ie))

     in = 1
     do i = 1, fem%element(ie)%style
        node => fem%element(ie)%node(i)%p
        do id = 1, fem%dof
           node%node_param%mass(id) = node%node_param%mass(id) + fem%element(ie)%m(in)
           node%node_param%c(id) = node%node_param%c(id) + fem%element(ie)%c(in,in)
           in = in + 1
        end do
     end do
  end do

  node => fem%nodelist%head
  do while(associated(node%next_node))
     node => node%next_node
     
     do id = 1, fem%dof
        node%node_param%inv_mc(id) = 1.d0/(node%node_param%mass(id) + 0.5d0*fem%dt*node%node_param%c(id))
     end do
  end do

end subroutine set_initial_matrix
