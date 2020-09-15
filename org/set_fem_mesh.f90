subroutine set_fem_mesh(fem)
  use module_fem
  implicit none
  type(fem_set) :: fem
  type(node_set), pointer :: node
  integer :: inode, style
  integer :: i, j
  
  do i = 1, fem%nelem
     style = fem%element(i)%style
     allocate(fem%element(i)%node(style))

     do j = 1, style
        inode = fem%element(i)%inode(j)
        call search_node(fem%nodelist,inode,fem%element(i)%node(j)%p)
     end do
  end do

end subroutine set_fem_mesh
