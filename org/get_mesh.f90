subroutine get_mesh(fem)
  use module_fem
  implicit none
  type(fem_set) :: fem
  type(node_param_set) :: node_param
  integer :: i

  open(11, file = "mesh/mesh.in")
  read(11,*) fem%nnode, fem%nelem, fem%dof
  write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++"
  write(*,*) " Read nodes   ", fem%nnode
  write(*,*) " Read elements", fem%nelem

  allocate(node_param%freedom(fem%dof))
  node_param%dof = fem%dof

  call init_nodelist(fem%nodelist)
  do i = 1, fem%nnode
     read(11,*) node_param%inode, node_param%x, node_param%z, node_param%freedom(1:fem%dof)
     call append_node(fem%nodelist,node_param)
  end do
!  call output_nodelist(fem%nodelist)
  
  allocate(fem%element(fem%nelem))

  do i = 1, fem%nelem
     read(11,*) fem%element(i)%ielem, fem%element(i)%style, fem%element(i)%vs, fem%element(i)%vp, fem%element(i)%rho

     fem%element(i)%dof =  fem%dof
     allocate(fem%element(i)%inode(fem%element(i)%style))
     read(11,*) fem%element(i)%inode(1:fem%element(i)%style), fem%element(i)%damping(1:4)
  end do

  close(11)

end subroutine get_mesh
