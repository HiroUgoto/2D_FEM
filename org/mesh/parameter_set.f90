module parameter_set
  implicit none

  ! Node object !
  type node_set
     integer :: dof
     real(8) :: x, z
     integer, allocatable :: freedom(:)
  end type node_set

  ! Element object !
  type element_set
     integer :: style
     integer :: inode(9)
     real(8) :: vs, vp, rho
     integer :: damping(4) 
  end type element_set

end module parameter_set
