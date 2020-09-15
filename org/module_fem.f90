module module_fem
  use module_node
  implicit none
  real(8), parameter :: pi = 3.14159265358979323846264d0


  type element_set
     integer :: ielem
     integer :: dof, style
     real(8) :: vs, vp, rho
     integer, allocatable :: inode(:)
     integer :: damping(4)
     type(node_pointer_set), allocatable :: node(:)
     real(8), allocatable :: m(:), k(:,:), c(:,:)
  end type element_set

  type src_set 
     real(8) :: x, z
     real(8) :: ds, tp, fp, dip, rake
     integer :: style
     type(node_pointer_set), allocatable :: node(:)
     real(8), allocatable :: coeff(:)
  end type src_set

  type site_set 
     integer :: style
     real(8) :: x, z
     type(node_pointer_set), allocatable :: node(:)
     real(8), allocatable :: disp(:), dispm(:), vel(:)
     real(8), allocatable :: coeff(:,:)
  end type site_set

  type fem_set
     type(nodelist_set) :: nodelist
     type(element_set), allocatable :: element(:)
     type(src_set), allocatable :: src(:)
     type(site_set), allocatable :: site(:)
     integer :: dof
     integer :: ntim, nnode, nelem
     integer :: nsrc, nsite
     real(8) :: dt, tim
  end type fem_set

end module module_fem
