module module_node
  implicit none

  type node_param_set
     integer :: dof
     integer :: inode
     real(8) :: x, z
     real(8), allocatable :: u(:), um(:)
     real(8), allocatable :: v(:)
     real(8), allocatable :: mass(:)
     real(8), allocatable :: inv_mc(:), c(:)
     real(8), allocatable :: force(:)
     integer, allocatable :: freedom(:)
  end type node_param_set

  type node_set
     type(node_param_set) :: node_param
     type(node_set), pointer :: next_node
  end type node_set

  type nodelist_set
     type(node_set), pointer :: head, end_node
     integer :: nnode
  end type nodelist_set

  type node_pointer_set
     type(node_set), pointer :: p
  end type node_pointer_set

!------------------------------------------------!
! # nodelist  
!     head => 1st node => .... => end_node       !
! 
!------------------------------------------------!
contains

  !----------------------------------------------!
  subroutine init_nodelist(nodelist)
    type(nodelist_set) :: nodelist

    nodelist%nnode = 0
    
    allocate(nodelist%head)
    allocate(nodelist%head%next_node)
    nullify(nodelist%head%next_node%next_node)

    nodelist%end_node => nodelist%head

  end subroutine init_nodelist

  !----------------------------------------------!
  subroutine erase_nodelist(nodelist)
    type(nodelist_set) :: nodelist
    type(node_set), pointer :: node, tmp_node

    if(associated(nodelist%head)) then
       node => nodelist%head
       do while(associated(node%next_node))
          tmp_node => node
          node => node%next_node
          deallocate(tmp_node)
       end do
       deallocate(node)
       
       nodelist%nnode = 0
    end if

  end subroutine erase_nodelist

  !----------------------------------------------!
  subroutine output_nodelist(nodelist)
    type(nodelist_set) :: nodelist
    type(node_set), pointer :: node

    node => nodelist%head
    do while(associated(node%next_node))
       node => node%next_node
       write(*,*) node%node_param%inode, node%node_param%x, node%node_param%freedom(1)
!       write(*,*) node%node_param%inode, node%node_param%mass(1)
    end do

  end subroutine output_nodelist

  !----------------------------------------------!
  subroutine append_node(nodelist,node_param)
    type(nodelist_set) :: nodelist
    type(node_param_set) :: node_param
    type(node_set), pointer :: new_node

    allocate(nodelist%end_node%next_node)
    new_node => nodelist%end_node%next_node

    new_node%node_param = node_param
    nullify(new_node%next_node)

    nodelist%end_node => new_node
    nodelist%nnode = nodelist%nnode + 1

  end subroutine append_node

  !----------------------------------------------!
  subroutine search_node(nodelist,inode,node_pointer)
    type(nodelist_set) :: nodelist
    type(node_set), pointer :: node, node_pointer 
    integer :: inode
    
    node => nodelist%head
    do while(associated(node%next_node))
       node => node%next_node
       
       if(node%node_param%inode == inode) then
          node_pointer => node
          return
       end if
    end do

  end subroutine search_node

end module module_node
