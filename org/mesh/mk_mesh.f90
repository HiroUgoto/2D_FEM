program mk_mesh
  use parameter_set
  implicit none
  integer :: nnode_x, nnode_z, nnode
  integer :: nelem_x, nelem_z, nelem
  real(8) :: area_x, area_z
  real(8) :: dx, dz
  real(8) :: dt, tim, dip, rake, tp, fp, ds
  real(8), allocatable :: sx(:), sz(:)
  real(8), allocatable :: sitex(:), sitez(:)
  integer(8), allocatable :: inode(:,:)
  type(node_set), allocatable :: node(:)
  type(element_set), allocatable :: element(:)
  integer :: nsrc, nsite
  integer :: dof
  integer :: i, k, in, ie, is

  !! Set Node !!

  dof = 3    ! 1: 2DSH, 2: 2DPSV, 3: 2DSH & 2DPSV

  area_x = 2500.d0  ![m]
  area_z = 2500.d0  ![m]

  nelem_x = 100
  nelem_z = 100

  nnode_x = 2*nelem_x + 1
  nnode_z = 2*nelem_z + 1
  nnode = nnode_x*nnode_z

  dx = area_x / dble(nelem_x)
  dz = area_z / dble(nelem_z)

  allocate(node(nnode))
  allocate(inode(nnode_x,nnode_z))

  in = 1
  do k = 1, nnode_z
     do i = 1, nnode_x
        inode(i,k) = in

        node(in)%x = 0.5d0*dx*(i-1)
        node(in)%z = 0.5d0*dz*(k-1)

        allocate(node(in)%freedom(dof))
        node(in)%freedom(1:dof) = 1

        in = in + 1
     end do
  end do

  !! Source parameters !!
  nsrc = 1

  allocate(sx(nsrc),sz(nsrc))

  sx(1) = 1250.d0
  sz(1) = 1000.d0

  dip = 0.d0
!  rake = 90.d0
  rake = 45.d0
  
  tp = 0.75d0
  fp = 2.d0
  ds = 10.d0

  !! Set calculation 
  dt = 0.0025d0
  tim = 3.d0

  !! Set observation
  nsite = 5

  allocate(sitex(nsite),sitez(nsite))

  do is = 1, nsite
     sitex(is) = 100.d0*dble(is-1) + sx(1)
     sitez(is) = 0.d0
  end do


  !! Set Mesh !!
  nelem = nelem_x*nelem_z
  allocate(element(nelem))
  
  ie = 1
  do k = 1, nelem_z
     do i = 1, nelem_x

        element(ie)%style = 9

        element(ie)%inode(1) = inode(2*i-1,2*k-1)
        element(ie)%inode(2) = inode(2*i+1,2*k-1)
        element(ie)%inode(3) = inode(2*i+1,2*k+1)
        element(ie)%inode(4) = inode(2*i-1,2*k+1)

        element(ie)%inode(5) = inode(2*i,2*k-1)
        element(ie)%inode(6) = inode(2*i+1,2*k)
        element(ie)%inode(7) = inode(2*i,2*k+1)
        element(ie)%inode(8) = inode(2*i-1,2*k)

        element(ie)%inode(9) = inode(2*i,2*k)

        element(ie)%vs  = 1000.d0
        element(ie)%vp  = 1730.d0
        element(ie)%rho = 2000.d0

        element(ie)%damping(1:4) = 0
        if(i == 1) element(ie)%damping(4) = 1
        if(i == nelem_x) element(ie)%damping(2) = 1
        if(k == nelem_z) element(ie)%damping(3) = 1
       
        ie = ie + 1
     end do
  end do


  !! Output meshfile !!
  open(11, file = "mesh.in")
  write(11,"(3i10)") nnode, nelem, dof

  do in = 1, nnode
     write(11,"(i10,2f20.10,2i10)") in, node(in)%x, node(in)%z, node(in)%freedom(1:dof)
  end do

  do ie = 1, nelem
     write(11,"(i10,i10,3f20.10)") ie, element(ie)%style, element(ie)%vs, element(ie)%vp, element(ie)%rho
     write(11,"(100i10)") element(ie)%inode(1:9), element(ie)%damping(1:4)
  end do

  close(11)

  !! Output solver file !!
  open(12, file = "solver.in")
  
  write(12,"(2f20.10)") dt, tim
  write(12,"(i10)") nsrc
  do is = 1, nsrc
     write(12,"(10f20.10)") sx(is), sz(is), ds, dip, rake, tp, fp
  end do

  write(12,"(i10)") nsite
  do is = 1, nsite
     write(12,"(2g20.10)") sitex(is), sitez(is)
  end do
  
  close(12)



end program mk_mesh


