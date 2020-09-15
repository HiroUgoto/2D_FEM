!-------------------------------------------------------------------!
subroutine mk_local_mass_matrix(element)
  use module_fem
  implicit none
  type(element_set) :: element
  type(node_set), pointer :: node
  real(8), allocatable :: xn(:,:)
  real(8), allocatable :: xg(:), wxg(:)
  real(8), allocatable :: n(:,:), m(:,:), mg(:,:)
  real(8) :: jacobi(2,2), det
  real(8) :: xi, zeta, wx, wz
  real(8) :: v, mass, m_tr
  integer :: dof, style
  integer :: ngauss
  integer :: i, j

  dof = element%dof
  style = element%style

  allocate(xn(style,2))

  do i = 1, style
     node => element%node(i)%p
     xn(i,1) = node%node_param%x 
     xn(i,2) = node%node_param%z
  end do
  
  ngauss = 5
  allocate(xg(ngauss),wxg(ngauss))
  call gauss_wgst(ngauss,xg,wxg)

  allocate(n(dof,style*dof))
  allocate(m(style*dof,style*dof))
  allocate(mg(style*dof,style*dof))

  v = 0.d0
  mg = 0.d0
  do i = 1, ngauss
     xi = xg(i) ;  wx = wxg(i)
     do j = 1, ngauss
        zeta = xg(j) ;  wz = wxg(j)

        call mk_n(dof,style,xi,zeta,n)
        call mk_m(dof,style,n,m)
        call mk_jacobi(style,xn,xi,zeta,det,jacobi)

        v = v + wx*wz*det
        mg = mg + m * wx*wz*det
     end do
  end do


  mass = element%rho * v

  m_tr = 0.d0
  do i = 1, dof*style
     m_tr = m_tr + mg(i,i)
  end do
  m_tr = m_tr / dble(dof)

  allocate(element%m(style*dof))
  do i = 1, dof*style
     element%m(i) = mg(i,i)* mass/m_tr
  end do

!  write(*,*)
!  write(*,"(100e12.4)") element%m(1:style*dof)

  deallocate(n,m,mg)

end subroutine mk_local_mass_matrix

!-------------------------------------------------------------------!
subroutine mk_local_stiffness_matrix(element)
  use module_fem
  implicit none
  type(element_set) :: element
  type(node_set), pointer :: node
  real(8), allocatable :: xn(:,:)
  real(8), allocatable :: xg(:), wxg(:)
  real(8), allocatable :: b(:,:), d(:,:), k(:,:)
  real(8) :: jacobi(2,2), det
  real(8) :: rmu, rlambda
  real(8) :: xi, zeta, wx, wz
  integer :: dof, dofs, style
  integer :: ngauss
  integer :: i, j

  dof = element%dof
  if(dof == 1) then
     dofs = 2
  else if(dof == 2) then
     dofs = 3
  else if(dof == 3) then
     dofs = 5
  end if

  style = element%style
  allocate(xn(style,2))

  do i = 1, style
     node => element%node(i)%p
     xn(i,1) = node%node_param%x 
     xn(i,2) = node%node_param%z
  end do

  ngauss = 5
  allocate(xg(ngauss),wxg(ngauss))
  call gauss_wgst(ngauss,xg,wxg)

  allocate(d(dofs,dofs))
  rmu = element%vs**2 * element%rho
  rlambda = element%vp**2 * element%rho - 2.d0*rmu
  call mk_d(dof,dofs,rmu,rlambda,d)

  allocate(b(dofs,style*dof))
  allocate(k(style*dof,style*dof))

  allocate(element%k(style*dof,style*dof))
  element%k = 0.d0

  do i = 1, ngauss
     xi = xg(i) ;  wx = wxg(i)
     do j = 1, ngauss
        zeta = xg(j) ;  wz = wxg(j)

        call mk_b(dof,dofs,style,xn,xi,zeta,b)
        call mk_k(dof,dofs,style,b,d,k)
        call mk_jacobi(style,xn,xi,zeta,det,jacobi)

        element%k = element%k + k * wx*wz*det

     end do
  end do

!  write(*,*)
!  do i = 1, style*dof
!     write(*,"(100e12.4)") element%k(i,1:style*dof)
!  end do

  deallocate(b,k)

end subroutine mk_local_stiffness_matrix

!-------------------------------------------------------------------!
subroutine mk_m(dof,style,n,m)
  implicit none
  integer :: dof, style
  real(8) :: n(dof,style*dof), m(style*dof,style*dof)

  m = matmul(transpose(n),n)

end subroutine mk_m

!-------------------------------------------------------------------!
subroutine mk_n(dof,style,xi,zeta,n)
  implicit none
  integer :: dof, style
  real(8) :: xi, zeta
  real(8) :: n(dof,style*dof), np(style)
  integer :: i, j, i1

  call shape_function_n(style,xi,zeta,np)

  n = 0.d0

  i1 = 1
  do i = 1, style
     do j = 1, dof
        n(j,i1) = np(i)
        i1 = i1 + 1
     end do
  end do

end subroutine mk_n

!-------------------------------------------------------------------!
subroutine mk_k(dof,dofs,style,b,d,k)
  implicit none
  integer :: dof, dofs, style
  real(8) :: b(dofs,style*dof), d(dofs,dofs)
  real(8) :: k(style*dof,style*dof)
  real(8) :: btd(style*dof,dofs)

  btd = matmul(transpose(b),d)
  k = matmul(btd,b)

end subroutine mk_k

!-------------------------------------------------------------------!
subroutine mk_d(dof,dofs,rmu,rlambda,d)
  implicit none
  integer :: dof, dofs
  real(8) :: rmu, rlambda
  real(8) :: d(dofs,dofs)

  d = 0.d0

  if(dof == 1) then
     d(1,1) = rmu
     d(2,2) = rmu

  else if(dof == 2) then
     d(1,1) = rlambda + 2.d0*rmu ;  d(1,2) = rlambda 
     d(2,1) = rlambda ;  d(2,2) = rlambda + 2.d0*rmu 
     d(3,3) = rmu 

  else if(dof == 3) then
     d(1,1) = rlambda + 2.d0*rmu ;  d(1,2) = rlambda 
     d(2,1) = rlambda ;  d(2,2) = rlambda + 2.d0*rmu 
     d(3,3) = rmu 

     d(4,4) = rmu
     d(5,5) = rmu
     
  end if

end subroutine mk_d

!-------------------------------------------------------------------!
subroutine mk_b(dof,dofs,style,xn,xi,zeta,b)
  implicit none
  integer :: dof, dofs, style
  real(8) :: xn(style,2)
  real(8) :: xi, zeta
  real(8) :: b(dofs,style*dof)
  real(8) :: dn(style,2)
  integer :: i, i1, i2, i3

  call mk_dn(dof,style,xn,xi,zeta,dn)

  b = 0.d0
  if(dof == 1) then
     do i = 1, style
        b(1,i) = dn(i,1)
        b(2,i) = dn(i,2)
     end do

  else if(dof == 2) then
     do i = 1, style
        i1 = (i-1)*2+1 ;  i2 = (i-1)*2+2

        b(1,i1) = dn(i,1) ;  b(1,i2) =    0.d0
        b(2,i1) =    0.d0 ;  b(2,i2) = dn(i,2)
        b(3,i1) = dn(i,2) ;  b(3,i2) = dn(i,1)
     end do

  else if(dof == 3) then
     b = 0.d0
     do i = 1, style
        i1 = (i-1)*3+1 ;  i2 = (i-1)*3+2 ;  i3 = (i-1)*3+3

        b(1,i1) = dn(i,1) ;  b(1,i2) =    0.d0
        b(2,i1) =    0.d0 ;  b(2,i2) = dn(i,2)
        b(3,i1) = dn(i,2) ;  b(3,i2) = dn(i,1)

        b(4,i3) = dn(i,1)
        b(5,i3) = dn(i,2)
     end do

  end if

end subroutine mk_b

!-------------------------------------------------------------------!
subroutine mk_dn(dof,style,xn,xi,zeta,dn)
  implicit none
  integer :: dof, style
  real(8) :: xn(style,2)
  real(8) :: xi, zeta
  real(8) :: dn(style,2)
  real(8) :: det_inv, jacobi_inv(2,2), dnp(style,2)

  call mk_inv_jacobi(style,xn,xi,zeta,det_inv,jacobi_inv)
  call shape_function_dn(style,xi,zeta,dnp)

  dn = matmul(dnp,jacobi_inv)

end subroutine mk_dn

!-------------------------------------------------------------------!
subroutine mk_inv_jacobi(style,xn,xi,zeta,det_inv,jacobi_inv)
  implicit none
  integer :: style
  real(8) :: xn(style,2)
  real(8) :: xi, zeta
  real(8) :: det, det_inv
  real(8) :: jacobi(2,2), jacobi_inv(2,2)

  call mk_jacobi(style,xn,xi,zeta,det,jacobi)

  jacobi_inv(1,1) =  jacobi(2,2)/det
  jacobi_inv(1,2) = -jacobi(1,2)/det

  jacobi_inv(2,1) = -jacobi(2,1)/det
  jacobi_inv(2,2) =  jacobi(1,1)/det

  det_inv = 1.d0/det

end subroutine mk_inv_jacobi

!-------------------------------------------------------------------!
subroutine mk_jacobi(style,xn,xi,zeta,det,jacobi)
  implicit none
  integer :: style
  real(8) :: xn(style,2)
  real(8) :: xi, zeta
  real(8) :: det, jacobi(2,2)
  real(8) :: dnp(style,2)
  
  call shape_function_dn(style,xi,zeta,dnp)

  jacobi = matmul(transpose(xn),dnp)
  det = jacobi(1,1)*jacobi(2,2) - jacobi(1,2)*jacobi(2,1)

end subroutine mk_jacobi
