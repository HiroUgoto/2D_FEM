subroutine check_inside(x,z,element,xi,zeta,ind)
  use module_fem
  implicit none
  type(element_set) :: element
  type(node_set), pointer :: node
  real(8) :: x, z
  real(8) :: xi, zeta
  real(8), allocatable :: xn(:,:)
  integer :: style, ind
  integer :: i
  
  style = 4

  allocate(xn(style,2))
  do i = 1, style
     node => element%node(i)%p
     xn(i,1) = node%node_param%x 
     xn(i,2) = node%node_param%z
  end do

  call x_to_xi(x,z,style,xn,xi,zeta)
  if((-1.d0 <= xi).and.(xi <= 1.d0).and. &
       (-1.d0 <= zeta).and.(zeta <= 1.d0)) then
     ind = 1
  else
     ind = 0
  end if

  deallocate(xn)

end subroutine check_inside


!-----------------------------------------------------------------------------!
subroutine x_to_xi(xref,zref,style,xn,xi,zeta)
  implicit none
  integer :: style
  real(8) :: xref, zref, xi, zeta
  real(8) :: xn(style,2)
  real(8) :: Jv(2), dJv_inv(2,2)
  real(8) :: norm, xi0, zeta0
  integer :: itr

  xi = 0.d0
  zeta = 0.d0
  
  xi0 = xi
  zeta0 = zeta

  do 
     call mk_J(xi0,zeta0,style,xn,xref,zref,Jv)
     call mk_dJ(xi0,zeta0,style,xn,xref,zref,dJv_inv)

     xi   = xi0    - dJv_inv(1,1)*Jv(1) - dJv_inv(1,2)*Jv(2)
     zeta = zeta0  - dJv_inv(2,1)*Jv(1) - dJv_inv(2,2)*Jv(2)

     norm = (xi0-xi)**2.d0 + (zeta0-zeta)**2.d0
!     write(*,*) xi, zeta, norm
     if(norm < 1.d-10) exit
     xi0 = xi
     zeta0 = zeta
  end do

end subroutine x_to_xi

!-----------------------------------------------------------------------------!
subroutine mk_J(xi,zeta,style,xn,xref,zref,Jv)
  implicit none
  integer :: style
  real(8) :: xi, zeta, xn(style,2)
  real(8) :: xref, zref, x, z
  real(8) :: Jv(2)
  real(8) :: dnp(style,2)
  real(8) :: dx(2,2)

  call shape_function_dn(style,xi,zeta,dnp)
  call xi_to_x(xi,zeta,style,xn,x,z)

  dx = matmul(transpose(dnp),xn)

  Jv(1) = -(xref-x)*dx(1,1) - (zref-z)*dx(1,2)
  Jv(2) = -(xref-x)*dx(2,1) - (zref-z)*dx(2,2)

end subroutine mk_J

!-----------------------------------------------------------------------------!
subroutine mk_dJ(xi,zeta,style,xn,xref,zref,dJv_inv)
  implicit none
  integer :: style
  real(8) :: xi, zeta, xn(style,2)
  real(8) :: xref, zref, x, z
  real(8) :: dJv(2,2), dJv_inv(2,2)
  real(8) :: ddnp(style,3), dnp(style,2)
  real(8) :: dx(2,2), ddx(2), det

  call shape_function_ddn(style,xi,zeta,ddnp)
  call shape_function_dn(style,xi,zeta,dnp)
  call xi_to_x(xi,zeta,style,xn,x,z)

  dx = matmul(transpose(dnp),xn)
  ddx = matmul(transpose(xn),ddnp(1:style,3))

  dJv(1,1) = dx(1,1)**2.d0 + dx(1,2)**2.d0
  dJv(2,2) = dx(2,1)**2.d0 + dx(2,2)**2.d0

  dJv(1,2) = dx(1,1)*dx(2,1) + dx(1,2)*dx(2,2) - (xref-x)*ddx(1) - (zref-z)*ddx(2)
  dJv(2,1) = dx(1,1)*dx(2,1) + dx(1,2)*dx(2,2) - (xref-x)*ddx(1) - (zref-z)*ddx(2)

  det = dJv(1,1)*dJv(2,2) - dJv(1,2)*dJv(2,1)
  
  dJv_inv(1,1) =  dJv(2,2)/det
  dJv_inv(2,2) =  dJv(1,1)/det

  dJv_inv(1,2) = -dJv(1,2)/det
  dJv_inv(2,1) = -dJv(2,1)/det

end subroutine mk_dJ

!-----------------------------------------------------------------------------!
subroutine xi_to_x(xi,zeta,style,xn,x,z)
  implicit none
  integer :: style
  real(8) :: xi, zeta, np(style)
  real(8) :: xn(style,2), xyzref(2), x, z

  call shape_function_n(style,xi,zeta,np)

  xyzref = matmul(transpose(xn),np)

  x = xyzref(1)
  z = xyzref(2)

end subroutine xi_to_x

