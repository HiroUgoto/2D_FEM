subroutine mk_local_damping_matrix(element)
  use module_fem
  implicit none
  type(element_set) :: element
  real(8), allocatable :: xg(:), wxg(:)
  real(8), allocatable :: xi_set(:), zeta_set(:)
  real(8), allocatable :: q(:,:), n(:,:), imp(:,:), nqn(:,:)
  real(8) :: xn2(2,2)
  real(8) :: length
  integer :: dof, style
  integer :: ngauss
  integer :: i, ig

  dof = element%dof
  style = element%style

  ngauss = 5
  allocate(xg(ngauss),wxg(ngauss))
  call gauss_wgst(ngauss,xg,wxg)

  allocate(xi_set(ngauss),zeta_set(ngauss))

  allocate(element%c(style*dof,style*dof))
  element%c = 0.d0

  do i = 1, 4
     if(element%damping(i) == 1) then
        call set_xz_on_edge(element,i,xn2,ngauss,xg,xi_set,zeta_set)

        allocate(q(dof,dof),imp(dof,dof))
        allocate(n(dof,dof*style),nqn(dof*style,dof*style))
        
        call mk_q(dof,xn2,q,length)
        call mk_imp(dof,element%vs,element%vp,element%rho,imp)

        do ig = 1, ngauss
           call mk_n(dof,style,xi_set(ig),zeta_set(ig),n)
           call mk_nqn(dof,style,n,q,imp,nqn)
           element%c = element%c + nqn * wxg(ig) * length/2.d0
        end do

        deallocate(q,imp,n,nqn)
     end if
  end do
  
!  write(*,*)
!  do i = 1, style*dof
!     write(*,"(100e12.4)") element%c(i,1:style*dof)
!  end do

end subroutine mk_local_damping_matrix

!-------------------------------------------------------------------!
subroutine set_xz_on_edge(element,dn,xn2,ngauss,xg,xi_set,zeta_set)
  use module_fem
  implicit none
  type(element_set) :: element
  type(node_set), pointer :: node1, node2
  real(8) :: xn2(2,2)
  integer :: dn
  integer :: ngauss
  real(8) :: xg(ngauss), xi_set(ngauss), zeta_set(ngauss)


  if(dn == 1) then
     node1 => element%node(1)%p
     node2 => element%node(2)%p

     xi_set   =    xg
     zeta_set = -1.d0

  else if(dn == 2) then
     node1 => element%node(2)%p
     node2 => element%node(3)%p

     xi_set   =  1.d0
     zeta_set =    xg

  else if(dn == 3) then
     node1 => element%node(3)%p
     node2 => element%node(4)%p

     xi_set   =    xg
     zeta_set =  1.d0

  else if(dn == 4) then
     node1 => element%node(4)%p
     node2 => element%node(1)%p

     xi_set   = -1.d0
     zeta_set =    xg

  end if

  xn2(1,1) = node1%node_param%x
  xn2(1,2) = node1%node_param%z

  xn2(2,1) = node2%node_param%x
  xn2(2,2) = node2%node_param%z

!  write(*,*) dn, node1%node_param%inode, node2%node_param%inode
!  write(*,*) xi_set
!  write(*,*) zeta_set

end subroutine set_xz_on_edge

!-------------------------------------------------------------------!
subroutine mk_nqn(dof,style,n,q,imp,nqn)
  implicit none
  integer :: dof, style
  real(8) :: n(dof,dof*style), q(dof,dof), imp(dof,dof)
  real(8) :: nqn(dof*style,dof*style)
  real(8) :: qiq(dof,dof)

  qiq = matmul(matmul(transpose(q),imp),q)
  nqn = matmul(matmul(transpose(n),qiq),n)
  
end subroutine mk_nqn

!-------------------------------------------------------------------!
subroutine mk_q(dof,xn,q,length)
  implicit none
  integer :: dof
  real(8) :: xn(2,2), q(dof,dof), length
  real(8) :: n(2), t(2)
  real(8) :: norm_n, norm_t
  
  t(1) = xn(2,1) - xn(1,1)
  t(2) = xn(2,2) - xn(1,2)

  n(1) =  t(2)
  n(2) = -t(1)

  norm_n = sqrt(n(1)*n(1) + n(2)*n(2))
  norm_t = sqrt(t(1)*t(1) + t(2)*t(2))

  length = norm_t

  n(1:2) = n(1:2)/norm_n
  t(1:2) = t(1:2)/norm_t

  if(dof == 1) then
     q(1,1) = 1.d0
  else if(dof == 2) then
     q(1,1) = n(1) ; q(1,2) = n(2) 
     q(2,1) = t(1) ; q(2,2) = t(2) 
  else if(dof == 3) then
     q(1,1) = n(1) ; q(1,2) = n(2) ; q(1,3) = 0.d0 
     q(2,1) = t(1) ; q(2,2) = t(2) ; q(2,3) = 0.d0
     q(3,1) = 0.d0 ; q(3,2) = 0.d0 ; q(3,3) = 1.d0
  end if

end subroutine mk_q

!-------------------------------------------------------------------!
subroutine mk_imp(dof,vs,vp,rho,imp)
  implicit none
  integer :: dof
  real(8) :: vs, vp, rho
  real(8) :: imp(dof,dof)

  if(dof == 1) then
     imp(1,1) = rho*vs
  else if(dof == 2) then
     imp(1,1) = rho*vp ; imp(1,2) =   0.d0 
     imp(2,1) =   0.d0 ; imp(2,2) = rho*vs 
  else if(dof == 3) then
     imp(1,1) = rho*vp ; imp(1,2) =   0.d0 ; imp(1,3) =   0.d0
     imp(2,1) =   0.d0 ; imp(2,2) = rho*vs ; imp(2,3) =   0.d0
     imp(3,1) =   0.d0 ; imp(3,2) =   0.d0 ; imp(3,3) = rho*vs
  end if

end subroutine mk_imp
