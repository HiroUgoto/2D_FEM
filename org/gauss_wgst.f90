subroutine gauss_wgst(n,x,w)
  implicit none
  real(8), parameter :: eps = 1.d-14
  real(8), parameter :: pi = 3.141592653589793238426433832795d0 
  integer :: n
  real(8) :: x(n), w(n)
  real(8) :: p1, p2, p3, pp, z, z1, xl, xm
  integer :: i, m, j

  m = (n+1)/2
  xm = 0.d0
  xl = 1.d0

  do i = 0, m-1
     z = cos(pi*dble(i+0.75d0)/(dble(n)+0.5d0))
     
     do
        p1 = 1.d0
        p2 = 0.d0
        do j = 0, n-1
           p3 = p2
           p2 = p1
           p1 = ((2.d0*j+1.d0)*z*p2-dble(j)*p3)/dble(j+1)
        end do
     
        pp = dble(n)*(z*p1-p2)/(z*z-1.d0)
        z1 = z
        z = z1-p1/pp
        if(abs(z-z1) < eps) exit
     end do
        
     x(i+1) = xm - xl*z
     x(n-i) = xm + xl*z
     w(i+1) = 2.d0*xl/((1.d0-z*z)*pp*pp)
     w(n-i) = w(i+1)
  end do

end subroutine gauss_wgst
