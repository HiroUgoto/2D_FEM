subroutine shape_function_n(style,xi,zeta,np)
  implicit none
  integer :: style
  real(8) :: xi, zeta
  real(8) :: np(style)

  if(style == 4) then
     np(1) = (1.d0 - xi)*(1.d0 - zeta) / 4.d0
     np(2) = (1.d0 + xi)*(1.d0 - zeta) / 4.d0
     np(3) = (1.d0 + xi)*(1.d0 + zeta) / 4.d0
     np(4) = (1.d0 - xi)*(1.d0 + zeta) / 4.d0
     
  else if(style == 8) then
     np(1) = (1.d0 - xi)*(1.d0 - zeta)*(-1.d0-xi-zeta) / 4.d0
     np(2) = (1.d0 + xi)*(1.d0 - zeta)*(-1.d0+xi-zeta) / 4.d0
     np(3) = (1.d0 + xi)*(1.d0 + zeta)*(-1.d0+xi+zeta) / 4.d0
     np(4) = (1.d0 - xi)*(1.d0 + zeta)*(-1.d0-xi+zeta) / 4.d0

     np(5) = (1.d0 - xi*xi)*(1.d0 - zeta) / 2.d0
     np(6) = (1.d0 + xi)*(1.d0 - zeta*zeta) / 2.d0
     np(7) = (1.d0 - xi*xi)*(1.d0 + zeta) / 2.d0
     np(8) = (1.d0 - xi)*(1.d0 - zeta*zeta) / 2.d0

  else if(style == 9) then
     np(1) =  (1.d0 - xi)*(1.d0 - zeta)*xi*zeta / 4.d0
     np(2) = -(1.d0 + xi)*(1.d0 - zeta)*xi*zeta / 4.d0
     np(3) =  (1.d0 + xi)*(1.d0 + zeta)*xi*zeta / 4.d0
     np(4) = -(1.d0 - xi)*(1.d0 + zeta)*xi*zeta / 4.d0
     
     np(5) = -(1.d0 - xi*xi)*zeta*(1.d0 - zeta) / 2.d0
     np(6) =  (1.d0 + xi)*xi*(1.d0 - zeta*zeta) / 2.d0
     np(7) =  (1.d0 - xi*xi)*zeta*(1.d0 + zeta) / 2.d0
     np(8) = -(1.d0 - xi)*xi*(1.d0 - zeta*zeta) / 2.d0

     np(9) = (1.d0-xi*xi)*(1.d0-zeta*zeta)

  else if(style == 12) then
     np(1) = (1.d0 - xi)*(1.d0 - zeta)*(-10.d0+9.d0*xi**2+9.d0*zeta**2) / 32.d0
     np(2) = (1.d0 + xi)*(1.d0 - zeta)*(-10.d0+9.d0*xi**2+9.d0*zeta**2) / 32.d0
     np(3) = (1.d0 + xi)*(1.d0 + zeta)*(-10.d0+9.d0*xi**2+9.d0*zeta**2) / 32.d0
     np(4) = (1.d0 - xi)*(1.d0 + zeta)*(-10.d0+9.d0*xi**2+9.d0*zeta**2) / 32.d0

     np(5) = (1.d0-3.d0*xi)*(1.d0 - xi**2)*(1.d0 - zeta) * 9.d0/32.d0
     np(6) = (1.d0+3.d0*xi)*(1.d0 - xi**2)*(1.d0 - zeta) * 9.d0/32.d0

     np(7) = (1.d0 + xi)*(1.d0-3.d0*zeta)*(1.d0 - zeta**2) * 9.d0/32.d0
     np(8) = (1.d0 + xi)*(1.d0+3.d0*zeta)*(1.d0 - zeta**2) * 9.d0/32.d0

     np(9)  = (1.d0+3.d0*xi)*(1.d0 - xi**2)*(1.d0 + zeta) * 9.d0/32.d0
     np(10) = (1.d0-3.d0*xi)*(1.d0 - xi**2)*(1.d0 + zeta) * 9.d0/32.d0
     
     np(11) = (1.d0 - xi)*(1.d0+3.d0*zeta)*(1.d0 - zeta**2) * 9.d0/32.d0
     np(12) = (1.d0 - xi)*(1.d0-3.d0*zeta)*(1.d0 - zeta**2) * 9.d0/32.d0

  else
     write(*,*) "    Error: Undefined element style in shape_function_n"
     stop
  end if

end subroutine shape_function_n

!-----------------------------------------------------------------------------!
subroutine shape_function_dn(style,xi,zeta,dnp)
  implicit none
  integer :: style
  real(8) :: xi, zeta
  real(8) :: dnp(style,2)

  if(style == 4) then
     dnp(1,1) = -(1.d0 - zeta) / 4.d0
     dnp(1,2) = -(1.d0 -   xi) / 4.d0
     
     dnp(2,1) =  (1.d0 - zeta) / 4.d0
     dnp(2,2) = -(1.d0 +   xi) / 4.d0
     
     dnp(3,1) =  (1.d0 + zeta) / 4.d0
     dnp(3,2) =  (1.d0 +   xi) / 4.d0
     
     dnp(4,1) = -(1.d0 + zeta) / 4.d0
     dnp(4,2) =  (1.d0 -   xi) / 4.d0
     
  else if(style == 8) then
     dnp(1,1) = (1.d0 - zeta)*(2.d0*xi+zeta) / 4.d0
     dnp(1,2) = (1.d0 -   xi)*(xi+2.d0*zeta) / 4.d0

     dnp(2,1) = (1.d0 - zeta)*(2.d0*xi-zeta) / 4.d0
     dnp(2,2) = -(1.d0 +  xi)*(xi-2.d0*zeta) / 4.d0

     dnp(3,1) = (1.d0 + zeta)*(2.d0*xi+zeta) / 4.d0
     dnp(3,2) = (1.d0 +   xi)*(xi+2.d0*zeta) / 4.d0

     dnp(4,1) = (1.d0 + zeta)*(2.d0*xi-zeta) / 4.d0
     dnp(4,2) = -(1.d0 -  xi)*(xi-2.d0*zeta) / 4.d0

     dnp(5,1) = -xi*(1.d0 - zeta)
     dnp(5,2) = (xi**2-1.d0) / 2.d0

     dnp(6,1) = (1.d0 - zeta**2) / 2.d0
     dnp(6,2) = -(1.d0 + xi)*zeta

     dnp(7,1) = -xi*(1.d0 + zeta)
     dnp(7,2) = (1.d0 - xi**2) / 2.d0

     dnp(8,1) = -(1.d0 - zeta**2) / 2.d0
     dnp(8,2) = -(1.d0 - xi)*zeta

  else if(style == 9) then
     dnp(1,1) =  (2.d0*xi-1.d0)*(zeta-1.d0)*zeta / 4.d0
     dnp(1,2) =  (xi-1.d0)*xi*(2.d0*zeta-1.d0) / 4.d0

     dnp(2,1) =  (2.d0*xi+1.d0)*(zeta-1.d0)*zeta / 4.d0
     dnp(2,2) =  (xi+1.d0)*xi*(2.d0*zeta-1.d0) / 4.d0

     dnp(3,1) =  (2.d0*xi+1.d0)*(zeta+1.d0)*zeta / 4.d0
     dnp(3,2) =  (xi+1.d0)*xi*(2.d0*zeta+1.d0) / 4.d0

     dnp(4,1) =  (2.d0*xi-1.d0)*(zeta+1.d0)*zeta / 4.d0
     dnp(4,2) =  (xi-1.d0)*xi*(2.d0*zeta+1.d0) / 4.d0

     dnp(5,1) =  xi*(1.d0-zeta)*zeta 
     dnp(5,2) =  (1.d0-xi*xi)*(2.d0*zeta-1.d0) / 2.d0

     dnp(6,1) =  (2.d0*xi+1.d0)*(1.d0-zeta*zeta) / 2.d0
     dnp(6,2) =  -xi*(1.d0+xi)*zeta 

     dnp(7,1) =  -xi*(1.d0+zeta)*zeta 
     dnp(7,2) =  (1.d0-xi*xi)*(2.d0*zeta+1.d0) / 2.d0

     dnp(8,1) =  (2.d0*xi-1.d0)*(1.d0-zeta*zeta) / 2.d0
     dnp(8,2) =  xi*(1.d0-xi)*zeta 

     dnp(9,1) =  -2.d0*xi*(1.d0-zeta*zeta)
     dnp(9,2) =  -2.d0*(1.d0-xi*xi)*zeta 

  else
     write(*,*) "    Error: Undefined element style in shape_function_dn"
     stop
  end if

end subroutine shape_function_dn

!-----------------------------------------------------------------------------!
subroutine shape_function_ddn(style,xi,zeta,ddnp)
  implicit none
  integer :: style
  real(8) :: xi, zeta
  real(8) :: ddnp(style,3)

  if(style == 4) then
     ddnp(1,1) = 0.d0
     ddnp(1,2) = 0.d0
     ddnp(1,3) = 1.d0 / 4.d0

     ddnp(2,1) = 0.d0
     ddnp(2,2) = 0.d0
     ddnp(2,3) = -1.d0 / 4.d0
     
     ddnp(3,1) = 0.d0
     ddnp(3,2) = 0.d0
     ddnp(3,3) = 1.d0 / 4.d0
     
     ddnp(4,1) = 0.d0
     ddnp(4,2) = 0.d0
     ddnp(4,3) = -1.d0 / 4.d0
     
  else if(style == 8) then
     ddnp(1,1) =  (1.d0 - zeta) / 2.d0
     ddnp(1,2) =  (1.d0 -   xi) / 2.d0
     ddnp(1,3) =  (-2.d0*xi - 2.d0*zeta + 1.d0)  / 4.d0

     ddnp(2,1) =  (1.d0 - zeta) / 2.d0
     ddnp(2,2) =  (1.d0 +   xi) / 2.d0
     ddnp(2,3) =  (-2.d0*xi + 2.d0*zeta - 1.d0)  / 4.d0

     ddnp(3,1) =  (1.d0 + zeta) / 2.d0
     ddnp(3,2) =  (1.d0 +   xi) / 2.d0
     ddnp(3,3) =  ( 2.d0*xi + 2.d0*zeta + 1.d0)  / 4.d0

     ddnp(4,1) =  (1.d0 + zeta) / 2.d0
     ddnp(4,2) =  (1.d0 -   xi) / 2.d0
     ddnp(4,3) =  (-2.d0*xi + 2.d0*zeta - 1.d0)  / 4.d0

     ddnp(5,1) =  -(1.d0 - zeta)
     ddnp(5,2) =   0.d0
     ddnp(5,3) =   xi

     ddnp(6,1) =   0.d0
     ddnp(6,2) =  -(1.d0 + xi)
     ddnp(6,3) =  -zeta

     ddnp(7,1) =  -(1.d0 + zeta)
     ddnp(7,2) =   0.d0
     ddnp(7,3) =  -xi

     ddnp(8,1) =   0.d0
     ddnp(8,2) =  -(1.d0 - xi)
     ddnp(8,3) =   zeta

  else
     write(*,*) "    Error: Undefined element style in shape_function_ddn"
     stop
  end if

end subroutine shape_function_ddn
