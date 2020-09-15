subroutine set_source(fem)
  use module_fem
  implicit none
  type(fem_set) :: fem
  real(8), allocatable :: xn(:,:)
  real(8), allocatable :: moment(:), b(:,:)
  real(8) :: xi, zeta
  real(8) :: rmu, ds, dip, rake, m0
  integer :: dof, dofs, style
  integer :: is, ie, ind
  integer :: i

  dof = fem%dof
  if(dof == 1) then
     dofs = 2
  else if(dof == 2) then
     dofs = 3
  else if(dof == 3) then
     dofs = 5
  end if

  allocate(moment(dofs))

  write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++"

  do is = 1, fem%nsrc
     do ie = 1, fem%nelem
        call check_inside(fem%src(is)%x,fem%src(is)%z,fem%element(ie),xi,zeta,ind)
        if(ind == 1) then
           style = fem%element(ie)%style
           fem%src(is)%style = style
           allocate(fem%src(is)%node(style))

           allocate(xn(style,2))
           do i = 1, style
              fem%src(is)%node(i)%p => fem%element(ie)%node(i)%p              
              xn(i,1) = fem%src(is)%node(i)%p%node_param%x 
              xn(i,2) = fem%src(is)%node(i)%p%node_param%z
           end do

           rmu = fem%element(ie)%vs**2.d0 * fem%element(ie)%rho

           ds = fem%src(is)%ds
           dip  = fem%src(is)%dip * pi/180.d0
           rake = fem%src(is)%rake * pi/180.d0

           m0 = ds*rmu

           allocate(b(dofs,style*dof))

           if(dof == 1) then
              moment(1) = -m0*sin(dip)*cos(rake)
              moment(2) = -m0*cos(dip)*cos(rake)
           else if(dof == 2) then
              moment(1) = -m0*sin(2*dip)*sin(rake)
              moment(2) =  m0*sin(2*dip)*sin(rake)
              moment(3) = -m0*cos(2*dip)*sin(rake)
           else if(dof == 3) then
              moment(1) = -m0*sin(2*dip)*sin(rake)
              moment(2) =  m0*sin(2*dip)*sin(rake)
              moment(3) = -m0*cos(2*dip)*sin(rake)              
              moment(4) = -m0*sin(dip)*cos(rake)
              moment(5) = -m0*cos(dip)*cos(rake)
           end if

           call mk_b(dof,dofs,style,xn,xi,zeta,b)

           allocate(fem%src(is)%coeff(style*dof))
           fem%src(is)%coeff(1:style*dof) = matmul(transpose(b),moment)
           
           write(*,*) " Source element No.", ie
           write(*,"(2x,10e12.3)") moment

           deallocate(xn,b)

           exit
        end if
     end do
  end do

end subroutine set_source
