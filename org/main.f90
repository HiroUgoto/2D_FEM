program main
  use module_fem
  implicit none
  type(fem_set) :: fem
  real(8) :: tim, tp, fp, slip
  integer :: it, is

  call get_mesh(fem)
  call get_solver(fem)

  call set_fem_mesh(fem)
  call set_initial_condition(fem)
  call set_initial_matrix(fem)

  call set_source(fem)
  call set_site(fem)

  write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++"

  open(21, file = "result/vel_x.dat")
  open(22, file = "result/vel_z.dat")
  open(23, file = "result/vel_y.dat")

  tp = fem%src(1)%tp
  fp = fem%src(1)%fp

  do it = 1, fem%ntim
     tim = it*fem%dt

     call ricker_func(tim,tp,fp,slip)
     call update_time(fem,fem%dt,slip)

     write(*,"(f12.7,10f15.7)") tim, fem%site(1)%vel(1), fem%site(1)%vel(3)
     write(21,"(f12.7,10f15.7)") tim, (fem%site(is)%vel(1), is = 1, fem%nsite)
     write(22,"(f12.7,10f15.7)") tim, (fem%site(is)%vel(2), is = 1, fem%nsite)
     write(23,"(f12.7,10f15.7)") tim, (fem%site(is)%vel(3), is = 1, fem%nsite)
  end do
  close(21)
  close(22)
  close(23)

  call erase_nodelist(fem%nodelist)

end program main


!-----------------------------------------------------------------------------!
subroutine ricker_func(tim,tp,fp,a)
  implicit none
  real(8) :: tim, t1, a
  real(8) :: tp, fp, pi

  pi = 2.d0*acos(0.d0)
  t1 = ((tim-tp)*pi*fp)**2.d0

  a = (2.d0*t1-1.d0)*exp(-t1)

end subroutine ricker_func
