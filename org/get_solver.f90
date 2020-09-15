subroutine get_solver(fem)
  use module_fem
  implicit none
  type(fem_set) :: fem
  integer :: is

  open(11, file = "mesh/solver.in")
  read(11,*) fem%dt, fem%tim
  fem%ntim = int(fem%tim/fem%dt)

  read(11,*) fem%nsrc

  write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++"
  write(*,*) " Read sources   ", fem%nsrc

  allocate(fem%src(fem%nsrc))
  do is = 1, fem%nsrc
     read(11,*) fem%src(is)%x, fem%src(is)%z, fem%src(is)%ds, fem%src(is)%dip, fem%src(is)%rake, &
          fem%src(is)%tp, fem%src(is)%fp
     write(*,"(i5,10f12.4)") is, fem%src(is)%x, fem%src(is)%z
  end do

  read(11,*) fem%nsite

  write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++"
  write(*,*) " Read stations  ", fem%nsite
  
  allocate(fem%site(fem%nsite))
  do is = 1, fem%nsite
     read(11,*) fem%site(is)%x, fem%site(is)%z
     write(*,"(i5,10f12.4)") is, fem%site(is)%x, fem%site(is)%z
  end do

  close(11)

end subroutine get_solver
