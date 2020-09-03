PROGRAM MAIN  

  implicit none

  ! REAL(8)    :: A,B
  integer :: num_args, ix
  character(len=12), dimension(:), allocatable :: args

  integer :: doy
  real(8) :: glat, glon, kp, ut
  real(8) :: w(2), mw(2)

  ! A = 1
  ! B = 2
  
  num_args = command_argument_count()
  allocate(args(num_args))  ! I've omitted checking the return status of the allocation 
  
  if (num_args .lt. 5) ERROR STOP "testhead DOY UT GLAT GLON KP"
     

  do ix = 1, num_args
     call get_command_argument(ix,args(ix))
     ! now parse the argument as you wish
  end do
     
  read(args(1),*) doy
  read(args(2),*) ut
  read(args(3),*) glat
  read(args(4),*) glon
  read(args(5),*) kp
  ! doy = int(args(1))
  ! ut = args(2)
  ! glat = args(3)
  ! glon = args(4)
  ! kp = args(5)
  
  print *,'doy:',doy
  print *,'ut:',ut
  print *,'glat:',glat
  print *,'glon:',glon
  print *,'kp:',kp
  
  call hltwim(doy, ut, glat, glon, kp, w, mw)

  print '(I5,f8.1,f6.1,3(f12.3,f10.3))', doy, ut, kp, glat, glon, w, mw
  

  ! PRINT*, A+B, 
END PROGRAM MAIN
