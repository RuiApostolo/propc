 ! real kind parametrs
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)
! integer, parameter :: qp = selected_real_kind(33, 4931)

! real arrays
real(dp),dimension(:,:,:),allocatable::array
real(dp),dimension(:,:),allocatable::hold
real(dp),dimension(3,2)::Volume
real(dp),dimension(:),allocatable::timestep

! real
real(dp) :: Lx, Ly, Lz, Hx, Hy, Hz, lmin, lmax

! integers
integer(sp)::l,ierror,Natom,Nstep,NMol,mol,qpoints,inputf=20,debugf=90
integer(sp)::Columns,StepMax,IgnoreFirst,Molsize,MolStart,StepOutput

! logicals
logical :: qtrig=.true.
logical :: debug,b_rg,b_ree,b_pq,b_trj

! strings
character(LEN=80),dimension(4)::headertext
character::Inputfile*64,Outputprefix*64,userinput*1

! New clock variables
! going to ignore the array, it gives user and system time in seconds, in this order.
real, dimension(2) :: tarray 
real(sp) ::  beg_cpu_time, end_cpu_time, cpu_time_now, cpu_time_last
character(LEN=20)::tempw, tsl, tss, eta