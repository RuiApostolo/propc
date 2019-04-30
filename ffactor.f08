! Code for calculating the formfactor of two polymer molecules
program radius
! From the group of Prof Philip Camp
! subroutines by Julien Sindt, Rui Apóstolo in 2016
! the file was modified and adapted by Joanna Faulds in 2019
! double precision added by Rui Apóstolo in 2019
! University of Edinburgh
! csv_file from flibs 0.9 - flibs.sourceforge.net/
  use csv_file
  use omp_lib
  implicit none

  ! real kind parametrs
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)
! integer, parameter :: qp = selected_real_kind(33, 4931)


  real(dp),dimension(:,:),allocatable::array
  real(dp),dimension(3,2)::Volume
  real(dp),dimension(:),allocatable::timestep
  real(dp),dimension(3)::hold
  character(LEN=80),dimension(4)::headertext
  integer(sp)::row,Natom,Nstep,stepmax,input,molsize,IgnoreFirst
  character(LEN=50)::inputfile
  real :: Lx, Ly, Lz, Hx, Hy, Hz
  logical :: qtrig=.true.
  input=22
  inputfile = '../dump.polymer.lammpstrj'
  open(unit=43, file='ff.csv', action='write', status='replace')
  call arraysize
  allocate(timestep(stepmax))
  open(input,file=inputfile)
  allocate(array(molsize,row))
  write(6,*) "Ignoring first:                   ",IgnoreFirst," timesteps"
  do Nstep=1,IgnoreFirst
      if (modulo(Nstep,100) == 0) write(6,*) Nstep
      call readheader
      call readdata
      ! call formfactor(1,molsize)
  end do
  write(6,*) "Now calculating timestep number:"
  write(6,*) Nstep
  do Nstep=IgnoreFirst+1,stepmax
      if (modulo(Nstep,20) == 0) write(6,*) Nstep
      call readheader
      call readdata
      call formfactor(1,molsize)
  end do

contains

  !*******************************************
  ! Retrieve matrix dimensions - Julien Sindt
  !*******************************************

  subroutine arraysize
    ! Defines size of array in file
    ! uncomment if inputfile name is different, check general variables
    ! write(6,*) "Please enter file name"
    ! read*,inputfile
    ! write(6,*) "How many steps?"
    ! set the number of frames to be read in
    ! stepmax = 5001
    stepmax = 2520
    ! write(6,*) "How many columns?"
    ! read*,row, set the number of rows
    row = 5
    ! write(6,*) "How many atoms?", set the number of atoms in the polymer molecule
    molsize = 5896
    ! write(6,*) "Ignore first X timesteps"
    IgnoreFirst = 2500
  end subroutine arraysize

  !*****************************
  ! Read headers - Julien Sindt
  !*****************************

  subroutine readheader
    !Subroutine to read the box dimension at every subsequent timestep
    read(input, '(A)') headertext(1) ! ITEM: TIMESTEP
    read(input,*) timestep(Nstep)
    read(input, '(A)') headertext(2) ! ITEM: NUMBER OF ATOMS
    read(input,*) Natom
    read(input, '(A)') headertext(3) ! ITEM: BOX BOUNDS pp pp pp
    read(input,*) Volume(1,1), Volume(1,2) ! Xmin and Xmax
    Lx = Volume(1,2) - Volume(1,1)
    Hx = Lx/2
    read(input,*) Volume(2,1), Volume(2,2) ! Ymin and Ymax
    Ly = Volume(2,2) - Volume(2,1)
    Hy = Ly/2
    read(input,*) Volume(3,1), Volume(3,2) ! Zmin and Zmax
    Lz = Volume(3,2) - Volume(3,1)
    Hz = Lz/2
    read(input, '(A)') headertext(4) ! ITEM: ATOMS id type x y z
   end subroutine readheader


  !***************************************
  ! Read data and tabulate - Julien Sindt
  !***************************************

  subroutine readdata
    !Read data from file and assign it to an array
    integer(sp)::i,j
    do i=1,molsize
      read(input,*) (array(i,j),j=1,row)
    end do
    do i=molsize+1,Natom
      read(input,*)
    end do
  end subroutine readdata

  !********************************************
  ! Molecule Rebuild - Written by Rui Apóstolo
  ! rebuilds molecule so it does not clip
  ! only valid for one sorted molecule
  ! It can't be used for two molecules.
  !********************************************

  subroutine molrebuild
    integer(sp)::i
    ! to avoid clipping of first atom
    if (array(1,3) >= hold(1)+Hx) then
      array(1,3) = array(1,3) - Lx
      else if (array(1,3) <= hold(1)-Hx) then
      array(1,3) = array(1,3) + Lx
    end if
    hold(1) = array(1,3)
    
    if (array(1,4) >= hold(2)+Hy) then
      array(1,4) = array(1,4) - Ly
      else if (array(1,4) <= hold(2)-Hy) then
      array(1,4) = array(1,4) + Ly
    end if
    hold(2) = array(1,4)
    
    if (array(1,5) >= hold(3)+Hz) then
      array(1,5) = array(1,5) - Lz
      else if (array(1,5) <= hold(3)-Hz) then
      array(1,5) = array(1,5) + Lz
    end if
    hold(3) = array(1,5)
    
    do i=2,molsize
      array(i,3) = array(i,3) - array(i-1,3)
      array(i,4) = array(i,4) - array(i-1,4)
      array(i,5) = array(i,5) - array(i-1,5)
      array(i,3) = array(i,3) - Lx*anint(array(i,3)/Lx)
      array(i,4) = array(i,4) - Ly*anint(array(i,4)/Ly)
      array(i,5) = array(i,5) - Lz*anint(array(i,5)/Lz)
      array(i,3) = array(i-1,3) + array(i,3)
      array(i,4) = array(i-1,4) + array(i,4)
      array(i,5) = array(i-1,5) + array(i,5)
    end do
  end subroutine molrebuild


  !********************************************
  ! Rg calculation - Written by Rui Apóstolo
  ! formula from Soft Matter, 2009, 5, 637–645
  !********************************************

  subroutine smrg(lower,upper)
    real(dp) xi, yi, zi, rdiff, rtot
    integer(sp) j,k,lower,upper
    rdiff = 0.0_dp
    do j=lower,upper-1
      do k=j+1,upper
        xi = array(j,3) - array(k,3)
        yi = array(j,4) - array(k,4)
        zi = array(j,5) - array(k,5)
        rdiff = rdiff + xi**2 + yi**2 + zi**2
      end do
    end do
    rtot = 1.0_dp*sqrt(rdiff/(molsize**2))
    write(42,*) Nstep, rtot
  end subroutine smrg


  !**************************************************
  ! Form Factor calculation - Written by Rui Apóstolo
  ! appended by Joanna Faulds to account for the minimum image convention.
  !**************************************************

  subroutine formfactor(lower,upper)
    real(dp) xj, yj, zj, qdiff, lmin, lmax, dummy_variable
    real(dp),dimension(:),allocatable::qvalues,pvalues,q
    integer(sp) j,k,m,lower,upper,qpoints
    ! write(6,*) "hello"
    ! write(6,*) size(q)
    qdiff = 0.0_dp
    ! if (Nstep == 1) then
        lmin = -3.0_dp
        lmax = 1.0_dp
        qpoints = 200
        allocate(qvalues(0:qpoints-1))
        allocate(pvalues(0:qpoints-1))
        allocate(q(0:qpoints-1))
        q = 0.0_dp
        do m=0,qpoints-1
          ! q(m) = 10.0**((1.0*((abs(lmin)+abs(lmax))/(1.0*qpoints))*m)+lmin)
          q(m) = 10.0_dp**(lmin + real(m)/real(qpoints) * (lmax-lmin))
        end do
      if (qtrig .eqv. .true.) then
        call csv_write(43,q,.true.)
        qtrig = .false.
      end if
    ! end if
    qvalues = 0.0_dp
    pvalues = 0.0_dp

    !$OMP PARALLEL DO            &
    !$OMP SCHEDULE(STATIC)      &
    !$OMP DEFAULT(SHARED)        &
    !$OMP PRIVATE(m, dummy_variable, j, k, xj, yj, zj, qdiff) 
    

    do m = 0, qpoints-1
      ! write(6,*) m, q(m)
      dummy_variable = real(upper,dp)
      do j = lower,upper-1
        do k = j+1,upper
            xj    = array(j,3) - array(k,3)
            xj= xj - Lx*anint(xj/Lx)
            yj    = array(j,4) - array(k,4)
            yj= yj - Ly*anint(yj/Ly)
            zj    = array(j,5) - array(k,5)
            zj= zj - Lz*anint(zj/Lz)
            qdiff = 1.0_dp*sqrt(xj**2 + yj**2 + zj**2)
            dummy_variable  = dummy_variable + 2.0_dp * sin(q(m)*qdiff)/(q(m)*qdiff)
        end do
      end do
      qvalues(m) = dummy_variable
    end do
    
    !$OMP END PARALLEL DO

    do m=0,qpoints-1
      pvalues(m)=qvalues(m)/(upper**2)
    end do
    call csv_write(43,pvalues,.true.)
  end subroutine formfactor

  !*******************************************
  ! Data output - Written by Rui Apóstolo
  ! rebuilds lammpstrj with only pol molecule
  !*******************************************

  subroutine output
    integer(sp)::i,j
    write(50, '(A)') trim(headertext(1)) ! ITEM: TIMESTEP
    write(50,*) timestep(Nstep)
    write(50, '(A)') trim(headertext(2)) ! ITEM: NUMBER OF ATOMS
    write(50,*) molsize
    write(50, '(A)') trim(headertext(3)) ! ITEM: BOX BOUNDS pp pp pp
    write(50,*) Volume(1,1), Volume(1,2) ! Xmin and Xmax
    write(50,*) Volume(2,1), Volume(2,2) ! Ymin and Ymax
    write(50,*) Volume(3,1), Volume(3,2) ! Zmin and Zmax
    write(50, '(A)') trim(headertext(4)) ! ITEM: ATOMS id type x y z
    do i=1,molsize
      write(50,*) int(array(i,1)), int(array(i,2)), (array(i,j),j=3,row)
    end do
  end subroutine output

end program radius
