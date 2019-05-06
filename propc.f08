! Code for calculating the formfactor of two polymer molecules
program propc
! From the group of Prof Philip Camp
! subroutines by Julien Sindt, Rui Apóstolo in 2016
! the file was modified and adapted by Joanna Faulds in 2019
! double precision added by Rui Apóstolo in 2019
! University of Edinburgh
use omp_lib
use M_time, only: ud,du
use, intrinsic :: iso_fortran_env
implicit none
include 'variables.inc'

debug = .false.
if (debug .eqv. .true.) then
  open(unit=debugf, file='debug.log', action='write', status='replace')
  write(6,*) "DEBUG MORE ON - CHECK DEBUG.LOG"
end if


! timer base on date_and_time intrinsic function
call timer(ttime,tbeg)
tlast=tbeg

! read from stin so that more than one file can be used
read(5,*)
read(5,'(a)',iostat=ierror)  Inputfile
  if (ierror /= 0) call systemexit("Input file")
read(5,'(a)',iostat=ierror)  Outputprefix
  if (ierror /= 0) call systemexit("Output file")
read(5,'(L6)',iostat=ierror)      b_rg
  if (ierror /= 0) call systemexit("Calculate Rg")
read(5,'(L6)',iostat=ierror)      b_ree
  if (ierror /= 0) call systemexit("Calculate Ree")
read(5,'(L6)',iostat=ierror)      b_pq
  if (ierror /= 0) call systemexit("Calculate p(q)")
read(5,'(L6)',iostat=ierror)      b_trj
  if (ierror /= 0) call systemexit("Calculate output trajectory")
read(5,'(I2)',iostat=ierror)      Columns
  if (ierror /= 0) call systemexit("Columns")
read(5,'(I6)',iostat=ierror)      StepMax
  if (ierror /= 0) call systemexit("StepMax")
read(5,'(I6)',iostat=ierror)      IgnoreFirst
  if (ierror /= 0) call systemexit("IgnoreFirst")
read(5,'(I6)',iostat=ierror)      NMol
  if (ierror /= 0) call systemexit("NMol")
read(5,'(I6)',iostat=ierror)      MolSize
  if (ierror /= 0) call systemexit("MolSize")
read(5,'(I6)',iostat=ierror)      MolStart
  if (ierror /= 0) call systemexit("MolStart")
read(5,'(I6)',iostat=ierror)      StepOutput
  if (ierror /= 0) call systemexit("StepOutput")
read(5,'(E6.0)',iostat=ierror)      lmin
  if (ierror /= 0) call systemexit("lmin")
read(5,'(E6.0)',iostat=ierror)      lmax
  if (ierror /= 0) call systemexit("lmax")
read(5,'(I6)',iostat=ierror)      qpoints
  if (ierror /= 0) call systemexit("qpoints")
read(5,'(a)',iostat=ierror)
  if (ierror >= 0) call systemexit("Too many arguments in inputfile")

! Input sanity checks
if (Inputfile == "") then
  write(6,*) "No name for input file given."
  ! see /usr/include/sysexits.h for error codes
  call EXIT(65)
end if
if (Columns < 5) then
  write(6,*) "'Columns' line in params.in must be an integer bigger than 4. X, Y and Z should be on rows 3, 4 and 5, respectively."
  call EXIT(65)
end if
if (StepMax < 0) then
  write(6,*) "'StepMax' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (StepMax < 1001) then
  write(6,*) StepMax," steps is an unusually low number for StepMax, make sure it's correct."
end if
if (StepMax-IgnoreFirst > 19999) then
  write(6,*) StepMax-IgnoreFirst," steps is an unusually high number of steps, make sure it's correct."
end if
if (StepMax-IgnoreFirst < 1) then
  write(6,*)      "StepMax - IgnoreFirst  must be an integer bigger than 0."
  call EXIT(65)
else if (StepMax-IgnoreFirst < 4999) then
  write(6,*) StepMax-IgnoreFirst," steps is an unusually low number of steps, make sure it's correct."
else if (IgnoreFirst < 0) then
  write(6,*) "'IgnoreFirst' line in params.in must be an integer bigger than or equal to 0."
  call EXIT(65)
end if
if (NMol < 1) then
  write(6,*) "'NMol' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (NMol > 10) then
  write(6,*) NMol," molecules is an unusually high number for NMol, make sure it's correct."
end if
if (MolSize < 1) then
  write(6,*) "'MolSize' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (MolSize < 101) then
  write(6,*) MolSize," atoms is an unusually low number for MolSize, make sure it's correct."
else if (MolSize > 1999) then
  write(6,*) MolSize," atoms is an unusually high number for MolSize, make sure it's correct."
end if
if (MolStart < 1) then
  write(6,*) "'MolStart' line in params.in must be an integer bigger than 0."
  call EXIT(65)
end if
if (StepOutput < 1) then
  write(6,*) "'StepOutput' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (StepOutput > 200) then
  write(6,*) StepOutput," is a large number of steps between console echoes. It might take a while to show."
end if
if (lmin < -3.0) then
  write(6,*) lmin,"is a very low number for the low exponent of the form factor, make sure it's correct."
  call EXIT(65)
else if (lmin > -1.0) then
  write(6,*) lmin,"is a very high number for the low exponent of the form factor, make sure it's correct."
end if
if (lmax <= lmin) then
  write(6,*) "lmax must be higher than lmin"
  call EXIT(65)
else if (lmax < -0.0) then
  write(6,*) lmax,"is a very low number for the high exponent of the form factor, make sure it's correct."
  call EXIT(65)
else if (lmax > 2.0) then
  write(6,*) lmax,"is a very high number for the high exponent of the form factor, make sure it's correct."
end if
if (qpoints < 1) then
  write(6,*) "'p(q)_num_points' line in params.in must be an integer bigger than 0."
else if (qpoints < 100) then
  write(6,*) qpoints,"is a very low number for the number of q points, make sure it's correct."
  call EXIT(65)
else if (qpoints > 200) then
  write(6,*) qpoints,"is a very high number for the number of q points, make sure it's correct."
end if

! ouput parameters
write(6,*)
write(6,*)      "Input file chosen:                      ",Inputfile
write(6,*)
call decwrite(b_rg,"Rg")
call decwrite(b_ree,"Ree")
call decwrite(b_pq,"p(q)")
call decwrite(b_trj,"Edited trj file")
write(6,*)      "Number of Rows in input file:           ",trim(i2str(Columns))
write(6,*)      "Last step number:                       ",trim(i2str(StepMax))
write(6,*)      "Ignoring first:                         ",trim(i2str(IgnoreFirst))," timesteps"
write(6,*)      "Molecule Size:                          ",trim(i2str(MolSize))
write(6,*)      "Type number for first atom in Molecule: ",trim(i2str(MolStart))
write(6,*)      "Status will be echoed on console every: ",trim(i2str(StepOutput))," timesteps"
write(6,*)      "Form factor lower exponent:             ",trim(r2str(lmin))
write(6,*)      "Form factor higher exponent:            ",trim(r2str(lmax))
write(6,*)      "Form factor number of q points:         ",trim(i2str(qpoints))
write(6,*)

if (debug .eqv. .true.) write(debugf,*) "Input done"

! allocate arrays
allocate(timestep(stepmax))
allocate(array(NMol,MolSize,Columns))
allocate(hold(NMol,3))
if (debug .eqv. .true.) write(debugf,*) "allocated variables"
! open input file
open(inputf,file=Inputfile,action='read',status='old')
! open output files
if (b_rg  .eqv. .true.) then
  open(unit=41, file='rg.dat', action='write', status='new')
  write(41,"(A)",advance='no') 'Nstep'
  do l=1,NMol-1
    write(41,"(A,I0)",advance='no') ' rg', l
  end do
  write(41,"(A,I0)") ' rg', NMol
end if
if (b_ree .eqv. .true.) then
  open(unit=42, file='ree.dat', action='write', status='new')
  write(42,"(A)",advance='no') 'Nstep'
  do l=1,NMol-1
    write(42,"(A,I0)",advance='no') ' ree', l
  end do
  write(42,"(A,I0)") ' ree', NMol
end if
if (b_pq  .eqv. .true.) then
  open(unit=43, file='pq.dat', action='write', status='new')
end if
if (b_trj .eqv. .true.) then
  open(unit=44, file='processed.lammpstrj', action='write', status='new')
end if


if (debug .eqv. .true.) write(debugf,*) "opened files"





write(6,*)      "Ignoring first:                         ",trim(i2str(IgnoreFirst))," timesteps"
write(6,*)

if (debug .eqv. .true.) write(debugf,*) "starting skipping"
do Nstep=1,IgnoreFirst
  if (debug .eqv. .true.) write(debugf,*) Nstep
  if (modulo(Nstep,100) == 0) write(6,*) Nstep
  call skipdata
end do
call ETIME(tarray,cpu_time_now)
call texttime(int(cpu_time_now-cpu_time_last),tsl)
write(6, *) "Ignored   ",IgnoreFirst,"  timesteps, took ",tsl
write(6,*)
write(6,*)
cpu_time_last = cpu_time_now
write(6,*) "Now calculating timestep number:"
write(6,*)
write(6,*) Nstep



do Nstep=IgnoreFirst+1,stepmax
  if (debug .eqv. .true.) write(debugf,*) "main loop"
  if (modulo(Nstep,StepOutput) == 0) call timer(Nstep)
    ! call ETIME(tarray,cpu_time_now) ! new clock
    ! call texttime(int(cpu_time_now-cpu_time_last),tsl)
    ! call texttime(int(cpu_time_now-beg_cpu_time),tss)
    ! call texttime(int((cpu_time_now-cpu_time_last)*(((stepmax-Nstep)/100)+1)),eta)
    ! write(6,*) "Timestep: ", Nstep, &
    !            "  This iteration: ", tsl, &
    !            "Running for: ", tss, &
    !            "ETA: ", eta
    ! cpu_time_last = cpu_time_now
  ! end if
  call readheader
  call readdata(MolStart,NMol,MolSize,Natom,Columns)
  if (Nstep==IgnoreFirst+1) then
    if (debug .eqv. .true.) write(debugf,*) "started hold array"
    do mol=1,NMol
      hold(mol,1) = array(mol,1,3)
      hold(mol,2) = array(mol,1,4)
      hold(mol,3) = array(mol,1,5)
    end do
  end if
  call molrebuild
  if (b_rg  .eqv. .true.) call rg(1,molsize,NMol)
  if (b_ree .eqv. .true.) call ree(1,molsize,NMol)
  if (b_pq  .eqv. .true.) call formfactor(1,molsize,NMol)
  if (b_trj .eqv. .true.) call outputtrj(1,molsize,NMol)
end do

! End of program clock call
call ETIME(tarray,end_cpu_time)
call texttime(int(end_cpu_time - beg_cpu_time), tempw)
write(6,*) "Total Time: ", ADJUSTL(tempw)


contains


  !**************************!
  ! Skip data - Rui Apóstolo !
  !**************************!

subroutine skipdata
  implicit none
  integer(sp) :: i
  read(inputf,*)
  read(inputf,*)
  read(inputf,*)
  read(inputf,*) Natom
  if (debug .eqv. .true.) write(debugf,*) "Skipped TS: ", Nstep," Natom: ",Natom
  do i=1,Natom+5
    read(inputf,*)
  end do
end subroutine skipdata


!*****************************!
! Read headers - Julien Sindt !
!*****************************!

subroutine readheader
  implicit none
  if (debug .eqv. .true.) write(debugf,*) "started readheader"
  !Subroutine to read the box dimension at every subsequent timestep
  read(inputf, '(A)') headertext(1) ! ITEM: TIMESTEP
  read(inputf,*) timestep(Nstep)
  read(inputf, '(A)') headertext(2) ! ITEM: NUMBER OF ATOMS
  read(inputf,*) Natom
  read(inputf, '(A)') headertext(3) ! ITEM: BOX BOUNDS pp pp pp
  read(inputf,*) Volume(1,1), Volume(1,2) ! Xmin and Xmax
  Lx = Volume(1,2) - Volume(1,1)
  Hx = Lx/2
  read(inputf,*) Volume(2,1), Volume(2,2) ! Ymin and Ymax
  Ly = Volume(2,2) - Volume(2,1)
  Hy = Ly/2
  read(inputf,*) Volume(3,1), Volume(3,2) ! Zmin and Zmax
  Lz = Volume(3,2) - Volume(3,1)
  Hz = Lz/2
  read(inputf, '(A)') headertext(4) ! ITEM: ATOMS id type x y z
  if (debug .eqv. .true.) then
    write(debugf,*) "readdata"
    write(debugf,*) headertext(1)
    write(debugf,*) timestep(Nstep)
    write(debugf,*) headertext(2)
    write(debugf,*) Natom
    write(debugf,*) headertext(3)
    write(debugf,*) Volume(1,1), Volume(1,2)
    write(debugf,*) Volume(2,1), Volume(2,2)
    write(debugf,*) Volume(3,1), Volume(3,2)
    write(debugf,*) headertext(4)
    write(debugf,*) "Lx,y,z: ",Lx,Ly,Lz
    write(debugf,*) "Hx,y,z: ",Hx,Hy,Hz
  end if
end subroutine readheader


!***************************************!
! Read data and tabulate - Rui Apóstolo !
!***************************************!

subroutine readdata(mstart,nmols,msize,atoms,col)
  implicit none
  !Read data from file and assign it to an array
  integer(sp),intent(in)::nmols,msize,atoms,col,mstart
  integer(sp)::i,j,k,n
  real(dp),dimension(col)::dummy
  if (debug .eqv. .true.) write(debugf,*) "started readdata"
  if (debug .eqv. .true.) write(debugf,*) (atoms-nmols*(msize-1))
  n=1
  k=1
  readloop: do i=1,(atoms-nmols*(msize-1))
    read(inputf,*) (dummy(j),j=1,col)
    if (debug .eqv. .true.) write(debugf,*) i, dummy
    if (int(dummy(2)) == mstart) then
      if (debug .eqv. .true.) write(debugf,*) "went to readdata molecule loop"
      array(n,k,1:col) = dummy(1:col)
      ! if (debug .eqv. .true.) write(debugf,*) (array(n,k,j),j=1,col)
        ! if (debug .eqv. .true.) write(debugf,*) (array(n,k,j),j=1,col)
      molloop: do k=2,msize
        read(inputf,*) (array(n,k,j),j=1,col)
        ! if (debug .eqv. .true.) write(debugf,*) (array(n,k,j),j=1,col)
      end do molloop
      ! if (debug .eqv. .true.) then
        ! do k=1,msize
          ! write(debugf,*) (array(n,k,j),j=1,col)
        ! end do
      ! end if
      n = n+1
      if (debug .eqv. .true.) write(debugf,*) "finishing readdata molecule loop"
    end if
  end do readloop
end subroutine readdata

!***********************************************!
! Molecule Rebuild - Written by Rui Apóstolo    !
! rebuilds molecule so it does not clip         !
! All atoms in same molecule must be contiguous !
!***********************************************!

subroutine molrebuild
  implicit none
  integer(sp)::i,n
  if (debug .eqv. .true.) write(debugf,*) "started molrebuild"
  do n=1,NMol
  ! to avoid clipping of first atom
    if (array(n,1,3) >= hold(n,1)+Hx) then
      array(n,1,3) = array(n,1,3) - Lx
      else if (array(n,1,3) <= hold(n,1)-Hx) then
      array(n,1,3) = array(n,1,3) + Lx
    end if
    hold(n,1) = array(n,1,3)
    
    if (array(n,1,4) >= hold(n,2)+Hy) then
      array(n,1,4) = array(n,1,4) - Ly
      else if (array(n,1,4) <= hold(n,2)-Hy) then
      array(n,1,4) = array(n,1,4) + Ly
    end if
    hold(n,2) = array(n,1,4)
    
    if (array(n,1,5) >= hold(n,3)+Hz) then
      array(n,1,5) = array(n,1,5) - Lz
      else if (array(n,1,5) <= hold(n,3)-Hz) then
      array(n,1,5) = array(n,1,5) + Lz
    end if
    hold(n,3) = array(n,1,5)
    
    do i=2,molsize
      array(n,i,3) = array(n,i,3) - array(n,i-1,3)
      array(n,i,4) = array(n,i,4) - array(n,i-1,4)
      array(n,i,5) = array(n,i,5) - array(n,i-1,5)
      array(n,i,3) = array(n,i,3) - Lx*anint(array(n,i,3)/Lx)
      array(n,i,4) = array(n,i,4) - Ly*anint(array(n,i,4)/Ly)
      array(n,i,5) = array(n,i,5) - Lz*anint(array(n,i,5)/Lz)
      array(n,i,3) = array(n,i-1,3) + array(n,i,3)
      array(n,i,4) = array(n,i-1,4) + array(n,i,4)
      array(n,i,5) = array(n,i-1,5) + array(n,i,5)
    end do
  end do
end subroutine molrebuild

!********************************************
! Rg calculation - Written by Rui Apóstolo
! formula from Soft Matter, 2009, 5, 637–645
!********************************************

subroutine rg(lower,upper,nmols)
  implicit none
  integer(sp),intent(in)::lower,upper,nmols
  real(dp):: xi, yi, zi, rdiff, rtot
  integer(sp):: j,k,n
  if (debug .eqv. .true.) write(debugf,*) "started rg"
  write(41,"(I0,A)",advance="no") Nstep, ' '
  do n=1,nmols
    rdiff = 0.0_dp
      do j=lower,upper-1
        do k=j+1,upper
          xi = array(n,j,3) - array(n,k,3)
          yi = array(n,j,4) - array(n,k,4)
          zi = array(n,j,5) - array(n,k,5)
          rdiff = rdiff + xi**2.0_dp + yi**2.0_dp + zi**2.0_dp
        end do
      end do
    rtot = 1.0_dp*sqrt(rdiff/(molsize**2.0_dp))
    write(41,"(f0.6,a)",advance="no") rtot, ' '
  end do
  write(41,"(A)") ' '
end subroutine rg

!*************************************!
! Calculate Ree -- Rui Apóstolo, 2018 !
!*************************************!

subroutine ree(lower,upper,nmols)
  implicit none
  integer(sp),intent(in)::lower,upper,nmols
  real(dp) :: reetot,xi,yi,zi
  integer(sp) :: i
  if (debug .eqv. .true.) write(debugf,*) "started ree"
  write(42,"(i0,a)",advance="no") Nstep, ' '
  do i=1,nmols
    xi = array(i,upper,3) - array(i,lower,3)
    yi = array(i,upper,4) - array(i,lower,4)
    zi = array(i,upper,5) - array(i,lower,5)
    reetot = 1.0_dp*sqrt(xi**2.0_dp + yi**2.0_dp + zi**2.0_dp)
    write(42,"(f0.6,a)",advance="no") reetot, ' '
  end do
  write(42,"(A)") ' '
end subroutine ree

!***********************************************************************!
! Form Factor calculation - Written by Rui Apóstolo                     !
! appended by Joanna Faulds to account for the minimum image convention.!
!***********************************************************************!

subroutine formfactor(lower,upper,nmols)
  implicit none
  integer(sp),intent(in)::lower,upper,nmols
  real(dp):: xj, yj, zj, qdiff, dummy_variable
  real(dp),dimension(:),allocatable::qvalues,pvalues,q
  integer(sp):: j,k,m,n,o
  if (debug .eqv. .true.) write(debugf,*) "started formfactor"
  allocate(qvalues(0:qpoints-1))
  allocate(pvalues(0:qpoints-1))
  allocate(q(0:qpoints-1))
  qdiff = 0.0_dp
  qvalues = 0.0_dp
  pvalues = 0.0_dp
  q = 0.0_dp
  do m=0,qpoints-1
    ! q(m) = 10.0**((1.0*((abs(lmin)+abs(lmax))/(1.0*qpoints))*m)+lmin)
    q(m) = 10.0_dp**(lmin + real(m)/real(qpoints) * (lmax-lmin))
  end do
  if (qtrig .eqv. .true.) then
    ! call csv_write(43,q,.true.)
    write(43,*) q
    qtrig = .false.
  end if


  !$OMP PARALLEL DO            &
  !$OMP SCHEDULE(STATIC)      &
  !$OMP DEFAULT(SHARED)        &
  !$OMP PRIVATE(m, dummy_variable, j, k, xj, yj, zj, qdiff,o,n) 

  do m = 0, qpoints-1
    ! write(6,*) m, q(m)
    dummy_variable = real(upper*nmols,dp)
    if (debug .eqv. .true.) write(debugf,*) dummy_variable
    do n=1,nmols
      do o=n,nmols
        if (n==o) then
          do j = lower,upper-1 
            do k = j+1,upper
                xj    = array(n,j,3) - array(o,k,3)
                xj= xj - Lx*anint(xj/Lx)
                yj    = array(n,j,4) - array(o,k,4)
                yj= yj - Ly*anint(yj/Ly)
                zj    = array(n,j,5) - array(o,k,5)
                zj= zj - Lz*anint(zj/Lz)
                qdiff = 1.0_dp*sqrt(xj**2 + yj**2 + zj**2)
                dummy_variable  = dummy_variable + 2.0_dp * sin(q(m)*qdiff)/(q(m)*qdiff)
            end do
          end do
        else
          do j = lower,upper
            do k = lower,upper
                xj    = array(n,j,3) - array(o,k,3)
                xj= xj - Lx*anint(xj/Lx)
                yj    = array(n,j,4) - array(o,k,4)
                yj= yj - Ly*anint(yj/Ly)
                zj    = array(n,j,5) - array(o,k,5)
                zj= zj - Lz*anint(zj/Lz)
                qdiff = 1.0_dp*sqrt(xj**2 + yj**2 + zj**2)
                dummy_variable  = dummy_variable + 2.0_dp * sin(q(m)*qdiff)/(q(m)*qdiff)
            end do
          end do
        end if
      end do
    end do
    qvalues(m) = dummy_variable
  end do
  !$OMP END PARALLEL DO

  do m=0,qpoints-1
    pvalues(m)=qvalues(m)/((upper*NMol)**2)
  end do
  ! call csv_write(43,pvalues,.true.)
  write(43,*) pvalues
end subroutine formfactor

!********************************************!
! Data output - Written by Rui Apóstolo      !
! rebuilds lammpstrj with selected molecules !
!********************************************!

subroutine outputtrj(lower,upper,nmols)
  implicit none
  integer(sp),intent(in)::lower,upper,nmols
  integer(sp)::i,j,n
  if (debug .eqv. .true.) write(debugf,*) "started outputtrj"
  write(44, '(A)') trim(headertext(1)) ! ITEM: TIMESTEP
  write(44,*) timestep(Nstep)
  write(44, '(A)') trim(headertext(2)) ! ITEM: NUMBER OF ATOMS
  write(44,*) molsize
  write(44, '(A)') trim(headertext(3)) ! ITEM: BOX BOUNDS pp pp pp
  write(44,*) Volume(1,1), Volume(1,2) ! Xmin and Xmax
  write(44,*) Volume(2,1), Volume(2,2) ! Ymin and Ymax
  write(44,*) Volume(3,1), Volume(3,2) ! Zmin and Zmax
  write(44, '(A)') trim(headertext(4)) ! ITEM: ATOMS id type x y z
  do n=1,nmols
    do i=lower,upper
      write(44,*) int(array(n,i,1)), int(array(n,i,2)), (array(n,i,j),j=3,Columns)
    end do
  end do
end subroutine outputtrj


!******************************************************!
! Exits program due to data input error - Rui Apóstolo !
!******************************************************!

subroutine systemexit(error)
  implicit none
  character(len=*),intent(in) :: error
    write(6,*) "Program exit due to error on data input in variable: ",error
    ! see /usr/include/sysexits.h for error codes
    call EXIT(65)
end subroutine systemexit


!***********************************************************************!
! Writes to console if module is going to be used or not - Rui Apóstolo !
!***********************************************************************!

subroutine decwrite(bool,module)
  implicit none
  character(len=*),intent(in) :: module
  logical,intent(in)::bool
  if (bool .eqv. .true. ) write(6,*) "Calculating                             ",module
  if (bool .eqv. .false.) write(6,*) "Skipping                                ",module
end subroutine


!*********************************************!
! Converts integers to strings - Rui Apóstolo !
!*********************************************!

function i2str(k) result(str)
  implicit none
    integer(sp), intent(in) :: k
    character*30::str
    write (str, *) k
    str = adjustl(str)
end function i2str


!******************************************!
! Converts reals to strings - Rui Apóstolo !
!******************************************!

function r2str(k) result(str)
  implicit none
    real(dp), intent(in) :: k
    character*30::str
    write (str, *) k
    str = adjustl(str)
end function r2str


!******************************************************!
! Time conversion - Written by Rui Apóstolo            !
! Converts seconds to days, hours, minutes and seconds !
!******************************************************!

subroutine texttime(time,ttime)
  implicit none
  integer(sp)::d,s,m,h
  real(dp)::ms
  real(dp),intent(in)::time
  character(LEN=20),intent(out)::ttime
  character(LEN=20)::ttemp
  ttime = ""
  ttemp = ""

  s = 0
  m = 0
  h = 0
  ! 123days 45h 67 min 89 sec 000ms
  ms = module(time,1.0)
  s = modulo(int(time-ms),60)
  m = modulo(int(time-s-ms),3600)
  h = modulo(int(time-m*60-s-ms),86400)
  d = floor(int(time)/86400)

  if (d==0) then
    write(ttemp, '(I4.1, A, I2.2, A, I2.2, A)') d,"d ",h, ":", m, ":", s, ".",ms
  else
    write(ttemp, '(I4.1, A, I2.2, A, I2.2, A)') h, ":", m, ":", s, ".",ms
  end if
  write(ttime, '(A)') trim(ADJUSTL(ttemp))

end subroutine texttime







subroutine timer(steps)
  ! requires subroutine texttime
  intrinsic none
  integer(sp),intent(in)::steps
  ! call cpu_time(cpu_time_now) ! old clock
  call ETIME(tarray,cpu_time_now) ! new clock
  call texttime(int(cpu_time_now-cpu_time_last),tsl)
  call texttime(int(cpu_time_now-beg_cpu_time),tss)
  call texttime(int((cpu_time_now-cpu_time_last)*(((stepmax-steps)/100)+1)),eta)
  write(6,*) "Timestep: ", steps, &
             "  This iteration: ", tsl, &
             "Running for: ", tss, &
             "ETA: ", eta
  cpu_time_last = cpu_time_now
end subroutine timer



!********************************************************************!
! Timer function - Written by Rui Apóstolo                           !
! Calculates run time based on call date_and_time intrinsic function !
!********************************************************************!

subroutine timer(text,newtime,begtime,nsteps,dstep,endtime)
  implicit none
  character(LEN=20),intent(out)::text
  integer(kind=4),dimension(8),intent(out)::newtime
  integer(kind=4),dimension(8),intent(inout),optional::begtime,endtime
  integer(sp),intent(in),optional::nsteps,dstep
  integer(kind=4),dimension(8)::tarray
  character(LEN=20)::temp
  ! call date_and_time(date, time, zone, values)
  ! date, CHARACTER*8, Output, Date, in form CCYYMMDD, where CCYY is the four-digit year, MM the two-digit month, and DD the two-digit day of the month. For example: 19980709 
  ! time, CHARACTER*10, Output, The current time, in the form hhmmss.sss, where hh is the hour, mm minutes, and ss.sss seconds and milliseconds. 
  ! zone, CHARACTER*5, Output, The time difference with respect to UTC, expressed in hours and minutes, in the form hhmm 
  ! values,INTEGER*4 VALUES(8), Output,An integer array of 8 elements described below.
  ! values (1): The year, as a 4-digit integer.
  ! values (2): The month, as an integer from 1 to 12.
  ! values (3): The day of the month, as an integer from 1 to 31.
  ! values (4): The time difference, in minutes, with respect to UTC.
  ! values (5): The hour of the day, as an integer from 1 to 23.
  ! values (6): The minutes of the hour, as an integer from 1 to 59.
  ! values (7): The seconds of the minute, as an integer from 0 to 60.
  ! values (8): The milliseconds of the second, in range 0 to 999.
  text = ""
  temp = ""
  if (present(endtime)) then
    call date_and_time(VALUES=newtime)
    tarray=newtime-endtime

    if (tarray(1:3) > 1) write(6,*) "test successful"




    write(text, '(A)') trim(ADJUSTL(temp))

  else if (present(begtime)) then
    tarray=newtime
    call date_and_time(VALUES=newtime)
    tarray = newtime - tarray



! just use epoch, this doesn't work



    if (tarray(5) == 0) then
      if (tarray(6) == 0) then
        write(temp, '(I2.2, A, I3.3,A)') tarray(7), '.', tarray(8), " seconds"
      else
        write(temp, '(I2.2, A, I2.2, A,I3.3,A)') tarray(6), " m ", tarray(7), '.', tarray(8), " s"
      end if
    else
      write(temp, '(I4.1, A, I2.2, A, I2.2, A,I3.3,A)') tarray(5), " h ", tarray(6), " m ", tarray(7), '.', tarray(8), " s"
    end if
    write(text, '(A)') trim(ADJUSTL(temp))
  else
    call date_and_time(VALUES=newtime)
  end if


end subroutine timer




end program propc
