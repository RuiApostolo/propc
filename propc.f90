! Code for calculating the formfactor of two polymer molecules
program propc
! From the group of Prof. Philip Camp
! original subroutines by Julien Sindt, Rui Apóstolo in 2016
! minimum image convention implemented by Joanna Faulds in 2019
! double precision added by Rui Apóstolo in 2019
! average vs total Rg added by Rui Apóstolo in 2021
! weighted form factor added by Rui Apóstolo in 2021
! University of Edinburgh
! version 2.3 -- 2021/11/30
! changelog from 2.2:
!  * added a new, optional weighted form factor calculation. To use it, you need a weights.in file
!  * added check for inputfile, and error message if it doesn't exist.
!  * added checks for outputfiles, and error messages if they already exist.
!  * incremented version number to params.in
!  * added version number to command-line log
use omp_lib
use, intrinsic :: iso_fortran_env
implicit none
include 'variables.f90'

write(6,*)
write(6,*) "propc version: ", trim(version)
write(6,*)

debug = .true.
if (debug .eqv. .true.) then
  open(unit=debugf, file='debug.log', action='write', status='replace')
  write(6,*) "DEBUG MODE ON - CHECK DEBUG.LOG"
  write(debugf,*) "propc version: ", trim(version)
end if

! Clock start
call tick(beg_time)

! read from stin so that more than one file can be used
read(5,*)
read(5,'(a)',iostat=ierror)  Inputfile
  if (debug .eqv. .true.) write(debugf,*) "Inputfile", Inputfile
  if (ierror /= 0) call systemexit("Input file")
read(5,'(a)',iostat=ierror)  Outputprefix
  if (debug .eqv. .true.) write(debugf,*) "Outputprefix", Outputprefix
  if (ierror /= 0) call systemexit("Output file")
read(5,'(L6)',iostat=ierror)      b_rg
  if (debug .eqv. .true.) write(debugf,*) "b_rg", b_rg
  if (ierror /= 0) call systemexit("Calculate Rg")
read(5,'(L6)',iostat=ierror)      b_rg_ind
  if (debug .eqv. .true.) write(debugf,*) "b_rg_ind", b_rg_ind
  if (ierror /= 0) call systemexit("One Rg value per molecule? (t= 1 Rg per molecule, f= Rg of all atoms, e.g., micelle)")
read(5,'(L6)',iostat=ierror)      b_ree
  if (debug .eqv. .true.) write(debugf,*) "b_ree", b_ree
  if (ierror /= 0) call systemexit("Calculate Ree")
read(5,'(L6)',iostat=ierror)      b_pq
  if (debug .eqv. .true.) write(debugf,*) "b_pq", b_pq
  if (ierror /= 0) call systemexit("Calculate p(q)")
read(5,'(L6)',iostat=ierror)      b_w_pq
  if (debug .eqv. .true.) write(debugf,*) "b_w_pq", b_w_pq
  if (ierror /= 0) call systemexit("Weight p(q)")
read(5,'(L6)',iostat=ierror)      b_pq_ind
  if (debug .eqv. .true.) write(debugf,*) "b_pq_ind", b_pq_ind
  if (ierror /= 0) call systemexit("Calculate individual p(q)")
read(5,'(L6)',iostat=ierror)      b_trj
  if (debug .eqv. .true.) write(debugf,*) "b_trj", b_trj
  if (ierror /= 0) call systemexit("Calculate output trajectory")
read(5,'(I2)',iostat=ierror)      Columns
  if (debug .eqv. .true.) write(debugf,*) "Columns", Columns
  if (ierror /= 0) call systemexit("Columns")
read(5,'(I6)',iostat=ierror)      StepMax
  if (debug .eqv. .true.) write(debugf,*) "StepMax", StepMax
  if (ierror /= 0) call systemexit("StepMax")
read(5,'(I6)',iostat=ierror)      IgnoreFirst
  if (debug .eqv. .true.) write(debugf,*) "IgnoreFirst", IgnoreFirst
  if (ierror /= 0) call systemexit("IgnoreFirst")
read(5,'(I6)',iostat=ierror)      NMol
  if (debug .eqv. .true.) write(debugf,*) "NMol", NMol
  if (ierror /= 0) call systemexit("NMol")
read(5,'(I6)',iostat=ierror)      MolSize
  if (debug .eqv. .true.) write(debugf,*) "MolSize", MolSize
  if (ierror /= 0) call systemexit("MolSize")
read(5,'(I6)',iostat=ierror)      MolStart
  if (debug .eqv. .true.) write(debugf,*) "MolStart", MolStart
  if (ierror /= 0) call systemexit("MolStartType")
read(5,'(I6)',iostat=ierror)      ReeFirst
  if (debug .eqv. .true.) write(debugf,*) "ReeFirst", ReeFirst
  if (ierror /= 0) call systemexit("ReeFirstAtom")
read(5,'(I6)',iostat=ierror)      ReeLast
  if (debug .eqv. .true.) write(debugf,*) "ReeLast", ReeLast
  if (ierror /= 0) call systemexit("ReeLastAtom")
read(5,'(I6)',iostat=ierror)      StepOutput
  if (debug .eqv. .true.) write(debugf,*) "StepOutput", StepOutput
  if (ierror /= 0) call systemexit("StepOutput")
read(5,'(E6.0)',iostat=ierror)      lmin
  if (debug .eqv. .true.) write(debugf,*) "lmin", lmin
  if (ierror /= 0) call systemexit("lmin")
read(5,'(E6.0)',iostat=ierror)      lmax
  if (debug .eqv. .true.) write(debugf,*) "lmax", lmax
  if (ierror /= 0) call systemexit("lmax")
read(5,'(I6)',iostat=ierror)      qpoints
  if (debug .eqv. .true.) write(debugf,*) "qpoints", qpoints
  if (ierror /= 0) call systemexit("qpoints")
read(5,'(a)',iostat=ierror)
  if (ierror >= 0) call systemexit("Too many arguments in inputfile")


if (debug .eqv. .true.) then
  write(debugf,*) 
  write(debugf,*) 
  write(debugf,*) "b_pq", b_pq
  write(debugf,*) "b_pq_ind", b_pq_ind
  write(debugf,*) b_pq .and. b_pq_ind
end if



! Input sanity checks
write(6,*) "Checking parameters:"
! see /usr/include/sysexits.h for error codes
if (Inputfile == "") then
  write(6,*) "No name for input file given."
  call_exit = .true.
end if
if ((b_rg_ind .eqv. .true.) .and. (b_rg .eqv. .false.)) then
  write(6,*) "Calculating the Rg for each molecule needs the Rg to be calculated at all. Changed b_rg to true."
  b_rg = .true.
end if
if ((b_pq_ind .eqv. .true.) .and. (b_pq .eqv. .false.)) then
  write(6,*) "Individual P(q) requires normal P(q) calculation. Changed b_pq to true."
  b_pq = .true.
end if
if ((b_pq_ind .eqv. .true.) .and. (NMol < 2)) then
  write(6,*) "Individual P(q) only valid for systems with more than one molecule. Changed b_ind_q to false."
  b_pq_ind = .false.
end if
if (Columns < 5) then
  write(6,*) "'Columns' line in params.in must be an integer bigger than 4. X, Y and Z should be on rows 3, 4 and 5, respectively."
  call_exit = .true.
end if
if (StepMax < 0) then
  write(6,*) "'StepMax' line in params.in must be an integer bigger than 0."
  call_exit = .true.
else if (StepMax < 1001) then
  write(6,*) StepMax," steps is an unusually low number for StepMax, make sure it's correct."
end if
if (StepMax-IgnoreFirst > 19999) then
  write(6,*) StepMax-IgnoreFirst," steps is an unusually high number of steps, make sure it's correct."
end if
if (StepMax-IgnoreFirst < 1) then
  write(6,*)      "StepMax - IgnoreFirst  must be an integer bigger than 0."
  call_exit = .true.
else if (StepMax-IgnoreFirst < 4999) then
  write(6,*) StepMax-IgnoreFirst," steps is an unusually low number of steps, make sure it's correct."
else if (IgnoreFirst < 0) then
  write(6,*) "'IgnoreFirst' line in params.in must be an integer bigger than or equal to 0."
  call_exit = .true.
end if
if (NMol < 1) then
  write(6,*) "'NMol' line in params.in must be an integer bigger than 0."
  call_exit = .true.
else if (NMol > 10) then
  write(6,*) NMol," molecules is an unusually high number for NMol, make sure it's correct."
end if
if (MolSize < 1) then
  write(6,*) "'MolSize' line in params.in must be an integer bigger than 0."
  call_exit = .true.
else if (MolSize > 1999) then
  write(6,*) MolSize," atoms is an unusually high number for MolSize, make sure it's correct."
end if
if (MolStart < 1) then
  write(6,*) "'MolStart' line in params.in must be an integer bigger than 0."
  call_exit = .true.
end if
if (ReeFirst < 0) then
  write(6,*) "'ReeFirstAtom' line in params.in must be an integer bigger than or equal to 0."
  call_exit = .true.
else if (ReeFirst > MolSize) then
  write(6,*) "ReeFirstAtom can't be outside the target molecule, ReeFirstAtom <= MolSize."
  call_exit = .true.
end if
if (ReeLast < 0) then
  write(6,*) "'ReeLastAtom' line in params.in must be an integer bigger than or equal to 0."
  call_exit = .true.
else if (ReeLast > MolSize) then
  write(6,*) "ReeLastAtom can't be outside the target molecule, ReeLastAtom <= MolSize."
  call_exit = .true.
end if
if (StepOutput < 1) then
  write(6,*) "'StepOutput' line in params.in must be an integer bigger than 0."
  call_exit = .true.
else if (StepOutput > 200) then
  write(6,*) StepOutput," is a large number of steps between console echoes. It might take a while to show."
end if
if (lmin < -3.0) then
  write(6,*) lmin,"is a very low number for the low exponent of the form factor, make sure it's correct."
  call_exit = .true.
else if (lmin > -1.0) then
  write(6,*) lmin,"is a very high number for the low exponent of the form factor, make sure it's correct."
end if
if (lmax <= lmin) then
  write(6,*) "lmax must be higher than lmin"
  call_exit = .true.
else if (lmax < -0.0) then
  write(6,*) lmax,"is a very low number for the high exponent of the form factor, make sure it's correct."
  call_exit = .true.
else if (lmax > 2.0) then
  write(6,*) lmax,"is a very high number for the high exponent of the form factor, make sure it's correct."
end if
if (qpoints < 1) then
  write(6,*) "'p(q)_num_points' line in params.in must be an integer bigger than 0."
else if (qpoints < 100) then
  write(6,*) qpoints,"is a very low number for the number of q points, make sure it's correct."
  call_exit = .true.
else if (qpoints > 200) then
  write(6,*) qpoints,"is a very high number for the number of q points, make sure it's correct."
end if

! file checks
write(6,*) ""
write(6,*) "Checking input and output files:"
inquire(file=Inputfile, exist=input_exists)
if (input_exists .eqv. .false.) then
    write(6,*) "      Input file ",trim(Inputfile)," doesn't exist. Exiting."
    call_exit = .true.
end if
if (b_w_pq .eqv. .true.) then
  inquire(file="weights.in", exist=weight_exists)
  if (weight_exists .eqv. .false.) then
      write(6,*) "      Input file weights.in doesn't exist. Exiting."
      call_exit = .true.
  end if
end if

if (b_rg .eqv. .true.) then
  call check_file("rg.dat")
  if (b_rg_ind .eqv. .true.) then
    call check_file("rg_per_molecule.dat")
  end if
end if

if (b_ree .eqv. .true.) then
  call check_file("ree.dat")
end if

if (b_pq .eqv. .true.) then
  call check_file("pq.dat")
  if (b_pq_ind .eqv. .true.) then
    call check_file("pq_ind.dat")
    call check_file("pq_all.dat")
    call check_file("pq_diff.dat")
  end if
end if

if (b_trj .eqv. .true.) then
  call check_file(trim(adjustl(Outputprefix//"processed.lammpstrj")))
end if

! exit because of input errors
if (call_exit .eqv. .true.) then
  call EXIT(65)
end if

! ouput parameters
write(6,*)
write(6,*)      "Input file chosen:                      ",Inputfile
write(6,*)
call decwrite(b_rg,"Rg")
if (b_rg_ind .eqv. .true.) then
  call decwrite(b_rg_ind,"Rg for each molecule")
else
  call decwrite(.true., "Rg for all atoms as one unit")
end if
call decwrite(b_ree,"Ree")
call decwrite(b_pq,"p(q)")
call decwrite(b_w_pq,"weighted p(q)")
call decwrite(b_pq_ind,"individual p(q)")
call decwrite(b_trj,"Edited trj file")
write(6,*)      "Number of Rows in input file:            ",trim(i2str(Columns))
write(6,*)      "Last step number:                        ",trim(i2str(StepMax))
write(6,*)      "Ignoring first:                          ",trim(i2str(IgnoreFirst))," timesteps"
write(6,*)      "Molecule Size:                           ",trim(i2str(MolSize))
write(6,*)      "Number of molecules                      ",trim(i2str(NMol))
write(6,*)      "Type number for first atom in Molecule:  ",trim(i2str(MolStart))
if (ReeFirst == 0) ReeFirst = 1
write(6,*)      "Index of first atom for Ree calculation: ",trim(i2str(ReeFirst))
if (ReeLast == 0) ReeLast = MolSize
write(6,*)      "Index of last atom for Ree calculation:  ",trim(i2str(ReeLast))
write(6,*)      "Status will be echoed on console every:  ",trim(i2str(StepOutput))," timesteps"
write(6,*)      "Form factor lower exponent:              ",trim(r2str(lmin))
write(6,*)      "Form factor higher exponent:             ",trim(r2str(lmax))
write(6,*)      "Form factor number of q points:          ",trim(i2str(qpoints))
write(6,*)

if (debug .eqv. .true.) write(debugf,*) "Input done"

! allocate arrays
allocate(timestep(stepmax))
allocate(array(Columns,MolSize,NMol))
allocate(hold(NMol,3))
if (b_pq .eqv. .true.) then
  allocate(total_pq(0:qpoints-1))
  total_pq = 0.0_dp
end if
if (b_pq_ind .eqv. .true.) then
  allocate(ind_pq(0:qpoints-1))
  allocate(diff_pq(0:qpoints-1))
  ind_pq = 0.0_dp
  diff_pq = 0.0_dp
end if
if (debug .eqv. .true.) write(debugf,*) "allocated variables"
! open input file
open(inputf,file=Inputfile,action='read',status='old')
if (b_w_pq .eqv. .true.) call read_weights("weights.in")
if ((b_w_pq .eqv. .false.) .and. (input_exists .eqv. .true.)) then
  ! get max atom_type
  if (debug .eqv. .true.) write(debugf,*) "Setting up default weights"
  Nstep = 1
  call readheader
  call readdata(MolStart,NMol,MolSize,Natom,Columns)
  if (debug .eqv. .true.) write(debugf,*) "default weights loop - data read successfully"
  Atom_Max_g = maxval(maxval(array(2,:,:), dim=1),dim=1)
  if (debug .eqv. .true.) write(debugf,*) "max atom type", Atom_Max_g
  ! allocate weights
  allocate(weights(Atom_Max_g))
  weights = 0.0_dp
  do l=1,molsize
    do l2=1,Nmol
      weights(int(array(2,l,l2))) = 1.0_dp
    end do
  end do
  ! cleanup
  deallocate(array)
  rewind(inputf)
  allocate(array(Columns,MolSize,NMol))
end if
! open output files
if (b_rg  .eqv. .true.) then
  if (b_rg_ind .eqv. .true.) then
    open(unit=41, file='rg_per_molecule.dat', action='write', status='new')
    write(41,"(A)",advance='no') 'Nstep'
    do l=1,NMol-1
      write(41,"(A,I0)",advance='no') ' rg', l
    end do
    write(41,"(A,I0)") ' rg', NMol
  else
    open(unit=41, file='rg.dat', action='write', status='new')
    write(41,"(A)") 'Nstep rg_total'
  end if
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
if (b_pq_ind  .eqv. .true.) then
  open(unit=45, file='pq_ind.dat', action='write', status='new')
  open(unit=46, file='pq_all.dat', action='write', status='new')
  open(unit=47, file='pq_diff.dat', action='write', status='new')
end if
if (b_trj .eqv. .true.) then
  open(unit=44, file=trim(adjustl(Outputprefix//"processed.lammpstrj")), action='write', status='new')
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
! call ETIME(tarray,time_now)
call sleep(1)
! write(6,*) tock(beg_time)
call texttime(tock(beg_time),tsl)
call tick(ready_time)
call tick(last_time)
write(6, *) "Ignored   ",IgnoreFirst,"  timesteps, took ",tsl
write(6,*)
write(6,*)
write(6,*) "Timestep:    ", &
           "This iteration:                 ", &
           "Running for:                    ", &
           "ETA:"


do Nstep=IgnoreFirst+1,stepmax
  if (debug .eqv. .true.) write(debugf,*) "main loop"
  if (modulo(Nstep,StepOutput) == 0) call timer(Nstep)
  call readheader
  call readdata(MolStart,NMol,MolSize,Natom,Columns)
  if (Nstep==IgnoreFirst+1) then
    if (debug .eqv. .true.) write(debugf,*) "started hold array"
    do mol=1,NMol
      hold(mol,1) = array(3,1,mol)
      hold(mol,2) = array(4,1,mol)
      hold(mol,3) = array(5,1,mol)
    end do
    call calc_weights
    call formfactorsetup(1,molsize,Nmol)
    if (debug .eqv. .true.) write(debugf,*) "tot_weight: ", tot_weight
  end if
  call molrebuild
  if ((b_rg .eqv. .true.) .and. (b_rg_ind .eqv. .true.)) call rg_ind(1,molsize,NMol)
  if ((b_rg .eqv. .true.) .and. (b_rg_ind .eqv. .false.)) call rg_tot(1,molsize,NMol)
  if (b_ree .eqv. .true.) call ree(ReeFirst,ReeLast,NMol)
  if ((b_pq .eqv. .true.) .and. (b_pq_ind .eqv. .false.)) call formfactor(1,molsize,NMol,total_pq)
  if ((b_pq .eqv. .true.) .and. (b_pq_ind .eqv. .true.)) call formfactor(1,molsize,NMol,total_pq,ind_pq,diff_pq)
  if (b_trj .eqv. .true.) call outputtrj(1,molsize,NMol)
end do

if (b_pq .eqv. .true.) then
  do l=0,qpoints-1
    write(43,*) 10.0_dp**(lmin + real(l)/real(qpoints) * (lmax-lmin)), total_pq(l)
  end do
end if

if (b_pq_ind .eqv. .true.) then
  write(46,*) "#q tot_pq ind_pq tot_div_ind_pq diff_pq"
  do l=0,qpoints-1
    write(45,*) 10.0_dp**(lmin + real(l)/real(qpoints) * (lmax-lmin)), ind_pq(l)
    write(46,*) 10.0_dp**(lmin + real(l)/real(qpoints) * (lmax-lmin)), total_pq(l), ind_pq(l), (total_pq(l)/ind_pq(l)), diff_pq(l)
    write(47,*) 10.0_dp**(lmin + real(l)/real(qpoints) * (lmax-lmin)), diff_pq(l)
  end do
end if

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
  readloop: do i=1,(atoms-nmols*(msize-1))
    read(inputf,*) (dummy(j),j=1,col)
    if (debug .eqv. .true.) write(debugf,*) i, dummy
    if (int(dummy(2)) == mstart) then
      k=1
      if (debug .eqv. .true.) write(debugf,*) "went to readdata molecule loop"
      array(1:col,k,n) = dummy(1:col)
      ! if (debug .eqv. .true.) write(debugf,*) (array(n,k,j),j=1,col)
        ! if (debug .eqv. .true.) write(debugf,*) (array(n,k,j),j=1,col)
      molloop: do k=2,msize
        read(inputf,*) (array(j,k,n),j=1,col)
        ! if (debug .eqv. .true.) write(debugf,*) (array(n,k,j),j=1,col)
      end do molloop
      ! if (debug .eqv. .true.) then
        ! do k=1,msize
          ! write(debugf,*) (array(n,k,j),j=1,col)
        ! end do
      ! end if
      n = n + 1
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
    if (array(3,1,n) >= hold(n,1)+Hx) then
      array(3,1,n) = array(3,1,n) - Lx
    else if (array(3,1,n) <= hold(n,1)-Hx) then
      array(3,1,n) = array(3,1,n) + Lx
    end if
    hold(n,1) = array(3,1,n)
    
    if (array(4,1,n) >= hold(n,2)+Hy) then
      array(4,1,n) = array(4,1,n) - Ly
    else if (array(4,1,n) <= hold(n,2)-Hy) then
      array(4,1,n) = array(4,1,n) + Ly
    end if
    hold(n,2) = array(4,1,n)
    
    if (array(5,1,n) >= hold(n,3)+Hz) then
      array(5,1,n) = array(5,1,n) - Lz
    else if (array(5,1,n) <= hold(n,3)-Hz) then
      array(5,1,n) = array(5,1,n) + Lz
    end if
    hold(n,3) = array(5,1,n)
    
    do i=2,molsize
      array(3,i,n) = array(3,i,n) - array(3,i-1,n)
      array(4,i,n) = array(4,i,n) - array(4,i-1,n)
      array(5,i,n) = array(5,i,n) - array(5,i-1,n)
      array(3,i,n) = array(3,i,n) - Lx*anint(array(3,i,n)/Lx)
      array(4,i,n) = array(4,i,n) - Ly*anint(array(4,i,n)/Ly)
      array(5,i,n) = array(5,i,n) - Lz*anint(array(5,i,n)/Lz)
      array(3,i,n) = array(3,i-1,n) + array(3,i,n)
      array(4,i,n) = array(4,i-1,n) + array(4,i,n)
      array(5,i,n) = array(5,i-1,n) + array(5,i,n)
    end do
  end do
end subroutine molrebuild


!*************************************************************!
! Average Rg calculation - averages Rg from several molecules !
! Written by Rui Apóstolo 2016                                !
! formula (13) from Soft Matter, 2009, 5, 637–645             !
!*************************************************************!

subroutine rg_ind(lower,upper,nmols)
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
        xi = array(3,j,n) - array(3,k,n)
        yi = array(4,j,n) - array(4,k,n)
        zi = array(5,j,n) - array(5,k,n)
        rdiff = rdiff + xi**2.0_dp + yi**2.0_dp + zi**2.0_dp
      end do
    end do
    rtot = 1.0_dp*sqrt(rdiff/(molsize**2.0_dp))
    write(41,"(f0.6,a)",advance="no") rtot, ' '
  end do
  write(41,"(A)") ' '
end subroutine rg_ind


!*********************************************************************!
! Total Rg calculation - calculates one Rg for the clump of molecules !
! Adapted from the rg_ind subroutine by Rui Apóstolo in 2021          !
! formula (13) from Soft Matter, 2009, 5, 637–645                     !
!*********************************************************************!

subroutine rg_tot(lower,upper,nmols)
  implicit none
  integer(sp),intent(in)::lower,upper,nmols
  real(dp):: xi, yi, zi, rdiff, rtot
  integer(sp):: j,k,n,m
  if (debug .eqv. .true.) write(debugf,*) "started rg"
  write(41,"(I0,A)",advance="no") Nstep, ' '
  rdiff = 0.0_dp
  !$OMP PARALLEL DO            &
  !$OMP SCHEDULE(AUTO)      &
  !$OMP DEFAULT(none)        &
  !$OMP SHARED(array, nmols, lower, upper, Lx, Ly, Lz) &
  !$OMP PRIVATE(j, k, xi, yi, zi, n, m) &
  !$OMP REDUCTION(+:rdiff)
  omol: do n=1,nmols
    imol: do m=n,nmols
      if (n == m) then
        oatom1: do j=lower,upper-1
          iatom1: do k=j+1,upper
            xi = array(3,j,m) - array(3,k,n)
            xi = xi - real(Lx*anint(xi/Lx),dp)
            yi = array(4,j,m) - array(4,k,n)
            yi = yi - real(Lx*anint(yi/Lx),dp)
            zi = array(5,j,m) - array(5,k,n)
            zi = zi - real(Lx*anint(zi/Lx),dp)
            rdiff = rdiff + xi**2.0_dp + yi**2.0_dp + zi**2.0_dp
          end do iatom1
        end do oatom1
      else
        oatom2: do j=lower,upper
          iatom2: do k=lower,upper
            xi = array(3,j,m) - array(3,k,n)
            xi = xi - real(Lx*anint(xi/Lx),dp)
            yi = array(4,j,m) - array(4,k,n)
            yi = yi - real(Lx*anint(yi/Lx),dp)
            zi = array(5,j,m) - array(5,k,n)
            zi = zi - real(Lx*anint(zi/Lx),dp)
            rdiff = rdiff + xi**2.0_dp + yi**2.0_dp + zi**2.0_dp
          end do iatom2
        end do oatom2
      end if 
    end do imol
  end do omol
  !$OMP END PARALLEL DO
  rtot = 1.0_dp*sqrt(rdiff/((Nmol*molsize)**2.0_dp))
  write(41,"(f0.6,a)") rtot
end subroutine rg_tot

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
    xi = array(3,upper,i) - array(3,lower,i)
    yi = array(4,upper,i) - array(4,lower,i)
    zi = array(5,upper,i) - array(5,lower,i)
    reetot = 1.0_dp*sqrt(xi**2.0_dp + yi**2.0_dp + zi**2.0_dp)
    write(42,"(f0.6,a)",advance="no") reetot, ' '
  end do
  write(42,"(A)") ' '
end subroutine ree


!***********************************************************************!
! Form Factor setup - Written by Rui Apóstolo                           !
! Subroutine call necessary before any call to formfactor.              !
!***********************************************************************!
subroutine formfactorsetup(lower, upper, nmols)
  implicit none
  integer(sp),intent(in)::lower,upper,nmols
  integer(sp) :: a, m
  if (debug .eqv. .true.) write(debugf,*) "started formfactor setup"
  allocate(q(0:qpoints-1))
  q = 0.0_dp
  do m=0,qpoints-1
    ! q(m) = 10.0**((1.0*((abs(lmin)+abs(lmax))/(1.0*qpoints))*m)+lmin)
    q(m) = 10.0_dp**(lmin + real(m,dp)/real(qpoints,dp) * (lmax-lmin))
  end do
  start_qvalue = 0.0_dp
  start_indqvalue = 0.0_dp
  do a = lower, upper
    start_indqvalue = start_indqvalue + weights(int(array(2,a,1)))**2.0_dp
  end do
  start_qvalue = start_indqvalue * real(NMols, dp)
  if (debug .eqv. .true.) write(debugf,*) "start_qvalue: ", start_qvalue
  if (debug .eqv. .true.) write(debugf,*) "start_indqvalue: ", start_indqvalue
  if (debug .eqv. .true.) write(debugf,*) "finished formfactor setup"
end subroutine formfactorsetup


!***********************************************************************!
! Form Factor calculation - Written by Rui Apóstolo                     !
! appended by Joanna Faulds to account for the minimum image convention.!
! call formfactorsetup once before calling this subroutine              !
!***********************************************************************!

subroutine formfactor(lower,upper,nmols,t_pq,i_pq,d_pq)
  implicit none
  integer(sp),intent(in)::lower,upper,nmols
  real(dp),dimension(:),allocatable,intent(inout)::t_pq
  real(dp),dimension(:),allocatable,intent(inout),optional::i_pq,d_pq
  real(dp):: xj, yj, zj, qdiff, weight
  real(dp),dimension(:),allocatable::qvalues,pvalues,ind_qvalues,ind_pvalues,diff_qvalues
  integer(sp):: j,k,n,o
  if (debug .eqv. .true.) write(debugf,*) "started formfactor"
  allocate(qvalues(0:qpoints-1))
  allocate(pvalues(0:qpoints-1))
  if (present(i_pq)) then
    if (present(d_pq)) then
      allocate(ind_qvalues(0:qpoints-1))
      allocate(ind_pvalues(0:qpoints-1))
      allocate(diff_qvalues(0:qpoints-1))
    else
      call systemexit("missing d_pq")
    end if
  end if
  qdiff = 0.0_dp
  pvalues = 0.0_dp

  ! if (b_pq .eqv. .true.)     qvalues = real(upper*nmols,dp)
  if (b_pq .eqv. .true.)     qvalues = start_qvalue
  ! if (b_pq_ind .eqv. .true.) ind_qvalues = real(upper,dp)
  if (b_pq_ind .eqv. .true.) ind_qvalues = start_indqvalue
  if (b_pq_ind .eqv. .true.) diff_qvalues = 0.0_dp

if (present(i_pq) .and. present(d_pq)) then
  !$OMP PARALLEL DO            &
  !$OMP SCHEDULE(AUTO)      &
  !$OMP DEFAULT(none)        &
  !$OMP SHARED(array,q, Lx, Ly, Lz, nmols, lower, upper, weights) &
  !$OMP PRIVATE(j, k, xj, yj, zj, qdiff,o,n,weight) &
  !$OMP REDUCTION(+:qvalues,ind_qvalues,diff_qvalues)
    ! write(6,*) m, q(m)
    ! if (debug .eqv. .true.) write(debugf,*) dummy_variable
    do n=1,nmols
      do o=n,nmols
        if (n==o) then
          do j = lower,upper-1 
            do k = j+1,upper
                xj = array(3,j,n) - array(3,k,o)
                xj = xj - real(Lx*anint(xj/Lx),dp)
                yj = array(4,j,n) - array(4,k,o)
                yj = yj - real(Ly*anint(yj/Ly),dp)
                zj = array(5,j,n) - array(5,k,o)
                zj = zj - real(Lz*anint(zj/Lz),dp)
                qdiff = 1.0_dp*dsqrt(xj**2.0_dp + yj**2.0_dp + zj**2.0_dp)
                weight = weights(int(array(2,j,n))) * weights(int(array(2,k,o)))
                ! tot_weight = tot_weight + weights(int(array(2,j,n))) + weights(int(array(2,k,o)))
                qvalues  = qvalues + 2.0_dp * weight * sin(q*qdiff)/(q*qdiff)
                ! qvalues  = qvalues + 2.0_dp * sin(q*qdiff)/(q*qdiff)
                ind_qvalues = ind_qvalues + 2.0_dp * weight * sin(q*qdiff)/(q*qdiff)
            end do
          end do
        else
          do j = lower,upper
            do k = lower,upper
                xj = array(3,j,n) - array(3,k,o)
                xj = xj - real(Lx*anint(xj/Lx),dp)
                yj = array(4,j,n) - array(4,k,o)
                yj = yj - real(Ly*anint(yj/Ly),dp)
                zj = array(5,j,n) - array(5,k,o)
                zj = zj - real(Lz*anint(zj/Lz),dp)
                qdiff = 1.0_dp*dsqrt(xj**2.0_dp + yj**2.0_dp + zj**2.0_dp)
                weight = weights(int(array(2,j,n))) * weights(int(array(2,k,o)))
                ! tot_weight = tot_weight + weights(int(array(2,j,n))) + weights(int(array(2,k,o)))
                qvalues  = qvalues + 2.0_dp * weight * sin(q*qdiff)/(q*qdiff)
                ! qvalues  = qvalues + 2.0_dp * sin(q*qdiff)/(q*qdiff)
                diff_qvalues = diff_qvalues + 2.0_dp * weight * sin(q*qdiff)/(q*qdiff)
            end do
          end do
        end if
      end do
    end do
  !$OMP END PARALLEL DO
  ind_pvalues = ind_qvalues / real(tot_weight/nmols,dp)
  ! diff_pvalues = diff_qvalues / real(nmols*(upper**2),dp)
  pvalues= qvalues / real(tot_weight,dp)
  t_pq = t_pq + pvalues/real(StepMax-IgnoreFirst,dp)
  i_pq = i_pq + ind_pvalues/real(StepMax-IgnoreFirst,dp)
  d_pq = d_pq + diff_qvalues/real(StepMax-IgnoreFirst,dp)
else
  !$OMP PARALLEL DO            &
  !$OMP SCHEDULE(AUTO)      &
  !$OMP DEFAULT(none)        &
  !$OMP SHARED(array,q, Lx, Ly, Lz, nmols, lower, upper, weights) &
  !$OMP PRIVATE(j, k, xj, yj, zj, qdiff,o,n,weight) &
  !$OMP REDUCTION(+:qvalues)
    ! write(6,*) m, q(m)
    ! if (debug .eqv. .true.) write(debugf,*) dummy_variable
    do n=1,nmols
      do o=n,nmols
        if (n==o) then
          do j = lower,upper-1 
            do k = j+1,upper
                xj    = array(3,j,n) - array(3,k,o)
                xj= xj - real(Lx*anint(xj/Lx),dp)
                yj    = array(4,j,n) - array(4,k,o)
                yj= yj - real(Ly*anint(yj/Ly),dp)
                zj    = array(5,j,n) - array(5,k,o)
                zj= zj - real(Lz*anint(zj/Lz),dp)
                qdiff = 1.0_dp*dsqrt(xj**2.0_dp + yj**2.0_dp + zj**2.0_dp)
                weight = weights(int(array(2,j,n))) * weights(int(array(2,k,o)))
                ! tot_weight = tot_weight + weights(int(array(2,j,n))) + weights(int(array(2,k,o)))
                qvalues  = qvalues + 2.0_dp * weight * sin(q*qdiff)/(q*qdiff)
                ! qvalues  = qvalues + 2.0_dp * sin(q*qdiff)/(q*qdiff)
            end do
          end do
        else
          do j = lower,upper
            do k = lower,upper
                xj    = array(3,j,n) - array(3,k,o)
                xj= xj - real(Lx*anint(xj/Lx),dp)
                yj    = array(4,j,n) - array(4,k,o)
                yj= yj - real(Ly*anint(yj/Ly),dp)
                zj    = array(5,j,n) - array(5,k,o)
                zj= zj - real(Lz*anint(zj/Lz),dp)
                qdiff = 1.0_dp*dsqrt(xj**2.0_dp + yj**2.0_dp + zj**2.0_dp)
                weight = weights(int(array(2,j,n))) * weights(int(array(2,k,o)))
                ! tot_weight = tot_weight + weights(int(array(2,j,n))) + weights(int(array(2,k,o)))
                qvalues  = qvalues + 2.0_dp * weight * sin(q*qdiff)/(q*qdiff)
                ! qvalues  = qvalues + 2.0_dp * sin(q*qdiff)/(q*qdiff)
            end do
          end do
        end if
      end do
    end do
  !$OMP END PARALLEL DO
  pvalues= qvalues / real(tot_weight,dp)
  t_pq = t_pq + pvalues/real(StepMax-IgnoreFirst,dp)
end if
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
      write(44,*) int(array(1,i,n)), int(array(2,i,n)), (array(j,i,n),j=3,Columns)
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
  if (bool .eqv. .true. ) write(6,*) "Calculating                              ",module
  if (bool .eqv. .false.) write(6,*) "Skipping                                 ",module
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


!***************************************!
! Calculate total weight - Rui Apóstolo !
!***************************************!

subroutine calc_weights
  implicit none
  integer(sp)::j

  tot_weight = 0.0_dp
  do j = 1, molsize
    ! tot_weight = tot_weight + (weights(int(array(2,j,n))))**2.0_dp
    tot_weight = tot_weight + weights(int(array(2,j,1)))
  end do
  tot_weight = (real(NMol, dp) * tot_weight)**2.0_dp

end subroutine calc_weights


!************************************************!
! Time conversion - Written by Rui Apóstolo      !
! Converts seconds to hours, minutes and seconds !
!************************************************!

subroutine texttime(seconds,ttime)
  implicit none
  integer(kind=8)::seconds
  integer::s,m,h,r
  character(LEN=32)::ttime,ttemp
  ttime = ""
  ttemp = ""
  s = 0
  r = 0
  m = 0
  h = 0
  ! 1234 h 56 min 78 sec
  s = modulo(seconds,60)
  r = (seconds - s)/60
  m = modulo(r,60)
  h = (r - m)/60
  if (h == 0) then
    if (m == 0) then
      write(ttemp, '(I2.2, A)') s, " seconds "
    else
      write(ttemp, '(I2.2, A, I2.2, A)') m, " min ", s, " sec "
    end if
  else
    write(ttemp, '(I4.1, A, I2.2, A, I2.2, A)') h, " h ", m, " min ", s, " sec "
  end if
  write(ttime, '(A)') trim(ADJUSTL(ttemp))
end subroutine texttime


!**********************************************************************************!
! Wall clock timing subroutine                                                     !
! from: stackoverflow.com/questions/5083051/timing-a-fortran-multithreaded-program !
! implemented by Rui Apóstolo in 2021                                              !
! *********************************************************************************!

subroutine tick(t)
    integer(kind=8),intent(OUT)::t

    call system_clock(t)
end subroutine tick

! returns time in seconds from now to time described by t
function tock(t) result(r)
    integer(kind=8),intent(in)::t
    integer(kind=8)::now,clock_rate,r

    call system_clock(now,clock_rate)
    r =real(now - t)/real(clock_rate)
end function tock

subroutine timer(steps)
  ! requires subroutine texttime
  integer(sp),intent(in)::steps
  time_since_last = tock(last_time)
  time_since_ready = tock(ready_time)
  call texttime(time_since_last,tsl)
  call texttime(time_since_ready,tss)
  call texttime(int(((real(StepMax-IgnoreFirst, dp)/(steps-IgnoreFirst))*(time_since_ready))-time_since_ready,kind=8),eta)
  write(6,*) steps, " ", tsl, tss, eta
  call tick(last_time)
end subroutine timer

!********************************************
! Exits program by user choice - Rui Apóstolo
!********************************************

subroutine userexit(answer)
  character,intent(in) :: answer*1
  select case(answer)
    case("y")
      write(6,*) "Program continuing by your choice."
    case("Y")
      write(6,*) "Program continuing by your choice."
    case default
      write(6,*) "Program exit by your choice."
      call EXIT(0)
  end select
end subroutine userexit


!******************************************
! Checks if the file exists -- Rui Apóstolo
!******************************************

subroutine check_file(filename)
  implicit none
  character(len=*), intent(in) :: filename
  logical :: file_exists=.false.
  inquire(file=filename, exist=file_exists)
  if (file_exists .eqv. .true.) then
    write(6,*) "      File ",trim(filename)," already exists. Rename or delete before running the program."
    call_exit = .true.
  end if
end subroutine check_file


!**********************************
! Read weights file -- Rui Apóstolo
!**********************************

subroutine read_weights(filename)
  implicit none
  character(len=*), intent(in) :: filename
  integer::io,n,atom_max,temp1,i
  real::temp2
  ! open file
  open(weightsf,file=filename,action='read',status='old')
  ! check max atype
  n = 0
  atom_max = 0
  do
    read(weightsf, *, iostat=io) temp1, temp2
    if (io /= 0) exit
    n = n + 1
    if (temp1 > atom_max) then
      atom_max = temp1
    end if
  end do
  rewind(weightsf)
  ! read file
  allocate(weights(atom_max))
  weights = 0.0_dp
  do i = 1, n
    read(weightsf, *, iostat=io) temp1, weights(temp1)
  end do
end subroutine read_weights

end program propc
