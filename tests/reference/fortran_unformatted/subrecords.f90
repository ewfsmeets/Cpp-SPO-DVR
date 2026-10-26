! Writes subrecords.bin: gfortran -fmax-subrecord-length=16 subrecords.f90 && ./a.out
! (records of 4, 8 and 40 bytes; the last one is split into subrecords of 16, 16 and 8 bytes)
program rec
  implicit none
  integer :: i
  real(8) :: a(5)
  logical :: l = .true.
  a = (/ (dble(i), i = 1, 5) /)
  open(9, file='subrecords.bin', form='unformatted', status='replace')
  write(9) 7
  write(9) l, 42
  write(9) a
  close(9)
end program
