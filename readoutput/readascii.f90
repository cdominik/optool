
program readasciifortran
  implicit none
  character*80 :: filename
  character*160 :: string
  integer :: iformat,nlam,nmu,ilam,iline,iang,nang
  logical :: scat
  double precision              :: dummy4(1:4),dummy6(1:6),dum
  double precision, allocatable :: lam(:),kabs(:),ksca(:),gsca(:)
  double precision, allocatable :: ang(:)
  double precision, allocatable :: f11(:,:),f12(:,:),f22(:,:),f33(:,:),f34(:,:),f44(:,:)


  !
  ! Set the file name.  This code can read both files with and without
  ! the scattering matrix, by looking at the value of iformat.
  !
  filename = "dustkapscatmat.dat"
  !
  ! Open the file 
  !
  open(unit=1,file=filename,status='old')
  !
  ! Skip all lines at the beginning of file that start with one of the
  ! symbols ";", "#" or "!". These are meant as comments. 300 such lines
  ! is the maximum. 
  !
  do iline=1,300
     read(1,'(A158)',end=10) string
     if (string(1:1).ne.'#') goto 20
  enddo
  stop 8901
10 continue
  write(*,*) 'ERROR while reading file ',trim(filename)
  write(*,*) 'Did not find data before end of file.'
  stop 8902
20 continue
  !
  ! First read the format number. Do so from the above read string.
  !
  read(string,*) iformat
  if(iformat.eq.3) then
     scat = .false.
  else if(iformat.eq.0) then
     scat = .true.
  else
     write(*,*) 'ERROR: Format number ',iformat,' of file ',trim(filename),&
          ' not known.'
     stop
  endif
  !
  ! Read the number of wavelengths
  !
  read(1,*) nlam
  !
  ! Read the number of scattering angles
  !
  if (scat) then
     read(1,*) nang
  endif
  !
  ! Allocate arrays
  !
  allocate(lam(nlam),kabs(nlam),ksca(nlam),gsca(nlam))
  !
  ! Read the frequencywavelength and the absorption and scattering opacities.
  ! Also read the gfactor.
  !
  do ilam=1,nlam
     read(1,*) dummy4
     lam(ilam)  = dummy4(1) / 1d4
     kabs(ilam) = dummy4(2)
     ksca(ilam) = dummy4(3)
     gsca(ilam) = dummy4(4)
  enddo
  if (scat) then
     allocate(ang(nang))
     allocate(f11(nlam,nang),f12(nlam,nang),f22(nlam,nang))
     allocate(f33(nlam,nang),f34(nlam,nang),f44(nlam,nang))
     !
     ! Read the angular grid
     !
     do iang=1,nang
        read(1,*) dum
        ang(iang) = dum
     enddo
     !
     ! Read the full scattering matrix. Each line contains
     ! the 11, 12, 22, 33, 34, 44 elements of the matrix
     !
     do ilam=1,nlam
        do iang=1,nang
           read(1,*) dummy6
           f11(ilam,iang) = dummy6(1)
           f12(ilam,iang) = dummy6(2)
           f22(ilam,iang) = dummy6(3)
           f33(ilam,iang) = dummy6(4)
           f34(ilam,iang) = dummy6(5)
           f44(ilam,iang) = dummy6(6)
        enddo
     enddo
  endif
  !
  ! Close the file
  ! 
  close(1)
  !
  ! Deallocate arrays
  !
  print *,lam(1),kabs(1),ksca(1),gsca(1)
  print *,f11(1,1),f44(1,1)
  print *,f11(nlam,nang),f44(nlam,nang)
  deallocate(lam,kabs,ksca,gsca)
  if (scat) then
     deallocate(ang)
     deallocate(f11,f12,f22,f33,f34,f44)  
  endif
end program readasciifortran

