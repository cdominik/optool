program readoptoolfits
  implicit none
  character*500 input
  character*80 comment,errmessage
  character*30 errtext
  integer*4 :: status,unit,blocksize,nfound,group,readwrite
  integer*4 :: firstpix,npixels,hdutype
  real*8  :: nullval,rho_gr,a1,a2,a3
  logical*4 :: anynull
  integer*4, dimension(4) :: naxes

  real*8,allocatable :: array(:,:),matrix(:,:,:)
  real*8,allocatable :: lam(:),kabs(:),ksca(:),kext(:)
  real*8,allocatable :: ang(:),F11(:,:),F12(:,:),F22(:,:),F33(:,:),F34(:,:),F44(:,:)

  integer i,j,ia,iT,iread,nlam,nang,isize
  real*8 l0,l1,tot,tot2,theta,asym,Pmax,HG,asym2,wasym2
  real rho_av

  ! Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  input = "dustkappa.fits"

  call ftgiou (unit,status)
  ! Open file
  readwrite=0
  call ftopen(unit,input,readwrite,blocksize,status)
  if (status /= 0) then
     print *,"Error reading particle file: ", trim(input)
     stop
  endif
  group=1
  firstpix=1
  nullval=-999
  
  !------------------------------------------------------------------------
  ! HDU0 : opacities
  !------------------------------------------------------------------------
  ! Check dimensions
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
  
  nlam = naxes(1)
  
  npixels=naxes(1)*naxes(2)
  
  ! Read model info
  
!  call ftgkye(unit,'a1',a1,comment,status)
!  call ftgkye(unit,'a2',a2,comment,status)
!  call ftgkye(unit,'a3',a3,comment,status)
!  call ftgkye(unit,'density',rho_av,comment,status)
  ! read_image
  allocate(array(nlam,4))
  
  call ftgpvd(unit,group,firstpix,npixels,nullval,array,anynull,status)

  !------------------------------------------------------------------------
  ! HDU 1: matrix
  !------------------------------------------------------------------------
  !  move to next hdu
  call ftmrhd(unit,1,hdutype,status)	
  
  ! Check dimensions
  call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
  
  npixels=naxes(1)*naxes(2)*naxes(3)
  nang = naxes(3)
  
  ! read_image
  allocate(matrix(nlam,6,180))
  
  call ftgpvd(unit,group,firstpix,npixels,nullval,matrix,anynull,status)
    
  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)
  
  !  Check for any error, and if so print out error messages
  !  Get the text string which describes the error
  if (status > 0) then
     call ftgerr(status,errtext)
     print *,"ERROR: error in reading fits file",errtext
     
     !  Read and print out all the error messages on the FITSIO stack
     call ftgmsg(errmessage)
     do while (errmessage .ne. ' ')
        print *,errmessage
        call ftgmsg(errmessage)
     end do
  endif

  !------------------------------------------------------------------------
  ! Copy the read-in values into the proper arrays
  !------------------------------------------------------------------------

  allocate(lam(nlam))
  allocate(Kabs(nlam))
  allocate(Ksca(nlam))
  allocate(Kext(nlam))
  allocate(ang(nang))

  lam(1:nlam)  = array(1:nlam,1)
  kext(1:nlam) = array(1:nlam,2)
  kabs(1:nlam) = array(1:nlam,3)
  ksca(1:nlam) = array(1:nlam,4)

  print *,nlam
  allocate(F11(nlam,nang))
  allocate(F12(nlam,nang))
  allocate(F22(nlam,nang))
  allocate(F33(nlam,nang))
  allocate(F34(nlam,nang))
  allocate(F44(nlam,nang))

  
  F11(1:nlam,1:nang)  = matrix(1:nlam,1,1:nang)
  F12(1:nlam,1:nang)  = matrix(1:nlam,2,1:nang)
  F22(1:nlam,1:nang)  = matrix(1:nlam,3,1:nang)
  F33(1:nlam,1:nang)  = matrix(1:nlam,4,1:nang)
  F34(1:nlam,1:nang)  = matrix(1:nlam,5,1:nang)
  F44(1:nlam,1:nang)  = matrix(1:nlam,6,1:nang)

  deallocate(array)
  deallocate(matrix)
  deallocate(lam,kext,kabs,ksca)
  deallocate(F11,F12,F22,F33,F34,F44)
  
end program readoptoolfits



