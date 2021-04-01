program ZalesakIC

implicit none

double precision :: x, y, r2, dx(2)
integer, parameter :: nfine=50
double precision, allocatable:: VOF(:,:,:)
integer :: i, j, k, ii, jj, nx, ny, nz


	!Create mesh

	nx = 50
	ny = 50
	nz = 1

	allocate(VOF(0:nx-1,0:ny-1,0:nz-1))

	VOF = 0.0d0

	dx(1) = 1.0/nx
	dx(2) = 1.0/ny

	open(unit=10,file='xvec.txt')	
	do i = 0, nx-1
		write(10,*)-0.5d0 + (i+0.5d0)*dx(1)
	end do
	close(10)

	open(unit=10,file='yvec.txt')	
	do j = 0, ny-1
		write(10,*)-0.5d0 + (j+0.5d0)*dx(2)
	end do
	close(10)

    open(unit=10,file='zalesakdiskIC.txt')

    ! limits of i, j, k to be lo(1), hi(1) etc. as it is in PeleLM
	do k=0,0 ! Since it is basically 2D
		do j= 0, ny-1
			y = -0.5d0 + j*dx(2)   ! nodal value of y. In PeleLM this will be float(j)*delta(2)+domnlo(2)
			do i= 0, nx-1
				x = -0.5d0 + i*dx(1)  ! nodal value of x. In PeleLM this will be float(i)*delta(1)+domnlo(1)
				! Get Zalesak distance on submesh
				do jj=1,nfine
					do ii=1,nfine
						if (init_zalesak((/x+(ii-1+0.5d0)*dx(1)/nfine,&
                  	&             y+(jj-1+0.5d0)*dx(2)/nfine,0.0d0/),&
                    &           (/0.0d0,0.25d0,0.0d0/),0.15d0,0.05d0,0.25d0).gt.0.0d0) then
                  		VOF(i,j,k)=VOF(i,j,k)+1.0d0/(nfine**2)
						end if
					end do
				end do
				write(10,*)i, j, VOF(i,j,k)
			end do
		end do
	end do

	close(10)

contains
   
  
  ! Zalesak disk initialization
  function init_zalesak(xyz,center,radius,width,height)
    implicit none
    double precision :: init_zalesak
    double precision, dimension(3), intent(in) :: xyz,center
    double precision, intent(in) :: radius,width,height
    double precision :: c,b,b1,b2,h1,h2
    c = radius-sqrt(sum((xyz-center)**2))
    b1 = center(1)-0.5d0*width
    b2 = center(1)+0.5d0*width
    h1 = center(2)-radius*cos(asin(0.5d0*width/radius))
    h2 = center(2)-radius+height
    if     (c>=0.0d0.and.xyz(1)<=b1.and.xyz(2)<=h2) then
       b = b1-xyz(1)
       init_zalesak = min(c,b)
    elseif (c>=0.0d0.and.xyz(1)>=b2.and.xyz(2)<=h2) then
       b = xyz(1)-b2
       init_zalesak = min(c,b)
    elseif (c>=0.0d0.and.xyz(1)>=b1.and.xyz(1)<=b2.and.xyz(2)>=h2) then
       b = xyz(2)-h2
       init_zalesak = min(c,b)
    elseif (c>=0.0d0.and.xyz(1)<=b1.and.xyz(2)>=h2) then
       b = sqrt(sum((xyz-(/b1,h2,0.0d0/))**2))
       init_zalesak = min(c,b)
    elseif (c>=0.0d0.and.xyz(1)>=b2.and.xyz(2)>=h2) then
       b = sqrt(sum((xyz-(/b2,h2,0.0d0/))**2))
       init_zalesak = min(c,b)
    elseif (xyz(1)>=b1.and.xyz(1)<=b2.and.xyz(2)<=h2.and.xyz(2)>=h1) then
       init_zalesak = -min(abs(xyz(1)-b1),abs(xyz(1)-b2),abs(xyz(2)-h2))
    elseif (xyz(1)>=b1.and.xyz(1)<=b2.and.xyz(2)<=h1) then
       init_zalesak = -min(sqrt(sum((xyz-(/b1,h1,0.0d0/))**2)),sqrt(sum((xyz-(/b2,h1,0.0d0/))**2)))
    else
       init_zalesak = c
    endif
  end function init_zalesak


end program ZalesakIC
