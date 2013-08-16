! This module is intended to hold the variables corresponding to
!  the time and spatial dimensions

      module wrf_dims_vars

      implicit none

      integer :: mix,mjx,mkx
      integer :: mixic,mjxic,mkxic
! 3D
      real, allocatable :: u(:,:,:),v(:,:,:),w(:,:,:),
     & ph(:,:,:),phb(:,:,:),gpttot(:,:,:),gph(:,:,:),
     & gphhalf(:,:,:),qvapor(:,:,:),theta(:,:,:),
     & p(:,:,:),pb(:,:,:),rhodry(:,:,:),
     & totpres(:,:,:),temp(:,:,:),qrain(:,:,:),
     & qsnow(:,:,:),qgraup(:,:,:)
! 2D
      real, allocatable :: mu(:,:),mub(:,:),surfptotal(:,:),
     & q2(:,:),t2(:,:),th2(:,:),psfc(:,:),slp(:,:),
     & hgt(:,:),u10(:,:),v10(:,:),windspeed_100m(:,:),
     & windspeed_200m(:,:),rainc(:,:),rainnc(:,:)
! 1D
      real, allocatable :: znw(:),znu(:)
! Non-WRF
! x,y,z,t,var(,mem)
      real, allocatable :: initial3d(:,:,:,:,:,:),
     & initmean3d(:,:,:,:,:),
     & variance3d(:,:,:,:,:),covar3d(:,:,:,:,:),sens3d(:,:,:,:,:),
     & corr3d(:,:,:,:,:),stdsens3d(:,:,:,:,:),tstat3d(:,:,:,:,:)
! x,y,t,var(,mem)
      real, allocatable :: initial2d(:,:,:,:,:),initmean2d(:,:,:,:),
     & variance2d(:,:,:,:),covar2d(:,:,:,:),sens2d(:,:,:,:),
     & corr2d(:,:,:,:),stdsens2d(:,:,:,:),tstat2d(:,:,:,:)
      real, allocatable :: response(:,:,:,:),rp(:)
      real :: rpmean,variance_rp

      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GET_DIMENSIONS(infile,permissr,iunit1,mix,mjx,mkx)
      ! subroutine retrieves x,y,z dimensions sizes from WRF file
      use netcdf
      implicit none
      character(len=100) :: infile
      integer :: rcode,permissr,iunit1,mix,mjx,mkx
      call open_file(infile,permissr,iunit1)
      rcode=nf_get_att_int(iunit1,nf_global,
     &           'WEST-EAST_GRID_DIMENSION',mix)
      rcode=nf_get_att_int(iunit1,nf_global,
     &           'SOUTH-NORTH_GRID_DIMENSION',mjx)
      rcode=nf_get_att_int(iunit1,nf_global,
     &           'BOTTOM-TOP_GRID_DIMENSION',mkx)
      call close_file(iunit1)
      END SUBROUTINE GET_DIMENSIONS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module wrf_dims_vars
