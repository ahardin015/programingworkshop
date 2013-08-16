      program vort_calc 
      USE NETCDF
      USE WRF_TOOLS
      use wrf_dims_vars
      use conversions_mod
      implicit none
      character(len=100) :: infile,outfile ! name of in/output files
      integer :: st_time,end_time,ntimes
      real, allocatable :: vort(:,:),ws(:,:,:)
      real, allocatable :: uwind10(:,:),vwind10(:,:)
      logical :: debug,msg ! flags for message printing
      integer :: iunit=10,iunit1,iunit2 ! unit numbers for namelist,infiles and outfile
      integer :: permissw=3,permissr=2 ! write/read permissions
      integer :: rcode,ncid ! integers used for writing netcdf 
      integer :: tdim,ydim,xdim,varid
      integer :: t,tmod,k
      real, parameter :: miss=-9999.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      namelist /vort_input/ infile,outfile,st_time,end_time,msg,debug
      open(iunit,file='vort_input.nml', status='old')
      read(iunit, nml = vort_input)
      close(iunit)

!-----------------------------------------------------------------------

! IC CALCULATIONS
      if (msg) print*, 'Doing ICs'
      call get_dimensions(infile,permissr,iunit1,mix,mjx,mkx)
      if (msg) print*, '  WRF Domain dimensions of file', mix,mjx,mkx
      ntimes=end_time-st_time+1
      if (msg) print*, 'Doing ', ntimes, 'total times'
!      allocate(u10(mix-1,mjx-1))
!      allocate(v10(mix-1,mjx-1))
      allocate(ws(mix-1,mjx-1,ntimes))
      if (msg) print*, '  Arrays allocated'
! READ DATA
      do t=1,ntimes
         tmod = st_time+t-1
         allocate(u10(mix-1,mjx-1))
         allocate(v10(mix-1,mjx-1))
         if (msg) print*, 'Doing time ',t,'index',tmod
         call open_file(infile,permissr,iunit1)
         call get_variable2d(iunit1,'U10',mix-1,mjx-1,tmod,u10)
         call get_variable2d(iunit1,'V10',mix-1,mjx-1,tmod,v10)
         call close_file(iunit1)
         if (msg) print*, '   Finished reading time',tmod
! WIND CALCULATIONS
         ws(:,:,t)=sqrt(u10**2+v10**2)
         deallocate(u10)
         deallocate(v10)
      enddo
      print*,size(ws)
        ! J loop then I loop to find vorticity
      if (msg) print*, 'calcualtions done'

! WRITE OUTPUT

      ! Create the netcdf file with output name
      rcode=nf_create(outfile,NF_CLOBBER,ncid)      
      if (debug) print*, 'Created output file', outfile
      !  Define dimensions
      rcode=nf_def_dim(ncid,'Time',ntimes,tdim)
      rcode=nf_def_dim(ncid,'south_north',mjx-1,ydim)
      rcode=nf_def_dim(ncid,'west_east',mix-1,xdim)
      if (debug) print*, 'Dimensions defined'
      
      rcode=nf_def_var(ncid,'WS10',NF_FLOAT,3,(/xdim,ydim,tdim/),varid)
      rcode=nf_put_att_text(ncid,varid,'description',13,'10m_windspeed')
      rcode=nf_put_att_text(ncid,varid,'units',3,'m/s')
!      rcode=nf_put_att_real(ncid,varid,'_FillValue',NF_FLOAT,1,miss)
      if (debug) print*, 'output variables defined'



      do k=1,ntimes
         if (debug) print*, 'writing data for index', k
         call write_variable2d(ncid,'WS10',mix-1,mjx-1,k,ws(:,:,k))
      enddo
      call close_file(ncid)


      if (msg) print*,'Done with main sensitivity program'
      END PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




