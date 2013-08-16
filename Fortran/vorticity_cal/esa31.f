      program v31_ESA 
      USE NETCDF
      USE WRF_TOOLS
      use wrf_dims_vars
      use conversions_mod
      implicit none
      integer :: varid
      integer :: ensnum ! number of ensemble members
      character(len=50) :: infilebaseinit,infilebaseresp ! base name of infiles; ens. # on end
      integer,parameter :: max_initvars=5,max_initlevels=10
      integer :: num_initvars
      integer :: num_initlevels
      integer :: num_ictimes
      character(len=10) :: initvar2d(max_initvars)
      character(len=10) :: initvar3d(max_initvars)
      character(len=10) :: initvar(2,max_initvars)
      character(len=40) :: initdescrip2d(max_initvars)
      character(len=40) :: initdescrip3d(max_initvars)
      character(len=40) :: initdescrip(2,max_initvars)
      character(len=5) :: initunits2d(max_initvars)
      character(len=5) :: initunits3d(max_initvars)
      character(len=5)  :: initunits(2,max_initvars)
      character(len=10) :: ic_ftype,rp_ftype
      character(len=100) :: outfile
      character(len=70) :: outvar2d,outvar3d
      integer :: timeinitial,tindex ! forecast hour and time index of initial variable
      integer :: starttime,endtime ! fhr (index) times
      integer :: plt_int ! Output interval of the sens file
      integer :: model_files_per_hour ! Output interval (num per hour) of wrf model
      real :: initlevels(max_initlevels) ! vertical levels of initial
      character(len=1) :: initvertunit,respvertunit ! vertical unit for interpolation
      integer :: num_rpvars
      character(len=10) :: resp_type, responsevar ! response variable name
      character(len=40) :: respdescrip
      character(len=5) :: respunits
      integer :: t_resp_st, t_resp_end ! WRF time index of response
      integer :: num_rptimes ! Time window over which rp funct is found
      real :: responselevel1,responselevel2 ! vertical levels of response
      integer :: respx,respy ! WRF coordinates of point/center of area for response
      integer :: dx,dy ! dimensions of area for response: (0,0) for a point
!      character(len=10) :: avgminmax ! avg, min, or max of resp area
      logical :: scatter
      integer :: scatterx,scattery
      logical :: debug,msg ! flags for message printing
      integer :: permissw=3,permissr=2 ! write/read permissions
      integer :: iunit=10,iunit1,iunit2 ! unit numbers for namelist,infiles and outfile
      integer :: i,j,k,pass,mem,func,arrx,arry! loop counters
      character(len=100) :: infile
      character(len=200) :: outvar
      character(len=4) :: en,fhr ! string with ensemble member number
      integer :: rcode ! netcdf error flag
      real :: rpfunction
      character(len=100) :: infile_name
      real,parameter :: miss=-9999. ! flag that lets us know when a variable is missing; typically
                             ! will be used when an interpolation has an underground result

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      namelist /sens_input/ ensnum,infilebaseinit,infilebaseresp, 
     & ic_ftype,rp_ftype,
     & outfile,starttime,endtime,plt_int,model_files_per_hour,
     & num_initvars,
     & initvar3d,initdescrip3d,initunits3d,initvar2d,initdescrip2d,
     & initunits2d,num_initlevels,initlevels,initvertunit,
     & resp_type,t_resp_st,t_resp_end,responsevar,respdescrip,
     & respunits,responselevel1,responselevel2,respvertunit,respx,
     & respy,dx,dy,scatter,scatterx,scattery,msg,debug
      open(iunit,file='sens_input.nml', status='old')
      read(iunit, nml = sens_input)
      close(iunit)

      initvar(1,:)=initvar2d
      initvar(2,:)=initvar3d
      initdescrip(1,:)=initdescrip2d
      initdescrip(2,:)=initdescrip3d
      initunits(1,:)=initunits2d
      initunits(2,:)=initunits3d
!-----------------------------------------------------------------------

! IC CALCULATIONS
      if (msg) print*, 'Doing ICs'
      infile=infile_name(infilebaseinit,1)
      call get_dimensions(infile,permissr,iunit1,mix,mjx,mkx)
      if (msg) print*, '  WRF Domain dimensions of file', mix,mjx,mkx
      num_ictimes=(endtime-starttime)/plt_int+1
      if (msg) print*,'num_ictimes',num_ictimes
      if (msg) print*,'IC vars ', initvar
      call allocate_arrays(1,num_ictimes,num_initvars,
     &    ensnum,num_initlevels)
      initmean2d=0.
      initmean3d=0.
      do mem=1,ensnum ! loop over ensemble
        if (msg) print*,'Working on member',mem
        infile=infile_name(infilebaseinit,mem)
        if (debug) print*,'  ',trim(infile)
!       do range of times from namelist
        do timeinitial=starttime,endtime,plt_int 
         tindex=(timeinitial-starttime)/plt_int+1
!        READ DATA
          call get_variables(infile,permissr,iunit1,timeinitial+1,
     & ic_ftype)
          call variable_processing()
          do i=1,num_initvars
            ! Surface vars
            call set_function(initvar2d(i),initial2d(i,tindex,:,:,
     &            mem),initlevels(1),initlevels(2),initvertunit)
            ! Vertical levels
            do k=1,num_initlevels
              call set_function(initvar3d(i),initial3d(i,tindex,
     &         :,:,k,mem),initlevels(k),initlevels(1),initvertunit)
            enddo ! vert levels
          enddo ! IC vars
        enddo ! times
      if (debug) then
        print*,'initial2d',initial2d(:,1,1,1,mem)
        print*,'initial3d',initial3d(:,1,1,1,1,mem)
      endif
      ! Do sum
      where (initial3d(:,:,:,:,:,mem) .ne. miss .and. 
     &       initmean3d .ne. miss)
        initmean3d=initmean3d+initial3d(:,:,:,:,:,mem)
       else where
        initmean3d=miss
      end where
      if (debug) print*
      enddo ! mem
      call deallocate_arrays(2) ! WRF vars
      if (msg) print*
      initmean2d=sum(initial2d,5)/ensnum
      where (initmean3d .ne. miss)
        initmean3d=initmean3d/ensnum
      end where
      if (debug) then
        print*,'initmean2d',initmean2d(:,1,1,1)
        print*,'initmean3d',initmean3d(:,1,1,1,1)
      endif

      if (msg) then
        print*, 'Doing IC variance'
        print*
      endif
      call ic_variance(ensnum)
      mixic=mix
      mjxic=mjx
!-------------------------------------------------------------------------------
! RESPONSE FUNCTION CALCULATIONS
      if (msg) then
        print*,'Doing response: '
        print*,'   ',trim(resp_type),' of ',responsevar
        print*,'   Response window  ',t_resp_st, t_resp_end
        print*,'   gridpt',respx,respy,'+-',dx,dy
      endif
!   GET DIMENSION
      infile=infile_name(infilebaseresp,1)
      if (debug) print*,trim(infile)
      call get_dimensions(infile,permissr,iunit1,mix,mjx,mkx)
      if (msg) print*, '  WRF Domain dimensions of file', mix,mjx,mkx
! CALCULATE RESPONSE TIME WINDOW- between t_resp_st and t_resp_end
      num_rptimes=(t_resp_end-t_resp_st)+1
      if (debug) print*, 'num_rptimes',num_rptimes
!   ALLOCATE SPACE
      if (debug) print*,'Allocating arrays'
      call allocate_arrays(2,num_rptimes,1,ensnum,num_initlevels)
      rpmean=0.
      variance_rp=0.
!   READ THE DATA
      do mem=1,ensnum ! pass 1: mean
        if (msg) print*,'Doing member',mem
        if (debug) print*,'  Reading data'
        infile=infile_name(infilebaseresp,mem)
        if (debug) print*,'  ',trim(infile)
        do timeinitial=t_resp_st,t_resp_end,1
          tindex=(timeinitial-t_resp_st)/1+1
          if(debug) print*, 'Doing time', timeinitial, tindex
          call get_variables(infile,permissr,iunit1,timeinitial,
     & rp_ftype)
          call variable_processing() ! Do some variable conversions
          call set_function(responsevar,response(tindex,:,:,mem),
     &    responselevel1,responselevel2,respvertunit)
          if (debug) print*,'  resp', response(tindex,respx,respy,mem)
        enddo ! loop over resp time window
        if (debug) print*, '       Got all data, now doing RESP'
        rp(mem)=rpfunction(response(:,respx,respy,mem),num_rptimes,
     &   resp_type)
          if (debug) print*,'Member',mem
          if (debug) print*,'time series',response(:,respx,respy,mem)
          if (debug) print*,'rp value',rp(mem)
          if (debug) print*,'maxval',maxval(response(:,respx,respy,mem))
          if (debug) print*,'maxloc',maxloc(response(:,respx,respy,mem))
          if (debug) print*,'  '
      enddo ! loop over mems
      rpmean=sum(rp)/ensnum
      call deallocate_arrays(1)
      print*, ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Revert back to IC dimensions
      mix=mixic
      mjx=mjxic
      call allocate_arrays(3,num_ictimes,num_initvars,
     &    ensnum,num_initlevels)
      if (debug) print*,'Doing ensemble statistics'
      call ensemble_stats(ensnum)
      if (msg) then
        print*,'Response function mean and variance'
        print*,rpmean,variance_rp
        print*,'IC mean'
        print*,initmean2d(:,1,1,1)
        print*,'IC variance'
        print*,variance2d(:,1,1,1)
        print*,'Sensitivity'
        print*, sens2d(:,1,1,1)
        print*,'Standardized Sensitivity'
        print*,stdsens2d(:,1,1,1)
        print*,'Correlation coefficient'
        print*,corr2d(:,1,1,1)
        print*,'T Stat'
        print*,tstat2d(:,1,1,1)
      endif

! WRITE OUTPUT
      call create_sens_netcdf(outfile,infilebaseinit,
     &  responsevar,resp_type,
     &  respdescrip,respunits,num_initvars,initvar,
     &  initdescrip,initunits,respx,respy,dx,dy,t_resp_st,
     &  t_resp_end,
     &  starttime,endtime,plt_int,
     &  model_files_per_hour,num_initlevels,initlevels,
     &  initvertunit,respvertunit)
      if (msg) then
        print*
        print*, 'Editing output file: ', trim(outfile)
      endif
      call open_file(outfile,permissw,iunit2)
      do j=1,num_initvars
      do timeinitial=starttime,endtime,plt_int
        tindex=(timeinitial-starttime)/plt_int+1
      if (debug) then
        print*,'fhr',timeinitial
        print*,'time index',tindex
        print*,'vars         ',initvar(:,j)
      endif

      if (.true.) then
      outvar='SENS_'//trim(resp_type)//'_'//trim(responsevar)
     &       //'_'//trim(initvar2d(j))
      call write_variable2d(iunit2,trim(outvar),mix-1,mjx-1,tindex,
     &      sens2d(j,tindex,:,:))
      outvar='SENS_'//trim(resp_type)//'_'//trim(responsevar)
     &       //'_'//trim(initvar3d(j))
      call write_variable3d(iunit2,trim(outvar),mix-1,mjx-1,
     &      num_initlevels,tindex,sens3d(j,tindex,:,:,:))
      outvar='STDSENS_'//trim(resp_type)//'_'//trim(responsevar)
     &       //'_'//trim(initvar2d(j))
      call write_variable2d(iunit2,trim(outvar),mix-1,mjx-1,tindex,
     &      stdsens2d(j,tindex,:,:))
      outvar='STDSENS_'//trim(resp_type)//'_'//trim(responsevar)
     &       //'_'//trim(initvar3d(j))
      call write_variable3d(iunit2,trim(outvar),mix-1,mjx-1,
     &      num_initlevels,tindex,stdsens3d(j,tindex,:,:,:))
      outvar='MEAN_'//trim(initvar2d(j))
      call write_variable2d(iunit2,trim(outvar),mix-1,mjx-1,tindex,
     &      initmean2d(j,tindex,:,:))
      outvar='MEAN_'//trim(initvar3d(j))
      call write_variable3d(iunit2,trim(outvar),mix-1,mjx-1,
     &      num_initlevels,tindex,initmean3d(j,tindex,:,:,:))
      outvar='VARIANCE_'//trim(initvar2d(j))
      call write_variable2d(iunit2,trim(outvar),mix-1,mjx-1,tindex,
     &      variance2d(j,tindex,:,:))
      outvar='VARIANCE_'//trim(initvar3d(j))
      call write_variable3d(iunit2,trim(outvar),mix-1,mjx-1,
     &      num_initlevels,tindex,variance3d(j,tindex,:,:,:))
      outvar='COVARIANCE_'//trim(responsevar)//'_'//trim(initvar2d(j))
      call write_variable2d(iunit2,trim(outvar),mix-1,mjx-1,tindex,
     &      covar2d(j,tindex,:,:))
      outvar='COVARIANCE_'//trim(responsevar)//'_'//trim(initvar3d(j))
      call write_variable3d(iunit2,trim(outvar),mix-1,mjx-1,
     &      num_initlevels,tindex,covar3d(j,tindex,:,:,:))

      endif
      enddo ! times
      enddo ! ic vars
      call close_file(iunit2)

      if (scatter) then
        call write_scatter_data(ensnum,initvar,initial2d(1,1,scatterx,
     &        scattery,:),initial3d(1,1,scatterx,scattery,1,:),
     &        responsevar,rp,resp_type)
      endif

      call deallocate_arrays(3)
      call deallocate_arrays(4)
      if (msg) print*,'Done with main sensitivity program'
      END PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE WRITE_SCATTER_DATA(ensnum,icvar,
     &   ic2d,ic3d,rpvar,rp,resp_type)
      implicit none
      integer :: ensnum,mem,d
      integer,parameter :: iunit=10
      character(len=10) :: icvar(2,1),rpvar,resp_type
      real :: ic2d(ensnum),ic3d(ensnum),ic(2,ensnum),rp(ensnum)
      ic(1,:)=ic2d
      ic(2,:)=ic3d
      open(iunit,file='scatter.txt',status='unknown')
      do d=1,2
        write(iunit,*) icvar(d,1),',',resp_type,' of ',rpvar
        do mem=1,ensnum
          write(iunit,*) ic(d,mem),',',rp(mem)
        enddo
      enddo
      close(iunit)
      END SUBROUTINE WRITE_SCATTER_DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE CREATE_SENS_NETCDF(esensfile,infilebase,
     &   rpvar,rptype,
     &   rpdescrip,rpunits,num_icvars,icvars,icdescrip,icunits,
     &   respx,respy,dx,dy,rpstart,rpend,
     &   starttime,endtime,plt_int, model_files_per_hour,
     &   zsize,levels,icvertunit,rpvertunit)
      use netcdf
      use wrf_dims_vars
      implicit none
      character(len=60) :: esensfile
      character(len=20) :: wrftemplate,infilebase,infile_name
      integer :: zsize,numtimes
      real, parameter :: miss=-9999.0
      integer :: rcode,ncid,varid
      integer :: xdim,ydim,zdim,tdim
      integer,allocatable :: dims(:)
      real :: levels(zsize),lon(mix-1,mjx-1),lat(mix-1,mjx-1)
      integer :: respx,respy,dx,dy,rpstart,rpend
      integer :: starttime,endtime,plt_int,model_files_per_hour
      integer :: num_icvars
      character(len=10) :: rpvar,rptype,icvars(2,num_icvars)
      character(len=50) :: varname
      character(len=40) :: icdescrip(2,num_icvars)
      character(len=40) :: rpdescrip
      character(len=100) :: descrip
      character(len=5) :: rpunits,icunits(2,num_icvars)
      character(len=10)  :: units
      character(len=1) :: icvertunit,rpvertunit
      integer :: d,i,j,k
      numtimes=(endtime-starttime)/plt_int+1
      wrftemplate=infile_name(infilebase,rpstart,1)
      ! Create netcdf file. Overwrite if it already exists.
      rcode=nf_create(esensfile,NF_CLOBBER,ncid)
      ! Define dimensions
      rcode=nf_def_dim(ncid,'Time',numtimes,tdim)
      rcode=nf_def_dim(ncid,'bottom_top',zsize,zdim)
      rcode=nf_def_dim(ncid,'south_north',mjx-1,ydim)
      rcode=nf_def_dim(ncid,'west_east',mix-1,xdim)
      rcode=nf_def_var(ncid,'LEVELS',NF_FLOAT,1,zdim,varid)
      rcode=nf_put_att_text(ncid,varid,'description',15,
     &     'Vertical levels')
      rcode=nf_put_att_text(ncid,varid,'vertunit',1,icvertunit)
      rcode=nf_put_att_real(ncid,varid,'_FillValue',NF_FLOAT,1,miss)
      rcode=nf_def_var(ncid,'LON',NF_FLOAT,2,(/xdim,ydim/),varid)
      rcode=nf_put_att_text(ncid,varid,'description',9,'Longitude')
      rcode=nf_put_att_text(ncid,varid,'units',3,'deg')
      rcode=nf_put_att_real(ncid,varid,'_FillValue',NF_FLOAT,1,miss)
      rcode=nf_def_var(ncid,'LAT',NF_FLOAT,2,(/xdim,ydim/),varid)
      rcode=nf_put_att_text(ncid,varid,'description',8,'Latitude')
      rcode=nf_put_att_text(ncid,varid,'units',3,'deg')
      rcode=nf_put_att_real(ncid,varid,'_FillValue',NF_FLOAT,1,miss)
      do d=2,3
        allocate(dims(d+1))
        dims(1)=xdim
        dims(2)=ydim
        if (d==2) then
          dims(3)=tdim
         else if (d==3) then
          dims(3)=zdim
          dims(4)=tdim
        endif
        do i=1,num_icvars
          do k=1,7
            if (k==1) then ! Raw Sens
              varname='SENS_'//trim(rptype)//'_'//trim(rpvar)//'_'//
     &               trim(icvars(d-1,i))
              descrip='Sensitivity of '//trim(rptype)//' of '
     &              //trim(rpdescrip)//' to '//trim(icdescrip(d-1,i))
              units=trim(rpunits)//'/'//trim(icunits(d-1,i))
             else if (k==2) then ! Std Sens
              varname='STDSENS_'//trim(rptype)//'_'//trim(rpvar)//'_'//
     &           trim(icvars(d-1,i))
              descrip='Standardized Sensitivity of '//
     &           trim(rptype)//' of '//trim(rpdescrip)//
     &           ' to '//trim(icdescrip(d-1,i))
              units=trim(rpunits)
             else if (k==3) then ! IC Mean
              varname='MEAN_'//trim(icvars(d-1,i))
              descrip='Mean '//trim(icdescrip(d-1,i))
              units=trim(icunits(d-1,i))
             else if (k==4) then ! IC Variance
              varname='VARIANCE_'//trim(icvars(d-1,i))
              descrip='Variance of '//trim(icdescrip(d-1,i))
              units=trim(icunits(d-1,i))//'^2'
             else if (k==5) then ! Correlation Coefficient
              varname='CORR_'//trim(rpvar)//'_'//trim(icvars(d-1,i))
              descrip='Correlation of '//trim(rpdescrip)//' and '//
     &               trim(icdescrip(d-1,i))
              units=' '
             else if (k==6) then ! T Stat
              varname='TSTAT_'//trim(rpvar)//'_'//trim(icvars(d-1,i))
              descrip='T-statistic for sensitivity of '//trim(rpdescrip)
     &               //' to '//trim(icdescrip(d-1,i))
              units=' '
             else if (k==7) then ! Covariance
              varname='COVARIANCE_'//trim(rpvar)//'_'//
     &           trim(icvars(d-1,i))
              descrip='Covariance of '//trim(rpdescrip)//' to '//
     &           trim(icdescrip(d-1,i))
              units=' '
            endif
            rcode=nf_def_var(ncid,trim(varname),NF_FLOAT,d+1,dims,varid)
            rcode=nf_put_att_text(ncid,varid,'description',
     &           len_trim(descrip),descrip)
            rcode=nf_put_att_text(ncid,varid,'units',
     &           len_trim(units),units)
            rcode=nf_put_att_real(ncid,varid,'_FillValue',
     &           NF_FLOAT,1,miss)
          enddo
        enddo
      deallocate(dims)
      enddo
      rcode=nf_put_att_text(ncid,NF_GLOBAL,'TITLE',20,
     &     'Ensemble Sensitivity')
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'RP_START',NF_INT,1,rpstart)
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'RP_END',NF_INT,1,rpend)
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'IC_START_TIME',NF_INT,
     &     1,starttime)
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'IC_END_TIME',NF_INT,
     &     1,endtime)
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'SENS_PLOT_INT',NF_INT,
     &     1,plt_int)
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'MODEL_FILES_PER_HOUR',NF_INT,
     &     1,model_files_per_hour)
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'RESPX',NF_INT,1,respx)
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'RESPY',NF_INT,1,respy)
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'DX',NF_INT,1,dx)
      rcode=nf_put_att_int(ncid,NF_GLOBAL,'DY',NF_INT,1,dy)
      rcode=nf_enddef(ncid)
      rcode=nf_inq_varid(ncid,'LEVELS',varid)
      rcode=nf_put_vara_real(ncid,varid,1,zsize,levels)
      rcode=nf_close(ncid)
      ! Now get stuff from pre-existing WRF file
      rcode=nf_open(wrftemplate,NF_NOWRITE,ncid)
       rcode=nf_inq_varid(ncid,'XLAT',varid)
       rcode=nf_get_vara_real(ncid,varid,(/1,1,1/),
     &     (/mix-1,mjx-1,1/),lat)
       rcode=nf_inq_varid(ncid,'XLONG',varid)
       rcode=nf_get_vara_real(ncid,varid,(/1,1,1/),
     &     (/mix-1,mjx-1,1/),lon)
      rcode=nf_close(ncid)
      rcode=nf_open(esensfile,NF_WRITE,ncid)
       rcode=nf_inq_varid(ncid,'LAT',varid)
       rcode=nf_put_vara_real(ncid,varid,(/1,1/),(/mix-1,mjx-1/),lat)
       rcode=nf_inq_varid(ncid,'LON',varid)
       rcode=nf_put_vara_real(ncid,varid,(/1,1/),(/mix-1,mjx-1/),lon)
      rcode=nf_close(ncid)
      END SUBROUTINE CREATE_SENS_NETCDF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION INFILE_NAME(infilebase,mem)
      ! This function creates a string of the form
      ! "infilebase_f{fcsthr}.{mem}"
      implicit none
      character(len=50) :: infilebase
      character(len=5) :: mem_str
      character(len=100) :: infile_name
      integer :: mem
      write(unit=mem_str, fmt='(I4)') mem
      mem_str=adjustl(mem_str)
      infile_name=trim(infilebase)//'.'//trim(mem_str)
      END FUNCTION INFILE_NAME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE IC_VARIANCE(ensnum)
      ! subroutine that calculates IC variance for 2D and 3D
      ! inputs: ensemble size
      ! output: variance
      use wrf_dims_vars
      implicit none
      integer :: ensnum,mem
      real, parameter :: miss=-9999.
      variance2d=0.
      variance3d=0.
      do mem=1,ensnum
        variance2d=variance2d+(initial2d(:,:,:,:,mem)-initmean2d)**2
        where (initmean3d .ne. miss)
          variance3d=variance3d+(initial3d(:,:,:,:,:,mem)-initmean3d)**2
         else where
          variance3d=miss
        end where
      enddo
      variance2d=variance2d/(ensnum-1.)
      where (variance3d .ne. miss)
        variance3d=variance3d/(ensnum-1.)
      end where
      END SUBROUTINE IC_VARIANCE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ENSEMBLE_STATS(ensnum)
      ! subroutine that calculates the covariances for a specific response function
      ! inputs: index for specific response function, ensemble size
      ! outputs: covariances
      use wrf_dims_vars
      implicit none
      integer :: rpvar,ensnum,mem
      real, parameter :: miss=-9999.
      covar2d=0.
      covar3d=0.
      corr2d=0.
      sens2d=0.
      stdsens2d=0.
      tstat2d=0.
      corr3d=0.
      sens3d=0.
      stdsens3d=0.
      tstat3d=0.

      do mem=1,ensnum
       ! Covariance
        covar2d=covar2d+(initial2d(:,:,:,:,mem)-initmean2d)*
     &                   (rp(mem)-rpmean)
        where (initmean3d .ne. miss)
          covar3d=covar3d+(initial3d(:,:,:,:,:,mem)-initmean3d)*
     &                     (rp(mem)-rpmean)
         else where
          covar3d=miss
        end where
      enddo ! mem
      covar2d=covar2d/(ensnum-1.)
      where (covar3d .ne. miss)
        covar3d=covar3d/(ensnum-1.)
      endwhere

      print*,'rp(:): ',rp
      print*,'rpmean: ',rpmean
      variance_rp=sum((rp(:)-rpmean)**2)/(ensnum-1.)
      print*,'Variance of response function: ',variance_rp
      print*, ' '

      ! Correlation coefficient, Sensitivity, Standardized Sensitivity
      corr2d=covar2d/sqrt(variance2d*variance_rp)
      sens2d=covar2d/variance2d
      stdsens2d=covar2d/sqrt(variance2d)
      where (variance3d .ne. miss)
        corr3d=covar3d/sqrt(variance3d*variance_rp)
        sens3d=covar3d/variance3d
        stdsens3d=covar3d/sqrt(variance3d)
!        print*,' did SENS CALCULATIONS'
!        print*, maxval(variance3d), 'variance'
!        print*, maxval(covar3d), 'covar'
!        print*, maxval(sens3d), 'sens'
!        print*, minval(variance3d), 'variance'
!        print*, minval(covar3d), 'covar'
!        print*, minval(sens3d), 'sens'
       else where
        corr3d=miss
        sens3d=miss
        stdsens3d=miss
      end where
      do mem=1,ensnum
        tstat2d=tstat2d+(rp(mem)-rpmean-sens2d*
     &          (initial2d(:,:,:,:,mem)-initmean2d))**2
        where (sens3d .ne. miss)
          tstat3d=tstat3d+(rp(mem)-rpmean-sens3d*
     &          (initial3d(:,:,:,:,:,mem)-initmean3d))**2
         else where
          tstat3d=-9999.
        end where
      enddo
      tstat2d=sens2d/sqrt(tstat2d/(ensnum-2.)/variance2d/(ensnum-1.))
      where (tstat3d .ne. miss)
        tstat3d=sens3d/sqrt(tstat3d/(ensnum-2.)/variance3d/(ensnum-1.))
      end where

      END SUBROUTINE ENSEMBLE_STATS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ALLOCATE_ARRAYS(func,times,vars,ensnum,levels)
      ! subroutine intended to allocate arrays with
      !  WRF dimensions mix,mjx,mkx, and the appropriate
      !  number of times and ensemble members
      use wrf_dims_vars
      implicit none
      integer :: func,times,vars,ensnum,levels
      if (func==1 .or. func==2) then
! 1D
      allocate(znw(mkx))
      allocate(znu(mkx-1))
! 3D
      allocate(u(mix,mjx-1,mkx-1))
      allocate(v(mix-1,mjx,mkx-1))
      allocate(w(mix-1,mjx-1,mkx-1))
      allocate(ph(mix-1,mjx-1,mkx))
      allocate(phb(mix-1,mjx-1,mkx))
      allocate(gpttot(mix-1,mjx-1,mkx))
      allocate(gph(mix-1,mjx-1,mkx))
      allocate(gphhalf(mix-1,mjx-1,mkx-1))
      allocate(qvapor(mix-1,mjx-1,mkx-1))
      allocate(qrain(mix-1,mjx-1,mkx-1))
      allocate(qsnow(mix-1,mjx-1,mkx-1))
      allocate(qgraup(mix-1,mjx-1,mkx-1))
      allocate(theta(mix-1,mjx-1,mkx-1))
      allocate(p(mix-1,mjx-1,mkx-1))
      allocate(pb(mix-1,mjx-1,mkx-1))
!      allocate(p_hyd(mix-1,mjx-1,mkx-1))
      allocate(rhodry(mix-1,mjx-1,mkx-1))
      allocate(totpres(mix-1,mjx-1,mkx-1))
      allocate(temp(mix-1,mjx-1,mkx-1))
! 2D
      allocate(mu(mix-1,mjx-1))
      allocate(mub(mix-1,mjx-1))
      allocate(surfptotal(mix-1,mjx-1))
      allocate(q2(mix-1,mjx-1))
      allocate(t2(mix-1,mjx-1))
      allocate(th2(mix-1,mjx-1))
      allocate(psfc(mix-1,mjx-1))
      allocate(slp(mix-1,mjx-1))
      allocate(hgt(mix-1,mjx-1))
      allocate(u10(mix-1,mjx-1))
      allocate(v10(mix-1,mjx-1))
      allocate(windspeed_100m(mix-1,mjx-1))
      allocate(windspeed_200m(mix-1,mjx-1))
      allocate(rainc(mix-1,mjx-1))
!      allocate(rainsh(mix-1,mjx-1))
      allocate(rainnc(mix-1,mjx-1))
      endif
      if (func==1) then
        allocate(initial2d(vars,times,mix-1,mjx-1,ensnum))
        allocate(initmean2d(vars,times,mix-1,mjx-1))
        allocate(variance2d(vars,times,mix-1,mjx-1))
        allocate(covar2d(vars,times,mix-1,mjx-1))
        allocate(initial3d(vars,times,mix-1,mjx-1,levels,ensnum))
        allocate(initmean3d(vars,times,mix-1,mjx-1,levels))
        allocate(variance3d(vars,times,mix-1,mjx-1,levels))
        allocate(covar3d(vars,times,mix-1,mjx-1,levels))
       else if (func==2) then
        allocate(response(times,mix-1,mjx-1,ensnum))
        allocate(rp(ensnum))
       else if (func==3) then
        allocate(corr2d(vars,times,mix-1,mjx-1))
        allocate(sens2d(vars,times,mix-1,mjx-1))
        allocate(stdsens2d(vars,times,mix-1,mjx-1))
        allocate(tstat2d(vars,times,mix-1,mjx-1))
        allocate(corr3d(vars,times,mix-1,mjx-1,levels))
        allocate(sens3d(vars,times,mix-1,mjx-1,levels))
        allocate(stdsens3d(vars,times,mix-1,mjx-1,levels))
        allocate(tstat3d(vars,times,mix-1,mjx-1,levels))
      endif
      END SUBROUTINE ALLOCATE_ARRAYS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE DEALLOCATE_ARRAYS(func)
       ! Subroutine intended to deallocate arrays.
       ! This is particularly necessary when doing cross-grid
       ! sensitivity where the response and IC domains have
       ! different dimensions.
      use wrf_dims_vars
      implicit none
      integer :: func ! 1 is IC, 2 is response, 3 is for at end of program
      if (func==1 .or. func==2) then
! 1D
      deallocate(znw)
      deallocate(znu)
! 3D
      deallocate(u)
      deallocate(v)
      deallocate(w)
      deallocate(ph)
      deallocate(phb)
      deallocate(gpttot)
      deallocate(gph)
      deallocate(gphhalf)
      deallocate(qvapor)
      deallocate(qrain)
      deallocate(qsnow)
      deallocate(qgraup)
      deallocate(theta)
      deallocate(p)
      deallocate(pb)
!      deallocate(p_hyd)
      deallocate(rhodry)
      deallocate(totpres)
      deallocate(temp)
! 2D
      deallocate(mu)
      deallocate(mub)
      deallocate(surfptotal)
      deallocate(q2)
      deallocate(t2)
      deallocate(th2)
      deallocate(psfc)
      deallocate(slp)
      deallocate(hgt)
      deallocate(u10)
      deallocate(v10)
      deallocate(windspeed_100m)
      deallocate(windspeed_200m)
      deallocate(rainc)
!      deallocate(rainsh)
      deallocate(rainnc)
      endif
      if (func==3) then
        deallocate(response)
        deallocate(initial2d)
        deallocate(initial3d)
        deallocate(rp)
        deallocate(initmean2d)
        deallocate(variance2d)
        deallocate(covar2d)
        deallocate(sens2d)
        deallocate(corr2d)
        deallocate(stdsens2d)
        deallocate(tstat2d)
        deallocate(initmean3d)
        deallocate(variance3d)
        deallocate(covar3d)
        deallocate(sens3d)
        deallocate(corr3d)
        deallocate(stdsens3d)
        deallocate(tstat3d)
      endif
      END SUBROUTINE DEALLOCATE_ARRAYS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GET_VARIABLES(filename,permission,unitnum,time,
     & filetype)
       ! Subroutine intended to extract specific variables from the wrf output
       !  Inputs are the filename, permission, unit #, time, and member
      use wrf_dims_vars
      use NETCDF
      implicit none 
      integer :: unitnum,time,permission
      character(len=67) :: filename
      character(len=10) :: filetype
      call open_file(filename,permission,unitnum)
      if (filetype == 'wrfout') then
        call get_variable3d(unitnum,'U',mix-1,mjx-1,mkx-1,time,u)
        call get_variable3d(unitnum,'V',mix-1,mjx-1,mkx-1,time,
     & v)!(:,:,:,time))
        call get_variable3d(unitnum,'W',mix-1,mjx-1,mkx-1,time,
     & w)!(:,:,:,time))
        call get_variable3d(unitnum,'PH',mix-1,mjx-1,mkx,time,
     & ph)!(:,:,:,time))
        call get_variable3d(unitnum,'PHB',mix-1,mjx-1,mkx,time,
     & phb)!(:,:,:,time))
        call get_variable3d(unitnum,'P',mix-1,mjx-1,mkx-1,time,
     & p)!(:,:,:,time))
        call get_variable3d(unitnum,'PB',mix-1,mjx-1,mkx-1,time,
     & pb)!(:,:,:,time))
        call get_variable3d(unitnum,'QVAPOR',mix-1,mjx-1,mkx-1,time,
     & qvapor)!(:,:,:,time))
        call get_variable3d(unitnum,'QRAIN',mix-1,mjx-1,mkx-1,time,
     & qrain)!(:,:,:,time))
        call get_variable3d(unitnum,'QSNOW',mix-1,mjx-1,mkx-1,time,
     & qsnow)!(:,:,:,time))
        call get_variable3d(unitnum,'QGRAUP',mix-1,mjx-1,mkx-1,time,
     & qgraup)!(:,:,:,time))
        call get_variable3d(unitnum,'T',mix-1,mjx-1,mkx-1,time,
     & theta)!(:,:,:,time))
        call get_variable2d(unitnum,'MU',mix-1,mjx-1,time,
     & mu)!(:,:,time))
        call get_variable2d(unitnum,'MUB',mix-1,mjx-1,time,
     & mub)!(:,:,time))
        call get_variable2d(unitnum,'T2',mix-1,mjx-1,time,
     & t2)!(:,:,time))
        call get_variable2d(unitnum,'Q2',mix-1,mjx-1,time,
     & q2)!(:,:,time))
        call get_variable2d(unitnum,'HGT',mix-1,mjx-1,time,
     & hgt)!(:,:,time))
        call get_variable2d(unitnum,'U10',mix-1,mjx-1,time,
     & u10)!(:,:,time))
        call get_variable2d(unitnum,'V10',mix-1,mjx-1,time,
     & v10)!(:,:,time))
        call get_variable2d(unitnum,'TH2',mix-1,mjx-1,time,
     & th2)!(:,:,time))
        call get_variable2d(unitnum,'RAINC',mix-1,mjx-1,time,
     & rainc)!(:,:,time))
        call get_variable2d(unitnum,'RAINNC',mix-1,mjx-1,time,
     & rainnc)!(:,:,time))
        call get_variable1d(unitnum,'ZNW',mkx,1,znw)!(:,time))
        call get_variable1d(unitnum,'ZNU',mkx-1,1,znu)!(:,time))
      else if (filetype == 'TTnc') then
        call get_variable2d(unitnum,'windspeed_100m',mix-1,mjx-1,time,
     & windspeed_100m)!(:,:,time))
        call get_variable2d(unitnum,'windspeed_200m',mix-1,mjx-1,time,
     & windspeed_200m)!(:,:,time))
      endif
      call close_file(unitnum)
      END SUBROUTINE GET_VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE INTERPOLATE_TO_LEVEL(varble,level,vertunit,
     &   outarray)
       ! Subroutine intended to take a 3D array and vertically interpolate to
       !  a specified pressure or height level
       ! var--variable to interpolate
       ! level--vertical level
       ! vertunit--string with vertical unit of level ('P', 'H', or 'M')
       ! outarray--output 2D interpolated array
      use wrf_tools
      use wrf_dims_vars
      implicit none
      character(len=1) :: vertunit
      real :: level
      real :: varble(mix-1,mjx-1,mkx-1)
      real :: outarray(mix-1,mjx-1)
      integer :: i,j,k
      if (vertunit=='P') then
        do i=1,mix-1
          do j=1,mjx-1
            outarray(i,j)=interp_pres(varble(i,j,:),
     &        totpres(i,j,:),level*100.,mkx-1)
          enddo
        enddo
       else if (vertunit=='H') then
        do i=1,mix-1
          do j=1,mjx-1
            outarray(i,j)=interp_hght(varble(i,j,:),gphhalf(i,j,:),
     &                       level,mkx-1)
          enddo
        enddo
       else if (vertunit=='M') then
        outarray=varble(:,:,int(level))
      endif
      end subroutine interpolate_to_level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION RPFUNCTION(resp,num_rptimes,rtype)
       ! Function is intented to find the maximum wind speed within
       ! a time window and return either the peak wind speed or the
       ! time of the peak wind speed.

       ! mix,mjx,mkx--the dimensions of the array
      use wrf_dims_vars
      implicit none
      integer :: num_rptimes, i, max_time
      real :: resp(num_rptimes)
      real :: rptemp,temp_slope,max_slope
      real :: rpfunction
      character(len=10) :: rtype,var 
      rptemp=0.0
      temp_slope=0.0
      max_slope=0.0
      if (rtype=='PEAK_SPD') rptemp=maxval(resp)
      if (rtype=='PEAK_TM') rptemp=maxloc(resp,DIM =1)
      if (rtype=='MIN_SPD') rptemp=minval(resp)
      if (rtype=='MIN_TM') rptemp=minloc(resp,DIM =1)
      if (rtype=='RAMP_SLOPE') then
        do i=1,num_rptimes-1
          temp_slope= resp(i+1)-resp(i)
          if (temp_slope > max_slope) then 
               max_slope = temp_slope
          endif
        enddo
        rptemp = max_slope
      endif       
      if (rtype=='RAMP_TIME') then
        do i=1,num_rptimes-1
          temp_slope= resp(i+1)-resp(i)
          if (temp_slope > max_slope) then
               max_slope = temp_slope
               max_time = i
          endif
        enddo
        rptemp = real(max_time)
      endif
      if (rtype=='DOWN_SLOPE') then
        do i=1,num_rptimes-1
          temp_slope= resp(i)-resp(i+1)
          if (temp_slope > max_slope) then
               max_slope = temp_slope
          endif
        enddo
        rptemp = max_slope
      endif
      if (rtype=='DOWN_TIME') then
        do i=1,num_rptimes-1
          temp_slope= resp(i)-resp(i+1)
          if (temp_slope > max_slope) then
               max_slope = temp_slope
               max_time = i
          endif
        enddo
        rptemp = real(max_time)
      endif
      rpfunction=rptemp
      END FUNCTION RPFUNCTION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SET_FUNCTION(varname,var,level1,level2,vertunit)
      use wrf_dims_vars
      use conversions_mod
      implicit none
      character(len=10) :: varname
      character(len=1) :: vertunit
      real :: level1,level2
      real :: var(mix-1,mjx-1),vartemp1(mix-1,mjx-1),
     &  vartemp2(mix-1,mjx-1),vartemp3d(mix-1,mjx-1,mkx-1)
      integer :: i,j,k,t
      if (varname=='U') call interpolate_to_level(u,
     &  level1,vertunit,var)

      if (varname=='V') call interpolate_to_level(v,
     &  level1,vertunit,var)
          
      if (varname=='W') call interpolate_to_level(w,
     &  level1,vertunit,var)

      if (varname=='MAXW') var=maxval(w,3)

      if (varname=='WIND') then
        call interpolate_to_level(u,level1,vertunit,
     &      vartemp1)
        call interpolate_to_level(v,level1,vertunit,
     &      vartemp2)
        var=sqrt(vartemp1**2+vartemp2**2)
      endif
      
      if (varname=='GPT') call interpolate_to_level(gpttot,
     &  level1,vertunit,var)
      if (varname=='GPH') call interpolate_to_level(gphhalf,
     &  level1,vertunit,var)

      if (varname=='QVAPOR') call interpolate_to_level(qvapor
     &  ,level1,vertunit,var)
      
      if (varname=='THETA') call interpolate_to_level(theta,
     &  level1,vertunit,var)
      if (varname=='T') call interpolate_to_level(temp,
     &level1,vertunit,var)
      if (varname=='TD') then
        ! Dewpoint temperature in Kelvin
        vartemp3d=243.5/(17.67/log((qvapor*
     &      totpres/100./(qvapor+0.622))/6.112)-1.)+273.15
        call interpolate_to_level(vartemp3d,level1,vertunit,var)
      endif

      if (varname=='SHEAR') then
        !Do u
        call interpolate_to_level(u,level1,vertunit,
     &      vartemp1)
        call interpolate_to_level(u,level2,vertunit,
     &      vartemp2)
        var=(vartemp2-vartemp1)**2
        !Do v
        call interpolate_to_level(v,level1,vertunit,
     &      vartemp1)
        call interpolate_to_level(v,level2,vertunit,
     &      vartemp2)
        var=sqrt(var+(vartemp2-vartemp1)**2)
      endif

      if (varname=='T2') var=t2
      if (varname=='TH2') var=th2
      if (varname=='Q2') var=q2
      if (varname=='TD2') then
        ! Dewpoint temperature in Kelvin
        vartemp1=q2*psfc/100./(q2+0.622) ! vapor pressure in mb
        var=243.5/(17.67/log(vartemp1/6.112)-1.)+273.15
      endif
      if (varname=='SLP') var=slp
      if (varname=='U10') var=u10
      if (varname=='V10') var=v10
      if (varname=='WIND100') var=windspeed_100m
      if (varname=='WIND200') var=windspeed_200m
      if (varname=='WIND10') var=sqrt(u10**2+v10**2)
      if (varname=='RAINC') var=rainc
      if (varname=='RAINNC') var=rainnc
      if (varname=='RAINTOT') var=rainc+rainnc
      if (varname=='DBZ') then
        call dbz_calc(qvapor,qrain,
     &      qsnow,
     &      qgraup,temp,totpres,
     &      vartemp3d,mix,mjx,mkx,1,1,1,1)
        call interpolate_to_level(vartemp3d,level1,vertunit,var)
      endif
      if (varname=='MDBZ') then
        call dbz_calc(qvapor,qrain,
     &      qsnow,
     &      qgraup,temp,totpres,
     &      vartemp3d,mix,mjx,mkx,1,1,1,1)
        do j=1,mjx-1
          do i=1,mix-1
            var(i,j)=maxval(vartemp3d(i,j,:))
          enddo
        enddo
      endif
      if (varname=='Z') then
        call dbz_calc(qvapor,qrain,
     &      qsnow,
     &      qgraup,temp,totpres,
     &      vartemp3d,mix,mjx,mkx,1,1,1,1)
        do j=1,mjx-1
          do i=1,mix-1
            var(i,j)=maxval(vartemp3d(i,j,:))
          enddo
        enddo
        var=10.**(var/10.)
      endif

!      print*,var(1,1)
      END SUBROUTINE SET_FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
