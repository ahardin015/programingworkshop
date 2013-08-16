!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   module netcdf - module that contains several subroutines used to 
!                   read and write netcdf data
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!  To compile module:
!  ifort -tpp7 -cm -axW -align -u -convert big_endian -xW -Vaxlib -c
!  module_netcdf.f -o module_netcdf.o
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      MODULE NETCDF

        include 'netcdf.inc'

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   close_file - subroutine that closes a netcdf file
!
!      fid - integer file id
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine close_file(fid)

      integer, intent(in) :: fid

      integer :: rcode

      rcode = nf_close(fid)
      if ( rcode .ne. 0 ) then
         call netcdf_error('File close', 'close', rcode, 0, 0, 0, 0)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_bdy_lev - subroutine that reads one vertical level of WRF 
!                 boundary data
!
!      fid - integer file id
!      var - character name of variable
!      ixy - number of grid points in x or y
!       iz - number of grid points in z
!      wid - number of grid points in width
!       it - time to get
!      dat - data read from file
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_bdy_lev(fid, var, ixy, iz, wid, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ixy, iz, wid, it 
      real, intent(out)             :: dat(ixy,wid)

      integer :: istart(4), iend(4), varid, rcode

      ! get variable id, set dimensions and write data
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ixy
      istart(2) = iz
      iend(2) = 1
      istart(3) = 1
      iend(3) = wid
      istart(4) = it
      iend(4) = 1

      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode .ne. 0) then
        call netcdf_error('get_bdy_lev', var, rcode, ixy, iz, wid, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_dimlen - function that reading the length of a dimension
!
!       fid - integer file id
!   dimname - name of dimension
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function get_dimlen(fid, dimname)

      character (len=*), intent(in) :: dimname
      integer, intent(in)           :: fid

      integer :: did, rcode, get_dimlen

      rcode = nf_inq_dimid(fid, dimname, did)
      rcode = nf_inq_dimlen(fid, did, get_dimlen)

      return
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_varble_attr - function that gets the value of an attribute 
!                     of a variable within a netcdf file. (real)
!
!      fid - integer file id
!     vnam - name of variable
!     anam - name of attribute
!
!     created June 2005 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function get_varble_attr(fid, vnam, anam)

      character (len=*), intent(in) :: vnam, anam
      integer, intent(in)           :: fid

      integer :: rcode, varid
      real :: get_varble_attr

      rcode = nf_inq_varid(fid, vnam, varid)
      rcode = nf_get_att_real(fid, varid, anam, get_varble_attr)

      return
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable1d - subroutine that reads netcdf data for a 1    
!                    dimensional quantity
!
!     fid - integer file id
!     var - name of the variable to read
!      iz - z dimension size
!      it - time to read the data from
!     dat - array containing data that has been read
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable1d(fid, var, iz, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, iz, it
      real, intent(out)             :: dat(iz)

      integer :: istart(2), iend(2), varid, rcode
      
      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = iz
      istart(2) = it
      iend(2) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable1d', var, rcode, iz, 0, 0, it)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable1d_local - subroutine that reads netcdf data for a 1    
!                          dimensional quantity over an interval
!
!       fid - integer file id
!       var - name of the variable to read
!       izs - z dimension start
!       ize - z dimension end
!        it - time to read the data from
!       dat - array containing data that has been read
!
!     created June 2004 Sebastien Dirren, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable1d_local(fid, var, izs, ize, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, izs, ize, it
      real, intent(out)             :: dat(ize-izs+1)

      integer :: istart(2), iend(2), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = izs
      iend(1) = ize-izs+1
      istart(2) = it
      iend(2) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        print*,'izs', istart(1)
        call netcdf_error('get_variable1d_local', var, rcode, 
     &                         iend(1), 0, 0, it) 
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable1d_int - subroutine that reads netcdf integer data 
!                        for a 1 dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to read
!       iz - z dimension size
!       it - time to read the data from
!      dat - array containing data that has been read
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable1d_int(fid, var, iz, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, iz, it
      integer, intent(out)          :: dat(iz)

      integer :: istart(2), iend(2), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = iz
      istart(2) = it
      iend(2) = 1

      ! read the dat
      rcode = nf_get_vara_int(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable1d_int', var, rcode, iz, 0, 0,it)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable2d - subroutine that reads netcdf data for a 2    
!                    dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to read
!       i1 - x dimension size
!       i2 - y dimension size
!       it - time to read the data from
!      dat - array containing data that has been read
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable2d(fid, var, i1, i2, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, i1, i2, it
      real, intent(out)             :: dat(i1,i2)

      integer :: istart(3), iend(3), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = i1
      istart(2) = 1
      iend(2) = i2
      istart(3) = it
      iend(3) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable2d', var, rcode, i1, i2, 0, it)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable2d_local - subroutine that reads netcdf data for a 2    
!                          dimensional quantity over a patch
!
!     fid - integer file id
!     var - name of the variable to read
!     i1s - x dimension start
!     i1e - x dimension end
!     i2s - x dimension start
!     i2e - x dimension end
!      it - time to read the data from
!     dat - array containing data that has been read
!
!     created June 2004 Sebastien Dirren, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable2d_local(fid, var, 
     &                           i1s, i1e, i2s, i2e, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, i1s, i1e, i2s, i2e, it
      real, intent(out)             :: dat(i1e-i1s+1,i2e-i2s+1)

      integer :: istart(3), iend(3), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = i1s
      iend(1) = i1e-i1s+1
      istart(2) = i2s
      iend(2) = i2e-i2s+1
      istart(3) = it
      iend(3) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        print*,'i1s i2s', istart(1:2)
        call netcdf_error('get_variable2d_local', var, rcode, 
     &                             iend(1), iend(2), 0, it)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable3d - subroutine that reads netcdf data for a 3 
!                    dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to read
!       i1 - x dimension size
!       i2 - y dimension size
!       i3 - z dimension size
!       it - time to read the data from
!      dat - array containing data that has been read
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable3d(fid, var, i1, i2, i3, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, i1, i2, i3, it
      real, intent(out)             :: dat(i1,i2,i3)

      integer :: istart(4), iend(4), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = i1
      istart(2) = 1
      iend(2) = i2
      istart(3) = 1
      iend(3) = i3
      istart(4) = it
      iend(4) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable3d', var, rcode, i1, i2, i3, it)
      endif
 
      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable3d_local - subroutine that reads netcdf data for a
!                          patch of a 3 dimensional quantity
!
!       fid - integer file id
!       var - name of the variable to read
!       i1s - x dimension start
!       i1e - x dimension end
!       i2s - y dimension start
!       i2e - y dimension end
!       i3s - x dimension start
!       i3e - x dimension end
!        it - time to read the data from
!       dat - array containing data that has been read
!
!     created June 2004 Sebastien Dirren, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable3d_local(fid, var, 
     &                          i1s, i1e, i2s, i2e, i3s, i3e, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid,i1s,i1e,i2s,i2e,i3s,i3e,it
      real, intent(out)           :: dat(i1e-i1s+1,i2e-i2s+1,i3e-i3s+1)

      integer :: istart(4), iend(4), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = i1s
      iend(1) = i1e-i1s+1
      istart(2) = i2s
      iend(2) = i2e-i2s+1
      istart(3) = i3s
      iend(3) = i3e-i3s+1
      istart(4) = it
      iend(4) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        print*,'i1s i2s i3s', istart(1:3)
        call netcdf_error('get_variable3d_local', var, rcode, 
     &                       iend(1), iend(2), iend(3), it)
      endif
 
      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable_lev - subroutine that reads netcdf data for one vertical 
!                      level of a 3-D quantity
!
!      fid - integer file id
!      var - name of the variable to read
!       ix - x dimension size
!       iy - y dimension size
!       iz - z level
!       it - time to read the data from
!      dat - array containing data that has been read
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable_lev(fid, var, ix, iy, iz, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ix, iy, iz, it
      real, intent(out)             :: dat(ix,iy)

      integer :: istart(4), iend(4), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = ix
      istart(2) = 1
      iend(2) = iy
      istart(3) = iz
      iend(3) = 1
      istart(4) = it
      iend(4) = 1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable_lev', var, rcode, ix, iy, iz,it)
      endif
 
      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable_vec - subroutine that reads netcdf data for a vector    
!                      of data
!
!       fid - integer file id
!       var - name of the variable to read
!        iz - length of the vector of data
!       dat - array containing data that has been read
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable_vec(fid, var, iz, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, iz
      real, intent(out)             :: dat(iz)

      integer :: istart, iend, varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart = 1
      iend = iz

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable_vec', var, rcode, iz, 0, 0, 0)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable_veci - subroutine that reads netcdf data for a vector    
!                       of integer data
!
!       fid - integer file id
!       var - name of the variable to read
!        iz - length of the vector of data
!       dat - array containing data that has been read
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable_veci(fid, var, iz, dat)

      character (len=*), intent(in)  :: var
      integer, intent(in)            :: fid, iz
      integer, intent(out)           :: dat(iz)

      integer :: istart, iend, varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart = 1
      iend = iz

      ! read the data
      rcode = nf_get_vara_int(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable_vec', var, rcode, iz, 0, 0, 0)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_variable_vec_local - subroutine that reads netcdf data for a     
!                            vector of data over a limited area
!
!       fid - integer file id
!       var - name of the variable to read
!    istart - index of beginning of data to read
!        ie - index of end of data to read
!       dat - array containing data that has been read
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_variable_vec_local(fid, var, istart, ie, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, istart, ie
      real, intent(out)             :: dat(ie-istart+1)

      integer :: iend, varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      iend = ie-istart+1

      ! read the data
      rcode = nf_get_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('get_variable_vec_local',var,rcode,istart,
     &                            ie,0,0)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_text - subroutine that reads netcdf text data for 1 dimensional 
!              quantity 
!
!      fid - integer file id
!      var - name of the variable to read
!       i1 - length of string to read
!       it - time to read the data from
!  outtext - array containing data that has been read
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_text(fid, var, i1, it, outtxt)

      character (len=*), intent(in)   :: var
      integer, intent(in)             :: fid, i1, it
      character (len=i1), intent(out) :: outtxt

      integer :: istart(2), iend(2), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = i1
      istart(2) = it
      iend(2) = 1

      ! read the dat
      rcode = nf_get_vara_text(fid, varid, istart, iend, outtxt)

      if (rcode.ne.0) then 
        call netcdf_error('get_text', var, rcode, i1, 0, 0, it)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   netcdf_error - subroutine that prints an error and stops a program 
!                  if a netcdf routine throws an error
!
!     place - program where the netcdf error happens
!       var - variable of error
!        rc - integer code of error
!        ix - dimension 1 of error
!        iy - dimension 2 of error
!        iz - dimension 3 of error
!        it - time dimension of error
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine netcdf_error(place, var, rc, ix, iy, iz, it)

      character (len=*), intent(in) :: place, var
      integer, intent(in)           :: rc, ix, iy, iz, it

      write(6,*) 'NETCDF ERROR !!!!!!!!'
      write(6,*)            
      write(6,*) 'Error occurs at ',place
      write(6,*) 'With variable ',var
      write(6,*) 'Netcdf message ',nf_strerror(rc)
      write(6,*) 'ix iy iz it',ix,iy,iz,it

      open(25, file='error_netcdf',status='unknown')
      write(25,*) place; close(25)

      stop
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   open_file - subroutine that opens a netcdf file 
!
!   filename - name of file to open
!    permiss - permissions of file to open
!        fid - integer file id
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine open_file(filename, permiss, fid)

      character (len=*), intent(in) :: filename
      integer, intent(in)           :: permiss
      integer, intent(out)          :: fid

      integer :: rcode

      rcode = nf_open(filename, permiss, fid)
      if ( rcode .ne. 0 ) then
         call netcdf_error(filename, 'open', rcode, 0, 0, 0, 0)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_bdy_lev - subroutine that writes one vertical level of 
!                   boundary data
!
!     fid - integer file id
!     var - name of variable
!     ixy - number of grid points in x or y
!      iz - vertical level to write out
!     wid - width of boundary data
!      it - time to write
!     dat - data to write out
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_bdy_lev(fid, var, ixy, iz, wid, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ixy, iz, wid, it
      real, intent(in)              :: dat(ixy,wid)

      integer :: istart(4), iend(4), varid, rcode

      ! get variable id, set dimensions and write data
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ixy
      istart(2) = iz
      iend(2) = 1
      istart(3) = 1
      iend(3) = wid
      istart(4) = it
      iend(4) = 1

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)

      if (rcode .ne. 0) then
        call netcdf_error('write_bdy_lev', var, rcode, ixy, iz, wid, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_text - subroutine that writes netcdf text data for 1 
!                dimensional quantity 
!
!      fid - integer file id
!      var - name of the variable to write
!       i1 - length of string to write
!       it - time to write the data from
!  outtext - array containing data to be written
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_text(fid, var, i1, it, outtxt)

      character (len=*), intent(in)   :: var
      integer, intent(in)             :: fid, i1, it
      character (len=i1), intent(in ) :: outtxt

      integer :: istart(2), iend(2), varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      istart(1) = 1
      iend(1) = i1
      istart(2) = it
      iend(2) = 1

      ! read the dat
      rcode = nf_put_vara_text(fid, varid, istart, iend, outtxt)

      if (rcode.ne.0) then 
        call netcdf_error('write_text', var, rcode, i1, 0, 0, it)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable1d - subroutine that writes netcdf dat for a 1
!                          dimensional quantity
!
!       fid - integer file id
!       var - name of the variable to write
!        ix - x dimension size
!        it - time to write out data
!       dat - array containing dat to write
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable1d(fid, var, ix, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ix, it
      real, intent(in)              :: dat(ix)

      integer :: istart(2), iend(2), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ix
      istart(2) = it
      iend(2) = 1

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable1d', var, rcode, ix, 0, 0, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable1d_int - subroutine that writes netcdf dat for a 1
!                          dimensional integer quantity
!
!     fid - integer file id
!     var - name of the variable to write
!      ix - x dimension size
!      it - time to write out data
!     dat - array containing dat to write
!
!     created June 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable1d_int(fid, var, ix, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ix, it, dat(ix)
      
      integer :: istart(2), iend(2), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ix
      istart(2) = it
      iend(2) = 1

      rcode = nf_put_vara_int(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable1d_int', var, rcode, ix, 0, 
     &                         0, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable2d - subroutine that writes netcdf dat for a 2
!                      dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to write
!       i1 - x dimension size
!       i2 - y dimension size
!       it - time to write out data
!      dat - array containing dat to write
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable2d(fid, var, i1, i2, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, i1, i2, it
      real, intent(in)              :: dat(i1,i2)

      integer :: istart(3), iend(3), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)
    
      print*, 'in write_variable2d module'

      istart(1) = 1
      iend(1) = i1
      istart(2) = 1
      iend(2) = i2
      istart(3) = it
      iend(3) = 1
      print*, 'checkpoint 1'

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      print*, 'checkpoint2'
      print*, rcode
      if (rcode.ne.0) then 
        call netcdf_error('write_variable2d', var, rcode, i1, i2, 0, it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable3d - subroutine that writes netcdf data for a 3 
!                      dimensional quantity
!
!      fid - integer file id
!      var - name of the variable to write
!       i1 - x dimension size
!       i2 - y dimension size
!       i3 - z dimension size
!       it - time to write data
!      dat - array containing data to write
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable3d(fid, var, i1, i2, i3, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, i1, i2, i3, it
      real, intent(in)              :: dat(i1,i2,i3)

      integer :: istart(4), iend(4), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = i1
      istart(2) = 1
      iend(2) = i2
      istart(3) = 1
      iend(3) = i3
      istart(4) = it
      iend(4) = 1

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable3d', var, rcode, i1, i2, i3,it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_var_lev - subroutine that writes netcdf data for one vertical 
!                   level of a 3D field
!
!     fid - integer file id
!     var - name of the variable to write
!      ix - x dimension size
!      iy - y dimension size
!      iz - vertical level to write
!      it - time to write out data
!     dat - array containing data to write
!
!     created June 2003 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_var_lev(fid, var, ix, iy, iz, it, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, ix, iy, iz, it
      real, intent(in)              :: dat(ix,iy)

      integer :: istart(4), iend(4), varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart(1) = 1
      iend(1) = ix
      istart(2) = 1
      iend(2) = iy
      istart(3) = iz
      iend(3) = 1
      istart(4) = it
      iend(4) = 1

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_var_lev', var, rcode, ix, iy, iz,it)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable_vec -  subroutine that writes netcdf data for a 
!                         vector of data
!
!      fid - integer file id
!      var - name of the variable to write
!       iz - x dimension size
!      dat - array containing dat to write
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable_vec(fid, var, iz, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, iz
      real, intent(in)              :: dat(iz)

      integer :: istart, iend, varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart = 1
      iend   = iz

      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable_vec', var, rcode, iz, 0, 0, 0)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable_veci -  subroutine that writes netcdf data for a 
!                          vector of integer data
!
!      fid - integer file id
!      var - name of the variable to write
!       iz - x dimension size
!      dat - array containing dat to write
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable_veci(fid, var, iz, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, iz, dat(iz)

      integer :: istart, iend, varid, rcode

      ! get variable id, set dimensions and write dat
      rcode = nf_inq_varid(fid, var, varid)

      istart = 1
      iend   = iz

      rcode = nf_put_vara_int(fid, varid, istart, iend, dat)
      
      if (rcode.ne.0) then 
        call netcdf_error('write_variable_vec', var, rcode, iz, 0, 0, 0)
      endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable_vec_local - subroutine that writes netcdf data for a     
!                              vector of data over a limited area
!
!       fid - integer file id
!       var - name of the variable to read
!    istart - index of beginning of data to read
!        ie - index of end of data to read
!       dat - array containing data that has been read
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable_vec_local(fid, var, istart, ie, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, istart, ie
      real, intent(in)              :: dat(ie-istart+1)

      integer :: iend, varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      iend = ie-istart+1

      ! read the data
      rcode = nf_put_vara_real(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('write_variable_vec_local',var,rcode,istart,
     &                            ie,0,0)
      endif

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write_variable_veci_local - subroutine that writes integer netcdf 
!                               data for a vector of data over a 
!                               limited area
!
!       fid - integer file id
!       var - name of the variable to read
!    istart - index of beginning of data to read
!        ie - index of end of data to read
!       dat - integer array containing data that has been read
!
!     created Oct. 2004 Ryan Torn, U. Washington
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_variable_veci_local(fid, var, istart, ie, dat)

      character (len=*), intent(in) :: var
      integer, intent(in)           :: fid, istart, ie, dat(ie-istart+1)

      integer :: iend, varid, rcode

      ! get the variable id
      rcode = nf_inq_varid(fid, var, varid)

      ! set the dimension arrays for reading
      iend = ie-istart+1

      ! read the data
      rcode = nf_put_vara_int(fid, varid, istart, iend, dat)

      if (rcode.ne.0) then 
        call netcdf_error('write_variable_veci_local',var,rcode,istart,
     &                            ie,0,0)
      endif

      return
      end subroutine

      END MODULE NETCDF
