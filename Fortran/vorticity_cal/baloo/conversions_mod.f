      module conversions_mod
! The wrf_dims_vars is where I will define the arrays
!  for the usual variables
      use wrf_dims_vars
      implicit none
      
      real, parameter :: basethet=300.
      real, parameter :: rdgas=287.
      real, parameter :: rvgas=461.6
      real, parameter :: preffer=100000.
      real, parameter :: kap=2./7.
      real, parameter :: grav=9.81
      real, parameter :: prefcont=preffer**(-(kap))
      real, parameter :: kapdiv=1./(1.-kap)
      contains
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine unstagger_winds()
! Subroutine to unstagger either u, v, and w winds
      use wrf_dims_vars
      implicit none
      integer :: i,j,k
      do k=1,mkx-1
      do j=1,mjx-1
      do i=1,mix-1
        u(i,j,k)=(u(i,j,k)+u(i+1,j,k))/2.
        v(i,j,k)=(v(i,j,k)+v(i,j+1,k))/2.
        w(i,j,k)=(w(i,j,k)+w(i,j,k+1))/2.
      enddo
      enddo
      enddo
      end subroutine unstagger_winds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gpot_and_gph()
! Subroutine to add geopotential perturbation to the base state geopotential
! Also calculates geopotential height from geopotential
      use wrf_dims_vars
      use wrf_tools
      implicit none
      integer i,j,k
      gpttot=ph+phb
! Geopotential height
      gph=gpttot/grav
! Geopotential height on half levels
      do i=1,mix-1
        do j=1,mjx-1
          call destag_zstag(znu,znw,mkx-1,gph(i,j,:),
     &     gphhalf(i,j,:))
        enddo
      enddo
      end subroutine gpot_and_gph
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sfc_pressure()
! Subroutine to calculate total surface pressure, dry density, total pressure
!  at a level, temperature, and psfc enroute to slp.
      use wrf_dims_vars
      use wrf_tools
      implicit none
      integer :: i,j,k
      real :: levee

!      do i=1,mix-1
!        do j=1,mjx-1
!          surfptotal(i,j)=mu(i,j)+mub(i,j)
!          do k=1,mkx-1
!            levee=znw(k+1)-znw(k)
!             rhodry(i,j,k)=-((levee)/(gpttot(i,j,k+1)-gpttot(i,j,k)))*
!     &         surfptotal(i,j)
!             totpres(i,j,k)=(theta(i,j,k)*prefcont*rhodry(i,j,k)*rdgas*
!     &         (1.+(rvgas/rdgas)*qvapor(i,j,k)))**kapdiv
!             temp(i,j,k)=theta(i,j,k)*((preffer/totpres(i,j,k))**-kap)
!      enddo ! end k loop

      surfptotal=mu+mub
      do k=1,mkx-1
        levee=znw(k+1)-znw(k)
        rhodry(:,:,k)=-((levee)/(gpttot(:,:,k+1)-
     &    gpttot(:,:,k)))*surfptotal
      enddo ! end k loop
      totpres=(theta*prefcont*rhodry*rdgas*
     &    (1.+(rvgas/rdgas)*qvapor))**kapdiv
      temp=theta*((preffer/totpres)**(-kap))

      psfc=totpres(:,:,1)*(2.71828**((gphhalf(:,:,1)-hgt)/
     &    (29.3*temp(:,:,1)*(1.+qvapor(:,:,1)))))
      do j=1,mjx-1
      do i=1,mix-1
        slp(i,j)=slp_standard_atmos(t2(i,j),
     &    psfc(i,j),q2(i,j),hgt(i,j))
      enddo
      enddo
!      print*,t2(200,80),psfc(200,80),q2(200,80),hgt(200,80)
      end subroutine sfc_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine variable_processing()
      ! subroutine calls other variable processing subroutines
      !   unstagger winds
      !   convert to full pot. temp
      !   calculate geopotential and gph
      !   calculate sfc pressure
      use wrf_dims_vars
      implicit none
      call unstagger_winds()
      theta=theta+basethet
      call gpot_and_gph()
      call sfc_pressure()
      end subroutine variable_processing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dbz_calc(qvp,qra,qsn,qgr,tmk,pres,dbz,
     &  mix,mjx,mkx,iliqskin,in0r,in0s,in0g)

! THIS SUBROUTINE IS COPIED FROM THE RIP CALCULATION OF SIMULATED
! REFLECTIVITY.  A FEW MODIFICATIONS WERE MADE, BUT NO EQUATIONS WERE
! CHANGED.

c     This routine computes equivalent reflectivity factor (in dBZ) at
c     each model grid point.  In calculating Ze, the RIP algorithm makes
c     assumptions consistent with those made in an early version
c     (ca. 1996) of the bulk mixed-phase microphysical scheme in the MM5
c     model (i.e., the scheme known as "Resiner-2").  For each species:
c
c     1. Particles are assumed to be spheres of constant density.  The
c     densities of rain drops, snow particles, and graupel particles are
c     taken to be rho_r = rho_l = 1000 kg m^-3, rho_s = 100 kg m^-3, and
c     rho_g = 400 kg m^-3, respectively. (l refers to the density of
c     liquid water.)
c
c     2. The size distribution (in terms of the actual diameter of the
c     particles, rather than the melted diameter or the equivalent solid
c     ice sphere diameter) is assumed to follow an exponential
c     distribution of the form N(D) = N_0 * exp( lambda*D ).
c
c     3. If in0X=0, the intercept parameter is assumed constant (as in
c     early Reisner-2), with values of 8x10^6, 2x10^7, and 4x10^6 m^-4,
c     for rain, snow, and graupel, respectively.  Various choices of
c     in0X are available (or can be added).  Currently, in0X=1 gives the
c     variable intercept for each species that is consistent with
c     Thompson, Rasmussen, and Manning (2004, Monthly Weather Review,
c     Vol. 132, No. 2, pp. 519-542.)
c
c     4. If iliqskin=1, frozen particles that are at a temperature above
c     freezing are assumed to scatter as a liquid particle.
c
c     More information on the derivation of simulated reflectivity in RIP
c     can be found in Stoelinga (2005, unpublished write-up).  Contact
c     Mark Stoelinga (stoeling@atmos.washington.edu) for a copy.

      implicit none

      integer :: i,j,k
      integer :: mix,mjx,mkx
      real :: qvp(mix-1,mjx-1,mkx-1),qra(mix-1,mjx-1,mkx-1),
     & qsn(mix-1,mjx-1,mkx-1),
     & qgr(mix-1,mjx-1,mkx-1),tmk(mix-1,mjx-1,mkx-1),
     & pres(mix-1,mjx-1,mkx-1),
     & dbz(mix-1,mjx-1,mkx-1)
      integer :: iliqskin,in0r,in0s,in0g

      real, parameter :: pi=4.*atan(1.),rgas=287.04,eps=0.622,
     & rhowat=1000.,
     & rn0_r=8.e6,
     & rn0_s=2.e7,
     & rn0_g=4.e6,
     & celkel=273.15,
     & gamma_seven=720.,
     & alpha=0.224,
     & rho_r=rhowat,
     & rho_s=100.,
     & rho_g=400.,
     & gon=5.e7,
     & ron_min=8.e6,
     & ron=8.e6,
     & ron2=1.e10,
     & son=2.e7,
     & r1=1.e-15,
     & ron_qr0=1.e-4

      real :: factor_r,factor_s,factorb_s,factor_g,factorb_g,
     & ronv,sonv,gonv,ron_delqr0,
     & ron_const1r,ron_const2r
      real :: rhoair,temp_c,z_e,t_virtual

      factor_r=gamma_seven*1.e18*(1./(pi*rho_r))**1.75
      factor_s=gamma_seven*1.e18*(1./(pi*rho_s))**1.75
     & *(rho_s/rhowat)**2 * alpha
      factor_g=gamma_seven*1.e18*(1./(pi*rho_g))**1.75
     & *(rho_g/rhowat)**2 * alpha
      ron_delqr0=0.25*ron_qr0
      ron_const1r=(ron2-ron_min)*0.5
      ron_const2r=(ron2+ron_min)*0.5

      do k=1,mkx-1
      do j=1,mjx-1
      do i=1,mix-1
!       if (qgr(i,j,k)>0) print*, qgr(i,j,k)
!       if (mod(i,50)==0 .and. mod(j,50)==0) then
!         print*,'Pres',pres(i,j,k),'Temp',tmk(i,j,k)
!         print*,'qvapor',qvp(i,j,k),'qrain',qra(i,j,k)
!         print*,'qsnow',qsn(i,j,k),'qgraup',qgr(i,j,k)
!       endif

        t_virtual=tmk(i,j,k)*(1.+qvp(i,j,k)/eps)/(1.+qvp(i,j,k))
        rhoair=pres(i,j,k)*1./(rgas*t_virtual)

        if (iliqskin==1 .and. tmk(i,j,k)>celkel) then
          factorb_s=factor_s/alpha
          factorb_g=factor_g/alpha
         else
          factorb_s=factor_s
          factorb_g=factor_g
        endif

        if (in0s==1) then
          temp_c=min(-0.001,tmk(i,j,k)-celkel)
          sonv=min(2.0e8,2.0e6*exp(-0.12*temp_c))
         else
          sonv=rn0_s
        endif
        if (in0g==1) then
          gonv=gon
          if (qgr(i,j,k)>r1) then
            gonv=2.38*(pi*rho_g/(rhoair*qgr(i,j,k)))**0.92
            gonv=max(1.e4,min(gonv,gon))
          endif
         else
          gonv=rn0_g
        endif
        if (in0r==1) then
          ronv=ron2
          if (qra(i,j,k)>r1) then
            ronv=ron_const1r*tanh((ron_qr0-qra(i,j,k))
     &       /ron_delqr0)+ron_const2r
          endif
         else
          ronv=rn0_r
        endif
        z_e=factor_r*(rhoair*qra(i,j,k))**1.75/ronv**0.75 +
     &    factorb_s*(rhoair*qsn(i,j,k))**1.75/sonv**0.75 +
     &    factorb_g*(rhoair*qgr(i,j,k))**1.75/gonv**0.75

!        if (mod(i,50)==0 .and. mod(j,50)==0) then
!          print*,z_e
!        endif
        z_e=max(z_e,0.001)
        dbz(i,j,k)=10.*log10(z_e)
!        if (mod(i,50)==0 .and. mod(j,50)==0) then
!          print*,'rhoair',rhoair
!          print*,i,j,k,'dbz',dbz(i,j,k)
!          print*
!        endif

      enddo
      enddo
      enddo

      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module
