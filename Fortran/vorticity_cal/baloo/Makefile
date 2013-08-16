# Makefile to compile sensitivity modules and main program
#   Created by Chris Bednarczyk, Texas Tech University, 2013

F90 = ifort
#F90 = /lustre/work/apps/gcc-4.4.7/bin/gfortran

#FFLAGS = -O2
FFLAGS = -O0

SRC = .
BIN = .
LIBCDF = $(NETCDF)/lib
INCCDF = $(NETCDF)/include
NETCDFCOMP = -L$(LIBCDF) -lnetcdf -I$(INCCDF)
OBJMOD = $(SRC)/module_netcdf.o $(SRC)/module_wrf_tools.o $(SRC)/module_map_utils.o $(SRC)/conversions_mod.o $(SRC)/wrfdimsvars_mod.o

esens :
	$(F90) -c module_netcdf.f
	$(F90) -c wrfdimsvars_mod.f
	$(F90) -c module_map_utils.f
	$(F90) -c module_wrf_tools.f
	$(F90) -c conversions_mod.f
	$(F90) $(FFLAGS) -o $(BIN)/made_file $(OBJMOD) vort_calc.f -limf $(NETCDFCOMP)
