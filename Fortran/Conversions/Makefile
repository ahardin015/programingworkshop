all: TestConvert

TestConvert: TestConversions.o TemperatureConversions.o
	gfortran -o TestConvert TestConversions.o TemperatureConversions.o

TestConversions.o: TestConversions.f90 TemperatureConversions.o
	gfortran -c TestConversions.f90

TemperatureConversions.o:
	gfortran -c TemperatureConversions.f90

clean:
	rm -rf *.o *.mod *.exe TestConvert
