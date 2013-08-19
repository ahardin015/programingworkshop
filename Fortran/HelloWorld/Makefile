all: HelloWorld

HelloWorld: HelloWorld.o
	gfortran -o HelloWorld HelloWorld.o

HelloWorld.o: HelloWorld.f90
	gfortran -c HelloWorld.f90

clean:
	rm -rf *.o *.exe HelloWorld
