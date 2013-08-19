# Basic Factorial Makefile example
#
# Written for inclusion as part of the Texas Tech University 
# Atmospheric Sciences Basic Programming Workshop
#
# Copyright 2013, Timothy Sliwinski
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

all: factorial factorial2 myfact

factorial: factorial.o
	gfortran factorial.o -o factorial

factorial.o: factorial.f90
	gfortran -c factorial.f90

factorial2: factorial2.o
	gfortran factorial2.o -o factorial2

factorial2.o: factorial2.f90
	gfortran -c factorial2.f90

myfact: fact.o myfactorial.o
	gfortran fact.o myfactorial.o -o myfact

fact.o: fact.f90
	gfortran -c fact.f90

myfactorial.o: myfactorial.f90
	gfortran -c myfactorial.f90

clean: 
	rm -rf *.o *.exe factorial factorial2 myfact
