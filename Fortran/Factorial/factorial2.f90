! Basic Factorial Calculation using Functions FORTRAN programming example
!
! Written for inclusion as part of the Texas Tech University 
! Atmospheric Sciences Basic Programming Workshop
!
! Copyright 2013, Timothy Sliwinski
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

PROGRAM Factorial

IMPLICIT NONE

INTEGER::num,result
INTEGER::Fact           ! Must declare function variable

num = 4

result = Fact(num) ! Function calls require "=" operator

PRINT*, result

END PROGRAM Factorial

!===========================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION: Fact
! USE: Computes and returns the factorial of the input number
! IN: Integer NUM (Mandatory)
! OUT: Return value (Fact) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION Fact(NUM)

INTEGER::NUM
INTEGER::i

Fact   = NUM
i      = NUM-1

DO WHILE (i>1)
	Fact=Fact*i
	i=i-1
END DO

RETURN  ! Required at the end of every function, otherwise no value is returned.
	! Value returned is the value in the variable with the same name as the
	! function (e.g in this function, that is the "Fact" variable)

END FUNCTION Fact






