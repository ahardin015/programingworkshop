HELLO WORLD!

=========================

This is the first FORTRAN program everyone becomes familiar with, and yet, nearly the most useless for anything other than demonstrating the 
basics of how a FORTRAN program is laid out and formatted.

The first thing you should notice is the PROGRAM declaration line. Every FORTRAN program begins with this declaration as it identifies the 
start of the program and also the name of the program. At the end of the program, you can see the PROGRAM statement again, but this time 
prefaced by the END keyword. Everything between the first PROGRAM declaration line and the END PROGRAM line is the body of the program where 
variables can be declared, operations can be performed, and functions and subroutines can be called. 

The most simple program that will compile and run is a two line program comprised of only the PROGRAM and END PROGRAM declarations. Hello 
World is therefore but a small step above that. It is important to note that FORTRAN uses the PROGRAM and END PROGRAM "block" layout for 
subroutine, function, loop, and module blocks using SUBROUTINE/END SUBROUTINE, FUNCTION/END FUNCTION, DO/END DO, and MODULE/END MODULE, 
respectively, so it is good to be familiar with that setup now. We will see these other blocks in later exmples, but for now, lets return to 
the Hello World example.

After we have declared our program and named it, the next step to getting useful output is putting the line "IMPLICIT NONE". In the early 
days of FORTRAN, certain variables would be assigned a type (INTEGER, REAL, DOUBLE PRECISION, CHARACTER, etc) based on the variable's name. 
However, this often lead to a lot of confusion during the debug process as variable names could not be very descriptive. To preserve 
backwards compatibility with these older programs, implicit typing is still supported by modern FORTRAN compilers, however it is 
considered bad practice to use implicit typing anymore so you should ALWAYS include the line "IMPLICIT NONE" to disable implicit typing of 
variables.

Once IMPLICIT NONE is declared, the next section would normally be the place to declare variables to be used in the program. Variables 
MUST be declared at the start of the program before they can be used. However, for the case of the Hello World program, there are no 
variables to be declared, and therefore, we move right into the execution portion of the program with the PRINT* line. In later programs, we 
will declare variables to use, and will revisit this topic then. 

The only portion of the Hello World that actually executes an operation is the PRINT* line. PRINT* is a basic operation that displays 
unformatted output text to the command line. In this case, we are simply printing the string "Hello World!". The two quotation marks (" ") 
signify that the characters between them are to be displayed as part of a string of characters. The values of variables can also be printed 
to the terminal by simply putting the variable name after the PRINT* without any quotation marks. However, to do so, you should also notice 
the comma after the PRINT* keyword. Commas are used to seperate different types of data to be displayed so that PRINT* can display them 
properly. For instance, if we wanted to print the value of a variable called "num" after printing "Hello World!", then we would change the line from:

	PRINT*, "Hello World!"

to instead be:

	PRINT*, "Hello World!", num

This way, PRINT* will interpret num as a variable who's value we wish to print and not a string of characters we want literally printed to 
the screen.

Compiling and Running Hello World

To compile a FORTRAN program, you first need to check that there is a FORTRAN compiler available on the system you are using. On Texas Tech's 
main Linux-based supercomputer, Hrothgar, you have the choice to use either the intel fortran compiler (ifort) or the open-source GNU fortran 
compiler (gfortran). For the examples here, we will opt to use the gfortran compiler, however what is covered should be easily reproducable 
with the Intel compiler as well. 

To compile our hello world example, we only need to issue one command as follows:

	gfortran HelloWorld.f90 -o HelloWorld

This command takes care of the full process of converting the source file (HelloWorld.f90) into an executable (HelloWorld). If we dissect 
this command into it's parts, there are three main portions. The first is the command itself which simply calls the compiler (gfortran). 
Next, we have to tell the compiler which file to compile, so we put the filename for our HelloWorld program, HelloWorld.f90. Finally, we have 
the option of telling the compiler what we want the result program to be called using the -o (dash lowercase letter O) option. If we do not 
use this -o option, the compiler will automatically name the new program a.out (or, a.exe depending on your system). Since a.out is not a 
very descriptive name and having descriptive names is an important part of having easily understandable code, it is recommended to always use 
the -o option to rename the resulting executable. 

After you run the above command, you should have your first executable. All that is left to do it to run it and make sure it works as you 
expect. To run this new program, simply type:

	./HelloWorld

and you should see "Hello World!" appear in the terminal.

Compiling Using Makefiles

Unix and Linux systems are often distributed with a utility called "make" with reads in a special file called a Makefile with helps to make 
the compilation process simpler. Each of the example programs we'll use will also have a Makefile included to mae compilation easier. To 
compile all of the program(s) in a file at once, you simply need to type "make" at the terminal, and make will read through the dependencies 
in the Makefile to create a properly configured executable program. 

If you look at the Makefile for the HelloWorld program, you'll see another file referenced that is not currently in your directory, but has 
the same HelloWorld name with a ".o" extension. This is an object file that is created when gfortran compiles a source code with the -c 
option. An object file is NOT an executable file, but rather a compiled form of the code in a source file that is ready to be linked with 
other .o files to create a final executable. To try this process on your own, type the following two commands in the order shown:

	gfortran -c HelloWorld.f90             (This compiles the source code into the object file with a .o extension)
	gfortran HelloWorld.o -o HelloWorld    (This takes the object file, and creates the executable named HelloWorld)

This process of creating object files is useful when multiple source files are involved, but for now, with a single source program like 
HelloWorld, it is a little overkill.

If you decide you want to clean up your directory to the pre-compiled state with the source code file and Makefile only, you can type "make 
clean", and the .o and executable files will be cleaned up by the make utility.
