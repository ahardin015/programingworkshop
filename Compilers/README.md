Compilers
==================

We use several compilers here. GNU and Intel are most common.
---

GNU
---
gcc gfortran g++

http://gcc.gnu.org/onlinedocs/gcc-4.8.1/gcc/

Example
---

gcc -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -m32 -std=gnu99 -O2 -g -o editor editor.o -I. -I${DREADER_PATH}/inc -L. -L${DREADER_PATH}/lib -lm -ldl -lc

Compiler flags
---

-g - adds debugging support with gdb

-Wall - print all warnings about weird looking code (Good practice to use this)

-L - path to libraries

-I - path to include files

-l - link to a library

-O - optimization 0,1,2,3

-o name - create the program with this name. If no -o is given a.out is usually used.

-c - compiles into an object. Usefull when needing to compile against several files.

--version - prints version info. Usefull when compiling against netcdf or other.

-v - lots of info

--help

Intel
---

icc ifort icpc

http://software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/cpp/win/

Compiler flags
---

Many are similar to above

-x compiles instructions for a processor type. Improve performance on intel cpus. SIMD optimizations.

-m compiles instructions for a processor type.


If problems occur, many options to try to fix.

