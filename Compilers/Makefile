# Makefile for radar software.
# System executable definitions

#DREADER_PATH = /home/curtis/research/software/dreader
DREADER_PATH = /Users/aereinha/Downloads/dreader/dreader
MAKE = make
CD = /bin/cd -f
CP = /bin/cp -f
MV = /bin/mv -f
RM = /bin/rm -f

# Compiler and archiver options
CC = gcc
#CFLAGS = -DHAVE_BYTESWAP_H -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -m32 -std=gnu99
CFLAGS = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -m32 -std=gnu99
#COPTS = -g
COPTS = -O2 -g
AR = ar -rscv

# Libraries
STDLIBS = -lm -ldl -lc
DREADERLIBS = -ldorade -lpiraq -lfileio -ldorade -lpiraq -lalgor -lfileio -lalgor
LIBS =  ${DREADERLIBS} ${STDLIBS}

# Directories
VPATH = src/algor:src/dorade:src/fileio:src/piraq:src/drivers
INCDIRS = -I. -I${DREADER_PATH}/inc
LIBDIRS = -L. -L${DREADER_PATH}/lib

# Compilation rules
.c.o:
	$(CC) $(CFLAGS) $(COPTS) $(INCDIRS) -c $<

# Object file lists
ALGOR_OBJS = algorithms.o coordinates.o interpolate.o procmatrix.o procmodel.o qualitycontrol.o statistics.o
DORADE_OBJS = procdorade.o readdorade.o writedorade.o
FILEIO_OBJS = allocate.o bswap.o cmdlineargs.o kml.o printout.o readfile.o readparam.o writefile.o
PIRAQ_OBJS = readpiraq.o writepiraq.o
DREADER_OBJS = ${ALGOR_OBJS} ${FILEIO_OBJS} ${DORADE_OBJS} ${PIRAQ_OBJS}

# Header file lists
DEFINE_HEADS = calcstructs.h doradedefs.h doradestructs.h modelstructs.h piraqdefs.h piraqstructs.h piraqxdefs.h
ALGOR_HEADS = algorithms.h coordinates.h interpolate.h procmatrix.h procmodel.h qualitycontrol.h statistics.h
DORADE_HEADS = procdorade.h readdorade.h writedorade.h
FILEIO_HEADS = allocate.h bswap.h cmdlineargs.h kml.h printout.h readfile.h readparam.h writefile.h
PIRAQ_HEADS = readpiraq.h writepiraq.h
DREADER_HEADS = ${DEFINE_HEADS} ${ALGOR_HEADS} ${FILEIO_HEADS} ${DORADE_HEADS} ${PIRAQ_HEADS}

# Targets
libs: algorlib doradelib fileiolib piraqlib
progs: editor extractor parser plotter rewriter tormodel translator
all: libs progs clean
debug: libs progs

algorlib: $(ALGOR_OBJS)
	$(AR) libalgor.a $(ALGOR_OBJS)
	$(MV) libalgor.a $(DREADER_PATH)/lib
	$(RM) $(ALGOR_OBJS)

doradelib: $(DORADE_OBJS)
	$(AR) libdorade.a $(DORADE_OBJS)
	$(MV) libdorade.a $(DREADER_PATH)/lib
	$(RM) $(DORADE_OBJS)

fileiolib: $(FILEIO_OBJS)
	$(AR) libfileio.a $(FILEIO_OBJS)
	$(MV) libfileio.a $(DREADER_PATH)/lib
	$(RM) $(FILEIO_OBJS)

piraqlib: $(PIRAQ_OBJS)
	$(AR) libpiraq.a $(PIRAQ_OBJS)
	$(MV) libpiraq.a $(DREADER_PATH)/lib
	$(RM) $(PIRAQ_OBJS)

editor: editor.o
	$(RM) $@
	$(CC) $(CFLAGS) $(COPTS) -o editor editor.o $(LIBDIRS) $(LIBS)
	$(MV) $@ $(DREADER_PATH)/bin

extractor: extractor.o
	$(RM) $@
	$(CC) $(CFLAGS) $(COPTS) -o extractor extractor.o $(LIBDIRS) $(LIBS)
	$(MV) $@ $(DREADER_PATH)/bin

parser: parser.o
	$(RM) $@
	$(CC) $(CFLAGS) $(COPTS) -o parser parser.o $(LIBDIRS) $(LIBS)
	$(MV) $@ $(DREADER_PATH)/bin

plotter: plotter.o
	$(RM) $@
	$(CC) $(CFLAGS) $(COPTS) -o plotter plotter.o $(LIBDIRS) $(LIBS)
	$(MV) $@ $(DREADER_PATH)/bin

rewriter: rewriter.o
	$(RM) $@
	$(CC) $(CFLAGS) $(COPTS) -o rewriter rewriter.o $(LIBDIRS) $(LIBS)
	$(MV) $@ $(DREADER_PATH)/bin

tormodel: tormodel.o
	$(RM) $@
	$(CC) $(CFLAGS) $(COPTS) -o tormodel tormodel.o $(LIBDIRS) $(LIBS)
	$(MV) $@ $(DREADER_PATH)/bin

translator: translator.o
	$(RM) $@
	$(CC) $(CFLAGS) $(COPTS) -o translator translator.o $(LIBDIRS) $(LIBS)
	$(MV) $@ $(DREADER_PATH)/bin

clean:
	$(RM) *.o
