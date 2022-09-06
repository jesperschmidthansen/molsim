# Makefile.in for molsim

# Macros

CC          =  gcc
CFLAGS      =  -g -O2 -fgnu89-inline -std=c99 -fPIC -W -Wall -Wextra -Wpedantic -pedantic
RANLIB      =  ranlib
AR          =  ar 
COMPLEX     =  -DCOMPLEX
OMP         =  -fopenmp
OBJECTS     =  sepmisc.o seputil.o separray.o sepinit.o sepprfrc.o sepintgr.o \
	       sepret.o sepmol.o sepcoulomb.o sepsampler.o sepomp.o	


# Compiling molsim
molsim.mex: molsim.c libsep.a _sep_lattice.o _sep_sfg.o
	mkoctfile --mex -DOCTAVE -Ofast -W -pedantic -c -v -Iinclude task.c
	mkoctfile --mex -DOCTAVE -Ofast -W -pedantic -c -v -Iinclude molsim.c
	mkoctfile --mex libsep.a task.o molsim.o _sep_lattice.o _sep_sfg.o -o molsim
	mv molsim.mex ../inst/
	cp ../tests/*.m ../inst/
	cp ../tests/*.xyz ../inst/
	cp ../tests/*.top ../inst/

# Building library
libsep.a: $(OBJECTS)
	$(AR) r libsep.a $(OBJECTS) 
	$(RANLIB) libsep.a

%.o:source/%.c 
	$(CC) $(COMPLEX) $(OMP) -c $(CFLAGS) -W -Wall -Iinclude/ $<


# Compiling tools
_sep_lattice.o: tools/_sep_lattice.c
	$(CC) $(CFLAGS) $(OMP) -c tools/_sep_lattice.c

_sep_sfg.o: tools/_sep_sfg.c
	$(CC) $(CFLAGS) $(OMP) -c tools/_sep_sfg.c

# Cleaning up
clean: 
	rm -f libsep.a *.o
	rm -f source/*~
	rm -f include/*~
	rm -f *~
	rm -f $(TOOLS)
	rm -f ../inst/molsim.mex

