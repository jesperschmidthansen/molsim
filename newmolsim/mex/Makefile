
all=lj.mex leapfrog.mex

all: $(all)

lj.mex: lj.c
	mkoctfile --mex lj.c 

leapfrog.mex: leapfrog.c
	mkoctfile --mex leapfrog.c

clean:
	rm -f *.mex
	touch *.c
