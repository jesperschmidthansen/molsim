
all=lj.mex leapfrog.mex neighblist.mex

all: $(all)

lj.mex: lj.c
	mkoctfile -Wall --mex lj.c 

leapfrog.mex: leapfrog.c
	mkoctfile -Wall --mex leapfrog.c

neighblist.mex: neighblist.c
	mkoctfile -Wall --mex neighblist.c

clean:
	rm -f *.mex
	touch *.c
