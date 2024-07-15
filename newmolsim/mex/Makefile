
all=lj.mex leapfrog.mex neighblist.mex calcdr2.mex

all: $(all)

calcdr2.mex: calcdr2.c
	mkoctfile -Wall --mex calcdr2.c

lj.mex: lj.c
	mkoctfile -Wall --mex lj.c 

leapfrog.mex: leapfrog.c
	mkoctfile -Wall --mex leapfrog.c

neighblist.mex: neighblist.c
	mkoctfile -Wall --mex neighblist.c

clean:
	rm -f *.mex
	touch *.c
