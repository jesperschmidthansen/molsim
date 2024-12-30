
CC = mkoctfile
FLAGS = -Wall -Wextra -Ofast --mex

srcs = $(wildcard *.c)
mexs = $(srcs:.c=.mex)

all: $(mexs)

%.mex: %.c
	$(CC) $(FLAGS) $< -o $@

clean:
	rm -f *.mex
	touch *.c
