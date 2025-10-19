
CC = mkoctfile
FLAGS += -Wall -Wextra -Ofast 

srcm = $(wildcard *.c)
mexs = $(srcm:.c=.mex)

srco = $(wildcard *.cpp)
octs = $(srco:.cpp=.oct)


all: $(mexs) $(octs)

%.mex: %.c
	$(CC) $(FLAGS) --mex $< -o $@

%.oct: %.cpp
	$(CC) $(FLAGS) $< -o $@

clean:
	rm -f *.mex *.oct
	touch *.c *.cpp

