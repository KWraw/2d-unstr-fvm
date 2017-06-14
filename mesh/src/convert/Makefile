CC    =  gcc
CFLAG =  -g -lm

OBJECTS = $(patsubst %.c,%.o,$(wildcard *.c))
HEADERS = grid.h

all : convert

.c.o :
	$(CC) $(CFLAG) -c $<

convert : $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $^ $(CFLAG)

clean :
	rm -f *.o *.vtk linadv

.PHONY : clean
