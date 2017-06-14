CC    =  gcc
CFLAG =  -g -lm

OBJECTS = $(patsubst %.c,%.o,$(wildcard *.c))
HEADERS = grid.h timer.h cloption.h

all : rect_trian_mesh

.c.o :
	$(CC) $(CFLAG) -c $<

rect_trian_mesh : $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $^ $(CFLAG)

clean :
	rm -f *.o *.vtk rect_trian_mesh

.PHONY : clean
