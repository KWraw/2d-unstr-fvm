OBJECTS = 2d-unstr-fvm-adv

CC = g++
CCFLAGS = -ansi -W -Wall -Wno-write-strings -g -std=c++11
LOCFLAGS =

SRCS = main.cc boundary.cc grid.cc def.cc recon.cc solver.cc test_suite.cc solution.cc
HEADERS = boundary.h grid.h def.h recon.h solver.h test_suite.h solution.h
OBJS = $(SRCS:.cc=.o)

all: $(OBJECTS)

%.o: %.cc %.h
	$(CC) $(CCFLAGS) $(LOCFLAGS) -c $< -o $@

%.o: %.cc
	$(CC) $(CCFLAGS) $(LOCFLAGS) -c $< -o $@

2d-unstr-fvm-adv: $(OBJS) $(HEADERS)
	$(CC) $(CCFLAGS) $(LOCFLAGS) $(OBJS) -o $@

clean:
	rm -rf $(OBJS) $(OBJECTS) *~
