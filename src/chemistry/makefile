# Build the chemistry.

CC = gcc

CCFLAGS = -DUNIX -DTHREADS -O3

all: parameters.o atom.o body.o molecule.o reaction.o thermal.o chemistry.o

parameters.o: parameters.hpp parameters.cpp
	$(CC) $(CCFLAGS) -c parameters.cpp

atom.o: atom.hpp atom.cpp body.hpp parameters.hpp
	$(CC) $(CCFLAGS) -c atom.cpp
	
body.o: body.hpp body.cpp parameters.hpp
	$(CC) $(CCFLAGS) -c body.cpp

molecule.o: molecule.hpp molecule.cpp
	$(CC) $(CCFLAGS) -c molecule.cpp

reaction.o: reaction.hpp reaction.cpp
	$(CC) $(CCFLAGS) -c reaction.cpp

thermal.o: thermal.hpp thermal.cpp parameters.hpp
	$(CC) $(CCFLAGS) -c thermal.cpp

chemistry.o: chemistry.hpp chemistry.cpp atom.hpp body.hpp thermal.hpp parameters.hpp
	$(CC) $(CCFLAGS) -c chemistry.cpp

clean:
	/bin/rm -f *.o
