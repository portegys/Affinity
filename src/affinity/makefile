# Affinity artificial chemistry makefile.

CC = gcc

CCFLAGS = -DUNIX -DTHREADS -I../../include -O3
LINKLIBS = -L../../lib -lsxmlgui -lglpng -lglut -lGLU -lGL -lm -lpthread -lstdc++

all: external ../../bin/affinity ../../bin/evolve_affinity

../../bin/affinity: affinity.h affinity.cpp \
    ../chemistry/atom.o ../chemistry/chemistry.o ../chemistry/parameters.o \
    ../chemistry/body.o ../chemistry/molecule.o ../chemistry/thermal.o \
    ../utility/baseObject.o ../utility/frustum.o  ../utility/octree.o \
    ../utility/camera.o ../utility/gettime.o ../utility/quaternion.o \
    ../utility/fileio.o ../utility/log.o ../utility/random.o \
    ../utility/frameRate.o ../utility/md5.o ../utility/spacial.o \
    ../../lib/libsxmlgui.a ../../lib/libglpng.a
	$(CC) $(CCFLAGS) -o ../../bin/affinity -DAFFINITY_MAIN affinity.cpp \
        ../chemistry/atom.o ../chemistry/chemistry.o ../chemistry/parameters.o \
        ../chemistry/body.o ../chemistry/molecule.o ../chemistry/thermal.o \
        ../utility/baseObject.o ../utility/frustum.o  ../utility/octree.o \
        ../utility/camera.o ../utility/gettime.o ../utility/quaternion.o \
        ../utility/fileio.o ../utility/log.o ../utility/random.o \
        ../utility/frameRate.o ../utility/md5.o ../utility/spacial.o \
        $(LINKLIBS)

../../bin/evolve_affinity: evolveAffinity.cpp affinity.h affinity.cpp \
    ../chemistry/atom.o ../chemistry/chemistry.o ../chemistry/parameters.o \
    ../chemistry/body.o ../chemistry/molecule.o ../chemistry/thermal.o \
    ../utility/baseObject.o ../utility/frustum.o  ../utility/octree.o \
    ../utility/camera.o ../utility/gettime.o ../utility/quaternion.o \
    ../utility/fileio.o ../utility/log.o ../utility/random.o \
    ../utility/frameRate.o ../utility/md5.o ../utility/spacial.o \
    ../../lib/libsxmlgui.a ../../lib/libglpng.a
	$(CC) $(CCFLAGS) -o ../../bin/evolve_affinity evolveAffinity.cpp affinity.cpp \
        ../chemistry/atom.o ../chemistry/chemistry.o ../chemistry/parameters.o \
        ../chemistry/body.o ../chemistry/molecule.o ../chemistry/thermal.o \
        ../utility/baseObject.o ../utility/frustum.o  ../utility/octree.o \
        ../utility/camera.o ../utility/gettime.o ../utility/quaternion.o \
        ../utility/fileio.o ../utility/log.o ../utility/random.o \
        ../utility/frameRate.o ../utility/md5.o ../utility/spacial.o \
        $(LINKLIBS)

external:
	@(cd ../gui; make)
	@(cd ../utility; make)
	@(cd ../chemistry; make)

clean:
	/bin/rm -f *.o
	@(cd ../gui; make clean)
	@(cd ../utility; make clean)
	@(cd ../chemistry; make clean)



