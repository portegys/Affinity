# Common makefile

CC = gcc

CCFLAGS = -DUNIX -O3

all: baseObject.o camera.o fileio.o frameRate.o frustum.o gettime.o log.o md5.o octree.o quaternion.o random.o spacial.o

baseObject.o: baseObject.hpp baseObject.cpp spacial.hpp
	$(CC) $(CCFLAGS) -c baseObject.cpp

camera.o: camera.hpp camera.cpp baseObject.hpp spacial.hpp
	$(CC) $(CCFLAGS) -c camera.cpp

fileio.o: fileio.h fileio.cpp
	$(CC) $(CCFLAGS) -c fileio.cpp

frameRate.o: frameRate.hpp frameRate.cpp
	$(CC) $(CCFLAGS) -c frameRate.cpp

frustum.o: frustum.hpp frustum.cpp
	$(CC) $(CCFLAGS) -c frustum.cpp

gettime.o: gettime.h gettime.cpp
	$(CC) $(CCFLAGS) -c gettime.cpp
	
log.o: log.cpp log.hpp
	$(CC) $(CCFLAGS) -c log.cpp

md5.o: md5.cpp md5.h
	$(CC) $(CCFLAGS) -c md5.cpp

octree.o: octree.hpp octree.cpp
	$(CC) $(CCFLAGS) -c octree.cpp

quaternion.o: quaternion.cpp quaternion.hpp
	$(CC) $(CCFLAGS) -c quaternion.cpp

random.o: random.hpp random.cpp
	$(CC) $(CCFLAGS) -c random.cpp

spacial.o: spacial.hpp spacial.cpp
	$(CC) $(CCFLAGS) -c spacial.cpp

clean:
	/bin/rm -f *.o
