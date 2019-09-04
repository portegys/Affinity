Affinity artificial chemistry

An architecture for an artificial chemistry featuring 3D
continuous physics and and chemical covalent bonding.

The "Affinity" artificial chemistry system is a search for interesting (e.g. auto-catalytic) chemical reactions and systems. 
It is loosely based on the VSEPR (Valence Shell Electron Pair Repulsion) chemistry model. Affinity has a number of tunable 
parameters, e.g. electron mass, orbital radii, etc.

Required packages:

OpenGL graphics and the GLUT package are required to build and
run the program. They can be obtained from the Mesa project:
www.mesa3d.org. The UNIX version also requires the gcc compiler,
the make command, and the bash shell. The Windows version requires
the Microsoft Visual Studio 2015 (or later) IDE.

To build:
UNIX: make
Windows: use VS solution.

To run:
Run the executables in the bin folder.

Usage:

affinity
      [-cycles <number of cycles>]
      [-numAtoms <number of atoms>]
      [-numThreads <number of threads (default=1)>]
      [-vesselRadius <vessel radius>]
      [-thermal <radius>,<x>,<y>,<z>,<temperature>] (multiple option)
      [-randomSeed <random seed>]
      [-load <load file name>]
      [-save <save file name>]
      [-dump <molecule dump file name> ('-' for stdout)]
      [-noGraphics (turn off graphics)]
      [-statsFreq <statistics gather frequency (0=never, 1=default)>]
      
e.g., affinity -numAtoms 10 -thermal 3.0,0.0,0.0,0.0,5.0

To run evolve_affinity chemistry parameter evolution program,
see usage in evolveAffinity.cpp file.



