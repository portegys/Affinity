/*
 * This software is provided under the terms of the GNU General
 * Public License as published by the Free Software Foundation.
 *
 * Copyright (c) 2006-2007 Tom Portegys, All Rights Reserved.
 * Permission to use, copy, modify, and distribute this software
 * and its documentation for NON-COMMERCIAL purposes and without
 * fee is hereby granted provided that this copyright notice
 * appears in all copies.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.
 */

/*
 * Chemistry for atoms interacting with charge and valence forces.
 */

#ifndef __CHEMISTRY__
#define __CHEMISTRY__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <errno.h>
#include <vector>
#ifdef THREADS
#include <pthread.h>
#endif
#include "parameters.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "thermal.hpp"
#include "../utility/random.hpp"
#include "../utility/octree.hpp"

namespace affinity
{
class Chemistry
{
public:

   // Parameters.
   Parameters *parameters;

   // Atoms.
   vector<Atom *> atoms;

   // Atom ID factory.
   int atomIDfactory;

   // Vessel size.
   float vesselRadius;

   // Random numbers.
   RANDOM randomSeed;
   Random *randomizer;

   // Atomic bodies.
   vector<OctObject *> bodies;
   Octree              *bodyTracker;

   // Thermal objects.
   vector<Thermal *> thermals;

   // Constructor.
#ifdef THREADS
   Chemistry(float vesselRadius, RANDOM randomSeed, int numThreads);
#else
   Chemistry(float vesselRadius, RANDOM randomSeed);
#endif

   // Initialize chemistry.
   void init(int numAtoms);

   // Destructor.
   ~Chemistry();

   // Clear chemistry.
   void clear();

   // Create atom.
   Atom *createAtom(int protons, int id = -1);

   // Add atom and return ID.
   int addAtom(Atom *);

   // Remove atom.
   void removeAtom(int id);

   // Get atom by ID.
   Atom *getAtom(int id);

   // Create thermal object.
   Thermal *createThermal(float radius, Vector& position, float temperature);

   // Add thermal object.
   void addThermal(Thermal *);

   // Update system.
   void update();

   // Bond updated?
   bool bondUpdate;

   // Mark and count atoms in molecule.
   void clearAtomMarks();
   void markMolecule(Atom *atom, vector<int>& atomCounts, int mark);

   // Generate molecules.
   vector<Molecule *> molecules;
   void generateMolecules();

   // Load and save atoms.
   void load(FILE *fp);
   void save(FILE *fp);
   void import(FILE *fp);

   // Molecule detectors.
   void getMoleculeStats(int& num, int& numClosed,
                         int& numTypes, int& numClosedTypes,
                         float& aveSize, float& aveClosedSize);

#if (O2_MOLECULES)
   int countO2();
#endif
#if (H2O_MOLECULES)
   int countH2O();
#endif
#if (CO2_MOLECULES)
   int countCO2();
#endif

private:

   void update(int threadNum);

#ifdef THREADS
   pthread_barrier_t updateBarrier;
   pthread_mutex_t   updateMutex;
   pthread_t         *threads;
   int               numThreads;
   struct ThreadInfo
   {
      Chemistry *chemistry;
      int       threadNum;
   };
   static void *updateThread(void *threadInfo);

   vector<OctObject *> moveList;
   bool                terminate;
#endif
};
}
#endif
