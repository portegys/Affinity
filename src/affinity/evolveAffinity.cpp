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

/**
 * Evolve Affinity chemistry.
 */

#ifdef WIN32
#include <windows.h>
#endif
#ifdef UNIX
#include <errno.h>
#endif
#include "affinity.h"
#include "../utility/log.hpp"

// Usage.
char *Usage[] =
{
   (char *)"To run:",
   (char *)"  evolve_affinity",
   (char *)"      -generations <evolution generations>",
   (char *)"      -cycles (cycles per run)",
#ifdef THREADS
   (char *)"      [-numThreads <number of threads (default=1)>]",
#endif
   (char *)"      [-input <evolution input file name> (for run continuation)]",
   (char *)"      -output <evolution output file name>",
   (char *)"      [-randomSeed <random seed> (for new run)]",
   (char *)"      [-logfile <log file name>]",
   (char *)"To unpack chemistries in affinity format:",
   (char *)"  evolve_affinity",
   (char *)"      -unpack",
   (char *)"      [-prefix <chemistry file name prefix>]",
   (char *)"      -input <evolution input file name>",
   NULL
};

// Evolvable parameter ranges.
#define MIN_VESSEL_RADIUS                    5.0f
#define MAX_VESSEL_RADIUS                    25.0f
#define MIN_ATOMS                            20
#define MAX_ATOMS                            100
#define MIN_MIN_ATOM_INITIAL_FORCE           0.0f
#define MAX_MIN_ATOM_INITIAL_FORCE           0.5f
#define MIN_MAX_ATOM_INITIAL_FORCE           1.0f
#define MAX_MAX_ATOM_INITIAL_FORCE           25.0f
#define MIN_MIN_NUCLEUS_PROTONS              1
#define MAX_MIN_NUCLEUS_PROTONS              1
#define MIN_MAX_NUCLEUS_PROTONS              Parameters::DEFAULT_MAX_NUCLEUS_PROTONS
#define MAX_MAX_NUCLEUS_PROTONS              Parameters::DEFAULT_MAX_NUCLEUS_PROTONS
#define MIN_PROTON_MASS                      0.1f
#define MAX_PROTON_MASS                      10.0f
#define MIN_ELECTRON_MASS                    0.01f
#define MAX_ELECTRON_MASS                    10.0f
#define MIN_PROTON_CHARGE                    0.01f
#define MAX_PROTON_CHARGE                    5.0f
#define MIN_ELECTRON_CHARGE                  -0.01f
#define MAX_ELECTRON_CHARGE                  -5.0f
#define MIN_CHARGE_GAUSSIAN_SPREAD           0.5f
#define MAX_CHARGE_GAUSSIAN_SPREAD           10.0f
#define MIN_NUCLEUS_BODY_RADIUS              0.1f
#define MAX_NUCLEUS_BODY_RADIUS              1.0f
#define MIN_ORBITAL_BODY_RADIUS              0.1f
#define MAX_ORBITAL_BODY_RADIUS              1.0f
#define MIN_MAX_BODY_RANGE                   0.1f
#define MAX_MAX_BODY_RANGE                   20.0f
#define MIN_BOND_LENGTH                      0.1f
#define MAX_BOND_LENGTH                      5.0f
#define MIN_BOND_STIFFNESS                   0.1f
#define MAX_BOND_STIFFNESS                   5.0f
#define MIN_BOND_DAMPER                      0.001f
#define MAX_BOND_DAMPER                      0.05f
#define MIN_NUCLEAR_REPULSION_STIFFNESS      0.1f
#define MAX_NUCLEAR_REPULSION_STIFFNESS      4.0f
#define MIN_COVALENT_BONDING_RANGE           0.5f
#define MAX_COVALENT_BONDING_RANGE           3.0f
#define MIN_MIN_COVALENT_BOND_FORCE          0.1f
#define MAX_MIN_COVALENT_BOND_FORCE          4.0f
#define MIN_COVALENT_BOND_STIFFNESS_SCALE    0.1f
#define MAX_COVALENT_BOND_STIFFNESS_SCALE    2.0f
#define MIN_THERMALS                         0
#define MAX_THERMALS                         5
#define MIN_MIN_THERMAL_RADIUS               0.1f
#define MAX_MIN_THERMAL_RADIUS               0.5f
#define MIN_MAX_THERMAL_RADIUS               1.0f
#define MAX_MAX_THERMAL_RADIUS               5.0f
#define MIN_MIN_THERMAL_TEMPERATURE          0.0f
#define MAX_MIN_THERMAL_TEMPERATURE          1.0f
#define MIN_MAX_TEMPERATURE                  10.0f
#define MAX_MAX_TEMPERATURE                  50.0f
#define MIN_UPDATE_STEP                      0.01f
#define MAX_UPDATE_STEP                      1.0f

// Evolution parameters.
#define FIT_POPULATION_SIZE                  20
#define NUM_MUTANTS                          10
#define NUM_OFFSPRING                        10
#define POPULATION_SIZE                      (FIT_POPULATION_SIZE + NUM_MUTANTS + NUM_OFFSPRING)
#define MUTATION_RATE                        0.1
#define DEFAULT_CYCLES                       500
#define EVOLVE_LOGGING                       LOG_TO_FILE
#define DEFAULT_EVOLVE_LOG_FILE_NAME         "evolve.log"
#define SAVE_FREQUENCY                       1

// Print usage.
void printUsage()
{
   int i;

   if ((EVOLVE_LOGGING == LOG_TO_FILE) || (EVOLVE_LOGGING == NO_LOG))
   {
      for (i = 0; Usage[i] != NULL; i++)
      {
         fprintf(stderr, "%s\n", Usage[i]);
      }
   }
   for (i = 0; Usage[i] != NULL; i++)
   {
      sprintf(Log::messageBuf, Usage[i]);
      Log::logInformation();
   }
}


// Evolution generations.
int Generations;

// Population file names.
char *InputFileName;
char *OutputFileName;

// Unpacking mode.
bool Unpack;
char *UnpackPrefix;

// Population member.
class Member
{
public:

   Chemistry *chemistry;
   float     fitness;
   int       generation;

   // Constructors.
   Member()
   {
      chemistry  = NULL;
      fitness    = 0.0f;
      generation = -1;
   }


   Member(int generation)
   {
      RANDOM s;
      int    i, n;
      Vector p;
      float  r, t, v;

      fitness          = 0.0f;
      this->generation = generation;
      s = Randomizer->RAND();
      Randomizer->RAND_PUSH();
      Randomizer->SRAND(s);
      n = Randomizer->RAND_CHOICE(MAX_ATOMS - MIN_ATOMS + 1) + MIN_ATOMS;
      if (generation == 0)
      {
         v = DEFAULT_VESSEL_RADIUS;
#ifdef THREADS
         chemistry = new Chemistry(v, s, NumThreads);
#else
         chemistry = new Chemistry(v, s);
#endif
         assert(chemistry != NULL);
         chemistry->init(n);
      }
      else
      {
         v = Randomizer->RAND_INTERVAL(MIN_VESSEL_RADIUS, MAX_VESSEL_RADIUS);
#ifdef THREADS
         chemistry = new Chemistry(v, s, NumThreads);
#else
         chemistry = new Chemistry(v, s);
#endif
         assert(chemistry != NULL);
         chemistry->parameters->MIN_ATOM_INITIAL_FORCE =
            Randomizer->RAND_INTERVAL(MIN_MIN_ATOM_INITIAL_FORCE, MAX_MIN_ATOM_INITIAL_FORCE);
         chemistry->parameters->MAX_ATOM_INITIAL_FORCE =
            Randomizer->RAND_INTERVAL(MIN_MAX_ATOM_INITIAL_FORCE, MAX_MAX_ATOM_INITIAL_FORCE);
         chemistry->parameters->MIN_NUCLEUS_PROTONS =
            Randomizer->RAND_CHOICE(MAX_MIN_NUCLEUS_PROTONS - MIN_MIN_NUCLEUS_PROTONS + 1) + MIN_MIN_NUCLEUS_PROTONS;
         chemistry->parameters->MAX_NUCLEUS_PROTONS =
            Randomizer->RAND_CHOICE(MAX_MAX_NUCLEUS_PROTONS - MIN_MAX_NUCLEUS_PROTONS + 1) + MIN_MAX_NUCLEUS_PROTONS;
         chemistry->parameters->PROTON_MASS =
            Randomizer->RAND_INTERVAL(MIN_PROTON_MASS, MAX_PROTON_MASS);
         chemistry->parameters->ELECTRON_MASS =
            Randomizer->RAND_INTERVAL(MIN_ELECTRON_MASS, MAX_ELECTRON_MASS);
         chemistry->parameters->PROTON_CHARGE =
            Randomizer->RAND_INTERVAL(MIN_PROTON_CHARGE, MAX_PROTON_CHARGE);
         chemistry->parameters->ELECTRON_CHARGE =
            Randomizer->RAND_INTERVAL(MIN_ELECTRON_CHARGE, MAX_ELECTRON_CHARGE);
         chemistry->parameters->CHARGE_GAUSSIAN_SPREAD =
            Randomizer->RAND_INTERVAL(MIN_CHARGE_GAUSSIAN_SPREAD, MAX_CHARGE_GAUSSIAN_SPREAD);
         chemistry->parameters->NUCLEUS_BODY_RADIUS =
            Randomizer->RAND_INTERVAL(MIN_NUCLEUS_BODY_RADIUS, MAX_NUCLEUS_BODY_RADIUS);
         chemistry->parameters->ORBITAL_BODY_RADIUS =
            Randomizer->RAND_INTERVAL(MIN_ORBITAL_BODY_RADIUS, MAX_ORBITAL_BODY_RADIUS);
         chemistry->parameters->MAX_BODY_RANGE =
            Randomizer->RAND_INTERVAL(MIN_MAX_BODY_RANGE, MAX_MAX_BODY_RANGE);
         chemistry->parameters->BOND_LENGTH =
            Randomizer->RAND_INTERVAL(MIN_BOND_LENGTH, MAX_BOND_LENGTH);
         chemistry->parameters->BOND_STIFFNESS =
            Randomizer->RAND_INTERVAL(MIN_BOND_STIFFNESS, MAX_BOND_STIFFNESS);
         chemistry->parameters->BOND_DAMPER =
            Randomizer->RAND_INTERVAL(MIN_BOND_DAMPER, MAX_BOND_DAMPER);
         chemistry->parameters->NUCLEAR_REPULSION_STIFFNESS =
            Randomizer->RAND_INTERVAL(MIN_NUCLEAR_REPULSION_STIFFNESS, MAX_NUCLEAR_REPULSION_STIFFNESS);
         chemistry->parameters->COVALENT_BONDING_RANGE =
            Randomizer->RAND_INTERVAL(MIN_COVALENT_BONDING_RANGE, MAX_COVALENT_BONDING_RANGE);
         chemistry->parameters->MIN_COVALENT_BOND_FORCE =
            Randomizer->RAND_INTERVAL(MIN_MIN_COVALENT_BOND_FORCE, MAX_MIN_COVALENT_BOND_FORCE);
         chemistry->parameters->COVALENT_BOND_STIFFNESS_SCALE =
            Randomizer->RAND_INTERVAL(MIN_COVALENT_BOND_STIFFNESS_SCALE, MAX_COVALENT_BOND_STIFFNESS_SCALE);
         chemistry->parameters->MIN_THERMAL_RADIUS =
            Randomizer->RAND_INTERVAL(MIN_MIN_THERMAL_RADIUS, MAX_MIN_THERMAL_RADIUS);
         chemistry->parameters->MAX_THERMAL_RADIUS =
            Randomizer->RAND_INTERVAL(MIN_MAX_THERMAL_RADIUS, MAX_MAX_THERMAL_RADIUS);
         chemistry->parameters->MIN_THERMAL_TEMPERATURE =
            Randomizer->RAND_INTERVAL(MIN_MIN_THERMAL_TEMPERATURE, MAX_MIN_THERMAL_TEMPERATURE);
         chemistry->parameters->MAX_TEMPERATURE =
            Randomizer->RAND_INTERVAL(MIN_MAX_TEMPERATURE, MAX_MAX_TEMPERATURE);
         chemistry->parameters->UPDATE_STEP =
            Randomizer->RAND_INTERVAL(MIN_UPDATE_STEP, MAX_UPDATE_STEP);
         chemistry->init(n);
#if ( ORGANIC_MOLECULES )
         // Install thermals at top and bottom of vessel.
         n = 2;
#else
         n = Randomizer->RAND_CHOICE(MAX_THERMALS - MIN_THERMALS + 1) + MIN_THERMALS;
#endif
         for (i = 0; i < n; i++)
         {
            r = Randomizer->RAND_INTERVAL(chemistry->parameters->MIN_THERMAL_RADIUS,
                                          chemistry->parameters->MAX_THERMAL_RADIUS);
            t = Randomizer->RAND_INTERVAL(chemistry->parameters->MIN_THERMAL_TEMPERATURE,
                                          chemistry->parameters->MAX_TEMPERATURE);
#if (ORGANIC_MOLECULES)
            if (n == 0)
            {
               p.Zero();
               p.y = v - r;
            }
            else
            {
               p.Zero();
               p.y = r - v;
            }
#else
            p.x = Randomizer->RAND_INTERVAL(-1.0f, 1.0f);
            p.y = Randomizer->RAND_INTERVAL(-1.0f, 1.0f);
            p.z = Randomizer->RAND_INTERVAL(-1.0f, 1.0f);
            p.Normalize(Randomizer->RAND_INTERVAL(0.0f, v - r));
#endif
            chemistry->createThermal(r, p, t);
         }
      }
      Randomizer->RAND_POP();
   }


   // Destructor.
   ~Member()
   {
      if (chemistry != NULL) { delete chemistry; }
   }


   // Evaluate.
   void evaluate()
   {
      for (int i = 0; i < Cycles; i++)
      {
         chemistry->update();
      }
      fitness = 0.0f;

#if (O2_MOLECULES)
      // Count O2 molecules.
      int o2 = chemistry->countO2();
      if (o2 > 0)
      {
         fitness += 100.0f + (float)(o2 - 1);
      }
#endif

#if (H2O_MOLECULES)
      // Count H2O molecules.
      int h2o = chemistry->countH2O();
      if (h2o > 0)
      {
         fitness += 100.0f + (float)(h2o - 1);
      }
#endif

#if (CO2_MOLECULES)
      // Count CO2 molecules.
      int co2 = chemistry->countCO2();
      if (co2 > 0)
      {
         fitness += 100.0f + (float)(co2 - 1);
      }
#endif

#if (ORGANIC_MOLECULES)
      int   num, numClosed, numTypes, numClosedTypes;
      float aveSize, aveClosedSize;
      chemistry->getMoleculeStats(num, numClosed, numTypes,
                                  numClosedTypes, aveSize, aveClosedSize);
      if (numClosed > 0)
      {
         fitness = ((float)numClosed * aveClosedSize) /
                   (float)numClosedTypes;
      }
#endif

#if (!H2O_MOLECULES && !O2_MOLECULES && !CO2_MOLECULES && !ORGANIC_MOLECULES)
      // Fitness is number of different-sized molecules.
      int         i, i2, j, j2, k, c;
      vector<int> atomCounts;
      vector<int> moleculeCounts;
      chemistry->clearAtomMarks();
      moleculeCounts.resize(chemistry->atoms.size());
      for (i = 0, i2 = moleculeCounts.size(); i < i2; i++) { moleculeCounts[i] = 0; }
      atomCounts.resize(chemistry->parameters->MAX_NUCLEUS_PROTONS);
      for (i = k = 0, i2 = chemistry->atoms.size(); i < i2; i++)
      {
         for (j = 0, j2 = atomCounts.size(); j < j2; j++)
         {
            atomCounts[j] = 0;
         }
         if (chemistry->atoms[i]->mark == -1)
         {
            chemistry->markMolecule(chemistry->atoms[i], atomCounts, k);
            for (j = c = 0, j2 = atomCounts.size(); j < j2; j++)
            {
               c += atomCounts[j];
            }
            moleculeCounts[c]++;
            k++;
         }
      }
      for (j = c = 0, j2 = moleculeCounts.size(); j < j2; j++)
      {
         if (moleculeCounts[j] > 0) { c++; }
      }
      fitness += (float)c;
#endif
   }


   // Create mutation of member.
   void mutate(Member *member)
   {
      RANDOM s;
      int    i, n;
      Vector p;
      float  r, t, v;

      assert(member->chemistry != NULL);
      fitness    = 0.0f;
      generation = member->generation + 1;
      s          = Randomizer->RAND();
      Randomizer->RAND_PUSH();
      Randomizer->SRAND(s);
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         v = Randomizer->RAND_INTERVAL(MIN_VESSEL_RADIUS, MAX_VESSEL_RADIUS);
      }
      else
      {
         v = member->chemistry->vesselRadius;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         n = Randomizer->RAND_CHOICE(MAX_ATOMS - MIN_ATOMS + 1) + MIN_ATOMS;
      }
      else
      {
         n = (int)member->chemistry->atoms.size();
      }
      if (chemistry != NULL) { delete chemistry; }
#ifdef THREADS
      chemistry = new Chemistry(v, s, NumThreads);
#else
      chemistry = new Chemistry(v, s);
#endif
      assert(chemistry != NULL);
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MIN_ATOM_INITIAL_FORCE =
            Randomizer->RAND_INTERVAL(MIN_MIN_ATOM_INITIAL_FORCE, MAX_MIN_ATOM_INITIAL_FORCE);
      }
      else
      {
         chemistry->parameters->MIN_ATOM_INITIAL_FORCE =
            member->chemistry->parameters->MIN_ATOM_INITIAL_FORCE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MAX_ATOM_INITIAL_FORCE =
            Randomizer->RAND_INTERVAL(MIN_MAX_ATOM_INITIAL_FORCE, MAX_MAX_ATOM_INITIAL_FORCE);
      }
      else
      {
         chemistry->parameters->MAX_ATOM_INITIAL_FORCE =
            member->chemistry->parameters->MAX_ATOM_INITIAL_FORCE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MIN_NUCLEUS_PROTONS =
            Randomizer->RAND_CHOICE(MAX_MIN_NUCLEUS_PROTONS - MIN_MIN_NUCLEUS_PROTONS + 1) + MIN_MIN_NUCLEUS_PROTONS;
      }
      else
      {
         chemistry->parameters->MIN_NUCLEUS_PROTONS =
            member->chemistry->parameters->MIN_NUCLEUS_PROTONS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MAX_NUCLEUS_PROTONS =
            Randomizer->RAND_CHOICE(MAX_MAX_NUCLEUS_PROTONS - MIN_MAX_NUCLEUS_PROTONS + 1) + MIN_MAX_NUCLEUS_PROTONS;
      }
      else
      {
         chemistry->parameters->MAX_NUCLEUS_PROTONS =
            member->chemistry->parameters->MAX_NUCLEUS_PROTONS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->PROTON_MASS =
            Randomizer->RAND_INTERVAL(MIN_PROTON_MASS, MAX_PROTON_MASS);
      }
      else
      {
         chemistry->parameters->PROTON_MASS =
            member->chemistry->parameters->PROTON_MASS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->ELECTRON_MASS =
            Randomizer->RAND_INTERVAL(MIN_ELECTRON_MASS, MAX_ELECTRON_MASS);
      }
      else
      {
         chemistry->parameters->ELECTRON_MASS =
            member->chemistry->parameters->ELECTRON_MASS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->PROTON_CHARGE =
            Randomizer->RAND_INTERVAL(MIN_PROTON_CHARGE, MAX_PROTON_CHARGE);
      }
      else
      {
         chemistry->parameters->PROTON_CHARGE =
            member->chemistry->parameters->PROTON_CHARGE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->ELECTRON_CHARGE =
            Randomizer->RAND_INTERVAL(MIN_ELECTRON_CHARGE, MAX_ELECTRON_CHARGE);
      }
      else
      {
         chemistry->parameters->ELECTRON_CHARGE =
            member->chemistry->parameters->ELECTRON_CHARGE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->CHARGE_GAUSSIAN_SPREAD =
            Randomizer->RAND_INTERVAL(MIN_CHARGE_GAUSSIAN_SPREAD, MAX_CHARGE_GAUSSIAN_SPREAD);
      }
      else
      {
         chemistry->parameters->CHARGE_GAUSSIAN_SPREAD =
            member->chemistry->parameters->CHARGE_GAUSSIAN_SPREAD;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->NUCLEUS_BODY_RADIUS =
            Randomizer->RAND_INTERVAL(MIN_NUCLEUS_BODY_RADIUS, MAX_NUCLEUS_BODY_RADIUS);
      }
      else
      {
         chemistry->parameters->NUCLEUS_BODY_RADIUS =
            member->chemistry->parameters->NUCLEUS_BODY_RADIUS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->ORBITAL_BODY_RADIUS =
            Randomizer->RAND_INTERVAL(MIN_ORBITAL_BODY_RADIUS, MAX_ORBITAL_BODY_RADIUS);
      }
      else
      {
         chemistry->parameters->ORBITAL_BODY_RADIUS =
            member->chemistry->parameters->ORBITAL_BODY_RADIUS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MAX_BODY_RANGE =
            Randomizer->RAND_INTERVAL(MIN_MAX_BODY_RANGE, MAX_MAX_BODY_RANGE);
      }
      else
      {
         chemistry->parameters->MAX_BODY_RANGE =
            member->chemistry->parameters->MAX_BODY_RANGE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->BOND_LENGTH =
            Randomizer->RAND_INTERVAL(MIN_BOND_LENGTH, MAX_BOND_LENGTH);
      }
      else
      {
         chemistry->parameters->BOND_LENGTH =
            member->chemistry->parameters->BOND_LENGTH;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->BOND_STIFFNESS =
            Randomizer->RAND_INTERVAL(MIN_BOND_STIFFNESS, MAX_BOND_STIFFNESS);
      }
      else
      {
         chemistry->parameters->BOND_STIFFNESS =
            member->chemistry->parameters->BOND_STIFFNESS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->BOND_DAMPER =
            Randomizer->RAND_INTERVAL(MIN_BOND_DAMPER, MAX_BOND_DAMPER);
      }
      else
      {
         chemistry->parameters->BOND_DAMPER =
            member->chemistry->parameters->BOND_DAMPER;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->NUCLEAR_REPULSION_STIFFNESS =
            Randomizer->RAND_INTERVAL(MIN_NUCLEAR_REPULSION_STIFFNESS, MAX_NUCLEAR_REPULSION_STIFFNESS);
      }
      else
      {
         chemistry->parameters->NUCLEAR_REPULSION_STIFFNESS =
            member->chemistry->parameters->NUCLEAR_REPULSION_STIFFNESS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->COVALENT_BONDING_RANGE =
            Randomizer->RAND_INTERVAL(MIN_COVALENT_BONDING_RANGE, MAX_COVALENT_BONDING_RANGE);
      }
      else
      {
         chemistry->parameters->COVALENT_BONDING_RANGE =
            member->chemistry->parameters->COVALENT_BONDING_RANGE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MIN_COVALENT_BOND_FORCE =
            Randomizer->RAND_INTERVAL(MIN_MIN_COVALENT_BOND_FORCE, MAX_MIN_COVALENT_BOND_FORCE);
      }
      else
      {
         chemistry->parameters->MIN_COVALENT_BOND_FORCE =
            member->chemistry->parameters->MIN_COVALENT_BOND_FORCE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->COVALENT_BOND_STIFFNESS_SCALE =
            Randomizer->RAND_INTERVAL(MIN_COVALENT_BOND_STIFFNESS_SCALE, MAX_COVALENT_BOND_STIFFNESS_SCALE);
      }
      else
      {
         chemistry->parameters->COVALENT_BOND_STIFFNESS_SCALE =
            member->chemistry->parameters->COVALENT_BOND_STIFFNESS_SCALE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MIN_THERMAL_RADIUS =
            Randomizer->RAND_INTERVAL(MIN_MIN_THERMAL_RADIUS, MAX_MIN_THERMAL_RADIUS);
      }
      else
      {
         chemistry->parameters->MIN_THERMAL_RADIUS =
            member->chemistry->parameters->MIN_THERMAL_RADIUS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MAX_THERMAL_RADIUS =
            Randomizer->RAND_INTERVAL(MIN_MAX_THERMAL_RADIUS, MAX_MAX_THERMAL_RADIUS);
      }
      else
      {
         chemistry->parameters->MAX_THERMAL_RADIUS =
            member->chemistry->parameters->MAX_THERMAL_RADIUS;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MIN_THERMAL_TEMPERATURE =
            Randomizer->RAND_INTERVAL(MIN_MIN_THERMAL_TEMPERATURE, MAX_MIN_THERMAL_TEMPERATURE);
      }
      else
      {
         chemistry->parameters->MIN_THERMAL_TEMPERATURE =
            member->chemistry->parameters->MIN_THERMAL_TEMPERATURE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->MAX_TEMPERATURE =
            Randomizer->RAND_INTERVAL(MIN_MAX_TEMPERATURE, MAX_MAX_TEMPERATURE);
      }
      else
      {
         chemistry->parameters->MAX_TEMPERATURE =
            member->chemistry->parameters->MAX_TEMPERATURE;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->UPDATE_STEP =
            Randomizer->RAND_INTERVAL(MIN_UPDATE_STEP, MAX_UPDATE_STEP);
      }
      else
      {
         chemistry->parameters->UPDATE_STEP =
            member->chemistry->parameters->UPDATE_STEP;
      }
      chemistry->init(n);
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         n = Randomizer->RAND_CHOICE(MAX_THERMALS - MIN_THERMALS + 1) + MIN_THERMALS;
      }
      else
      {
         n = (int)member->chemistry->thermals.size();
      }
      for (i = 0; i < n; i++)
      {
         r = Randomizer->RAND_INTERVAL(chemistry->parameters->MIN_THERMAL_RADIUS,
                                       chemistry->parameters->MAX_THERMAL_RADIUS);
         t = Randomizer->RAND_INTERVAL(chemistry->parameters->MIN_THERMAL_TEMPERATURE,
                                       chemistry->parameters->MAX_TEMPERATURE);
         p.x = Randomizer->RAND_INTERVAL(-1.0f, 1.0f);
         p.y = Randomizer->RAND_INTERVAL(-1.0f, 1.0f);
         p.z = Randomizer->RAND_INTERVAL(-1.0f, 1.0f);
         p.Normalize(Randomizer->RAND_INTERVAL(0.0f, v - r));
         chemistry->createThermal(r, p, t);
      }
      Randomizer->RAND_POP();
   }


   // Mate members.
   void mate(Member *member1, Member *member2)
   {
      RANDOM s;
      int    i, n;
      Vector p;
      float  r, t, v;

      assert(member1->chemistry != NULL && member2->chemistry != NULL);
      fitness = 0.0f;
      if (member1->generation > member2->generation)
      {
         generation = member1->generation + 1;
      }
      else
      {
         generation = member2->generation + 1;
      }
      s = Randomizer->RAND();
      Randomizer->RAND_PUSH();
      Randomizer->SRAND(s);
      if (Randomizer->RAND_BOOL())
      {
         v = member1->chemistry->vesselRadius;
      }
      else
      {
         v = member2->chemistry->vesselRadius;
      }
      if (Randomizer->RAND_BOOL())
      {
         n = (int)member1->chemistry->atoms.size();
      }
      else
      {
         n = (int)member2->chemistry->atoms.size();
      }
      if (chemistry != NULL) { delete chemistry; }
#ifdef THREADS
      chemistry = new Chemistry(v, s, NumThreads);
#else
      chemistry = new Chemistry(v, s);
#endif
      assert(chemistry != NULL);
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MIN_ATOM_INITIAL_FORCE =
            member1->chemistry->parameters->MIN_ATOM_INITIAL_FORCE;
      }
      else
      {
         chemistry->parameters->MIN_ATOM_INITIAL_FORCE =
            member2->chemistry->parameters->MIN_ATOM_INITIAL_FORCE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MAX_ATOM_INITIAL_FORCE =
            member1->chemistry->parameters->MAX_ATOM_INITIAL_FORCE;
      }
      else
      {
         chemistry->parameters->MAX_ATOM_INITIAL_FORCE =
            member2->chemistry->parameters->MAX_ATOM_INITIAL_FORCE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MIN_NUCLEUS_PROTONS =
            member1->chemistry->parameters->MIN_NUCLEUS_PROTONS;
      }
      else
      {
         chemistry->parameters->MIN_NUCLEUS_PROTONS =
            member2->chemistry->parameters->MIN_NUCLEUS_PROTONS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MAX_NUCLEUS_PROTONS =
            member1->chemistry->parameters->MAX_NUCLEUS_PROTONS;
      }
      else
      {
         chemistry->parameters->MAX_NUCLEUS_PROTONS =
            member2->chemistry->parameters->MAX_NUCLEUS_PROTONS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->PROTON_MASS =
            member1->chemistry->parameters->PROTON_MASS;
      }
      else
      {
         chemistry->parameters->PROTON_MASS =
            member2->chemistry->parameters->PROTON_MASS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->ELECTRON_MASS =
            member1->chemistry->parameters->ELECTRON_MASS;
      }
      else
      {
         chemistry->parameters->ELECTRON_MASS =
            member2->chemistry->parameters->ELECTRON_MASS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->PROTON_CHARGE =
            member1->chemistry->parameters->PROTON_CHARGE;
      }
      else
      {
         chemistry->parameters->PROTON_CHARGE =
            member2->chemistry->parameters->PROTON_CHARGE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->ELECTRON_CHARGE =
            member1->chemistry->parameters->ELECTRON_CHARGE;
      }
      else
      {
         chemistry->parameters->ELECTRON_CHARGE =
            member2->chemistry->parameters->ELECTRON_CHARGE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->CHARGE_GAUSSIAN_SPREAD =
            member1->chemistry->parameters->CHARGE_GAUSSIAN_SPREAD;
      }
      else
      {
         chemistry->parameters->CHARGE_GAUSSIAN_SPREAD =
            member2->chemistry->parameters->CHARGE_GAUSSIAN_SPREAD;
      }
      if (Randomizer->RAND_CHANCE(MUTATION_RATE))
      {
         chemistry->parameters->NUCLEUS_BODY_RADIUS =
            member1->chemistry->parameters->NUCLEUS_BODY_RADIUS;
      }
      else
      {
         chemistry->parameters->NUCLEUS_BODY_RADIUS =
            member2->chemistry->parameters->NUCLEUS_BODY_RADIUS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->ORBITAL_BODY_RADIUS =
            member1->chemistry->parameters->ORBITAL_BODY_RADIUS;
      }
      else
      {
         chemistry->parameters->ORBITAL_BODY_RADIUS =
            member2->chemistry->parameters->ORBITAL_BODY_RADIUS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MAX_BODY_RANGE =
            member1->chemistry->parameters->MAX_BODY_RANGE;
      }
      else
      {
         chemistry->parameters->MAX_BODY_RANGE =
            member2->chemistry->parameters->MAX_BODY_RANGE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->BOND_LENGTH =
            member1->chemistry->parameters->BOND_LENGTH;
      }
      else
      {
         chemistry->parameters->BOND_LENGTH =
            member2->chemistry->parameters->BOND_LENGTH;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->BOND_STIFFNESS =
            member1->chemistry->parameters->BOND_STIFFNESS;
      }
      else
      {
         chemistry->parameters->BOND_STIFFNESS =
            member2->chemistry->parameters->BOND_STIFFNESS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->BOND_DAMPER =
            member1->chemistry->parameters->BOND_DAMPER;
      }
      else
      {
         chemistry->parameters->BOND_DAMPER =
            member2->chemistry->parameters->BOND_DAMPER;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->NUCLEAR_REPULSION_STIFFNESS =
            member1->chemistry->parameters->NUCLEAR_REPULSION_STIFFNESS;
      }
      else
      {
         chemistry->parameters->NUCLEAR_REPULSION_STIFFNESS =
            member2->chemistry->parameters->NUCLEAR_REPULSION_STIFFNESS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->COVALENT_BONDING_RANGE =
            member1->chemistry->parameters->COVALENT_BONDING_RANGE;
      }
      else
      {
         chemistry->parameters->COVALENT_BONDING_RANGE =
            member2->chemistry->parameters->COVALENT_BONDING_RANGE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MIN_COVALENT_BOND_FORCE =
            member1->chemistry->parameters->MIN_COVALENT_BOND_FORCE;
      }
      else
      {
         chemistry->parameters->MIN_COVALENT_BOND_FORCE =
            member2->chemistry->parameters->MIN_COVALENT_BOND_FORCE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->COVALENT_BOND_STIFFNESS_SCALE =
            member1->chemistry->parameters->COVALENT_BOND_STIFFNESS_SCALE;
      }
      else
      {
         chemistry->parameters->COVALENT_BOND_STIFFNESS_SCALE =
            member2->chemistry->parameters->COVALENT_BOND_STIFFNESS_SCALE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MIN_THERMAL_RADIUS =
            member1->chemistry->parameters->MIN_THERMAL_RADIUS;
      }
      else
      {
         chemistry->parameters->MIN_THERMAL_RADIUS =
            member2->chemistry->parameters->MIN_THERMAL_RADIUS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MAX_THERMAL_RADIUS =
            member1->chemistry->parameters->MAX_THERMAL_RADIUS;
      }
      else
      {
         chemistry->parameters->MAX_THERMAL_RADIUS =
            member2->chemistry->parameters->MAX_THERMAL_RADIUS;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MIN_THERMAL_TEMPERATURE =
            member1->chemistry->parameters->MIN_THERMAL_TEMPERATURE;
      }
      else
      {
         chemistry->parameters->MIN_THERMAL_TEMPERATURE =
            member2->chemistry->parameters->MIN_THERMAL_TEMPERATURE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->MAX_TEMPERATURE =
            member1->chemistry->parameters->MAX_TEMPERATURE;
      }
      else
      {
         chemistry->parameters->MAX_TEMPERATURE =
            member2->chemistry->parameters->MAX_TEMPERATURE;
      }
      if (Randomizer->RAND_BOOL())
      {
         chemistry->parameters->UPDATE_STEP =
            member1->chemistry->parameters->UPDATE_STEP;
      }
      else
      {
         chemistry->parameters->UPDATE_STEP =
            member2->chemistry->parameters->UPDATE_STEP;
      }
      chemistry->init(n);
      if (Randomizer->RAND_BOOL())
      {
         n = (int)member1->chemistry->thermals.size();
      }
      else
      {
         n = (int)member2->chemistry->thermals.size();
      }
      for (i = 0; i < n; i++)
      {
         r = Randomizer->RAND_INTERVAL(chemistry->parameters->MIN_THERMAL_RADIUS,
                                       chemistry->parameters->MAX_THERMAL_RADIUS);
         t = Randomizer->RAND_INTERVAL(chemistry->parameters->MIN_THERMAL_TEMPERATURE,
                                       chemistry->parameters->MAX_TEMPERATURE);
         p.x = Randomizer->RAND_INTERVAL(-1.0f, 1.0f);
         p.y = Randomizer->RAND_INTERVAL(-1.0f, 1.0f);
         p.z = Randomizer->RAND_INTERVAL(-1.0f, 1.0f);
         p.Normalize(Randomizer->RAND_INTERVAL(0.0f, v - r));
         chemistry->createThermal(r, p, t);
      }
      Randomizer->RAND_POP();
   }


   // Load.
   void load(FILE *fp)
   {
      if (chemistry == NULL)
      {
#ifdef THREADS
         chemistry = new Chemistry(DEFAULT_VESSEL_RADIUS, (RANDOM)0, NumThreads);
#else
         chemistry = new Chemistry(DEFAULT_VESSEL_RADIUS, (RANDOM)0);
#endif
         assert(chemistry != NULL);
      }
      chemistry->load(fp);
      FREAD_FLOAT(&fitness, fp);
      FREAD_INT(&generation, fp);
   }


   // Save.
   void save(FILE *fp)
   {
      assert(chemistry != NULL);
      chemistry->save(fp);
      FWRITE_FLOAT(&fitness, fp);
      FWRITE_INT(&generation, fp);
   }


   // Print.
   void print()
   {
      sprintf(Log::messageBuf, "Fitness = %f, Generation = %d", fitness, generation);
      Log::logInformation();
      sprintf(Log::messageBuf, "RANDOM_SEED = %d", chemistry->randomSeed);
      Log::logInformation();
      sprintf(Log::messageBuf, "VESSEL_RADIUS = %f", chemistry->vesselRadius);
      Log::logInformation();
      sprintf(Log::messageBuf, "NUM_ATOMS = %d", chemistry->atoms.size());
      Log::logInformation();
      sprintf(Log::messageBuf, "MIN_ATOM_INITIAL_FORCE = %f", chemistry->parameters->MIN_ATOM_INITIAL_FORCE);
      Log::logInformation();
      sprintf(Log::messageBuf, "MAX_ATOM_INITIAL_FORCE = %f", chemistry->parameters->MAX_ATOM_INITIAL_FORCE);
      Log::logInformation();
      sprintf(Log::messageBuf, "MIN_NUCLEUS_PROTONS = %d", chemistry->parameters->MIN_NUCLEUS_PROTONS);
      Log::logInformation();
      sprintf(Log::messageBuf, "MAX_NUCLEUS_PROTONS = %d", chemistry->parameters->MAX_NUCLEUS_PROTONS);
      Log::logInformation();
      sprintf(Log::messageBuf, "PROTON_MASS = %f", chemistry->parameters->PROTON_MASS);
      Log::logInformation();
      sprintf(Log::messageBuf, "ELECTRON_MASS = %f", chemistry->parameters->ELECTRON_MASS);
      Log::logInformation();
      sprintf(Log::messageBuf, "PROTON_CHARGE = %f", chemistry->parameters->PROTON_CHARGE);
      Log::logInformation();
      sprintf(Log::messageBuf, "ELECTRON_CHARGE = %f", chemistry->parameters->ELECTRON_CHARGE);
      Log::logInformation();
      sprintf(Log::messageBuf, "CHARGE_GAUSSIAN_SPREAD = %f", chemistry->parameters->CHARGE_GAUSSIAN_SPREAD);
      Log::logInformation();
      sprintf(Log::messageBuf, "NUCLEUS_BODY_RADIUS = %f", chemistry->parameters->NUCLEUS_BODY_RADIUS);
      Log::logInformation();
      sprintf(Log::messageBuf, "ORBITAL_BODY_RADIUS = %f", chemistry->parameters->ORBITAL_BODY_RADIUS);
      Log::logInformation();
      sprintf(Log::messageBuf, "MAX_BODY_RANGE = %f", chemistry->parameters->MAX_BODY_RANGE);
      Log::logInformation();
      sprintf(Log::messageBuf, "BOND_LENGTH = %f", chemistry->parameters->BOND_LENGTH);
      Log::logInformation();
      sprintf(Log::messageBuf, "BOND_STIFFNESS = %f", chemistry->parameters->BOND_STIFFNESS);
      Log::logInformation();
      sprintf(Log::messageBuf, "BOND_DAMPER = %f", chemistry->parameters->BOND_DAMPER);
      Log::logInformation();
      sprintf(Log::messageBuf, "NUCLEAR_REPULSION_STIFFNESS = %f", chemistry->parameters->NUCLEAR_REPULSION_STIFFNESS);
      Log::logInformation();
      sprintf(Log::messageBuf, "COVALENT_BONDING_RANGE = %f", chemistry->parameters->COVALENT_BONDING_RANGE);
      Log::logInformation();
      sprintf(Log::messageBuf, "MIN_COVALENT_BOND_FORCE = %f", chemistry->parameters->MIN_COVALENT_BOND_FORCE);
      Log::logInformation();
      sprintf(Log::messageBuf, "COVALENT_BOND_STIFFNESS_SCALE = %f", chemistry->parameters->COVALENT_BOND_STIFFNESS_SCALE);
      Log::logInformation();
      sprintf(Log::messageBuf, "NUM_THERMALS = %d", chemistry->thermals.size());
      Log::logInformation();
      sprintf(Log::messageBuf, "MIN_THERMAL_RADIUS = %f", chemistry->parameters->MIN_THERMAL_RADIUS);
      Log::logInformation();
      sprintf(Log::messageBuf, "MAX_THERMAL_RADIUS = %f", chemistry->parameters->MAX_THERMAL_RADIUS);
      Log::logInformation();
      sprintf(Log::messageBuf, "MIN_THERMAL_TEMPERATURE = %f", chemistry->parameters->MIN_THERMAL_TEMPERATURE);
      Log::logInformation();
      sprintf(Log::messageBuf, "MAX_TEMPERATURE = %f", chemistry->parameters->MAX_TEMPERATURE);
      Log::logInformation();
      sprintf(Log::messageBuf, "UPDATE_STEP = %f", chemistry->parameters->UPDATE_STEP);
      Log::logInformation();
   }
};

// Population.
Member *Population[POPULATION_SIZE];

// Start/end functions.
void logParameters();
void loadPopulation(FILE *fp);
void savePopulation(FILE *fp);

// Evolve functions.
void evolve(), evaluate(), prune(), mutate(), mate();

int
main(int argc, char *argv[])
{
   int  i;
   FILE *fp;
   bool gotlogopt, gotprefix;

   // Logging.
   Log::LOGGING_FLAG = EVOLVE_LOGGING;
   Log::setLogFileName((char *)DEFAULT_EVOLVE_LOG_FILE_NAME);

   // Parse arguments.
   Generations   = -1;
   InputFileName = OutputFileName = NULL;
   RandomSeed    = INVALID_RANDOM;
   gotlogopt     = gotprefix = false;
   Unpack        = false;
   UnpackPrefix  = (char *)"affinity";

   for (i = 1; i < argc; i++)
   {
      if (strcmp(argv[i], "-generations") == 0)
      {
         i++;
         if (i >= argc)
         {
            printUsage();
            exit(1);
         }
         Generations = atoi(argv[i]);
         if (Generations < 0)
         {
            printUsage();
            exit(1);
         }
         continue;
      }

      if (strcmp(argv[i], "-cycles") == 0)
      {
         i++;
         if (i >= argc)
         {
            printUsage();
            exit(1);
         }
         if ((Cycles = atoi(argv[i])) < 0)
         {
            printUsage();
            exit(1);
         }
         continue;
      }

#ifdef THREADS
      if (strcmp(argv[i], "-numThreads") == 0)
      {
         i++;
         if (i >= argc)
         {
            printUsage();
            exit(1);
         }
         if ((NumThreads = atoi(argv[i])) < 1)
         {
            printUsage();
            exit(1);
         }
         continue;
      }
#endif

      if (strcmp(argv[i], "-input") == 0)
      {
         i++;
         if (i >= argc)
         {
            printUsage();
            exit(1);
         }
         InputFileName = argv[i];
         continue;
      }

      if (strcmp(argv[i], "-output") == 0)
      {
         i++;
         if (i >= argc)
         {
            printUsage();
            exit(1);
         }
         OutputFileName = argv[i];
         continue;
      }

      if (strcmp(argv[i], "-randomSeed") == 0)
      {
         i++;
         if (i >= argc)
         {
            printUsage();
            exit(1);
         }
         RandomSeed = (RANDOM)atoi(argv[i]);
         continue;
      }

      if (strcmp(argv[i], "-logfile") == 0)
      {
         i++;
         if (i >= argc)
         {
            printUsage();
            exit(1);
         }
         Log::setLogFileName(argv[i]);
         gotlogopt = true;
         continue;
      }

      if (strcmp(argv[i], "-unpack") == 0)
      {
         Unpack = true;
         continue;
      }

      if (strcmp(argv[i], "-prefix") == 0)
      {
         i++;
         if (i >= argc)
         {
            printUsage();
            exit(1);
         }
         UnpackPrefix = argv[i];
         gotprefix    = true;
         continue;
      }

      printUsage();
      exit(1);
   }

   if (Unpack)
   {
      if ((InputFileName == NULL) || (OutputFileName != NULL) ||
          (Cycles >= 0) || (RandomSeed != INVALID_RANDOM) || gotlogopt)
      {
         printUsage();
         exit(1);
      }
   }
   else
   {
      if ((OutputFileName == NULL) || (Cycles < 0) || gotprefix)
      {
         printUsage();
         exit(1);
      }
   }

   // Seed random numbers.
   if (InputFileName != NULL)
   {
      if (RandomSeed != INVALID_RANDOM)
      {
         printUsage();
         exit(1);
      }
      if ((fp = fopen(InputFileName, "r")) == NULL)
      {
         sprintf(Log::messageBuf, "Cannot load population file %s", InputFileName);
         Log::logError();
         exit(1);
      }
      FREAD_LONG(&RandomSeed, fp);
      Randomizer = new Random(RandomSeed);
      assert(Randomizer != NULL);
      Randomizer->RAND_LOAD(fp);
   }
   else
   {
      fp = NULL;
      if (RandomSeed == INVALID_RANDOM)
      {
         RandomSeed = (RANDOM)time(NULL);
      }
      Randomizer = new Random(RandomSeed);
      assert(Randomizer != NULL);
      Randomizer->SRAND(RandomSeed);
   }

   // Log run parameters.
   if (!Unpack)
   {
      Log::logInformation((char *)"Initializing evolve:");
      sprintf(Log::messageBuf, "generations = %d", Generations);
      Log::logInformation();
      if (InputFileName != NULL)
      {
         sprintf(Log::messageBuf, "input = %s", InputFileName);
         Log::logInformation();
      }
      sprintf(Log::messageBuf, "output = %s", OutputFileName);
      Log::logInformation();
      logParameters();
   }

   // Create population.
   if (fp == NULL)
   {
      for (i = 0; i < POPULATION_SIZE; i++)
      {
         Population[i] = new Member(0);
         assert(Population[i] != NULL);
      }
   }
   else
   {
      // Continue run.
      loadPopulation(fp);
      fclose(fp);
   }

   // Unpacking chemistries?
   if (Unpack)
   {
      char  fileName[200];
      int   w;
      float v;
      bool  b;
      for (int i = 0; i < POPULATION_SIZE; i++)
      {
         sprintf(fileName, "%s%d.chem", UnpackPrefix, i);
         if ((fp = fopen(fileName, "w")) == NULL)
         {
            fprintf(stderr, "Cannot save to file %s\n", fileName);
            exit(1);
         }
         FWRITE_LONG(&RandomSeed, fp);
         Randomizer->RAND_SAVE(fp);
         w = WINDOW_WIDTH;
         FWRITE_INT(&w, fp);
         w = WINDOW_HEIGHT;
         FWRITE_INT(&w, fp);
         v = 0.0f;
         FWRITE_FLOAT(&v, fp);
         FWRITE_FLOAT(&v, fp);
         FWRITE_FLOAT(&v, fp);
         w = 0;
         FWRITE_INT(&w, fp);
         b = false;
         FWRITE_BOOL(&b, fp);
         FWRITE_BOOL(&b, fp);
         w = -1;
         FWRITE_INT(&w, fp);
         FWRITE_INT(&w, fp);
         Population[i]->chemistry->save(fp);
         fclose(fp);
      }
   }
   else
   {
      // Evolution loop.
      Log::logInformation((char *)"Begin evolve:");
      for (i = 0; i < Generations; i++)
      {
         sprintf(Log::messageBuf, "Generation = %d", i);
         Log::logInformation();
         evolve();

         // Save population?
         if ((i % SAVE_FREQUENCY) == 0)
         {
            if ((fp = fopen(OutputFileName, "w")) == NULL)
            {
               sprintf(Log::messageBuf, "Cannot save to population file %s", OutputFileName);
               Log::logError();
               exit(1);
            }
            FWRITE_LONG(&RandomSeed, fp);
            Randomizer->RAND_SAVE(fp);
            savePopulation(fp);
            fclose(fp);
         }
      }

      // Save population.
      if ((fp = fopen(OutputFileName, "w")) == NULL)
      {
         sprintf(Log::messageBuf, "Cannot save to population file %s", OutputFileName);
         Log::logError();
         exit(1);
      }
      FWRITE_LONG(&RandomSeed, fp);
      Randomizer->RAND_SAVE(fp);
      savePopulation(fp);
      fclose(fp);
   }

   // Release memory.
   for (i = 0; i < POPULATION_SIZE; i++)
   {
      if (Population[i] != NULL)
      {
         delete Population[i];
         Population[i] = NULL;
      }
   }

   // Close log.
   if (!Unpack)
   {
      Log::logInformation((char *)"End evolve");
      Log::close();
   }

   return(0);
}


// Log run parameters.
void logParameters()
{
   Log::logInformation((char *)"Evolve Parameters:");
   sprintf(Log::messageBuf, "RandomSeed = %d", RandomSeed);
   Log::logInformation();
   sprintf(Log::messageBuf, "FIT_POPULATION_SIZE = %d", FIT_POPULATION_SIZE);
   Log::logInformation();
   sprintf(Log::messageBuf, "NUM_MUTANTS = %d", NUM_MUTANTS);
   Log::logInformation();
   sprintf(Log::messageBuf, "NUM_OFFSPRING = %d", NUM_OFFSPRING);
   Log::logInformation();
   sprintf(Log::messageBuf, "CYCLES = %d", Cycles);
   Log::logInformation();
   Log::logInformation((char *)"Evolvable Parameters:");
   sprintf(Log::messageBuf, "MIN_VESSEL_RADIUS = %f", MIN_VESSEL_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_VESSEL_RADIUS = %f", MAX_VESSEL_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_ATOMS = %d", MIN_ATOMS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_ATOMS = %d", MAX_ATOMS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MIN_ATOM_INITIAL_FORCE = %f", MIN_MIN_ATOM_INITIAL_FORCE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MIN_ATOM_INITIAL_FORCE = %f", MAX_MIN_ATOM_INITIAL_FORCE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MAX_ATOM_INITIAL_FORCE = %f", MIN_MAX_ATOM_INITIAL_FORCE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MAX_ATOM_INITIAL_FORCE = %f", MAX_MAX_ATOM_INITIAL_FORCE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MIN_NUCLEUS_PROTONS = %d", MIN_MIN_NUCLEUS_PROTONS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MIN_NUCLEUS_PROTONS = %d", MAX_MIN_NUCLEUS_PROTONS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MAX_NUCLEUS_PROTONS = %d", MIN_MAX_NUCLEUS_PROTONS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MAX_NUCLEUS_PROTONS = %d", MAX_MAX_NUCLEUS_PROTONS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_PROTON_MASS = %f", MIN_PROTON_MASS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_PROTON_MASS = %f", MAX_PROTON_MASS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_ELECTRON_MASS = %f", MIN_ELECTRON_MASS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_ELECTRON_MASS = %f", MAX_ELECTRON_MASS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_PROTON_CHARGE = %f", MIN_PROTON_CHARGE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_PROTON_CHARGE = %f", MAX_PROTON_CHARGE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_ELECTRON_CHARGE = %f", MIN_ELECTRON_CHARGE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_ELECTRON_CHARGE = %f", MAX_ELECTRON_CHARGE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_CHARGE_GAUSSIAN_SPREAD = %f", MIN_CHARGE_GAUSSIAN_SPREAD);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_CHARGE_GAUSSIAN_SPREAD = %f", MAX_CHARGE_GAUSSIAN_SPREAD);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_NUCLEUS_BODY_RADIUS = %f", MIN_NUCLEUS_BODY_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_NUCLEUS_BODY_RADIUS = %f", MAX_NUCLEUS_BODY_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_ORBITAL_BODY_RADIUS = %f", MIN_ORBITAL_BODY_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_ORBITAL_BODY_RADIUS = %f", MAX_ORBITAL_BODY_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MAX_BODY_RANGE = %f", MIN_MAX_BODY_RANGE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MAX_BODY_RANGE = %f", MAX_MAX_BODY_RANGE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_BOND_LENGTH = %f", MIN_BOND_LENGTH);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_BOND_LENGTH = %f", MAX_BOND_LENGTH);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_BOND_STIFFNESS = %f", MIN_BOND_STIFFNESS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_BOND_STIFFNESS = %f", MAX_BOND_STIFFNESS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_BOND_DAMPER = %f", MIN_BOND_DAMPER);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_BOND_DAMPER = %f", MAX_BOND_DAMPER);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_COVALENT_BONDING_RANGE = %f", MIN_COVALENT_BONDING_RANGE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_COVALENT_BONDING_RANGE = %f", MAX_COVALENT_BONDING_RANGE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MIN_COVALENT_BOND_FORCE = %f", MIN_MIN_COVALENT_BOND_FORCE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MIN_COVALENT_BOND_FORCE = %f", MAX_MIN_COVALENT_BOND_FORCE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_COVALENT_BOND_STIFFNESS_SCALE = %f", MIN_COVALENT_BOND_STIFFNESS_SCALE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_COVALENT_BOND_STIFFNESS_SCALE = %f", MAX_COVALENT_BOND_STIFFNESS_SCALE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_THERMALS = %d", MIN_THERMALS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_THERMALS = %d", MAX_THERMALS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MIN_THERMAL_RADIUS = %f", MIN_MIN_THERMAL_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MIN_THERMAL_RADIUS = %f", MAX_MIN_THERMAL_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MAX_THERMAL_RADIUS = %f", MIN_MAX_THERMAL_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MAX_THERMAL_RADIUS = %f", MAX_MAX_THERMAL_RADIUS);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MIN_THERMAL_TEMPERATURE = %f", MIN_MIN_THERMAL_TEMPERATURE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MIN_THERMAL_TEMPERATURE = %f", MAX_MIN_THERMAL_TEMPERATURE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_MAX_TEMPERATURE = %f", MIN_MAX_TEMPERATURE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_MAX_TEMPERATURE = %f", MAX_MAX_TEMPERATURE);
   Log::logInformation();
   sprintf(Log::messageBuf, "MIN_UPDATE_STEP = %f", MIN_UPDATE_STEP);
   Log::logInformation();
   sprintf(Log::messageBuf, "MAX_UPDATE_STEP = %f", MAX_UPDATE_STEP);
   Log::logInformation();
}


// Load evolution population.
void loadPopulation(FILE *fp)
{
   for (int i = 0; i < POPULATION_SIZE; i++)
   {
      Population[i] = new Member();
      assert(Population[i] != NULL);
      Population[i]->load(fp);
   }
}


// Save evolution population.
void savePopulation(FILE *fp)
{
   for (int i = 0; i < POPULATION_SIZE; i++)
   {
      Population[i]->save(fp);
   }
}


// Evolution generation.
void evolve()
{
   // Evaluate member fitness.
   evaluate();

   // Prune unfit members.
   prune();

   // Create new members by mutation.
   mutate();

   // Create new members by mating.
   mate();
}


// Evaluate member fitnesses.
void evaluate()
{
   Log::logInformation((char *)"Evaluate:");

   for (int i = 0; i < POPULATION_SIZE; i++)
   {
      Population[i]->evaluate();
      sprintf(Log::messageBuf, "  Member = %d, Fitness = %f, Generation = %d",
              i, Population[i]->fitness, Population[i]->generation);
      Log::logInformation();
   }
}


// Prune unfit members.
void prune()
{
   double max;
   int    i, j, m;
   Member *member;
   Member *fitPopulation[FIT_POPULATION_SIZE];

   Log::logInformation((char *)"Select:");
   for (i = 0; i < FIT_POPULATION_SIZE; i++)
   {
      m = -1;
      for (j = 0; j < POPULATION_SIZE; j++)
      {
         member = Population[j];
         if (member == NULL)
         {
            continue;
         }
         if ((m == -1) || (member->fitness > max))
         {
            m   = j;
            max = member->fitness;
         }
      }
      member           = Population[m];
      Population[m]    = NULL;
      fitPopulation[i] = member;
      sprintf(Log::messageBuf, "  Fitness = %f, Generation = %d",
              member->fitness, member->generation);
      Log::logInformation();
   }
   for (i = 0; i < POPULATION_SIZE; i++)
   {
      if (Population[i] != NULL)
      {
         delete Population[i];
         Population[i] = NULL;
      }
   }
   for (i = 0; i < FIT_POPULATION_SIZE; i++)
   {
      Population[i] = fitPopulation[i];
   }
}


// Mutate members.
void mutate()
{
   int    i, j;
   Member *member, *mutant;

   Log::logInformation((char *)"Mutate:");
   for (i = 0; i < NUM_MUTANTS; i++)
   {
      // Select a fit member to mutate.
      j      = Randomizer->RAND_CHOICE(FIT_POPULATION_SIZE);
      member = Population[j];

      // Create mutant member.
      mutant = new Member();
      assert(mutant != NULL);
      mutant->mutate(member);
      Population[FIT_POPULATION_SIZE + i] = mutant;
      sprintf(Log::messageBuf, "  Member = %d -> Member = %d", j, FIT_POPULATION_SIZE + i);
      Log::logInformation();
   }
}


// Produce offspring by mating parents.
void mate()
{
   int    i, j, k;
   Member *member1, *member2, *offspring;

   Log::logInformation((char *)"Mate:");
   if (FIT_POPULATION_SIZE < 2)
   {
      return;
   }
   for (i = 0; i < NUM_OFFSPRING; i++)
   {
      // Select a pair of fit members to mate.
      j       = Randomizer->RAND_CHOICE(FIT_POPULATION_SIZE);
      member1 = Population[j];
      while ((k = Randomizer->RAND_CHOICE(FIT_POPULATION_SIZE)) == j)
      {
      }
      member2 = Population[k];

      // Create offspring.
      offspring = new Member();
      assert(offspring != NULL);
      offspring->mate(member1, member2);
      Population[FIT_POPULATION_SIZE + NUM_MUTANTS + i] = offspring;
      sprintf(Log::messageBuf, "  Members = %d, %d -> Member = %d",
              j, k, FIT_POPULATION_SIZE + NUM_MUTANTS + i);
      Log::logInformation();
   }
}
