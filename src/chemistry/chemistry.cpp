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

#include "chemistry.hpp"
using namespace affinity;

// Constructor.
#ifdef THREADS
Chemistry::Chemistry(float vesselRadius, RANDOM randomSeed, int numThreads)
#else
Chemistry::Chemistry(float vesselRadius, RANDOM randomSeed)
#endif
{
   parameters = new Parameters();
   assert(parameters != NULL);
   this->vesselRadius = vesselRadius;
   this->randomSeed   = randomSeed;
   randomizer         = NULL;
   bodyTracker        = NULL;
   bondUpdate         = false;

#ifdef THREADS
   // Start additional chemistry update threads.
   assert(numThreads > 0);
   terminate        = false;
   this->numThreads = numThreads;
   if (numThreads > 1)
   {
      if (pthread_barrier_init(&updateBarrier, NULL, numThreads) != 0)
      {
         fprintf(stderr, "pthread_barrier_init failed, errno=%d\n", errno);
         exit(1);
      }
      if (pthread_mutex_init(&updateMutex, NULL) != 0)
      {
         fprintf(stderr, "pthread_mutex_init failed, errno=%d\n", errno);
         exit(1);
      }
      threads = new pthread_t[numThreads - 1];
      assert(threads != NULL);
      struct ThreadInfo *info;
      for (int i = 0; i < numThreads - 1; i++)
      {
         info = new struct ThreadInfo;
         assert(info != NULL);
         info->chemistry = this;
         info->threadNum = i + 1;
         if (pthread_create(&threads[i], NULL, updateThread, (void *)info) != 0)
         {
            fprintf(stderr, "pthread_create failed, errno=%d\n", errno);
            exit(1);
         }
      }
   }
#endif
}


// Initialize chemistry.
void Chemistry::init(int numAtoms)
{
   int i, n;

   clear();
   randomizer = new Random(randomSeed);
   assert(randomizer != NULL);
   bodyTracker = new Octree(0.0f, 0.0f, 0.0f,
                            vesselRadius * 1.5f, parameters->BOND_LENGTH);
   assert(bodyTracker != NULL);
   for (i = 0; i < numAtoms; i++)
   {
      n = -1;
#if (H2O_MOLECULES || O2_MOLECULES)
      if (randomizer->RAND_BOOL())
      {
         n = 6;
      }
      else
      {
         n = 1;
      }
#if (CO2_MOLECULES)
      if (randomizer->RAND_CHOICE(4) == 0)
      {
         n = 4;
      }
#endif
#else
#if (CO2_MOLECULES)
      if (randomizer->RAND_CHOICE(3) == 0)
      {
         n = 4;
      }
      else
      {
         n = 6;
      }
#endif
#endif
#if (ORGANIC_MOLECULES)
      switch (randomizer->RAND_CHOICE(6))
      {
      case 0:
      case 1:
         n = 1;
         break;

      case 2:
      case 3:
         n = 4;
         break;

      case 4:
         n = 5;
         break;

      case 5:
         n = 6;
         break;
      }
#endif
      if (n == -1)
      {
         n = randomizer->RAND_CHOICE(
            parameters->MAX_NUCLEUS_PROTONS -
            parameters->MIN_NUCLEUS_PROTONS + 1) +
             parameters->MIN_NUCLEUS_PROTONS;
      }
      createAtom(n);
   }
}


// Destructor.
Chemistry::~Chemistry()
{
#ifdef THREADS
   // Terminate threads.
   if (numThreads > 1)
   {
      // Unblock threads waiting on update barrier.
      terminate = true;
      update();
      for (int i = 0; i < numThreads - 1; i++)
      {
         pthread_join(threads[i], NULL);
         pthread_detach(threads[i]);
      }
      pthread_mutex_destroy(&updateMutex);
      pthread_barrier_destroy(&updateBarrier);
      delete threads;
   }
#endif
   clear();
   if (parameters != NULL)
   {
      delete parameters;
   }
}


// Clear chemistry.
void Chemistry::clear()
{
   int i, j;

   atomIDfactory = 0;
   for (i = 0, j = (int)molecules.size(); i < j; i++)
   {
      if (molecules[i] != NULL)
      {
         delete molecules[i];
      }
   }
   molecules.clear();
   bodies.clear();
   if (bodyTracker != NULL)
   {
      delete bodyTracker;
   }
   bodyTracker = NULL;
   for (i = 0, j = (int)atoms.size(); i < j; i++)
   {
      if (atoms[i] != NULL)
      {
         delete atoms[i];
      }
   }
   atoms.clear();
   for (i = 0, j = (int)thermals.size(); i < j; i++)
   {
      if (thermals[i] != NULL)
      {
         delete thermals[i];
      }
   }
   thermals.clear();
   if (randomizer != NULL)
   {
      delete randomizer;
   }
   randomizer = NULL;
}


// Create atom and add to system.
Atom *Chemistry::createAtom(int protons, int id)
{
   int       s, s2, o, o2;
   Atom      *atom;
   OctObject *b;

   if (bodyTracker == NULL)
   {
      init(0);
   }
   if (id == -1)
   {
      atom = new Atom(parameters, atomIDfactory, protons, randomizer);
      assert(atom != NULL);
      atomIDfactory++;
   }
   else
   {
      atom = new Atom(parameters, id, protons, randomizer);
      assert(atom != NULL);
   }
   atoms.push_back(atom);
   atom->nucleus.position.x = (float)randomizer->RAND_INTERVAL(-1.0f, 1.0f);
   atom->nucleus.position.y = (float)randomizer->RAND_INTERVAL(-1.0f, 1.0f);
   atom->nucleus.position.z = (float)randomizer->RAND_INTERVAL(-1.0f, 1.0f);
   atom->nucleus.position.Normalize((float)randomizer->RAND_INTERVAL(0.0f,
                                                              vesselRadius - ((atom->shells.size() * parameters->BOND_LENGTH) +
                                                                              atom->nucleus.radius)));
   atom->nucleus.forces.x = (float)randomizer->RAND_INTERVAL(-1.0f, 1.0f);
   atom->nucleus.forces.y = (float)randomizer->RAND_INTERVAL(-1.0f, 1.0f);
   atom->nucleus.forces.z = (float)randomizer->RAND_INTERVAL(-1.0f, 1.0f);
   atom->nucleus.forces.Normalize(
	   (float)randomizer->RAND_INTERVAL(parameters->MIN_ATOM_INITIAL_FORCE,
                                parameters->MAX_ATOM_INITIAL_FORCE));
   b = new OctObject(atom->nucleus.position, (void *)&atom->nucleus);
   assert(b != NULL);
   bodies.push_back(b);
   bodyTracker->insert(b);
   for (s = 0, s2 = (int)atom->shells.size(); s < s2; s++)
   {
      for (o = 0, o2 = (int)atom->shells[s].orbitals.size(); o < o2; o++)
      {
         atom->shells[s].orbitals[o].position += atom->nucleus.position;
         b = new OctObject(atom->shells[s].orbitals[o].position,
                           (void *)&atom->shells[s].orbitals[o]);
         assert(b != NULL);
         bodies.push_back(b);
         bodyTracker->insert(b);
      }
   }
   return(atom);
}


// Add atom to system.
int Chemistry::addAtom(Atom *atom)
{
   int       s, s2, o, o2;
   OctObject *b;

   if (bodyTracker == NULL)
   {
      init(0);
   }
   atom->setID(atomIDfactory);
   atomIDfactory++;
   atom->setParameters(parameters);
   atoms.push_back(atom);
   b = new OctObject(atom->nucleus.position, (void *)&atom->nucleus);
   assert(b != NULL);
   bodies.push_back(b);
   bodyTracker->insert(b);
   for (s = 0, s2 = (int)atom->shells.size(); s < s2; s++)
   {
      for (o = 0, o2 = (int)atom->shells[s].orbitals.size(); o < o2; o++)
      {
         b = new OctObject(atom->shells[s].orbitals[o].position,
                           (void *)&atom->shells[s].orbitals[o]);
         assert(b != NULL);
         bodies.push_back(b);
         bodyTracker->insert(b);
      }
   }
   return(atom->getID());
}


// Remove atom from system.
void Chemistry::removeAtom(int id)
{
   int  i, j;
   Body *body;

   vector<Atom *>      tmpAtoms;
   vector<OctObject *> tmpBodies;

   if (bodyTracker == NULL)
   {
      init(0);
   }
   for (i = 0, j = (int)bodies.size(); i < j; i++)
   {
      body = (Body *)bodies[i]->client;
      if (body->id == id)
      {
         bodyTracker->remove(bodies[i]);
      }
      else
      {
         tmpBodies.push_back(bodies[i]);
         if ((body->covalentBody != NULL) && (body->covalentBody->id == id))
         {
            body->covalentBody->covalentBody = NULL;
            body->covalentBody = NULL;
         }
      }
   }
   bodies.clear();
   for (i = 0, j = (int)tmpBodies.size(); i < j; i++)
   {
      bodies.push_back(tmpBodies[i]);
   }
   for (i = 0, j = (int)atoms.size(); i < j; i++)
   {
      if (atoms[i]->getID() != id)
      {
         tmpAtoms.push_back(atoms[i]);
      }
      else
      {
         delete atoms[i];
      }
   }
   atoms.clear();
   for (i = 0, j = (int)tmpAtoms.size(); i < j; i++)
   {
      atoms.push_back(tmpAtoms[i]);
   }
}


// Get atom by ID.
Atom *Chemistry::getAtom(int id)
{
   int i, j;

   if (bodyTracker == NULL)
   {
      init(0);
   }
   for (i = 0, j = (int)atoms.size(); i < j; i++)
   {
      if (atoms[i]->getID() == id)
      {
         return(atoms[i]);
      }
   }
   return(NULL);
}


// Create thermal object and add to system.
Thermal *Chemistry::createThermal(float radius, Vector& position, float temperature)
{
   Thermal *thermal;

   if (bodyTracker == NULL)
   {
      init(0);
   }
   thermal = new Thermal(parameters, radius, position, temperature);
   assert(thermal != NULL);
   thermals.push_back(thermal);
   return(thermal);
}


// Add thermal object to system.
void Chemistry::addThermal(Thermal *thermal)
{
   if (bodyTracker == NULL)
   {
      init(0);
   }
   thermal->parameters = parameters;
   thermals.push_back(thermal);
}


// Update system.
void Chemistry::update()
{
   // Update with base thread.
   update(0);
}


void Chemistry::update(int threadNum)
{
   int    i, i2, j, j2, k;
   float  b, d, s;
   Vector x, f, n, v, p, m;
   Body   *b1, *b2;

   list<OctObject *>           searchList;
   list<OctObject *>::iterator searchItr;
   enum {
      BODY_SHIFT_TRIES = 10
   };
#ifdef THREADS
   vector<Body *> unbond;
   vector<Body *> bond1, bond2;
#endif

   // Ensure initialization.
   if ((threadNum == 0) && (bodyTracker == NULL))
   {
      init(0);
   }
   bondUpdate = false;

#ifdef THREADS
   // Synchronize threads.
   if (numThreads > 1)
   {
      i = pthread_barrier_wait(&updateBarrier);
      if ((i != PTHREAD_BARRIER_SERIAL_THREAD) && (i != 0))
      {
         fprintf(stderr, "pthread_barrier_wait failed, errno=%d\n", errno);
         exit(1);
      }
      if (terminate)
      {
         if (threadNum == 0)
         {
            return;
         }
         pthread_exit(NULL);
      }
   }
#endif

   // Break over-extended bonds.
   for (i = 0, i2 = (int)bodies.size(); i < i2; i++)
   {
#ifdef THREADS
      // Divide work among threads.
      if ((i % numThreads) != threadNum)
      {
         continue;
      }
#endif

      b1 = (Body *)bodies[i]->client;
      if ((b1->covalentBody != NULL) && (b1->id < b1->covalentBody->id))
      {
         x = b1->covalentBody->position - b1->position;
         d = x.Magnitude();
         if (d > parameters->COVALENT_BONDING_RANGE)
         {
#ifdef THREADS
            if (numThreads > 1)
            {
               unbond.push_back(b1);
            }
            else
            {
               b1->covalentBody->covalentBody = NULL;
               b1->covalentBody = NULL;
               bondUpdate       = true;
            }
#else
            b1->covalentBody->covalentBody = NULL;
            b1->covalentBody = NULL;
            bondUpdate       = true;
#endif
         }
      }
   }
#ifdef THREADS
   // Re-group threads before performing actual unbonding.
   if (numThreads > 1)
   {
      pthread_barrier_wait(&updateBarrier);
      for (i = 0, i2 = (int)unbond.size(); i < i2; i++)
      {
         b1 = unbond[i];
         b1->covalentBody->covalentBody = NULL;
         b1->covalentBody = NULL;
         bondUpdate       = true;
      }
   }
#endif

   // Do covalent bonding:
   // A covalent bond is a 0-length spring connecting orbitals.
   // The stiffness of the spring is proportional to the covalent
   // bonding force. A bond forms when valence orbitals draw
   // within a certain distance of each other.
#ifdef THREADS
   if (numThreads > 1)
   {
      pthread_barrier_wait(&updateBarrier);
   }
#endif
   for (i = 0, i2 = (int)bodies.size(); i < i2; i++)
   {
#ifdef THREADS
      // Divide work among threads.
      if ((i % numThreads) != threadNum)
      {
         continue;
      }
#endif

      // Bond to nearby orbitals.
      b1 = (Body *)bodies[i]->client;
      searchList.clear();
      bodyTracker->search(b1->position, parameters->MAX_BODY_RANGE, searchList);
      for (searchItr = searchList.begin();
           searchItr != searchList.end(); searchItr++)
      {
         b2 = (Body *)(*searchItr)->client;
         if (b1->id == b2->id)
         {
            continue;
         }
         if (!b1->hasValence || !b2->hasValence)
         {
            continue;
         }
         x = b2->position - b1->position;
         d = x.Magnitude();
         if (d > parameters->COVALENT_BONDING_RANGE)
         {
            continue;
         }
         b = b1->getCovalentForce(b2);
         if ((b1->covalentBody != NULL) && (b <= b1->getCovalentForce(b1->covalentBody)))
         {
            continue;
         }
         if ((b2->covalentBody != NULL) && (b <= b2->getCovalentForce(b2->covalentBody)))
         {
            continue;
         }
         if (b1->covalentBody != NULL)
         {
            b1->covalentBody->covalentBody = NULL;
            b1->covalentBody = NULL;
         }
         if (b2->covalentBody != NULL)
         {
            b2->covalentBody->covalentBody = NULL;
            b2->covalentBody = NULL;
         }
#ifdef THREADS
         if (numThreads > 1)
         {
            bond1.push_back(b1);
            bond2.push_back(b2);
         }
         else
         {
            b1->covalentBody = b2;
            b2->covalentBody = b1;
            bondUpdate       = true;
         }
#else
         b1->covalentBody = b2;
         b2->covalentBody = b1;
         bondUpdate       = true;
#endif
      }
   }
#ifdef THREADS
   // Re-group threads before performing actual bonding.
   if (numThreads > 1)
   {
      pthread_barrier_wait(&updateBarrier);
      for (i = 0, i2 = (int)bond1.size(); i < i2; i++)
      {
         b1 = bond1[i];
         b2 = bond2[i];

         // Serialize to prevent duplicate bonds.
         pthread_mutex_lock(&updateMutex);
         if ((b1->covalentBody == NULL) && (b2->covalentBody == NULL))
         {
            b1->covalentBody = b2;
            b2->covalentBody = b1;
            bondUpdate       = true;
         }
         pthread_mutex_unlock(&updateMutex);
      }
   }
#endif

   // Do charge forces.
#ifdef THREADS
   if (numThreads > 1)
   {
      pthread_barrier_wait(&updateBarrier);
   }
#endif
   for (i = 0, i2 = (int)bodies.size(); i < i2; i++)
   {
#ifdef THREADS
      if ((i % numThreads) != threadNum)
      {
         continue;
      }
#endif
      b1 = (Body *)bodies[i]->client;
      searchList.clear();
      bodyTracker->search(b1->position, parameters->MAX_BODY_RANGE, searchList);
      for (searchItr = searchList.begin();
           searchItr != searchList.end(); searchItr++)
      {
         b2 = (Body *)(*searchItr)->client;
         if (b1 == b2)
         {
            continue;
         }
         if ((b1->id == b2->id) && (b1->shell != b2->shell))
         {
            continue;
         }
         if (b1->covalentBody == b2)
         {
            continue;
         }

         // Prevent overlapping bodies.
         for (k = 0; k < BODY_SHIFT_TRIES; k++)
         {
            x = b2->position - b1->position;
            d = x.Magnitude();
            if (d > tol)
            {
               break;
            }
            d = parameters->ORBITAL_BODY_RADIUS * 0.1f;
            b1->position.x += (float)randomizer->RAND_INTERVAL(-d, d);
            b1->position.y += (float)randomizer->RAND_INTERVAL(-d, d);
            b1->position.z += (float)randomizer->RAND_INTERVAL(-d, d);
            b2->position.x += (float)randomizer->RAND_INTERVAL(-d, d);
            b2->position.y += (float)randomizer->RAND_INTERVAL(-d, d);
            b2->position.z += (float)randomizer->RAND_INTERVAL(-d, d);
         }
         x.Normalize();

         // Charge force: gaussian with max=charge product.
         f = x * (b1->charge * b2->charge) *
             (float)exp(-(double)((d * d) /
                                  (parameters->CHARGE_GAUSSIAN_SPREAD *
                                   parameters->CHARGE_GAUSSIAN_SPREAD)));
         b1->forces -= f;
         b2->forces += f;
      }
   }

   // Do nuclear repulsion forces:
   // A nucleus repulses "foreign" bodies within its outer shell.
   for (i = 0, i2 = (int)atoms.size(); i < i2; i++)
   {
#ifdef THREADS
      if ((i % numThreads) != threadNum)
      {
         continue;
      }
#endif
      b1 = &atoms[i]->nucleus;
      searchList.clear();
      bodyTracker->search(b1->position, parameters->MAX_BODY_RANGE, searchList);
      for (searchItr = searchList.begin();
           searchItr != searchList.end(); searchItr++)
      {
         b2 = (Body *)(*searchItr)->client;
         if (b1->id == b2->id)
         {
            continue;
         }
         x = b2->position - b1->position;
         d = x.Magnitude() - (parameters->BOND_LENGTH *
                              (float)(atoms[i]->shells.size() + 1));
         if (d > 0.0f)
         {
            continue;
         }
         x.Normalize();
         f = (-parameters->NUCLEAR_REPULSION_STIFFNESS *
              (float)atoms[i]->number * d * x);
         b1->forces -= f;
         b2->forces += f;
      }
   }

   // Update orbital bond forces.
   for (i = 0, i2 = (int)atoms.size(); i < i2; i++)
   {
#ifdef THREADS
      if ((i % numThreads) != threadNum)
      {
         continue;
      }
#endif
      atoms[i]->updateOrbitalBonds();
   }

   // Update covalent bond forces.
   for (i = 0, i2 = (int)bodies.size(); i < i2; i++)
   {
#ifdef THREADS
      if ((i % numThreads) != threadNum)
      {
         continue;
      }
#endif
      b1 = (Body *)bodies[i]->client;
      if ((b1->covalentBody != NULL) && (b1->id < b1->covalentBody->id))
      {
         b1->updateCovalentBond();
      }
   }

   // Update atom velocities and positions.
   for (i = 0, i2 = (int)atoms.size(); i < i2; i++)
   {
#ifdef THREADS
      if ((i % numThreads) != threadNum)
      {
         continue;
      }
#endif
      atoms[i]->update(parameters->UPDATE_STEP);
   }

   // Contain bodies inside vessel.
   for (i = 0, i2 = (int)bodies.size(); i < i2; i++)
   {
#ifdef THREADS
      if ((i % numThreads) != threadNum)
      {
         continue;
      }
#endif
      b1 = (Body *)bodies[i]->client;
      d  = b1->position.Magnitude();
      if (d > vesselRadius * 1.25f)
      {
         b1->position.Normalize(vesselRadius - b1->radius);
      }
      p = b1->position + (b1->velocity * parameters->UPDATE_STEP);
      if ((d >= (vesselRadius - b1->radius)) && (d < p.Magnitude()))
      {
         if (d > tol)
         {
            n = -b1->position;
            n.Normalize();
            v            = b1->velocity;
            b1->velocity = v - (2.0f * (v * n) * n);
         }
         continue;
      }

      // Handle collisions with thermal objects.
      // Bodies tend to achieve a speed equal to the temperature of the thermal.
      for (j = 0, j2 = (int)thermals.size(); j < j2; j++)
      {
         n = b1->position - thermals[j]->position;
         d = n.Magnitude();
         if ((d <= (thermals[j]->radius + b1->radius)) &&
             ((p - thermals[j]->position).Magnitude() < d))
         {
            if (d > tol)
            {
               n.Normalize();
               v            = b1->velocity;
               b1->velocity = v - (2.0f * (v * n) * n);
               s            = b1->velocity.Magnitude();
               s           += (thermals[j]->temperature - s) / 2.0f;
               b1->velocity.Normalize(s);
            }
            break;
         }
      }
   }

#ifdef THREADS
   // Update body tracker.
   if (numThreads > 1)
   {
      if (threadNum == 0)
      {
         moveList.clear();
      }
      pthread_barrier_wait(&updateBarrier);
   }
#endif
   for (i = 0, i2 = (int)bodies.size(); i < i2; i++)
   {
#ifdef THREADS
      if ((i % numThreads) != threadNum)
      {
         continue;
      }
#endif
      b1 = (Body *)bodies[i]->client;
#ifdef THREADS
      if (numThreads > 1)
      {
         bodies[i]->position = b1->position;

         // Save object if move causes a tracker change.
         if (!bodies[i]->isInside(bodies[i]->node))
         {
            pthread_mutex_lock(&updateMutex);
            moveList.push_back(bodies[i]);
            pthread_mutex_unlock(&updateMutex);
         }
      }
      else
      {
         bodies[i]->move(b1->position);
      }
#else
      bodies[i]->move(b1->position);
#endif
   }
#ifdef THREADS
   // Update tracker with moved objects.
   if (numThreads > 1)
   {
      pthread_barrier_wait(&updateBarrier);
      if (threadNum == 0)
      {
         for (i = 0, i2 = (int)moveList.size(); i < i2; i++)
         {
            b1 = (Body *)moveList[i]->client;
            moveList[i]->move(b1->position);
         }
      }
   }
#endif

#ifdef THREADS
   // Re-group threads.
   if (numThreads > 1)
   {
      pthread_barrier_wait(&updateBarrier);
   }
#endif
}


#ifdef THREADS
// Chemistry update thread.
void *Chemistry::updateThread(void *arg)
{
   struct ThreadInfo *info      = (struct ThreadInfo *)arg;
   Chemistry         *chemistry = info->chemistry;
   int               threadNum  = info->threadNum;

   delete info;
   while (true)
   {
      chemistry->update(threadNum);
   }
   return(NULL);
}


#endif

// Clear atom marks.
void Chemistry::clearAtomMarks()
{
   int i, i2;

   for (i = 0, i2 = (int)atoms.size(); i < i2; i++)
   {
      atoms[i]->mark = -1;
   }
}


// Mark molecule.
void Chemistry::markMolecule(Atom *atom, vector<int>& atomCounts, int mark)
{
   int  s, s2, o, o2;
   Body *body;

   if (atom->mark != -1)
   {
      return;
   }
   atom->mark = mark;
   atomCounts[atom->number]++;
   for (s = 0, s2 = (int)atom->shells.size(); s < s2; s++)
   {
      for (o = 0, o2 = (int)atom->shells[s].orbitals.size(); o < o2; o++)
      {
         if ((body = atom->shells[s].orbitals[o].covalentBody) != NULL)
         {
            markMolecule(getAtom(body->id), atomCounts, mark);
         }
      }
   }
}


// Generate molecules.
void Chemistry::generateMolecules()
{
   int i, i2, j, j2, k;

   vector<int> atomCounts;
   Molecule    *molecule;

   for (i = 0, i2 = (int)molecules.size(); i < i2; i++)
   {
      if (molecules[i] != NULL)
      {
         delete molecules[i];
      }
   }
   molecules.clear();
   clearAtomMarks();
   atomCounts.resize(parameters->MAX_NUCLEUS_PROTONS + 1);
   for (i = k = 0, i2 = (int)atoms.size(); i < i2; i++)
   {
      if (atoms[i]->mark == -1)
      {
         for (j = 0, j2 = (int)atomCounts.size(); j < j2; j++)
         {
            atomCounts[j] = 0;
         }
         markMolecule(atoms[i], atomCounts, k);
         k++;
         molecule = new Molecule(this, atoms[i]);
         assert(molecule != NULL);
         molecules.push_back(molecule);
      }
   }
}


// Load chemistry.
void Chemistry::load(FILE *fp)
{
   int       i, j, p, q, id, id2, s, s2, o, o2;
   OctObject *b;
   Atom      *atom;
   Body      *b1, *b2;
   Thermal   *thermal;

   init(0);
   parameters->load(fp);
   FREAD_INT(&atomIDfactory, fp);
   FREAD_FLOAT(&vesselRadius, fp);
   randomizer->RAND_LOAD(fp);
   FREAD_INT(&j, fp);
   for (i = 0; i < j; i++)
   {
      atom = new Atom(parameters);
      assert(atom != NULL);
      atom->load(fp);
      atoms.push_back(atom);
      b = new OctObject(atom->nucleus.position, (void *)&atom->nucleus);
      assert(b != NULL);
      bodies.push_back(b);
      bodyTracker->insert(b);
      for (s = 0, s2 = (int)atom->shells.size(); s < s2; s++)
      {
         for (o = 0, o2 = (int)atom->shells[s].orbitals.size(); o < o2; o++)
         {
            b = new OctObject(atom->shells[s].orbitals[o].position,
                              (void *)&atom->shells[s].orbitals[o]);
            assert(b != NULL);
            bodies.push_back(b);
            bodyTracker->insert(b);
         }
      }
   }
   FREAD_INT(&j, fp);
   for (i = 0; i < j; i++)
   {
      FREAD_INT(&id, fp);
      FREAD_INT(&s, fp);
      FREAD_INT(&o, fp);
      for (p = 0, q = (int)bodies.size(); p < q; p++)
      {
         b1 = (Body *)bodies[p]->client;
         if ((b1->id == id) && (b1->shell == s) && (b1->orbital == o))
         {
            break;
         }
      }
      assert(p < q);
      FREAD_INT(&id2, fp);
      FREAD_INT(&s2, fp);
      FREAD_INT(&o2, fp);
      for (p = 0, q = (int)bodies.size(); p < q; p++)
      {
         b2 = (Body *)bodies[p]->client;
         if ((b2->id == id2) && (b2->shell == s2) && (b2->orbital == o2))
         {
            break;
         }
      }
      assert(p < q);
      b1->covalentBody = b2;
      b2->covalentBody = b1;
   }
   FREAD_INT(&j, fp);
   for (i = 0; i < j; i++)
   {
      thermal = new Thermal(parameters);
      assert(thermal != NULL);
      thermal->load(fp);
      thermals.push_back(thermal);
   }
}


// Save chemistry.
void Chemistry::save(FILE *fp)
{
   int     i, j, k;
   Atom    *atom;
   Body    *body;
   Thermal *thermal;

   if (bodyTracker == NULL)
   {
      init(0);
   }
   parameters->save(fp);
   FWRITE_INT(&atomIDfactory, fp);
   FWRITE_FLOAT(&vesselRadius, fp);
   randomizer->RAND_SAVE(fp);
   j = (int)atoms.size();
   FWRITE_INT(&j, fp);
   for (i = 0; i < j; i++)
   {
      atom = atoms[i];
      atom->save(fp);
   }
   for (i = k = 0, j = (int)bodies.size(); i < j; i++)
   {
      body = (Body *)bodies[i]->client;
      if ((body->covalentBody != NULL) && (body->id < body->covalentBody->id))
      {
         k++;
      }
   }
   FWRITE_INT(&k, fp);
   for (i = 0, j = (int)bodies.size(); i < j; i++)
   {
      body = (Body *)bodies[i]->client;
      if ((body->covalentBody != NULL) && (body->id < body->covalentBody->id))
      {
         FWRITE_INT(&body->id, fp);
         FWRITE_INT(&body->shell, fp);
         FWRITE_INT(&body->orbital, fp);
         FWRITE_INT(&body->covalentBody->id, fp);
         FWRITE_INT(&body->covalentBody->shell, fp);
         FWRITE_INT(&body->covalentBody->orbital, fp);
      }
   }
   j = (int)thermals.size();
   FWRITE_INT(&j, fp);
   for (i = 0; i < j; i++)
   {
      thermal = thermals[i];
      thermal->save(fp);
   }
}


// Import chemistry into current one.
void Chemistry::import(FILE *fp)
{
   int       i, j;
   Chemistry *chemistry;
   Atom      *atom;
   Thermal   *thermal;

#ifdef THREADS
   chemistry = new Chemistry(vesselRadius, randomSeed, 1);
#else
   chemistry = new Chemistry(vesselRadius, randomSeed);
#endif
   assert(chemistry != NULL);
   chemistry->load(fp);

   for (i = 0, j = (int)chemistry->atoms.size(); i < j; i++)
   {
      atom = chemistry->atoms[i];
      addAtom(atom);
   }

   for (i = 0, j = (int)chemistry->thermals.size(); i < j; i++)
   {
      thermal = chemistry->thermals[i];
      addThermal(thermal);
   }

   chemistry->atoms.clear();
   chemistry->thermals.clear();
   delete chemistry;
}


// Get molecule statistics.
void Chemistry::getMoleculeStats(int& num, int& numClosed,
                                 int& numTypes, int& numClosedTypes,
                                 float& aveSize, float& aveClosedSize)
{
   int      i, j, j2, s, cs;
   Molecule *m;

   vector<Molecule *> u, cu;

   generateMolecules();
   s = cs = numClosed = numTypes = numClosedTypes = 0;
   for (i = 0, num = (int)molecules.size(); i < num; i++)
   {
      m  = molecules[i];
      s += m->size();
      for (j = 0, j2 = (int)u.size(); j < j2; j++)
      {
         if (m->equals(u[j]))
         {
            break;
         }
      }
      if (j == j2)
      {
         numTypes++;
         u.push_back(m);
      }
      if (m->isClosed())
      {
         numClosed++;
         cs += m->size();
         for (j = 0, j2 = (int)cu.size(); j < j2; j++)
         {
            if (m->equals(cu[j]))
            {
               break;
            }
         }
         if (j == j2)
         {
            numClosedTypes++;
            cu.push_back(m);
         }
      }
   }
   if (num > 0)
   {
      aveSize = (float)s / (float)num;
   }
   else
   {
      aveSize = 0.0f;
   }
   if (numClosed > 0)
   {
      aveClosedSize = (float)cs / (float)numClosed;
   }
   else
   {
      aveClosedSize = 0.0f;
   }
}


#if (O2_MOLECULES)
// Count O2 molecules.
int Chemistry::countO2()
{
   int  i, i2, j, j2, k, o2;
   Atom *atom, *atom1, *atom2;
   Body *body1, *body2;

   for (i = o2 = 0, i2 = atoms.size(); i < i2; i++)
   {
      // Oxygen?
      atom = atoms[i];
      if (atom->number != 6)
      {
         continue;
      }
      body1 = body2 = NULL;
      for (j = 0, j2 = atom->shells[0].orbitals.size(); j < j2; j++)
      {
         if (body1 == NULL)
         {
            body1 = &atom->shells[0].orbitals[j];
            body1 = body1->covalentBody;
         }
         else if (body2 == NULL)
         {
            body2 = &atom->shells[0].orbitals[j];
            body2 = body2->covalentBody;
         }
         else
         {
            break;
         }
      }
      if ((body1 == NULL) || (body2 == NULL))
      {
         continue;
      }
      for (k = 0; k < i2; k++)
      {
         if (atoms[k]->getID() == body1->id)
         {
            atom1 = atoms[k];
         }
         if (atoms[k]->getID() == body2->id)
         {
            atom2 = atoms[k];
         }
      }
      if ((atom1->number == 6) && (atom2->number == 6) &&
          (atom1->getID() == atom2->getID()) &&
          (atom->getID() > atom1->getID()))
      {
         o2++;
      }
   }
   return(o2);
}


#endif

#if (H2O_MOLECULES)
// Count H2O molecules
int Chemistry::countH2O()
{
   int  i, i2, j, j2, k, h2o;
   Atom *atom, *atom1, *atom2;
   Body *body1, *body2;

   h2o = 0;
   for (i = 0, i2 = atoms.size(); i < i2; i++)
   {
      atom = atoms[i];
      if (atom->number != 6)
      {
         continue;                                // oxygen?
      }
      body1 = body2 = NULL;
      for (j = 0, j2 = atom->shells[0].orbitals.size(); j < j2; j++)
      {
         if (body1 == NULL)
         {
            body1 = &atom->shells[0].orbitals[j];
            body1 = body1->covalentBody;
         }
         else if (body2 == NULL)
         {
            body2 = &atom->shells[0].orbitals[j];
            body2 = body2->covalentBody;
         }
         else
         {
            break;
         }
      }
      if ((body1 == NULL) || (body2 == NULL))
      {
         continue;
      }
      for (k = 0; k < i2; k++)
      {
         if (atoms[k]->getID() == body1->id)
         {
            atom1 = atoms[k];
         }
         if (atoms[k]->getID() == body2->id)
         {
            atom2 = atoms[k];
         }
      }
      if ((atom1->number == 1) && (atom2->number == 1))
      {
         h2o++;
      }
   }
   return(h2o);
}


#endif

#if (CO2_MOLECULES)
// Count CO2 molecules.
int Chemistry::countCO2()
{
   int  i, i2, j, j2, k, co2;
   Atom *atom, *atom1, *atom2, *atom3, *atom4;
   Body *body1, *body2, *body3, *body4;

   co2 = 0;
   for (i = 0, i2 = atoms.size(); i < i2; i++)
   {
      atom = atoms[i];
      if (atom->number != 4)
      {
         continue;                                // carbon?
      }
      body1 = body2 = body3 = body4 = NULL;
      for (j = 0, j2 = atom->shells[0].orbitals.size(); j < j2; j++)
      {
         if (body1 == NULL)
         {
            body1 = &atom->shells[0].orbitals[j];
            body1 = body1->covalentBody;
         }
         else if (body2 == NULL)
         {
            body2 = &atom->shells[0].orbitals[j];
            body2 = body2->covalentBody;
         }
         else if (body3 == NULL)
         {
            body3 = &atom->shells[0].orbitals[j];
            body3 = body3->covalentBody;
         }
         else if (body4 == NULL)
         {
            body4 = &atom->shells[0].orbitals[j];
            body4 = body4->covalentBody;
         }
         else
         {
            break;
         }
      }
      if ((body1 == NULL) || (body2 == NULL) ||
          (body3 == NULL) || (body4 == NULL))
      {
         continue;
      }
      for (k = 0; k < i2; k++)
      {
         if (atoms[k]->getID() == body1->id)
         {
            atom1 = atoms[k];
         }
         if (atoms[k]->getID() == body2->id)
         {
            atom2 = atoms[k];
         }
         if (atoms[k]->getID() == body3->id)
         {
            atom3 = atoms[k];
         }
         if (atoms[k]->getID() == body4->id)
         {
            atom4 = atoms[k];
         }
      }
      if ((atom1->number == 6) && (atom2->number == 6) &&
          (atom3->number == 6) && (atom4->number == 6))
      {
         if (((atom1 == atom2) && (atom3 == atom4)) ||
             ((atom1 == atom3) && (atom2 == atom4)) ||
             ((atom1 == atom4) && (atom2 == atom3)))
         {
            co2++;
         }
      }
   }
   return(co2);
}


#endif
