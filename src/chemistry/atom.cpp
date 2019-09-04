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
 * Atom.
 */

#include "atom.hpp"
using namespace affinity;

enum {
   ATOM_DEPLOYMENT_STEPS = 100, BODY_SHIFT_TRIES = 10
};

// Constructor.
Atom::Atom(Parameters *parameters, int id, int protons, Random *randomizer)
{
   int    i, j, k, n;
   float  d;
   Vector x, f;
   Body   *b1, *b2;

   vector<Body *> bodies;

   // Create nucleus.
   this->parameters = parameters;
   assert(protons >= parameters->MIN_NUCLEUS_PROTONS &&
          protons <= parameters->MAX_NUCLEUS_PROTONS);
   number             = protons;
   nucleus.parameters = parameters;
   nucleus.id         = id;
   nucleus.mass       = (float)number * parameters->PROTON_MASS;
   nucleus.radius     = parameters->NUCLEUS_BODY_RADIUS;
   nucleus.charge     = parameters->PROTON_CHARGE;
   mark = -1;

   // Generate atom color based on atomic number.
   generateColor();

   // Create orbital shells.
   if (number <= 8)
   {
      shells.resize(1);
      shells[0].number = 0;
   }
   else
   {
      shells.resize(2);
      shells[0].number = 0;
      shells[1].number = 1;
   }
   switch (number)
   {
   case 1:
      shells[0].orbitals.resize(1);
      shells[0].orbitals[0].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[0].charge = parameters->ELECTRON_CHARGE;
      getValence(shells[0].orbitals[0].valence[0],
                 shells[0].orbitals[0].valence[1]);
      shells[0].orbitals[0].hasValence = true;
      break;

   case 2:
      shells[0].orbitals.resize(2);
      shells[0].orbitals[0].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[0].charge = parameters->ELECTRON_CHARGE;
      shells[0].orbitals[1].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[1].charge = parameters->ELECTRON_CHARGE;
      getValence(shells[0].orbitals[0].valence[0],
                 shells[0].orbitals[0].valence[1]);
      getValence(shells[0].orbitals[1].valence[0],
                 shells[0].orbitals[1].valence[1]);
      shells[0].orbitals[0].hasValence    =
         shells[0].orbitals[1].hasValence = true;
      break;

   case 3:
      shells[0].orbitals.resize(3);
      shells[0].orbitals[0].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[0].charge = parameters->ELECTRON_CHARGE;
      shells[0].orbitals[1].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[1].charge = parameters->ELECTRON_CHARGE;
      shells[0].orbitals[2].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[2].charge = parameters->ELECTRON_CHARGE;
      getValence(shells[0].orbitals[0].valence[0],
                 shells[0].orbitals[0].valence[1]);
      getValence(shells[0].orbitals[1].valence[0],
                 shells[0].orbitals[1].valence[1]);
      getValence(shells[0].orbitals[2].valence[0],
                 shells[0].orbitals[2].valence[1]);
      shells[0].orbitals[0].hasValence       =
         shells[0].orbitals[1].hasValence    =
            shells[0].orbitals[2].hasValence = true;
      break;

   case 4:
      shells[0].orbitals.resize(4);
      shells[0].orbitals[0].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[0].charge = parameters->ELECTRON_CHARGE;
      shells[0].orbitals[1].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[1].charge = parameters->ELECTRON_CHARGE;
      shells[0].orbitals[2].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[2].charge = parameters->ELECTRON_CHARGE;
      shells[0].orbitals[3].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[3].charge = parameters->ELECTRON_CHARGE;
      getValence(shells[0].orbitals[0].valence[0],
                 shells[0].orbitals[0].valence[1]);
      getValence(shells[0].orbitals[1].valence[0],
                 shells[0].orbitals[1].valence[1]);
      getValence(shells[0].orbitals[2].valence[0],
                 shells[0].orbitals[2].valence[1]);
      getValence(shells[0].orbitals[3].valence[0],
                 shells[0].orbitals[3].valence[1]);
      shells[0].orbitals[0].hasValence          =
         shells[0].orbitals[1].hasValence       =
            shells[0].orbitals[2].hasValence    =
               shells[0].orbitals[3].hasValence = true;
      break;

   case 5:
      shells[0].orbitals.resize(4);
      shells[0].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[1].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[1].charge = parameters->ELECTRON_CHARGE;
      shells[0].orbitals[2].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[2].charge = parameters->ELECTRON_CHARGE;
      shells[0].orbitals[3].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[3].charge = parameters->ELECTRON_CHARGE;
      getValence(shells[0].orbitals[1].valence[0],
                 shells[0].orbitals[1].valence[1]);
      getValence(shells[0].orbitals[2].valence[0],
                 shells[0].orbitals[2].valence[1]);
      getValence(shells[0].orbitals[3].valence[0],
                 shells[0].orbitals[3].valence[1]);
      shells[0].orbitals[1].hasValence       =
         shells[0].orbitals[2].hasValence    =
            shells[0].orbitals[3].hasValence = true;
      break;

   case 6:
      shells[0].orbitals.resize(4);
      shells[0].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[1].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[1].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[2].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[2].charge = parameters->ELECTRON_CHARGE;
      shells[0].orbitals[3].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[3].charge = parameters->ELECTRON_CHARGE;
      getValence(shells[0].orbitals[2].valence[0],
                 shells[0].orbitals[2].valence[1]);
      getValence(shells[0].orbitals[3].valence[0],
                 shells[0].orbitals[3].valence[1]);
      shells[0].orbitals[2].hasValence    =
         shells[0].orbitals[3].hasValence = true;
      break;

   case 7:
      shells[0].orbitals.resize(4);
      shells[0].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[1].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[1].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[2].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[2].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[3].mass   = parameters->ELECTRON_MASS;
      shells[0].orbitals[3].charge = parameters->ELECTRON_CHARGE;
      getValence(shells[0].orbitals[3].valence[0],
                 shells[0].orbitals[3].valence[1]);
      shells[0].orbitals[3].hasValence = true;
      break;

   case 8:
      shells[0].orbitals.resize(4);
      shells[0].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[1].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[1].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[2].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[2].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[3].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[3].charge = parameters->ELECTRON_CHARGE * 2.0f;
      break;
   }
   if (number > 8)
   {
      shells[0].orbitals.resize(4);
      shells[0].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[1].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[1].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[2].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[2].charge = parameters->ELECTRON_CHARGE * 2.0f;
      shells[0].orbitals[3].mass   = parameters->ELECTRON_MASS * 2.0f;
      shells[0].orbitals[3].charge = parameters->ELECTRON_CHARGE * 2.0f;
      switch (number)
      {
      case 9:
         shells[1].orbitals.resize(1);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[0].valence[0],
                    shells[1].orbitals[0].valence[1]);
         shells[1].orbitals[0].hasValence = true;
         break;

      case 10:
         shells[1].orbitals.resize(2);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[0].valence[0],
                    shells[1].orbitals[0].valence[1]);
         getValence(shells[1].orbitals[1].valence[0],
                    shells[1].orbitals[1].valence[1]);
         shells[1].orbitals[0].hasValence    =
            shells[1].orbitals[1].hasValence = true;
         break;

      case 11:
         shells[1].orbitals.resize(3);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[0].valence[0],
                    shells[1].orbitals[0].valence[1]);
         getValence(shells[1].orbitals[1].valence[0],
                    shells[1].orbitals[1].valence[1]);
         getValence(shells[1].orbitals[2].valence[0],
                    shells[1].orbitals[2].valence[1]);
         shells[1].orbitals[0].hasValence       =
            shells[1].orbitals[1].hasValence    =
               shells[1].orbitals[2].hasValence = true;
         break;

      case 12:
         shells[1].orbitals.resize(4);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[3].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[3].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[0].valence[0],
                    shells[1].orbitals[0].valence[1]);
         getValence(shells[1].orbitals[1].valence[0],
                    shells[1].orbitals[1].valence[1]);
         getValence(shells[1].orbitals[2].valence[0],
                    shells[1].orbitals[2].valence[1]);
         getValence(shells[1].orbitals[3].valence[0],
                    shells[1].orbitals[3].valence[1]);
         shells[1].orbitals[0].hasValence          =
            shells[1].orbitals[1].hasValence       =
               shells[1].orbitals[2].hasValence    =
                  shells[1].orbitals[3].hasValence = true;
         break;

      case 13:
         shells[1].orbitals.resize(5);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[3].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[3].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[4].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[4].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[0].valence[0],
                    shells[1].orbitals[0].valence[1]);
         getValence(shells[1].orbitals[1].valence[0],
                    shells[1].orbitals[1].valence[1]);
         getValence(shells[1].orbitals[2].valence[0],
                    shells[1].orbitals[2].valence[1]);
         getValence(shells[1].orbitals[3].valence[0],
                    shells[1].orbitals[3].valence[1]);
         getValence(shells[1].orbitals[4].valence[0],
                    shells[1].orbitals[4].valence[1]);
         shells[1].orbitals[0].hasValence             =
            shells[1].orbitals[1].hasValence          =
               shells[1].orbitals[2].hasValence       =
                  shells[1].orbitals[3].hasValence    =
                     shells[1].orbitals[4].hasValence = true;
         break;

      case 14:
         shells[1].orbitals.resize(6);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[3].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[3].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[4].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[4].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[5].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[5].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[0].valence[0],
                    shells[1].orbitals[0].valence[1]);
         getValence(shells[1].orbitals[1].valence[0],
                    shells[1].orbitals[1].valence[1]);
         getValence(shells[1].orbitals[2].valence[0],
                    shells[1].orbitals[2].valence[1]);
         getValence(shells[1].orbitals[3].valence[0],
                    shells[1].orbitals[3].valence[1]);
         getValence(shells[1].orbitals[4].valence[0],
                    shells[1].orbitals[4].valence[1]);
         getValence(shells[1].orbitals[5].valence[0],
                    shells[1].orbitals[5].valence[1]);
         shells[1].orbitals[0].hasValence                =
            shells[1].orbitals[1].hasValence             =
               shells[1].orbitals[2].hasValence          =
                  shells[1].orbitals[3].hasValence       =
                     shells[1].orbitals[4].hasValence    =
                        shells[1].orbitals[5].hasValence = true;
         break;

      case 15:
         shells[1].orbitals.resize(6);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[3].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[3].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[4].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[4].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[5].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[5].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[1].valence[0],
                    shells[1].orbitals[1].valence[1]);
         getValence(shells[1].orbitals[2].valence[0],
                    shells[1].orbitals[2].valence[1]);
         getValence(shells[1].orbitals[3].valence[0],
                    shells[1].orbitals[3].valence[1]);
         getValence(shells[1].orbitals[4].valence[0],
                    shells[1].orbitals[4].valence[1]);
         getValence(shells[1].orbitals[5].valence[0],
                    shells[1].orbitals[5].valence[1]);
         shells[1].orbitals[1].hasValence             =
            shells[1].orbitals[2].hasValence          =
               shells[1].orbitals[3].hasValence       =
                  shells[1].orbitals[4].hasValence    =
                     shells[1].orbitals[5].hasValence = true;
         break;

      case 16:
         shells[1].orbitals.resize(6);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[3].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[3].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[4].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[4].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[5].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[5].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[2].valence[0],
                    shells[1].orbitals[2].valence[1]);
         getValence(shells[1].orbitals[3].valence[0],
                    shells[1].orbitals[3].valence[1]);
         getValence(shells[1].orbitals[4].valence[0],
                    shells[1].orbitals[4].valence[1]);
         getValence(shells[1].orbitals[5].valence[0],
                    shells[1].orbitals[5].valence[1]);
         shells[1].orbitals[2].hasValence          =
            shells[1].orbitals[3].hasValence       =
               shells[1].orbitals[4].hasValence    =
                  shells[1].orbitals[5].hasValence = true;
         break;

      case 17:
         shells[1].orbitals.resize(6);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[3].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[3].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[4].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[4].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[5].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[5].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[3].valence[0],
                    shells[1].orbitals[3].valence[1]);
         getValence(shells[1].orbitals[4].valence[0],
                    shells[1].orbitals[4].valence[1]);
         getValence(shells[1].orbitals[5].valence[0],
                    shells[1].orbitals[5].valence[1]);
         shells[1].orbitals[3].hasValence       =
            shells[1].orbitals[4].hasValence    =
               shells[1].orbitals[5].hasValence = true;
         break;

      case 18:
         shells[1].orbitals.resize(6);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[3].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[3].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[4].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[4].charge = parameters->ELECTRON_CHARGE;
         shells[1].orbitals[5].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[5].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[4].valence[0],
                    shells[1].orbitals[4].valence[1]);
         getValence(shells[1].orbitals[5].valence[0],
                    shells[1].orbitals[5].valence[1]);
         shells[1].orbitals[4].hasValence    =
            shells[1].orbitals[5].hasValence = true;
         break;

      case 19:
         shells[1].orbitals.resize(6);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[3].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[3].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[4].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[4].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[5].mass   = parameters->ELECTRON_MASS;
         shells[1].orbitals[5].charge = parameters->ELECTRON_CHARGE;
         getValence(shells[1].orbitals[5].valence[0],
                    shells[1].orbitals[5].valence[1]);
         shells[1].orbitals[5].hasValence = true;
         break;

      case 20:
         shells[1].orbitals.resize(6);
         shells[1].orbitals[0].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[0].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[1].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[1].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[2].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[2].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[3].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[3].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[4].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[4].charge = parameters->ELECTRON_CHARGE * 2.0f;
         shells[1].orbitals[5].mass   = parameters->ELECTRON_MASS * 2.0f;
         shells[1].orbitals[5].charge = parameters->ELECTRON_CHARGE * 2.0f;
         break;
      }
   }

   // Finish initialization and deploy the orbitals around the nucleus.
   bodies.push_back(&nucleus);
   for (i = 0; i < (int)shells.size(); i++)
   {
      for (j = 0; j < (int)shells[i].orbitals.size(); j++)
      {
         b1             = &shells[i].orbitals[j];
         b1->parameters = parameters;
         b1->id         = id;
         b1->radius     = parameters->ORBITAL_BODY_RADIUS;
         b1->shell      = i;
         b1->orbital    = j;
         d = parameters->BOND_LENGTH * (float)(i + 1);
         do
         {
            b1->position.x = (float)randomizer->RAND_INTERVAL(-1.0f, 1.0f);
            b1->position.y = (float)randomizer->RAND_INTERVAL(-1.0f, 1.0f);
            b1->position.z = (float)randomizer->RAND_INTERVAL(-1.0f, 1.0f);
         } while (b1->position.Magnitude() < tol);
         b1->position.Normalize(d);
         bodies.push_back(b1);
      }
   }
   for (n = 0; n < ATOM_DEPLOYMENT_STEPS; n++)
   {
      for (i = 0; i < (int)bodies.size(); i++)
      {
         b1 = bodies[i];
         for (j = i + 1; j < (int)bodies.size(); j++)
         {
            b2 = bodies[j];
            if ((b1->id == b2->id) && (b1->shell != b2->shell))
            {
               continue;
            }
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
            f = x * (b1->charge * b2->charge) *
                (float)exp(-(double)((d * d) /
                                     (parameters->CHARGE_GAUSSIAN_SPREAD *
                                      parameters->CHARGE_GAUSSIAN_SPREAD)));
            b1->forces -= f;
            b2->forces += f;
         }
      }
      updateOrbitalBonds();
      for (i = 0; i < (int)bodies.size(); i++)
      {
         b1 = bodies[i];
         b1->update(parameters->UPDATE_STEP);
      }
   }
   for (i = 0; i < (int)shells.size(); i++)
   {
      for (j = 0; j < (int)shells[i].orbitals.size(); j++)
      {
         b1            = &shells[i].orbitals[j];
         b1->position -= nucleus.position;
         b1->velocity.Zero();
      }
   }
   nucleus.position.Zero();
   nucleus.velocity.Zero();
}


// Default constructor.
Atom::Atom(Parameters *parameters)
{
   this->parameters = parameters;
   number           = 0;
   for (int i = 0; i < 4; i++)
   {
      color[i] = 0.0f;
   }
   mark = -1;
}


// Destructor.
Atom::~Atom()
{
   int i, j;

   for (i = 0, j = (int)shells.size(); i < j; i++)
   {
      shells[i].orbitals.clear();
   }
   shells.clear();
}


// Get ID.
int Atom::getID()
{
   return(nucleus.id);
}


// Set ID.
void Atom::setID(int id)
{
   int  i, i2, j, j2;
   Body *b;

   nucleus.id = id;
   for (i = 0, i2 = (int)shells.size(); i < i2; i++)
   {
      for (j = 0, j2 = (int)shells[i].orbitals.size(); j < j2; j++)
      {
         b     = &shells[i].orbitals[j];
         b->id = id;
      }
   }
}


// Get position.
Vector& Atom::getPosition()
{
   return(nucleus.getPosition());
}


// Set position.
void Atom::setPosition(Vector& position)
{
   int    i, i2, j, j2;
   Body   *b;
   Vector d, p;

   d = position - nucleus.position;
   nucleus.setPosition(position);
   for (i = 0, i2 = (int)shells.size(); i < i2; i++)
   {
      for (j = 0, j2 = (int)shells[i].orbitals.size(); j < j2; j++)
      {
         b = &shells[i].orbitals[j];
         p = b->getPosition() + d;
         b->setPosition(p);
      }
   }
}


// Get velocity.
Vector& Atom::getVelocity()
{
   return(nucleus.getVelocity());
}


// Set velocity.
void Atom::setVelocity(Vector& velocity)
{
   int  i, i2, j, j2;
   Body *b;

   nucleus.setVelocity(velocity);
   for (i = 0, i2 = (int)shells.size(); i < i2; i++)
   {
      for (j = 0, j2 = (int)shells[i].orbitals.size(); j < j2; j++)
      {
         b = &shells[i].orbitals[j];
         b->setVelocity(velocity);
      }
   }
}


// Get parameters.
Parameters *Atom::getParameters()
{
   return(parameters);
}


// Set parameters.
void Atom::setParameters(Parameters *p)
{
   int  i, i2, j, j2;
   Body *b;

   parameters = nucleus.parameters = p;
   for (i = 0, i2 = (int)shells.size(); i < i2; i++)
   {
      for (j = 0, j2 = (int)shells[i].orbitals.size(); j < j2; j++)
      {
         b             = &shells[i].orbitals[j];
         b->parameters = p;
      }
   }
}


// Update atom.
void Atom::update()
{
   update(parameters->UPDATE_STEP);
}


void Atom::update(float step)
{
   int s, s2, o, o2;

   nucleus.update(step);
   for (s = 0, s2 = (int)shells.size(); s < s2; s++)
   {
      for (o = 0, o2 = (int)shells[s].orbitals.size(); o < o2; o++)
      {
         shells[s].orbitals[o].update(step);
      }
   }
}


// Update nucleus-orbital bond forces.
void Atom::updateOrbitalBonds()
{
   int    s, s2, o, o2;
   float  d;
   Vector x, v, f;

   for (s = 0, s2 = (int)shells.size(); s < s2; s++)
   {
      for (o = 0, o2 = (int)shells[s].orbitals.size(); o < o2; o++)
      {
         // Use spring equation.
         x = shells[s].orbitals[o].position - nucleus.position;
         d = x.Magnitude() - (parameters->BOND_LENGTH * (float)(s + 1));
         x.Normalize();
         v = shells[s].orbitals[o].velocity - nucleus.velocity;
         f = (-parameters->BOND_STIFFNESS * d * x) +
             (-parameters->BOND_DAMPER * (v * x) * x);
         nucleus.forces -= f;
         shells[s].orbitals[o].forces += f;
      }
   }
}


// Get orbital valence.
// out = # "surplus" electrons in outer shell
// in = # surplus holes.
void Atom::getValence(float& out, float& in)
{
   assert(parameters->MAX_NUCLEUS_PROTONS == 20);
   switch (number)
   {
   case 1:
      out = 1.0f / 1.0f;
      in  = 7.0f / 1.0f;
      break;

   case 2:
      out = 2.0f / 2.0f;
      in  = 6.0f / 2.0f;
      break;

   case 3:
      out = 3.0f / 3.0f;
      in  = 5.0f / 3.0f;
      break;

   case 4:
      out = 4.0f / 4.0f;
      in  = 4.0f / 4.0f;
      break;

   case 5:
      out = 5.0f / 3.0f;
      in  = 3.0f / 3.0f;
      break;

   case 6:
      out = 6.0f / 2.0f;
      in  = 2.0f / 2.0f;
      break;

   case 7:
      out = 7.0f / 1.0f;
      in  = 1.0f / 1.0f;
      break;

   case 8:
      out = 0.0f;
      in  = 0.0f;
      break;

   case 9:
      out = 1.0f / 1.0f;
      in  = 11.0f / 1.0f;
      break;

   case 10:
      out = 2.0f / 2.0f;
      in  = 10.0f / 2.0f;
      break;

   case 11:
      out = 3.0f / 3.0f;
      in  = 9.0f / 3.0f;
      break;

   case 12:
      out = 4.0f / 4.0f;
      in  = 8.0f / 4.0f;
      break;

   case 13:
      out = 5.0f / 5.0f;
      in  = 7.0f / 5.0f;
      break;

   case 14:
      out = 6.0f / 6.0f;
      in  = 6.0f / 6.0f;
      break;

   case 15:
      out = 7.0f / 5.0f;
      in  = 5.0f / 5.0f;
      break;

   case 16:
      out = 8.0f / 4.0f;
      in  = 4.0f / 4.0f;
      break;

   case 17:
      out = 9.0f / 3.0f;
      in  = 3.0f / 3.0f;
      break;

   case 18:
      out = 10.0f / 2.0f;
      in  = 2.0f / 2.0f;
      break;

   case 19:
      out = 11.0f / 1.0f;
      in  = 1.0f / 1.0f;
      break;

   case 20:
      out = 0.0f;
      in  = 0.0f;
      break;

   default:
      out = 0.0f;
      in  = 0.0f;
      break;
   }
}


// Generate atom color based on atomic number.
void Atom::generateColor()
{
   int i, j, k, r, g, b;

   assert(number >= parameters->MIN_NUCLEUS_PROTONS &&
          number <= parameters->MAX_NUCLEUS_PROTONS);
   color[3] = 1.0f;
   k        = (number + 1) % 3;
   j        = (number + 1) % 2;
   for (i = r = g = b = 0; i < number; i++)
   {
      if (j != k)
      {
         switch (j)
         {
         case 0:
            r++;
            break;

         case 1:
            g++;
            break;

         case 2:
            b++;
            break;
         }
      }
      j = (j + 1) % 3;
   }
   switch (k)
   {
   case 0:
      switch (number % 4)
      {
      case 0:
         break;

      case 1:
         r = -r;
         break;

      case 2:
         g = -g;
         break;

      case 3:
         r = -r;
         g = -g;
         break;
      }
      break;

   case 1:
      switch (number % 4)
      {
      case 0:
         break;

      case 1:
         g = -g;
         break;

      case 2:
         b = -b;
         break;

      case 3:
         g = -g;
         b = -b;
         break;
      }
      break;

   case 2:
      switch (number % 4)
      {
      case 0:
         break;

      case 1:
         b = -b;
         break;

      case 2:
         r = -r;
         break;

      case 3:
         b = -b;
         r = -r;
         break;
      }
      break;
   }
   color[0] = (float)r * 3.0f / (float)parameters->MAX_NUCLEUS_PROTONS;
   color[1] = (float)g * 3.0f / (float)parameters->MAX_NUCLEUS_PROTONS;
   color[2] = (float)b * 3.0f / (float)parameters->MAX_NUCLEUS_PROTONS;
   for (i = 0; i < 3; i++)
   {
      if (color[i] < 0.0f)
      {
         color[i] = 1.0f - color[i];
      }
      color[i] = (color[i] * 0.8f) + 0.2f;
      if (color[i] > 1.0f)
      {
         color[i] = 1.0f;
      }
   }
}


// Draw atom.
// When naming atom, combine shell and orbital into name.
void Atom::draw(bool showOrbitals, bool name)
{
   int     i, i2, j, j2;
   Body    *body;
   GLfloat grey[4] = { 0.5f, 0.5f, 0.5f, 1.0f };

   glMaterialfv(GL_FRONT, GL_AMBIENT, color);
   glLineWidth(2.0f);
   if (showOrbitals)
   {
      // Draw nucleus.
      glPushMatrix();
      if (name)
      {
         glPushName(nucleus.id << 7);
      }
      glTranslatef(nucleus.position.x, nucleus.position.y, nucleus.position.z);
      glutSolidSphere(nucleus.radius, 20, 20);
      if (name)
      {
         glPopName();
      }
      glPopMatrix();

      // Draw orbitals.
      for (i = 0, i2 = (int)shells.size(); i < i2; i++)
      {
         for (j = 0, j2 = (int)shells[i].orbitals.size(); j < j2; j++)
         {
            glPushMatrix();
            glBegin(GL_LINES);
            glVertex3f(nucleus.position.x, nucleus.position.y, nucleus.position.z);
            glVertex3f(shells[i].orbitals[j].position.x,
                       shells[i].orbitals[j].position.y, shells[i].orbitals[j].position.z);
            glEnd();
            if (name)
            {
               glPushName(nucleus.id << 7 | 0x1 << 6 | i << 4 | j);
            }
            glTranslatef(shells[i].orbitals[j].position.x,
                         shells[i].orbitals[j].position.y, shells[i].orbitals[j].position.z);
            if (shells[i].orbitals[j].hasValence)
            {
               glutSolidTorus(shells[i].orbitals[j].radius / 2.0f, shells[i].orbitals[j].radius, 20, 20);
            }
            else
            {
               glutSolidSphere(shells[i].orbitals[j].radius, 20, 20);
            }
            if (name)
            {
               glPopName();
            }
            glPopMatrix();
         }
      }
   }
   else
   {
      glPushMatrix();
      if (name)
      {
         glPushName(nucleus.id << 7);
      }
      glTranslatef(nucleus.position.x, nucleus.position.y, nucleus.position.z);
      if ((shells.size() > 0) && (shells[0].orbitals.size() > 0))
      {
         glutSolidSphere((shells.size() * parameters->BOND_LENGTH) +
                         (shells[0].orbitals[0].radius * 0.5f), 20, 20);
      }
      else
      {
         glutSolidSphere((shells.size() * parameters->BOND_LENGTH) +
                         (nucleus.radius * 0.5f), 20, 20);
      }
      if (name)
      {
         glPopName();
      }
      glPopMatrix();
   }

   // Draw covalent bonds.
   glMaterialfv(GL_FRONT, GL_AMBIENT, grey);
   for (i = 0, i2 = (int)shells.size(); i < i2; i++)
   {
      for (j = 0, j2 = (int)shells[i].orbitals.size(); j < j2; j++)
      {
         body = &shells[i].orbitals[j];
         if (body->covalentBody != NULL)
         {
            glBegin(GL_LINES);
            glVertex3f(body->position.x, body->position.y, body->position.z);
            glVertex3f(body->covalentBody->position.x,
                       body->covalentBody->position.y, body->covalentBody->position.z);
            glEnd();
         }
      }
   }
   glLineWidth(1.0f);
}


// Highlight the atom.
void Atom::highlight(bool showOrbitals, int shell, int orbital)
{
   GLfloat on[4]  = { 0.1f, 0.1f, 0.1f, 1.0f };
   GLfloat off[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

   glMaterialfv(GL_FRONT, GL_EMISSION, on);
   glLineWidth(2.0f);
   if (showOrbitals)
   {
      // Hightlight nucleus?
      if ((shell == -1) || (orbital == -1))
      {
         glPushMatrix();
         glTranslatef(nucleus.position.x, nucleus.position.y, nucleus.position.z);
         glutWireSphere(nucleus.radius * 2.0f, 5, 5);
         glPopMatrix();
      }
      else
      {
         // Highlight orbital.
         glPushMatrix();
         glTranslatef(shells[shell].orbitals[orbital].position.x,
                      shells[shell].orbitals[orbital].position.y, shells[shell].orbitals[orbital].position.z);
         glutWireSphere(shells[shell].orbitals[orbital].radius * 2.0f, 5, 5);
         glPopMatrix();
      }
   }
   else
   {
      glPushMatrix();
      glTranslatef(nucleus.position.x, nucleus.position.y, nucleus.position.z);
      if ((shells.size() > 0) && (shells[0].orbitals.size() > 0))
      {
         glutWireSphere(((shells.size() * parameters->BOND_LENGTH) +
                         (shells[0].orbitals[0].radius * 0.5f)) * 2.0f, 5, 5);
      }
      else
      {
         glutWireSphere(((shells.size() * parameters->BOND_LENGTH) +
                         (nucleus.radius * 0.5f)) * 2.0f, 5, 5);
      }
      glPopMatrix();
   }
   glMaterialfv(GL_FRONT, GL_EMISSION, off);
   glLineWidth(1.0f);
}


// Compare atoms by atomic number.
bool Atom::equals(Atom *atom)
{
   return(number == atom->number);
}


bool Atom::operator==(Atom& atom)
{
   return(equals(&atom));
}


bool Atom::operator!=(Atom& atom)
{
   return(!(*this == atom));
}


// Load atom.
void Atom::load(FILE *fp)
{
   int i, j, p, q;

   FREAD_INT(&number, fp);
   nucleus.parameters = parameters;
   generateColor();
   nucleus.load(fp);
   for (i = 0, j = (int)shells.size(); i < j; i++)
   {
      shells[i].orbitals.clear();
   }
   shells.clear();
   FREAD_INT(&i, fp);
   shells.resize(i);
   for (i = 0, j = (int)shells.size(); i < j; i++)
   {
      FREAD_INT(&shells[i].number, fp);
      FREAD_INT(&p, fp);
      shells[i].orbitals.resize(p);
      for (p = 0, q = (int)shells[i].orbitals.size(); p < q; p++)
      {
         shells[i].orbitals[p].parameters = parameters;
         shells[i].orbitals[p].load(fp);
      }
   }
}


// Save atom.
void Atom::save(FILE *fp)
{
   int i, j, p, q;

   FWRITE_INT(&number, fp);
   nucleus.save(fp);
   j = (int)shells.size();
   FWRITE_INT(&j, fp);
   for (i = 0; i < j; i++)
   {
      FWRITE_INT(&shells[i].number, fp);
      q = (int)shells[i].orbitals.size();
      FWRITE_INT(&q, fp);
      for (p = 0; p < q; p++)
      {
         shells[i].orbitals[p].save(fp);
      }
   }
}
