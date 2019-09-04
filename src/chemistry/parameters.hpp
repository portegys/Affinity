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
 * Affinity chemistry parameters.
 */

#ifndef __PARAMETERS__
#define __PARAMETERS__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../utility/fileio.h"

// Define chemistry-specific compilation directives here.

// Create O2 molecules.
#define O2_MOLECULES         0

// Create H2O molecules.
#define H2O_MOLECULES        0

// Create CO2 molecules.
#define CO2_MOLECULES        0

// Create "organic" molecules.
#define ORGANIC_MOLECULES    1

namespace affinity
{
class Parameters
{
public:

   static const float DEFAULT_MIN_ATOM_INITIAL_FORCE;
   float              MIN_ATOM_INITIAL_FORCE;
   static const float DEFAULT_MAX_ATOM_INITIAL_FORCE;
   float              MAX_ATOM_INITIAL_FORCE;
   static const int   DEFAULT_MIN_NUCLEUS_PROTONS;
   int                MIN_NUCLEUS_PROTONS;
   static const int   DEFAULT_MAX_NUCLEUS_PROTONS;
   int                MAX_NUCLEUS_PROTONS;
   static const float DEFAULT_PROTON_MASS;
   float              PROTON_MASS;
   static const float DEFAULT_ELECTRON_MASS;
   float              ELECTRON_MASS;
   static const float DEFAULT_PROTON_CHARGE;
   float              PROTON_CHARGE;
   static const float DEFAULT_ELECTRON_CHARGE;
   float              ELECTRON_CHARGE;
   static const float DEFAULT_CHARGE_GAUSSIAN_SPREAD;
   float              CHARGE_GAUSSIAN_SPREAD;
   static const float DEFAULT_NUCLEUS_BODY_RADIUS;
   float              NUCLEUS_BODY_RADIUS;
   static const float DEFAULT_ORBITAL_BODY_RADIUS;
   float              ORBITAL_BODY_RADIUS;
   static const float DEFAULT_MAX_BODY_RANGE;
   float              MAX_BODY_RANGE;
   static const float DEFAULT_BOND_LENGTH;
   float              BOND_LENGTH;
   static const float DEFAULT_BOND_STIFFNESS;
   float              BOND_STIFFNESS;
   static const float DEFAULT_BOND_DAMPER;
   float              BOND_DAMPER;
   static const float DEFAULT_NUCLEAR_REPULSION_STIFFNESS;
   float              NUCLEAR_REPULSION_STIFFNESS;
   static const float DEFAULT_COVALENT_BONDING_RANGE;
   float              COVALENT_BONDING_RANGE;
   static const float DEFAULT_MIN_COVALENT_BOND_FORCE;
   float              MIN_COVALENT_BOND_FORCE;
   static const float DEFAULT_COVALENT_BOND_STIFFNESS_SCALE;
   float              COVALENT_BOND_STIFFNESS_SCALE;
   static const float DEFAULT_MIN_THERMAL_RADIUS;
   float              MIN_THERMAL_RADIUS;
   static const float DEFAULT_MAX_THERMAL_RADIUS;
   float              MAX_THERMAL_RADIUS;
   static const float DEFAULT_MIN_THERMAL_TEMPERATURE;
   float              MIN_THERMAL_TEMPERATURE;
   static const float DEFAULT_MAX_TEMPERATURE;
   float              MAX_TEMPERATURE;
   static const float DEFAULT_UPDATE_STEP;
   float              UPDATE_STEP;

   // Constructor.
   Parameters();

   // Load and save parameters.
   void load(FILE *fp);
   void save(FILE *fp);
};
}
#endif
