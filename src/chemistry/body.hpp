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
 * Body.
 * A physical object that can act and be acted upon by charge and valence forces.
 */

#ifndef __BODY__
#define __BODY__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "parameters.hpp"
#include "../utility/vector.hpp"
#include "../utility/fileio.h"

namespace affinity
{
// Body.
class Body
{
public:

   // Parameters.
   Parameters *parameters;

   int    id;                                     // atom id
   int    shell;                                  // shell (-1=nucleus)
   int    orbital;                                // orbital (-1=nucleus)
   float  mass;                                   // mass
   float  radius;                                 // radius
   float  charge;                                 // charge (+|-)
   float  valence[2];                             // valence (export/import)
   bool   hasValence;                             // has valence?
   Body   *covalentBody;                          // covalent bonded body
   Vector position;                               // position
   Vector velocity;                               // velocity
   Vector forces;                                 // impinging forces

   // Constructors.
   Body(Parameters *parameters = NULL);
   Body(Parameters *parameters, float mass, float radius, float charge,
        float *valence, bool hasValence, Vector& position, Vector& velocity);

   // Get and set position.
   Vector& getPosition();
   void setPosition(Vector& position);

   // Get and set velocity.
   Vector& getVelocity();
   void setVelocity(Vector& velocity);

   // Update body.
   void update(float step);

   // Get covalent bonding force with other body.
   float getCovalentForce(Body *);

   // Update covalent bond forces.
   void updateCovalentBond();

   // Load and save body.
   void load(FILE *fp);
   void save(FILE *fp);
};
}
#endif
