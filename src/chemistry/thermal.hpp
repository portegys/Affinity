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
 * Thermal object.
 */

#ifndef __THERMAL__
#define __THERMAL__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <GL/glut.h>
#include "parameters.hpp"
#include "../utility/vector.hpp"
#include "../utility/fileio.h"

namespace affinity
{
// Thermal object.
// Colliding bodies tend to achieve a speed
// equal to the temperature of the thermal.
class Thermal
{
public:

   // Parameters.
   Parameters *parameters;

   float  radius;                                 // radius
   Vector position;                               // position
   float  temperature;                            // temperature

   // Constructors.
   Thermal(Parameters *parameters = NULL);
   Thermal(Parameters *parameters, float radius,
           Vector& position, float temperature);

   // Draw thermal.
   void draw();

   // Load and save thermal.
   void load(FILE *fp);
   void save(FILE *fp);
};
}
#endif
