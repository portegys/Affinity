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

#include "thermal.hpp"
using namespace affinity;

// Constructors.
Thermal::Thermal(Parameters *parameters)
{
   this->parameters = parameters;
   radius           = temperature = 0.0f;
}


Thermal::Thermal(Parameters *parameters, float radius,
                 Vector& position, float temperature)
{
   this->parameters = parameters;
   assert(radius >= parameters->MIN_THERMAL_RADIUS &&
          radius <= parameters->MAX_THERMAL_RADIUS);
   assert(temperature >= parameters->MIN_THERMAL_TEMPERATURE &&
          temperature <= parameters->MAX_TEMPERATURE);
   this->radius      = radius;
   this->position    = position;
   this->temperature = temperature;
}


// Draw thermal.
void Thermal::draw()
{
   GLfloat color[4];

   color[0] = color[1] = color[2] =
                            temperature / parameters->MAX_TEMPERATURE;
   color[3] = 1.0f;
   glMaterialfv(GL_FRONT, GL_AMBIENT, color);
   glPushMatrix();
   glTranslatef(position.x, position.y, position.z);
   glutSolidSphere(radius, 50, 50);
   glPopMatrix();
}


// Load thermal.
void Thermal::load(FILE *fp)
{
   FREAD_FLOAT(&radius, fp);
   FREAD_FLOAT(&position.x, fp);
   FREAD_FLOAT(&position.y, fp);
   FREAD_FLOAT(&position.z, fp);
   FREAD_FLOAT(&temperature, fp);
}


// Save thermal.
void Thermal::save(FILE *fp)
{
   FWRITE_FLOAT(&radius, fp);
   FWRITE_FLOAT(&position.x, fp);
   FWRITE_FLOAT(&position.y, fp);
   FWRITE_FLOAT(&position.z, fp);
   FWRITE_FLOAT(&temperature, fp);
}
