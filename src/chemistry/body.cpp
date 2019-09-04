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
 */

#include "body.hpp"
#include "atom.hpp"
using namespace affinity;

// Constructors.
Body::Body(Parameters *parameters)
{
   this->parameters = parameters;
   id           = shell = orbital = -1;
   mass         = radius = charge = valence[0] = valence[1] = 0.0f;
   hasValence   = false;
   covalentBody = NULL;
}


Body::Body(Parameters *parameters, float mass, float radius, float charge,
           float *valence, bool hasValence, Vector& position, Vector& velocity)
{
   this->parameters = parameters;
   id               = shell = orbital = -1;
   this->mass       = mass;
   this->radius     = radius;
   this->charge     = charge;
   this->valence[0] = valence[0];
   this->valence[1] = valence[1];
   this->hasValence = hasValence;
   covalentBody     = NULL;
   this->position   = position;
   this->velocity   = velocity;
}


// Get position.
Vector& Body::getPosition()
{
   return(position);
}


// Set position.
void Body::setPosition(Vector& position)
{
   this->position = position;
}


// Get velocity.
Vector& Body::getVelocity()
{
   return(velocity);
}


// Set velocity.
void Body::setVelocity(Vector& velocity)
{
   this->velocity = velocity;
}


// Update body.
void Body::update(float step)
{
   velocity += forces / mass;
   if (velocity.Magnitude() > parameters->MAX_TEMPERATURE)
   {
      velocity.Normalize(parameters->MAX_TEMPERATURE);
   }
   position += velocity * step;
   forces.Zero();
}


// Get covalent bonding force with other body.
float Body::getCovalentForce(Body *body)
{
   if (!hasValence || !body->hasValence)
   {
      return(0.0f);
   }
#ifdef NEVER
   return(min(valence[0], body->valence[1]) +
          min(body->valence[0], valence[1]));
#endif
   return(((fabs((float)valence[0] - (float)body->valence[0]) +
            fabs((float)valence[1] - (float)body->valence[1])) / 2.0f) +
          parameters->MIN_COVALENT_BOND_FORCE);
}


// Update covalent bond forces.
void Body::updateCovalentBond()
{
   float  d;
   Vector x, v, f;

   // Use spring equation.
   if (covalentBody == NULL)
   {
      return;
   }
   x = covalentBody->position - position;
   d = x.Magnitude();
   if (d < tol)
   {
      return;
   }
   x.Normalize();
   v = covalentBody->velocity - velocity;
   f = (-getCovalentForce(covalentBody) *
        parameters->COVALENT_BOND_STIFFNESS_SCALE * d * x);
   forces -= f;
   covalentBody->forces += f;
}


// Load body.
void Body::load(FILE *fp)
{
   FREAD_INT(&id, fp);
   FREAD_INT(&shell, fp);
   FREAD_INT(&orbital, fp);
   FREAD_FLOAT(&mass, fp);
   FREAD_FLOAT(&radius, fp);
   FREAD_FLOAT(&charge, fp);
   FREAD_FLOAT(&valence[0], fp);
   FREAD_FLOAT(&valence[1], fp);
   FREAD_BOOL(&hasValence, fp);
   FREAD_FLOAT(&position.x, fp);
   FREAD_FLOAT(&position.y, fp);
   FREAD_FLOAT(&position.z, fp);
   FREAD_FLOAT(&velocity.x, fp);
   FREAD_FLOAT(&velocity.y, fp);
   FREAD_FLOAT(&velocity.z, fp);
}


// Save body.
void Body::save(FILE *fp)
{
   FWRITE_INT(&id, fp);
   FWRITE_INT(&shell, fp);
   FWRITE_INT(&orbital, fp);
   FWRITE_FLOAT(&mass, fp);
   FWRITE_FLOAT(&radius, fp);
   FWRITE_FLOAT(&charge, fp);
   FWRITE_FLOAT(&valence[0], fp);
   FWRITE_FLOAT(&valence[1], fp);
   FWRITE_BOOL(&hasValence, fp);
   FWRITE_FLOAT(&position.x, fp);
   FWRITE_FLOAT(&position.y, fp);
   FWRITE_FLOAT(&position.z, fp);
   FWRITE_FLOAT(&velocity.x, fp);
   FWRITE_FLOAT(&velocity.y, fp);
   FWRITE_FLOAT(&velocity.z, fp);
}
