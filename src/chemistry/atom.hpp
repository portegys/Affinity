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
 * A 3D structure consisting of a nucleus of protons surrounded by
 * a number of electrons shells, each consisting of a number of electron
 * orbitals. The orbitals are connected to the nucleus by electric
 * charge and spring-like forces. Atoms may bond through both charge
 * and valence forces. Valence forces allow electron sharing to achieve
 * integral shell occupancy.
 */

#ifndef __ATOM__
#define __ATOM__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <GL/glut.h>
#include "body.hpp"
#include "../utility/random.hpp"
using namespace std;

namespace affinity
{
class Atom
{
public:

   // Parameters.
   Parameters *parameters;

   // Electron shell.
   class Shell
   {
public:

      int          number;                        // shell number
      vector<Body> orbitals;                      // electron orbitals
   };

   int           number;                          // number of protons
   Body          nucleus;                         // nucleus
   vector<Shell> shells;                          // electron shells
   GLfloat       color[4];                        // color
   int           mark;                            // mark

   // Constructors.
   Atom(Parameters *parameters = NULL);
   Atom(Parameters *parameters, int id, int protons,
        Random *randomizer);

   // Destructor.
   ~Atom();

   // Get and set ID.
   int getID();
   void setID(int id);

   // Get and set position.
   Vector& getPosition();
   void setPosition(Vector& position);

   // Get and set velocity.
   Vector& getVelocity();
   void setVelocity(Vector& velocity);

   // Get and set parameters.
   Parameters *getParameters();
   void setParameters(Parameters *);

   // Update atom.
   void update();
   void update(float step);

   // Update nucleus-orbital bond forces.
   void updateOrbitalBonds();

   // Get orbital valence.
   void getValence(float& out, float& in);

   // Generate color.
   void generateColor();

   // Draw atom.
   void draw(bool showOrbitals, bool name = false);

   // Highlight the atom.
   void highlight(bool showOrbitals, int shell, int orbital);

   // Compare atoms by atomic number.
   bool equals(Atom *atom);
   bool operator==(Atom& atom);
   bool operator!=(Atom& atom);

   // Load and save atom.
   void load(FILE *fp);
   void save(FILE *fp);
};
}
#endif
