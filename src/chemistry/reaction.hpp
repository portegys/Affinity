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
 * Reaction: a recursive molecule decomposition into bonded components.
 * At each level, the molecule is divided by the weakest bonds.
 * This can be used as a "recipe" to create the molecule from atoms.
 */

#ifndef __REACTION__
#define __REACTION__

#include <stdio.h>
#include <stdlib.h>
#include <pair>
#include <vector>
#include <assert.h>
#include "chemistry.hpp"
using namespace std;

namespace affinity
{
class Reaction
{
public:

   // The molecule.
   Chemistry *molecule;

   // Components.

   vector<Chemistry *> components;

   // Constructor.
   Reaction(Molecule *molecule);

   // Destructor.
   ~Reaction();

   // Load and save reaction.
   void load(FILE *fp);
   void save(FILE *fp);

   // Print.
   void print(FILE *fp = stdout);
};
}
#endif
