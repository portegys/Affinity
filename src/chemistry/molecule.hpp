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
 * Molecule: a bonded graph of atoms.
 */

#ifndef __MOLECULE__
#define __MOLECULE__

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "atom.hpp"
#include "../utility/md5.h"
using namespace std;

namespace affinity
{
class Chemistry;

// Molecule: a bonded graph of atoms.
class Molecule
{
public:

   // Molecule atom IDs.
   vector<int> atomIDs;

   // Unique molecule code.
   unsigned char code[MD5_SIZE];

   // Constructor.
   Molecule(Chemistry *chemistry, Atom *atom);

   // Destructor.
   ~Molecule();

   // Get molecule size (number of atoms).
   int size();

   // Molecule contains atom/body?
   bool contains(Atom *atom);
   bool contains(Body *body);

   // Molecule is "closed" (all bonds connected)?
   bool isClosed();

   // Compare molecules by code.
   bool equals(Molecule *molecule);
   bool operator==(Molecule& molecule);
   bool operator!=(Molecule& molecule);

   // Print.
   void print(FILE *fp = stdout);

   // Chemistry.
   Chemistry *chemistry;

   // Get component atom IDs.
   void getIDs(Atom *, vector<int>& ids);

   // Less-than comparison of atoms by atomic number.
   static bool ltcmpAtoms(Atom *a, Atom *b);

   // Atom with its subtree of atoms.
   class AtomTree
   {
public:

      Chemistry          *chemistry;
      Atom               *atom;
      int                parentBonds;
      unsigned char      code[MD5_SIZE];
      vector<AtomTree *> children;
      AtomTree(Chemistry *, Atom *atom = NULL);
      ~AtomTree();
      void generateCode(bool hashNumbers = true);

      bool expanded;
      void expand(vector<int>& path);

private:

      static bool ltcmpCode(AtomTree *a, AtomTree *b);
   };
};
}
#endif
