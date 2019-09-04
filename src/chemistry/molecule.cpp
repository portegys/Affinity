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
 * Molecule.
 */

#include "molecule.hpp"
#include "chemistry.hpp"
using namespace affinity;

// Constructor.
// Analyze molecule given a component atom.
Molecule::Molecule(Chemistry *chemistry, Atom *atom)
{
   int i, i2;

   vector<Atom *> atoms;
   AtomTree       *root, *child;

   // Save chemistry.
   this->chemistry = chemistry;

   // Get component atom IDs.
   getIDs(atom, atomIDs);

   // Sort atoms by atomic number.
   for (i = 0, i2 = (int)atomIDs.size(); i < i2; i++)
   {
      atoms.push_back(chemistry->getAtom(atomIDs[i]));
   }
   sort(atoms.begin(), atoms.end(), ltcmpAtoms);

   // Generate code for molecule using atom trees for each atom.
   root = new AtomTree(chemistry);
   assert(root != NULL);
   for (i = 0, i2 = (int)atoms.size(); i < i2; i++)
   {
      child = new AtomTree(chemistry, atoms[i]);
      assert(child != NULL);
      root->children.push_back(child);
   }
   root->generateCode();
   memcpy(code, root->code, MD5_SIZE);
   delete root;
}


// Destructor.
Molecule::~Molecule()
{
   atomIDs.clear();
}


// Get molecule size (number of atoms).
int Molecule::size()
{
   return((int)atomIDs.size());
}


// Molecule contains atom?
bool Molecule::contains(Atom *atom)
{
   return(contains(&atom->nucleus));
}


// Molecule contains body?
bool Molecule::contains(Body *body)
{
   int i, j;

   for (i = 0, j = (int)atomIDs.size(); i < j &&
        atomIDs[i] <= body->id; i++)
   {
      if (atomIDs[i] == body->id)
      {
         return(true);
      }
   }
   return(false);
}


// Molecule is "closed" (all bonds connected)?
bool Molecule::isClosed()
{
   int  i, i2, j, j2, s;
   Atom *atom;

   for (i = 0, i2 = (int)atomIDs.size(); i < i2; i++)
   {
      atom = chemistry->getAtom(atomIDs[i]);
      s    = (int)atom->shells.size() - 1;
      for (j = 0, j2 = (int)atom->shells[s].orbitals.size(); j < j2; j++)
      {
         if (atom->shells[s].orbitals[j].hasValence &&
             (atom->shells[s].orbitals[j].covalentBody == NULL))
         {
            return(false);
         }
      }
   }
   return(true);
}


// Compare molecules by code.
bool Molecule::equals(Molecule *molecule)
{
   if (memcmp(code, molecule->code, MD5_SIZE) == 0)
   {
      return(true);
   }
   else
   {
      return(false);
   }
}


bool Molecule::operator==(Molecule& molecule)
{
   return(equals(&molecule));
}


bool Molecule::operator!=(Molecule& molecule)
{
   return(!(*this == molecule));
}


// Print molecule.
void Molecule::print(FILE *fp)
{
   int i, j;

   fprintf(fp, "Atom ID/Numbers: ");
   for (i = 0, j = (int)atomIDs.size(); i < j; i++)
   {
      fprintf(fp, "%d/%d ", atomIDs[i], chemistry->getAtom(atomIDs[i])->number);
   }
   fprintf(fp, "\n");
   fprintf(fp, "Code: ");
   for (i = 0; i < MD5_SIZE; i++)
   {
      fprintf(fp, "%d ", code[i]);
   }
   fprintf(fp, "\n");
}


// Get component atom IDs.
void Molecule::getIDs(Atom *atom, vector<int>& ids)
{
   int  i, i2, j, j2;
   Body *body;

   for (i = 0, i2 = (int)ids.size(); i < i2; i++)
   {
      if (ids[i] == atom->getID())
      {
         return;
      }
   }
   ids.push_back(atom->getID());
   sort(ids.begin(), ids.end());

   // Get connected atom IDs.
   i = (int)atom->shells.size() - 1;
   for (j = 0, j2 = (int)atom->shells[i].orbitals.size(); j < j2; j++)
   {
      if ((body = atom->shells[i].orbitals[j].covalentBody) != NULL)
      {
         getIDs(chemistry->getAtom(body->id), ids);
      }
   }
}


// Less-than comparison atoms by atomic number.
bool Molecule::ltcmpAtoms(Atom *a, Atom *b)
{
   return(a->number < b->number);
}


// Atom tree constructor.
Molecule::AtomTree::AtomTree(Chemistry *chemistry, Atom *atom)
{
   this->chemistry = chemistry;
   this->atom      = atom;
   if (atom == NULL)
   {
      parentBonds = 0;
   }
   else
   {
      parentBonds = 1;
   }
   expanded = false;
}


// Atom tree destructor.
Molecule::AtomTree::~AtomTree()
{
   int i, i2;

   for (i = 0, i2 = (int)children.size(); i < i2; i++)
   {
      delete children[i];
   }
   children.clear();
}


// Generate code.
void Molecule::AtomTree::generateCode(bool hashNumbers)
{
   int               i, j, s;
   struct MD5Context md5c;
   unsigned char     *input;

   // Root?
   s = (int)children.size();
   if (atom == NULL)
   {
      // Iteratively expand children until either fully
      // expanded or all have distinct codes.
      vector<int> path;
      while (true)
      {
         // Generate codes without labels to encode structure at this point.
         for (i = 0; i < s; i++)
         {
            children[i]->generateCode(false);
         }
         sort(children.begin(), children.end(), ltcmpCode);
         for (i = 0; i < s - 1; i++)
         {
            if (!children[i]->expanded && (ltcmpCode(children[i], children[i + 1]) == 0))
            {
               break;
            }
         }
         if (i == s - 1)
         {
            break;
         }
         for (j = i + 2; j < s && ltcmpCode(children[i], children[j]) == 0; j++)
         {
         }
         for ( ; i < j; i++)
         {
            path.clear();
            children[i]->expand(path);
         }
      }
   }

   // Recursively generate code.
   for (i = 0; i < s; i++)
   {
      children[i]->generateCode();
   }
   sort(children.begin(), children.end(), ltcmpCode);
   MD5Init(&md5c);
   j = 0;
   if (atom != NULL)
   {
      if (hashNumbers)
      {
         j = 2;
      }
      else
      {
         j = 1;
      }
   }
   input = new unsigned char[(s * MD5_SIZE) + j];
   assert(input != NULL);
   if (atom != NULL)
   {
      if (hashNumbers)
      {
         input[0] = (unsigned char)atom->number;
         input[1] = (unsigned char)parentBonds;
      }
      else
      {
         input[0] = (unsigned char)parentBonds;
      }
   }
   for (i = 0; i < s; i++)
   {
      memcpy(&input[j], children[i]->code, MD5_SIZE);
      j += MD5_SIZE;
   }
   MD5Update(&md5c, input, j);
   MD5Final(code, &md5c);
   delete input;
}


// Expand tree one level deeper.
void Molecule::AtomTree::expand(vector<int>& path)
{
   int      i, i2, j, j2, k, id;
   Body     *body;
   AtomTree *child;

   vector<int>        path2;
   vector<AtomTree *> tmpTree;

   if (expanded)
   {
      return;
   }
   expanded = true;
   if (children.size() == 0)
   {
      // Expand this level.
      i = (int)atom->shells.size() - 1;
      for (j = 0, j2 = (int)atom->shells[i].orbitals.size(); j < j2; j++)
      {
         if ((body = atom->shells[i].orbitals[j].covalentBody) != NULL)
         {
            child = new AtomTree(chemistry, chemistry->getAtom(body->id));
            assert(child != NULL);
            children.push_back(child);
         }
      }

      // Consolidate multiple bonds.
      for (i = 0, i2 = (int)children.size(); i < i2; i++)
      {
         if (children[i] == NULL)
         {
            continue;
         }
         for (j = i + 1, k = 0; j < i2; j++)
         {
            if (children[j] == NULL)
            {
               continue;
            }
            if (children[i]->atom->getID() == children[j]->atom->getID())
            {
               k++;
               delete children[j];
               children[j] = NULL;
            }
         }
         children[i]->parentBonds += k;
      }
      for (i = 0; i < i2; i++)
      {
         if (children[i] != NULL)
         {
            tmpTree.push_back(children[i]);
         }
      }
      children.clear();
      for (i = 0, i2 = (int)tmpTree.size(); i < i2; i++)
      {
         children.push_back(tmpTree[i]);
      }
      if (children.size() > 0)
      {
         expanded = false;
      }
   }
   else
   {
      // Expand deeper.
      id = atom->getID();
      for (i = 0, i2 = (int)path.size(); i < i2; i++)
      {
         if (id == path[i])
         {
            break;
         }
      }
      if (i == i2)
      {
         path.push_back(id);
         for (i = 0, i2 = (int)children.size(); i < i2; i++)
         {
            id = children[i]->atom->getID();
            path2.clear();
            for (j = 0, j2 = (int)path.size(); j < j2; j++)
            {
               path2.push_back(path[j]);
            }
            children[i]->expand(path2);
            if (!children[i]->expanded)
            {
               expanded = false;
            }
         }
      }
   }
}


// Less-than comparison by code.
bool Molecule::AtomTree::ltcmpCode(AtomTree *a, AtomTree *b)
{
   if (memcmp(a->code, b->code, MD5_SIZE) < 0)
   {
      return(true);
   }
   else
   {
      return(false);
   }
}
