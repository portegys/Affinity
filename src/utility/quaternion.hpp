//***************************************************************************//
//* File Name: quaternion.hpp                                               *//
//* Author:    Tom Portegys, portegys@ilstu.edu                             *//
//* Date Made: 07/25/02                                                     *//
//* File Desc: Class declaration representing a quaternion.                 *//
//* Rev. Date:                                                              *//
//* Rev. Desc:                                                              *//
//*                                                                         *//
//***************************************************************************//

#ifndef __QUATERNION_HPP__
#define __QUATERNION_HPP__

#include <math.h>
#include "vector.hpp"

class cQuaternion
{
public:

   // Quaternion.
   float quat[4];

   // Constructors.
   cQuaternion()
   {
      quat[0] = 0.0;
      quat[1] = 0.0;
      quat[2] = 0.0;
      quat[3] = 1.0;
   }


   cQuaternion(float q[4])
   {
      quat[0] = q[0];
      quat[1] = q[1];
      quat[2] = q[2];
      quat[3] = q[3];
   }


   // Vector operations.
   void vzero(float *);
   void vset(float *, float, float, float);
   void vsub(const float *, const float *, float *);
   void vcopy(const float *, float *);
   void vcross(const float *, const float *, float *);
   float vlength(const float *);
   void vscale(float *, float);
   void vnormal(float *);
   float vdot(const float *, const float *);
   void vadd(const float *, const float *, float *);

   // Clear.
   void clear()
   {
      quat[0] = 0.0;
      quat[1] = 0.0;
      quat[2] = 0.0;
      quat[3] = 1.0;
   }


   // Add quaternions.
   void add_quats(cQuaternion& q1, cQuaternion& q2, cQuaternion& dest);

   // Multiply quaternions.
   void mult_quats(cQuaternion& q1, cQuaternion& q2, cQuaternion& dest);

   // Normalize a quaternion.
   void normalize_quat();

   // Build a rotation matrix given a quaternion rotation.
   void build_rotmatrix(float m[4][4]);

   // Load an axis-angle rotation into quaternion.
   void loadRotation(float angle, float *axis);

   // Merge an axis-angle rotation into quaternion.
   void mergeRotation(float angle, float *axis);

   // Make quaternion from Euler angles.
   void makeQFromEulerAngles(float pitch, float yaw, float roll);

   // Make Euler angles from quaternion.
   void makeEulerAnglesFromQ(float& pitch, float& yaw, float& roll);

private:

   int NormalCount;                               // Normalization counter.
};
#endif
