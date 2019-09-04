//***************************************************************************//
//* File Name: spacial.cpp                                                  *//
//* Author:    Tom Portegys, portegys@ilstu.edu                             *//
//* Date Made: 07/25/02                                                     *//
//* File Desc: Class implementation details representing spacial            *//
//*            properties: rotation, translation, scale, and speed.         *//
//* Rev. Date:                                                              *//
//* Rev. Desc:                                                              *//
//*                                                                         *//
//***************************************************************************//

#include "spacial.hpp"

// Update rotation and translation state.
void cSpacial::update()
{
   cQuaternion xq, yq, zq, q1, q2;
   GLfloat     v[3];

   v[0] = 1.0f;
   v[1] = 0.0f;
   v[2] = 0.0f;
   xq.loadRotation(pitch * (float)DIV_PI_180, v);
   v[0] = 0.0f;
   v[1] = 1.0f;
   v[2] = 0.0f;
   yq.loadRotation(yaw * (float)DIV_PI_180, v);
   qcalc->mult_quats(xq, yq, q1);
   v[0] = 0.0f;
   v[1] = 0.0f;
   v[2] = 1.0f;
   zq.loadRotation(roll * (float)DIV_PI_180, v);
   qcalc->mult_quats(q1, zq, q2);
   q1.quat[0] = qcalc->quat[0];
   q1.quat[1] = qcalc->quat[1];
   q1.quat[2] = qcalc->quat[2];
   q1.quat[3] = qcalc->quat[3];
   qcalc->mult_quats(q1, q2, *qcalc);
   qcalc->build_rotmatrix(rotmatrix);
   v[0] = rotmatrix[2][0];
   v[1] = rotmatrix[2][1];
   v[2] = rotmatrix[2][2];
   normalize(v);
   x += (v[0] * speed * speedFactor);
   y += (v[1] * speed * speedFactor);
   z += (v[2] * speed * speedFactor);
}


// Get billboard (face toward) rotation to target point given in local coordinates.
// Return axis and angle for rotation to accomplish billboard.
// rotation[0-2]=axis, rotation[3]=angle
void cSpacial::getBillboard(GLfloat *target, GLfloat *rotation)
{
   GLfloat forward[3];

   // Check for coincidence.
   for (int i = 0; i < 4; i++)
   {
      rotation[i] = 0.0;
   }
   if ((target[0] == 0.0) && (target[1] == 0.0) && (target[2] == 0.0))
   {
      return;
   }

   // Find the rotation from the forward vector to the target.
   getForward(forward);
   getBillboard(target, forward, rotation);
}


// Get billboard rotation from source vector to target vector.
// Return axis and angle for rotation to accomplish billboard.
// rotation[0-2]=axis, rotation[3]=angle
void cSpacial::getBillboard(GLfloat *target, GLfloat *source, GLfloat *rotation)
{
   Vector  v1, v2, v3;
   GLfloat d;

   // Check for invalid condition.
   for (int i = 0; i < 4; i++)
   {
      rotation[i] = 0.0;
   }
   if ((fabs(target[0]) < tol) && (fabs(target[1]) < tol) && (fabs(target[2]) < tol))
   {
      return;
   }
   if ((fabs(source[0]) < tol) && (fabs(source[1]) < tol) && (fabs(source[2]) < tol))
   {
      return;
   }

   // The axis to rotate about is the cross product
   // of the forward vector and the vector to the target.
   v1.x = target[0];
   v1.y = target[1];
   v1.z = target[2];
   v1.Normalize();
   v2.x = source[0];
   v2.y = source[1];
   v2.z = source[2];
   v2.Normalize();
   v3 = v1 ^ v2;
   v3.Normalize();
   rotation[0] = v3.x;
   rotation[1] = v3.y;
   rotation[2] = v3.z;

   // The angle to rotate is the dot product of the vectors.
   d = v1 * v2;
   if (d > 0.9999)
   {
      rotation[3] = 0.0;
   }
   else
   {
      rotation[3] = RadiansToDegrees(acos(d));
   }
}
