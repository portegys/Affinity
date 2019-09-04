//***************************************************************************//
//* File Name: spacial.hpp                                                  *//
//* Author:    Tom Portegys, portegys@ilstu.edu                             *//
//* Date Made: 07/25/02                                                     *//
//* File Desc: Class declaration representing spacial properties:           *//
//*             rotation, translation, scale, and speed.                    *//
//* Rev. Date:                                                              *//
//* Rev. Desc:                                                              *//
//*                                                                         *//
//***************************************************************************//

#ifndef __SPACIAL_HPP__
#define __SPACIAL_HPP__

#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <assert.h>
#include "quaternion.hpp"
#include "matrix.hpp"
#include "vector.hpp"

#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD    std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<GLfloat>   Matrix;
#else
typedef matrix            Matrix;
#endif

// Pi and radians
#define M_PI          3.14159265358979323846
#define DIV_PI_180    .01745329251
#define DIV_180_PI    57.29577951

class cSpacial
{
public:

   // Rotation.
   GLfloat pitch, yaw, roll;                      // Rotation rates.
   GLfloat rotmatrix[4][4];
   class cQuaternion * qcalc;

   // Translation.
   GLfloat x, y, z;

   // Scale.
   GLfloat scale;

   // Speed.
   GLfloat speed;
   GLfloat speedFactor;

   // Constructors.
   cSpacial()
   {
      qcalc = new cQuaternion();
      assert(qcalc != NULL);
      clear();
   }


   cSpacial(GLfloat pitch, GLfloat yaw, GLfloat roll,
            GLfloat x, GLfloat y, GLfloat z, GLfloat scale, GLfloat speed)
   {
      qcalc = new cQuaternion();
      assert(qcalc != NULL);
      initialize(pitch, yaw, roll, x, y, z, scale, speed);
   }


   void clear()
   {
      initialize(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);
   }


   void initialize(GLfloat pitch, GLfloat yaw, GLfloat roll,
                   GLfloat x, GLfloat y, GLfloat z, GLfloat scale, GLfloat speed,
                   cQuaternion *qcalc = NULL)
   {
      if (qcalc == NULL)
      {
         this->qcalc->clear();
      }
      else
      {
         for (int i = 0; i < 4; i++) { this->qcalc->quat[i] = qcalc->quat[i]; }
      }
      this->qcalc->build_rotmatrix(rotmatrix);
      this->pitch = pitch;
      this->yaw   = yaw;
      this->roll  = roll;
      this->x     = x;
      this->y     = y;
      this->z     = z;
      this->scale = scale;
      this->speed = speed;
      speedFactor = 1.0f;
   }


   // Destructor.
   ~cSpacial() { delete qcalc; }

   // Update rotation and translation state.
   void update();

   // Get direction vectors.
   void getRight(GLfloat *v)
   {
      v[0] = rotmatrix[0][0];
      v[1] = rotmatrix[0][1];
      v[2] = rotmatrix[0][2];
   }


   void getUp(GLfloat *v)
   {
      v[0] = rotmatrix[1][0];
      v[1] = rotmatrix[1][1];
      v[2] = rotmatrix[1][2];
   }


   void getForward(GLfloat *v)
   {
      v[0] = rotmatrix[2][0];
      v[1] = rotmatrix[2][1];
      v[2] = rotmatrix[2][2];
   }


   // Get model transformation matrix.
   void getModelTransform(GLfloat *matrix)
   {
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glTranslatef(x, y, z);
      qcalc->build_rotmatrix(rotmatrix);
      glMultMatrixf(&rotmatrix[0][0]);
      glScalef(scale, scale, scale);
      glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
      glPopMatrix();
   }


   // Get world coordinates from local.
   void localToWorld(GLfloat *local, GLfloat *world)
   {
      int     i, j;
      GLfloat m[16];
      Matrix  x(4, 4), p(4, 1), t(4, 1);

      getModelTransform(m);
      for (i = 0; i < 4; i++)
      {
         for (j = 0; j < 4; j++)
         {
            x(i, j) = m[(j * 4) + i];
         }
      }
      p(0, 0)  = local[0];
      p(1, 0)  = local[1];
      p(2, 0)  = local[2];
      p(3, 0)  = 1.0;
      t        = x * p;
      world[0] = t(0, 0);
      world[1] = t(1, 0);
      world[2] = t(2, 0);
   }


   // Transform local point.
   void transformPoint(GLfloat *point)
   {
      localToWorld(point, point);
   }


   // Inverse transform local point.
   void inverseTransformPoint(GLfloat *point)
   {
      int     i, j;
      GLfloat m[16];
      Matrix  x(4, 4), y(4, 4), p(4, 1), t(4, 1);

      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();
      getModelTransform(m);
      for (i = 0; i < 4; i++)
      {
         for (j = 0; j < 4; j++)
         {
            x(i, j) = m[(j * 4) + i];
         }
      }
      y        = !x;
      p(0, 0)  = point[0];
      p(1, 0)  = point[1];
      p(2, 0)  = point[2];
      p(3, 0)  = 1.0;
      t        = y * p;
      point[0] = t(0, 0);
      point[1] = t(1, 0);
      point[2] = t(2, 0);
      glPopMatrix();
   }


   // Normalize vector.
   inline static void normalize(GLfloat *v)
   {
      GLfloat d = (sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2])));

      v[0] = v[0] / d;
      v[1] = v[1] / d;
      v[2] = v[2] / d;
   }


   // Get billboard rotation to given target point.
   void getBillboard(GLfloat *target, GLfloat *rotation);

   // Get billboard rotation from source to target vectors.
   void getBillboard(GLfloat *target, GLfloat *source, GLfloat *rotation);

   // Load an axis-angle rotation into quaternion.
   void loadRotation(GLfloat angle, GLfloat *axis)
   {
      qcalc->loadRotation(DegreesToRadians(angle), axis);
      build_rotmatrix();
   }


   // Merge an axis-angle rotation into quaternion.
   void mergeRotation(GLfloat angle, GLfloat *axis)
   {
      qcalc->mergeRotation(DegreesToRadians(angle), axis);
      build_rotmatrix();
   }


   // Build rotation matrix from quaternion.
   void build_rotmatrix()
   {
      qcalc->build_rotmatrix(rotmatrix);
   }


   // Point-to-point distance.
   static GLfloat pointDistance(GLfloat *p1, GLfloat *p2)
   {
      GLfloat dx = p1[0] - p2[0];
      dx *= dx;

      GLfloat dy = p1[1] - p2[1];
      dy *= dy;
      GLfloat dz = p1[2] - p2[2];
      dz *= dz;
      return(sqrt(dx + dy + dz));
   }
};
#endif
