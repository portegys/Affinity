//***************************************************************************//
//* File Name: baseObject.hpp                                               *//
//* Author:    Tom Portegys, portegys@ilstu.edu                             *//
//* Date Made: 07/25/02                                                     *//
//* File Desc: Class declaration representing a basic object.               *//
//* Rev. Date:                                                              *//
//* Rev. Desc:                                                              *//
//*                                                                         *//
//***************************************************************************//

#ifndef __BASE_OBJECT_HPP__
#define __BASE_OBJECT_HPP__

#ifdef WIN32
#include <windows.h>
#endif
#include <GL/glut.h>
#include "spacial.hpp"
#include "fileio.h"

class BaseObject
{
public:

   // Constructor.
   BaseObject()
   {
      m_spacial = new cSpacial();
      assert(m_spacial != NULL);
      m_pitch = m_yaw = m_roll = 0.0f;
      m_speed = 0.0f;
   }


   // Destructor.
   ~BaseObject()
   {
      delete m_spacial;
   }


   // Rotations.
   GLfloat getPitch() { return(m_pitch); }
   GLfloat getYaw()   { return(m_yaw); }
   GLfloat getRoll()  { return(m_roll); }
   void setPitch(GLfloat n)
   {
      m_spacial->pitch = (n - m_pitch);
      m_spacial->update();
      m_spacial->pitch = 0.0f;
      m_pitch          = n;
   }


   void setYaw(GLfloat n)
   {
      m_spacial->yaw = (n - m_yaw);
      m_spacial->update();
      m_spacial->yaw = 0.0f;
      m_yaw          = n;
   }


   void setRoll(GLfloat n)
   {
      m_spacial->roll = (n - m_roll);
      m_spacial->update();
      m_spacial->roll = 0.0f;
      m_roll          = n;
   }


   void addPitch(GLfloat n)
   {
      m_spacial->pitch = n;
      m_spacial->update();
      m_spacial->pitch = 0.0f;
      m_pitch          = n + m_pitch;
   }


   void addYaw(GLfloat n)
   {
      m_spacial->yaw = n;
      m_spacial->update();
      m_spacial->yaw = 0.0f;
      m_yaw          = n + m_yaw;
   }


   void addRoll(GLfloat n)
   {
      m_spacial->roll = n;
      m_spacial->update();
      m_spacial->roll = 0.0f;
      m_roll          = n + m_roll;
   }


   // Get direction vectors.
   void getRight(GLfloat *v)   { m_spacial->getRight(v); }
   void getUp(GLfloat *v)      { m_spacial->getUp(v); }
   void getForward(GLfloat *v) { m_spacial->getForward(v); }

   // Position.
   void getPosition(GLfloat *v)
   {
      v[0] = m_spacial->x;
      v[1] = m_spacial->y;
      v[2] = m_spacial->z;
   }


   void setPosition(GLfloat *v)
   {
      m_spacial->x = v[0];
      m_spacial->y = v[1];
      m_spacial->z = v[2];
      m_spacial->update();
   }


   // Scale.
   GLfloat getScale() { return(m_spacial->scale); }
   void setScale(GLfloat n) { m_spacial->scale = n; }

   // Speed.
   GLfloat getSpeed() { return(m_speed); }
   void setSpeed(GLfloat speed)
   {
      m_speed = speed;
   }


   void addSpeed(GLfloat speed)
   {
      m_speed += speed;
      if (m_speed < 0.0) { m_speed = 0.0f; }
   }


   GLfloat getSpeedFactor() { return(m_spacial->speedFactor); }
   void setSpeedFactor(GLfloat n)
   {
      m_spacial->speedFactor = n;
      if (m_spacial->speedFactor < 0.0) { m_spacial->speedFactor = 0.0; }
   }


   // Set spacial state.
   cSpacial *getSpacial() { return(m_spacial); }

   // Set spacial state.
   void setSpacial(GLfloat pitch, GLfloat yaw, GLfloat roll,
                   GLfloat speed, cSpacial *spacial)
   {
      m_spacial->initialize(0.0f, 0.0f, 0.0f,
                            spacial->x, spacial->y, spacial->z,
                            spacial->scale, 0.0f, spacial->qcalc);
      m_pitch = pitch;
      m_yaw   = yaw;
      m_roll  = roll;
      m_speed = speed;
   }


   // Reset spacial state.
   void resetSpacial()
   {
      m_spacial->clear();
      m_pitch = m_yaw = m_roll = 0.0f;
      m_speed = 0.0f;
   }


   // Update.
   void update()
   {
      m_spacial->speed = m_speed;
      m_spacial->update();
      m_spacial->speed = 0.0;
   }


   // Get model transformation matrix.
   void getModelTransform(GLfloat *matrix)
   {
      m_spacial->getModelTransform(matrix);
   }


   // Get world coordinates from local.
   void localToWorld(GLfloat *local, GLfloat *world)
   {
      m_spacial->localToWorld(local, world);
   }


   // Transform local point.
   void transformPoint(GLfloat *point)
   {
      m_spacial->transformPoint(point);
   }


   // Inverse transform local point.
   void inverseTransformPoint(GLfloat *point)
   {
      m_spacial->inverseTransformPoint(point);
   }


   // Get billboard (face toward) rotation to target point given in local coordinates.
   // Return axis and angle for rotation to accomplish billboard.
   // rotation[0-2]=axis, rotation[3]=angle
   void getBillboard(GLfloat *target, GLfloat *rotation)
   {
      m_spacial->getBillboard(target, rotation);
   }


   // Get billboard (face toward) rotation from source to target vectors.
   // Return axis and angle for rotation to accomplish billboard.
   // rotation[0-2]=axis, rotation[3]=angle
   void getBillboard(GLfloat *target, GLfloat *source, GLfloat *rotation)
   {
      m_spacial->getBillboard(target, source, rotation);
   }


   // Load an axis-angle rotation into quaternion.
   void loadRotation(GLfloat angle, GLfloat *axis)
   {
      m_spacial->loadRotation(angle, axis);
   }


   // Merge an axis-angle rotation into quaternion.
   void mergeRotation(GLfloat angle, GLfloat *axis)
   {
      m_spacial->mergeRotation(angle, axis);
   }


   // Draw axes.
   void drawAxes(GLfloat span = 1.0f);

   // Set block vertices.
   static void setBlockVertices(Vector *vertices, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);

   // Draw a block.
   static void drawBlock(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);

   // Load object.
   void load(char *filename);
   void load(FILE *fp);

   // Save object.
   void save(char *filename);
   void save(FILE *fp);

protected:

   // Spacial properties.
   cSpacial *m_spacial;
   GLfloat  m_pitch, m_yaw, m_roll;
   GLfloat  m_speed;
};
#endif
