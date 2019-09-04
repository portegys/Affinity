//***************************************************************************//
//* File Name: frameRate.hpp                                                *//
//*    Author: Chris McBride chris_a_mcbride@hotmail.com                    *//
//* Date Made: 04/06/02                                                     *//
//* File Desc: Frame rate counter for any project, windows specific.        *//
//* Rev. Date: 11/26/02                                                     *//
//* Rev. Desc: Added frame rate independence and UNIX functionality (TEP)   *//
//*                                                                         *//
//***************************************************************************//
#ifndef __FRAMERATE_HPP__
#define __FRAMERATE_HPP__

#include "gettime.h"

class FrameRate
{
public:

   // Frame rate recalculation frequency (secs).
   enum { FRAME_RECALC_FREQUENCY=2 };

   float        targetFPS;                        // Target frames per second (FPS).
   float        FPS;                              // Current FPS.
   float        speedFactor;                      // Frame rate independence speed factor.
   static float initialSpeedFactor;
   static float maxSpeedFactor;

   FrameRate(float targetFPS)
   {
      this->targetFPS = targetFPS;
      FPS             = targetFPS;
      speedFactor     = initialSpeedFactor;
      frameCount      = 0;
      lastTime        = gettime();
   }


   // Update: call per frame.
   void update()
   {
      TIME currentTime, delta;

      // Count the frame.
      frameCount++;

      // Get the time delta.
      currentTime = gettime();
      delta       = (currentTime - lastTime) / 1000;

      // Time to recalculate frame rate?
      if (delta >= FRAME_RECALC_FREQUENCY)
      {
         // Calculate new values.
         FPS = (float)frameCount / (float)delta;
         if (FPS > 0.0f)
         {
            speedFactor = targetFPS / FPS;
         }
         else
         {
            speedFactor = 0.0f;
         }
         if (speedFactor > maxSpeedFactor) { speedFactor = maxSpeedFactor; }
         frameCount = 0;
         lastTime   = currentTime;
      }
   }


   // Reset.
   void reset()
   {
      FPS         = targetFPS;
      speedFactor = 1.0;
      frameCount  = 0;
      lastTime    = gettime();
   }


private:

   int  frameCount;
   TIME lastTime;
};
#endif
