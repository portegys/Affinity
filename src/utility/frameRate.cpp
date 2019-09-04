//***************************************************************************//
//* File Name: frameRate.cpp                                                *//
//*    Author: Chris McBride chris_a_mcbride@hotmail.com                    *//
//* Date Made: 04/06/02                                                     *//
//* File Desc: Frame rate counter for any project, windows specific.        *//
//* Rev. Date: 11/26/02                                                     *//
//* Rev. Desc: Added frame rate independence and UNIX functionality (TEP)   *//
//*                                                                         *//
//***************************************************************************//

#include "frameRate.hpp"

// Initial speed factor.
float FrameRate::initialSpeedFactor = 1.0f;

// Maximum speed factor.
float FrameRate::maxSpeedFactor = 5.0f;
