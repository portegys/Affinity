/*
 * Get time in milliseconds since the initial call.
 */

#ifndef __GETTIME__
#define __GETTIME__

#ifdef UNIX
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#else
#include <windows.h>
#endif
#include <assert.h>

typedef unsigned long   TIME;
#define INVALID_TIME    0xffffffffUL
TIME gettime();
#endif
