// Load and store functions.

#ifndef __FILEIO__
#define __FILEIO__

#ifdef WIN32
#include <windows.h>
#endif
#ifdef UNIX
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// Find paths to files.
char *getResourcePath(char *file);
char *getDataPath(char *file);
char *getPath(char *dir, char *file);

// File I/O.
#define FREAD_INT        myfreadInt
#define FWRITE_INT       myfwriteInt
#define FREAD_LONG       myfreadLong
#define FWRITE_LONG      myfwriteLong
#define FREAD_FLOAT      myfreadFloat
#define FWRITE_FLOAT     myfwriteFloat
#define FREAD_DOUBLE     myfreadDouble
#define FWRITE_DOUBLE    myfwriteDouble
#define FREAD_BOOL       myfreadBool
#define FWRITE_BOOL      myfwriteBool
#define FREAD_CHAR       myfreadChar
#define FWRITE_CHAR      myfwriteChar
#define FREAD_BYTES      myfreadBytes
#define FWRITE_BYTES     myfwriteBytes
#define FREAD_STRING     myfreadString
#define FWRITE_STRING    myfwriteString
void myfreadInt(int *, FILE *);
void myfwriteInt(int *, FILE *);
void myfreadLong(unsigned long *, FILE *);
void myfwriteLong(unsigned long *, FILE *);
void myfreadFloat(float *, FILE *);
void myfwriteFloat(float *, FILE *);
void myfreadDouble(double *, FILE *);
void myfwriteDouble(double *, FILE *);
void myfreadBool(bool *, FILE *);
void myfwriteBool(bool *, FILE *);
void myfreadChar(unsigned char *, FILE *);
void myfwriteChar(unsigned char *, FILE *);
void myfreadBytes(unsigned char *, int size, FILE *);
void myfwriteBytes(unsigned char *, int size, FILE *);
void myfreadString(char *, int size, FILE *);
void myfwriteString(char *, int size, FILE *);
#endif
