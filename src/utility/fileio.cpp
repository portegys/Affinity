// Load and save in binary or readable format.

#include "fileio.h"

// Get path to resource file.
// Application responsible for freeing memory.
char *getResourcePath(char *file)
{
   return(getPath((char *)"resource", file));
}


// Get path to data file.
// Application responsible for freeing memory.
char *getDataPath(char *file)
{
   return(getPath((char *)"data", file));
}


// Get path to file within directory.
// Application responsible for freeing memory.
char *getPath(char *dir, char *file)
{
   char *home, *path;

   // Sanity check.
   if ((file == NULL) || (*file == '\0'))
   {
      return(NULL);
   }

   // Fixed path?
   if ((file[0] == '/') || (file[0] == '.'))
   {
      if ((path = (char *)malloc(strlen(file) + 1)) == NULL)
      {
         fprintf(stderr, "getPath: cannot malloc path memory\n");
         exit(1);
      }
      strcpy(path, file);
      return(path);
   }

   // Check AFFINITY_HOME directory path environment variable.
   if ((home = getenv("AFFINITY_HOME")) != NULL)
   {
      if ((path = (char *)malloc(strlen(home) + strlen(dir) + 2)) == NULL)
      {
         fprintf(stderr, "getPath: cannot malloc path memory\n");
         exit(1);
      }
      sprintf(path, "%s/%s", home, dir);
#ifdef UNIX
      if (access(path, F_OK) != -1)
#else
      if (GetFileAttributes(path) != 0xffffffff)
#endif
      {
         // Add the file.
         free(path);
         if ((path = (char *)malloc(strlen(home) +
                                    strlen(dir) + strlen(file) + 3)) == NULL)
         {
            fprintf(stderr, "getPath: cannot malloc path memory\n");
            exit(1);
         }
         sprintf(path, "%s/%s/%s", home, dir, file);
         return(path);
      }
      else
      {
         free(path);
      }
   }

   // Try relative paths.
   if ((path = (char *)malloc(strlen(dir) + 4)) == NULL)
   {
      fprintf(stderr, "getPath: cannot malloc path memory\n");
      exit(1);
   }
   sprintf(path, "../%s", dir);
#ifdef UNIX
   if (access(path, F_OK) != -1)
#else
   if (GetFileAttributes(path) != 0xffffffff)
#endif
   {
      free(path);
      if ((path = (char *)malloc(strlen(dir) + strlen(file) + 5)) == NULL)
      {
         fprintf(stderr, "getPath: cannot malloc path memory\n");
         exit(1);
      }
      sprintf(path, "../%s/%s", dir, file);
      return(path);
   }
   else
   {
      free(path);
   }
   if ((path = (char *)malloc(strlen(dir) + 7)) == NULL)
   {
      fprintf(stderr, "getPath: cannot malloc path memory\n");
      exit(1);
   }
   sprintf(path, "../../%s", dir);
#ifdef UNIX
   if (access(path, F_OK) != -1)
#else
   if (GetFileAttributes(path) != 0xffffffff)
#endif
   {
      free(path);
      if ((path = (char *)malloc(strlen(dir) + strlen(file) + 8)) == NULL)
      {
         fprintf(stderr, "getPath: cannot malloc path memory\n");
         exit(1);
      }
      sprintf(path, "../../%s/%s", dir, file);
      return(path);
   }
   else
   {
      free(path);
   }

   // Default to input file.
   if ((path = (char *)malloc(strlen(file) + 1)) == NULL)
   {
      fprintf(stderr, "getPath: cannot malloc path memory\n");
      exit(1);
   }
   strcpy(path, file);
   return(path);
}


void myfreadInt(int *i, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fread(i, sizeof(int), 1, fp) == 1);
#else
   fscanf(fp, "%d", i);
#endif
}


void myfwriteInt(int *i, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fwrite(i, sizeof(int), 1, fp) == 1);
#else
   fprintf(fp, "%d\n", *i);
#endif
}


void myfreadLong(unsigned long *l, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fread(l, sizeof(unsigned long), 1, fp) == 1);
#else
   fscanf(fp, "%d", l);
#endif
}


void myfwriteLong(unsigned long *l, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fwrite(l, sizeof(unsigned long), 1, fp) == 1);
#else
   fprintf(fp, "%d\n", *l);
#endif
}


void myfreadFloat(float *f, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fread(f, sizeof(float), 1, fp) == 1);
#else
   char buf[100];
   fscanf(fp, "%s", buf);
   *f = (float)atof(buf);
#endif
}


void myfwriteFloat(float *f, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fwrite(f, sizeof(float), 1, fp) == 1);
#else
   fprintf(fp, "%f\n", *f);
#endif
}


void myfreadDouble(double *d, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fread(d, sizeof(double), 1, fp) == 1);
#else
   char buf[100];
   fscanf(fp, "%s", buf);
   *d = (double)atof(buf);
#endif
}


void myfwriteDouble(double *d, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fwrite(d, sizeof(double), 1, fp) == 1);
#else
   fprintf(fp, "%f\n", *d);
#endif
}


void myfreadBool(bool *b, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fread(b, sizeof(bool), 1, fp) == 1);
#else
   int v;
   fscanf(fp, "%d", &v);
   if (v == 1)
   {
      *b = true;
   }
   else
   {
      *b = false;
   }
#endif
}


void myfwriteBool(bool *b, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fwrite(b, sizeof(bool), 1, fp) == 1);
#else
   if (*b)
   {
      fprintf(fp, "1\n");
   }
   else
   {
      fprintf(fp, "0\n");
   }
#endif
}


void myfreadChar(unsigned char *c, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fread(c, sizeof(unsigned char), 1, fp) == 1);
#else
   fscanf(fp, "%c", c);
#endif
}


void myfwriteChar(unsigned char *c, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fwrite(c, sizeof(unsigned char), 1, fp) == 1);
#else
   fprintf(fp, "%c\n", *c);
#endif
}


void myfreadBytes(unsigned char *bytes, int size, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fread(bytes, size, 1, fp) == 1);
#else
   int           len  = (2 * size) + 1;
   unsigned char *buf = new unsigned char[len];
   assert(buf != NULL);
   fscanf(fp, "%s", buf);
   int i, j, d1, d2;
   for (i = 0; i < size; i++)
   {
      j = 2 * i;
      if ((buf[j] >= '0') && (buf[j] <= '9'))
      {
         d1 = buf[j] - '0';
      }
      else
      {
         d1 = buf[j] - 'a' + 10;
      }
      j++;
      if ((buf[j] >= '0') && (buf[j] <= '9'))
      {
         d2 = buf[j] - '0';
      }
      else
      {
         d2 = buf[j] - 'a' + 10;
      }
      bytes[i] = (d1 * 16) + d2;
   }
   delete [] buf;
#endif
}


void myfwriteBytes(unsigned char *bytes, int size, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fwrite(bytes, size, 1, fp) == 1);
#else
   int  len  = (2 * size) + 1;
   char *buf = new char[len];
   assert(buf != NULL);
   for (int i = 0; i < size; i++)
   {
      sprintf(&buf[2 * i], "%02x", bytes[i]);
   }
   buf[len - 1] = '\0';
   fprintf(fp, "%s\n", buf);
   delete [] buf;
#endif
}


void myfreadString(char *str, int size, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fread(str, size, 1, fp) == 1);
#else
   char *buf = new char[size + 1];
   assert(buf != NULL);
   fscanf(fp, "%s", buf);
   for (int i = 0; i < size; i++)
   {
      if (buf[i] == '`')
      {
         buf[i] = ' ';
      }
   }
   strncpy(str, buf, size);
   str[size - 1] = '\0';
   delete buf;
#endif
}


void myfwriteString(char *str, int size, FILE *fp)
{
#ifdef BINARY_FILE_FORMAT
   assert(fwrite(str, size, 1, fp) == 1);
#else
   size = (int)strlen(str) + 1;
   char *buf = new char[size];
   assert(buf != NULL);
   for (int i = 0; i < size; i++)
   {
      buf[i] = str[i];
      if (str[i] == ' ')
      {
         buf[i] = '`';
      }
   }
   fprintf(fp, "%s\n", buf);
   delete buf;
#endif
}
