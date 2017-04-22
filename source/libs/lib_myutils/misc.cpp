//
// ATUS2 - The ATUS2 package is atom interferometer Toolbox developed at ZARM
// (CENTER OF APPLIED SPACE TECHNOLOGY AND MICROGRAVITY), Germany. This project is
// founded by the DLR Agentur (Deutsche Luft und Raumfahrt Agentur). Grant numbers:
// 50WM0942, 50WM1042, 50WM1342.
// Copyright (C) 2017 Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
//
// This file is part of ATUS2.
//
// ATUS2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ATUS2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ATUS2.  If not, see <http://www.gnu.org/licenses/>.
//


#include <string>
#include <limits>
#include <cmath>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

unsigned getbits( unsigned x, unsigned p, unsigned n)
{
  return ((x >> (p+1-n)) & ~(~0 << n));
}

double sign( double val )
{
  return copysign(1,val);
}

double Heaviside( double x )
{
  return 0.5*(1.0+sign(x));
}

double rect( double x, double xl, double xr )
{
  return Heaviside(x-xl)-Heaviside(x-xr);
}

float FIX_FLOAT(float arg)
{
  float retval;
  char *mem  = reinterpret_cast<char *>(&arg);
  char *mem2 = reinterpret_cast<char *>(&retval);

  mem2[0] = mem[3];
  mem2[1] = mem[2];
  mem2[2] = mem[1];
  mem2[3] = mem[0];
  return retval;
}

float FIX_INT(int arg)
{
  int retval;
  char *mem  = reinterpret_cast<char *>(&arg);
  char *mem2 = reinterpret_cast<char *>(&retval);

  mem2[0] = mem[3];
  mem2[1] = mem[2];
  mem2[2] = mem[1];
  mem2[3] = mem[0];
  return retval;
}

unsigned short FIX_SHORT(unsigned short arg)
{
  unsigned short retval;
  char *mem  = reinterpret_cast<char *>(&arg);
  char *mem2 = reinterpret_cast<char *>(&retval);

  mem2[0] = mem[1];
  mem2[1] = mem[0];
  return retval;
}

void Find_Oldest_File_in_Dir( string path, string pattern, string &retval, int &nosubfolder )
{
  DIR *pDir;
  struct dirent *pDirEnt;
  struct stat fileinfo;
  time_t reftime;

  double timedif_old = std::numeric_limits<double>::max();
  double timedif;
  string filename;

  nosubfolder = 0;

  time(&reftime);
  printf("reftime: %s", ctime(&reftime));

  retval.empty();
  pDir = opendir( path.c_str() );
  pDirEnt = readdir( pDir );
  while ( pDirEnt != nullptr )
  {
    filename = pDirEnt->d_name;
    if ( stat( pDirEnt->d_name, &fileinfo ) == -1 )
    {
      pDirEnt = readdir( pDir );
      continue;
    }

    /*
        printf("filename: %s - ", pDirEnt->d_name);

        switch (fileinfo.st_mode & S_IFMT)
        {
          case S_IFBLK:  printf("block device\n");            break;
          case S_IFCHR:  printf("character device\n");        break;
          case S_IFDIR:  printf("directory\n");               break;
          case S_IFIFO:  printf("FIFO/pipe\n");               break;
          case S_IFLNK:  printf("symlink\n");                 break;
          case S_IFREG:  printf("regular file\n");            break;
          case S_IFSOCK: printf("socket\n");                  break;
          default:       printf("unknown?\n");                break;
        }

        printf(" device: %lld\n",                       fileinfo.st_dev);
        printf(" inode: %ld\n",                         fileinfo.st_ino);
        printf(" protection: %o\n",                     fileinfo.st_mode);
        printf(" number of hard links: %d\n",           fileinfo.st_nlink);
        printf(" user ID of owner: %d\n",               fileinfo.st_uid);
        printf(" group ID of owner: %d\n",              fileinfo.st_gid);
        printf(" device type (if inode device): %lld\n",fileinfo.st_rdev);
        printf(" total size, in bytes: %ld\n",          fileinfo.st_size);
        printf(" blocksize for filesystem I/O: %ld\n",  fileinfo.st_blksize);
        printf(" number of blocks allocated: %ld\n",    fileinfo.st_blocks);
        printf(" time of last access: %ld : %s",        fileinfo.st_atime, ctime(&fileinfo.st_atime));
        printf(" time of last modification: %ld : %s",  fileinfo.st_mtime, ctime(&fileinfo.st_mtime));
        printf(" time of last change: %ld : %s",        fileinfo.st_ctime, ctime(&fileinfo.st_ctime));
    */

    if ( (fileinfo.st_mode & S_IFMT) == S_IFDIR )
    {
      nosubfolder++;
    }
    if ( (fileinfo.st_mode & S_IFMT) != S_IFREG )
    {
      pDirEnt = readdir( pDir );
      continue;
    }

    if ( filename.size() > pattern.size() )
    {
      if ( filename.compare( filename.size()-pattern.size(), pattern.size(), pattern ) == 0 )
      {
        timedif = difftime(reftime,fileinfo.st_mtime);
        if ( fabs(timedif) < fabs(timedif_old) )
        {
          retval = filename;
          timedif_old = timedif;
        }

        //printf("Last status change:       %s", ctime(&fileinfo.st_ctime));
        //printf("Last file access:         %s", ctime(&fileinfo.st_atime));
        //printf("Last file modification:   %s", ctime(&fileinfo.st_mtime));
      }
    }
    pDirEnt = readdir( pDir );
  }

  nosubfolder = nosubfolder-2;

  closedir( pDir );
}

#if __APPLE__ && __MACH__
void sincos( double x, double *sinus, double *kosinus )
{
  *sinus = sin(x);
  *kosinus = cos(x);
}

void sincosf( float x, float *sinus, float *kosinus )
{
  *sinus = sinf(x);
  *kosinus = cosf(x);
}
#endif
