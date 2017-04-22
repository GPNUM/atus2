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

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <cmath>
#include <omp.h>
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"

using namespace std;

template <class T> bool Read_Data( const char *filename, generic_header *header, T *&field )
{
  FILE *fh = fopen( filename, "r" );

  if ( fh == nullptr ) return false;

  fseek( fh, sizeof(generic_header), SEEK_SET );
  size_t Nges = header->nDimX * header->nDimY * header->nDimZ;

  fread( (void *)field, sizeof(T), Nges, fh );
  fclose(fh);
  return true;
}

int main(int argc, char *argv[])
{
  FILE *fh = nullptr;

  int no_of_threads = 4;
  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );
  fftw_init_threads();
  fftw_plan_with_nthreads( no_of_threads );
  omp_set_num_threads( no_of_threads );

  generic_header header = {};

  if ( argc > 1 )
  {
    fh = fopen( argv[1], "r" );
    if ( fh == nullptr )
    {
      printf( "Could not open file %s.\n", argv[1] );
      exit(0);
    }

    fread( &header, sizeof(generic_header), 1, fh );

    printf( "### %s\n", argv[1] );
    printf( "# nDims    == %lld\n", header.nDims );
    printf( "# nDimX    == %lld\n", header.nDimX );
    printf( "# nDimY    == %lld\n", header.nDimY );
    printf( "# nDimZ    == %lld\n", header.nDimZ );
    printf( "# nDatatyp == %lld\n", header.nDatatyp );
    printf( "# bAtom    == %d\n", header.bAtom );
    printf( "# bComplex == %d\n", header.bComplex );
    printf( "# t        == %g\n", header.t );
    printf( "# dt       == %g\n", header.dt );
    printf( "# xMin     == %g\n", header.xMin );
    printf( "# xMax     == %g\n", header.xMax );
    printf( "# yMin     == %g\n", header.yMin );
    printf( "# yMax     == %g\n", header.yMax );
    printf( "# zMin     == %g\n", header.zMin );
    printf( "# zMax     == %g\n", header.zMax );
    printf( "# dx       == %g\n", header.dx );
    printf( "# dy       == %g\n", header.dy );
    printf( "# dz       == %g\n", header.dz );
    printf( "# dkx      == %g\n", header.dkx );
    printf( "# dky      == %g\n", header.dky );
    printf( "# dkz      == %g\n", header.dkz );
    fclose(fh);
  }

  if ( header.nDims == 2 ) header.nDimZ = 1;
  long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  char filename[255] = {};

  // 2D Stuff
  if ( header.nDims == 2 &&  header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    double *res_x = fftw_alloc_real( header.nDimX );
    double *res_y = fftw_alloc_real( header.nDimY );
    Fourier::cft_2d *ft1 = new Fourier::cft_2d(header);
    ft1->SetFix(true);

    fftw_complex *in = ft1->Getp2In();
    fftw_complex *out = ft1->Getp2Out();

    Read_Data<fftw_complex>( argv[1], &header, in );
    memset( reinterpret_cast<void *>(res_x), '0', sizeof(double)*header.nDimX );
    memset( reinterpret_cast<void *>(res_y), '0', sizeof(double)*header.nDimY );

    ft1->ft(-1);

    FILE *fh2 = fopen( "ana_mom.bin", "w" );
    fwrite( (void *)&header, sizeof(generic_header), 1, fh2 );
    fwrite( (void *)out, sizeof(fftw_complex), Nges, fh2 );
    fclose(fh2);


    for ( int j=0; j < header.nDimY; j++ )
    {
      for ( int i=0; i<header.nDimX; i++ )
      {
        int ij = j+header.nDimY*i;
        res_x[i] += (out[ij][0]*out[ij][0] + out[ij][1]*out[ij][1]);
      }
    }

    const int shift_x = header.nDimX/2;
    fh = fopen( "ana_mom_x.txt", "w" );
    for ( int l=0; l<header.nDimX; l++ )
    {
      fprintf( fh, "%g\t%g\n", double(l-shift_x)*header.dkx, res_x[l] / double(header.nDimY));
    }
    fclose(fh);

    for ( int i=0; i < header.nDimX; i++ )
    {
      for ( int j=0; j<header.nDimY; j++ )
      {
        int ij = j+header.nDimY*i;
        res_y[j] += (out[ij][0]*out[ij][0] + out[ij][1]*out[ij][1]);
      }
    }

    const int shift_y = header.nDimY/2;
    fh = fopen( "ana_mom_y.txt", "w" );
    for ( int l=0; l<header.nDimY; l++ )
    {
      fprintf( fh, "%g\t%g\n", double(l-shift_y)*header.dky, res_y[l] / double(header.nDimX));
    }
    fclose(fh);

    delete ft1;

    fftw_free(res_x);
    fftw_free(res_y);
  }

  if ( header.nDims == 1 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    double *res = fftw_alloc_real( header.nDimX );
    Fourier::cft_1d *ft1 = new Fourier::cft_1d( header );

    fftw_complex *in = ft1->Getp2In();
    fftw_complex *out = ft1->Getp2Out();

    Read_Data<fftw_complex>( argv[1], &header, in );
    bzero( (void *)res, sizeof(double)*header.nDimX);

    ft1->ft(-1);

    for ( int i=0; i<header.nDimX; i++ )
    {
      res[i] += (out[i][0]*out[i][0] + out[i][1]*out[i][1]);
    }


    const int shift_x = header.nDimX/2;
    fh = fopen( "ana_mom.txt", "w" );
    for ( int l=0; l<header.nDimX; l++ )
    {
      fprintf( fh, "%g\t%g\n", double(l-shift_x)*header.dkx, res[l] / double(header.nDimY));
    }
    fclose(fh);

    delete ft1;
    fftw_free(res);
  }

  // 2D Stuff
  if ( header.nDims == 2 && header.bComplex == 0 && header.nDatatyp == sizeof(double) )
  {
  }

  // 3D Stuff
  if ( header.nDims == 3 &&  header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    double *res_x = fftw_alloc_real( header.nDimX );
    double *res_y = fftw_alloc_real( header.nDimY );
    double *res_z = fftw_alloc_real( header.nDimZ );
    Fourier::cft_3d *ft1 = new Fourier::cft_3d( header );
    ft1->SetFix(true);

    fftw_complex *in = ft1->Getp2In();
    fftw_complex *out = ft1->Getp2Out();

    Read_Data<fftw_complex>( argv[1], &header, in );
    bzero( (void *)res_x, sizeof(double)*header.nDimX);
    bzero( (void *)res_y, sizeof(double)*header.nDimY);
    bzero( (void *)res_z, sizeof(double)*header.nDimZ);

    ft1->ft(-1);

    FILE *fh2 = fopen( "ana_mom.bin", "w" );
    fwrite( (void *)&header, sizeof(generic_header), 1, fh2 );
    fwrite( (void *)out, sizeof(fftw_complex), Nges, fh2 );
    fclose(fh2);

    for ( int k=0; k<header.nDimZ; k++ )
    {
      for ( int j=0; j<header.nDimY; j++ )
      {
        for ( int i=0; i<header.nDimX; i++ )
        {
          int ijk = k+header.nDimZ*(j+header.nDimY*i);
          res_x[i] += (out[ijk][0]*out[ijk][0] + out[ijk][1]*out[ijk][1]);
        }
      }
    }

    const int shift_x = header.nDimX/2;
    fh = fopen( "ana_mom_x.txt", "w" );
    for ( int l=0; l<header.nDimX; l++ )
    {
      fprintf( fh, "%g\t%g\n", double(l-shift_x)*header.dkx, res_x[l] / double(header.nDimY));
    }
    fclose(fh);

    for ( int k=0; k<header.nDimZ; k++ )
    {
      for ( int i=0; i<header.nDimX; i++ )
      {
        for ( int j=0; j<header.nDimY; j++ )
        {
          int ijk = k+header.nDimZ*(j+header.nDimY*i);
          res_y[j] += (out[ijk][0]*out[ijk][0] + out[ijk][1]*out[ijk][1]);
        }
      }
    }


    const int shift_y = header.nDimY/2;
    fh = fopen( "ana_mom_y.txt", "w" );
    for ( int l=0; l<header.nDimY; l++ )
    {
      fprintf( fh, "%g\t%g\n", double(l-shift_y)*header.dky, res_y[l] / double(header.nDimX));
    }
    fclose(fh);

    delete ft1;

    fftw_free(res_x);
    fftw_free(res_y);
  }

  if ( header.nDims == 3 && header.bComplex == 0 && header.nDatatyp == sizeof(double) )
  {
  }

  fftw_cleanup_threads();
  return EXIT_SUCCESS;
}
