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

/******************************************************************************
Extrahiert Schnitte aus 3D Daten
*******************************************************************************
Aufrufe:
    slice_3d Dateiname x y z
    slice_3d Dateiname x x 512
    slice_3d Dateiname x 512 x
    slice_3d Dateiname 512 x x
******************************************************************************/

#include "fftw3.h"
#include "my_structs.h"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <omp.h>

void extract_x_const( fftw_complex *Phi, fftw_complex *slice, generic_header *header, const int x0 )
{
  const long long dimY = header->nDimY;
  const long long dimZ = header->nDimZ;

  #pragma omp parallel for collapse(2)
  for ( long long j=0; j<dimY; j++ )
  {
    for ( long long k=0; k<dimZ; k++)
    {
      long long ij  = k+dimZ*j;
      long long ijk = k+dimZ*(j+dimY*x0);

      slice[ij][0] = Phi[ijk][0];
      slice[ij][1] = Phi[ijk][1];
    }
  }
}

void extract_y_const( fftw_complex *Phi, fftw_complex *slice, generic_header *header, const int y0 )
{
  const long long dimX = header->nDimX;
  const long long dimY = header->nDimY;
  const long long dimZ = header->nDimZ;

  #pragma omp parallel for collapse(2)
  for ( long long i=0; i<dimX; i++ )
  {
    for ( long long k=0; k<dimZ; k++)
    {
      long long ij  = k+dimZ*i;
      long long ijk = k+dimZ*(y0+dimY*i);

      slice[ij][0] = Phi[ijk][0];
      slice[ij][1] = Phi[ijk][1];
    }
  }
}

void extract_z_const( fftw_complex *Phi, fftw_complex *slice, generic_header *header, const int z0 )
{
  const long long dimX = header->nDimX;
  const long long dimY = header->nDimY;
  const long long dimZ = header->nDimZ;

  #pragma omp parallel for collapse(2)
  for ( long long i=0; i<dimX; i++ )
  {
    for ( long long j=0; j<dimY; j++)
    {
      long long ij  = j+dimY*i;
      long long ijk = z0+dimZ*(j+dimY*i);

      slice[ij][0] = Phi[ijk][0];
      slice[ij][1] = Phi[ijk][1];
    }
  }
}

int main( int argc, const char *argv[])
{
  if ( argc != 5 )
  {
    printf( "slice_3d filename x x idx\n" );
    printf( "slice_3d filename x idx x\n" );
    printf( "slice_3d filename idx x x\n" );
    return 0;
  }

  FILE *fh = nullptr;

  int no_of_threads = 4;
  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );
  omp_set_num_threads( no_of_threads );

  generic_header header_3d = {};
  generic_header header_2d = {};

  fh = fopen( argv[1], "r" );
  if ( fh == nullptr )
  {
    printf( "Could not open file %s.\n", argv[1] );
    exit(0);
  }

  fread( &header_3d, sizeof(generic_header), 1, fh );
  memcpy( &header_2d, &header_3d, sizeof(generic_header) );
  header_2d.nDims = 2;
  header_2d.nDimZ = 1;

  printf( "### %s\n", argv[1] );
  printf( "# nDims    == %lld\n", header_3d.nDims );
  printf( "# nDimX    == %lld\n", header_3d.nDimX );
  printf( "# nDimY    == %lld\n", header_3d.nDimY );
  printf( "# nDimZ    == %lld\n", header_3d.nDimZ );
  printf( "# nDatatyp == %lld\n", header_3d.nDatatyp );
  printf( "# bAtom    == %d\n", header_3d.bAtom );
  printf( "# bComplex == %d\n", header_3d.bComplex );
  printf( "# t        == %g\n", header_3d.t );
  printf( "# dt       == %g\n", header_3d.dt );
  printf( "# xMin     == %g\n", header_3d.xMin );
  printf( "# xMax     == %g\n", header_3d.xMax );
  printf( "# yMin     == %g\n", header_3d.yMin );
  printf( "# yMax     == %g\n", header_3d.yMax );
  printf( "# zMin     == %g\n", header_3d.zMin );
  printf( "# zMax     == %g\n", header_3d.zMax );
  printf( "# dx       == %g\n", header_3d.dx );
  printf( "# dy       == %g\n", header_3d.dy );
  printf( "# dz       == %g\n", header_3d.dz );
  printf( "# dkx      == %g\n", header_3d.dkx );
  printf( "# dky      == %g\n", header_3d.dky );
  printf( "# dkz      == %g\n", header_3d.dkz );

  char filename[255];
  long nBytes=0;

  if ( header_3d.nDims == 3 && header_3d.bComplex == 1 && header_3d.bAtom == 1 && header_3d.nDatatyp == sizeof(fftw_complex) )
  {
    // Lesen von Psi[ijk]
    fseek( fh, sizeof(generic_header), SEEK_SET );

    size_t Nges = header_3d.nDimX*header_3d.nDimY*header_3d.nDimZ;

    fftw_complex *field = (fftw_complex *)fftw_alloc_complex( header_3d.nDatatyp*Nges );
    fftw_complex *slice = nullptr;

    fread( (void *)field, header_3d.nDatatyp, Nges, fh );
    fclose(fh);

    // x const
    if ( strcmp(argv[3],"x") == 0 && strcmp(argv[4],"x") == 0 )
    {
      int x0 = atoi( argv[2] );
      if ( x0 > header_3d.nDimX || x0 < 0 )
      {
        printf( "index out of bounds\n" );
        free( field );
        exit(0);
      }
      nBytes = header_3d.nDatatyp*header_3d.nDimY*header_3d.nDimZ;
      slice = (fftw_complex *)fftw_malloc( nBytes );
      extract_x_const( field, slice, &header_3d, x0 );
      sprintf( filename, "slice_%.3g__%d_x_x.bin", header_3d.t, x0 );
    }

    // y const
    if ( strcmp(argv[2],"x") == 0 && strcmp(argv[4],"x") == 0 )
    {
      int y0 = atoi( argv[3] );
      if ( y0 > header_3d.nDimY || y0 < 0 )
      {
        printf( "index out of bounds\n" );
        free( field );
        exit(0);
      }
      nBytes = header_3d.nDatatyp*header_3d.nDimX*header_3d.nDimZ;
      slice = (fftw_complex *)fftw_malloc( nBytes );
      extract_y_const( field, slice, &header_3d, y0 );
      sprintf( filename, "slice_%.3g__x_%d_x.bin", header_3d.t, y0 );
    }

    // z const
    if ( strcmp(argv[2],"x") == 0 && strcmp(argv[3],"x") == 0 )
    {
      int z0 = atoi( argv[4] );
      if ( z0 > header_3d.nDimZ || z0 < 0 )
      {
        printf( "index out of bounds\n" );
        free( field );
        exit(0);
      }
      nBytes = header_3d.nDatatyp*header_3d.nDimX*header_3d.nDimY;
      slice  = (fftw_complex *)fftw_malloc( nBytes );
      extract_z_const( field, slice, &header_3d, z0 );
      sprintf( filename, "slice_%.3g__x_x_%d.bin", header_3d.t, z0 );
    }
    fclose(fh);

    FILE *fh2 = fopen( filename, "w" );
    fwrite( &header_2d, sizeof(generic_header),1,fh2);
    fwrite( slice, nBytes, 1, fh2 );
    fclose(fh2);

    free( slice );
    free( field );
  }
}
