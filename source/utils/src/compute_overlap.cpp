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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "fftw3.h"
#include "my_structs.h"

using namespace std;

bool Read( const char *filename, const generic_header &header, fftw_complex *field )
{
  FILE *fh = fopen( filename, "r" );

  if ( fh == nullptr ) return false;

  fseek( fh, sizeof(generic_header), SEEK_SET );
  size_t Nges = header.nDimX * header.nDimY * header.nDimZ;

  fread( (void *)field, sizeof(fftw_complex), Nges, fh );
  fclose(fh);
  return true;
}

void Scale( const generic_header &header, const double fak, fftw_complex *u )
{
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  #pragma omp parallel for
  for ( long long l=0; l<Nges; l++ )
  {
    u[l][0] *= fak;
    u[l][1] *= fak;
  }
}

double Particle_Number( const generic_header &header, fftw_complex *u )
{
  double retval=0;
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  #pragma omp parallel for reduction(+:retval)
  for ( long long l=0; l<Nges; l++ )
  {
    retval += u[l][0]*u[l][0] + u[l][1]*u[l][1];
  }
  return retval;
}

double Overlap( const generic_header &header, fftw_complex *u, fftw_complex *v )
{
  double retval=0;
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  #pragma omp parallel for reduction(+:retval)
  for ( long long l=0; l<Nges; l++ )
  {
    retval += u[l][0]*v[l][0] + u[l][1]*v[l][1];
  }
  return retval;
}

int main(int argc, char *argv[])
{
  if ( argc != 3 )
  {

    return EXIT_SUCCESS;
  }

  FILE *fh = nullptr;
  FILE *fh_2 = nullptr;

  int no_of_threads = 4;
  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );
  omp_set_num_threads( no_of_threads );

  generic_header header = {};
  generic_header header_2 = {};

  fh = fopen( argv[1], "r" );
  if ( fh == nullptr )
  {
    printf( "Could not open file %s.\n", argv[1] );
    exit(0);
  }

  fh_2 = fopen( argv[2], "r" );
  if ( fh_2 == nullptr )
  {
    printf( "Could not open file %s.\n", argv[2] );
    exit(0);
  }

  fread( &header, sizeof(generic_header), 1, fh );
  fread( &header_2, sizeof(generic_header), 1, fh_2 );

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

  if ( header.nDims == 1 )
  {
    header.nDimY=1;
    header.nDimZ=1;
  };
  if ( header.nDims == 2 )
  {
    header.nDimZ=1;
  };

  long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  fftw_complex *wf1 = fftw_alloc_complex( Nges );
  fftw_complex *wf2 = fftw_alloc_complex( Nges );

  Read( argv[1], header, wf1 );
  Read( argv[2], header, wf2 );

  double N1 = Particle_Number( header, wf1 );
  double N2 = Particle_Number( header, wf2 );

  double overlap;

  //1D Stuff
  if ( header.nDims == 1 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    N1 *= header.dx;
    N2 *= header.dx;

    printf( "N1 = %g, N2 = %g\n", N1, N2 );

    Scale( header, 1.0/sqrt(N1), wf1 );
    Scale( header, 1.0/sqrt(N2), wf2 );

    overlap = Overlap( header, wf1, wf2 );
    overlap *= header.dx;
  }

  // 2D Stuff
  if ( header.nDims == 2 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    N1 *= header.dx * header.dy;
    N2 *= header.dx * header.dy;

    printf( "N1 = %g, N2 = %g\n", N1, N2 );

    Scale( header, 1/sqrt(N1), wf1 );
    Scale( header, 1/sqrt(N2), wf2 );

    overlap = Overlap( header, wf1, wf2 );
    overlap *= header.dx * header.dy;
  }

  // 3D Stuff
  if ( header.nDims == 3 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    N1 *= header.dx * header.dy * header.dz;
    N2 *= header.dx * header.dy * header.dz;

    printf( "N1 = %g, N2 = %g\n", N1, N2 );

    Scale( header, 1/sqrt(N1), wf1 );
    Scale( header, 1/sqrt(N2), wf2 );

    overlap = Overlap( header, wf1, wf2 );
    overlap *= header.dx * header.dy * header.dz;
  }

  printf( "overlap == %g\n", overlap );

  fftw_free(wf1);
  fftw_free(wf2);
  return EXIT_SUCCESS;
}
