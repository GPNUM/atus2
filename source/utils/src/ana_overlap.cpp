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
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"
#include "fftw3.h"

using namespace std;

bool Read( const char *filename, const generic_header &header, double * field )
{
  FILE *fh = fopen( filename, "r" );

  if ( fh == nullptr ) return false;

  fseek( fh, sizeof(generic_header), SEEK_SET );
  size_t Nges = header.nDimX * header.nDimY * header.nDimZ;

  fread( (void *)field, 2*sizeof(double), Nges, fh );
  fclose(fh);
  return true;
}

int main(int argc, char *argv[])
{
  FILE *fh = nullptr;

  double PN;

  int no_of_threads = 4;
  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );
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
    // printf( "# dkx      == %g\n", header.dkx );
    // printf( "# dky      == %g\n", header.dky );
    // printf( "# dkz      == %g\n", header.dkz );
    fclose(fh);
  }

  if ( header.nDims == 1 )
  {
    header.nDimY=1;
    header.nDimZ=1;
    header.dy=1;
    header.dz=1;
  };
  if ( header.nDims == 2 )
  {
    header.nDimZ=1;
    header.dz=1;
  };

  long long Ntot = header.nDimX*header.nDimY*header.nDimZ;
  double * wf1 = new double[2*Ntot];
  double * wf2 = new double[2*Ntot];
  
  Read( argv[1], header, wf1 );
  Read( argv[2], header, wf2 );

  double re=0, im=0, N1=0, N2=0, m1=0, m2=0;
  #pragma omp parallel for reduction(+:re,im,N1,N2) reduction(max:m1,m2)
  for ( long long i=0; i<2*Ntot; i+=2 )
  {
    re += (wf1[i]*wf2[i] + wf1[i+1]*wf2[i+1]);
    im += (wf1[i]*wf2[i+1] - wf1[i+1]*wf2[i]);
    N1 += (wf1[i]*wf1[i] + wf1[i+1]*wf1[i+1]);
    N2 += (wf2[i]*wf2[i] + wf2[i+1]*wf2[i+1]);
    m1 = std::max( m1, wf1[i]*wf1[i] + wf1[i+1]*wf1[i+1] );
    m2 = std::max( m2, wf2[i]*wf2[i] + wf2[i+1]*wf2[i+1] );
  }
  
  re *= header.dx*header.dy*header.dz;
  im *= header.dx*header.dy*header.dz;
  N1 *= header.dx*header.dy*header.dz;
  N2 *= header.dx*header.dy*header.dz;
  double overlap = (re*re + im*im)/sqrt(N1)/sqrt(N2);

  std::cout << "overlap = " << overlap << std::endl;
  std::cout << "N1 = " << N1 << ", N2 = " << N2 << std::endl;
  std::cout << "max density = " << m1/N1 << ", " << m2/N2 << std::endl;

  delete [] wf1;
  delete [] wf2;
  
  return EXIT_SUCCESS;
}
