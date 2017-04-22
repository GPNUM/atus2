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

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "CPoint.h"
#include "cft_2d.h"

/**
 * \brief Setup testfunction
 */
void fkt1( Fourier::cft_2d &cft )
{
  CPoint<2> x;

  fftw_complex *data = cft.Getp2In();

  if ( data == nullptr )
  {
    std::cerr << "Critical error occurred: data == nullptr" << std::endl;
    throw;
  }

  #pragma omp parallel for private(x)
  for ( int64_t l=0; l<cft.Get_Dim_RS(); l++ )
  {
    x = cft.Get_x(l);
    data[l][0] = exp(-(x*x));
    data[l][1] = 0.0;
  }
}

int main(int argc, char *argv[])
{
  generic_header header = {};
  header.nDatatyp = sizeof(fftw_complex);
  header.bComplex = true;
  header.nDims = 2;
  header.nDimX = 512;
  header.nDimY = 512;
  header.nDimZ = 1;
  header.xMin = -5.0;
  header.xMax = -header.xMin;
  header.yMin = -5.0;
  header.yMax = -header.yMin;
  const double LX = fabs( header.xMax-header.xMin );
  const double LY = fabs( header.yMax-header.yMin );
  header.dkx = 2*M_PI/LX;
  header.dky = 2*M_PI/LY;
  header.dx = LX/header.nDimX;
  header.dy = LY/header.nDimY;
  header.nself_and_data = sizeof(generic_header) + header.nDatatyp * header.nDimX * header.nDimY * header.nDimZ;

  fftw_init_threads();
  fftw_plan_with_nthreads( 4 );
  omp_set_num_threads( 4 );

  printf( "omp_get_num_threads  == %d\n", omp_get_num_threads());
  printf( "omp_get_max_threads  == %d\n", omp_get_max_threads());
  printf( "dimX                 == %lld\n", header.nDimX );
  printf( "dimY                 == %lld\n", header.nDimY );
  printf( "Anzahl Stuetzstellen == %lld\n", header.nDimX*header.nDimY );
  printf( "Benoetigter Speicher == %.4f MB\n", float(2*header.nDimX*header.nDimY*sizeof(double)/1024/1024) );
  printf( "dx                   == %g\n", header.dx );
  printf( "dy                   == %g\n", header.dy );
  printf( "dkx                  == %g\n", header.dkx );
  printf( "dky                  == %g\n", header.dky );

  Fourier::cft_2d cft( header );

  fkt1( cft );

  cft.save( "f0.bin" );
  cft.ft(-1);

  CPoint<2> k;
  fftw_complex *data = cft.Getp2Out();
  #pragma omp parallel for private(k)
  for ( int64_t l=0; l<cft.Get_Dim_FS(); l++ )
  {
    k = cft.Get_k(l);
    double fak = -(k*k);
    data[l][0] *= fak;
    data[l][1] *= fak;
  }

  cft.ft(1);
  cft.save( "f1.bin" );

  fftw_cleanup_threads();
  return EXIT_SUCCESS;
}
