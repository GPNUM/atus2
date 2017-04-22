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
#include <cmath>
#include "CPoint.h"
#include "rft_1d.h"

/**
 * \brief Setup testfunction
 */
void fkt1( Fourier::rft_1d &rft )
{
  CPoint<1> x;
  double *data = rft.Getp2InReal();

  if ( data == nullptr )
  {
    std::cerr << "Critical error occurred: data == nullptr" << std::endl;
    throw;
  }

  #pragma omp parallel for private(x)
  for ( int64_t l=0; l<rft.Get_Dim_RS(); l++ )
  {
    x = rft.Get_x(l);
    data[l] = exp(-(x*x));
  }
}

int main(int argc, char *argv[])
{
  generic_header header = {};
  header.nDims = 1;
  header.nDimX = 1024;
  header.nDimY = 1;
  header.nDimZ = 1;
  header.xMin = -10.0;
  header.xMax = -header.xMin;
  const double LX = fabs( header.xMax-header.xMin );
  header.dkx = 2*M_PI/LX;
  header.dx = LX/header.nDimX;
  header.nself_and_data = sizeof(generic_header) + header.nDatatyp * header.nDimX * header.nDimY * header.nDimZ;

  printf( "dx  = %g\n", header.dx );
  printf( "dkx = %g\n", header.dkx );

  Fourier::rft_1d rft( header );

  fkt1( rft );
  rft.save("f0.bin");
  rft.Diff_x();
  //rft.Diff_xx();
  rft.save("f1.bin");

  return EXIT_SUCCESS;
}

