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
#include <cstring>
#include <cmath>
#include <deque>
#include <vector>
#include <fstream>
#include <algorithm>
#include "noise.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

using namespace std;

int main(int argc, char *argv[])
{
  const int fak = 4;

  generic_header header = {};

  header.nself    = sizeof(generic_header);
  header.nDatatyp = sizeof(fftw_complex);
  header.nDims    = 1;
  header.nDimX    = 256;
  header.nDimY    = 1;
  header.nDimZ    = 1;
  header.dt       = 0.1;
  header.xMin     = -10;
  header.xMax     = -header.xMin;
  header.dx       = (header.xMax-header.xMin) / double(header.nDimX);
  header.dkx      = 2.0*M_PI/fabs(header.xMax-header.xMin);
  header.nself_and_data = header.nself + (header.nDimX*header.nDimY*header.nDimZ)*header.nDatatyp;

  generic_header header2 = header;
  header2.nDimX = fak * header.nDimX;
  header2.dx = header2.dx/double(fak);


  Fourier::CNoise<Fourier::rft_1d,1> noise( header );
  noise.color_noise(2);
  noise.save("noise.bin");
  noise.ft(-1);

  Fourier::rft_1d interpol( header2 );

  std::cout << "dimX == " << noise.Get_Dim_X() << std::endl;
  std::cout << "dimX == " << interpol.Get_Dim_X() << std::endl;
  std::cout << "dimX == " << noise.Get_red_Dim() << std::endl;
  std::cout << "dimX == " << interpol.Get_red_Dim() << std::endl;

  // extract chunk
  fftw_complex * noise_out = noise.Getp2Out();
  fftw_complex * interpol_out = interpol.Getp2Out();

  memcpy( reinterpret_cast<void*>(interpol_out), reinterpret_cast<void*>(noise_out), sizeof(fftw_complex) * noise.Get_red_Dim() );

  interpol.save( "inoise_ft.bin" );
  interpol.ft(1);
  interpol.save( "inoise.bin" );

  fftw_cleanup_threads();
return EXIT_SUCCESS;
}
