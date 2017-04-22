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

//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
#include "noise.h"

int main(int argc, char *argv[])
{
  const int N = 18*1024;

  generic_header header = {};

  header.nself    = sizeof(generic_header);
  header.nDatatyp = sizeof(fftw_complex);
  header.nDims    = 1;
  header.nDimX    = N;
  header.nDimY    = 1;
  header.nDimZ    = 1;
  header.bAtom    = 1;
  header.bComplex = 1;
  header.xMin     = -160.0;
  header.xMax     = -header.xMin;
  header.dx       = fabs( header.xMax-header.xMin )/double(header.nDimX);
  header.dkx      = 2.0*M_PI/fabs(header.xMax-header.xMin);
  header.nself_and_data = header.nself + header.nDimX*header.nDatatyp;

  Fourier::CNoise<Fourier::cft_1d,1> noise( header );
  noise.white_noise();
  noise.save( "white_noise_1d.bin" );
  noise.color_noise(2);
  noise.save( "color_noise_1d.bin" );
  noise.custom_noise();
  noise.save( "custom_noise_1d.bin" );
  /*
    header.nDims    = 2;
    header.nDimY    = N;
    header.t        = 0.0;
    header.yMin     = -5.0;
    header.yMax     = -header.xMin;
    header.dy       = fabs( header.yMax-header.yMin )/double(header.nDimY);
    header.dky      = 2.0*M_PI/fabs(header.yMax-header.yMin);
    header.nself_and_data = header.nself + (header.nDimX + header.nDimY)*header.nDatatyp;

    Fourier::CNoise<Fourier::cft_2d,2> noise2( header );
    noise2.color_noise(2);
    noise2.save( "color_noise_2d.bin" );
    noise2.white_noise();
    noise2.save( "white_noise_2d.bin" );
    noise2.custom_noise();
    noise2.save( "custom_noise_2d.bin" );
  */
}
