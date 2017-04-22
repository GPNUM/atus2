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

/** Calculates the PSD
 *
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <omp.h>
#include "my_structs.h"
#include "noise3_2d.h"
#include "ParameterHandler.h"

int main(int argc, char *argv[])
{

  if( argc < 2 )
    {
      printf( "No signal binary file specified.\n" );
      return EXIT_FAILURE;
    }

  generic_header header;
  ifstream fsignal(argv[1], ifstream::binary );
  fsignal.read( (char*)&header, sizeof(generic_header));

  int NX = header.nDimX;
  int NY = header.nDimY;

  Fourier::cft_2d ft(header);
  fftw_complex* signal = ft.Getp2In();

  fsignal.read( (char*)signal, sizeof(fftw_complex)*NX*NY );
  fsignal.close();

  ft.ft(-1);
  fftw_complex* fourier = ft.Getp2Out();

  double dx = header.dx;
  double dy = header.dy;

  int ij = 0;
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      ij = j + NY*i;
      signal[ij][0] = pow(dx*dy*fourier[ij][0],2);
    }
  }

  char* bin_header = reinterpret_cast<char*>(&header);
  char* bin_signal = reinterpret_cast<char*>(signal);

  ofstream file1( "psd.bin", ofstream::binary );
  file1.write( bin_header, sizeof(generic_header) );
  file1.write( bin_signal, NX*NY*sizeof(fftw_complex) );
  file1.close();

}
