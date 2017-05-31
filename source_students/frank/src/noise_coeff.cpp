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
#include "rft_1d.h"

int main(int argc, char *argv[])
{
  if( argc < 2 )
    {
      printf( "No signal binary file specified.\n" );
      return EXIT_FAILURE;
    }

  generic_header header;
  ifstream fsignal(argv[1], ifstream::binary );
  if (fsignal.fail()) {
    std::cout << "File not found: " << argv[1] << std::endl;
    abort();
  }
  fsignal.read( (char*)&header, sizeof(generic_header));
  long NT = header.nDimX;
  long NX = header.nDimY;

  std::cout << NT << "\t" << NX << "\t" << std::endl;

  generic_header ft_header = header;
  ft_header.nDims = 1;
  ft_header.nDimX = ft_header.nDimY;
  ft_header.dx = ft_header.dy;
  ft_header.xMin = ft_header.yMin;
  ft_header.xMax = ft_header.yMax;
  ft_header.dkx = ft_header.dky;
  ft_header.nDimY = 1;
  ft_header.dy = 1;
  ft_header.yMin = 0;
  ft_header.yMax = 0;
  ft_header.dky = 1;

  Fourier::rft_1d noiseft(ft_header);

  double *signal_real;
  double *noise = noiseft.Getp2InReal();
  vector<double> m_noise;
  vector<double> m_dx_noise;
  vector<double> m_dx2_noise;
  vector<double> out;
  m_noise.resize(NX);
  m_dx_noise.resize(NX);
  m_dx2_noise.resize(NX);
  out.resize(NX);

  // Local functions inserting the prefactor term for 2nd, 1st, 0th derivative
  auto d2 = [&](int i){ return (1.0 - m_noise[i]); };
  auto d1 = [&](int i){ return (-m_dx_noise[i] +
                                0.5*(1.0 - m_noise[i])
                                *m_dx_noise[i]); };
  auto d0 = [&](int i){ return ((1.0 - m_noise[i])
                                *(-0.25*m_dx2_noise[i]/(1.0 + m_noise[i])
                                  + 5.0/16.0 * pow(m_dx_noise[i]
                                                   / (1.0 + m_noise[i]),2))
                                + 0.25*pow(m_dx_noise[i], 2) / (1.0 + m_noise[i])
                                - 0.125 * (1.0 - m_noise[i])*pow(m_dx_noise[i], 2)
                                / (1 + m_noise[i])); };

  char* bin_header = reinterpret_cast<char*>(&header);

  ofstream of_d1( "d1.bin", ofstream::binary );
  ofstream of_d0( "d0.bin", ofstream::binary );
  of_d1.write( bin_header, sizeof(generic_header) );
  of_d0.write( bin_header, sizeof(generic_header) );

  char* bin_signal;
  bin_signal = reinterpret_cast<char*>(out.data());

  for (int64_t i = 0; i < NT; i++) {
    fsignal.read( (char*)noise, sizeof(double)*NX );
    for (int64_t j = 0; j < NX; j++) {
      m_noise[j] = noise[j];
    }
    noiseft.Diff_x();
    for (int64_t j = 0; j < NX; j++) {
      m_dx_noise[j] = noise[j];
    }
    noiseft.Diff_x();
    for (int64_t j = 0; j < NX; j++) {
      m_dx2_noise[j] = noise[j];
    }
    for (int64_t j = 0; j < NX; j++) {
      out[j] = d1(j);
    }
    of_d1.write( bin_signal, NX*sizeof(double) );
    for (int64_t j = 0; j < NX; j++) {
      out[j] = d0(j);
    }
    of_d0.write( bin_signal, NX*sizeof(double) );

  }

}
