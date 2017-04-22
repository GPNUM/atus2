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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
      printf( "No parameter xml file specified.\n" );
      return EXIT_FAILURE;
    }

  int no_of_threads = 4;
  char* envstr = getenv( "MY_NO_OF_THREADS" );
  if( envstr != NULL ) no_of_threads = atoi( envstr );
  printf("Number of threads %i\n", no_of_threads);

  omp_set_num_threads( no_of_threads );

  generic_header header;
  std::ifstream fdata(argv[1], ifstream::binary );
  fdata.read( (char*)&header, sizeof(generic_header));

  int N = 1000;
  int prob_dist[N] = {};
  double nMin = -1.0;
  if (argc > 2) {
    nMin = atof(argv[2]);
  }
  double nMax = -nMin;
  double dn = (nMax-nMin)/N;
  std::cout << "nDimX \t" << header.nDimX << std::endl;
  std::cout << "nDimY \t" << header.nDimY << std::endl;
  printf("N %i\n", N);
  printf("dn %f\n", dn);
  printf("nMin %f\n", nMin);
  printf("nMax %f\n", nMax);

  double *data_real;
  fftw_complex *data_complex;
  if (header.bComplex) {
    data_complex = new fftw_complex[header.nDimY];
    for (int i = 0; i < header.nDimX; i++) {
      fdata.read( (char*)data_complex, sizeof(fftw_complex)*header.nDimY );
      for (int j = 0; j <header.nDimY; j++) {
        prob_dist[static_cast<int>((data_complex[j][0]-nMin)/dn)] += 1;
      }
    }
  } else {
    data_real = new double[header.nDimY];
    for (int i = 0; i < header.nDimX; i++) {
      fdata.read( (char*)data_real, sizeof(double)*header.nDimY );
      for (int j = 0; j <header.nDimY; j++) {
        prob_dist[static_cast<int>((data_real[j]-nMin)/dn)] += 1;
      }
    }
  }

  ofstream file1( "prob_dist.txt");
  for (int i = 0; i < N; i++) {
    file1 << i*dn+nMin << "\t" << prob_dist[i]/static_cast<double>(header.nDimX*header.nDimY) << std::endl;
  }
  fdata.close();

  if (header.bComplex) {
    delete[] data_complex;
  } else {
    delete[] data_real;
  }

}
