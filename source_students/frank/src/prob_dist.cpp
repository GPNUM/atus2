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

using namespace std;

int main(int argc, char *argv[])
{
  if( argc < 2 )
    {
      cout << "Error: No data file specified." << endl;
      cout << "Usage: " << argv[0] << " data.bin [N]" << endl;
      return EXIT_FAILURE;
    }

  int no_of_threads = 4;
  char* envstr = getenv( "MY_NO_OF_THREADS" );
  if( envstr != NULL ) no_of_threads = atoi( envstr );
  cout << "Number of threads: " << no_of_threads << endl;

  omp_set_num_threads( no_of_threads );

  generic_header header;
  ifstream fdata(argv[1], ifstream::binary );
  if (fdata.fail()) {
    cout << "File not found: " << argv[1] << endl;
    abort();
  }
  fdata.read( (char*)&header, sizeof(generic_header));

  double *data = new double[header.nDimX*header.nDimY];
  if (data == nullptr) {
    exit(EXIT_FAILURE);
  }
  fdata.read( (char*)data, sizeof(double)*header.nDimX*header.nDimY );

  int N = 1000;
  if (argc > 2) {
    N = atoi(argv[2]);
  }
  cout << "N: " << N << endl;

  double max_val = 0.0;
  #pragma omp parallel for reduction(max:max_val)
  for (int64_t i = 0; i < header.nDimX*header.nDimY; ++i) {
    if (abs(data[i]) > max_val) {
      max_val = abs(data[i]);
    }
  }

  cout << "Maximum: " << max_val << endl;
  double min_val = -max_val;
  double dn = 2*max_val/N;

  int64_t prob_dist[N] = {};

  #pragma omp parallel
  {
    int64_t priv_prob_dist[N] = {};
    #pragma omp for
    for (int64_t i = 0; i < header.nDimX*header.nDimY; i++) {
      priv_prob_dist[static_cast<int>((data[i]-min_val)/dn)] += 1;
    }
    #pragma omp critical
    {
      for (int64_t i = 0; i < N; ++i) {
        prob_dist[i] += priv_prob_dist[i];
      }
    }
  }

  ofstream file1( "prob_dist.txt");
  for (int i = 0; i < N; i++) {
    file1 << min_val+i*dn << "\t" << prob_dist[i]/static_cast<double>(header.nDimX*header.nDimY) << endl;
  }
  fdata.close();

  delete[] data;

}
