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
    cout <<  "No Number specified." << endl;
    return EXIT_FAILURE;
  }

  int no_of_threads = 4;
  char* envstr = getenv( "MY_NO_OF_THREADS" );

  if( envstr != NULL ) no_of_threads = atoi( envstr );
  cout << "Number of threads " << no_of_threads << endl;

  omp_set_num_threads( no_of_threads );
  ofstream fout("stat.txt");
  char foo[1000];
  int N = stoi(argv[1]);
  int dn = 100;
  if (argc > 2) {
    dn = atoi(argv[2]);
  }

  fout << "# x mean variance sigma N" << endl;
  for (int j = 1; j <= N; j++) {
    sprintf(foo, "%i.000_1.bin", j*dn);
    cout << foo << endl;
    generic_header header;
    ifstream fdata(foo, ifstream::binary );
    if (fdata.fail()) {
      cout << "File not found: " << foo << endl;
      exit(EXIT_FAILURE);
    }
    fdata.read( (char*)&header, sizeof(generic_header));

    int NX = header.nDimX;
    fftw_complex data[NX];
    fdata.read( (char*)data, NX*sizeof(fftw_complex));

    double sum = pow(data[0][0],2) + pow(data[0][1],2)*header.dx/3.0;
    double mean = header.xMin*sum;
    double variance = header.xMin*header.xMin*sum;
    for (int i = 1; i < NX-1; i++) {
      double value = pow(data[i][0],2) + pow(data[i][1],2);
      double x = header.xMin + i*header.dx;
      sum += ((i%2)? 4 : 2)*value*header.dx/3.0;
      mean += ((i%2)? 4 : 2)*x*value*header.dx/3.0;
      variance += ((i%2)? 4 : 2)*x*x*value*header.dx/3.0;
    }
    double val = pow(data[N-1][0],2) + pow(data[N-1][1],2)*header.dx/3.0;
    sum += val;
    mean += header.xMax*val;
    variance += header.xMax*val;
    mean /= sum;
    variance /= sum;

    fout << j*dn << "\t" << mean << "\t" << variance << "\t" << sqrt(abs(variance)) << "\t" << sum << endl;
  }
  cout << "stat.txt written." << endl;
}
