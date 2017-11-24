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
#include "fftw3.h"

using namespace std;

int main(int argc, char *argv[])
{

  if (argc < 2) {
    cout << "Missing arguments" << endl;
    return EXIT_FAILURE;
  }
  cout << "Initialization...";
  cout.flush();
  generic_header header = {};

  cout << "File: " << argv[1] << endl;
  ifstream fdata(argv[1], ifstream::binary );
  if (fdata.fail()) {
    cout << "File not found: " << argv[1] << endl;
    abort();
  }
  fdata.read( (char*)&header, sizeof(generic_header));

  if (header.nDims != 2) {
    cout << "nDims == 2 required." << endl;
    abort();
  }

  if (header.bComplex) {
    cout << "Only real data supported" << endl;
    abort();
  }

  cout << "Done." << endl;
  cout << "Reading Data...";
  cout.flush();
  double* data = new double[header.nDimX*header.nDimY];
  assert(data != nullptr);
  fdata.read( (char*)data, header.nDimX*header.nDimY*sizeof(double));
  cout << "Done." << endl;

  cout << header.nDimX << " " << header.nDimY << endl;
  for (int64_t i = 1; i < header.nDimX; ++i) {
    for (int64_t j = 0; j < header.nDimY; ++j) {
      data[j+i*header.nDimY] += data[j+(i-1)*header.nDimY];
    }
  }
  cout << "Write to file...";
  cout.flush();
  char* bin_header = reinterpret_cast<char*>(&header);
  ofstream file1("integral-fun.bin", ofstream::binary);
  file1.write(bin_header, sizeof(generic_header));
  file1.write(reinterpret_cast<char*>(data), header.nDimX*header.nDimY*sizeof(double));
  cout << "Done." << endl;


  return 0;
}
