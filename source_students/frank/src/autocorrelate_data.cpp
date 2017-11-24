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
#include <fstream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <deque>
#include <vector>
#include "noise.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "my_structs.h"
#include "ParameterHandler.h"

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

  Fourier::rft_2d ft(header);

  cout << "Done." << endl;
  cout << "Reading Data...";
  cout.flush();
  double* data = ft.Getp2InReal();
  assert(data != nullptr);
  fdata.read( (char*)data, header.nDimX*header.nDimY*sizeof(double));

  cout << "Done." << endl;

  cout << "Fourier Transformation...";
  cout.flush();
  ft.ft(-1);
  ft.save("ft-data.bin", false);
  cout << "Done." << endl;

  fftw_complex* fourier = ft.Getp2Out();
  assert(fourier != nullptr);

  cout << "Square data...";
  cout.flush();
  #pragma omp parallel for
  for (int64_t i = 0; i < ft.Get_Dim_FS(); ++i) {
    fourier[i][0] = pow(fourier[i][0],2) + pow(fourier[i][1],2);
    fourier[i][1] = 0.0;
  }
  cout << "Done." << endl;

  cout << "Write to file...";
  generic_header psd_header = header;
  psd_header.nDimX = ft.Get_Dim_X();
  psd_header.nDimY = ft.Get_red_Dim();
  psd_header.fs = 1;
  cout.flush();
  psd_header.bComplex = 1;
  char* bin_psd_header = reinterpret_cast<char*>(&psd_header);
  ofstream file2("power-spectral-density.bin", ofstream::binary);
  file2.write(bin_psd_header, sizeof(generic_header));
  file2.write(reinterpret_cast<char*>(fourier), sizeof(fftw_complex)*ft.Get_Dim_FS());
  cout << "Done." << endl;

  cout << "Backtransformation...";
  cout.flush();
  ft.ft(1);
  cout << "Done." << endl;

  cout << "Write to file...";
  cout.flush();

  char* bin_header = reinterpret_cast<char*>(&header);
  ofstream file1("autocorrelation.bin", ofstream::binary);
  file1.write(bin_header, sizeof(generic_header));
  file1.write(reinterpret_cast<char*>(data), header.nDimX*header.nDimY*sizeof(double));
  cout << "Done." << endl;

  return 0;
}
