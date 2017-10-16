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
#include <istream>
#include <ostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <omp.h>
#include "my_structs.h"
#include "ParameterHandler.h"
#include "fftw3.h"

using namespace std;

int main(int argc, char *argv[])
{
  if ( argc < 3) {
    std::cout << "Error: Arguments missing.\n Try: " << argv[0] << " data1.bin data2.bin" << std::endl;
    return EXIT_FAILURE;
  }

  generic_header header1, header2;
  ifstream fdata1(argv[1], ifstream::binary);
  if (fdata1.fail()) {
    std::cout << "File not found: " << argv[1] << std::endl;
    abort();
  }

  ifstream fdata2(argv[2], ifstream::binary);
  if (fdata2.fail()) {
    std::cout << "File not found: " << argv[2] << std::endl;
    abort();
  }

  fdata1.read((char*)&header1, sizeof(generic_header));
  fdata2.read((char*)&header2, sizeof(generic_header));

  if (header1.nDims != 1) {
    std::cout << "nDims > 1 not yet implemented...\n";
    abort();
  }
  if (header1.nDims != header2.nDims) {
    std::cout << "nDims not equal size: " << header1.nDims << " != "
              << header2.nDims << std::endl;
    abort();
  }
  if (header1.nDimX != header2.nDimX) {
    std::cout << "nDimX not equal size: " << header1.nDimX << " != "
              << header2.nDimX << std::endl;
    abort();
  }
  if (header1.nDimY != header2.nDimY) {
    std::cout << "nDimY not equal size: " << header1.nDimY << " != "
              << header2.nDimY << std::endl;
    abort();
  }
  if (header1.nDimZ != header2.nDimZ) {
    std::cout << "nDimZ not equal size: " << header1.nDimZ << " != "
              << header2.nDimZ << std::endl;
    abort();
  }

  int64_t nDimXYZ = header1.nDimX*header1.nDimY*header1.nDimZ;
  fftw_complex *data1, *data2;
  data1 = fftw_alloc_complex(nDimXYZ);
  if (data1 == nullptr) {
    std::cout << "Allocation failed. Size wanted: " << sizeof(fftw_complex)*nDimXYZ
              << std::endl;
    abort();
  }
  data2 = fftw_alloc_complex(nDimXYZ);
  if (data2 == nullptr) {
    std::cout << "Allocation failed. Size wanted: " << sizeof(fftw_complex)*nDimXYZ
              << std::endl;
    abort();
  }

  fdata1.read( (char*)data1, sizeof(fftw_complex)*nDimXYZ );
  fdata2.read( (char*)data2, sizeof(fftw_complex)*nDimXYZ );

  double min_err = 1.0e-6;
  if (argc == 4) {
    min_err = atof(argv[3]);
  }

  double max_diff = 0.0;
  int64_t max_i = 0;
  for (int64_t i = 0; i < nDimXYZ; ++i) {
    double diff_real  = fabs(data1[i][0] - data2[i][0]);
    double diff_imag  = fabs(data1[i][1] - data2[i][1]);
    double max_val = max(diff_real, diff_imag);
    if (max_val > max_diff) {
      max_diff = max_val;
      max_i = i;
    }
  }
  if (max_diff > min_err) {
    std::cout << "Not equal! (diff > " << min_err << ")" << std::endl;
    std::cout << "max_i = " << max_i << " max_diff: " << max_diff << std::endl;
    return -1;
  }

  fdata1.close();
  fdata2.close();
  fftw_free(data1);
  fftw_free(data2);

  std::cout << "Equal! (diff < " << min_err << ")" << std::endl;
  std::cout << "max_i = " << max_i << " max_diff: " << max_diff << std::endl;
  return 0;
}
