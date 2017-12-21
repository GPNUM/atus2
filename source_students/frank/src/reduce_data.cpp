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

using namespace std;

long pixel = 512;

int main(int argc, char *argv[])
{
  if( argc < 2 )
  {
    cout << "Reduces 2d array to " << pixel << "x" << pixel << " array." << endl;
    cout << "Usage: reduce_data data.bin [option] [x-offset] [y-offset]" << endl;
    cout << "Error: No signal binary file specified." << endl;
    return EXIT_FAILURE;
  }

  generic_header header;
  ifstream fsignal(argv[1], ifstream::binary );
  if (fsignal.fail()) {
    cout << "File not found: " << argv[1] << endl;
    exit(EXIT_FAILURE);
  }
  fsignal.read( (char*)&header, sizeof(generic_header));
  long NX = header.nDimX;
  long NY = header.nDimY;

  if ( (NX < pixel) || (NY < pixel) ) {
    cout << "Data already small" << endl;
    return 0;
  }

  long rNX = NX/pixel;
  long rNY = NY/pixel;

  header.nDimX = pixel;
  header.nDimY = pixel;
  header.dx *= rNX;
  header.dy *= rNY;

  if (argc == 2) {
    cout << "Normal Mode" << endl;
  } else {
    cout << "Alternative Mode" << endl;
  }
  cout << NX << "\t" << NY << "\t" << endl;
  cout << rNX << "\t" << rNY << "\t" << endl;

  double *signal_real;
  double *rsignal_real;
  fftw_complex *signal_complex;
  fftw_complex *rsignal_complex;

  if (header.bComplex) {
    signal_complex = fftw_alloc_complex( NX * NY );
    rsignal_complex = fftw_alloc_complex( pixel * pixel );
    fsignal.read( (char*)signal_complex, sizeof(fftw_complex)*NX*NY );
    fsignal.close();
  } else {
    signal_real = new double[ NX * NY ];
    rsignal_real = new double[ pixel * pixel ];
    fsignal.read( (char*)signal_real, sizeof(double)*NX*NY );
    fsignal.close();
  }

  long ij;
  long rij;
  if (argc == 2) {
    if (header.bComplex) {
      for (long i = 0; i < pixel; i++) {
        for (long j = 0; j < pixel; j++) {
          rij = j + i*pixel;
          ij = j*rNY + i*NY*rNX;
          rsignal_complex[rij][0] = signal_complex[ij][0];
        }
      }
    } else {
      for (long i = 0; i < pixel; i++) {
        for (long j = 0; j < pixel; j++) {
          rij = j + i*pixel;
          ij = j*rNY + i*NY*rNX;
          rsignal_real[rij] = signal_real[ij];
        }
      }
    }

  } else {
    int64_t NX_offset = 0;
    int64_t NY_offset = 0;
    if (argc > 4) {
      NX_offset = atoi(argv[3]);
      NY_offset = atoi(argv[4]);
      assert(NX_offset < (NX-pixel));
      assert(NY_offset < (NY-pixel));
      cout << "New offsets" << endl;
    }
    if (header.bComplex) {
      for (long i = 0; i < pixel; i++) {
        for (long j = 0; j < pixel; j++) {
          rij = j + i*pixel;
          ij = (j+NY_offset) + (i+NX_offset)*NY;
          rsignal_complex[rij][0] = signal_complex[ij][0];
        }
      }
    } else {
      for (long i = 0; i < pixel; i++) {
        for (long j = 0; j < pixel; j++) {
          rij = j + i*pixel;
          ij = (j+NY_offset) + (i+NX_offset)*NY;
          rsignal_real[rij] = signal_real[ij];
        }
      }
    }
  }

  char* bin_header = reinterpret_cast<char*>(&header);
  string foo = string(argv[1]);
  ofstream file1;
  if (argc == 2) {
    file1.open( "red_"+foo, ofstream::binary );
  } else {
    file1.open( "red_alt_"+foo, ofstream::binary );
  }
  file1.write( bin_header, sizeof(generic_header) );
  char* bin_signal;
  if (header.bComplex) {
    bin_signal = reinterpret_cast<char*>(rsignal_complex);
    file1.write( bin_signal, pixel*pixel*sizeof(fftw_complex) );
  } else {
    bin_signal = reinterpret_cast<char*>(rsignal_real);
    file1.write( bin_signal, pixel*pixel*sizeof(double) );
  }

  file1.close();
  if (header.bComplex) {
    fftw_free(signal_complex);
    fftw_free(rsignal_complex);
  } else {
    delete[] signal_real;
    delete[] rsignal_real;
  }

}
