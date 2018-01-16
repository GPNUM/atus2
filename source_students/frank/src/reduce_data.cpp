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


int main(int argc, char *argv[])
{
  if( argc < 2 )
  {
    cout << "Reduces 2d array to rNX x rNY array." << endl;
    cout << "Usage: "<< argv[0] << " data.bin [di=1] [dj=1] [rioff=0] [rjoff=0] [rNX=512] [rNY=512]" << endl;
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

  int64_t NX = header.nDimX;
  int64_t NY = header.nDimY;
  int64_t rNX = 512;
  int64_t rNY = 512;
  int64_t rioff = 0;
  int64_t rjoff = 0;
  int64_t di = 1;
  int64_t dj = 1;

  if (argc > 2) {
    di = atoi(argv[2]);
  }
  if (argc > 3) {
    dj = atoi(argv[3]);
  }
  if (argc > 4) {
    rioff = atoi(argv[4]);
  }
  if (argc > 5) {
    rjoff = atoi(argv[5]);
  }
  if (argc > 6) {
    rNX = atoi(argv[6]);
  }
  if (argc > 7) {
    rNY = atoi(argv[7]);
  }

  cout << "NX=" << NX << " NY=" << NY << " rNX=" << rNX << " rNY=" << rNY << " rioff=" << rioff << " rjoff=" << rjoff << " di=" << di << " dj=" << dj<<endl;
  assert(rNX <= NX);
  assert(rNY <= NY);
  assert(((rNX+rioff)*di) <= NX);
  assert(((rNY+rjoff)*dj) <= NY);

  header.nDimX = rNX;
  header.xMin += rioff*header.dx;
  header.dx *= di;
  header.xMax = header.xMin+header.dx*rNX;

  header.nDimY = rNY;
  header.yMin += rjoff*header.dy;
  header.dy *= dj;
  header.yMax = header.yMin+header.dy*rNY;

  double *signal_real;
  double *rsignal_real;
  fftw_complex *signal_complex;
  fftw_complex *rsignal_complex;

  if (header.bComplex) {
    signal_complex = fftw_alloc_complex( NX * NY );
    rsignal_complex = fftw_alloc_complex( rNX * rNY );
    fsignal.read( (char*)signal_complex, sizeof(fftw_complex)*NX*NY );
    fsignal.close();
  } else {
    signal_real = new double[ NX * NY ];
    rsignal_real = new double[ rNX * rNY ];
    fsignal.read( (char*)signal_real, sizeof(double)*NX*NY );
    fsignal.close();
  }

  if (header.bComplex) {
    for (int64_t ri = 0; ri < rNX; ri++) {
      for (int64_t rj = 0; rj < rNY; rj++) {
        int64_t rij = rj + ri*rNY;
        int64_t ij = (rj*dj+rjoff) + (ri*di+rioff)*NY;
        rsignal_complex[rij][0] = signal_complex[ij][0];
        rsignal_complex[rij][1] = signal_complex[ij][1];
      }
    }
  } else {
    for (int64_t ri = 0; ri < rNX; ri++) {
      for (int64_t rj = 0; rj < rNY; rj++) {
        int64_t rij = rj + ri*rNY;
        int64_t ij = (rj*dj+rjoff) + (ri*di+rioff)*NY;
        rsignal_real[rij] = signal_real[ij];
      }
    }
  }

  char* bin_header = reinterpret_cast<char*>(&header);
  string file = string(argv[1]);
  string file_out = "red_"+to_string(di)+"_"+to_string(dj)+"_"+to_string(rioff)+"_"+to_string(rjoff)+"_"+to_string(rNX)+"_"+to_string(rNY)+"_"+file;
  ofstream fout;
  fout.open( file_out, ofstream::binary );

  fout.write( bin_header, sizeof(generic_header) );
  char* bin_signal;
  if (header.bComplex) {
    bin_signal = reinterpret_cast<char*>(rsignal_complex);
    fout.write( bin_signal, rNX*rNY*sizeof(fftw_complex) );
  } else {
    bin_signal = reinterpret_cast<char*>(rsignal_real);
    fout.write( bin_signal, rNX*rNY*sizeof(double) );
  }

  fout.close();
  if (header.bComplex) {
    fftw_free(signal_complex);
    fftw_free(rsignal_complex);
  } else {
    delete[] signal_real;
    delete[] rsignal_real;
  }

}
