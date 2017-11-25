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

#include "ParameterHandler.h"
#include <iostream>
#include <fstream>

using namespace std;

const int fak = 8;

int main(int argc, char *argv[])
{
  if( argc < 3 )
  {
    cout << "Error: Too few arguments\n"
         << "Parameter xml file and wavefunction .bin needed.\n"
         << "eg. " << argv[0] << " params.xml inf_0.000_1.bin " << endl;
    return EXIT_FAILURE;
  }

  ParameterHandler params(argv[1]); // XML file
  double dt = 1000.0;
  double duration = 0.0;
  for (auto seq_item : params.m_sequence) {
    if( seq_item.name != "set_momentum" ) { // Ignore set momentum step
      if (seq_item.dt < dt) {
        dt = seq_item.dt;
      }
      duration += seq_item.duration.front();
    }
  }
  cout << "Duration: " << duration << endl;
  cout << "dt: " << dt << endl;
  assert(dt > 0.0);
  assert(duration > 0.0);

  generic_header header = {};
  {
    ifstream fheader( argv[2], ifstream::binary );
    if (fheader.fail()) {
      cout << "Error: File not found: " << argv[2] << endl;
      exit(EXIT_FAILURE);
    }
    fheader.read( (char*)&header, sizeof(generic_header));
    fheader.close();
  }
  assert(header.nDims > 0);
  assert(header.nDimX > 0);

  header.nself    = sizeof(generic_header);
  header.nDatatyp = sizeof(double);
  header.bComplex = 0;
  header.nDims    = 2;
  header.nDimZ    = 1;
  header.nDimY    = header.nDimX/fak;
  header.nDimX    = duration/dt/fak;
  header.yMin     = header.xMin;
  header.yMax     = header.xMax;
  header.xMin     = 0;
  header.xMax     = duration;
  header.dy       = header.dx*fak;
  header.dky      = 2.0*M_PI/fabs(header.yMax-header.yMin);
  header.dx       = fabs( header.xMax-header.xMin )/double(header.nDimX);
  header.dkx      = 2.0*M_PI/fabs(header.xMax-header.xMin);
  header.nself_and_data = header.nself + (header.nDimX*header.nDimY*header.nDimZ)*header.nDatatyp;

  double *data = new double[ header.nDimX * header.nDimY ];

  #pragma omp parallel for collapse(2)
  for (int64_t i = 0; i < header.nDimX; ++i) {
    for (int64_t j = 0; j < header.nDimY; ++j) {
      int64_t ij = j+i*header.nDimY;
      double x = header.xMin + i*header.dx;
      double y = header.yMin + j*header.dy;
      data[ij] = sin(x/10.0)*sin(y/10.0);
    }
  }

  ofstream file1( "sine.bin", ofstream::binary );
  file1.write( reinterpret_cast<char*>(&header), sizeof(generic_header) );
  file1.write( reinterpret_cast<char*>(data), sizeof(double)*header.nDimX*header.nDimY );

  delete[] data;

  return 0;
}
