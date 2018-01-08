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

#include <iomanip>
#include <sstream>
#include <fstream>
#include <omp.h>
#include "my_structs.h"
#include "ParameterHandler.h"
#include "fftw3.h"

using namespace std;

int main(int argc, char *argv[])
{
  if( argc < 2 )
  {
    cout << "Error: No data binary file specified." << endl;
    cout << "Usage: " << argv[0] << " Seq_1_1.bin [_N postfix =1]" << endl;
    return EXIT_FAILURE;
  }

  string postfix = "1";
  if (argc > 2) {
    postfix = string(argv[2]);
  }

  generic_header header;
  ifstream fdata(argv[1], ifstream::binary );
  if (fdata.fail()) {
    cout << "File not found: " << argv[1] << endl;
    exit(EXIT_FAILURE);
  }

  // Unpack wavefunctions from sequence pack
  while (true) {
    fdata.read( (char*)&header, sizeof(generic_header));
    if (fdata.eof()) {
      break;
    }

    // Create Output file
    stringstream stream;
    stream << fixed << setprecision(3) << header.t;
    string file = stream.str() + "_" + postfix + ".bin";

    ofstream fout(file, ofstream::binary);
    if (fout.fail()) {
      cout << "Could not create file : " << file  << endl;
      exit(EXIT_FAILURE);
    }

    fout.write( (char*)&header, sizeof(generic_header) );

    int64_t N = header.nDimX*header.nDimY*header.nDimZ;
    if (header.bComplex) {
      fftw_complex data[N];
      fdata.read( (char*)data, N*sizeof(fftw_complex) );
      fout.write( (char*)data, N*sizeof(fftw_complex) );
    } else {
      double data[N];
      fdata.read( (char*)data, N*sizeof(double) );
      fout.write( (char*)data, N*sizeof(double) );
    }
  }
}
