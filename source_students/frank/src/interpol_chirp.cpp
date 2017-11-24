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
#include "rft_1d.h"
#include "ParameterHandler.h"

using namespace std;

int main(int argc, char *argv[])
{

  if( argc < 2 )
  {
    cout << "No Chirp file specified." << endl;
    return EXIT_FAILURE;
  }

  string sgarbage;
  double garbage;
  ifstream fdata_count(argv[1]);
  if (fdata_count.fail()) {
    cout << "File not found: " << argv[1] << endl;
    exit(EXIT_FAILURE);
  }
  int count =0;
  getline(fdata_count, sgarbage);
  getline(fdata_count, sgarbage);
  while (fdata_count.eof() == 0) {
    count++;
    getline(fdata_count, sgarbage);
  }

  generic_header header, header_interpol;
  header.nDims = 1;
  header.nDimX = count;
  header.nDimY = 1;
  header.nDimZ = 1;
  header.dx = 2.0*M_PI/header.nDimX;
  header.xMin = 0;
  header.xMax = 2.0*M_PI-header.dx;
  header.dkx = 1;

  cout << header.nDimX << endl;
  double expansion = 100.0;
  header_interpol = header;
  header_interpol.nDimX *= expansion;
  header_interpol.dx /= expansion;

  Fourier::rft_1d ft(header);
  Fourier::rft_1d interpol_ft(header_interpol);
  double* data = ft.Getp2InReal();

  ifstream fdata(argv[1]);
  if (fdata.fail()) {
    cout << "File not found: " << argv[1] << endl;
    exit(EXIT_FAILURE);
  }
  getline(fdata, sgarbage);
  for (int i = 0; i < header.nDimX; i++) {
    fdata >> garbage;
    fdata >> garbage;
    fdata >> data[i];
    cout << data[i] << endl;
    fdata >> garbage;
    fdata >> garbage;
  }

  ft.ft(-1);
  fftw_complex* fourier = ft.Getp2Out();
  fftw_complex* interpolft_out = interpol_ft.Getp2Out();

  memset( reinterpret_cast<void*>(interpolft_out), 0, sizeof(fftw_complex)*interpol_ft.Get_Dim_FS() );
  for( int i=0; i<ft.Get_Dim_FS(); i++ )
  {
    interpolft_out[i][0] = fourier[i][0];
    interpolft_out[i][1] = fourier[i][1];
  }

  interpol_ft.ft(1);

  double* interpolft_in = interpol_ft.Getp2InReal();

  ofstream file1("Chirp_interpol.txt");
  if (file1.fail()) {
    cout << "Error opening file: " << "Chirp_interpol.txt" << endl;
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < interpol_ft.Get_Dim_X(); i++) {
    file1 << header_interpol.xMin+header_interpol.dx*i << "\t" << interpolft_in[i] << endl;
  }

}
