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

/** Generates 2D correlated random numbers and its 1st and 2nd derivatives (DEPRECATED version, use noise_gen instead)
 *
 * noise_generator <param.xml> [<duration> [<dt>]]
 */
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

int main(int argc, char *argv[])
{
  printf("WARNING: Program deprecated. Please use noise_gen.\n");
  if( argc < 2 )
    {
      printf( "No parameter xml file specified.\n" );
      return EXIT_FAILURE;
    }

  int no_of_threads = 4;
  char* envstr = getenv( "MY_NO_OF_THREADS" );
  if( envstr != nullptr ) no_of_threads = atoi( envstr );
  printf("Number of threads %i\n", no_of_threads);

  fftw_init_threads();
  fftw_plan_with_nthreads( no_of_threads );
  omp_set_num_threads( no_of_threads );

  ParameterHandler params(argv[1]); // XML file
  generic_header header;
// ifstream fheader(params.Get_simulation("FILENAME"), ifstream::binary );
// fheader.read( (char*)&header, sizeof(generic_header));
// fheader.close();

  double dt = 0;
  double duration = 0;
  for (auto seq_item : params.m_sequence) {
    if( seq_item.name != "set_momentum" ) { // Ignore set momentum step
      dt = seq_item.dt;
      duration += seq_item.duration.front();
    }
  }

  if (argc > 2) {
    duration = stod(argv[2]);
  }
  if (argc > 3) {
    dt = stod(argv[3]);
  }

  printf("dt, duration %f %f\n", dt, duration);

  long NT = duration/dt;
  long NY = params.Get_NX();

  // Time component is saved in x-dimension component
  header.nself    = sizeof(generic_header);
  header.nDatatyp = sizeof(fftw_complex);
  header.nself_and_data = sizeof(generic_header)+header.nDatatyp*NT*NY;
  header.nDims    = 2;
  header.dx       = dt;
  header.dkx      = 2.0*M_PI/(duration);
  header.nDimX    = NT;
  header.nDimY    = NY;
  header.nDimZ    = 1;
  header.yMin     = params.Get_xMin();
  header.yMax     = params.Get_xMax();
  header.t = NT;
  header.dt = dt;
  header.bAtom = 1;
  header.bComplex = 1;
  header.zMin = 1;
  header.zMax = 1;
  header.dz = 1;
  header.dkz = 1;
  header.dy       = fabs( header.yMax-header.yMin )/double(header.nDimY);
  header.dky      = 2.0*M_PI/fabs(header.yMax-header.yMin);
  header.xMin     = 0;
  header.xMax     = duration;
  header.nself_and_data = header.nself + (header.nDimX + header.nDimY)*header.nDatatyp;

  std::cout << "NY " << NY << std::endl;
  std::cout << "NT " << NT << std::endl;

  Fourier::CNoise2_2D noise_gen( header );
  if (noise_gen.Getp2In() == nullptr) {
    printf("Not Allocated...\n");
    exit(-1);
  }
  double exp_x;
  try {
    exp_x = params.Get_Constant("Noise_xCorr");
  } catch (...) {
    std::cout << "Using default for exp_x" << std::endl;
    exp_x = 2.0;
  }
  double exp_t;
  try {
    exp_t = params.Get_Constant("Noise_tCorr");
  } catch (...) {
    std::cout << "Using default for exp_t" << std::endl;
    exp_t = 2.0;
  }

  double sigma;
  try {
    sigma = params.Get_Constant("Noise_Sigma");
  } catch (...) {
    std::cout << "Using default for sigma" << std::endl;
    sigma = 1.0;
  }

  std::cout << exp_t << "\t" << exp_x << "\t"  << params.Get_Constant("Noise_Amplitude")<< std::endl;

  std::cout.flush();
  noise_gen.Do_Noise_Metric(1.0, sigma, exp_t, exp_x);
  fftw_complex* noise_data = noise_gen.Getp2In();

  char* bin_header = reinterpret_cast<char*>(&header);
  char* noise = reinterpret_cast<char*>(noise_data);

  ofstream file1( "Noise.bin", ofstream::binary );
  file1.write( bin_header, sizeof(generic_header) );
  file1.write( noise, NT*NY*sizeof(fftw_complex) );
  file1.close();
  std::cout << "Noise done." << std::endl;

  // noise_gen.Diff_y();

  // ofstream file2( "dxNoise.bin", ofstream::binary );
  // file2.write( bin_header, sizeof(generic_header) );
  // file2.write( noise, NT*NY*sizeof(fftw_complex) );
  // file2.close();
  // std::cout << "dxNoise done." << std::endl;

  // noise_gen.Diff_y();
  // ofstream file3( "dx2Noise.bin", ofstream::binary );
  // file3.write( bin_header, sizeof(generic_header) );
  // file3.write( noise, NT*NY*sizeof(fftw_complex) );
  // file3.close();
  // std::cout << "dx2Noise done." << std::endl;

}
