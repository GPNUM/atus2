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
#include <cmath>
#include <deque>
#include <vector>
#include <fstream>
#include "noise.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "my_structs.h"
#include "ParameterHandler.h"

using namespace std;

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
      printf( "Too few arguments\n" );
      printf( "Parameter xml file and wavefunction .bin needed.\n" );
      printf( "eg. noise_gen params.xml inf_0.000_1.bin\n" );
      return EXIT_FAILURE;
    }

  const int no_chunks = 16;
  const int fak = 4;

  int no_of_threads = 4;
  char* envstr = getenv( "MY_NO_OF_THREADS" );
  if( envstr != NULL ) no_of_threads = atoi( envstr );
  printf("Number of threads %i\n", no_of_threads);

  fftw_init_threads();
  fftw_plan_with_nthreads( no_of_threads );
  omp_set_num_threads( no_of_threads );

  ParameterHandler params(argv[1]); // XML file
  double dt = 1000;
  double duration = 0;
  for (auto seq_item : params.m_sequence) {
    if( seq_item.name != "set_momentum" ) { // Ignore set momentum step
      if (seq_item.dt < dt) {
        dt = seq_item.dt;
      }
      duration += seq_item.duration.front();
    }
  }
  printf("duration: %f\n", duration);
  printf("dt: %f\n", dt);
  assert(dt > 0.0);
  assert(duration > 0.0);
  int exp =ceil(log2(duration/dt));
  printf("exp %i\n", exp);
  duration = pow(2.0, ceil(log2(duration/dt)))*dt;
  printf("new duration: %f\n", duration);
  printf("lower duration: %f\n", pow(2.0, exp-1)*dt);

  generic_header header = {};
  {
    ifstream fheader(argv[2], ifstream::binary );
    fheader.read( (char*)&header, sizeof(generic_header));
    fheader.close();
  }

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

  // header for the chunk
  const double chunk_len =  fabs( header.xMax-header.xMin ) / double(no_chunks);
  const double chunk_dkx =  2.0*M_PI/fabs(chunk_len);

  printf( "dimX == %lld\n", header.nDimX );
  printf( "dimY == %lld\n", header.nDimY );
  printf( "dx  == %g\n", header.dx );
  printf( "dkx == %g\n", header.dkx );
  printf( "dy  == %g\n", header.dy );
  printf( "dky == %g\n", header.dky );

  double runtime_size = (header.nDimX*header.nDimY*sizeof(double)+header.nDimX*(header.nDimY/2 + 1)*sizeof(fftw_complex))/pow(1024,3);
  double end_size = header.nDimX*header.nDimY*sizeof(double)/pow(1024,3);
  printf("Memory: Runtime %f GB, End %f GB\n", runtime_size, end_size);

  Fourier::CNoise<Fourier::rft_2d,2> noise( header );
  noise.color_noise_custom(2,chunk_dkx);
  noise.save("Noise.bin");

  fftw_cleanup_threads();
  return EXIT_SUCCESS;
}
