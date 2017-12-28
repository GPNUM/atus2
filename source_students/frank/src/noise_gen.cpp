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

const int fak = 8;
const int no_of_chunks = fak*fak;

int main(int argc, char *argv[])
{
  if( argc < 2 )
    {
      cout << "Too few arguments\n"
           << argv[0] << " params.xml [inf_0.000_1.bin] [seed] [no_of_threads]" << endl;
      return EXIT_FAILURE;
    }

  string wavefun = "../inf_zero.bin";
  if (argc > 2) {
    wavefun = argv[2];
  }

  int64_t seed = time(nullptr);
  if (argc > 3) {
    cout << "User provided Seed!" << endl;
    seed = stol(argv[3]);
  }
  cout << "Seed: " << seed << endl;

  int no_of_threads = 4;
  char* envstr = getenv( "MY_NO_OF_THREADS" );
  if( envstr != NULL ) no_of_threads = atoi( envstr );
  if (argc > 4) {
    cout << "User provided Number of Threads!" << endl;
    no_of_threads = stoi(argv[4]);
  }
  cout << "Number of threads: " << no_of_threads << endl;

  fftw_init_threads();
  fftw_plan_with_nthreads( no_of_threads );
  omp_set_num_threads( no_of_threads );

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
  int exponent = ceil(log2(duration/dt));
  cout << "Exponent: " << exponent << endl;
  duration = pow(2.0, ceil(log2(duration/dt)))*dt;
  cout << "New Duration: " << duration << endl;
  cout << "Lower Duration: " << pow(2.0, exponent-1)*dt << endl;

  generic_header header = {};
  {
    ifstream fheader(wavefun, ifstream::binary );
    if (fheader.fail()) {
      cout << "File not found: " << wavefun << endl;
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

  // header for the chunk
  const double chunk_len =  fabs( header.xMax-header.xMin ) / double(no_of_chunks);
  const double chunk_dkx =  2.0*M_PI/fabs(chunk_len);

  const double runtime_size = (header.nDimX*header.nDimY*sizeof(double)+header.nDimX*(header.nDimY/2 + 1)*sizeof(fftw_complex))/pow(1024,3);
  const double end_size = header.nDimX*header.nDimY*sizeof(double)/pow(1024,3);
  cout << "Memory: Runtime "<< runtime_size << "GB, End " << end_size << "GB\n" << endl;

  const double chunk_x = fabs(header.xMax-header.xMin)/fabs(no_of_chunks)/10.0;
  const double chunk_y = fabs(header.yMax-header.yMin)/10.0;
  cout << "chunk_x: " << chunk_x << ", x " << chunk_x*header.dx << endl;
  cout << "chunk_y: " << chunk_y << ", y " << chunk_y*header.dy << endl;
//  vector<double> corr_length = {1000,1000};
  vector<double> corr_length = {chunk_x,chunk_y};

  cout << "corr_length_x\t" << corr_length[0] << endl;
  cout << "corr_length_y\t" << corr_length[1] << endl;

  vector<double> mink = {chunk_dkx, 0.0};

  Fourier::CNoise<Fourier::rft_2d,2> noise( header, seed );
  noise.color_noise_exp(2.0, mink, corr_length);
  noise.save("Noise.bin");

  ofstream log("noise_gen.log");
  log << "Seed\t" << seed << endl;
  log << "Threads\t" << no_of_threads << endl;
  log << "Expansion_X\t" << fak << endl;
  log << "Expansion_Y\t" << fak << endl;
  log << "Chunks\t" << no_of_chunks << endl;
  log << "corr_length_x\t" << corr_length[0] << endl;
  log << "corr_length_y\t" << corr_length[1] << endl;
  fftw_cleanup_threads();
  return EXIT_SUCCESS;
}
