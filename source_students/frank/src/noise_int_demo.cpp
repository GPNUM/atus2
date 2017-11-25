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

using namespace std;

int main(int argc, char *argv[])
{
  const int no_chunks = 2;
  const int fak = 4;

  generic_header header = {};

  header.nself    = sizeof(generic_header);
  header.nDatatyp = sizeof(fftw_complex);
  header.nDims    = 2;
  header.nDimX    = 4096; // Zeit
  header.nDimY    = 2048; // Ort
  header.nDimZ    = 1;
  header.xMin     = -100; // t_i
  header.xMax     = 100; // t_f
  header.dx       = fabs( header.xMax-header.xMin )/double(header.nDimX);
  header.dkx      = 2.0*M_PI/fabs(header.xMax-header.xMin);
  header.yMin     = -5;
  header.yMax     = -header.yMin;
  header.dy       = fabs( header.yMax-header.yMin )/double(header.nDimY);
  header.dky      = 2.0*M_PI/fabs(header.yMax-header.yMin);
  header.nself_and_data = header.nself + (header.nDimX*header.nDimY*header.nDimZ)*header.nDatatyp;

  // header for the chunk
  const double chunk_len =  fabs( header.xMax-header.xMin ) / double(no_chunks);
  generic_header header2 = header;
  header2.nDimX = header.nDimX/no_chunks;
  header2.xMax = header.xMin + chunk_len;
  header2.dkx = 2.0*M_PI/fabs( header2.xMax-header2.xMin );

  // header for the interpolated data
  generic_header header3 = header2;
  header3.nDimX = fak*header2.nDimX;
  header3.nDimY = fak*header2.nDimY;
  header3.dx = header3.dx/double(fak);
  header3.dy = header3.dy/double(fak);
  header3.nself_and_data = header.nself + (header3.nDimX*header3.nDimY*header3.nDimZ)*header.nDatatyp;

  printf( "dimX == %lld\n", header.nDimX );
  printf( "dimY == %lld\n", header.nDimY );
  printf( "dimX2 == %lld\n", header2.nDimX );
  printf( "dimY2 == %lld\n", header2.nDimY );
  printf( "dimX3 == %lld\n", header3.nDimX );
  printf( "dimY3 == %lld\n", header3.nDimY );
  printf( "dx  == %g\n", header.dx );
  printf( "dkx == %g\n", header.dkx );
  printf( "dy  == %g\n", header.dy );
  printf( "dky == %g\n", header.dky );
  printf( "dx2  == %g\n", header2.dx );
  printf( "dkx2 == %g\n", header2.dkx );
  printf( "dy2  == %g\n", header2.dy );
  printf( "dky2 == %g\n", header2.dky );
  printf( "dx3  == %g\n", header3.dx );
  printf( "dkx3 == %g\n", header3.dkx );
  printf( "dy3  == %g\n", header3.dy );
  printf( "dky3 == %g\n", header3.dky );

  Fourier::cft_2d ft( header2, false, true );

  Fourier::CNoise<Fourier::rft_2d,2> noise( header );
  std::vector<double> mink = {header2.dkx, 0.0};
  noise.color_noise_custom(2,mink);
  //noise.white_noise();
  noise.save("noise.bin");

  Fourier::rft_2d chunkft( header2 );
  Fourier::rft_2d interpolft( header3 );

  // extract chunk
  double * noise_in = noise.Getp2InReal();
  double * chunk_in = chunkft.Getp2InReal();
  double * interpol_in = interpolft.Getp2InReal();
  fftw_complex * chunk_out = chunkft.Getp2Out();
  fftw_complex * interpolft_out = interpolft.Getp2Out();

  const int64_t Nx = chunkft.Get_Dim_X();
  const int64_t Nyred = chunkft.Get_red_Dim();
  const int64_t shifti = interpolft.Get_Dim_X() - chunkft.Get_Dim_X();
  const int64_t Nynew = interpolft.Get_red_Dim();

  for( int i=0; i<no_chunks; i++ )
  {
    // copy chunk from the coarse grid
    int64_t offset = i*chunkft.Get_Dim_RS();
    const void * ptr = reinterpret_cast<void*>(&noise_in[offset]);
    memcpy( reinterpret_cast<void*>(chunk_in), ptr, sizeof(double) * chunkft.Get_Dim_RS() );
    chunkft.save( "chunk_" + std::to_string(i) + ".bin" );
    chunkft.ft(-1);

    memset( reinterpret_cast<void*>(interpolft_out), 0, sizeof(fftw_complex)*interpolft.Get_Dim_FS() );
    #pragma omp parallel for collapse(2)
    for( int i=0; i<Nx/2; i++ )
    {
      for( int j=0; j<Nyred; j++ )
      {
        interpolft_out[j+i*Nynew][0] = chunk_out[j+i*Nyred][0];
        interpolft_out[j+i*Nynew][1] = chunk_out[j+i*Nyred][1];
      }
    }

    #pragma omp parallel for collapse(2)
    for( int i=Nx/2; i<Nx; i++ )
    {
      for( int j=0; j<Nyred; j++ )
      {
        interpolft_out[j+(i+shifti)*Nynew][0] = chunk_out[j+i*Nyred][0];
        interpolft_out[j+(i+shifti)*Nynew][1] = chunk_out[j+i*Nyred][1];
      }
    }
    interpolft.save( "fichunk_" + std::to_string(i) + ".bin", false );
    interpolft.ft(1);
    interpolft.save( "ichunk_" + std::to_string(i) + ".bin" );

    double maxval1 = 0;
    #pragma omp parallel for reduction(max:maxval1)
    for( int64_t l=0; l<chunkft.Get_Dim_RS(); l++ )
    {
      if( fabs(chunk_in[l]) > maxval1 ) maxval1 = fabs(chunk_in[l]);
    }

    double maxval2 = 0;
    #pragma omp parallel for reduction(max:maxval2)
    for( int64_t l=0; l<interpolft.Get_Dim_RS(); l++ )
    {
      if( fabs(interpol_in[l]) > maxval2 ) maxval2 = fabs(interpol_in[l]);
    }

    std::cout << "max val chunk = " << maxval1 << std::endl;
    std::cout << "max val interpol = " << maxval2 << std::endl;
  }

  fftw_cleanup_threads();
return EXIT_SUCCESS;
}
