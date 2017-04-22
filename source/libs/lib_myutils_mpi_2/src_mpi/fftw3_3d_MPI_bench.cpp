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
#include <cmath>
#include "cft_3d_MPI.h"

using namespace std;

void fkt1( MPI::Fourier::cft_3d_MPI &ft )
{
  fftw_complex *data = ft.Get_p2_Data();
  ptrdiff_t N = ft.Get_loc_n_rs();

  CPoint<3> x;
  for ( ptrdiff_t l=0; l<N; l++ )
  {
    x = ft.Get_x( l );

    data[l][0] = exp(-0.5*(x*x));
    data[l][1] = 0.0;
  }
}

//--------------------------------------------------------------------------------
int main( int argc, char *argv[] )
{
  double ti, tf, tacc;

  const int repno = 3;
  const int N = 1024;
  char pname[MPI_MAX_PROCESSOR_NAME];

  int nprocs, myrank, len;

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Get_processor_name( pname, &len );
  fftw_mpi_init();

  printf( "(%d) - %s\n", myrank, pname );

  //*** Setze Dateikopf-Informationen *******************************************
  //*****************************************************************************
  generic_header header = {};

  header.nself    = sizeof(generic_header);
  header.nDatatyp = sizeof(fftw_complex);
  header.nDims    = 3;
  header.nDimX    = N/2;
  header.nDimY    = N;
  header.nDimZ    = N;
  header.bAtom    = 1;
  header.bComplex = 1;
  header.t        = 0.0;
  header.xMin     = -10.0;
  header.xMax     = -header.xMin;
  header.yMin     = -10.0;
  header.yMax     = -header.xMin;
  header.zMin     = -10.0;
  header.zMax     = -header.zMin;
  header.dx       = fabs( header.xMax-header.xMin )/double(header.nDimX);
  header.dkx      = 2.0*M_PI/fabs(header.xMax-header.xMin);
  header.dy       = fabs( header.yMax-header.yMin )/double(header.nDimY);
  header.dky      = 2.0*M_PI/fabs(header.yMax-header.yMin);
  header.dz       = fabs( header.zMax-header.zMin )/double(header.nDimZ);
  header.dkz      = 2.0*M_PI/fabs(header.zMax-header.zMin);
  header.nself_and_data = header.nself + (header.nDimX + header.nDimY + header.nDimZ)*header.nDatatyp;

  printf( "dx  == %g\n", header.dx );
  printf( "dy  == %g\n", header.dy );
  printf( "dkx == %g\n", header.dkx );
  printf( "dky == %g\n", header.dky );

  ti = MPI_Wtime();
  MPI::Fourier::cft_3d_MPI ft1 = MPI::Fourier::cft_3d_MPI( &header );
  tf = MPI_Wtime();
  fkt1( ft1 );

  printf( "(%d) init time %.10g\n", myrank, tf-ti );

  tacc = 0.0;
  for ( int l=0; l<repno; l++ )
  {
    ti = MPI_Wtime();
    ft1.ft(-1);
    ft1.ft(1);
    tf = MPI_Wtime();
    tacc += tf - ti;
    printf( "(%d) (%d) tf-ti == %.10g\n", myrank, l, tf-ti );
  }

  printf( "(%d) avg %.10g\n", myrank, tacc/double(repno) );

  MPI_Finalize();
  return 0;
}
