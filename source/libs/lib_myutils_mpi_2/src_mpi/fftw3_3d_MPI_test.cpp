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
#include "fftw3-mpi.h"
#include "my_structs.h"
#include "cft_3d_MPI.h"

using namespace std;

void fkt( MPI::Fourier::cft_3d_MPI &ft )
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
  const int N = 128;
  char pname[MPI_MAX_PROCESSOR_NAME];

  int nprocs, myrank, len;

  MPI_Init( &argc, &argv );
  fftw_mpi_init();

  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
  MPI_Get_processor_name( pname, &len );

  printf( "(%d) - %s\n", myrank, pname );

  //*** Setze Dateikopf-Informationen *******************************************
  //*****************************************************************************
  generic_header header = {};

  header.nself    = sizeof(generic_header);
  header.nDatatyp = sizeof(fftw_complex);
  header.nDims    = 3;
  header.nDimX    = N;
  header.nDimY    = N;
  header.nDimZ    = N;
  header.bAtom    = 1;
  header.bComplex = 1;
  header.t        = 0.0;
  header.xMin     = -100.0;
  header.xMax     = -header.xMin;
  header.yMin     = -100.0;
  header.yMax     = -header.xMin;
  header.zMin     = -100.0;
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
  printf( "dz  == %g\n", header.dz );
  printf( "dkx == %g\n", header.dkx );
  printf( "dky == %g\n", header.dky );
  printf( "dkz == %g\n", header.dkz );

  MPI::Fourier::cft_3d_MPI *ft1 = new MPI::Fourier::cft_3d_MPI( &header );
  fkt( *ft1 );

  ft1->Write_File( "gauss.bin" );
  ft1->D_x();
  ft1->Write_File( "D_x.bin" );
  fkt( *ft1 );
  ft1->D_y();
  ft1->Write_File( "D_y.bin" );
  fkt( *ft1 );
  ft1->D_z();
  ft1->Write_File( "D_z.bin" );

  delete ft1;
  MPI_Finalize();
  return 0;
}
