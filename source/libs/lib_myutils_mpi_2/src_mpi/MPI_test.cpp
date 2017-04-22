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

#include <fftw3-mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>


using namespace std;
int main(int argc, char **argv)
{
  const ptrdiff_t N0 = 4096, N1 = 4096;
  fftw_plan plan;
  fftw_complex *data;
  ptrdiff_t alloc_local, local_n0, local_0_start, i, j;
  //std::cout << "Test1" << std::endl;
  MPI_Init(&argc, &argv);

  fftw_mpi_init();
  //cout << "Test3" << endl;

  /* get local data size and allocate */
  alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
                                       &local_n0, &local_0_start);
  data = fftw_alloc_complex(alloc_local);
  //cout << "Test4" << endl;

  /* create plan for in-place forward DFT */
  plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
                              FFTW_FORWARD, FFTW_ESTIMATE);
  //cout << "Test5" << endl;

  /* initialize data to some function my_function(x,y) */
  for (i = 0; i < local_n0; ++i) for (j = 0; j < N1; ++j)
    {
      data[i*N1 + j][0] = 0.0;//my_function(local_0_start + i, j);
      data[i*N1 + j][1] = 0.0;
    }
  /* compute transforms, in-place, as many times as desired */

  fftw_execute(plan);
  fftw_destroy_plan(plan);
  //cout << "Test6" << endl;

  MPI_Finalize();
}