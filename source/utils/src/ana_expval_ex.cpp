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
#include <omp.h>
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"

using namespace std;

template<int dim, class T1>
void Expval_Momentum( const generic_header &header, const vector<CPoint<dim>> &ks, const vector<double> &thresholds, fftw_complex *wf, vector<CPoint<dim>> &all_retvals )
{
  all_retvals.clear();
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  assert( ks.size() == thresholds.size() );

  T1 ft(header);

  fftw_complex *in = ft.Getp2In();
  memcpy( in, wf, sizeof(fftw_complex)*Nges );
  ft.ft(-1);

  for ( int i=0; i<ks.size(); i++ )
  {
    CPoint<dim> k_ref = ks[i];
    CPoint<dim> k_diff;
    CPoint<dim> retval;
    double k_eps = thresholds[i];

    double *tmp=nullptr;
    double *tmp2=nullptr;
    double PN=0;
    #pragma omp parallel
    {
      const int nthreads = omp_get_num_threads();
      const int ithread = omp_get_thread_num();

      double den;
      CPoint<dim> k;
      fftw_complex *Psi = ft.Getp2In();

      #pragma omp single
      {
        tmp = new double[dim*nthreads];
        for (int i=0; i<(dim*nthreads); i++) tmp[i] = 0;
        tmp2 = new double[nthreads];
        for (int i=0; i<nthreads; i++) tmp2[i] = 0;
      }

      #pragma omp for
      for ( int l=0; l<Nges; l++ )
      {
        k = ft.Get_k(l);
        k_diff = k_ref - k;
        if ( sqrt(fabs(k_diff*k_diff)) > k_eps ) continue;
        den = (Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
        for (int i=0; i<dim; i++ )
          tmp[ithread*dim+i] += k[i]*den;
        tmp2[ithread] += den;
      }

      #pragma omp for
      for (int i=0; i<dim; i++)
        for (int t=0; t<nthreads; t++)
          retval[i] += tmp[dim*t + i];

      #pragma omp single
      for (int t=0; t<nthreads; t++)
        PN += tmp2[t];
    }

    for ( int i=0; i<dim; i++ )
      retval[i] /= PN;

    all_retvals.push_back(retval);
    delete [] tmp;
    delete [] tmp2;
  }
}

template<int dim, class T1>
void Expval_2nd_moment_Momentum(  const generic_header &header, const vector<CPoint<dim>> &moms, const vector<CPoint<dim>> &ks, const vector<double> &thresholds, fftw_complex *wf, vector<CPoint<dim>> &all_retvals  )
{
  all_retvals.clear();
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  assert( ks.size() == thresholds.size() );

  T1 ft(header);

  fftw_complex *in = ft.Getp2In();
  memcpy( in, wf, sizeof(fftw_complex)*Nges );
  ft.ft(-1);

  for ( int i=0; i<ks.size(); i++ )
  {
    CPoint<dim> k_ref = ks[i];
    CPoint<dim> k_diff;
    CPoint<dim> retval;
    CPoint<dim> mu = moms[i];
    double k_eps = thresholds[i];

    double *tmp;
    double *tmp2;
    double PN=0;
    #pragma omp parallel
    {
      const int nthreads = omp_get_num_threads();
      const int ithread = omp_get_thread_num();

      double den;
      CPoint<dim> k;
      fftw_complex *Psi = ft.Getp2In();

      #pragma omp single
      {
        tmp = new double[dim*nthreads];
        for (int i=0; i<(dim*nthreads); i++) tmp[i] = 0;
        tmp2 = new double[nthreads];
        for (int i=0; i<nthreads; i++) tmp2[i] = 0;
      }

      #pragma omp for
      for ( int l=0; l<Nges; l++ )
      {
        k = ft.Get_k(l);
        k_diff = k_ref - k;
        if ( sqrt(fabs(k_diff*k_diff)) > k_eps ) continue;

        den = (Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
        for (int i=0; i<dim; i++ )
          tmp[ithread*dim+i] += (k[i]-mu[i])*(k[i]-mu[i])*den;
        tmp2[ithread] += den;
      }

      #pragma omp for
      for (int i=0; i<dim; i++)
        for (int t=0; t<nthreads; t++)
          retval[i] += tmp[dim*t + i];

      #pragma omp single
      for (int t=0; t<nthreads; t++)
        PN += tmp2[t];
    }

    for ( int i=0; i<dim; i++ )
      retval[i] /= PN;

    all_retvals.push_back(retval);
    delete [] tmp;
    delete [] tmp2;
  }
}

int main(int argc, char *argv[])
{
  FILE *fh = nullptr;

  double PN;

  int no_of_threads = 4;
  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );
  omp_set_num_threads( no_of_threads );

  generic_header header;
  bzero( &header, sizeof(generic_header) );

  if ( argc > 1 )
  {
    fh = fopen( argv[1], "r" );
    if ( fh == nullptr )
    {
      printf( "Could not open file %s.\n", argv[1] );
      exit(0);
    }

    fread( &header, sizeof(generic_header), 1, fh );

    // printf( "### %s\n", argv[1] );
    // printf( "# nDims    == %lld\n", header.nDims );
    // printf( "# nDimX    == %lld\n", header.nDimX );
    // printf( "# nDimY    == %lld\n", header.nDimY );
    // printf( "# nDimZ    == %lld\n", header.nDimZ );
    // printf( "# nDatatyp == %lld\n", header.nDatatyp );
    // printf( "# bAtom    == %d\n", header.bAtom );
    // printf( "# bComplex == %d\n", header.bComplex );
    // printf( "# t        == %g\n", header.t );
    // printf( "# dt       == %g\n", header.dt );
    // printf( "# xMin     == %g\n", header.xMin );
    // printf( "# xMax     == %g\n", header.xMax );
    // printf( "# yMin     == %g\n", header.yMin );
    // printf( "# yMax     == %g\n", header.yMax );
    // printf( "# zMin     == %g\n", header.zMin );
    // printf( "# zMax     == %g\n", header.zMax );
    // printf( "# dx       == %g\n", header.dx );
    // printf( "# dy       == %g\n", header.dy );
    // printf( "# dz       == %g\n", header.dz );
    // printf( "# dkx      == %g\n", header.dkx );
    // printf( "# dky      == %g\n", header.dky );
    // printf( "# dkz      == %g\n", header.dkz );
  }

  if ( header.nDims == 1 )
  {
    header.nDimY=1;
    header.nDimZ=1;
  };
  if ( header.nDims == 2 )
  {
    header.nDimZ=1;
  };

  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;
  fftw_complex *wave_function = (fftw_complex *)fftw_malloc( sizeof(fftw_complex) *  Nges );

  fread( (void *)wave_function, sizeof(fftw_complex), Nges, fh );
  fclose(fh);


  //1D Stuff
  if ( header.nDims == 1 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    std::vector<CPoint<1>> ks;
    std::vector<double> thresholds;

    CPoint<1> pt1 = {0};
    CPoint<1> pt2 = {16};

    ks.push_back(pt1);
    ks.push_back(pt2);

    thresholds.push_back(5.0);
    thresholds.push_back(5.0);

    vector<CPoint<1>> first_moment_P;
    vector<CPoint<1>> second_moment_P;

    Expval_Momentum<1,Fourier::cft_1d>( header, ks, thresholds, wave_function, first_moment_P );
    Expval_2nd_moment_Momentum<1,Fourier::cft_1d>( header, first_moment_P, ks, thresholds, wave_function, second_moment_P );

    for ( auto i : first_moment_P )
    {
      cout << "\t" << i[0];
    }
    for ( auto i : second_moment_P )
    {
      cout << "\t" << sqrt(i[0]);
    }
    cout << endl;
  }

  fftw_free(wave_function);
  return EXIT_SUCCESS;
}
