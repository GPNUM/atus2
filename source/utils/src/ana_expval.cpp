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
#include <omp.h>
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"
#include "fftw3.h"

using namespace std;
using namespace Fourier;

bool Read( const char *filename, const generic_header &header, fftw_complex *field )
{
  FILE *fh = fopen( filename, "r" );

  if ( fh == nullptr ) return false;

  fseek( fh, sizeof(generic_header), SEEK_SET );
  size_t Nges = header.nDimX * header.nDimY * header.nDimZ;

  fread( (void *)field, sizeof(fftw_complex), Nges, fh );
  fclose(fh);
  return true;
}

template<int dim>
double Phase( const generic_header &header, cft_base<dim> *ft )
{
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  double tmp = 0;
  double retval = 0;

  fftw_complex *Psi = ft->Getp2In();

  for ( int l=0; l<Nges; l++ )
  {
    tmp += atan2(Psi[l][1],Psi[l][0])*(Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
    retval += (Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
  }
  tmp = tmp / retval;

  return tmp;
}

template<int dim>
double Particle_Number( const generic_header &header, cft_base<dim> *ft )
{
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  double retval=0;

  fftw_complex *Psi = ft->Getp2In();

  #pragma omp parallel for reduction(+:retval)
  for ( int l=0; l<Nges; l++ )
    retval += (Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);

  return retval;
}

template<int dim>
void Expval_Position( const generic_header &header, cft_base<dim> *ft, CPoint<dim> &retval )
{
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;
  double *tmp;
  #pragma omp parallel
  {
    const int nthreads = omp_get_num_threads();
    const int ithread = omp_get_thread_num();

    double den;
    CPoint<dim> x;
    fftw_complex *Psi = ft->Getp2In();

    #pragma omp single
    {
      tmp = new double[dim*nthreads];
      for (int i=0; i<(dim*nthreads); i++) tmp[i] = 0;
    }

    #pragma omp for
    for ( int l=0; l<Nges; l++ )
    {
      x = ft->Get_x(l);
      den = (Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
      for (int i=0; i<dim; i++ )
        tmp[ithread*dim+i] += x[i]*den;
    }

    #pragma omp for
    for (int i=0; i<dim; i++)
    {
      for (int t=0; t<nthreads; t++)
        retval[i] += tmp[dim*t + i];
    }
  }
  delete [] tmp;
}

template<int dim>
void Expval_2nd_moment( const generic_header &header, cft_base<dim> *ft, CPoint<dim> mu, CPoint<dim> &retval )
{
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;
  double *tmp;
  #pragma omp parallel
  {
    const int nthreads = omp_get_num_threads();
    const int ithread = omp_get_thread_num();

    double den;
    CPoint<dim> x;
    fftw_complex *Psi = ft->Getp2In();

    #pragma omp single
    {
      tmp = new double[dim*nthreads];
      for (int i=0; i<(dim*nthreads); i++) tmp[i] = 0;
    }

    #pragma omp for
    for ( int l=0; l<Nges; l++ )
    {
      x = ft->Get_x(l);
      den = (Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
      for (int i=0; i<dim; i++ )
        tmp[ithread*dim+i] += (x[i]-mu[i])*(x[i]-mu[i])*den;
    }

    #pragma omp for
    for (int i=0; i<dim; i++)
    {
      for (int t=0; t<nthreads; t++)
        retval[i] += tmp[dim*t + i];
    }
  }
  delete [] tmp;
}

template<int dim>
void Expval_Momentum( const generic_header &header, cft_base<dim> *ft, CPoint<dim> &retval )
{
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;
  double *tmp;
  #pragma omp parallel
  {
    const int nthreads = omp_get_num_threads();
    const int ithread = omp_get_thread_num();

    double den;
    CPoint<dim> k;
    fftw_complex *Psi = ft->Getp2In();

    #pragma omp single
    {
      tmp = new double[dim*nthreads];
      for (int i=0; i<(dim*nthreads); i++) tmp[i] = 0;
    }

    #pragma omp for
    for ( int l=0; l<Nges; l++ )
    {
      k = ft->Get_k(l);
      den = (Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
      for (int i=0; i<dim; i++ )
        tmp[ithread*dim+i] += k[i]*den;
    }

    #pragma omp for
    for (int i=0; i<dim; i++)
      for (int t=0; t<nthreads; t++)
        retval[i] += tmp[dim*t + i];
  }
  delete [] tmp;
}

template<int dim>
void Expval_2nd_moment_Momentum( const generic_header &header, cft_base<dim> *ft, CPoint<dim> mu, CPoint<dim> &retval )
{
  const long long Nges = header.nDimX * header.nDimY * header.nDimZ;
  double *tmp;
  #pragma omp parallel
  {
    const int nthreads = omp_get_num_threads();
    const int ithread = omp_get_thread_num();

    double den;
    CPoint<dim> k;
    fftw_complex *Psi = ft->Getp2In();

    #pragma omp single
    {
      tmp = new double[dim*nthreads];
      for (int i=0; i<(dim*nthreads); i++) tmp[i] = 0;
    }

    #pragma omp for
    for ( int l=0; l<Nges; l++ )
    {
      k = ft->Get_k(l);
      den = (Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
      for (int i=0; i<dim; i++ )
        tmp[ithread*dim+i] += (k[i]-mu[i])*(k[i]-mu[i])*den;
    }

    #pragma omp for
    for (int i=0; i<dim; i++)
      for (int t=0; t<nthreads; t++)
        retval[i] += tmp[dim*t + i];
  }
  delete [] tmp;
}

int main(int argc, char *argv[])
{
  FILE *fh = nullptr;

  double PN;

  int no_of_threads = 4;
  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );
  omp_set_num_threads( no_of_threads );

  generic_header header = {};

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
    fclose(fh);
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

  //1D Stuff
  if ( header.nDims == 1 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    CPoint<1> first_moment;
    CPoint<1> second_moment;
    CPoint<1> first_moment_P;
    CPoint<1> second_moment_P;

    Fourier::cft_1d *ft1 = new Fourier::cft_1d(header);
    Read( argv[1], header, ft1->Getp2In() );

    double phase = Phase( header, ft1 );

    PN = header.dx*Particle_Number( header, ft1 );
    Expval_Position( header, ft1, first_moment );
    first_moment *= (header.dx) / PN;
    Expval_2nd_moment( header, ft1, first_moment, second_moment );
    second_moment *= (header.dx) / PN;

    ft1->ft(-1);
    Expval_Momentum( header, ft1, first_moment_P );
    first_moment_P *= (header.dkx) / PN;
    Expval_2nd_moment_Momentum( header, ft1, first_moment_P, second_moment_P );
    second_moment_P *= (header.dkx) / PN;

    delete ft1;

    // cout << "N = " << PN << endl;
    // cout << "<X> = " << first_moment << endl;
    // cout << "<(X-<X>)^2> = " << second_moment << endl;
    // cout << "<P> = " << first_moment_P << endl;
    // cout << "<(P-<P>)^2> = " << second_moment_P << endl;
    cout << phase << "\t" << PN << "\t" << first_moment << "\t" << second_moment << "\t" << first_moment_P << "\t" << second_moment_P << endl;
  }

  // 2D Stuff
  if ( header.nDims == 2 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    CPoint<2> first_moment;
    CPoint<2> second_moment;
    CPoint<2> first_moment_P;
    CPoint<2> second_moment_P;

    Fourier::cft_2d *ft2 = new Fourier::cft_2d(header);
    Read( argv[1], header, ft2->Getp2In() );

    PN = header.dx*header.dy*Particle_Number( header, ft2 );
    Expval_Position( header, ft2, first_moment );
    first_moment *= (header.dx*header.dy) / PN;
    Expval_2nd_moment( header, ft2, first_moment, second_moment );
    second_moment *= (header.dx*header.dy) / PN;

    ft2->ft(-1);
    Expval_Momentum( header, ft2, first_moment_P );
    first_moment_P *= (header.dkx*header.dky) / PN;
    Expval_2nd_moment_Momentum( header, ft2, first_moment_P, second_moment_P );
    second_moment_P *= (header.dkx*header.dky) / PN;

    delete ft2;

    // cout << "N = " << PN << endl;
    // cout << "<X> = " << first_moment << endl;
    // cout << "<(X-<X>)^2> = " << second_moment << endl;
    // cout << "<P> = " << first_moment_P << endl;
    // cout << "<(P-<P>)^2> = " << second_moment_P << endl;
    cout << PN << "\t" << first_moment << "\t" << second_moment << "\t" << first_moment_P << "\t" << second_moment_P << endl;
  }

  // 3D Stuff
  if ( header.nDims == 3 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    CPoint<3> first_moment;
    CPoint<3> second_moment;
    CPoint<3> first_moment_P;
    CPoint<3> second_moment_P;

    Fourier::cft_3d *ft3 = new Fourier::cft_3d(header);
    Read( argv[1], header, ft3->Getp2In() );

    PN = header.dx*header.dy*header.dz*Particle_Number( header, ft3 );
    Expval_Position( header, ft3, first_moment );
    first_moment *= (header.dx*header.dy*header.dz) / PN;
    Expval_2nd_moment( header, ft3, first_moment, second_moment );
    second_moment *= (header.dx*header.dy*header.dz) / PN;

    ft3->ft(-1);
    Expval_Momentum( header, ft3, first_moment_P );
    first_moment_P *= (header.dkx*header.dky*header.dkz) / PN;
    Expval_2nd_moment_Momentum( header, ft3, first_moment_P, second_moment_P );
    second_moment_P *= (header.dkx*header.dky*header.dkz) / PN;

    delete ft3;

    // cout << "N = " << PN << endl;
    // cout << "<X> = " << first_moment << endl;
    // cout << "<(X-<X>)^2> = " << second_moment << endl;
    // cout << "<P> = " << first_moment_P << endl;
    // cout << "<(P-<P>)^2> = " << second_moment_P << endl;
    cout << PN << "\t" << first_moment << "\t" << second_moment << "\t" << first_moment_P << "\t" << second_moment_P << endl;
  }

  return EXIT_SUCCESS;
}
