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
#include <ctime>
#include <fstream>
#include <vector>
#include <sys/time.h>
#include <omp.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "fftw3.h"
#include "my_structs.h"
#include "anyoption.h"
#include "eigenfunctions_HO.h"

using namespace std;

struct my_params
{
  int Nx;
  int Ny;
  int Nz;
  double Lx;
  double Ly;
  double Lz;
  fftw_complex *wf;  // wave function
  generic_header header;
};

// v entspricht a_ijk
double my_f ( const gsl_vector *v, void *params )
{
  double retval=0, x, y, z;
  my_params *p = reinterpret_cast<my_params *>(params);

  fftw_complex *wf = p->wf;  // wf[ijk][0] real teil, wf[ijk][1] imag teil

  long long shift_x = p->header.nDimX/2;
  long long shift_y = p->header.nDimX/2;
  long long shift_z = p->header.nDimX/2;
  double dx = p->header.dx;
  double dy = p->header.dy;
  double dz = p->header.dz;


  for ( long long i=0; i<p->Nx; i++ )
  {
    x = double(i-shift_x)*dx;
    for ( long long j=0; j<p->Ny; j++ )
    {
      y = double(j-shift_y)*dy;

    }
  }
  /*
    double x, y;
    double *p = (double *)params;

    x = gsl_vector_get(v, 0);
    y = gsl_vector_get(v, 1);

    return p[2] * (x - p[0]) * (x - p[0]) +
             p[3] * (y - p[1]) * (y - p[1]) + p[4];*/
  return retval*retval;
}

// 3D Fall
void display_3D( fftw_complex *field, generic_header *header, int opt, const int i0, const int j0, const int k0 )
{
  /*
  const int dimX = header->nDimX;
  const int dimY = header->nDimY;
  const int dimZ = header->nDimZ;
  const double dx = header->dx;
  const double dy = header->dy;
  const double dz = header->dz;
  const double x0 = header->xMin - header->dx;
  const double y0 = header->yMin - header->dy;
  const double z0 = header->zMin - header->dz;

  int i, j, k, ijk;
  double x, y, z;
  */

  my_params params;
  params.header = *header;
  params.wf = field;
  params.Nx = 10;
  params.Ny = 10;
  params.Nz = 10;
  params.Lx = 1;
  params.Ly = 1;
  params.Lz = 1;


  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = nullptr;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  long long tot = 2*params.Nx*params.Ny*params.Nz;

  /* Starting point */
  gsl_vector *x = gsl_vector_alloc (tot);
  gsl_vector_set_all (x, 1.0);

  /* Set initial step sizes to 1 */
  gsl_vector *ss = gsl_vector_alloc (tot);
  gsl_vector_set_all (ss, 1.0);

  /* Initialize method and iterate */
  minex_func.n = params.Nx*params.Ny*params.Nz;
  minex_func.f = my_f;
  minex_func.params = (void *)(&params);

  s = gsl_multimin_fminimizer_alloc (T, tot);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status)
      break;

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1e-2);

    if (status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }

    /*
          printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
                  iter,
                  gsl_vector_get (s->x, 0),
                  gsl_vector_get (s->x, 1),
                  s->fval, size);
    */
  }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

// 2D Fall
void display_2D( fftw_complex *field, generic_header *header, int opt, const int i0, const int j0 )
{
  const int    dimX = header->nDimX;
  const int    dimY = header->nDimY;
  const double dx   = header->dx;
  const double dy   = header->dy;
  const double x0   = header->xMin - header->dx;
  const double y0   = header->yMin - header->dy;

  int i, j, ij;
  double x, y;


}

// 1D Fall
void display_1D( fftw_complex *field, generic_header *header )
{
}

int main(int argc, char *argv[])
{
  int i0=-1, j0=-1, k0=-1, options=0;

  string filename;

  AnyOption *opt = new AnyOption();
  opt->noPOSIX();
  //opt->setVerbose();
  //opt->autoUsagePrint(true);

  opt->addUsage( "" );
  opt->addUsage( "Usage: gpo2 [options] filename" );
  opt->addUsage( "" );
  opt->addUsage( " -h --help	Prints this help " );
  opt->addUsage( "" );

  opt->setFlag(  "help", 'h' );

  opt->processCommandArgs( argc, argv );

  if ( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) opt->printUsage();

  if ( opt->getArgc() != 0 ) filename = opt->getArgv(0);
  else opt->printUsage();

  delete opt;

  generic_header header;

  ifstream in;
  in.open(filename);
  if ( !in.is_open() ) return EXIT_FAILURE;
  in.read( (char *)&header, sizeof(generic_header) );

  /*
    int no_of_threads = 4;
    char* envstr = getenv( "MY_NO_OF_THREADS" );
    if( envstr != NULL ) no_of_threads = atoi( envstr );
    omp_set_num_threads( no_of_threads );
  */

  printf( "### %s\n", filename.c_str() );
  printf( "# nDims    == %lld\n", header.nDims );
  printf( "# nDimX    == %lld\n", header.nDimX );
  printf( "# nDimY    == %lld\n", header.nDimY );
  printf( "# nDimZ    == %lld\n", header.nDimZ );
  printf( "# nDatatyp == %lld\n", header.nDatatyp );
  printf( "# bAtom    == %d\n", header.bAtom );
  printf( "# bComplex == %d\n", header.bComplex );
  printf( "# t        == %g\n", header.t );
  printf( "# dt       == %g\n", header.dt );
  printf( "# xMin     == %g\n", header.xMin );
  printf( "# xMax     == %g\n", header.xMax );
  printf( "# yMin     == %g\n", header.yMin );
  printf( "# yMax     == %g\n", header.yMax );
  printf( "# zMin     == %g\n", header.zMin );
  printf( "# zMax     == %g\n", header.zMax );
  printf( "# dx       == %g\n", header.dx );
  printf( "# dy       == %g\n", header.dy );
  printf( "# dz       == %g\n", header.dz );
  printf( "# dkx      == %g\n", header.dkx );
  printf( "# dky      == %g\n", header.dky );
  printf( "# dkz      == %g\n", header.dkz );
  printf( "# i0       == %d\n", i0 );
  printf( "# j0       == %d\n", j0 );
  printf( "# k0       == %d\n", k0 );

  size_t total_no_bytes;

  /*
  // 1D Stuff
  if( header.nDims == 1 && header.bComplex == 1 )
  {
    total_no_bytes = sizeof(fftw_complex)*header.nDimX;
    fftw_complex* field = fftw_alloc_complex( header.nDimX );
    in.read( (char*)field, total_no_bytes );

    display_1D( field, &header );

    fftw_free(field);
  }
  */
  // 2D Stuff
  if ( header.nDims == 2 && header.bComplex == 1 )
  {
    total_no_bytes = sizeof(fftw_complex)*header.nDimX*header.nDimY;

    fftw_complex *field = fftw_alloc_complex(header.nDimX*header.nDimY);
    in.read( (char *)field, total_no_bytes );

    display_2D( field, &header, options, i0, j0 );

    fftw_free(field);
  }

  if ( header.nDims == 3 && header.bComplex == 1  )
  {
    total_no_bytes = sizeof(fftw_complex)*header.nDimX*header.nDimY*header.nDimZ;

    fftw_complex *field = fftw_alloc_complex(header.nDimX*header.nDimY*header.nDimZ);
    in.read( (char *)field, total_no_bytes );

    display_3D( field, &header, options, i0, j0, k0 );

    fftw_free(field);
  }

}
