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
#include <fstream>
#include "cxxopts.hpp"
#include "fftw3.h"
#include "my_structs.h"

extern double sign( const double val );

enum
{
  I_CONST = 0x01,
  J_CONST = 0x02,
  K_CONST = 0x04,
  RT = 0x08,
  IT = 0x10,
  PH = 0x20
};

using namespace std;

double mode1( double rt, double it )
{
  return rt;
}

double mode2( double rt, double it )
{
  return it;
}

double mode3( double rt, double it )
{
  return rt*rt+it*it;
}

double mode4( double rt, double it )
{
  return atan2(it,rt);
}

// 3D Fall
void display_3D( fftw_complex *field, generic_header *header, int opt, const int i0, const int j0, const int k0 )
{
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
  double (*fkt)(double, double);

  if ( opt & RT ) fkt = &mode1;
  else if ( opt & IT ) fkt = &mode2;
  else if ( opt & PH ) fkt = &mode4;
  else fkt = &mode3;

  if ( i0 < 0 || i0 > dimX ) opt &= ~I_CONST;
  if ( j0 < 0 || j0 > dimY ) opt &= ~J_CONST;
  if ( k0 < 0 || k0 > dimY ) opt &= ~K_CONST;

  if ( (opt & (I_CONST|J_CONST)) ==  (I_CONST|J_CONST) )
  {
    z = z0;
    for ( k=0; k<dimZ; k++ )
    {
      z += dz;
      ijk = k+dimZ*(j0+dimY*i0);
      printf( "%e\t%e\n", z, (*fkt)(field[ijk][0],field[ijk][1]) );
    }
    return;
  }

  if ( (opt & (I_CONST|K_CONST)) ==  (I_CONST|K_CONST) )
  {
    y = y0;
    for ( j=0; j<dimY; j++ )
    {
      y += dy;
      ijk = k0+dimZ*(j+dimY*i0);
      printf( "%e\t%e\n", y, (*fkt)(field[ijk][0],field[ijk][1]) );
    }
    return;
  }

  if ( (opt & (J_CONST|K_CONST)) ==  (J_CONST|K_CONST) )
  {
    x = x0;
    for ( i=0; i<dimX; i++ )
    {
      x += dx;
      ijk = k0+dimZ*(j0+dimY*i);
      printf( "%e\t%e\n", x, (*fkt)(field[ijk][0],field[ijk][1]) );
    }
    return;
  }

  if ( opt & K_CONST )
  {
    x = x0;
    y = y0;
    for ( i=0; i<dimX; i++ )
    {
      x += dx;
      y = y0;
      for ( j=0; j<dimY; j++ )
      {
        y += dy;
        ijk = k0+dimZ*(j+dimY*i);
        printf( "%e\t%e\t%e\n", x, y, (*fkt)(field[ijk][0],field[ijk][1]) );
      }
      printf( "\n" );
    }
    return;
  }

  if ( opt & J_CONST )
  {
    x = x0;
    z = z0;
    for ( i=0; i<dimX; i++ )
    {
      x += dx;
      z = z0;
      for ( k=0; k<dimZ; k++ )
      {
        z += dz;
        ijk = k+dimZ*(j0+dimY*i);
        printf( "%e\t%e\t%e\n", x, z, (*fkt)(field[ijk][0],field[ijk][1]) );
      }
      printf( "\n" );
    }
    return;
  }

  if ( opt & I_CONST )
  {
    y = y0;
    z = z0;
    for ( j=0; j<dimY; j++ )
    {
      y += dy;
      z = z0;
      for ( k=0; k<dimZ; k++ )
      {
        z += dz;
        ijk = k+dimZ*(j+dimY*i0);
        printf( "%e\t%e\t%e\n", y, z, (*fkt)(field[ijk][0],field[ijk][1]) );
      }
      printf( "\n" );
    }
    return;
  }
}

// 3D Fall
void display_3D( double *field, generic_header *header, int opt, const int i0, const int j0, const int k0 )
{
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

  if ( i0 < 0 || i0 > dimX ) opt &= ~I_CONST;
  if ( j0 < 0 || j0 > dimY ) opt &= ~J_CONST;
  if ( k0 < 0 || k0 > dimY ) opt &= ~K_CONST;

  if ( (opt & (I_CONST|J_CONST)) ==  (I_CONST|J_CONST) )
  {
    z = z0;
    for ( k=0; k<dimZ; k++ )
    {
      z += dz;
      ijk = k+dimZ*(j0+dimY*i0);
      printf( "%e\t%e\n", z, field[ijk] );
    }
    return;
  }

  if ( (opt & (I_CONST|K_CONST)) ==  (I_CONST|K_CONST) )
  {
    y = y0;
    for ( j=0; j<dimY; j++ )
    {
      y += dy;
      ijk = k0+dimZ*(j+dimY*i0);
      printf( "%e\t%e\n", y, field[ijk] );
    }
    return;
  }

  if ( (opt & (J_CONST|K_CONST)) ==  (J_CONST|K_CONST) )
  {
    x = x0;
    for ( i=0; i<dimX; i++ )
    {
      x += dx;
      ijk = k0+dimZ*(j0+dimY*i);
      printf( "%e\t%e\n", x, field[ijk] );
    }
    return;
  }

  if ( opt & K_CONST )
  {
    x = x0;
    y = y0;
    for ( i=0; i<dimX; i++ )
    {
      x += dx;
      y = y0;
      for ( j=0; j<dimY; j++ )
      {
        y += dy;
        ijk = k0+dimZ*(j+dimY*i);
        printf( "%e\t%e\t%e\n", x, y, field[ijk] );
      }
      printf( "\n" );
    }
    return;
  }

  if ( opt & J_CONST )
  {
    x = x0;
    z = z0;
    for ( i=0; i<dimX; i++ )
    {
      x += dx;
      z = z0;
      for ( k=0; k<dimZ; k++ )
      {
        z += dz;
        ijk = k+dimZ*(j0+dimY*i);
        printf( "%e\t%e\t%e\n", x, z, field[ijk] );
      }
      printf( "\n" );
    }
    return;
  }

  if ( opt & I_CONST )
  {
    y = y0;
    z = z0;
    for ( j=0; j<dimY; j++ )
    {
      y += dy;
      z = z0;
      for ( k=0; k<dimZ; k++ )
      {
        z += dz;
        ijk = k+dimZ*(j+dimY*i0);
        printf( "%e\t%e\t%e\n", y, z, field[ijk] );
      }
      printf( "\n" );
    }
    return;
  }
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
  double (*fkt)(double, double);
  if ( i0 < 0 || i0 > dimX ) opt &= ~I_CONST;
  if ( j0 < 0 || j0 > dimY ) opt &= ~J_CONST;

  if ( opt & RT ) fkt = &mode1;
  else if ( opt & IT ) fkt = &mode2;
  else if ( opt & PH ) fkt = &mode4;
  else fkt = &mode3;

  if ( opt & I_CONST )
  {
    x = x0 + double(i0+1)*dx;
    y = y0;
    for ( j=0; j<dimY; j++ )
    {
      y += dy;
      ij = j+dimY*i0;
      printf( "%e\t%e\n", y, (*fkt)(field[ij][0],field[ij][1]) );
    }
    return;
  }

  if ( opt & J_CONST )
  {
    x = x0;
    y = y0 + double(j0+1)*dy;
    for ( i=0; i<dimX; i++ )
    {
      x += dx;
      ij = j0+dimY*i;
      printf( "%e\t%e\n", x, (*fkt)(field[ij][0],field[ij][1]) );
    }
    return;
  }

  x = x0;
  for ( i=0; i<dimX; i++ )
  {
    x += dx;
    y  = y0;
    for ( j=0; j<dimY; j++ )
    {
      y += dy;
      ij = j+dimY*i;

      printf( "%g\t%g\t%.10e\n", x, y, (*fkt)(field[ij][0],field[ij][1]) );
    }
    printf( "\n" );
  }
}

// 2D Fall
void display_2D( double *field, generic_header *header, int opt, const int i0, const int j0 )
{
  const int dimX = header->nDimX;
  const int dimY = header->nDimY;
  const double dx = header->dx;
  const double dy = header->dy;
  const double x0 = header->xMin - header->dx;
  const double y0 = header->yMin - header->dy;

  int i, j, ij;
  double x, y;

  if ( i0 < 0 || i0 > dimX ) opt &= ~I_CONST;
  if ( j0 < 0 || j0 > dimY ) opt &= ~J_CONST;

  if ( opt & I_CONST )
  {
    x = x0 + double(i0+1)*dx;
    y = y0;
    for ( j=0; j<dimY; j++ )
    {
      y += dy;
      ij = j+dimY*i0;
      printf( "%e\t%e\n", y, field[ij] );
    }
    return;
  }

  if ( opt & J_CONST )
  {
    x = x0;
    y = y0 + double(j0+1)*dy;
    for ( i=0; i<dimX; i++ )
    {
      x += dx;
      ij = j0+dimY*i;
      printf( "%e\t%e\n", x, field[ij] );
    }
    return;
  }

  x = x0;
  for ( i=0; i<dimX; i++ )
  {
    x += dx;
    y  = y0;
    for ( j=0; j<dimY; j++ )
    {
      y += dy;
      ij = j+dimY*i;

      printf( "%g\t%g\t%.10e\n", x, y, field[ij] );
    }
    printf( "\n" );
  }
}

// 1D Fall
void display_1D( fftw_complex *field, generic_header *header, int opt )
{
  const int    dimX = header->nDimX;
  const double dx   = header->dx;
  double x   = header->xMin - header->dx;

  double (*fkt)(double, double);

  if ( opt & RT ) fkt = &mode1;
  else if ( opt & IT ) fkt = &mode2;
  else if ( opt & PH ) fkt = &mode4;
  else fkt = &mode3;

  for ( int i=0; i<dimX; i++ )
  {
      x += dx;
      printf( "%e\t%e\n", x, (*fkt)(field[i][0],field[i][1]) );
  }
}

// 1D Fall
void display_1D( double *field, generic_header *header )
{
  const int    dimX = header->nDimX;
  const double dx   = header->dx;
  double x    = header->xMin - header->dx;

  for ( int i=0; i<dimX; i++ )
  {
    x += dx;
    printf( "%e\t%e\n", x, field[i] );
  }
}

int main(int argc, char *argv[])
{
  int i0=-1, j0=-1, k0=-1, flags=0;

  cxxopts::Options options("gpo3", "Pipe program for gnuplot output");

  options.add_options()
  ("i,xi",  "Fix x coord. for slices at index i", cxxopts::value<int>()->default_value("-1") )
  ("j,yi",  "Fix y coord. for slices at index j", cxxopts::value<int>()->default_value("-1") )
  ("k,zi",  "Fix z coord. for slices at index k", cxxopts::value<int>()->default_value("-1") )
  ("a,re",  "output real part of a complex wave function", cxxopts::value<bool>()->default_value("false") )
  ("b,im",  "output imag part of a complex wave function", cxxopts::value<bool>()->default_value("false") )
  ("p,ph",  "output phase of a complex wave function", cxxopts::value<bool>()->default_value("false") )
  ("positional", "Positional arguments: these are the arguments that are entered without an option", cxxopts::value<std::vector<std::string>>())
  ("help","Print help")
  ;

  options.parse_positional({"positional"});
  auto result = options.parse(argc, argv);

  string filename;

  try
  {
    if (result.arguments().size() == 0)
    {
      std::cout << options.help({""}) << std::endl;
      return EXIT_FAILURE;
    }

    if( result.count("positional") > 0 )
    {
      filename = result["positional"].as<std::vector<std::string>>()[0];
    }
    else
    {
      std::cout << "error parsing options: missing file name" << std::endl;
      return EXIT_FAILURE;
    }

    i0 = result["i"].as<int>();
    j0 = result["j"].as<int>();
    k0 = result["k"].as<int>();
    if( i0 > 0 ) flags |= I_CONST;
    if( j0 > 0 ) flags |= J_CONST;
    if( k0 > 0 ) flags |= K_CONST;

    if( result["re"].as<bool>() ) flags |= RT;
    if( result["im"].as<bool>() ) flags |= IT;
    if( result["ph"].as<bool>() ) flags |= PH;
  }
  catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  generic_header header;

  ifstream in;
  in.open(filename.c_str());
  if ( !in.is_open() ) return EXIT_FAILURE;
  in.read( (char *)&header, sizeof(generic_header) );

  int no_of_threads = 4;
  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );
  omp_set_num_threads( no_of_threads );

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

  // 1D Stuff
  if ( header.nDims == 1 && header.bComplex == 1 )
  {
    total_no_bytes = sizeof(fftw_complex)*header.nDimX;
    fftw_complex *field = fftw_alloc_complex( header.nDimX );
    in.read( (char *)field, total_no_bytes );

    display_1D( field, &header, flags );

    fftw_free(field);
  }

  // 1D Stuff
  if ( header.nDims == 1 && header.bComplex == 0 )
  {
    total_no_bytes = sizeof(double)*header.nDimX;
    double *field = fftw_alloc_real( header.nDimX );
    in.read( (char *)field, total_no_bytes );

    display_1D( field, &header);

    fftw_free(field);
  }

  // 2D Stuff
  if ( header.nDims == 2 && header.bComplex == 1 )
  {
    total_no_bytes = sizeof(fftw_complex)*header.nDimX*header.nDimY;

    fftw_complex *field = fftw_alloc_complex(header.nDimX*header.nDimY);
    in.read( (char *)field, total_no_bytes );

    display_2D( field, &header, flags, i0, j0 );

    fftw_free(field);
  }

  // 2D Stuff
  if ( header.nDims == 2 && header.bComplex == 0 )
  {
    total_no_bytes = sizeof(double)*header.nDimX*header.nDimY;
    double *field = fftw_alloc_real(header.nDimX*header.nDimY);
    in.read( (char *)field, total_no_bytes );

    display_2D( field, &header, flags, i0, j0 );

    fftw_free(field);
  }

  if ( header.nDims == 3 && header.bComplex == 1  )
  {
    total_no_bytes = sizeof(fftw_complex)*header.nDimX*header.nDimY*header.nDimZ;

    fftw_complex *field = fftw_alloc_complex(header.nDimX*header.nDimY*header.nDimZ);
    in.read( (char *)field, total_no_bytes );

    display_3D( field, &header, flags, i0, j0, k0 );

    fftw_free(field);
  }

  if ( header.nDims == 3 && header.bComplex == 0  )
  {
    total_no_bytes = sizeof(double)*header.nDimX*header.nDimY*header.nDimZ;

    double *field = fftw_alloc_real(header.nDimX*header.nDimY*header.nDimZ);
    in.read( (char *)field, total_no_bytes );

    display_3D( field, &header, flags, i0, j0, k0 );

    fftw_free(field);
  }
}
