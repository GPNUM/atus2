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

using namespace std;

enum
{
  O_FIX = 0x01,
  O_FAK = 0x02
};

void inflate_3d_fix( fftw_complex *field_old, generic_header *header_old, fftw_complex *field_new, generic_header *header_new, const int fak )
{
  memcpy( header_new, header_old, sizeof(generic_header) );

  header_new->xMin *= double(fak);
  header_new->xMax *= double(fak);
  header_new->yMin *= double(fak);
  header_new->yMax *= double(fak);
  header_new->zMin *= double(fak);
  header_new->zMax *= double(fak);
  header_new->dx   *= double(fak);
  header_new->dy   *= double(fak);
  header_new->dz   *= double(fak);
  header_new->dkx  /= double(fak);
  header_new->dky  /= double(fak);
  header_new->dkz  /= double(fak);
  header_new->dt = 0.001;

  const int dimX = header_old->nDimX;
  const int dimY = header_old->nDimY;
  const int dimZ = header_old->nDimZ;

  memset( field_new, 0, dimX*dimY*dimZ*sizeof(fftw_complex) );

  int i, j, k, ijk, ijk2;

  const int offsetx = dimX*(fak-1)/2/fak;
  const int offsety = dimY*(fak-1)/2/fak;
  const int offsetz = dimZ*(fak-1)/2/fak;

  //printf( "offset = %d\n", offsetx );

  #pragma omp parallel for private(j,k,ijk,ijk2)
  for ( i=0; i<dimX; i+=fak )
  {
    for ( j=0; j<dimY; j+=fak )
    {
      for ( k=0; k<dimZ; k+=fak )
      {
        ijk  = k+dimZ*(j+dimY*i);
        ijk2 = offsetz+(k/fak)+dimZ*(offsety+(j/fak)+dimY*(offsetx+(i/fak)));
        memcpy( &field_new[ijk2], &field_old[ijk], sizeof(fftw_complex) );
      }
    }
  }
}

void inflate_3d( fftw_complex *field_old, generic_header *header_old, fftw_complex *field_new, generic_header *header_new, const int fak )
{
  memcpy( header_new, header_old, sizeof(generic_header) );

  header_new->nDimX *= fak;
  header_new->nDimY *= fak;
  header_new->nDimZ *= fak;
  header_new->xMin *= double(fak);
  header_new->xMax *= double(fak);
  header_new->yMin *= double(fak);
  header_new->yMax *= double(fak);
  header_new->zMin *= double(fak);
  header_new->zMax *= double(fak);
  header_new->dkx  /= double(fak);
  header_new->dky  /= double(fak);
  header_new->dkz  /= double(fak);
  header_new->dt = 0.001;

  const int dimX = header_old->nDimX;
  const int dimY = header_old->nDimY;
  const int dimZ = header_old->nDimZ;
  const int dimY_new = header_new->nDimY;
  const int dimZ_new = header_new->nDimZ;

  memset( field_new, 0, fak*fak*fak*dimX*dimY*dimZ*sizeof(fftw_complex) );

  int i, j, k, ijk, ijk2;

  const int offsetx = (header_new->nDimX-header_old->nDimX)/2;
  const int offsety = (header_new->nDimY-header_old->nDimY)/2;
  const int offsetz = (header_new->nDimZ-header_old->nDimZ)/2;

  #pragma omp parallel for private(j,k,ijk,ijk2)
  for ( i=0; i<dimX; i++ )
  {
    for ( j=0; j<dimY; j++ )
    {
      for ( k=0; k<dimZ; k++ )
      {
        ijk  = k+dimZ*(j+dimY*i);
        ijk2 = offsetz + k + dimZ_new*((offsety+j)+dimY_new*(offsetx+i));
        memcpy( &field_new[ijk2], &field_old[ijk], sizeof(fftw_complex) );
      }
    }
  }
}

void inflate_2d_fix( fftw_complex *field_old, generic_header *header_old, fftw_complex *field_new, generic_header *header_new, const int fak )
{
  memcpy( header_new, header_old, sizeof(generic_header) );

  header_new->xMin *= double(fak);
  header_new->xMax *= double(fak);
  header_new->yMin *= double(fak);
  header_new->yMax *= double(fak);
  header_new->dx   *= double(fak);
  header_new->dy   *= double(fak);
  header_new->dkx  /= double(fak);
  header_new->dky  /= double(fak);
  header_new->dt = 0.001;

  const int dimX = header_old->nDimX;
  const int dimY = header_old->nDimY;

  memset( field_new, 0, dimX*dimY*sizeof(fftw_complex) );

  int i, j, ij, ij2;

  const int offsetx = dimX*(fak-1)/2/fak;
  const int offsety = dimY*(fak-1)/2/fak;

  #pragma omp parallel for private(j,ij,ij2)
  for ( i=0; i<dimX; i+=fak )
  {
    for ( j=0; j<dimY; j+=fak )
    {
      ij  = j+dimY*i;
      ij2 = offsety+(j/fak)+dimY*(offsetx+(i/fak));
      memcpy( &field_new[ij2], &field_old[ij], sizeof(fftw_complex) );
    }
  }
}

void inflate_2d( fftw_complex *field_old, generic_header *header_old, fftw_complex *field_new, generic_header *header_new, const int fak )
{
  memcpy( header_new, header_old, sizeof(generic_header) );

  header_new->nDimX *= fak;
  header_new->nDimY *= fak;
  header_new->nDimZ = 1;
  header_new->xMin *= double(fak);
  header_new->xMax *= double(fak);
  header_new->yMin *= double(fak);
  header_new->yMax *= double(fak);
  header_new->dkx  /= double(fak);
  header_new->dky  /= double(fak);
  header_new->dt = 0.001;

  const int dimX = header_old->nDimX;
  const int dimY = header_old->nDimY;
  const int dimY_new = header_new->nDimY;

  memset( field_new, 0, fak*fak*dimX*dimY*sizeof(fftw_complex) );

  int i, j, ij, ij2;

  const int offsetx = (header_new->nDimX-header_old->nDimX)/2;
  const int offsety = (header_new->nDimY-header_old->nDimY)/2;

  #pragma omp parallel for private(j,ij,ij2)
  for ( i=0; i<dimX; i++ )
  {
    for ( j=0; j<dimY; j++ )
    {
      ij  = j+dimY*i;
      ij2 = (offsety+j)+dimY_new*(offsetx+i);
      memcpy( &field_new[ij2], &field_old[ij], sizeof(fftw_complex) );
    }
  }
}

void inflate_1d( fftw_complex *field_old, generic_header *header_old, fftw_complex *field_new, generic_header *header_new, const int fak )
{
  memcpy( header_new, header_old, sizeof(generic_header) );

  header_new->nDimX *= fak;
  header_new->nDimY = 1;
  header_new->nDimZ = 1;
  header_new->xMin *= double(fak);
  header_new->xMax *= double(fak);
  header_new->dkx  /= double(fak);
  header_new->dt = 0.001;

  const int  offsetx = (header_new->nDimX-header_old->nDimX)/2;

  memset( field_new, 0, header_new->nDimX*sizeof(fftw_complex) );
  memcpy( &field_new[offsetx], field_old, header_old->nDimX*sizeof(fftw_complex) );
}

int main(int argc, char *argv[])
{
  int fak=-1, flags=0;
  string filename, filename2;

  cxxopts::Options options("inflate_domain", "\nInflate the domain of the wave function.\n");

  options.add_options()
  ("f,fak",  "Fix x coord. for slices at index i", cxxopts::value<int>()->default_value("-1") )
  ("g,fix",  "Fix y coord. for slices at index j", cxxopts::value<int>()->default_value("-1") )
  ("positional", "Positional arguments: these are the arguments that are entered without an option", cxxopts::value<std::vector<std::string>>())
  ("help","Print help")
  ;
  
  options.parse_positional({"positional"});
  auto result = options.parse(argc, argv);

  try
  {
    if (result.count("") == 0)
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

    if( result["fak"].as<int>() != -1 )
    {
      fak = result["fak"].as<int>();
      flags |= O_FAK;
    }
    if( result["fix"].as<int>() != -1 )
    {
      fak = result["fix"].as<int>();
      flags |= O_FIX;
    }
    if( result["fak"].as<int>() == result["fix"].as<int>() )
    {
      cout << "Hmmm" << endl;
      return EXIT_FAILURE;
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }


/*
  opt->addUsage( "" );
  opt->addUsage( "--fak integer		" );
  opt->addUsage( "--fix			inflate with fixed number of points" );

*/

  int no_of_threads = 4;
  char *envstr = getenv("MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );
  omp_set_num_threads( no_of_threads );

  generic_header header_old, header_new;
  memset( (void *)&header_new, 0, sizeof(generic_header) );
  memset( (void *)&header_old, 0, sizeof(generic_header) );

  ifstream in( filename.c_str() );
  if ( in )
  {
    in.read( (char *)&header_old, sizeof(generic_header) );

    printf( "### %s\n", filename.c_str() );
    printf( "# nDims    == %lld\n", header_old.nDims );
    printf( "# nDimX    == %lld\n", header_old.nDimX );
    printf( "# nDimY    == %lld\n", header_old.nDimY );
    printf( "# nDimZ    == %lld\n", header_old.nDimZ );
    printf( "# nDatatyp == %lld\n", header_old.nDatatyp );
    printf( "# bAtom    == %d\n", header_old.bAtom );
    printf( "# bComplex == %d\n", header_old.bComplex );
    printf( "# t        == %g\n", header_old.t );
    printf( "# dt       == %g\n", header_old.dt );
    printf( "# xMin     == %g\n", header_old.xMin );
    printf( "# xMax     == %g\n", header_old.xMax );
    printf( "# yMin     == %g\n", header_old.yMin );
    printf( "# yMax     == %g\n", header_old.yMax );
    printf( "# zMin     == %g\n", header_old.zMin );
    printf( "# zMax     == %g\n", header_old.zMax );
    printf( "# dx       == %g\n", header_old.dx );
    printf( "# dy       == %g\n", header_old.dy );
    printf( "# dz       == %g\n", header_old.dz );
    printf( "# dkx      == %g\n", header_old.dkx );
    printf( "# dky      == %g\n", header_old.dky );
    printf( "# dkz      == %g\n", header_old.dkz );
  }
  else
  {
    printf( "Cannot open file: %s\n", filename.c_str() );
    return EXIT_FAILURE;
  }

  if ( header_old.nDims == 3 && header_old.bComplex == 1 && header_old.bAtom == 1 && header_old.nDatatyp == sizeof(fftw_complex) )
  {
    // Lesen von Psi[ijk]
    in.seekg( sizeof(generic_header), in.beg );
    size_t Nges = header_old.nDimX*header_old.nDimY*header_old.nDimZ;
    size_t Nges2 = 0;
    fftw_complex *field_old = (fftw_complex *)fftw_malloc( header_old.nDatatyp*Nges );
    fftw_complex *field_new = nullptr;
    in.read( (char *)field_old, header_old.nDatatyp*Nges );
    in.close();

    if ( flags & O_FIX )
    {
      Nges2 = Nges;
      field_new = (fftw_complex *)fftw_malloc( header_old.nDatatyp*Nges2 );
      inflate_3d_fix( field_old, &header_old, field_new, &header_new, fak );
    }
    else
    {
      Nges2 = fak*fak*fak*Nges;
      field_new = (fftw_complex *)fftw_malloc( header_old.nDatatyp*Nges2 );
      inflate_3d( field_old, &header_old, field_new, &header_new, fak );
    }

    filename2 = "inf_" + filename;
    ofstream out( filename2.c_str() );
    out.write( (char *)&header_new, sizeof(generic_header) );
    out.write( (char *)field_new, Nges2*sizeof(fftw_complex) );
    out.close();

    fftw_free(field_old);
    fftw_free(field_new);
  }

  if ( header_old.nDims == 2 && header_old.bComplex == 1 && header_old.bAtom == 1 && header_old.nDatatyp == sizeof(fftw_complex) )
  {
    // Lesen von Psi[ijk]
    in.seekg( sizeof(generic_header), in.beg );
    size_t Nges = header_old.nDimX*header_old.nDimY;
    size_t Nges2 = 0;
    fftw_complex *field_old = (fftw_complex *)fftw_malloc( header_old.nDatatyp*Nges );
    fftw_complex *field_new = nullptr;
    in.read( (char *)field_old, header_old.nDatatyp*Nges );
    in.close();

    if ( flags & O_FIX )
    {
      Nges2 = Nges;
      field_new = (fftw_complex *)fftw_malloc( header_old.nDatatyp*Nges2 );
      inflate_2d_fix( field_old, &header_old, field_new, &header_new, fak );
    }
    else
    {
      Nges2 = fak*fak*Nges;
      field_new = (fftw_complex *)fftw_malloc( header_old.nDatatyp*Nges2 );
      inflate_2d( field_old, &header_old, field_new, &header_new, fak );
    }

    filename2 = "inf_" + filename;
    ofstream out( filename2.c_str() );
    out.write( (char *)&header_new, sizeof(generic_header) );
    out.write( (char *)field_new, Nges2*sizeof(fftw_complex) );
    out.close();

    fftw_free(field_old);
    fftw_free(field_new);
  }

  if ( header_old.nDims == 1 && header_old.bComplex == 1 && header_old.bAtom == 1 && header_old.nDatatyp == sizeof(fftw_complex) )
  {
    // Lesen von Psi[ijk]
    in.seekg( sizeof(generic_header), in.beg );
    size_t Nges = header_old.nDimX;
    size_t Nges2 = 0;
    fftw_complex *field_old = (fftw_complex *)fftw_malloc( header_old.nDatatyp*Nges );
    fftw_complex *field_new = nullptr;
    in.read( (char *)field_old, header_old.nDatatyp*Nges );
    in.close();

    if ( flags & O_FIX )
    {
      // Nges2 = Nges;
      // field_new = (fftw_complex*)fftw_malloc( header_old.nDatatyp*Nges2 );
      // inflate_2d_fix( field_old, &header_old, field_new, &header_new, fak );
    }
    else
    {
      Nges2 = fak*Nges;
      field_new = (fftw_complex *)fftw_malloc( header_old.nDatatyp*Nges2 );
      inflate_1d( field_old, &header_old, field_new, &header_new, fak );
    }

    filename2 = "inf_" + filename;
    ofstream out( filename2.c_str() );
    out.write( (char *)&header_new, sizeof(generic_header) );
    out.write( (char *)field_new, Nges2*sizeof(fftw_complex) );
    out.close();

    fftw_free(field_old);
    fftw_free(field_new);
  }

  printf( "# new header\n" );
  printf( "# nDims    == %lld\n", header_new.nDims );
  printf( "# nDimX    == %lld\n", header_new.nDimX );
  printf( "# nDimY    == %lld\n", header_new.nDimY );
  printf( "# nDimZ    == %lld\n", header_new.nDimZ );
  printf( "# nDatatyp == %lld\n", header_new.nDatatyp );
  printf( "# bAtom    == %d\n", header_new.bAtom );
  printf( "# bComplex == %d\n", header_new.bComplex );
  printf( "# t        == %g\n", header_new.t );
  printf( "# dt       == %g\n", header_new.dt );
  printf( "# xMin     == %g\n", header_new.xMin );
  printf( "# xMax     == %g\n", header_new.xMax );
  printf( "# yMin     == %g\n", header_new.yMin );
  printf( "# yMax     == %g\n", header_new.yMax );
  printf( "# zMin     == %g\n", header_new.zMin );
  printf( "# zMax     == %g\n", header_new.zMax );
  printf( "# dx       == %g\n", header_new.dx );
  printf( "# dy       == %g\n", header_new.dy );
  printf( "# dz       == %g\n", header_new.dz );
  printf( "# dkx      == %g\n", header_new.dkx );
  printf( "# dky      == %g\n", header_new.dky );
  printf( "# dkz      == %g\n", header_new.dkz );

  return EXIT_SUCCESS;
}
