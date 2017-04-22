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

#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include "fftw3.h"
#include "my_structs.h"
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>

using namespace std;

template <class T> bool Read( const char *filename, generic_header *header, T *&field )
{
  FILE *fh = fopen( filename, "r" );

  if ( fh == nullptr ) return false;

  fseek( fh, sizeof(generic_header), SEEK_SET );
  size_t Nges = header->nDimX * header->nDimY * header->nDimZ;

  field = reinterpret_cast<T *>(fftw_malloc( sizeof(T)*Nges ));
  fread( (void *)field, sizeof(T), Nges, fh );
  fclose(fh);
  return true;
}

int main(int argc, char *argv[])
{
  FILE *fh = nullptr;

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

    printf( "### %s\n", argv[1] );
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
    fclose(fh);
  }

  double origin[] = {-double(header.nDimX-header.nDimX/2) *header.dx, -double(header.nDimY-header.nDimY/2) *header.dy, -double(header.nDimZ-header.nDimZ/2) *header.dz };

  if ( header.nDims == 2 )
  {
    header.nDimZ = 1;
    origin[2]=0;
    header.dz=1;
  };

  long long Nges = header.nDimX * header.nDimY * header.nDimZ;

  char filename[255];
  bzero( filename, sizeof(filename)/sizeof(char));
  int strl = strlen( argv[1] );
  memcpy( filename, argv[1], strl-4 );
  strcat( filename, ".vti" );

  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
  imageData->SetDimensions(header.nDimX,header.nDimY,header.nDimZ);

#if VTK_MAJOR_VERSION <= 5
  imageData->SetNumberOfScalarComponents(1);
  imageData->SetScalarTypeToDouble();
#else
  imageData->AllocateScalars(VTK_DOUBLE, 1);
#endif
  imageData->SetOrigin(origin);
  imageData->SetSpacing(header.dx,header.dy,header.dz);

  // 2D Stuff
  if ( header.nDims == 2 && header.bComplex == 0 && header.nDatatyp == sizeof(double) )
  {
    double *field = nullptr;
    Read( argv[1], &header, field );

    for ( long long i=0; i<header.nDimX; i++ )
    {
      for ( long long j=0; j<header.nDimY; j++)
      {
        double *pixel = static_cast<double *>(imageData->GetScalarPointer(i,j,0));
        long long ij = j+header.nDimY*i;
        pixel[0] = field[ij]*field[ij];
      }
    }
    fftw_free(field);
  }

  if ( header.nDims == 2 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    fftw_complex *field = nullptr;
    Read( argv[1], &header, field );

    for ( long long i=0; i<header.nDimX; i++ )
    {
      for ( long long j=0; j<header.nDimY; j++)
      {
        double *pixel = static_cast<double *>(imageData->GetScalarPointer(i,j,0));
        long long ij = j+header.nDimY*i;
        pixel[0] = field[ij][0]*field[ij][0]+field[ij][1]*field[ij][1];
      }
    }
    fftw_free(field);
  }

  // 3D Stuff
  if ( header.nDims == 3 && header.bComplex == 0 && header.nDatatyp == sizeof(double) )
  {
    double *field = nullptr;
    Read( argv[1], &header, field );

    for ( long long i=0; i<header.nDimX; i++ )
    {
      for ( long long j=0; j<header.nDimY; j++)
      {
        for ( long long k=0; k<header.nDimZ; k++)
        {
          double *pixel = static_cast<double *>(imageData->GetScalarPointer(i,j,k));
          long long ijk = k+header.nDimZ*(j+header.nDimY*i);
          pixel[0] = field[ijk]*field[ijk]+field[ijk]*field[ijk];
        }
      }
    }
    fftw_free(field);
  }

  if ( header.nDims == 3 && header.bComplex == 1 && header.nDatatyp == sizeof(fftw_complex) )
  {
    fftw_complex *field = nullptr;
    Read( argv[1], &header, field );

    for ( long long i=0; i<header.nDimX; i++ )
    {
      for ( long long j=0; j<header.nDimY; j++)
      {
        for ( long long k=0; k<header.nDimZ; k++)
        {
          double *pixel = static_cast<double *>(imageData->GetScalarPointer(i,j,k));
          long long ijk = k+header.nDimZ*(j+header.nDimY*i);
          pixel[0] = field[ijk][0]*field[ijk][0]+field[ijk][1]*field[ijk][1];
        }
      }
    }
    fftw_free(field);
  }

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName(filename);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInputConnection(imageData->GetProducerPort());
#else
  writer->SetInputData(imageData);
#endif
  writer->Write();

  return EXIT_SUCCESS;
}
