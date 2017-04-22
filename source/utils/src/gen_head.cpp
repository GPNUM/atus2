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
#include "my_structs.h"

int main(int argc, char *argv[])
{
  FILE *fh = nullptr;

  generic_header header = {};

  if ( argc > 1 )
  {
    fh = fopen( argv[1], "rb" );
    if ( fh == nullptr )
    {
      printf( "Could not open file %s.\n", argv[1] );
      exit(0);
    }

    fread( &header, sizeof(generic_header), 1, fh );
    fclose(fh);

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
  }
}
