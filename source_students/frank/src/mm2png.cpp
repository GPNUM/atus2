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

/******************************************************************************
Extrahiert Schnitte aus 2D Daten
*******************************************************************************
Aufrufe:
    
*******************************************************************************
..2013
*******************************************************************************
Želimir Marojević (zeli@zarm.uni-bremen.de)
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <omp.h>
#include <fstream>
#include "my_structs.h"
#include "pngwriter.h"
#include "lis.h"

using namespace std;

string removeExtension( const string& filename ) 
{
  size_t lastdot = filename.find_last_of(".");
  if (lastdot == string::npos) return filename;
return filename.substr(0, lastdot); 
}	

int main( int argc, char *argv[] )
{
//   if( argc != 4 ) 
//   {
//     return 0; 
//   }
 
  int is, ie, i, j, jj, ln, n, Mformat;
  int no_of_threads = 4;
  char* envstr = getenv( "MY_NO_OF_THREADS" );
  if( envstr != NULL ) no_of_threads = atoi( envstr );
  omp_set_num_threads( no_of_threads );

  LIS_MATRIX A;
  lis_initialize(&argc, &argv);  
  lis_matrix_create(LIS_COMM_WORLD,&A);
  lis_input_matrix( A, argv[1] );
  lis_matrix_get_range( A, &is, &ie );
  lis_matrix_get_size( A, &ln, &n );
  lis_matrix_get_type( A, &Mformat );
  
  printf( "matrix format == %d\n", Mformat );

  string filename = argv[1];
  string filename2 = removeExtension( filename );
  filename = filename2 + ".png";
  
  pngwriter pngfile( n, n, 1.0, filename.c_str() );
  for(i=0;i<n;i++)
  {
    for(j=A->ptr[i];j<A->ptr[i+1];j++)
    {
      jj = is + A->index[j]+1;
      pngfile.plot( jj, n-(is+i), 0.0, 0.0, 1.0 );
    }    
  }
  
  pngfile.close();
  
  // Cleanup
  lis_matrix_destroy( A );
  lis_finalize(); 
}
