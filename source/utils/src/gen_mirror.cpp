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
#include <fstream>
#include <cstring>
#include <array>
#include <algorithm>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "my_structs.h"
#include "noise3_2d.h"
#include "fftw3.h"
#include "zernike.h"

int main(int argc, char *argv[])
{
  const int N = 256;
  FILE *fh1=nullptr;

  //*** Setze Dateikopf-Informationen *******************************************
  //*****************************************************************************
  generic_header header = {};

  header.nself    = sizeof(generic_header);
  header.nDatatyp = sizeof(double);
  header.nDims    = 2;
  header.nDimX    = N;
  header.nDimY    = N;
  header.nDimZ    = 1;
  header.bAtom    = 1;
  header.bComplex = 0;
  header.t        = 0.0;
  header.xMin     = -50; //in \mu m
  header.xMax     = -header.xMin;
  header.yMin     = -50;
  header.yMax     = -header.yMin;
  header.dx       = fabs( header.xMax-header.xMin )/double(header.nDimX);
  header.dkx      = 2.0*M_PI/fabs(header.xMax-header.xMin);
  header.dy       = fabs( header.yMax-header.yMin )/double(header.nDimY);
  header.dky      = 2.0*M_PI/fabs(header.yMax-header.yMin);
  header.nself_and_data = header.nself + (header.nDimX + header.nDimY)*header.nDatatyp;

//   printf( "dy  == %g\n", header.dx );
//   printf( "dz  == %g\n", header.dy );
//   printf( "dky == %g\n", header.dkx );
//   printf( "dkz == %g\n", header.dky );

  double kx, ky;
  double lambda = 780.241; //nm
  int ij;

  Fourier::CNoise2_2D noise( header );

  zernike zern( header );
  zern.read_zernike();
  zern.rescale_zern(30);
  zern.print_zernike();

  double zern_errors[N*N];
  zern.calc_zernike( header, zern_errors );

  std::array<double,N *N> zern_test;

  for ( int i=0; i<header.nDimX; i++ )
  {
    for ( int j=0; j<header.nDimY; j++ )
    {
      ij  = j+N*i;
      zern_test[ij] = 2*M_PI*zern_errors[ij]/lambda;
    }
  }

  double sigma = 10.0;
  double rho = 0.1;
  int p = 3;

  if ( argc == 4 )
  {
    sigma = atof(argv[1]);
    rho = atof(argv[2]);
    p = atof(argv[3]);
  }

  //Generate Mirror
  noise.Do_Noise_Mirror( sigma, rho, p ); //sigma(nm)=10.0 rho=0.1 p=3

  std::array<double,N *N> phase;
  std::array<double,N *N> surface;
  std::array<double,N *N> full_mirror;


  for ( int i=0; i<header.nDimX; i++ )
  {
    for ( int j=0; j<header.nDimY; j++ )
    {
      ij  = j+N*i;
      surface[ij] = noise.Get_Val_re(i,j);
      phase[ij] = 2*M_PI*noise.Get_Val_re(i,j)/lambda;
      full_mirror[ij] = phase[ij] + zern_test[ij];
    }
  }

  //Mirror Analysis
  //-----------------------------------------------
  //Surface height distribution
  int n_of_bins = 100;
  double min = *std::min_element(surface.begin(),surface.end());
  double max = *std::max_element(surface.begin(),surface.end());

  std::cout << "Peak to Valley = " << max -min << std::endl;

  gsl_histogram *h = gsl_histogram_alloc(n_of_bins);
  gsl_histogram_set_ranges_uniform (h, min, max);

  for ( int l=0; l<header.nDimX*header.nDimY; l++ )
    gsl_histogram_increment(h, surface[l]);

  double height_sigma = gsl_histogram_sigma(h);
  std::cout << "Surface Height Distribution (sigma) = " << height_sigma << std::endl;

  //----------------------------------------------------
  //Autocorrelation
  generic_header header_ac = {};
  header_ac = header;
  header_ac.nDatatyp = sizeof(fftw_complex);
  header_ac.bComplex = 1;
  header_ac.nself_and_data = header.nself + (header.nDimX + header.nDimY)*header.nDatatyp;

  Fourier::cft_2d *ft = new Fourier::cft_2d( header_ac );
  ft->SetFix(true);

  fftw_complex *in = ft->Getp2In();
  fftw_complex *out = ft->Getp2Out();

  bzero(in,sizeof(fftw_complex)*N*N);
  for ( int i=0; i<N; i++ )
  {
    for ( int j=0; j<N; j++ )
    {
      ij  = j+N*i;
      in[ij][0]=surface[ij];
    }
  }

  ft->ft(-1);

  for ( int l=0; l<N*N; l++ )
  {
    in[l][0] = 2*M_PI*(in[l][0]*in[l][0]+in[l][1]*in[l][1]); //Faktor 2*pi durch doppelte Fouriertransformation
    in[l][1] = 0;
  }

  ft->ft(1);

  //----------------------------------------------------------
  //Interpolation + correlation length (1/e)
  double x[N], y[N], xi, yi;

  double lim = in[N/2+N*N/2][0]/2.71828182846;
  double eps = 10000000;
  double cor_l=1;
  for ( int i=0; i<N; i++ )
  {
    int j = N/2;
    ij = j+N*i;
    x[i]=(i-N/2)*header.dx;
    y[i]=in[ij][0];
  }

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
  gsl_spline_init (spline, x, y, N);

  fh1 = fopen("spline.txt","w");
  for (xi = x[0]; xi < x[N-1]; xi += 0.0001)
  {
    yi = gsl_spline_eval (spline, xi, acc);
    double tmp = lim - yi;
    if ( abs(tmp) < eps )
    {
      eps = abs(tmp);
      cor_l = abs(xi);
    }
    fprintf (fh1,"%g %g\n", xi, yi);
  }
  fclose(fh1);
  printf("Correlation Length = %g\n",cor_l);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  //-------------------------------------------------
  //Write
  char *header2 = reinterpret_cast<char *>(&header);
  char *header3 = reinterpret_cast<char *>(&header_ac);
  char *Mirror1 = reinterpret_cast<char *>(&phase);
  char *Mirror2 = reinterpret_cast<char *>(&surface);
  char *zernike = reinterpret_cast<char *>(&zern_test);
  char *full = reinterpret_cast<char *>(&full_mirror);
  char *korr = reinterpret_cast<char *>(in);

  ofstream file1( "Mirror.bin", ofstream::binary );
  file1.write( header2, sizeof(generic_header) );
  file1.write( Mirror1, N*N*sizeof(double) );
  file1.close();

  ofstream file2( "Surface.bin", ofstream::binary );
  file2.write( header2, sizeof(generic_header) );
  file2.write( Mirror2, N*N*sizeof(double) );
  file2.close();

  ofstream file4( "zernike.bin", ofstream::binary );
  file4.write( header2, sizeof(generic_header) );
  file4.write( zernike, N*N*sizeof(double) );
  file4.close();

  ofstream file5( "full.bin", ofstream::binary );
  file5.write( header2, sizeof(generic_header) );
  file5.write( full, N*N*sizeof(double) );
  file5.close();

  ofstream file3( "Korrelation.bin", ofstream::binary );
  file3.write( header3, sizeof(generic_header) );
  file3.write( korr, N*N*sizeof(fftw_complex) );
  file3.close();

  fh1 = fopen("surface_height.txt","w");
  gsl_histogram_fprintf(fh1,h,"%g", "%g");
  fclose(fh1);

  delete ft;
}
