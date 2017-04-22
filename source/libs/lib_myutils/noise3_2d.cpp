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
#include <ostream>
#include <fstream>
#include <omp.h>
#include "noise3_2d.h"

using namespace std;

namespace Fourier
{
  /**
   * \brief CNoise2_2D Constructor
   *
   * @param header Header Information
   */
  CNoise2_2D::CNoise2_2D( const generic_header &header ) : Fourier::cft_2d( header )
  {
    m_seed    = time(nullptr);
    srand( m_seed );
  }

  /**
   * \brief Generate Noise for mirror simulations
   *
   * @param sigma_0
   * @param rho
   * @param p
   */
  void CNoise2_2D::Do_Noise_Mirror( const double sigma_0, double rho, int p )
  {
    const int max_thrno = omp_get_max_threads();

    const gsl_rng_type *T;
    gsl_rng **r = new gsl_rng*[max_thrno];

    T = gsl_rng_default;

    for ( int i=0; i<max_thrno; i++ )
    {
      r[i] = gsl_rng_alloc(T);
      gsl_rng_set( r[i], rand() );
    }

    m_bfix=true;

    int ti, tj;
    double kx, ky, im, re, PSD;
    double A = m_dx*m_dim_x*m_dy*m_dim_y;//in Einheiten von mm^2, sodass PSD=nm^2 mm^2;
    double h0 = 0.0;
    //Generate anti-symmetric Phi and normalization h0
    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        int ij = j+m_dim_y*i;
        Get_k( i, j, ti, tj, kx, ky );
        h0 += 1/(1+pow(sqrt(kx*kx+ky*ky)/rho,p));
        double zahl = gsl_ran_flat( r[omp_get_thread_num()], 0, 1 );
        m_out[ij][0] = zahl;
        m_out[ij][1] = 0;
      }
    }

    ft(-1);

    FILE *fh = fopen( "PSD.txt", "w" );
    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        Get_k( i, j, ti, tj, kx, ky );
        int tij=tj+m_dim_y*ti;
        double phi = atan2(m_out[tij][1],m_out[tij][0]);
        sincos( phi, &im, &re );
        PSD = sigma_0*sigma_0*A/(h0*4*M_PI*M_PI)*(1/(1+pow(sqrt(kx*kx+ky*ky)/rho,p)));
        fprintf(fh,"%g\t%g\n",sqrt(kx*kx+ky*ky),PSD);
        m_out[tij][0] = sqrt(A*PSD)*re;
        m_out[tij][1] = sqrt(A*PSD)*im;
      }
      fprintf( fh, "\n" );
    }
    fclose(fh);

    ft(1);

    //Calculate Autocorrelation (Fourier Transform of PSD)
    //     std::array<double,N>;
    //     cft_1d* ft1 = new cft_1d( field, header );

    //Calculate Sigma from Wavefront Map /Roughness RMS
    double sigma = 0;
    double mean = 0;

    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        int ij = j+m_dim_y*i;
        mean+=m_in[ij][0];
      }
    }

    mean = mean / (m_dim_y*m_dim_x);

    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        int ij = j+m_dim_y*i;
        sigma+=(m_in[ij][0]-mean)*(m_in[ij][0]-mean);
      }
    }
    sigma=sqrt(sigma*(m_dx*m_dy/A));
    printf("Roughness (Sigma) == %g\n",sigma);

    for ( int i=0; i<max_thrno; i++ )
      gsl_rng_free(r[i]);

    delete [] r;
  }

  /**
   * \brief Calculate simple Power Spectral Density
   *
   * @param kx x-coordinate of wavevector k
   * @param ky y-coordinate of wavevector k
   * @param exponent Exponent of Power Spectral Density
   */
  double PSD(double kx, double ky, double expX, double expY)
  {
    return fabs(1.0 / (1 + pow(kx, expX) + pow(ky, expY)) );
  }

  /** Generate 2d correlated random numbers for metric noise
   *
   */
  void CNoise2_2D::Do_Noise_Metric(double max_noise, double sigma, double exp_x, double exp_y)
  {
    const int max_thrno = omp_get_max_threads();

    const gsl_rng_type *T;
    gsl_rng **r = new gsl_rng*[max_thrno];

    T = gsl_rng_default;

    for ( int i=0; i<max_thrno; i++ )
    {
      r[i] = gsl_rng_alloc(T);
      gsl_rng_set( r[i], rand() );
    }

    m_bfix=true;

    int ti, tj, ij, tij;
    double kx, ky;
    double A = m_dx*m_dim_x*m_dy*m_dim_y;//in Einheiten von mm^2, sodass PSD=nm^2 mm^2;

    std::cout << "A: " << A << std::endl;
    #pragma omp parallel for private(ij)
    for (int i=0; i<m_dim_x; i++ )
    {
      for (int j=0; j<m_dim_y; j++ )
      {
        ij = j+m_dim_y*i;
        m_in[ij][0] = gsl_rng_uniform( r[omp_get_thread_num()] );
        m_in[ij][1] = 0.0;
      }
    }
    double h0 = 0.0;
    for (int i=0; i<m_dim_x; i++ )
    {
      for (int j=0; j<m_dim_y; j++ )
      {
        Get_k( i, j, ti, tj, kx, ky );
        h0 += PSD(kx, ky, exp_x, exp_y);
      }
    }
    printf("Setup done\n");
    ft(-1);
    printf("Fourier done\n");

    #pragma omp parallel for private(ij,ti,tj,kx,ky,tij)
    for (int i=0; i<m_dim_x; i++ )
    {
      for (int j=0; j<m_dim_y; j++ )
      {
        ij = j+m_dim_y*i;
        Get_k( i, j, ti, tj, kx, ky );
        tij=tj+m_dim_y*ti;
        if ((kx==0) || (ky==0))
        {
          m_out[tij][0] = 0.0;
          m_out[tij][1] = 0.0;
        }
        else
        {
          m_out[tij][0] *= sqrt(A*A*pow(sigma,2)*PSD(kx, ky, exp_x, exp_y)/h0);
          m_out[tij][1] *= sqrt(A*A*pow(sigma,2)*PSD(kx, ky, exp_x, exp_y)/h0);
        }
      }
    }
    printf("PSD done\n");
    ft(1);
    printf("Backfourier done\n");
    printf("max_thrno %i\n", max_thrno);

    // Get maximum norm
    double max_val = 0.0;
    #pragma omp parallel for reduction(max:max_val) private(ij)
    for (int i = 0; i < m_dim_x; i++)
    {
      for (int j = 0; j < m_dim_y; j++)
      {
        ij = j+m_dim_y*i;
        max_val = max(max_val, fabs(m_in[ij][0]));
      }
    }

    // Rescaling of noise amplitude
    #pragma omp parallel for private(ij)
    for (int i = 0; i < m_dim_x; i++)
    {
      for (int j = 0; j < m_dim_y; j++)
      {
        ij = j+m_dim_y*i;
        m_in[ij][0] *= max_noise/max_val;
      }
    }
    std::cout << "dx*NX dy*NY " << m_dx *m_dim_x << "\t" << m_dy *m_dim_y << std::endl;
    std::cout << "sigma exp_x exp_y: "<< sigma << "\t" << exp_x << "\t"<< exp_y << std::endl;
    std::cout << "maxNoise: " << max_val << std::endl;
    std::cout << "h0: " << h0 << std::endl;
    printf("all done\n");
    for ( int i=0; i<max_thrno; i++ )
      gsl_rng_free(r[i]);
    printf("free done\n");
    delete [] r;
  }

  /**
   * \brief Get real value at index (i,j)
   *
   * Get real part of array element at (i,j)
   *
   * @param i Array Index in x-dimension
   * @param j Array Index in y-dimension
   */
  double CNoise2_2D::Get_Val_re( const int i, const int j )
  {
    return m_in[j+m_dim_y*i][0];
  }

  /**
   * \brief Get real part of value at index ij
   *
   * Get real part of array element at (ij)
   *
   * @param ij Array Index
   */
  double CNoise2_2D::Get_Val( const int ij )
  {
    return m_in[ij][0];
  }

  /**
   * \brief Get imag value at index (i,j)
   *
   * Get imag part of array element at (i,j)
   *
   * @param i Array Index in x-dimension
   * @param j Array Index in y-dimension
   */
  double CNoise2_2D::Get_Val_im( const int i, const int j )
  {
    return m_in[j+m_dim_y*i][1];
  }
}
