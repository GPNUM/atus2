/* * ATUS2 - The ATUS2 package is atom interferometer Toolbox developed at ZARM
 * (CENTER OF APPLIED SPACE TECHNOLOGY AND MICROGRAVITY), Germany. This project is
 * founded by the DLR Agentur (Deutsche Luft und Raumfahrt Agentur). Grant numbers:
 * 50WM0942, 50WM1042, 50WM1342.
 * Copyright (C) 2017 Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
 *
 * This file is part of ATUS2.
 *
 * ATUS2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ATUS2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ATUS2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <array>
#include <omp.h>
#include "fftw3.h"
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"
#include "rft_1d.h"
#include "rft_2d.h"
#include "rft_3d.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

namespace Fourier
{
  template<class T, int dim>
  class CNoise : public T
  {
  public:
    CNoise( const generic_header& header, const int64_t a = time(nullptr) ) : T(header) , m_r(omp_get_max_threads())
    {
      m_seed = a;
      srand( m_seed );

      const gsl_rng_type * TT = gsl_rng_default;
      for( auto& it : m_r )
      {
        it = gsl_rng_alloc(TT);
        gsl_rng_set( it, rand() );
      }
    }

    ~CNoise()
    {
      for( auto& it : m_r )
        gsl_rng_free(it);
    }

    void reset()
    {
      for( auto& it : m_r )
        gsl_rng_free(it);

      const gsl_rng_type * TT = gsl_rng_default;
      for( auto& it : m_r )
      {
        it = gsl_rng_alloc(TT);
        gsl_rng_set( it, rand() );
      }
    }

    void white_noise()
    {
      std::cout << "dim_fs " << this->m_dim_fs << std::endl;

      #pragma omp parallel for
      for( int64_t l=0; l<this->m_dim_fs; l++ )
      {
        m_out[l][0] = gsl_ran_flat (m_r[omp_get_thread_num()],-0.5,0.5);
        m_out[l][1] = 0;
      }

      this->ft(1);

      *this *= 1/max();
    };

    void custom_noise()
    {
      CPoint<dim> k;
      #pragma omp parallel for private(k)
      for( int64_t l=0; l<this->m_dim_fs; l++ )
      {
        k = this->Get_k(l);
        m_out[l][0] = 1/(1+(k*k)) * gsl_ran_flat (m_r[omp_get_thread_num()],-0.5,0.5);
        m_out[l][1] = 0;
      }
      this->ft(1);

      *this *= 1/max();
    };

    void color_noise( const double exponent=1 )
    {
      CPoint<dim> k;

      #pragma omp parallel for private(k)
      for( int64_t l=0; l<this->m_dim_fs; l++ )
      {
        k = this->Get_k(l);

        bool bzero = false;
        for( int i=0; i<dim; i++ )
        {
          bzero |= (k[i] == 0);
        }

        if( bzero )
        {
          m_out[l][0] = 0;
          m_out[l][1] = 0;
        }
        else
        {
          double tmp=1;
          for( int i=0; i<dim; i++ )
          {
            tmp += pow(k[i],fabs(exponent));
          }
          double tmp2 = 1/tmp;
          m_out[l][0] = gsl_ran_flat (m_r[omp_get_thread_num()],-0.5,0.5)*tmp2;
          m_out[l][1] = gsl_ran_flat (m_r[omp_get_thread_num()],-0.5,0.5)*tmp2;
        }
      }

      this->ft(1);

      *this *= 1/max();
    };

    void color_noise_custom( const double exponent, const std::vector<double> mink)
    {
      std::vector<double> corr_length(dim, 1.0);
      color_noise_exp( exponent, mink, corr_length );
    };

    template <typename ContType>
    void color_noise_exp( const double exponent, const ContType& mink, const ContType& corr_length)
    // Generate separate exp(-|tau|/t_corr) correlated Noise
    {
      CPoint<dim> k;

      // Generate Gaussian Distributed White Noise in real space
      double sigma = 0.2;
      #pragma omp parallel for private(k)
      for (int64_t i = 0; i < this->m_dim; ++i) {
        m_in_real[i] = gsl_ran_gaussian_ziggurat(m_r[omp_get_thread_num()], sigma);
      }
      this->ft(-1);

      // Multiply with desired Power Spectral Density in Fourier Space
      #pragma omp parallel for private(k)
      for( int64_t l=0; l<this->m_dim_fs; l++ )
      {
        k = this->Get_k(l);

        bool bzero = false;
        for( int i=0; i<dim; i++ )
        {
          bzero |= (k[i] == 0);
        }

        // Clears values below cutoff frequency.
        for (int i = 0; i < dim; ++i) {
          if( fabs(k[i]) < mink[i] ) // fabs() needed to keep negative frequencies
          {
            bzero = true;
          }
        }

        if( bzero )
        {
          m_out[l][0] = 0;
          m_out[l][1] = 0;
        }
        else
        {
          double psd = 1.0;
          for( int i=0; i<dim; i++ )
          {
            // exp(-|tau|/t_corr) correlated Noise
            psd *= 1.0/corr_length[i]/(pow(k[i], fabs(exponent)) + pow(1.0/corr_length[i], fabs(exponent)));
          }
          double sqrt_psd = sqrt(psd);
          m_out[l][0] *= sqrt_psd;
          m_out[l][1] *= sqrt_psd;
        }
      }

      this->ft(1);

      // Normalize to 1
      auto max_value = this->max();
      std::cout << "Maximum value: " << max_value << std::endl;
      assert(1.0 > fabs(max_value));
      *this *= 1.0/max_value;
    };

    CNoise& operator*=(const double s)
    {
      if( m_in_real == nullptr )
      {
        #pragma omp parallel for
        for( int64_t l=0; l<this->m_dim; l++ )
        {
          m_in[l][0] *= s;
          m_in[l][1] *= s;
        }
      }
      else
      {
        #pragma omp parallel for
        for( int64_t l=0; l<this->m_dim; l++ )
        {
          m_in_real[l] *= s;
        }
      }
      return *this;
    };

    double max()
    {
      double retval=0;
      if( m_in_real == nullptr )
      {
        #pragma omp parallel for reduction(max:retval)
        for( int64_t l=0; l<this->m_dim; l++ )
        {
          if( fabs(m_in[l][0]) > retval ) retval = fabs(m_in[l][0]);
          if( fabs(m_in[l][1]) > retval ) retval = fabs(m_in[l][1]);
        }
      }
      else
      {
        #pragma omp parallel for reduction(max:retval)
        for( int64_t l=0; l<this->m_dim; l++ )
        {
          if( fabs(m_in_real[l]) > retval ) retval = fabs(m_in_real[l]);
        }
      }
      return retval;
    }

    void set_seed( const double a )
    {
      m_seed = a;
      srand(a);
    };

    double get_seed()
    {
      return m_seed;
    };

    double get_val_re( const int64_t i )
    {
      if( m_in != nullptr )
        return m_in[i][0];
      else
        return m_in_real[i];
    };

    double get_val_im( const int64_t i )
    {
      if( m_in != nullptr )
        return m_in[i][0];
      else
        return 0;
    };

  protected:
    int64_t m_seed;
    std::vector<gsl_rng*> m_r;

    using T::m_in;
    using T::m_in_real;
    using T::m_out;
    using T::m_header;
  };
}
