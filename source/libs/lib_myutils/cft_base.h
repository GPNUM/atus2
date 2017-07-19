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

#include <fstream>
#include <cassert>
#include <cstring>
#include "fftw3.h"
#include <cmath>
#include "CPoint.h"
#include "my_structs.h"

#pragma once

namespace Fourier
{
  enum TYPE { REAL, COMPLEX };

  template <int dim>
  class cft_base
  {
  public:
    /**
    * \brief Constructor of cft_base
    *
    * @param header Header information to construct cft_base object
    * @param b Whether inplace transformation is done
    */
    cft_base( const generic_header& header, bool b=true, bool f=false, Fourier::TYPE t=Fourier::TYPE::COMPLEX ) : m_bInplace(b), m_bfix(f), m_type(t)
    {
      if( header.nDims != dim )
      {
        std::cerr << "Critical error: header.nDims does not match template parameter dim" << std::endl;
        throw;
      }

      Setup(header);

      if ( m_type == Fourier::TYPE::COMPLEX )
      {
        if( b )
        {
          m_in_real = nullptr;
          m_in  = fftw_alloc_complex( m_dim );
          assert(m_in != nullptr);
          m_out = m_in;
          std::memset( m_in, 0, m_dim*sizeof(fftw_complex));
        }
        else
        {
          m_in_real = nullptr;
          m_in  = fftw_alloc_complex( m_dim );
          assert(m_in != nullptr);
          m_out = fftw_alloc_complex( m_dim );
          assert(m_out != nullptr);
          std::memset( m_in, 0, m_dim*sizeof(fftw_complex));
          std::memset( m_out, 0, m_dim*sizeof(fftw_complex));
        }
      }
      else
      {
          m_in_real = fftw_alloc_real( m_dim );;
          assert(m_in_real != nullptr);
          m_in  = nullptr;
          m_out = fftw_alloc_complex( m_dim_fs );
          assert(m_out != nullptr);
          std::memset( m_in_real, 0, m_dim*sizeof(double));
          std::memset( m_out, 0, m_dim_fs*sizeof(fftw_complex));
      }
    }

    /**
    * \brief Deconstructor of cft_base
    */
    virtual ~cft_base()
    {
      fftw_destroy_plan( m_forwardPlan );
      fftw_destroy_plan( m_backwardPlan );

      if ( m_type == Fourier::TYPE::COMPLEX )
      {
        if( !m_bInplace )
        {
          fftw_free( m_in );
          fftw_free( m_out );
        }
        else
        {
          fftw_free( m_in );
        }
      }
      else
      {
        fftw_free( m_in_real );
        fftw_free( m_out );
      }
    }

    void SetFix( bool bval ) { m_bfix = bval; };

    void save( const std::string& filename, bool rs=true )
    {
      std::ofstream ofs(filename);
      generic_header header = m_header;

      if ( rs && ( m_type == Fourier::TYPE::COMPLEX ) )
      {
        header.nDatatyp = sizeof(fftw_complex);
        header.bComplex = true;
        header.fs = 0;
      }
      else if ( rs && ( m_type == Fourier::TYPE::REAL ) )
      {
        header.nDatatyp = sizeof(double);
        header.bComplex = false;
        header.fs = 0;
      }
      else if ( !rs )
      {
        header.nDatatyp = sizeof(fftw_complex);
        header.bComplex = true;
        header.fs = 1;
        switch( dim )
        {
          case 1: header.nDimX = (header.nDimX/2)+1;
          break;
          case 2: header.nDimY = (header.nDimY/2)+1;
          break;
          case 3: header.nDimZ = (header.nDimZ/2)+1;
          break;
        }
      }

      ofs.write(reinterpret_cast<char*>(&header),sizeof(generic_header));
      if( rs )
      {
        if ( m_type == Fourier::TYPE::COMPLEX )
        {
          ofs.write(reinterpret_cast<char*>(m_in),sizeof(fftw_complex)*this->m_dim);
        }
        else
        {
          ofs.write(reinterpret_cast<char*>(m_in_real),sizeof(double)*this->m_dim);
        }
      }
      else
      {
        ofs.write(reinterpret_cast<char*>(m_out),sizeof(fftw_complex)*this->m_dim_fs);
      }
    };

    virtual void ft(int)=0;
    virtual CPoint<dim> Get_k(const int64_t)=0;
    virtual CPoint<dim> Get_x(const int64_t)=0;

    double * Getp2InReal() { return m_in_real; }
    fftw_complex * Getp2In() { return m_in; }
    fftw_complex * Getp2Out() { return m_out; }

    int Get_Dim_X() { return m_dim_x; };
    int Get_Dim_Y() { return m_dim_y; };
    int Get_Dim_Z() { return m_dim_z; };
    int64_t Get_red_Dim() { return m_red_dim; };
    int64_t Get_Dim_RS() { return m_dim; }; /// total number of sampling points in real space
    int64_t Get_Dim_FS() { return m_dim_fs; }; /// total number of sampling points in fourier space
  protected:

    int m_dim_x; /// Number of sampling points in x-dimension
    int m_dim_y; /// Number of sampling points in y-dimension
    int m_dim_z; /// Number of sampling points in z-dimension
    int m_shift_x; /// Index shift in x-dimension
    int m_shift_y; /// Index shift in y-dimension
    int m_shift_z; /// Index shift in z-dimension
    int64_t m_red_dim;
    int64_t m_dim; /// Product of sampling points in real space for each spatial direction
    int64_t m_dim_fs; /// Product of sampling points in fourier space for each spatial direction
    int m_isign; /// Last Transformation direction

    bool m_bInplace; /// Whether inplace transformation is performed
    bool m_bfix; /// Whether Ordering is fixed
    Fourier::TYPE m_type; /// decides if we deal with r2c or c2c

    double m_dx; /// Stepsize in x-direction
    double m_dy; /// Stepsize in y-direction
    double m_dz; /// Stepsize in z-direction
    double m_dkx; /// Stepsize in kx-direction
    double m_dky; /// Stepsize in ky-direction
    double m_dkz; /// Stepsize in kz-direction

    double * m_in_real; ///
    fftw_complex * m_in; /// Input array in real space
    fftw_complex * m_out; /// Output array in fourier space

    fftw_plan m_forwardPlan; /// Plan for forward transformation
    fftw_plan m_backwardPlan; /// Plan for backward transformation

    generic_header m_header;
  private:
    /**
    * \brief Helper routine for setting up cft_base
    *
    * @param header Header information to construct cft object
    */
    void Setup( const generic_header& header )
    {
      m_header   = header;
      m_dim_x    = header.nDimX;
      m_dim_y    = ( header.nDimY == 0 ) ? (1) : (header.nDimY);
      m_dim_z    = ( header.nDimZ == 0 ) ? (1) : (header.nDimZ);
      m_dx       = header.dx;
      m_dy       = header.dy;
      m_dz       = header.dz;
      m_shift_x  = m_dim_x/2;
      m_shift_y  = m_dim_y/2;
      m_shift_z  = m_dim_z/2;
      m_dim      = m_dim_x*m_dim_y*m_dim_z;
      m_isign    = 0;
      m_dkx      = header.dkx;
      m_dky      = header.dky;
      m_dkz      = header.dkz;
      m_dim_fs   = m_dim;
      m_red_dim  = 0;

      assert( m_dim_x > 0 );
      assert( m_dim_y >= 0 );
      assert( m_dim_z >= 0 );

      if ( m_type == Fourier::TYPE::REAL )
      {
        switch( dim )
        {
          case 1: m_red_dim = m_dim_x/2 + 1;
                  m_dim_fs = m_red_dim;
          break;
          case 2: m_red_dim = m_dim_y/2 + 1;
                  m_dim_fs =  m_dim_x*m_red_dim;
          break;
          case 3: m_red_dim = m_dim_z/2 + 1;
                  m_dim_fs = m_dim_x*m_dim_y*m_red_dim;
          break;
        }
      }
    }
  };
} // end of namespace
