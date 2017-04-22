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


#ifndef __class_cft_3d_MPI__
#define __class_cft_3d_MPI__

#include "fftw3-mpi.h"
#include "cft_base_mpi.h"

namespace MPI { namespace Fourier
{
  class cft_3d_MPI : public cft_base_MPI<3>
  {
  public:
    cft_3d_MPI( generic_header* );
    virtual ~cft_3d_MPI() {};

    void Get_k( const ptrdiff_t, const ptrdiff_t, const ptrdiff_t, ptrdiff_t&, ptrdiff_t&, ptrdiff_t&, double&, double&, double& );
    void Get_kz( const ptrdiff_t, ptrdiff_t&, double& );
    double Get_kz( const ptrdiff_t );
    double Get_z( const ptrdiff_t ); 
    
    void Laplace();
    void D_x();
    void D_y();
    void D_z();
    void D_xx();
    void D_yy();
    void D_zz();
   
    void local_to_global_rs( const ptrdiff_t, ptrdiff_t&, ptrdiff_t&, ptrdiff_t& );
    void local_to_global_rs( const ptrdiff_t, std::vector<ptrdiff_t>& );
    void local_to_global_fs( const ptrdiff_t, ptrdiff_t&, ptrdiff_t&, ptrdiff_t& );
    void local_to_global_fs( const ptrdiff_t, std::vector<ptrdiff_t>& );
    void local_to_x( const ptrdiff_t, CPoint<3>& );
    void local_to_k( const ptrdiff_t, CPoint<3>& );
   
    void ft( int isign ); // -1 (forward) oder +1 (backward)
    
    CPoint<3> Get_k(const ptrdiff_t);
    CPoint<3> Get_x(const ptrdiff_t);
  };
} }
#endif
