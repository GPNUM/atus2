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

/*  Solving:
 *  |V11  Ω1    0 |
 *  |Ω1   V22   Ω2|
 *  |0    Ω2   V33|
 */

#include <iomanip>
#include <cmath>
#include <omp.h>
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"
#include "muParser.h"
#include "ParameterHandler.h"
#include "CRT_Base_IF.h"

using namespace std;

namespace RT_Solver
{
  template<class T, int dim>
  class Raman_single : public CRT_Base_IF<T,dim,3>
  {
  public:
    Raman_single( ParameterHandler * );
    virtual ~Raman_single() {};

  protected:
    bool run_custom_sequence( const sequence_item & );
  };

  template<class T, int dim>
  Raman_single<T,dim>::Raman_single( ParameterHandler *p ) : CRT_Base_IF<T,dim,3>( p )
  {
    CPoint<dim> pt1;
    CPoint<dim> pt2;
    pt2[0] = 2*this->laser_k[0];
    CPoint<dim> pt3;
    pt3[0] = -2*this->laser_k[0];

   this->m_rabi_momentum_list.push_back(pt1);
    this->m_rabi_momentum_list.push_back(pt2);
    this->m_rabi_momentum_list.push_back(pt3);
  }

  template<class T, int dim>
  bool Raman_single<T,dim>::run_custom_sequence( const sequence_item &item )
  {
    // return true if a custom sequence is found or else
    return false;
  }
}

int main( int argc, char *argv[] )
{
  if ( argc != 2 )
  {
    printf( "No parameter xml file specified.\n" );
    return EXIT_FAILURE;
  }

  ParameterHandler params(argv[1]);
  int dim=0;

  try
  {
    std::string tmp = params.Get_simulation("DIM");
    dim = std::stod(tmp);
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }

  int no_of_threads = 4;
  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );

  fftw_init_threads();
  fftw_plan_with_nthreads( no_of_threads );
  omp_set_num_threads( no_of_threads );

  try
  {
    if ( dim == 1 )
    {
      RT_Solver::Raman_single<Fourier::cft_1d,1> rtsol( &params );
      rtsol.run_sequence();
    }
    else if ( dim == 2 )
    {
      RT_Solver::Raman_single<Fourier::cft_2d,2> rtsol( &params );
      rtsol.run_sequence();
    }
    else if ( dim == 3 )
    {
      RT_Solver::Raman_single<Fourier::cft_3d,3> rtsol( &params );
      rtsol.run_sequence();
    }
    else
    {
      cout << "You have found a new dimension!" << endl;
    }
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }
  catch (std::string &str)
  {
    cout << str << endl;
  }

  fftw_cleanup_threads();
  return EXIT_SUCCESS;
}
