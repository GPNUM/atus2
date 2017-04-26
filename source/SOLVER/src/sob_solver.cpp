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
#include <cmath>
#include <ctime>
#include <iomanip>
#include <omp.h>
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"
#include "muParser.h"
#include "ParameterHandler.h"
#include "CSOB_Base.h"

using namespace std;

namespace SOB_Solver
{
  template<class T,int dim>
  class CSOB_Min : public CSOB_Base<T,dim,1>
  {
  public:
    CSOB_Min( ParameterHandler * );

    void Setup_Guess();
    void Setup_Potential();
  protected:
    using CSOB_Base<T,dim,1>::m_no_of_pts;
    using CSOB_Base<T,dim,1>::m_params;
    using CSOB_Base<T,dim,1>::m_fields;
    using CSOB_Base<T,dim,1>::m_Psi;
    using CSOB_Base<T,dim,1>::m_Potential;
  };

  template<class T,int dim>
  CSOB_Min<T,dim>::CSOB_Min( ParameterHandler *p ) : CSOB_Base<T,dim,1>( p )
  {
    Setup_Guess();
    Setup_Potential();
  }

  template<class T,int dim>
  void CSOB_Min<T,dim>::Setup_Potential()
  {
    #pragma omp parallel
    {
      CPoint<dim> x;
      try
      {
        std::string dimstr = std::to_string(dim);

        mu::Parser loc_mup;
        m_params->Setup_muParser( loc_mup );
        loc_mup.SetExpr(m_params->Get_simulation("POTENTIAL_" + dimstr + "D"));
        loc_mup.DefineVar("x", &x[0]);

        switch (dim)
        {
        case 2:
          loc_mup.DefineVar("y", &x[1]);
          break;
        case 3:
          loc_mup.DefineVar("y", &x[1]);
          loc_mup.DefineVar("z", &x[2]);
          break;
        }

        #pragma omp for
        for ( int l=0; l<m_no_of_pts; l++ )
        {
          x = m_fields[0]->Get_x(l);
          m_Potential[0][l] = loc_mup.Eval();
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
    }
  }

  template<class T,int dim>
  void CSOB_Min<T,dim>::Setup_Guess()
  {
    #pragma omp parallel
    {
      CPoint<dim> x;
      try
      {
        std::string dimstr = std::to_string(dim);

        mu::Parser loc_mup;
        m_params->Setup_muParser( loc_mup );
        loc_mup.SetExpr(m_params->Get_simulation("GUESS_" + dimstr + "D"));
        loc_mup.DefineVar("x", &x[0]);

        switch (dim)
        {
        case 2:
          loc_mup.DefineVar("y", &x[1]);
          break;
        case 3:
          loc_mup.DefineVar("y", &x[1]);
          loc_mup.DefineVar("z", &x[2]);
          break;
        }

        #pragma omp for
        for ( int l=0; l<m_no_of_pts; l++ )
        {
          x = m_fields[0]->Get_x(l);
          m_Psi[0][l][0] = loc_mup.Eval();
          m_Psi[0][l][1] = 0;
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
    }
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

  vector<string> filenames;
  filenames.push_back(params.Get_simulation("FILENAME"));

  try
  {
    if ( dim == 1 )
    {
      SOB_Solver::CSOB_Min<Fourier::cft_1d,1> sol( &params );
      sol.run();
      sol.Save(filenames,true); // true -> scale the wavefunction to N
      sol.Save_Zero();
    }
    else if ( dim == 2 )
    {
      SOB_Solver::CSOB_Min<Fourier::cft_2d,2> sol( &params );
      sol.run();
      sol.Save(filenames,true); // true -> scale the wavefunction to N
      sol.Save_Zero();
    }
    else if ( dim == 3 )
    {
      SOB_Solver::CSOB_Min<Fourier::cft_3d,3> sol( &params );
      sol.run();
      sol.Save(filenames,true); // true -> scale the wavefunction to N
      sol.Save_Zero();
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
