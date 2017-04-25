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
#include "cft_2d_MPI.h"
#include "cft_3d_MPI.h"
#include "muParser.h"
#include "ParameterHandler.h"
#include "CSOB_Base_mpi.h"

using namespace std;

namespace MPI
{
  namespace SOB_Solver
  {
    template<class T,int dim>
    class CSOB_Min : public CSOB_Base_MPI<T,dim,1>
    {
    public:
      CSOB_Min( ParameterHandler * );
      virtual ~CSOB_Min() {};

    protected:
      void Setup_Guess();
      void Setup_Potential();

      using CSOB_Base_MPI<T,dim,1>::m_Psi;
      using CSOB_Base_MPI<T,dim,1>::m_Potential;
      using CSOB_Base_MPI<T,dim,1>::m_params;
    };

    template<class T,int dim>
    CSOB_Min<T,dim>::CSOB_Min( ParameterHandler *p ) : CSOB_Base_MPI<T,dim,1>( p )
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
          mu::Parser loc_mup;
          m_params->Setup_muParser( loc_mup );
          loc_mup.SetExpr(m_params->Get_simulation("POTENTIAL"));
          loc_mup.DefineVar("x", &x[0]);
          loc_mup.DefineVar("y", &x[1]);

          if ( dim == 3 ) loc_mup.DefineVar("z", &x[2]);

          #pragma omp for
          for ( ptrdiff_t l=0; l<this->m_no_of_pts; l++ )
          {
            x = this->m_fields[0]->Get_x(l);
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
          mu::Parser loc_mup;
          this->m_params->Setup_muParser( loc_mup );
          loc_mup.SetExpr(m_params->Get_simulation("GUESS_2D"));
          loc_mup.DefineVar("x", &x[0]);
          loc_mup.DefineVar("y", &x[1]);

          if ( dim == 3 ) loc_mup.DefineVar("z", &x[2]);

          #pragma omp for
          for ( ptrdiff_t l=0; l<this->m_no_of_pts; l++ )
          {
            x = this->m_fields[0]->Get_x(l);
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

  MPI_Init(&argc, &argv);
  fftw_mpi_init();

  vector<string> filenames;
  filenames.push_back(params.Get_simulation("FILENAME"));

  try
  {
    if ( dim == 2 )
    {
      MPI::SOB_Solver::CSOB_Min<MPI::Fourier::cft_2d_MPI,2> sol( &params );
      sol.run();
      sol.Save(filenames,true); // true -> scale the wavefunction to N
    }
    else if ( dim == 3 )
    {
      MPI::SOB_Solver::CSOB_Min<MPI::Fourier::cft_3d_MPI,3> sol( &params );
      sol.run();
      sol.Save(filenames,true); // true -> scale the wavefunction to N
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

  MPI_Finalize();
  return EXIT_SUCCESS;
}
