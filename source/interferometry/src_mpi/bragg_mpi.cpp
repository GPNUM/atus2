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
#include "CRT_Base_IF_mpi.h"

using namespace std;

namespace MPI
{
  namespace RT_Solver
  {
    template<class T, int dim>
    class Bragg_single : public CRT_Base_IF_mpi<T,dim,2>
    {
    public:
      Bragg_single( ParameterHandler * );
      virtual ~Bragg_single() {};

    protected:
      void Do_Bragg_ad();
      static void Do_Bragg_ad_Wrapper(void *, sequence_item &);

      bool run_custom_sequence( const sequence_item & );

      using CRT_Base_mpi<T,dim,2>::m_map_stepfcts;
      using CRT_Base_mpi<T,dim,2>::MTime;
      using CRT_Base_IF_mpi<T,dim,2>::Amp;
      using CRT_Base_IF_mpi<T,dim,2>::laser_k;
      using CRT_Base_IF_mpi<T,dim,2>::laser_dk;
      using CRT_Base_IF_mpi<T,dim,2>::chirp;
      using CRT_Base_IF_mpi<T,dim,2>::chirp_rate;
      using CRT_Base_IF_mpi<T,dim,2>::laser_domh;
      using CRT_Base_IF_mpi<T,dim,2>::beta;
      using CRT_Base_IF_mpi<T,dim,2>::phase;
      using CRT_Base_IF_mpi<T,dim,2>::DeltaL;
    };

    template<class T, int dim>
    Bragg_single<T,dim>::Bragg_single( ParameterHandler *p ) : CRT_Base_IF_mpi<T,dim,2>( p )
    {
      this->m_fields[0]->Read_File( p->Get_simulation("FILENAME") );
      this->m_fields[1]->Read_File( p->Get_simulation("FILENAME_2") );

      m_map_stepfcts["bragg_ad"] = &Do_Bragg_ad_Wrapper;

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
    bool Bragg_single<T,dim>::run_custom_sequence( const sequence_item &item )
    {
      return false;
    }

    template<class T, int dim>
    void Bragg_single<T,dim>::Do_Bragg_ad_Wrapper ( void *ptr, sequence_item & /*seq*/ )
    {
      Bragg_single *self = static_cast<Bragg_single *>(ptr);
      self->Do_Bragg_ad();
    }

    template<class T, int dim>
    void Bragg_single<T,dim>::Do_Bragg_ad() // ad -> analytic diagonalization
    {
      MTime.enter_section(__FUNCTION__);

      fftw_complex *Psi_1 = this->m_fields[0]->Get_p2_Data();
      fftw_complex *Psi_2 = this->m_fields[1]->Get_p2_Data();

      const double dt = this->Get_dt();
      const double t1 = this->Get_t();
      fftw_complex O11, O12, O21, O22, gamma_p, gamma_m, eta;
      double re1, im1, tmp1, tmp2, V11, V22, Ep, Em, Omega, test;
      FILE *fh = nullptr;
      CPoint<dim> x;

      double dw = laser_domh[0]+chirp_rate[0]*t1;

      for ( int l=0; l<this->m_no_of_pts; l++ )
      {
        x = this->m_fields[0]->Get_x(l);

        tmp1 = Psi_1[l][0]*Psi_1[l][0]+Psi_1[l][1]*Psi_1[l][1];
        tmp2 = Psi_2[l][0]*Psi_2[l][0]+Psi_2[l][1]*Psi_2[l][1];

        V11 = this->m_gs[0]*tmp1+this->m_gs[1]*tmp2+(beta*x);
        V22 = this->m_gs[2]*tmp1+this->m_gs[3]*tmp2-DeltaL[1]+(beta*x);

        Omega = Amp[0]*cos(laser_k[0]*x[0]-dw*t1+phase[0]/2);

        if ( Omega == 0.0 )
        {
          sincos( -dt*V11, &im1, &re1 );
          tmp1 = Psi_1[l][0];
          Psi_1[l][0] = tmp1*re1-Psi_1[l][1]*im1;
          Psi_1[l][1] = tmp1*im1+Psi_1[l][1]*re1;

          sincos( -dt*V22, &im1, &re1 );
          tmp1 = Psi_2[l][0];
          Psi_2[l][0] = tmp1*re1-Psi_2[l][1]*im1;
          Psi_2[l][1] = tmp1*im1+Psi_2[l][1]*re1;
          continue;
        }

        sincos( -0.5*(laser_dk[0]*x[0]+phase[0]), &im1, &re1 );
        eta[0] = re1;
        eta[1] = im1;

        tmp1 = sqrt((V11-V22)*(V11-V22)+4.0*Omega*Omega);
        Ep = 0.5*(V11+V22+tmp1);
        Em = 0.5*(V11+V22-tmp1);

        sincos( -dt*Ep, &im1, &re1 );
        tmp1 = 1.0/fabs((V11-Ep)*(V11-Ep)+Omega*Omega);
        gamma_p[0] = tmp1*re1;
        gamma_p[1] = tmp1*im1;

        sincos( -dt*Em, &im1, &re1 );
        tmp1 = 1.0/fabs((V11-Em)*(V11-Em)+Omega*Omega);
        gamma_m[0] = tmp1*re1;
        gamma_m[1] = tmp1*im1;

        tmp1 = Omega*Omega;
        O11[0] = tmp1*(gamma_p[0]+gamma_m[0]);
        O11[1] = tmp1*(gamma_p[1]+gamma_m[1]);

        tmp1 = V11-Ep;
        tmp2 = V11-Em;

        O22[0] = tmp1*tmp1*gamma_p[0]+tmp2*tmp2*gamma_m[0];
        O22[1] = tmp1*tmp1*gamma_p[1]+tmp2*tmp2*gamma_m[1];

        O12[0] = -Omega*(tmp1*gamma_p[0]+tmp2*gamma_m[0]);
        O12[1] = -Omega*(tmp1*gamma_p[1]+tmp2*gamma_m[1]);

        O21[0] = O12[0];
        O21[1] = O12[1];

        tmp1 = O12[0];
        O12[0] = eta[0]*tmp1-eta[1]*O12[1];
        O12[1] = eta[1]*tmp1+eta[0]*O12[1];

        tmp1 = O21[0];
        O21[0] = eta[0]*tmp1+eta[1]*O21[1];
        O21[1] = eta[0]*O21[1]-eta[1]*tmp1;

        gamma_p[0] = Psi_1[l][0];
        gamma_p[1] = Psi_1[l][1];
        gamma_m[0] = Psi_2[l][0];
        gamma_m[1] = Psi_2[l][1];

        Psi_1[l][0] = (O11[0]*gamma_p[0]-O11[1]*gamma_p[1]) + (O12[0]*gamma_m[0]-O12[1]*gamma_m[1]);
        Psi_1[l][1] = (O11[0]*gamma_p[1]+O11[1]*gamma_p[0]) + (O12[0]*gamma_m[1]+O12[1]*gamma_m[0]);

        Psi_2[l][0] = (O21[0]*gamma_p[0]-O21[1]*gamma_p[1]) + (O22[0]*gamma_m[0]-O22[1]*gamma_m[1]);
        Psi_2[l][1] = (O21[0]*gamma_p[1]+O21[1]*gamma_p[0]) + (O22[0]*gamma_m[1]+O22[1]*gamma_m[0]);
      }
      MTime.exit_section(__FUNCTION__);
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

  try
  {
    if ( dim == 2 )
    {
      MPI::RT_Solver::Bragg_single<MPI::Fourier::cft_2d_MPI,2> rtsol( &params );
      rtsol.run_sequence();
    }
    else if ( dim == 3 )
    {
      MPI::RT_Solver::Bragg_single<MPI::Fourier::cft_3d_MPI,3> rtsol( &params );
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

  MPI_Finalize();
  return EXIT_SUCCESS;
}
