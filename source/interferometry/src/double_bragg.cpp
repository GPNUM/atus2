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
 *  |V22  Ω1    0 |
 *  |Ω1   V11   Ω2|
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
  class Bragg_double : public CRT_Base_IF<T,dim,3>
  {
  public:
    Bragg_double( ParameterHandler * );
    virtual ~Bragg_double() {};

  protected:
    static void Do_Double_Bragg_ad_Wrapper(void *, sequence_item &seq);
    void Do_Double_Bragg_ad();

    bool run_custom_sequence( const sequence_item & );

    using CRT_Base<T,dim,3>::m_map_stepfcts;
    using CRT_Base<T,dim,3>::m_gs;
    using CRT_Base_IF<T,dim,3>::DeltaL;
    using CRT_Base_IF<T,dim,3>::Amp;
    using CRT_Base_IF<T,dim,3>::laser_k;
    using CRT_Base_IF<T,dim,3>::laser_dk;
    using CRT_Base_IF<T,dim,3>::chirp;
    using CRT_Base_IF<T,dim,3>::chirp_rate;
    using CRT_Base_IF<T,dim,3>::laser_domh;
    using CRT_Base_IF<T,dim,3>::beta;
    using CRT_Base_IF<T,dim,3>::phase;
  };

  template<class T, int dim>
  Bragg_double<T,dim>::Bragg_double( ParameterHandler *p ) : CRT_Base_IF<T,dim,3>( p )
  {
    m_map_stepfcts["bragg_ad"] = &Do_Double_Bragg_ad_Wrapper;

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
  bool Bragg_double<T,dim>::run_custom_sequence( const sequence_item &item )
  {
    // return true if a custom sequence is found or else
    return false;
  }

  template<class T, int dim>
  void Bragg_double<T,dim>::Do_Double_Bragg_ad_Wrapper( void *ptr, sequence_item &seq )
  {
    Bragg_double *self = static_cast<Bragg_double *>(ptr);
    self->Do_Double_Bragg_ad();
  }

  template<class T, int dim>
  void Bragg_double<T,dim>::Do_Double_Bragg_ad()
  {
    #pragma omp parallel
    {
      fftw_complex *Psi_1 = this->m_fields[0]->Getp2In();
      fftw_complex *Psi_2 = this->m_fields[1]->Getp2In();
      fftw_complex *Psi_3 = this->m_fields[2]->Getp2In();

      const double dt = this->Get_dt();
      const double t1 = this->Get_t()+0.5*dt;
      CPoint<dim> x;

      // Pulseshapes in time
      double F = this->Amplitude_at_time();
      if( F < 0 ) F = 0;

      fftw_complex O11, O12, O21, O22, O13, O31, O33, O32, O23, gamma_1, gamma_2, gamma_3, eta;
      double re1, im1, tmp1, tmp2, tmp3, V11, V22, E1, E2, Omega_p, Omega_m, test, dw;
      int I, J;
      int i,j;

      std::cout << F << std::endl;
      #pragma omp for
      for ( int l=0; l<this->m_no_of_pts; l++ )
      {
        x = this->m_fields[0]->Get_x(l);

        tmp1 = Psi_1[l][0]*Psi_1[l][0]+Psi_1[l][1]*Psi_1[l][1];
        tmp2 = Psi_2[l][0]*Psi_2[l][0]+Psi_2[l][1]*Psi_2[l][1];
        tmp3 = Psi_3[l][0]*Psi_3[l][0]+Psi_3[l][1]*Psi_3[l][1];

        V11 = m_gs[0]*tmp1+m_gs[1]*tmp2+m_gs[2]*tmp3+beta*x;
        V22 = m_gs[3]*tmp1+m_gs[4]*tmp2+m_gs[5]*tmp3-DeltaL[1]+beta*x;

        double phase_error = phase[0];//+m_Mirror[jk];
        dw = laser_domh[0]+(chirp_rate[0])*t1;
        Omega_p = F*Amp[0]*cos(laser_k[0]*x[0]-dw*t1-phase[0]/2);
        dw = laser_domh[0]-(chirp_rate[0])*t1;
        Omega_m = F*Amp[0]*cos(-laser_k[0]*x[0]-dw*t1-phase[0]/2);

        if ((Omega_m == 0.0 ) && ( Omega_p == 0.0 ))
        {
          sincos( -dt*V11, &im1, &re1 );
          tmp1 = Psi_1[l][0];
          Psi_1[l][0] = tmp1*re1-Psi_1[l][1]*im1;
          Psi_1[l][1] = tmp1*im1+Psi_1[l][1]*re1;

          sincos( -dt*V22, &im1, &re1 );
          tmp1 = Psi_2[l][0];
          Psi_2[l][0] = tmp1*re1-Psi_2[l][1]*im1;
          Psi_2[l][1] = tmp1*im1+Psi_2[l][1]*re1;
          tmp1 = Psi_3[l][0];
          Psi_3[l][0] = tmp1*re1-Psi_3[l][1]*im1;
          Psi_3[l][1] = tmp1*im1+Psi_3[l][1]*re1;
          continue;
        }

        sincos( -0.5*laser_dk[0]*x[0], &im1, &re1 );
        eta[0] = re1;
        eta[1] = im1;

        tmp1 = sqrt((V11-V22)*(V11-V22)+4.0*(Omega_p*Omega_p+Omega_m*Omega_m));
        E1 = 0.5*(V11+V22+tmp1);
        E2 = 0.5*(V11+V22-tmp1);

        sincos( -dt*E1, &im1, &re1 );
        tmp1 = 1.0/fabs((V22-E1)*(V22-E1)+Omega_p*Omega_p+Omega_m*Omega_m);
        gamma_1[0] = tmp1*re1;
        gamma_1[1] = tmp1*im1;

        sincos( -dt*E2, &im1, &re1 );
        tmp1 = 1.0/fabs((V22-E2)*(V22-E2)+Omega_p*Omega_p+Omega_m*Omega_m);
        gamma_2[0] = tmp1*re1;
        gamma_2[1] = tmp1*im1;

        sincos( -dt*V22, &im1, &re1);
        tmp1 = 1.0/fabs(Omega_p*Omega_p+Omega_m*Omega_m);
        gamma_3[0] = tmp1*re1;
        gamma_3[1] = tmp1*im1;

        tmp1 = (V22-E2)*(V22-E2);
        tmp2 = (V22-E1)*(V22-E1);
        O11[0] = tmp1*gamma_2[0] + tmp2*gamma_1[0];
        O11[1] = tmp1*gamma_2[1] + tmp2*gamma_1[1];

        tmp1 = (V22-E2)*gamma_2[0] + (V22-E1)*gamma_1[0];
        tmp2 = (V22-E2)*gamma_2[1] + (V22-E1)*gamma_1[1];

        O12[0] = -Omega_p*tmp1;
        O12[1] = -Omega_p*tmp2;

        O13[0] = -Omega_m*tmp1;
        O13[1] = -Omega_m*tmp2;

        O21[0] = O12[0];
        O21[1] = O12[1];

        O31[0] = O13[0];
        O31[1] = O13[1];

        tmp1 = O12[0];
        O12[0] = eta[0]*tmp1-eta[1]*O12[1];
        O12[1] = eta[1]*tmp1+eta[0]*O12[1];

        tmp1 = O21[0];
        O21[0] = eta[0]*tmp1+eta[1]*O21[1];
        O21[1] = eta[0]*O21[1]-eta[1]*tmp1;

        tmp1 = O31[0];
        O31[0] = eta[0]*tmp1-eta[1]*O31[1];
        O31[1] = eta[1]*tmp1+eta[0]*O31[1];

        tmp1 = O13[0];
        O13[0] = eta[0]*tmp1+eta[1]*O13[1];
        O13[1] = eta[0]*O13[1]-eta[1]*tmp1;

        tmp1 = Omega_p*Omega_p;
        tmp2 = Omega_m*Omega_m;
        O22[0] = tmp1*(gamma_2[0]+gamma_1[0]) + tmp2*gamma_3[0];
        O22[1] = tmp1*(gamma_2[1]+gamma_1[1]) + tmp2*gamma_3[1];

        O33[0] = tmp2*(gamma_2[0]+gamma_1[0]) + tmp1*gamma_3[0];
        O33[1] = tmp2*(gamma_2[1]+gamma_1[1]) + tmp1*gamma_3[1];

        tmp1 = Omega_m*Omega_p;
        O23[0] = tmp1*(gamma_1[0]+gamma_2[0]-gamma_3[0]);
        O23[1] = tmp1*(gamma_1[1]+gamma_2[1]-gamma_3[1]);

        O32[0] = O23[0];
        O32[1] = O23[1];

        tmp1 = O23[0];
        O23[0] = (eta[0]*eta[0]+eta[1]*eta[1])*tmp1+2*eta[0]*eta[1]*O23[1];
        O23[1] = (eta[0]*eta[0]+eta[1]*eta[1])*O23[1]-2*eta[0]*eta[1]*tmp1;

        tmp1 = O32[0];
        O32[0] = (eta[0]*eta[0]-eta[1]*eta[1])*tmp1-2*eta[0]*eta[1]*O32[1];
        O32[1] = (eta[0]*eta[0]-eta[1]*eta[1])*O23[1]+2*eta[0]*eta[1]*tmp1;

        gamma_1[0] = Psi_1[l][0];
        gamma_1[1] = Psi_1[l][1];
        gamma_2[0] = Psi_2[l][0];
        gamma_2[1] = Psi_2[l][1];
        gamma_3[0] = Psi_3[l][0];
        gamma_3[1] = Psi_3[l][1];

        Psi_1[l][0] = (O11[0]*gamma_1[0]-O11[1]*gamma_1[1]) + (O12[0]*gamma_2[0]-O12[1]*gamma_2[1]) + (O13[0]*gamma_3[0]-O13[1]*gamma_3[1]);
        Psi_1[l][1] = (O11[0]*gamma_1[1]+O11[1]*gamma_1[0]) + (O12[0]*gamma_2[1]+O12[1]*gamma_2[0]) + (O13[0]*gamma_3[1]+O13[1]*gamma_3[0]);

        Psi_2[l][0] = (O21[0]*gamma_1[0]-O21[1]*gamma_1[1]) + (O22[0]*gamma_2[0]-O22[1]*gamma_2[1]) + (O23[0]*gamma_3[0]-O23[1]*gamma_3[1]);
        Psi_2[l][1] = (O21[0]*gamma_1[1]+O21[1]*gamma_1[0]) + (O22[0]*gamma_2[1]+O22[1]*gamma_2[0]) + (O23[0]*gamma_3[1]+O23[1]*gamma_3[0]);

        Psi_3[l][0] = (O31[0]*gamma_1[0]-O31[1]*gamma_1[1]) + (O32[0]*gamma_2[0]-O32[1]*gamma_2[1]) + (O33[0]*gamma_3[0]-O33[1]*gamma_3[1]);
        Psi_3[l][1] = (O31[0]*gamma_1[1]+O31[1]*gamma_1[0]) + (O32[0]*gamma_2[1]+O32[1]*gamma_2[0]) + (O33[0]*gamma_3[1]+O33[1]*gamma_3[0]);
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

  int no_of_threads = 10;
  char *envstr = getenv( "MY_NO_OF_THREADS" );
  if ( envstr != nullptr ) no_of_threads = atoi( envstr );

  fftw_init_threads();
  fftw_plan_with_nthreads( no_of_threads );
  omp_set_num_threads( no_of_threads );

  try
  {
    if ( dim == 1 )
    {
      RT_Solver::Bragg_double<Fourier::cft_1d,1> rtsol( &params );
      rtsol.run_sequence();
    }
    else if ( dim == 2 )
    {
      RT_Solver::Bragg_double<Fourier::cft_2d,2> rtsol( &params );
      rtsol.run_sequence();
    }
    else if ( dim == 3 )
    {
      RT_Solver::Bragg_double<Fourier::cft_3d,3> rtsol( &params );
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
