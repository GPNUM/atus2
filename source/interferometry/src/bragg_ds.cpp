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
 *  |V11  Ω1 |
 *  |Ω1   V22|
 *
 */

#include <cmath>
#include <iomanip>
#include <omp.h>
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"
#include "muParser.h"
#include "ParameterHandler.h"
#include "CRT_Base_IF_2.h"

using namespace std;

namespace RT_Solver_2
{
  template<class T, int dim>
  class Bragg_single : public CRT_Base_IF_2<T,dim,4>
  {
  public:
    Bragg_single( ParameterHandler * );
    virtual ~Bragg_single() {};

  protected:
    void Do_Bragg_ad(sequence_item &seq);
    static void Do_Bragg_ad_Wrapper(void *, sequence_item &seq);

    bool run_custom_sequence( const sequence_item &item );

    using CRT_Base_2<T,dim,4>::m_map_stepfcts;
    using CRT_Base_IF_2<T,dim,4>::DeltaL;
    using CRT_Base_IF_2<T,dim,4>::Amp;
    using CRT_Base_IF_2<T,dim,4>::Amp2;
    using CRT_Base_IF_2<T,dim,4>::laser_k;
    using CRT_Base_IF_2<T,dim,4>::laser_k2;
    using CRT_Base_IF_2<T,dim,4>::laser_dk;
    using CRT_Base_IF_2<T,dim,4>::laser_dk2;
    using CRT_Base_IF_2<T,dim,4>::laser_domh;
    using CRT_Base_IF_2<T,dim,4>::laser_domh2;
    using CRT_Base_IF_2<T,dim,4>::beta;
    using CRT_Base_IF_2<T,dim,4>::beta2;
    using CRT_Base_IF_2<T,dim,4>::phase;
    using CRT_Base_IF_2<T,dim,4>::phase2;
  };

  template<class T, int dim>
  Bragg_single<T,dim>::Bragg_single( ParameterHandler *p ) : CRT_Base_IF_2<T,dim,4>( p )
  {
    CPoint<dim> pt1;
    CPoint<dim> pt2;
    pt2[0] = 2*this->laser_k[0];
    CPoint<dim> pt3;
    pt3[0] = -2*this->laser_k[0];

    this->m_rabi_momentum_list.push_back(pt1);
    this->m_rabi_momentum_list.push_back(pt2);
    this->m_rabi_momentum_list.push_back(pt3);

    pt1[0] = 0;
    pt2[0] = 2*this->laser_k2[0];
    pt3[0] = -2*this->laser_k2[0];

    this->m_rabi_momentum_list2.push_back(pt1);
    this->m_rabi_momentum_list2.push_back(pt2);
    this->m_rabi_momentum_list2.push_back(pt3);

    m_map_stepfcts["bragg_ad"] = &Do_Bragg_ad_Wrapper;
  }

  template<class T, int dim>
  bool Bragg_single<T,dim>::run_custom_sequence( const sequence_item &item )
  {
    // return true if a custom sequence is found or else
    return false;
  }

  template<class T, int dim>
  void Bragg_single<T,dim>::Do_Bragg_ad_Wrapper ( void *ptr, sequence_item &seq )
  {
    Bragg_single *self = static_cast<Bragg_single *>(ptr);
    self->Do_Bragg_ad(seq);
  }

  template<class T, int dim>
  void Bragg_single<T,dim>::Do_Bragg_ad(sequence_item &seq) // ad -> analytic diagonalization
  {
    vector<fftw_complex *> Psi;
    for ( int i=0; i<4; i++ )
      Psi.push_back(this->m_fields[i]->Getp2In());

    const double dt = this->Get_dt();
    const double t1 = this->Get_t()+0.5*dt;
    int mode1=1, mode2=1;

    double time = t1-seq.time;

    //nicht richtig (nur für den ersten Puls)
    if ( time > seq.duration[0] )
      mode1 = 0;
    if ( time > seq.duration[1] )
      mode2 = 0;

    #pragma omp parallel
    {
      CPoint<dim> x;
      complex<double> M00, M01, M10, M11, eta, Psi_neu_0, Psi_neu_1, Psi_alt_0, Psi_alt_1, tmp, d0, d1;
      array<complex<double>,2> gamma;
      array<double,4> phi;
      double Ep, Em, Omega, Omega2, tmp1, V11, V22, re1, im1;

      #pragma omp for
      for ( int l=0; l<this->m_no_of_pts; l++ )
      {
        x = this->m_fields[0]->Get_x(l);

        for ( int i=0; i<2; i++ )
        {
          phi[i]=0;
          for ( int j=0; j<4; j++ )
            phi[i] += this->m_gs[j+4*i]*(Psi[j][l][0]*Psi[j][l][0] + Psi[j][l][1]*Psi[j][l][1]);
          phi[i] += beta[0]*x[0]-DeltaL[i];
        }

        for ( int i=2; i<4; i++ )
        {
          phi[i]=0;
          for ( int j=0; j<4; j++ )
            phi[i] += this->m_gs[j+4*i]*(Psi[j][l][0]*Psi[j][l][0] + Psi[j][l][1]*Psi[j][l][1]);
          phi[i] += beta2[0]*x[0]-DeltaL[i];
        }

        Omega  = mode1*Amp[0]*cos(laser_k[0]*x[0]-(laser_domh[0]+this->chirp_rate[0]*t1)*t1+phase[0]/2);
        Omega2 = mode2*Amp2[0]*cos(laser_k2[0]*x[0]-(laser_domh2[0]+this->chirp_rate2[0]*t1)*t1+phase2[0]/2);

        if ( Omega == 0 && Omega2 == 0 )
        {
          for ( int i=0; i<4; i++)
          {
            sincos( -phi[i]*dt, &im1, &re1 );
            tmp1 = Psi[i][l][0];
            Psi[i][l][0] = tmp1*re1-Psi[i][l][1]*im1;
            Psi[i][l][1] = tmp1*im1+Psi[i][l][1]*re1;
          }
          continue;
        }
        else if ( Omega == 0.0 )
        {
          for ( int i=0; i<2; i++)
          {
            sincos( -phi[i]*dt, &im1, &re1 );
            tmp1 = Psi[i][l][0];
            Psi[i][l][0] = tmp1*re1-Psi[i][l][1]*im1;
            Psi[i][l][1] = tmp1*im1+Psi[i][l][1]*re1;
          }

          V11 = phi[2];
          V22 = phi[3];

          eta = exp(complex<double>( 0, -0.5*(laser_dk2[0]*x[0]+phase2[0]) ));

          tmp1 = sqrt((V11-V22)*(V11-V22)+4.0*Omega2*Omega2);
          Ep = 0.5*(V11+V22+tmp1);
          Em = 0.5*(V11+V22-tmp1);
          d0 = V11-Ep;
          d1 = V11-Em;
          gamma[0] = exp(complex<double>( 0, -dt*Ep )) / fabs(d0*d0+Omega2*Omega2);
          gamma[1] = exp(complex<double>( 0, -dt*Em )) / fabs(d1*d1+Omega2*Omega2);

          tmp = d0*gamma[0]+d1*gamma[1];
          M00 = Omega2*Omega2*(gamma[0]+gamma[1]);
          M01 = -Omega2*eta*tmp;
          M10 = -Omega2*conj(eta)*tmp;
          M11 = d0*d0*gamma[0] + d1*d1*gamma[1];

          Psi_alt_0 = complex<double>(Psi[2][l][0], Psi[2][l][1]);
          Psi_alt_1 = complex<double>(Psi[3][l][0], Psi[3][l][1]);

          Psi_neu_0 = M00*Psi_alt_0 + M01*Psi_alt_1;
          Psi_neu_1 = M10*Psi_alt_0 + M11*Psi_alt_1;

          Psi[2][l][0] = real(Psi_neu_0);
          Psi[2][l][1] = imag(Psi_neu_0);
          Psi[3][l][0] = real(Psi_neu_1);
          Psi[3][l][1] = imag(Psi_neu_1);
          continue;
        }
        else if ( Omega2 == 0.0)
        {
          for ( int i=2; i<4; i++)
          {
            sincos( -phi[i]*dt, &im1, &re1 );
            tmp1 = Psi[i][l][0];
            Psi[i][l][0] = tmp1*re1-Psi[i][l][1]*im1;
            Psi[i][l][1] = tmp1*im1+Psi[i][l][1]*re1;
          }

          V11 = phi[0];
          V22 = phi[1];

          eta = exp(complex<double>( 0, -0.5*(laser_dk[0]*x[0]+phase[0]) ));

          tmp1 = sqrt((V11-V22)*(V11-V22)+4.0*Omega*Omega);
          Ep = 0.5*(V11+V22+tmp1);
          Em = 0.5*(V11+V22-tmp1);
          d0 = V11-Ep;
          d1 = V11-Em;
          gamma[0] = exp(complex<double>( 0, -dt*Ep )) / fabs(d0*d0+Omega*Omega);
          gamma[1] = exp(complex<double>( 0, -dt*Em )) / fabs(d1*d1+Omega*Omega);

          tmp = d0*gamma[0] + d1*gamma[1];
          M00 = Omega*Omega*(gamma[0]+gamma[1]);
          M01 = -Omega*eta*tmp;
          M10 = -Omega*conj(eta)*tmp;
          M11 = d0*d0*gamma[0] + d1*d1*gamma[1];

          Psi_alt_0 = complex<double>(Psi[0][l][0], Psi[0][l][1]);
          Psi_alt_1 = complex<double>(Psi[1][l][0], Psi[1][l][1]);

          Psi_neu_0 = M00*Psi_alt_0 + M01*Psi_alt_1;
          Psi_neu_1 = M10*Psi_alt_0 + M11*Psi_alt_1;

          Psi[0][l][0] = real(Psi_neu_0);
          Psi[0][l][1] = imag(Psi_neu_0);
          Psi[1][l][0] = real(Psi_neu_1);
          Psi[1][l][1] = imag(Psi_neu_1);
          continue;
        }

        // first species
        V11 = phi[0];
        V22 = phi[1];

        eta = exp(complex<double>( 0, -0.5*(laser_dk[0]*x[0]+phase[0]) ));

        tmp1 = sqrt((V11-V22)*(V11-V22)+4.0*Omega*Omega);
        Ep = 0.5*(V11+V22+tmp1);
        Em = 0.5*(V11+V22-tmp1);
        d0 = V11-Ep;
        d1 = V11-Em;
        gamma[0] = exp(complex<double>( 0, -dt*Ep )) / fabs(d0*d0+Omega*Omega);
        gamma[1] = exp(complex<double>( 0, -dt*Em )) / fabs(d1*d1+Omega*Omega);

        //cout << "o " <<  << sqrt((V11-Ep)*(V11-Ep)+4.0*Omega*Omega) << ", " << sqrt((V11-Ep)*(V11-Ep)+4.0*Omega*Omega) << endl;

        tmp = d0*gamma[0]+d1*gamma[1];
        M00 = Omega*Omega*(gamma[0]+gamma[1]);
        M01 = -Omega*eta*tmp;
        M10 = -Omega*conj(eta)*tmp;
        M11 = d0*d0*gamma[0] + d1*d1*gamma[1];

        Psi_alt_0 = complex<double>(Psi[0][l][0], Psi[0][l][1]);
        Psi_alt_1 = complex<double>(Psi[1][l][0], Psi[1][l][1]);

        Psi_neu_0 = M00*Psi_alt_0 + M01*Psi_alt_1;
        Psi_neu_1 = M10*Psi_alt_0 + M11*Psi_alt_1;

        Psi[0][l][0] = real(Psi_neu_0);
        Psi[0][l][1] = imag(Psi_neu_0);
        Psi[1][l][0] = real(Psi_neu_1);
        Psi[1][l][1] = imag(Psi_neu_1);

        // second species
        V11 = phi[2];
        V22 = phi[3];

        eta = exp(complex<double>( 0, -0.5*(laser_dk2[0]*x[0]+phase2[0]) ));

        tmp1 = sqrt((V11-V22)*(V11-V22)+4.0*Omega2*Omega2);
        Ep = 0.5*(V11+V22+tmp1);
        Em = 0.5*(V11+V22-tmp1);
        d0 = V11-Ep;
        d1 = V11-Em;
        gamma[0] = exp(complex<double>( 0, -dt*Ep )) / fabs(d0*d0+Omega2*Omega2);
        gamma[1] = exp(complex<double>( 0, -dt*Em )) / fabs(d1*d1+Omega2*Omega2);

        tmp = d0*gamma[0]+d1*gamma[1];
        M00 = Omega2*Omega2*(gamma[0]+gamma[1]);
        M01 = -Omega2*eta*tmp;
        M10 = -Omega2*conj(eta)*tmp;
        M11 = d0*d0*gamma[0] + d1*d1*gamma[1];

        Psi_alt_0 = complex<double>(Psi[2][l][0], Psi[2][l][1]);
        Psi_alt_1 = complex<double>(Psi[3][l][0], Psi[3][l][1]);

        Psi_neu_0 = M00*Psi_alt_0 + M01*Psi_alt_1;
        Psi_neu_1 = M10*Psi_alt_0 + M11*Psi_alt_1;

        Psi[2][l][0] = real(Psi_neu_0);
        Psi[2][l][1] = imag(Psi_neu_0);
        Psi[3][l][0] = real(Psi_neu_1);
        Psi[3][l][1] = imag(Psi_neu_1);
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

  try
  {
    if ( dim == 1 )
    {
      RT_Solver_2::Bragg_single<Fourier::cft_1d,1> rtsol( &params );
      rtsol.run_sequence();
    }
    else if ( dim == 2 )
    {
      RT_Solver_2::Bragg_single<Fourier::cft_2d,2> rtsol( &params );
      rtsol.run_sequence();
    }
    else if ( dim == 3 )
    {
      RT_Solver_2::Bragg_single<Fourier::cft_3d,3> rtsol( &params );
      rtsol.run_sequence();
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
