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

#include <ostream>
#include <fstream>
#include <string>
#include <cstring>
#include <array>

#include "strtk.hpp"
#include "CRT_Base_2.h"
#include "ParameterHandler.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_blas.h"

using namespace std;

#ifndef __class_CRT_Base_IF_2__
#define __class_CRT_Base_IF_2__

template <class T, int dim, int no_int_states>
class CRT_Base_IF_2 : public CRT_Base_2<T,dim,no_int_states>
{
public:
  CRT_Base_IF_2( ParameterHandler * );
  virtual ~CRT_Base_IF_2();

  void run_sequence();

protected:
  using CRT_Base_2<T,dim,no_int_states>::m_header;
  using CRT_Base_2<T,dim,no_int_states>::m_params;
  using CRT_Base_2<T,dim,no_int_states>::m_fields;
  using CRT_Base_2<T,dim,no_int_states>::m_custom_fct;
  using CRT_shared::m_no_of_pts;

  CPoint<dim> beta,beta2;
  std::array<double,no_int_states> DeltaL;
  std::array<double,2> Amp, Amp2, laser_k, laser_k2, laser_dk, laser_dk2, laser_domh, laser_domh2, phase, phase2, chirp_rate, chirp_rate2;

  static void Do_NL_Step_Wrapper(void *,sequence_item &seq);

  void Do_NL_Step();

  void UpdateParams();
  void Output_rabi_freq_list(string, string, const long long);
  void Output_chirps_list(string);

  virtual bool run_custom_sequence( const sequence_item & )=0;

  fftw_complex *m_workspace;

  double m_rabi_threshold;
  vector<CPoint<dim>> m_rabi_momentum_list; // this list contains the expected momenta
  vector<CPoint<dim>> m_rabi_momentum_list2;
  list<list<double>> m_rabi_freq_list; // is this list we store the rabi freq after each outer loop, if enabled
  list<list<double>> m_rabi_freq_list2;
  list<list<double>> m_chirps_list; // is this list we store the rabi freq of each last step during the phase sweep

  void compute_rabi_integrals();
};

template <class T, int dim, int no_int_states>
CRT_Base_IF_2<T,dim,no_int_states>::CRT_Base_IF_2( ParameterHandler *params ) : CRT_Base_2<T,dim,no_int_states>(params)
{
  this->m_map_stepfcts["freeprop"] = &Do_NL_Step_Wrapper;

  m_workspace=nullptr;

  UpdateParams();
}

template <class T, int dim, int no_int_states>
CRT_Base_IF_2<T,dim,no_int_states>::~CRT_Base_IF_2()
{
  if ( m_workspace != nullptr ) fftw_free(m_workspace);
}

template <class T, int dim, int no_int_states>
void CRT_Base_IF_2<T,dim,no_int_states>::UpdateParams()
{
  char s[100];
  for ( int i=0; i<dim; i++)
  {
    beta[i] = m_params->Get_VConstant("Beta",i);
    beta2[i] = m_params->Get_VConstant("Beta_2",i);
  }
  for ( int i=0; i<no_int_states; i++)
    DeltaL[i] = m_params->Get_VConstant("Delta_L",i);

  try
  {
    Amp[0] = m_params->Get_VConstant("Amp_1",0);
    Amp[1] = m_params->Get_VConstant("Amp_1",1);
    Amp2[0] = m_params->Get_VConstant("Amp_2",0);
    Amp2[1] = m_params->Get_VConstant("Amp_2",1);
    laser_k[0] = m_params->Get_Constant("laser_k");
    laser_k2[0] = m_params->Get_Constant("laser_k_2");
    laser_domh[0] = m_params->Get_Constant("laser_domh");
    laser_domh2[0] = m_params->Get_Constant("laser_domh_2");
    m_rabi_threshold = m_params->Get_Constant("rabi_threshold");
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }
  catch (std::string &str )
  {
    cout << str << endl;
  }
}

template <class T, int dim, int no_int_states>
void CRT_Base_IF_2<T,dim,no_int_states>::compute_rabi_integrals()
{
  int n = m_rabi_momentum_list.size();
  int m = m_rabi_momentum_list2.size();
  if ( n == 0 || m == 0)
  {
    std::cerr << "WARNING: Rabi momentum list is empty. Cannot compute rabi frquencies." << endl;
    return;
  }

  fftw_complex *psik = m_fields[0]->Getp2In();
  fftw_complex *psik2 = m_fields[no_int_states/2]->Getp2In();
  m_fields[0]->ft(-1);
  m_fields[no_int_states/2]->ft(-1);

  double res[n];
  double res2[m];
  memset(res,0,n*sizeof(double));
  memset(res2,0,m*sizeof(double));

  double *tmp;
  double *tmp2;
  #pragma omp parallel
  {
    const int nthreads = omp_get_num_threads();
    const int ithread = omp_get_thread_num();

    CPoint<dim> k1, d;

    #pragma omp single
    {
      tmp = new double[n*nthreads];
      for (int i=0; i<(n*nthreads); i++) tmp[i] = 0;
      tmp2 = new double[m*nthreads];
      for (int i=0; i<(m*nthreads); i++) tmp2[i] = 0;
    }

    #pragma omp for
    for ( int i=0; i<m_no_of_pts; i++ )
    {
      k1 = m_fields[0]->Get_k(i);
      int j=0;
      for ( auto k0 : m_rabi_momentum_list )
      {
        d = (k0)-k1;
        if ( sqrt(d*d) < m_rabi_threshold )
        {
          tmp[ithread*n+j] += psik[i][0]*psik[i][0] + psik[i][1]*psik[i][1];
          continue;
        }
        j++;
      }
      j=0;
      for ( auto k0 : m_rabi_momentum_list2 )
      {
        d = (k0)-k1;
        if ( sqrt(d*d) < m_rabi_threshold )
        {
          tmp2[ithread*m+j] += psik2[i][0]*psik2[i][0] + psik2[i][1]*psik2[i][1];
          continue;
        }
        j++;
      }
    }

    #pragma omp for
    for (int i=0; i<n; i++)
    {
      for (int t=0; t<nthreads; t++)
      {
        res[i] += tmp[n*t + i];
      }
    }
    #pragma omp for
    for (int i=0; i<m; i++)
    {
      for (int t=0; t<nthreads; t++)
      {
        res2[i] += tmp2[m*t + i];
      }
    }
  }
  delete [] tmp;
  delete [] tmp2;

  list<double> tmpvec;
  list<double> tmpvec2;

  for ( int i=0 ; i<n; i++ )
  {
    tmpvec.push_back(this->m_ar_k*res[i]);
  }

  for ( int i=0 ; i<m; i++ )
  {
    tmpvec2.push_back(this->m_ar_k*res2[i]);
  }

  m_rabi_freq_list.push_back(tmpvec);
  m_rabi_freq_list2.push_back(tmpvec2);

  m_fields[0]->ft(1);
  m_fields[no_int_states/2]->ft(1);
}

template <class T, int dim, int no_int_states>
void CRT_Base_IF_2<T,dim,no_int_states>::Output_rabi_freq_list( string filename, string filename2, const long long Nk )
{
  ofstream txtfile( filename );

  double t=0;

  txtfile << "# time \t";
  for ( auto i : m_rabi_momentum_list )
  {
    txtfile << i << "\t";
  }
  txtfile << endl;

  for ( auto i : m_rabi_freq_list )
  {
    txtfile << t << "\t";
    for ( auto j : i )
    {
      txtfile << j << "\t";
    }
    txtfile << endl;
    t+=this->m_header.dt*double(Nk);
  }

  ofstream txtfile2( filename2 );

  t=0;

  txtfile2 << "# time \t";
  for ( auto i : m_rabi_momentum_list2 )
  {
    txtfile2 << i << "\t";
  }
  txtfile2 << endl;

  for ( auto i : m_rabi_freq_list2 )
  {
    txtfile2 << t << "\t";
    for ( auto j : i )
    {
      txtfile2 << j << "\t";
    }
    txtfile2 << endl;
    t+=this->m_header.dt*double(Nk);
  }
}

template <class T, int dim, int no_int_states>
void CRT_Base_IF_2<T,dim,no_int_states>::Output_chirps_list( string filename )
{
  ofstream txtfile( filename );

  txtfile << "# Step \t Chirp \t ";
  for ( auto i : m_rabi_momentum_list )
  {
    txtfile << i << "\t";
  }
  txtfile << endl;

  for ( auto i : m_chirps_list )
  {
    for ( auto j : i )
    {
      txtfile << j << "\t";
    }
    txtfile << endl;
  }
}

template <class T, int dim, int no_int_states>
void CRT_Base_IF_2<T,dim,no_int_states>::Do_NL_Step_Wrapper ( void *ptr, sequence_item &seq )
{
  CRT_Base_IF_2<T,dim,no_int_states> *self = static_cast<CRT_Base_IF_2<T,dim,no_int_states>*>(ptr);
  self->Do_NL_Step();
}

template <class T, int dim, int no_int_states>
void CRT_Base_IF_2<T,dim,no_int_states>::Do_NL_Step()
{
  const double dt = -m_header.dt;
  double re1, im1, tmp1, phi[no_int_states];
  CPoint<dim> x;

  vector<fftw_complex *> Psi;
  for ( int i=0; i<no_int_states; i++ )
    Psi.push_back(m_fields[i]->Getp2In());

  for ( int l=0; l<this->m_no_of_pts; l++ )
  {
    for ( int i=0; i<no_int_states/2; i++ )
    {
      phi[i] = 0;
      for ( int j=0; j<no_int_states; j++ )
      {
        phi[i] += this->m_gs[j+no_int_states*i]*(Psi[j][l][0]*Psi[j][l][0] + Psi[j][l][1]*Psi[j][l][1]);
      }
      x = m_fields[0]->Get_x(l);
      phi[i] += beta[0]*x[0]-DeltaL[i]; //+ statt -?
      phi[i] *= dt;
    }

    for ( int i=no_int_states/2; i<no_int_states; i++ )
    {
      phi[i] = 0;
      for ( int j=0; j<no_int_states; j++ )
      {
        phi[i] += this->m_gs[j+no_int_states*i]*(Psi[j][l][0]*Psi[j][l][0] + Psi[j][l][1]*Psi[j][l][1]);
      }
      x = m_fields[0]->Get_x(l);
      phi[i] += beta2[0]*x[0]-DeltaL[i]; //+ statt -?
      phi[i] *= dt;
    }

    for ( int i=0; i<no_int_states; i++ )
    {
      sincos( phi[i], &im1, &re1 );

      tmp1 = Psi[i][l][0];
      Psi[i][l][0] = Psi[i][l][0]*re1 - Psi[i][l][1]*im1;
      Psi[i][l][1] = Psi[i][l][1]*re1 + tmp1*im1;
    }
  }
}


template <class T, int dim, int no_int_states>
void CRT_Base_IF_2<T,dim,no_int_states>::run_sequence()
{
  if ( m_fields.size() != no_int_states )
  {
    std::cerr << "Critical Error: m_fields.size() != no_int_states\n";
    exit(EXIT_FAILURE);
  }

  StepFunction step_fct=nullptr;
  StepFunction half_step_fct=nullptr;
  StepFunction full_step_fct=nullptr;

  char filename[1024];
  char filename2[1024];

  std::cout << "FYI: Found " << m_params->m_sequence.size() << " sequences." << std::endl;
  int nrm = m_rabi_momentum_list.size();
  if ( nrm == 0 )
  {
    std::cerr << "WARNING: Rabi momentum list is empty. Cannot compute rabi frquencies." << endl;
  }

  try
  {
    half_step_fct = this->m_map_stepfcts.at("half_step");
    full_step_fct = this->m_map_stepfcts.at("full_step");
  }
  catch (const std::out_of_range &oor)
  {
    std::cerr << "Critical Error: Invalid fct ptr to half_step or full_step ()" << oor.what() << ")\n";
    exit(EXIT_FAILURE);
  }

  int seq_counter=1;

  for ( auto seq : m_params->m_sequence )
  {
    if ( run_custom_sequence(seq) ) continue;

    if ( seq.name == "set_momentum" )
    {
      std::vector<std::string> vec;
      strtk::parse(seq.content,",",vec);
      CPoint<dim> P;

      assert( vec.size() > dim );
      assert( seq.comp <= dim );

      for ( int i=0; i<dim; i++ )
        P[i] = stod(vec[i]);

      this->Setup_Momentum( P, seq.comp );

      std::cout << "FYI: started new sequence " << seq.name << "\n";
      std::cout << "FYI: momentum set for component " << seq.comp << "\n";
      continue;
    }

    double max_duration = 0;
    for ( int i = 0; i < seq.duration.size(); i++)
      if (seq.duration[i] > max_duration )
        max_duration = seq.duration[i];

    int subN = int(max_duration / seq.dt);
    int Nk = seq.Nk;
    int Na = subN / seq.Nk;

    seq.time = this->Get_t();

    std::cout << "FYI: started new sequence " << seq.name << "\n";
    std::cout << "FYI: sequence no : " << seq_counter << "\n";
    std::cout << "FYI: duration    : " << max_duration << "\n";
    std::cout << "FYI: dt          : " << seq.dt << "\n";
    std::cout << "FYI: Na          : " << Na << "\n";
    std::cout << "FYI: Nk          : " << Nk << "\n";
    std::cout << "FYI: Na*Nk*dt    : " << double(Na*Nk)*seq.dt << "\n";

    if ( double(Na*Nk)*seq.dt != max_duration )
      std::cout << "FYI: double(Na*Nk)*seq.dt != max_duration\n";

    if ( this->Get_dt() != seq.dt )
      this->Set_dt(seq.dt);

    try
    {
      step_fct = this->m_map_stepfcts.at(seq.name);
    }
    catch (const std::out_of_range &oor)
    {
      std::cerr << "Critical Error: Invalid squence name " << seq.name << "\n(" << oor.what() << ")\n";
      exit(EXIT_FAILURE);
    }

    if ( seq.name == "freeprop" ) seq.no_of_chirps=1;
    double backup_t = m_header.t;

    m_chirps_list.clear();
    for ( int s=0; s<seq.no_of_chirps; ++s )
    {
      m_rabi_freq_list.clear();
      m_rabi_freq_list2.clear();

      chirp_rate[0] = 0.0;
      for ( int i=1; i<=Na; i++ )
      {
        (*half_step_fct)(this,seq);
        for ( int j=2; j<=Nk; j++ )
        {
          (*step_fct)(this,seq);
          (*full_step_fct)(this,seq);
        }
        (*step_fct)(this,seq);
        (*half_step_fct)(this,seq);

        std::cout << "t = " << to_string(m_header.t) << std::endl;

        if ( seq.output_freq == freq::each )
        {
          for ( int k=0; k<no_int_states; k++ )
          {
            sprintf( filename, "%.3f_%d.bin", this->Get_t(), k+1 );
            this->Save_Phi( filename, k );
          }
        }

        if ( seq.compute_pn_freq == freq::each )
        {
          for ( int c=0; c<no_int_states; c++ )
            std::cout << "N[" << c << "] = " << this->Get_Particle_Number(c) << std::endl;
        }

        if ( seq.rabi_output_freq == freq::each )
        {
          compute_rabi_integrals();
        }

        if ( seq.custom_freq == freq::each && m_custom_fct != nullptr )
        {
          (*m_custom_fct)(this,seq);
        }
      }

      if ( seq.output_freq == freq::last )
      {
        for ( int k=0; k<no_int_states; k++ )
        {
          sprintf( filename, "%.3f_%d.bin", this->Get_t(), k+1 );
          this->Save_Phi( filename, k );
        }
      }

      if ( seq.compute_pn_freq == freq::last )
      {
        for ( int c=0; c<no_int_states; c++ )
          std::cout << "N[" << c << "] = " << this->Get_Particle_Number(c) << std::endl;
      }

      if ( seq.custom_freq == freq::last && m_custom_fct != nullptr )
      {
        (*m_custom_fct)(this,seq);
      }

      if ( seq.rabi_output_freq == freq::last )
      {
        compute_rabi_integrals();
        cout << "Spezie 1" << endl;
        cout << "Rabi output for timestep ";
        for ( auto i : m_rabi_freq_list )
        {
          cout << this->Get_t() << "\n";
          for ( auto j : i ) cout << j << "\t";
          cout << endl;
        }

        cout << "Spezie 2" << endl;
        for ( auto i : m_rabi_freq_list2 )
        {
          cout << this->Get_t() << "\n";
          for ( auto j : i ) cout << j << "\t";
          cout << endl;
        }
      }

      if ( seq.rabi_output_freq == freq::each )
      {
        sprintf(filename, "Spezie_%d_Rabi_%d_%d.txt", 1,seq_counter, s );
        sprintf(filename2, "Spezie_%d_Rabi_%d_%d.txt", 2,seq_counter, s );

        Output_rabi_freq_list(filename,filename2,seq.Nk);
      }

      if ( seq.no_of_chirps > 1 && nrm > 0)
      {
        compute_rabi_integrals();
        list<double> tmp = m_rabi_freq_list.back();
        tmp.push_front(chirp_rate[0]);
        tmp.push_front(s);
        m_chirps_list.push_back(tmp);
      }

      // Loading files from last sequence
      if ( seq.no_of_chirps > 1 && s<seq.no_of_chirps-1 )
      {
        for ( int k=0; k<no_int_states; k++ )
        {
          sprintf( filename, "%.3f_%d.bin", backup_t, k+1 );
          ifstream file1( filename, ifstream::binary );
          file1.seekg( sizeof(generic_header), ifstream::beg );
          file1.read( (char *)m_fields[k]->Getp2In(), sizeof(fftw_complex)*m_no_of_pts );
        }
        this->m_header.t = backup_t;
      }
    }   // end of phase scan loop

    if ( seq.no_of_chirps > 1 )
    {
      sprintf(filename, "Chirp_%d.txt", seq_counter );
      Output_chirps_list(filename);
    }

    seq_counter++;
  } // end of sequence loop
}
#endif
