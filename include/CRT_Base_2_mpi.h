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
#include <array>
#include "fftw3.h"
#include <cmath>
#include "my_structs.h"
#include "CRT_shared_mpi.h"
#include "ParameterHandler.h"
#include "timer.h"

using namespace std;

#ifndef __class_CRT_Base_2_mpi__
#define __class_CRT_Base_2_mpi__

template <class T, int dim, int no_int_states>
class CRT_Base_2_mpi : public CRT_shared_mpi
{
public:
  CRT_Base_2_mpi( ParameterHandler * );
  virtual ~CRT_Base_2_mpi();

  void run_sequence();

  void Dump_2( ofstream & );

  const double Get_Particle_Number( const int comp=0 ) const;

  std::array<fftw_complex *,no_int_states> m_Psi_1;
  std::array<fftw_complex *,no_int_states> m_Psi_2;

  void Setup_Momentum( CPoint<dim>, const int comp=0 );
  void Expval_Position( CPoint<dim> &, const int comp=0 );
  void Expval_Momentum( CPoint<dim> &, const int comp=0 );
  void Save_Psi( std::string, const int comp=0 );
  Timer MTime;
protected:
  void Init();
  void Allocate();
  void LoadFiles();
  void Do_FT_Step_full();
  void Do_FT_Step_half();
  static void Do_FT_Step_full_Wrapper(void *,sequence_item &);
  static void Do_FT_Step_half_Wrapper(void *,sequence_item &);

  virtual bool run_custom_sequence( const sequence_item & )=0;

  ParameterHandler *m_params;

  CPoint<dim> m_alpha_1;
  CPoint<dim> m_alpha_2;

  fftw_complex *m_full_step_1;
  fftw_complex *m_half_step_1;
  fftw_complex *m_full_step_2;
  fftw_complex *m_half_step_2;

  std::array<double,no_int_states *no_int_states> m_gs;
  std::array<T *,no_int_states> m_fields;

  std::map<std::string,StepFunction> m_map_stepfcts;
  StepFunction m_custom_fct;
};

template <class T, int dim, int no_int_states>
CRT_Base_2_mpi<T,dim,no_int_states>::CRT_Base_2_mpi( ParameterHandler *params )
{
  m_params = params;

  Read_header( m_params->Get_simulation("FILENAME"),dim);
  assert( m_header.nDims == dim );

  Allocate();
  Init();
  LoadFiles();

  //m_map_stepfcts.insert ( std::pair<string,StepFunction>("half_step",&Do_FT_Step_half_Wrapper) );
  //m_map_stepfcts.insert ( std::pair<string,StepFunction>("full_step",&Do_FT_Step_full_Wrapper) );
  m_map_stepfcts["half_step"] = &Do_FT_Step_half_Wrapper;
  m_map_stepfcts["full_step"] = &Do_FT_Step_full_Wrapper;
  m_custom_fct=nullptr;

  string tmpstr;

  for ( int i=0; i<no_int_states; i++ )
  {
    tmpstr = "GS_" + to_string(i+1);
    for ( int j=0; j<no_int_states; j++ )
    {
      m_gs[j+no_int_states*i] = m_params->Get_VConstant( tmpstr, j );
    }
  }

  for ( int i=0; i<dim; i++ )
  {
    m_alpha_1[i] = m_params->Get_VConstant( "Alpha_1", i );
    m_alpha_2[i] = m_params->Get_VConstant( "Alpha_2", i );
  }
  m_header.dt = params->Get_dt();
}

template <class T, int dim, int no_int_states>
CRT_Base_2_mpi<T,dim,no_int_states>::~CRT_Base_2_mpi()
{
  for ( int i=0; i<no_int_states; i++ )
  {
    delete m_fields[i];
  }
  fftw_free( m_full_step_1 );
  fftw_free( m_half_step_1 );
  fftw_free( m_full_step_2 );
  fftw_free( m_half_step_2 );
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Allocate()
{
  for ( int i=0; i<no_int_states; i++ )
    m_fields[i] = new T(&m_header);

  m_full_step_1 = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*m_alloc );
  m_half_step_1 = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*m_alloc );
  m_full_step_2 = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*m_alloc );
  m_half_step_2 = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*m_alloc );
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::LoadFiles()
{
  m_fields[0]->Read_File( m_params->Get_simulation("FILENAME") );

  for ( int i=1; i<no_int_states; i++ )
  {
    string str = "FILENAME_" + to_string(i+1);
    m_fields[i]->Read_File( m_params->Get_simulation(str) );
  }
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Init()
{
  const double dt = -m_header.dt;
  double phi1, phi2;

  CPoint<dim> k;

  for ( int i=0; i<m_no_of_pts; i++ )
  {
    k = m_fields[0]->Get_k(i);

    phi1 = dt*(k.scale(m_alpha_1)*k);
    phi2 = dt*(k.scale(m_alpha_2)*k);
    m_half_step_1[i][0] = cos(0.5*phi1);
    m_half_step_1[i][1] = sin(0.5*phi1);
    m_full_step_1[i][0] = cos(phi1);
    m_full_step_1[i][1] = sin(phi1);
    m_half_step_2[i][0] = cos(0.5*phi2);
    m_half_step_2[i][1] = sin(0.5*phi2);
    m_full_step_2[i][0] = cos(phi2);
    m_full_step_2[i][1] = sin(phi2);
  }
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Do_FT_Step_full()
{
  MTime.enter_section("Do_FT_Step_full");
  for ( int c=0; c<no_int_states; c++ )
    m_fields[c]->ft(-1);

  double tmp1;
  for ( int i=0; i<no_int_states; i++ )
  {
    fftw_complex *Psi = m_fields[i]->Get_p2_Data();

    for ( int l=0; l<m_no_of_pts; l++ )
    {
      tmp1 = Psi[l][0];
      Psi[l][0] = Psi[l][0]*m_full_step_1[l][0] - Psi[l][1]*m_full_step_1[l][1];
      Psi[l][1] = Psi[l][1]*m_full_step_1[l][0] + tmp1*m_full_step_1[l][1];
    }
  }

  for ( int c=0; c<no_int_states; c++ )
    m_fields[c]->ft(1);

  m_header.t += m_header.dt;
  MTime.exit_section("Do_FT_Step_full");

}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Do_FT_Step_full_Wrapper ( void *ptr, sequence_item & /*item*/ )
{
  CRT_Base_2_mpi<T,dim,no_int_states> *self = static_cast<CRT_Base_2_mpi<T,dim,no_int_states>*>(ptr);
  self->Do_FT_Step_full();
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Do_FT_Step_half_Wrapper ( void *ptr, sequence_item & /*item*/ )
{
  CRT_Base_2_mpi<T,dim,no_int_states> *self = static_cast<CRT_Base_2_mpi<T,dim,no_int_states>*>(ptr);
  self->Do_FT_Step_half();
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Do_FT_Step_half()
{
  MTime.enter_section("Do_FT_Step_half");
  for ( int c=0; c<no_int_states; c++ )
    m_fields[c]->ft(-1);

  double tmp1;
  for ( int i=0; i<no_int_states; i++ )
  {
    fftw_complex *Psi = m_fields[i]->Get_p2_Data();

    for ( int l=0; l<m_no_of_pts; l++ )
    {
      tmp1 = Psi[l][0];
      Psi[l][0] = Psi[l][0]*m_half_step_1[l][0] - Psi[l][1]*m_half_step_1[l][1];
      Psi[l][1] = Psi[l][1]*m_half_step_1[l][0] + tmp1*m_half_step_1[l][1];
    }
  }

  for ( int c=0; c<no_int_states; c++ )
    m_fields[c]->ft(1);

  m_header.t += 0.5*m_header.dt;
  MTime.exit_section("Do_FT_Step_half");

}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Setup_Momentum( CPoint<dim> px, const int comp )
{
  if ( comp<0 || comp>no_int_states ) throw std::string("Error in " + string(__func__) + ": comp out of bounds\n");

  fftw_complex *Psi = m_fields[comp]->Get_p2_Data();

  CPoint<dim> x;

  double re, im, re2, im2;

  for ( int l=0; l<m_no_of_pts; l++ )
  {
    x = m_fields[0]->Get_x(l);
    sincos(px*x,&im,&re);

    re2 = Psi[l][0];
    im2 = Psi[l][1];
    Psi[l][0] = re2*re-im2*im;
    Psi[l][1] = re2*im+im2*re;
  }
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Expval_Position( CPoint<dim> &retval, const int comp )
{
  if ( comp<0 || comp>no_int_states ) throw std::string("Error in " + string(__func__) + ": comp out of bounds\n");

  fftw_complex *Psi = m_fields[comp]->Get_p2_Data();

  double lx[dim] = {};
  double n;

  CPoint<dim> x;

  for ( int l=0; l<m_no_of_pts; l++ )
  {
    x = m_fields[0]->Get_x(l);
    n = x*(Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
    for (int i=0; i<dim; i++ )
      lx[i] += x[i]*n;
  }

  for (int i=0; i<dim; i++ )
  {
    MPI_Allreduce(&retval[i],&retval[i],1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    retval[i] = m_ar*lx[i];
  }
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Expval_Momentum( CPoint<dim> &retval, const int comp )
{
  if ( comp<0 || comp>no_int_states ) throw std::string("Error in " + string(__func__) + ": comp out of bounds\n");

  fftw_complex *Psi = m_fields[comp]->Get_p2_Data();

  CPoint<dim> k;
  double tmp1;

  fftw_complex *m_Psi_fs = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*m_no_of_pts );
  std::memcpy(m_Psi_fs, Psi, sizeof(fftw_complex)*m_no_of_pts );
  for ( int i=0; i<dim; i++ )
  {
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      k = m_fields[0]->Get_k(l);

      tmp1 = Psi[l][0];
      Psi[l][0] = -k[i]*Psi[l][1];
      Psi[l][1] = k[i]*tmp1;
    }

    m_fields[comp]->ft(1);

    tmp1=0;
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      tmp1 += m_Psi_fs[l][0]*Psi[l][1] - m_Psi_fs[l][1]*Psi[l][0];
    }
    MPI_Allreduce(&retval[i],&retval[i],1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    retval[i] = m_ar*tmp1;
    std::memcpy(Psi, m_Psi_fs, sizeof(fftw_complex)*m_no_of_pts );
  }
  m_fields[comp]->ft(-1);
  fftw_free(m_Psi_fs);
}

template <class T, int dim, int no_int_states>
const double CRT_Base_2_mpi<T,dim,no_int_states>::Get_Particle_Number( const int comp ) const
{
  if ( comp<0 || comp>no_int_states ) throw std::string("Error in " + string(__func__) + ": comp out of bounds\n");

  fftw_complex *Psi = m_fields[comp]->Get_p2_Data();

  double retval=0.0;
  for ( int l=0; l<m_no_of_pts; l++ )
  {
    retval += (Psi[l][0]*Psi[l][0] + Psi[l][1]*Psi[l][1]);
  }
  MPI_Allreduce(&retval,&retval,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return m_ar*retval;
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::Save_Psi( std::string filename, const int comp )
{
  if ( comp<0 || comp>no_int_states ) throw std::string("Error in " + string(__func__) + ": comp out of bounds\n");

  fftw_complex *Psi = m_fields[comp]->Get_p2_Data();

  double *data = reinterpret_cast<double *>(Psi);

  MPI_Status status;
  MPI_File   fh;

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  MPI_Offset offset = sizeof(generic_header) + sizeof(fftw_complex)*rank*m_no_of_pts;

  MPI_File_open( MPI_COMM_WORLD, const_cast<char *>(filename.c_str()), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );

  if ( rank == 0 )
    MPI_File_write( fh, &m_header, sizeof(generic_header), MPI_BYTE, &status );

  MPI_File_write_at( fh, offset, data, 2*m_no_of_pts, MPI_DOUBLE, MPI_STATUS_IGNORE );
  MPI_File_close( &fh );
  MPI_Barrier( MPI_COMM_WORLD );
}

template <class T, int dim, int no_int_states>
void CRT_Base_2_mpi<T,dim,no_int_states>::run_sequence()
{
  if ( m_fields.size() != no_int_states )
  {
    std::cerr << "Critical Error: m_fields.size() != no_int_states\n";
    exit(EXIT_FAILURE);
  }

  //MTime.enter_section(__func__);

  StepFunction step_fct=nullptr;
  StepFunction half_step_fct=nullptr;
  StepFunction full_step_fct=nullptr;

  char filename[1024];
  char filename2[1024];

  std::cout << "FYI: Found " << m_params->m_sequence.size() << " sequences." << std::endl;

  try
  {
    half_step_fct = this->m_map_stepfcts.at("half_step");
    full_step_fct = this->m_map_stepfcts.at("full_step");
  }
  catch (const std::out_of_range &oor)
  {
    std::cerr << "Critical Error: Invalid fct ptr to half_step or full_step ()" << oor.what() << ')\n';
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
      std::cerr << "Critical Error: Invalid squence name " << seq.name << "\n(" << oor.what() << ')\n';
      exit(EXIT_FAILURE);
    }

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
          this->Save_Psi( filename, k );
        }
      }

      if ( seq.compute_pn_freq == freq::each )
      {
        for ( int c=0; c<no_int_states; c++ )
          std::cout << "N[" << c << "] = " << this->Get_Particle_Number(c) << std::endl;
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
        m_fields[k]->Write_File(filename);
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
    seq_counter++;
  } // end of sequence loop
}

// ofstream& operator<<( ofstream& stream, CRT_Base_2_mpi& obj );
// ostream& operator<<( ostream& stream, CRT_Base_2_mpi& obj );
#endif
