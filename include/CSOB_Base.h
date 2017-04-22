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

/*
 * arXiv:0906.3206v1 [quant-ph] 17 Jun 2009
 */

#include <ostream>
#include <fstream>
#include <string>
#include <cstring>
#include <array>

#include "CRT_shared.h"
#include "cft_base.h"
#include "ParameterHandler.h"

using namespace std;

#ifndef __class_CSOB_Base__
#define __class_CSOB_Base__

template <class T, int dim, int no_wf>
class CSOB_Base : public CRT_shared
{
public:
  CSOB_Base( ParameterHandler * );
  virtual ~CSOB_Base();

  double Get_Particle_Number(const int comp=0);

  void run();

  void Save( double *, std::string );
  void Save( fftw_complex *, std::string );
  void Save( std::vector<std::string>, bool=false );
  void Save_Zero();
  void Dump_2( ofstream & );
protected:
  double m_res_tot;
  double m_epsilon;
  double m_stepsize;

  void Init_Potential();
  void Compute_Laplace();
  void Compute_Sobolev_Gradient();
  void Compute_Sobolev_Psi();
  void Project_Sobolev_Gradient();
  void Compute_mu();
  void Compute_res();
  void Normalize_P_Sobolev_Gradient();
  void Renormalize_All_Psi();

  ParameterHandler *m_params;

  std::array<double,no_wf> m_mu;
  std::array<double,no_wf> m_N;
  std::array<CPoint<dim>,no_wf> m_alpha;
  std::array<double *,no_wf> m_operator_fs;
  std::array<double *,no_wf> m_Laplace_operator_fs;
  std::array<double *,no_wf> m_Potential;

  std::array<fftw_complex *,no_wf> m_Psi;
  std::array<fftw_complex *,no_wf> m_Laplace_Psi;
  std::array<fftw_complex *,no_wf> m_Psi_sob;

  void Init();

  std::array<double,dim> m_res;
  std::array<double,no_wf *no_wf> m_gs;
  std::array<T *,no_wf> m_fields;
  void Allocate();
};

/** Constructor
  * @params Pointer to a ParameterHandler object
  */
template <class T, int dim, int no_wf>
CSOB_Base<T,dim,no_wf>::CSOB_Base( ParameterHandler *params )
{
  m_params = params;

  m_header = {};
  m_header.nself = sizeof(generic_header);
  m_header.bAtom = 1;
  m_header.bComplex = 1;
  m_header.nDatatyp = sizeof(fftw_complex);
  m_header.nDims = dim;
  m_header.t    = 0;
  m_header.dt   = 0.01;

  switch ( dim )
  {
  case 1:
    m_header.nDimX = params->Get_NX();
    m_header.nDimY = 1;
    m_header.nDimZ = 1;
    m_header.xMax = params->Get_xMax();
    m_header.xMin = params->Get_xMin();
    m_header.dx   = fabs(m_header.xMax-m_header.xMin)/double(m_header.nDimX);
    m_header.dkx  = 2*M_PI/fabs(m_header.xMax-m_header.xMin);
    m_ar = m_header.dx;
    m_ar_k = m_header.dkx;
    m_no_of_pts = m_header.nDimX;
    m_shift_x = m_header.nDimX/2;
    m_no_of_pts_red = m_shift_x+1;
    break;
  case 2:
    m_header.nDimX = params->Get_NX();
    m_header.nDimY = params->Get_NY();
    m_header.nDimZ = 1;
    m_header.xMax = params->Get_xMax();
    m_header.xMin = params->Get_xMin();
    m_header.yMax = params->Get_yMax();
    m_header.yMin = params->Get_yMax();
    m_header.dx   = fabs(m_header.xMax-m_header.xMin)/double(m_header.nDimX);
    m_header.dy   = fabs(m_header.yMax-m_header.yMin)/double(m_header.nDimY);
    m_header.dkx  = 2*M_PI/fabs(m_header.xMax-m_header.xMin);
    m_header.dky  = 2*M_PI/fabs(m_header.yMax-m_header.yMin);
    m_ar = m_header.dx*m_header.dy;
    m_ar_k =  m_header.dkx*m_header.dky;
    m_no_of_pts = m_header.nDimX*m_header.nDimY;
    m_shift_x = m_header.nDimX/2;
    m_shift_y = m_header.nDimY/2;
    m_no_of_pts_red = m_header.nDimX*(m_shift_y+1);
    break;
  case 3:
    m_header.nDimX = params->Get_NX();
    m_header.nDimY = params->Get_NY();
    m_header.nDimZ = params->Get_NZ();
    m_header.xMax = params->Get_xMax();
    m_header.xMin = params->Get_xMin();
    m_header.yMax = params->Get_yMax();
    m_header.yMin = params->Get_yMax();
    m_header.zMax = params->Get_zMax();
    m_header.zMin = params->Get_zMin();
    m_header.dx   = fabs(m_header.xMax-m_header.xMin)/double(m_header.nDimX);
    m_header.dy   = fabs(m_header.yMax-m_header.yMin)/double(m_header.nDimY);
    m_header.dz   = fabs(m_header.zMax-m_header.zMin)/double(m_header.nDimZ);
    m_header.dkx  = 2*M_PI/fabs(m_header.xMax-m_header.xMin);
    m_header.dky  = 2*M_PI/fabs(m_header.yMax-m_header.yMin);
    m_header.dkz  = 2*M_PI/fabs(m_header.zMax-m_header.zMin);
    m_ar = m_header.dx*m_header.dy*m_header.dz;
    m_ar_k = m_header.dkx*m_header.dky*m_header.dkz;
    m_no_of_pts = m_header.nDimX*m_header.nDimY*m_header.nDimZ;
    m_shift_x = m_header.nDimX/2;
    m_shift_y = m_header.nDimY/2;
    m_shift_z = m_header.nDimZ/2;
    m_no_of_pts_red = m_header.nDimX*m_header.nDimY*(m_shift_z+1);
    break;
  }

  Allocate();
  Init();

  string tmpstr;

  for ( int i=0; i<no_wf; i++ )
  {
    tmpstr = "GS_" + to_string(i+1);
    for ( int j=0; j<no_wf; j++ )
    {
      m_gs[j+no_wf*i] = m_params->Get_VConstant( tmpstr, j ) * sqrt(m_params->Get_VConstant( "N", i )) * sqrt(m_params->Get_VConstant( "N", j ));
    }
  }

  for ( int i=0; i<no_wf; i++ )
  {
    tmpstr = "Alpha_" + to_string(i+1);
    m_alpha[i] = m_params->Get_VConstant( tmpstr, 0 );
  }
  m_header.dt = params->Get_dt();

  m_epsilon = 1e-10;
  m_stepsize = m_params->Get_Constant( "stepsize" );
}

template <class T, int dim, int no_wf>
CSOB_Base<T,dim,no_wf>::~CSOB_Base()
{
  for ( int i=0; i<no_wf; i++ )
  {
    delete m_fields[i];
    delete m_operator_fs[i];
    delete m_Laplace_operator_fs[i];
    delete m_Potential[i];
    delete m_Psi[i];
    delete m_Laplace_Psi[i];
    delete m_Psi_sob[i];
  }
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Allocate()
{
  for ( int i=0; i<no_wf; i++ )
  {
    m_operator_fs[i] = (double *)fftw_malloc( sizeof(double)*m_no_of_pts );
    m_Laplace_operator_fs[i] = (double *)fftw_malloc( sizeof(double)*m_no_of_pts );
    m_Potential[i] = (double *)fftw_malloc( sizeof(double)*m_no_of_pts );

    m_Psi[i] = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*m_no_of_pts );
    m_Laplace_Psi[i] = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*m_no_of_pts );
    m_Psi_sob[i] = (fftw_complex *)fftw_malloc( sizeof(fftw_complex)*m_no_of_pts );

    m_fields[i] = new T(m_header);
    m_fields[i]->SetFix(false);
  }
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Init()
{
  CPoint<dim> k;
  double phi;

  #pragma omp parallel for private(k,phi)
  for ( int j=0; j<no_wf; j++ )
  {
    double *op = m_operator_fs[j];
    double *op2 = m_Laplace_operator_fs[j];
    for ( int i=0; i<m_no_of_pts; i++ )
    {
      k = m_fields[0]->Get_k(i);
      phi = k.scale(m_alpha[j])*k;
      //phi = k*k;
      op[i] = 1.0/(1.0+phi);
      op2[i] = -phi;
    }
  }
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Compute_Laplace()
{
  fftw_complex *in = nullptr;
  fftw_complex *out = nullptr;
  fftw_complex *Psi = nullptr;
  fftw_complex *Laplace_Psi = nullptr;
  double *op = nullptr;

  // compute Laplace
  for ( int i=0; i<no_wf; i++ )
  {
    Psi = m_Psi[i];
    Laplace_Psi = m_Laplace_Psi[i];
    in = m_fields[i]->Getp2In();
    out = m_fields[i]->Getp2Out();
    op =  m_Laplace_operator_fs[i];

    memcpy( (void *)(in), (void *)(Psi), m_no_of_pts*sizeof(fftw_complex) );
    m_fields[i]->ft(-1);

    #pragma omp parallel for
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      out[l][0] *= op[l];
      out[l][1] *= op[l];
    }

    m_fields[i]->ft(1);
    memcpy( (void *)(Laplace_Psi), (void *)(in), m_no_of_pts*sizeof(fftw_complex) );
  }
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Compute_Sobolev_Psi()
{
  fftw_complex *in = nullptr;
  fftw_complex *out = nullptr;
  fftw_complex *Psi = nullptr;
  fftw_complex *Psi_sob = nullptr;
  double *op = nullptr;

  // solves (1-Laplace) Psi_sob = Psi
  for ( int i=0; i<no_wf; i++ )
  {
    Psi = m_Psi[i];
    Psi_sob = m_Psi_sob[i];
    in = m_fields[i]->Getp2In();
    out = m_fields[i]->Getp2Out();
    op = m_operator_fs[i];

    memcpy( (void *)(in), (void *)(Psi), m_no_of_pts*sizeof(fftw_complex) );
    m_fields[i]->ft(-1);

    #pragma omp parallel for
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      out[l][0] *= op[l];
      out[l][1] *= op[l];
    }

    m_fields[i]->ft(1);
    memcpy( (void *)(Psi_sob), (void *)(in), m_no_of_pts*sizeof(fftw_complex) );
  }
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Compute_Sobolev_Gradient()
{
  fftw_complex *in = nullptr;
  fftw_complex *out = nullptr;
  fftw_complex *Psi = nullptr;
  fftw_complex *Laplace_Psi = nullptr;
  double *op = nullptr;
  double *pot = nullptr;

  double NLpot;

  // solves (1-Laplace) Sobolev_Gradient = L2_Gradient
  for ( int i=0; i<no_wf; i++ )
  {
    // assemble the L2 gradient
    Psi = m_Psi[i];
    Laplace_Psi = m_Laplace_Psi[i];
    in = m_fields[i]->Getp2In();
    out = m_fields[i]->Getp2Out();
    op =  m_operator_fs[i];
    pot = m_Potential[i];

    #pragma omp parallel for
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      NLpot=0;
      for ( int j=0; j<no_wf; j++ )
        NLpot += m_gs[j+i*no_wf]*(m_Psi[j][l][0]*m_Psi[j][l][0] + m_Psi[j][l][1]*m_Psi[j][l][1]);

      in[l][0] = -Laplace_Psi[l][0] + (pot[l]+NLpot)*Psi[l][0];
      in[l][1] = -Laplace_Psi[l][1] + (pot[l]+NLpot)*Psi[l][1];
    }

    // compute the Sobolev gradient
    m_fields[i]->ft(-1);
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      out[l][0] *= op[l];
      out[l][1] *= op[l];
    }
    m_fields[i]->ft(1);
  }
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Project_Sobolev_Gradient()
{
  fftw_complex *Psi_sob = nullptr;
  fftw_complex *Psi = nullptr;
  fftw_complex *Sob_grad = nullptr;

  for ( int i=0; i<no_wf; i++ )
  {
    Sob_grad = m_fields[i]->Getp2In();
    Psi = m_Psi[i];
    Psi_sob = m_Psi_sob[i];

    double help[] = {0,0};

    #pragma omp parallel for
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      help[0] += (Psi[l][0]*Sob_grad[l][0] + Psi[l][1]*Sob_grad[l][1]);
      help[1] += (Psi[l][0]*Psi_sob[l][0] + Psi[l][1]*Psi_sob[l][1]);
    }

    double fak = help[0]/help[1];

    #pragma omp parallel for
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      Sob_grad[l][0] -= fak*Psi_sob[l][0];
      Sob_grad[l][1] -= fak*Psi_sob[l][1];
    }
  }
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Normalize_P_Sobolev_Gradient()
{
  for ( int i=0; i<no_wf; i++ )
  {
    fftw_complex *Sob_grad = m_fields[i]->Getp2In();
    double fak = 1.0/sqrt(m_res[i]);

    #pragma omp parallel for
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      Sob_grad[l][0] *= fak;
      Sob_grad[l][1] *= fak;
    }
  }
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Compute_res()
{
  fftw_complex *Sob_grad = nullptr;

  for ( int i=0; i<no_wf; i++ )
  {
    Sob_grad = m_fields[i]->Getp2In();

    double tmp=0;
    #pragma omp parallel for reduction(+:tmp)
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      tmp += (Sob_grad[l][0]*Sob_grad[l][0] + Sob_grad[l][1]*Sob_grad[l][1]);
    }
    m_res[i]=tmp;
  }

  for ( int i=0; i<no_wf; i++ )
  {
    m_res[i] *= m_ar;
  }

  m_res_tot=0;
  for ( auto i : m_res )
    m_res_tot += i*i;

  m_res_tot = sqrt(m_res_tot);
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Compute_mu()
{
  Compute_Laplace();

  fftw_complex *Psi = nullptr;
  fftw_complex *Laplace_Psi = nullptr;
  double *pot = nullptr;

  for ( int i=0; i<no_wf; i++ )
  {
    double mu=0;

    Psi = m_Psi[i];
    Laplace_Psi = m_Laplace_Psi[i];
    pot = m_Potential[i];

    #pragma omp parallel for reduction(+:mu)
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      double NLpot=0;
      for ( int j=0; j<no_wf; j++ )
        NLpot += m_gs[j+i*no_wf]*(m_Psi[j][l][0]*m_Psi[j][l][0] + m_Psi[j][l][1]*m_Psi[j][l][1]);

      mu += -(Psi[l][0]*Laplace_Psi[l][0]+Psi[l][1]*Laplace_Psi[l][1]) + (pot[l]+NLpot)*(Psi[l][0]*Psi[l][0]+Psi[l][1]*Psi[l][1]);
    }

    m_mu[i] = m_fields[i]->Get_Ar()*mu;
  }
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::run()
{
  Renormalize_All_Psi();

  int counter=0;
  do
  {
    Compute_Sobolev_Psi();

    // for( int i=0; i<no_wf; i++ )
    // {
    //   string filename = "sob_psi_" + to_string(i)+ "_" + to_string(counter) + ".bin";
    //   Save( m_Psi_sob[i], filename );
    // }

    Compute_Laplace();

    // for( int i=0; i<no_wf; i++ )
    // {
    //   string filename = "lap_psi_" + to_string(i) + "_" + to_string(counter) + ".bin";
    //   Save( m_Laplace_Psi[i], filename );
    // }

    Compute_Sobolev_Gradient();

    // for( int i=0; i<no_wf; i++ )
    // {
    //   string filename = "sob_grad_" + to_string(i) + "_" + to_string(counter) + ".bin";
    //   Save(  m_fields[i]->Getp2In(), filename );
    // }

    Project_Sobolev_Gradient();

    // for( int i=0; i<no_wf; i++ )
    // {
    //   string filename = "Psob_grad_" + to_string(i) + "_" + to_string(counter) + ".bin";
    //   Save( m_fields[i]->Getp2In(), filename );
    // }

    Compute_res();

    bool brenorm=false;
    for ( int i=0; i<no_wf; i++ )
    {
      m_N[i] = Get_Particle_Number(i);
      if ( fabs(m_N[i]-1) > 1e-5 ) brenorm = brenorm || true;
      cout << "N[" << i << "] = " << m_N[i] << endl;
    }

    printf( "brenorm = %s\n", brenorm ? "true" : "false");

    if ( brenorm )
      Renormalize_All_Psi();

    //Normalize_P_Sobolev_Gradient();

    for ( int i=0; i<no_wf; i++ )
    {
      fftw_complex *sobgrad = m_fields[i]->Getp2In();
      fftw_complex *Psi = m_Psi[i];

      #pragma omp parallel for
      for ( int l=0; l<m_no_of_pts; l++ )
      {
        Psi[l][0] -= m_stepsize*sobgrad[l][0];
        Psi[l][1] -= m_stepsize*sobgrad[l][1];
      }
    }

    // for( int i=0; i<no_wf; i++ )
    // {
    //   string filename = "psi_" + to_string(i) + "_" + to_string(counter) + ".bin";
    //   Save( m_Psi[i], filename );
    // }

    cout << "--- " << counter << endl;
    cout << "res " << m_res_tot << endl;
    counter++;
  }
  while ( m_res_tot > m_epsilon );
}

template <class T, int dim, int no_wf>
double CSOB_Base<T,dim,no_wf>::Get_Particle_Number( const int comp )
{
  if ( comp<0 || comp>no_wf ) throw std::string("Error in " + std::string(__func__) + ": comp out of bounds\n");

  fftw_complex *Psi = m_Psi[comp];
  double retval=0.0;
  #pragma omp parallel for reduction(+:retval)
  for ( int l=0; l<m_no_of_pts; l++ )
  {
    retval += (Psi[l][0]*Psi[l][0] + Psi[l][1]*Psi[l][1]);
  }
  return m_ar*retval;
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Renormalize_All_Psi()
{
  for ( int i=0; i<no_wf; i++ )
    m_N[i] = Get_Particle_Number( i );

  for ( int i=0; i<no_wf; i++ )
  {
    const double fak = 1.0/sqrt(m_N[i]);
    fftw_complex *Psi = m_Psi[i];

    #pragma omp parallel for
    for ( int l=0; l<m_no_of_pts; l++ )
    {
      Psi[l][0] *= fak;
      Psi[l][1] *= fak;
    }
  }

  for ( int i=0; i<no_wf; i++ )
    m_N[i] = 1;
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Save( double *data, std::string filename )
{
  generic_header header2 = m_header;
  header2.bComplex=false;
  header2.nDatatyp=sizeof(double);

  char *header = reinterpret_cast<char *>(&header2);
  char *Psi = reinterpret_cast<char *>(data);

  ofstream file1( filename, ofstream::binary );
  file1.write( header, sizeof(generic_header) );
  file1.write( Psi, m_no_of_pts*sizeof(fftw_complex) );
  file1.close();
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Save( fftw_complex *data, std::string filename )
{
  char *header = reinterpret_cast<char *>(&m_header);
  char *Psi = reinterpret_cast<char *>(data);

  ofstream file1( filename, ofstream::binary );
  file1.write( header, sizeof(generic_header) );
  file1.write( Psi, m_no_of_pts*sizeof(fftw_complex) );
  file1.close();
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Save( std::vector<std::string> filenames, bool rescale )
{
  if ( rescale )
  {
    for ( int i=0; i<no_wf; i++ )
    {
      fftw_complex *Psi = m_Psi[i];
      const double fak = sqrt(m_params->Get_VConstant("N",i));

      #pragma omp parallel for
      for ( int l=0; l<m_no_of_pts; l++ )
      {
        Psi[l][0] *= fak;
        Psi[l][1] *= fak;
      }
    }
  }
  for ( int i=0; i<no_wf; i++ )
    Save( m_Psi[i], filenames[i] );
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Save_Zero()
{
  memset( m_Psi_sob[0], 0, m_no_of_pts*sizeof(fftw_complex) );
  Save( m_Psi_sob[0], "zero.bin" );
}

template <class T, int dim, int no_wf>
void CSOB_Base<T,dim,no_wf>::Dump_2( ofstream &stream )
{
  stream << m_header.t << "\t";
  //stream << Get_N() << "\t";
}

template<class T, int dim, int no_wf>
ofstream &operator<<( ofstream &stream, CSOB_Base<T,dim,no_wf> &obj )
{
  obj.Dump_2(stream);
  return stream;
}

#endif
