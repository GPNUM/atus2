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

#include <iostream>
#include <fstream>
#include <iomanip>

#include <cstdio>
#include <cmath>
#include <complex>
#include <ctime>
#include <cstring>
#include <sys/time.h>
#include <omp.h>

#include "cft_1d.h"
#include "rft_1d.h"
#include "rft_2d.h"
#include "my_structs.h"
#include "muParser.h"
#include "ParameterHandler.h"
#include "CRT_Base_IF.h"
#include "fftw3.h"
#include "lis.h"

using namespace std;

namespace RT_Solver
{

  class Noise_Data
  {
  public:
    //! Default constructor
    Noise_Data(string filename, int no_chunks, int chunk_expansion);

    //! Copy constructor
    Noise_Data(const Noise_Data &other);

    //! Move constructor
    Noise_Data(Noise_Data &&other) noexcept;

    //! Destructor
    virtual ~Noise_Data() noexcept;

    //! Copy assignment operator
    Noise_Data& operator=(const Noise_Data &other);

    //! Move assignment operator
    Noise_Data& operator=(Noise_Data &&other) noexcept;

    const double* Get_Noise(int64_t NT);
    double dt;
  protected:
  private:
    Fourier::rft_2d *chunkft;
    Fourier::rft_2d *interpolft;
    generic_header source_header;
    generic_header chunk_header;
    generic_header interpol_header;
    int no_of_chunks;
    int expansion;
    ifstream fnoise;
    int chunk_size;
    int64_t chunk_bytes;
    int64_t current_chunk;

    void Get_Chunk(int64_t chunk);
  };

  Noise_Data::Noise_Data( string filename, int no_chunks, int chunk_expansion ){
    no_of_chunks = no_chunks;
    expansion = chunk_expansion;

    fnoise.open( filename, ifstream::binary );
    if (fnoise.fail()) {
      cout << "File not found: " << filename << endl;
      exit(EXIT_FAILURE);
    }
    fnoise.read( (char*)&source_header, sizeof(generic_header) );

    chunk_size = source_header.nDimX/no_of_chunks;
    cout << "noise_nDimX: " << source_header.nDimX << endl;
    cout << "no_of_chunks: " << no_of_chunks << endl;
    cout << "chunk_size: " << chunk_size << endl;

    chunk_header = source_header;
    chunk_header.nDimX = chunk_size;
    chunk_header.xMax = chunk_header.xMin + chunk_size*chunk_header.dx;
    chunk_header.dkx = 2.0*M_PI/fabs( chunk_header.xMax-chunk_header.xMin );
    chunk_header.nself_and_data = chunk_header.nself + (chunk_header.nDimX*chunk_header.nDimY*chunk_header.nDimZ)*chunk_header.nDatatyp;
    chunkft = new Fourier::rft_2d(chunk_header);

    chunk_bytes = chunk_header.nDimX*chunk_header.nDimY*sizeof(double);

    interpol_header = chunk_header;
    interpol_header.nDimX *= expansion;
    interpol_header.dx /= expansion;
    interpol_header.nDimY *= expansion;
    interpol_header.dy /= expansion;
    interpol_header.nself_and_data = interpol_header.nself + (interpol_header.nDimX*interpol_header.nDimY*interpol_header.nDimZ)*interpol_header.nDatatyp;
    interpolft = new Fourier::rft_2d(interpol_header);

    dt = interpol_header.dx;
    current_chunk = -1;

    cout << "source: nDimX: " << source_header.nDimX
         <<  " dx: " << source_header.dx
         << " nDimY: " << source_header.nDimY
         << " dy: " << source_header.dy << endl;
    cout << "chunk: nDimX: " << chunk_header.nDimX
         <<  " dx: " << chunk_header.dx
         << " nDimY: " << chunk_header.nDimY
         << " dy: " << chunk_header.dy << endl;
    cout << "interpol: nDimX: " << interpol_header.nDimX
         <<  " dx: " << interpol_header.dx
         << " nDimY: " << interpol_header.nDimY
         << " dy: " << interpol_header.dy << endl;
  }

  Noise_Data::~Noise_Data() {
    fnoise.close();
    delete chunkft;
    delete interpolft;
  }

  void Noise_Data::Get_Chunk(int64_t chunk) {
    double * chunk_in = chunkft->Getp2InReal();
    fftw_complex * chunk_out = chunkft->Getp2Out();
    fftw_complex * interpolft_out = interpolft->Getp2Out();

    const int64_t Nx = chunkft->Get_Dim_X();
    const int64_t Nyred = chunkft->Get_red_Dim();
    const int64_t shifti = interpolft->Get_Dim_X() - chunkft->Get_Dim_X();
    const int64_t Nynew = interpolft->Get_red_Dim();

    cout << endl << "New Chunk " << chunk;
    cout << " Old Chunk: " << current_chunk << endl;
    assert(chunk < no_of_chunks);

    // Read next chunk
    fnoise.seekg(sizeof(generic_header)+chunk*chunk_bytes);
    fnoise.read( (char*)chunk_in, chunk_bytes);
    // chunkft->save( "chunk_" + to_string(chunk) + ".bin" );

    // Expand
    chunkft->ft(-1);
    memset( reinterpret_cast<void*>(interpolft_out), 0, sizeof(fftw_complex)*interpolft->Get_Dim_FS() );
    #pragma omp parallel for collapse(2)
    for( int i=0; i<Nx/2; i++ )
    {
      for( int j=0; j<Nyred; j++ )
      {
        interpolft_out[j+i*Nynew][0] = chunk_out[j+i*Nyred][0];
        interpolft_out[j+i*Nynew][1] = chunk_out[j+i*Nyred][1];
      }
    }

    #pragma omp parallel for collapse(2)
    for( int i=Nx/2; i<Nx; i++ )
    {
      for( int j=0; j<Nyred; j++ )
      {
        interpolft_out[j+(i+shifti)*Nynew][0] = chunk_out[j+i*Nyred][0];
        interpolft_out[j+(i+shifti)*Nynew][1] = chunk_out[j+i*Nyred][1];
      }
    }
    // interpolft->save( "fichunk_" + to_string(chunk) + ".bin", false );
    interpolft->ft(1);
    // interpolft->save( "ichunk_" + to_string(chunk) + ".bin" );

    current_chunk = chunk;
  }

  const double* Noise_Data::Get_Noise(int64_t NT){
    const int64_t interpol_nx = interpolft->Get_Dim_X();
    const int64_t interpol_ny = interpolft->Get_Dim_Y();
    const int64_t interpol_size = interpol_nx*interpol_ny;
    int64_t demanded_chunk = NT/interpol_nx;

    if (demanded_chunk != current_chunk) {
      Get_Chunk(demanded_chunk);
    }

    int64_t offset = interpol_ny*(NT - current_chunk*interpol_nx);
    assert((offset < (interpol_nx*interpol_ny)) && (offset >= 0));

    double * interpol_in = interpolft->Getp2InReal();
    const double *noise_ptr = interpol_in + offset;

    return noise_ptr;
  }


  class CRT_Propagation_1D : public CRT_Base_IF<Fourier::cft_1d,1,2>
  {
  public:
    CRT_Propagation_1D( ParameterHandler* );
    ~CRT_Propagation_1D();
  protected:
    void Do_Bragg_ad();
    void Do_Single_Noise_Step_half(fftw_complex* psi);
    void Do_Noise_Step_half(sequence_item& seq);
    void Do_Noise_Step_full(sequence_item& seq);

    void init_cn_matrix(LIS_MATRIX *A, int N, double diag, double alpha, double dx);
    void set_cn_matrix_with_metric_noise(LIS_MATRIX A, int sign);
    void lis_csr_add_value(LIS_MATRIX A, int i, int j, LIS_SCALAR val);
    void lis_csr_set_value(LIS_MATRIX A, int i, int j, LIS_SCALAR val);
    static void Do_Noise_Step_half_Wrapper(void*, sequence_item& seq);
    static void Do_Noise_Step_full_Wrapper(void*, sequence_item& seq);
    static void Do_Bragg_ad_Wrapper(void*, sequence_item& seq);

    bool run_custom_sequence( const sequence_item& item );

    // lA * psi_n+1 = rA * psi_n = b
    LIS_MATRIX m_cn_lA, m_cn_rA;
    LIS_VECTOR m_x, m_b;
    LIS_SOLVER m_solver;
    double* m_x_backup;

    vector<double> m_noise;
    vector<double> m_dx_noise;
    vector<double> m_dx2_noise;
    Fourier::rft_1d *noiseft;
    const int noise_expansion = 8;
    const int no_of_chunks = 64;
    Noise_Data *noise_data;
    double m_max_noise;

    double alpha;

    ofstream m_oftotal;
    bool no_noise_run;
  };


  CRT_Propagation_1D::CRT_Propagation_1D( ParameterHandler* p ) : CRT_Base_IF( p )
  {
    no_noise_run = (fabs(m_params->Get_Constant("Noise_Amplitude")) < 1e-8);
    if (no_noise_run) {
      cout << "No Noise Run" << endl;
    }
    m_oftotal.open("total.txt");
    double dt = 0.0;
    double duration = 0.0;
    for (auto seq_item : m_params->m_sequence) {
      dt = seq_item.dt;
      duration += seq_item.duration.front();
      cout << "Chirps = " << seq_item.no_of_chirps << endl;
    }
    int NT = duration/dt;

    m_noise.resize(m_no_of_pts);
    m_dx_noise.resize(m_no_of_pts);
    m_dx2_noise.resize(m_no_of_pts);

    if (not no_noise_run) {
      noiseft = new Fourier::rft_1d(m_header);
      noise_data = new Noise_Data(m_params->Get_simulation("NOISE"), no_of_chunks, noise_expansion);
    }

    m_max_noise = m_params->Get_Constant("Noise_Amplitude");
    cout << "Max Noise: " << m_max_noise << endl;

    m_map_stepfcts["bragg_ad"] = &Do_Bragg_ad_Wrapper;
    m_map_stepfcts["half_step"] = &Do_Noise_Step_half_Wrapper;
    m_map_stepfcts["full_step"] = &Do_Noise_Step_full_Wrapper;

    CPoint<1> pt1 = {0};
    CPoint<1> pt2 = {2*this->laser_k[0]};
    CPoint<1> pt3 = {-2*this->laser_k[0]};

    this->m_rabi_momentum_list.push_back(pt1);
    this->m_rabi_momentum_list.push_back(pt2);
    this->m_rabi_momentum_list.push_back(pt3);

    int N = this->m_no_of_pts;
    lis_vector_create(0, &m_x);
    lis_vector_set_size(m_x, 0, 2*N);
    m_x_backup = m_x->value;
    lis_vector_create(0, &m_b);
    lis_vector_set_size(m_b, 0, 2*N);

    const double dx = m_header.dx;

    alpha = m_alpha[0]*0.5;
    cout << "m_alpha: " << m_alpha[0] << '\n';
    cout << "alpha: " << alpha << endl;
    init_cn_matrix(&m_cn_rA, N, 2.0/dt, -alpha, dx);
    init_cn_matrix(&m_cn_lA, N, 2.0/dt, alpha, dx);
    lis_solver_create(&m_solver);
    char *lis_options = (char*)"-i gmres -maxiter 10000";
    lis_solver_set_option(lis_options, m_solver );

    if (no_noise_run) {
      for (int i = 0; i < m_no_of_pts; i++) {
        m_noise[i] = 0.0;
        m_dx_noise[i] = 0.0;
        m_dx2_noise[i] = 0.0;
      }
      set_cn_matrix_with_metric_noise(m_cn_rA, -1);
      set_cn_matrix_with_metric_noise(m_cn_lA, +1);
    }

  }

  CRT_Propagation_1D::~CRT_Propagation_1D() {
    m_x->value = m_x_backup;
    lis_vector_destroy(m_x);
    lis_vector_destroy(m_b);
    lis_matrix_destroy(m_cn_rA);
    lis_matrix_destroy(m_cn_lA);
    lis_solver_destroy(m_solver);

    if (not no_noise_run) {
      delete noiseft;
      delete noise_data;
    }

    m_oftotal.close();
  }

  /**
   * \brief In this function one can add additional custom sequences.
   *
   * An example can be found in Bragg_1D_double. (LINK!)
   * @param item sequence_item for sequence name
   * @retval bool true if custom sequence is found
   * */
  bool CRT_Propagation_1D::run_custom_sequence( const sequence_item& item )
  {
    // return true if a custom sequence is found or else
    if( item.name == "pseudo_beamsplitter")
    {
      vector<string> vec;
      strtk::parse(item.content,",",vec);
      CPoint<1> P;
      CPoint<1> x;

      assert( vec.size() == 1 );
      P[0] = stod(vec[0]);

      fftw_complex* Psi = m_fields[0]->Getp2In();

      for( int l=0; l<m_no_of_pts; l++ )
      {
        double re, im;
        x = this->m_fields[0]->Get_x(l);
        sincos(P*x,&im,&re);

        double re2 = Psi[l][0];
        double im2 = Psi[l][1];

        Psi[l][0] = (re2*re-im2*im+re2)/sqrt(2);
        Psi[l][1] = (re2*im+im2*re+im2)/sqrt(2);
      }
      return true;
    }
    return false;
  }

  /** Initialize Crank-Nichelson matrix
   *
   * The matrix is stored in a lis_matrix structure. Call
   * lis_matrix_destroy to deallocate the lis_matrix. A
   * five-point-stencil difference scheme is used. The formula for the
   * Crank-Nichelson matrix is A = (diag - alpha*d^2/dx^2)
   *
   * \param A  Matrix to be set
   * \param diag  Diagonal matrix part before laplace operator
   * \param alpha  Factor before laplace operator
   * \param dx  Spatial difference dx
   *
   */
  void CRT_Propagation_1D::init_cn_matrix(LIS_MATRIX *A, int N, double diag, double alpha,
                                          double dx) {
    double beta = -alpha;
    lis_matrix_create(0, A);
    lis_matrix_set_size(*A, 0, 2*N);


    lis_matrix_set_value(LIS_INS_VALUE, 2*0, (2*0)-3+(2*N), -alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*0, (2*0)-1+(2*N), 16*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*0+1, (2*0+1)-5+(2*N), -beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*0+1, (2*0+1)-3+(2*N), 16*beta/(12.0*pow(dx,2)),*A);

    lis_matrix_set_value(LIS_INS_VALUE, 2*1, (2*1)-3+(2*N), -alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*1+1, (2*1+1)-5+(2*N), -beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))+5-(2*N), -alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)+3-(2*N), -beta/(12.0*pow(dx,2)),*A);

    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-1), (2*(N-1))+3-(2*N), 16*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-1), (2*(N-1))+5-(2*N), -alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-1)+1, (2*(N-1)+1)+1-(2*N), 16*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-1)+1, (2*(N-1)+1)+3-(2*N), -beta/(12.0*pow(dx,2)),*A);

    // real
    lis_matrix_set_value(LIS_INS_VALUE, 0, 0, diag,*A);
    lis_matrix_set_value(LIS_INS_VALUE, 0, 0+1, -30*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 0, 0+3, 16*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 0, 0+5, -alpha/(12.0*pow(dx,2)),*A);
    // imag
    lis_matrix_set_value(LIS_INS_VALUE, 1, 1-1, -30*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 1, 1, diag,*A);
    lis_matrix_set_value(LIS_INS_VALUE, 1, 1+1, 16*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 1, 1+3, -beta/(12.0*pow(dx,2)),*A);
    // real
    lis_matrix_set_value(LIS_INS_VALUE, 2, 2-1, 16*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2, 2  , diag,*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2, 2+1, -30*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2, 2+3, 16*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2, 2+5, -alpha/(12.0*pow(dx,2)),*A);
    // imag
    lis_matrix_set_value(LIS_INS_VALUE, 3, (3)-3, 16*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 3, (3)-1, -30*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 3, (3)  , diag,*A);
    lis_matrix_set_value(LIS_INS_VALUE, 3, (3)+1, 16*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 3, (3)+3, -beta/(12.0*pow(dx,2)),*A);

    for (int i = 2; i < (N-2); i++) {
      // real
      lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)-3, -alpha/(12.0*pow(dx,2)),*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)-1, 16*alpha/(12.0*pow(dx,2)),*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)  , diag,*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)+1, -30*alpha/(12.0*pow(dx,2)),*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)+3, 16*alpha/(12.0*pow(dx,2)),*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)+5, -alpha/(12.0*pow(dx,2)),*A);
      // imag
      lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)-5, -beta/(12.0*pow(dx,2)),*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)-3, 16*beta/(12.0*pow(dx,2)),*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)-1, -30*beta/(12.0*pow(dx,2)),*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)  , diag,*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)+1, 16*beta/(12.0*pow(dx,2)),*A);
      lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)+3, -beta/(12.0*pow(dx,2)),*A);
    }

    // real
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))-3, -alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))-1, 16*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))  , diag,*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))+1, -30*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))+3, 16*alpha/(12.0*pow(dx,2)),*A);
    // imag
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)-5, -beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)-3, 16*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)-1, -30*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)  , diag,*A);
    lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)+1, 16*beta/(12.0*pow(dx,2)),*A);
    // real
    lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1))-3, -alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1))-1, 16*alpha/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1)), diag,*A);
    lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1))+1, -30*alpha/(12.0*pow(dx,2)),*A);
    // imag
    lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1)-5, -beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1)-3, 16*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1)-1, -30*beta/(12.0*pow(dx,2)),*A);
    lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1), diag,*A);

    lis_matrix_set_type(*A, LIS_MATRIX_CSR);
    if (lis_matrix_assemble(*A)) {
      cout << "Error in lis_matrix_assemble" << endl;
      throw;
    }
  }

  /** Add value to existing element of a lis_csr_matrix. (Unsafe)
   *
   * This helper function add a value to an existing element of a csr
   * typed lis_matrix. If you want to set a new value use
   * lis_csr_set_value instead. The matrix element has to be set via
   * lis_matrix_set_value before lis_matrix is assembled! No checking
   * is done. Memory corruption if element was not previously set...
   *
   * \param A  Matrix to be changed
   * \param row  row of the value to be changed
   * \param col  col of the value to be changed
   * \param val  value to be added to element A(row, col)
   *
   */
  void CRT_Propagation_1D::lis_csr_add_value(LIS_MATRIX A, int row, int col,
                                             LIS_SCALAR val) {
    for (int i = A->ptr[row]; i < A->ptr[row+1]; i++) {
      if (col == A->index[i]) {
        A->value[i] += val;
        return;
      }
    }
  }

  /** Set new value to existing element of a lis_csr_matrix. (Unsafe)
   *
   * This helper function sets a new value to an existing element of a
   * csr typed lis_matrix. If you want to add a value to the old one
   * use lis_csr_add_value instead. The matrix element has to be set
   * via lis_matrix_set_value before lis_matrix is assembled! No
   * checking is done. Memory corruption if element was not previously
   * set...
   *
   * \param A  Matrix to be changed
   * \param row  row of the value to be changed
   * \param col  col of the value to be changed
   * \param val  value to be added to element A(row, col)
   *
   */
  void CRT_Propagation_1D::lis_csr_set_value(LIS_MATRIX A, int row, int col,
                                             LIS_SCALAR val) {
    for (int i = A->ptr[row]; i < A->ptr[row+1]; i++) {
      if (col == A->index[i]) {
        A->value[i] = val;
        return;
      }
    }
  }

  /** Set Crank-Nichelson Matrix with metric noise terms.
   *
   * A Five-Point-Stencil difference scheme is used.
   *
   * \param A  Matrix to be set
   * \param sign  sign of laplace operator
   *
   */
  void CRT_Propagation_1D::set_cn_matrix_with_metric_noise(LIS_MATRIX A, int sign) {
    const int N = this->m_no_of_pts;
    const double alpha = sign*0.5*m_alpha[0];
    const double beta = -alpha;
    const double dx = m_header.dx;

    // Five-Point Stencel prefactors for 1st and 2nd differential
    const double fps1[] = { 1.0/(12.0*dx),
                            -8.0/(12.0*dx),
                            0.0,
                            8.0/(12.0*dx),
                            -1.0/(12.0*dx) };
    const double fps2[] = { -1.0/(12.0*pow(dx,2)),
                            16.0/(12.0*pow(dx,2)),
                            -30.0/(12.0*pow(dx,2)),
                            16.0/(12.0*pow(dx,2)),
                            -1.0/(12.0*pow(dx,2)) };

    // Local functions inserting the prefactor term for 2nd, 1st, 0th derivative
    auto d2 = [&](int i){ return (1.0 - m_noise[i]); };
    auto d1 = [&](int i){ return (-m_dx_noise[i] +
                                  0.5*(1.0 - m_noise[i])
                                  *m_dx_noise[i]); };
    auto d0 = [&](int i){ return ((1.0 - m_noise[i])
                                  *(-0.25*m_dx2_noise[i]/(1.0 + m_noise[i])
                                    + 5.0/16.0 * pow(m_dx_noise[i]
                                                     / (1.0 + m_noise[i]),2))
                                  + 0.25*pow(m_dx_noise[i], 2) / (1.0 + m_noise[i])
                                  - 0.125 * (1.0 - m_noise[i])*pow(m_dx_noise[i], 2)
                                  / (1 + m_noise[i])); };

    // Wrap around
    {
      int i = 0;
      // real
      lis_csr_set_value(A, 2*i, (2*i)-3+(2*N),
                        alpha*(fps2[0]*(d2(i-2+N) + d2(i))/2.0
                               + fps1[0]*(d1(i-2+N) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i, (2*i)-1+(2*N),
                        alpha*(fps2[1]*(d2(i-1+N) + d2(i))/2.0
                               + fps1[1]*(d1(i-1+N) + d1(i))/2.0));
      // imag
      lis_csr_set_value(A, 2*i+1, (2*i+1)-5+(2*N),
                        beta*(fps2[0]*(d2(i-2+N) + d2(i))/2.0
                              + fps1[0]*(d1(i-2+N) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i+1, (2*i+1)-3+(2*N),
                        beta*(fps2[1]*(d2(i-1+N) + d2(i))/2.0
                              + fps1[1]*(d1(i-1+N) + d1(i))/2.0));

      i = 1;
      // real
      lis_csr_set_value(A, 2*i, (2*i)-3+(2*N),
                        alpha*(fps2[0]*(d2(i-2+N) + d2(i))/2.0
                               + fps1[0]*(d1(i-2+N) + d1(i))/2.0));
      // imag
      lis_csr_set_value(A, 2*i+1, (2*i+1)-5+(2*N),
                        beta*(fps2[0]*(d2(i-2+N) + d2(i))/2.0
                              + fps1[0]*(d1(i-2+N) + d1(i))/2.0));
      i = N-2;
      // real
      lis_csr_set_value(A, 2*i, (2*i)+5-(2*N),
                        alpha*(fps2[0]*(d2(i+2-N) + d2(i))/2.0
                               + fps1[3]*(d1(i+2-N) + d1(i))/2.0));
      // imag
      lis_csr_set_value(A, 2*i+1, (2*i+1)+3-(2*N),
                        beta*(fps2[0]*(d2(i+2-N) + d2(i))/2.0
                              + fps1[3]*(d1(i+2-N) + d1(i))/2.0));

      i = N-1;
      // real
      lis_csr_set_value(A, 2*i, (2*i)+3-(2*N),
                        alpha*(fps2[1]*(d2(i+1-N) + d2(i))/2.0
                               + fps1[2]*(d1(i+1-N) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i, (2*i)+5-(2*N),
                        alpha*(fps2[0]*(d2(i+2-N) + d2(i))/2.0
                               + fps1[3]*(d1(i+2-N) + d1(i))/2.0));
      // imag
      lis_csr_set_value(A, 2*i+1, (2*i+1)+1-(2*N),
                        beta*(fps2[1]*(d2(i+1-N) + d2(i))/2.0
                              + fps1[2]*(d1(i+1-N) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i+1, (2*i+1)+3-(2*N),
                        beta*(fps2[0]*(d2(i+2-N) + d2(i))/2.0
                              + fps1[3]*(d1(i+2-N) + d1(i))/2.0));
    }

    for (int i = 0; i < 2*N; i++) {
      lis_csr_set_value(A, i, i, 2.0/m_header.dt);
    }

    // Five-Point-Stencil difference scheme
    // real
    lis_csr_set_value(A, 0, 0+1,
                      alpha*(fps2[2]*d2(0)
                             + d0(0)));
    lis_csr_set_value(A, 0, 0+3,
                      alpha*(fps2[1]*(d2(0+1) + d2(0))/2.0
                             + fps1[2]*(d1(0+1) + d1(0))/2.0));
    lis_csr_set_value(A, 0, 0+5,
                      alpha*(fps2[0]*(d2(0+2) + d2(0))/2.0
                             + fps1[3]*(d1(0+2) + d1(0))/2.0));
    // imag
    lis_csr_set_value(A, 1, 1-1,
                      beta*(fps2[2]*d2(0)
                            + d0(0)));
    lis_csr_set_value(A, 1, 1+1,
                      beta*(fps2[1]*(d2(0+1) + d2(0))/2.0
                            + fps1[2]*(d1(0+1) + d1(0))/2.0));
    lis_csr_set_value(A, 1, 1+3,
                      beta*(fps2[0]*(d2(0+2) + d2(0))/2.0
                            + fps1[3]*(d1(0+2) + d1(0))/2.0));
    // real
    lis_csr_set_value(A, 2, 2-1,
                      alpha*(fps2[1]*(d2(1-1) + d2(1))/2.0
                             + fps1[1]*(d1(1-1) + d1(1))/2.0));
    lis_csr_set_value(A, 2, 2+1,
                      alpha*(fps2[2]*d2(1)
                             + d0(1)));
    lis_csr_set_value(A, 2, 2+3,
                      alpha*(fps2[1]*(d2(1+1) + d2(1))/2.0
                             + fps1[2]*(d1(1+1) + d1(1))/2.0));
    lis_csr_set_value(A, 2, 2+5,
                      alpha*(fps2[0]*(d2(1+2) + d2(1))/2.0
                             + fps1[3]*(d1(1+2) + d1(1))/2.0));
    // imag
    lis_csr_set_value(A, 3, 3-3,
                      beta*(fps2[1]*(d2(1-1) + d2(1))/2.0
                            + fps1[1]*(d1(1-1) + d1(1))/2.0));
    lis_csr_set_value(A, 3, 3-1,
                      beta*(fps2[2]*d2(1)
                            + d0(1)));
    lis_csr_set_value(A, 3, 3+1,
                      beta*(fps2[1]*(d2(1+1) + d2(1))/2.0
                            + fps1[2]*(d1(1+1) + d1(1))/2.0));
    lis_csr_set_value(A, 3, 3+3,
                      beta*(fps2[0]*(d2(1+2) + d2(1))/2.0
                            + fps1[3]*(d1(1+2) + d1(1))/2.0));

    #pragma omp parallel for
    for (int i = 2; i < (N-2); i++) {
      // real
      lis_csr_set_value(A, 2*i, (2*i)-3,
                        alpha*(fps2[0]*(d2(i-2) + d2(i))/2.0
                               + fps1[0]*(d1(i-2) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i, (2*i)-1,
                        alpha*(fps2[1]*(d2(i-1) + d2(i))/2.0
                               + fps1[1]*(d1(i-1) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i, (2*i)+1,
                        alpha*(fps2[2]*d2(i)
                               + d0(i)));
      lis_csr_set_value(A, 2*i, (2*i)+3,
                        alpha*(fps2[3]*(d2(i+1) + d2(i))/2.0
                               + fps1[3]*(d1(i+1) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i, (2*i)+5,
                        alpha*(fps2[4]*(d2(i+2) + d2(i))/2.0
                               + fps1[4]*(d1(i+2) + d1(i))/2.0));
      // imag
      lis_csr_set_value(A, 2*i+1, (2*i+1)-5,
                        beta*(fps2[0]*(d2(i-2) + d2(i))/2.0
                              + fps1[0]*(d1(i-2) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i+1, (2*i+1)-3,
                        beta*(fps2[1]*(d2(i-1) + d2(i))/2.0
                              + fps1[1]*(d1(i-1) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i+1, (2*i+1)-1,
                        beta*(fps2[2]*d2(i)
                              + d0(i)));
      lis_csr_set_value(A, 2*i+1, (2*i+1)+1,
                        beta*(fps2[3]*(d2(i+1) + d2(i))/2.0
                              + fps1[3]*(d1(i+1) + d1(i))/2.0));
      lis_csr_set_value(A, 2*i+1, (2*i+1)+3,
                        beta*(fps2[4]*(d2(i+2) + d2(i))/2.0
                              + fps1[4]*(d1(i+2) + d1(i))/2.0));
    }

    // real
    lis_csr_set_value(A, (2*(N-2)), (2*(N-2))-3,
                      alpha*(fps2[0]*(d2((N-2)-2) + d2((N-2)))/2.0
                             + fps1[0]*(d1((N-2)-2) + d1((N-2)))/2.0));
    lis_csr_set_value(A, (2*(N-2)), (2*(N-2))-1,
                      alpha*(fps2[1]*(d2((N-2)-1) + d2((N-2)))/2.0
                             + fps1[1]*(d1((N-2)-1) + d1((N-2)))/2.0));
    lis_csr_set_value(A, (2*(N-2)), (2*(N-2))+1,
                      alpha*(fps2[2]*d2(N-2)
                             + d0(N-2)));
    lis_csr_set_value(A, (2*(N-2)), (2*(N-2))+3,
                      alpha*(fps2[1]*(d2(N-2+1) + d2(N-2))/2.0
                             + fps1[2]*(d1(N-2+1) + d1(N-2))/2.0));
    // imag
    lis_csr_set_value(A, (2*(N-2)+1), (2*(N-2)+1)-5,
                      beta*(fps2[0]*(d2(N-2-2) + d2(N-2))/2.0
                            + fps1[0]*(d1(N-2-2) + d1(N-2))/2.0));
    lis_csr_set_value(A, (2*(N-2)+1), (2*(N-2)+1)-3,
                      beta*(fps2[1]*(d2(N-2-2) + d2(N-2))/2.0
                            + fps1[1]*(d1(N-2-2) + d1(N-2))/2.0));
    lis_csr_set_value(A, (2*(N-2)+1), (2*(N-2)+1)-1,
                      beta*(fps2[2]*d2(N-2)
                            + d0(N-2)));
    lis_csr_set_value(A, (2*(N-2)+1), (2*(N-2)+1)+1,
                      beta*(fps2[1]*(d2(N-2+1) + d2(N-2))/2.0
                            + fps1[2]*(d1(N-2+1) + d1(N-2))/2.0));
    // real
    lis_csr_set_value(A, (2*(N-1)), (2*(N-1))-3,
                      alpha*(fps2[0]*(d2(N-1-2) + d2(N-1))/2.0
                             + fps1[0]*(d1(N-1-2) + d1(N-1))/2.0));
    lis_csr_set_value(A, (2*(N-1)), (2*(N-1))-1,
                      alpha*(fps2[1]*(d2(N-1-1) + d2(N-1))/2.0
                             + fps1[1]*(d1(N-1-1) + d1(N-1))/2.0));
    lis_csr_set_value(A, (2*(N-1)), (2*(N-1))+1,
                      alpha*(fps2[2]*d2(N-1)
                             + d0(N-1)));
    // imag
    lis_csr_set_value(A, (2*(N-1)+1), (2*(N-1)+1)-5,
                      beta*(fps2[0]*(d2(N-1-2) + d2(N-1))/2.0
                            + fps1[0]*(d1(N-1-2) + d1(N-1))/2.0));
    lis_csr_set_value(A, (2*(N-1)+1), (2*(N-1)+1)-3,
                      beta*(fps2[1]*(d2(N-1-1) + d2(N-1))/2.0
                            + fps1[1]*(d1(N-1-1) + d1(N-1))/2.0));
    lis_csr_set_value(A, (2*(N-1)+1), (2*(N-1)+1)-1,
                      beta*(fps2[2]*d2(N-1)
                            + d0(N-1)));
  }

  /** Wrapper function for Do_Noise_Step_half.
   *
   * \param ptr  Function pointer to be set to Do_Noise_Step_half
   * \param seq  sequence (Not used)
   *
   */
  void CRT_Propagation_1D::Do_Noise_Step_half_Wrapper ( void* ptr,
                                                        sequence_item& seq )
  {
    CRT_Propagation_1D *self = static_cast<CRT_Propagation_1D*>(ptr);
    self->Do_Noise_Step_half(seq);
  }

  /** Half step with metric noise is performed for a single wavefunction.
   *
   * The Half step is perfomed by solving the Crank-Nichelson system
   * via LIS routines.
   *
   * \param psi  Wave function to be evolved
   *
   */
  void CRT_Propagation_1D::Do_Single_Noise_Step_half(fftw_complex* psi) {
    int N = this->m_no_of_pts;
    LIS_INT err,iter,iter_double,iter_quad;
    double time,itime,ptime,p_c_time,p_i_time;
    LIS_REAL resid;

    static int myiter = 0;
    m_x->value = reinterpret_cast<double*>(psi);
    err = lis_matvec(m_cn_rA, m_x, m_b);

    // if (myiter < 3) {
    //   char * my_argument = const_cast<char*> (("rmatrix-"+to_string(myiter)+".txt").c_str());
    //   lis_output_matrix(m_cn_rA, LIS_FMT_MM, my_argument);
    //   char * my_argument2 = const_cast<char*> (("lmatrix-"+to_string(myiter)+".txt").c_str());
    //   lis_output_matrix(m_cn_lA, LIS_FMT_MM, my_argument2);
    //   myiter++;
    // }

    if (err) {
      cout << "ERROR: lis_matvec" << endl;
      exit(EXIT_FAILURE);
    }

    for (int i = 0; i < 2*N; i++) {
      m_x->value[i] = 0.0;
    }

    err=lis_solve(m_cn_lA, m_b, m_x, m_solver);
    CHKERR(err);
    if (err) {
      cout << "ERROR: lis_solve" << endl;
      exit(EXIT_FAILURE);
    }

    lis_solver_get_iterex(m_solver,&iter,&iter_double,&iter_quad);
    lis_solver_get_timeex(m_solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
    lis_solver_get_residualnorm(m_solver,&resid);
    if (iter > 1000) {
      cout << "Error: Iter too high" << endl;
      cout << "Time: " << m_header.t << ", Iter: " << iter << endl;
      exit(EXIT_FAILURE);
    }
  }


  /** Half step with metric noise is performed on all wavefunctions.
   *
   */
  void CRT_Propagation_1D::Do_Noise_Step_half(sequence_item& seq) {
    bool update_noise = fmod(round(fabs(10*m_header.t/m_header.dt)), 10) == 0;
    // Sanity Check, bool should alternate
    static bool last_bool = false;
    assert(last_bool != update_noise);
    last_bool = update_noise;

    if ((not no_noise_run) && update_noise) // update noise every fullstep
    {
      for (int i = 0; i < m_no_of_pts; i++) {
        m_noise[i] = 0.0;
      }

      int64_t NT = static_cast<int64_t>(floor(m_header.t/noise_data->dt+0.5));
      int nSteps = static_cast<int>(floor(seq.dt/noise_data->dt+0.5));
      assert(nSteps > 0);

      for (int j = 0; j < nSteps; j++) {
        const double* noise = noise_data->Get_Noise(NT+j);
        for (int i = 0; i < m_no_of_pts; i++) {
          m_noise[i] += m_max_noise*noise[i]/nSteps;
        }
      }

      double *diff_noise = noiseft->Getp2InReal();
      for (int i = 0; i < m_no_of_pts; i++) {
        diff_noise[i] = m_noise[i];
      }
      noiseft->Diff_x();
      for (int i = 0; i < m_no_of_pts; i++) {
        m_dx_noise[i] = diff_noise[i];
      }
      noiseft->Diff_x();
      for (int i = 0; i < m_no_of_pts; i++) {
        m_dx2_noise[i] = diff_noise[i];
      }

      // set_cn_matrix_with_metric_noise(m_cn_rA, -1);
      // set_cn_matrix_with_metric_noise(m_cn_lA, +1);
    }
    if (update_noise) {
      set_cn_matrix_with_metric_noise(m_cn_rA, -1);
      set_cn_matrix_with_metric_noise(m_cn_lA, +1);
    }

    double total = 0;
    int no_int_states = 2;
    m_oftotal << setprecision(12);
    for( int c=0; c<no_int_states; c++ ) {
      double nParticles = this->Get_Particle_Number(c);
      m_oftotal << nParticles << "\t";
      total += nParticles;
    }
    m_oftotal << total << endl;

    Do_Single_Noise_Step_half(this->m_fields[0]->Getp2In());
    Do_Single_Noise_Step_half(this->m_fields[1]->Getp2In());

    m_header.t += 0.5*m_header.dt;
  }

  /** Wrapper function for Do_Noise_Step_full.
   *
   * \param ptr  Function pointer to be set to Do_Noise_Step_full
   * \param seq  sequence (not used)
   */
  void CRT_Propagation_1D::Do_Noise_Step_full_Wrapper ( void* ptr,
                                                        sequence_item& seq )
  {
    CRT_Propagation_1D *self = static_cast<CRT_Propagation_1D*>(ptr);
    self->Do_Noise_Step_full(seq);
  }

  /** Full step with metric noise is performed on all wavefunctions.
   *
   * The full step is perfomed by performing two consecutive half steps.
   *
   */
  void CRT_Propagation_1D::Do_Noise_Step_full(sequence_item& seq) {
    Do_Noise_Step_half(seq);
    Do_Noise_Step_half(seq);
  }

  /** Wrapper function for Do_Bragg_ad.
   *
   * \param ptr  Function pointer to be set to Do_Bragg_ad
   * \param seq  sequence (not used)
   */
  void CRT_Propagation_1D::Do_Bragg_ad_Wrapper ( void* ptr,
                                                 sequence_item& seq )
  {
    ignore = seq;
    CRT_Propagation_1D *self = static_cast<CRT_Propagation_1D*>(ptr);
    self->Do_Bragg_ad();
  }

  /**
   * \brief This function computes the laser-atom interaction by means of an analytical diagonalisation
   *
   * See XXX for further information about the method used for the diagonalisation
   * */
  void CRT_Propagation_1D::Do_Bragg_ad() // ad -> analytic diagonalization
  {
    /*
      vector<fftw_complex*> Psi;
      Psi.push_back(m_fields[0]->Getp2In());
      Psi.push_back(m_fields[1]->Getp2In());

      fftw_complex* Psi_1 = this->m_fields[0]->Getp2In();
      fftw_complex* Psi_2 = this->m_fields[1]->Getp2In();

      const double dt = this->Get_dt();
      const double t1 = this->Get_t()+0.5*dt;

      double re1, im1, tmp1, tmp2, V11, V22, Ep, Em, Omega;

      CPoint<1> x;
      array<complex<double>,2> gamma;
      complex<double> Psi_alt_0, Psi_alt_1, Psi_neu_0, Psi_neu_1, M00, M01, M10, M11, eta;

      #pragma omp parallel for private(x,gamma,re1,im1,tmp1,tmp2,V11,V22,Ep,Em,Omega,Psi_alt_0,Psi_alt_1,Psi_neu_0,Psi_neu_1,M00,M01,M10,M11,eta)
      for( int l=0; l<this->m_no_of_pts; l++ )
      {
      x = this->m_fields[0]->Get_x(l);

      tmp1 = Psi_1[l][0]*Psi_1[l][0]+Psi_1[l][1]*Psi_1[l][1];
      tmp2 = Psi_2[l][0]*Psi_2[l][0]+Psi_2[l][1]*Psi_2[l][1];

      V11 = this->m_gs[0]*tmp1+this->m_gs[1]*tmp2+beta*x;
      V22 = this->m_gs[2]*tmp1+this->m_gs[3]*tmp2-DeltaL[1]+beta*x;

      Omega = Amp[0]*cos(laser_k[0]*x[0]-(laser_domh[0]+chirp_rate[0]*t1)*t1+phase[0]/2);

      Psi_alt_0 = complex<double>( Psi[0][l][0], Psi[0][l][1]);
      Psi_alt_1 = complex<double>( Psi[1][l][0], Psi[1][l][1]);

      if( Omega == 0.0 )
      {
      Psi_neu_0 = exp(complex<double>( -dt*V11, 0)) * Psi_alt_0;
      Psi[0][l][0] = real(Psi_neu_0);
      Psi[0][l][1] = imag(Psi_neu_0);

      Psi_neu_1 = exp(complex<double>( -dt*V22, 0)) * Psi_alt_1;
      Psi[1][l][0] = real(Psi_neu_1);
      Psi[1][l][1] = imag(Psi_neu_1);
      continue;
      }

      eta = exp(complex<double>(0,-0.5*(laser_dk[0]*x[0]+phase[0])));

      tmp1 = sqrt((V11-V22)*(V11-V22)+4.0*Omega*Omega);
      Ep = 0.5*(V11+V22+tmp1);
      Em = 0.5*(V11+V22-tmp1);

      gamma[0] = exp(complex<double>( 0, -dt*Ep )) / fabs((V11-Ep)*(V11-Ep)+Omega*Omega);
      gamma[1] = exp(complex<double>( 0, -dt*Em )) / fabs((V11-Em)*(V11-Em)+Omega*Omega);

      tmp1 = V11-Ep;
      tmp2 = V11-Em;

      M00 = Omega*Omega*(gamma[0]+gamma[1]);
      M01 = -Omega*eta*(tmp1*gamma[0]+tmp2*gamma[1]);
      M10 = -Omega*conj(eta)*(tmp1*gamma[0]+tmp2*gamma[1]);
      M11 = tmp1*tmp1*gamma[0] + tmp2*tmp2*gamma[1];

      Psi_neu_0 = M00 * Psi_alt_0 + M01 * Psi_alt_1;
      Psi_neu_1 = M10 * Psi_alt_0 + M11 * Psi_alt_1;

      Psi[0][l][0] = real(Psi_neu_0);
      Psi[0][l][1] = imag(Psi_neu_0);
      Psi[1][l][0] = real(Psi_neu_1);
      Psi[1][l][1] = imag(Psi_neu_1);
      }
      }*/
    fftw_complex* Psi_1 = this->m_fields[0]->Getp2In();
    fftw_complex* Psi_2 = this->m_fields[1]->Getp2In();

    const double dt = this->Get_dt();
    const double t1 = this->Get_t()+0.5*dt;
    fftw_complex O11, O12, O21, O22, gamma_p, gamma_m, eta;
    double re1, im1, Ep, Em;

    CPoint<1> x;
    #pragma omp parallel for private(x,re1,im1,Ep,Em,O11,O12,O21,O22,gamma_p,gamma_m,eta)
    for( int l=0; l<this->m_no_of_pts; l++ )
    {
      x = this->m_fields[0]->Get_x(l);

      double tmp1 = Psi_1[l][0]*Psi_1[l][0]+Psi_1[l][1]*Psi_1[l][1];
      double tmp2 = Psi_2[l][0]*Psi_2[l][0]+Psi_2[l][1]*Psi_2[l][1];

      double V11 = this->m_gs[0]*tmp1+this->m_gs[1]*tmp2+beta*x;
      double V22 = this->m_gs[2]*tmp1+this->m_gs[3]*tmp2-DeltaL[1]+beta*x;

      double Omega = Amp[0]*cos(laser_k[0]*x[0]-(laser_domh[0]+chirp_rate[0]*t1)*t1+phase[0]/2);

      if( Omega == 0.0 )
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
  }
}

int main( int argc, char* argv[] )
{
  if( argc != 2 )
  {
    cout << "No parameter xml file specified." << endl;
    cout << argv[0] << " freeprop.xml" << endl;
    return EXIT_FAILURE;
  }

  ParameterHandler params(argv[1]);

  int no_of_threads = 4;
  char* envstr = getenv( "MY_NO_OF_THREADS" );
  if( envstr != nullptr ) no_of_threads = atoi( envstr );

  cout << "Number of threads: " << no_of_threads << endl;
  fftw_init_threads();
  fftw_plan_with_nthreads( no_of_threads );
  omp_set_num_threads( no_of_threads );

  lis_initialize(&argc, &argv);
  sequence_item seq;

  try
  {
    RT_Solver::CRT_Propagation_1D rtsol( &params );
    rtsol.run_sequence();
  }
  catch(mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }
  catch(string &str)
  {
    cout << str << endl;
  }

  fftw_cleanup_threads();
  lis_finalize();
  return EXIT_SUCCESS;
}
