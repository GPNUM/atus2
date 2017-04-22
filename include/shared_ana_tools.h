// External call: Memory = (no_of_momentum_states + 1) * sizeof(Psi)
// Internal call: Memory = (no_of_momentum_states) * sizeof(Psi)

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <dirent.h>
#include <cmath>

#include "fftw3.h"
#include "my_structs.h"
#include "CPoint.h"
#include "ParameterHandler.h"

//For double and CPoint
template <class T>
struct ana_data
{
  double time = 0;
  T data;
  std::vector<T> state_data;

  ana_data()
  {
    data = 0;
  };

  bool operator < (const ana_data &str) const
  {
    return (time < str.time);
  }
};

template <int dim, class T>
class Shared_Ana_Tools
{
public:
  Shared_Ana_Tools ( std::string, ParameterHandler * );
  Shared_Ana_Tools ( const generic_header &, ParameterHandler * );
  Shared_Ana_Tools ( ParameterHandler * );
  ~Shared_Ana_Tools();

  void Set_Wavefunction( T *, const double );

  void Write( int seq = 0 );
  void Run_Analysis();
  void Run_Analysis_in_Directory();

  const double &Get_Particle_Number ()
  {
    Particle_Number( false );
    return m_no_of_particles.back().data;
  }
  const std::vector<double> &Get_Particle_Number_States ()
  {
    Particle_Number( true );
    return m_no_of_particles.back().state_data;
  }

  const CPoint<dim> &Get_Expval_Position ()
  {
    Expval_Position( false );
    return m_exp_pos.back().data;
  }
  const std::vector<CPoint<dim>> &Get_Expval_Position_States ()
  {
    Expval_Position( true );
    return m_exp_pos.back().state_data;
  }

  const CPoint<dim> &Get_Expval_Momentum ()
  {
    Expval_Momentum( false );
    return m_exp_mom.back().data;
  }
  const std::vector<CPoint<dim>> &Get_Expval_Momentum_States ()
  {
    Expval_Momentum( true );
    return m_exp_mom.back().state_data;
  }

protected:
  //Function pointer
  typedef void (Shared_Ana_Tools::*fctptr)( const bool );

  //For storing data in a vector
  template <class P>
  using Data = std::vector<ana_data<P>>;

  void Read_File ( std::string );
  void Read_Header();
  void Init();
  void Init_fcts();

  template <class P>
  bool Data_Already_Calculated( Data<P> &, const bool );

  void Momentum_States_Pos( );

  void Particle_Number( const bool );
  void Expval_Position( const bool );
  void Expval_Momentum( const bool );
  void Save_Momentum_States_Pos( const bool );
  void Save_Momentum_States( const bool );
  void FWHM( const bool );
  void Phase( const bool );

  T *m_ft;
  generic_header m_header;

  long long m_no_of_pts;
  ParameterHandler *m_params;
  std::vector<analyze_item> m_options;

  std::vector<CPoint<dim>> m_rabi_momentum_list;
  std::map<std::string,fctptr> m_map_fcts;
  std::vector<T *>     m_momentum_states;

  double              m_ar,
                      m_ar_k,
                      m_rabi_threshold;
  Data<double>        m_no_of_particles;
  Data<CPoint<dim>>   m_exp_pos,
       m_exp_mom;
};

template<int dim, class T>
Shared_Ana_Tools<dim,T>::Shared_Ana_Tools( ParameterHandler *params ) : m_params(params)
{
  try
  {
    Read_File(m_params->Get_simulation("FILENAME"));
  }
  catch ( std::string &str )
  {
    std::cout << str << std::endl;
  }
  Init();
  Init_fcts();
}

template<int dim, class T>
Shared_Ana_Tools<dim,T>::Shared_Ana_Tools( const generic_header &header, ParameterHandler *params ) : m_params(params), m_header(header)
{
  Init();
  Init_fcts();
}

template<int dim, class T>
Shared_Ana_Tools<dim,T>::Shared_Ana_Tools( std::string filename, ParameterHandler *params ) : m_params(params)
{
  Read_File( filename );

  Init();
  Init_fcts();
}

template <int dim, class T>
Shared_Ana_Tools <dim,T>::~Shared_Ana_Tools()
{
  delete m_ft;
  for ( auto state : m_momentum_states)
    delete state;
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Set_Wavefunction( T *ft, const double time )
{
  m_ft = ft;
  m_header.t = time;

  for ( auto state : m_momentum_states)
    delete state;
  m_momentum_states.clear();
}

template <int dim, class T>
template <typename P>
bool Shared_Ana_Tools <dim,T>::Data_Already_Calculated( Data<P> &vec, const bool option )
{
  ana_data<P> tmp;
  tmp.time = m_header.t;

  if ( vec.empty() )
  {
    vec.push_back(tmp);
  }

  if ( vec.back().time != m_header.t )
  {
    vec.push_back(tmp);
  }

  if ( vec.back().data == 0 && option == false)
    return false;
  if ( vec.back().state_data.empty() && option == true)
    return false;

  return true;
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Write( int seq )
{
  if ( ! m_no_of_particles.empty() )
  {
    std::sort( m_no_of_particles.begin(), m_no_of_particles.end() );

    std::string file;
    file = "Particle_Number_" + std::to_string(seq) + ".txt";
    std::ofstream myfile( file );
    for ( auto it : m_no_of_particles )
    {
      myfile << it.time << "\t" << it.data << "\t";
      for ( auto it2 : it.state_data )
        myfile << it2 << "\t";
      myfile << "\n";
    }
    myfile.close();
  }

  if ( !m_exp_pos.empty() )
  {
    std::sort( m_exp_pos.begin(), m_exp_pos.end() );

    std::string file;
    file = "Expectation_Pos_" + std::to_string(seq) + ".txt";
    std::ofstream myfile( file );
    for ( auto it : m_exp_pos )
    {
      myfile << it.time << "\t" << it.data << "\t";
      for ( auto it2 : it.state_data )
        myfile << it2 << "\t";
      myfile << "\n";
    }
    myfile.close();
  }

  if ( !m_exp_mom.empty() )
  {
    std::sort( m_exp_mom.begin(), m_exp_mom.end() );

    std::string file;
    file = "Expectation_Mom_" + std::to_string(seq) + ".txt";
    std::ofstream myfile( file );
    for ( auto it : m_exp_mom )
    {
      myfile << it.time << "\t" << it.data << "\t";
      for ( auto it2 : it.state_data )
        myfile << it2 << "\t";
      myfile << "\n";
    }
    myfile.close();
  }

}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Run_Analysis_in_Directory()
{
  DIR *dir;
  struct dirent *files;
  dir=opendir(".");
  while ( (files = readdir(dir)) )
  {
    std::string file = files->d_name;

    if ( file.size() < 10 ) continue;
    if ( file.compare( file.size()-6, 6, "_1.bin" ) == 0 )
    {
      for ( auto state : m_momentum_states)
        delete state;
      m_momentum_states.clear();
      Read_File(file);
      Run_Analysis();
    }
  }
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Read_File( std::string filename )
{

  std::ifstream myfile (filename, std::ifstream::binary);
  myfile.read ((char *) &m_header, sizeof(generic_header));

  m_ft = new T (m_header);
  double m_no_of_pts = m_header.nDimX*m_header.nDimY*m_header.nDimZ;

  myfile.seekg( sizeof(generic_header));
  myfile.read((char *)m_ft->Getp2In(),sizeof(fftw_complex)*m_no_of_pts);
  myfile.close();
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Init()
{
  m_no_of_pts = m_header.nDimX*m_header.nDimY*m_header.nDimZ;

  switch ( dim )
  {
  case 1:
    m_ar   = m_header.dx;
    m_ar_k = m_header.dkx;
    break;
  case 2:
    m_ar   = m_header.dx * m_header.dy;
    m_ar_k = m_header.dkx * m_header.dky;
    break;
  case 3:
    m_ar   = m_header.dx * m_header.dy * m_header.dz;
    m_ar_k = m_header.dkx * m_header.dky * m_header.dkz;
    break;
  }
  try
  {
    m_rabi_threshold = m_params->Get_Constant("rabi_threshold");

    int p = 1;
    std::vector<double> vec;
    CPoint<dim> state;
    while (m_params->Get_MomentumStates(vec, p))
    {
      if ( vec.size() != dim ) throw std::string("Error: Wrong dimension of momentum states");
      state = vec;
      m_rabi_momentum_list.push_back(state);
      p++;
    }
  }
  catch ( std::string &str )
  {
    std::cout << str << std::endl;
  }
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Init_fcts()
{
  m_map_fcts["particle_number"] = &Shared_Ana_Tools::Particle_Number;
  m_map_fcts["expval_position"] = &Shared_Ana_Tools::Expval_Position;
  m_map_fcts["expval_momentum"] = &Shared_Ana_Tools::Expval_Momentum;
  m_map_fcts["mom_states_pos"]  = &Shared_Ana_Tools::Save_Momentum_States_Pos;
  m_map_fcts["mom_states"]      = &Shared_Ana_Tools::Save_Momentum_States;
  m_map_fcts["fwhm"]            = &Shared_Ana_Tools::FWHM;
  m_map_fcts["phase"]           = &Shared_Ana_Tools::Phase;
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Momentum_States_Pos( )
{
  if ( ! m_momentum_states.empty() ) return;

  if ( m_rabi_momentum_list.empty() )
  {
    std::cout << "Error: No momentum states specified." << std::endl;
    throw;
  }


  for ( auto pos_of_states : m_rabi_momentum_list )
  {
    fftw_complex *Psi;
    CPoint<dim> k, d;
    m_momentum_states.push_back( new T(m_header) );

    auto it = --m_momentum_states.end();
    Psi = (*it)->Getp2In();
    memcpy( Psi, m_ft->Getp2In(), m_no_of_pts*sizeof(fftw_complex) );
    (*it)->ft(-1);

    for ( int i=0; i<m_no_of_pts; i++ )
    {
      k = m_ft->Get_k(i);
      d = pos_of_states - k;
      if ( sqrt(d*d) > m_rabi_threshold )
      {
        Psi[i][0] = 0;
        Psi[i][1] = 0;
      }
    }
    (*it)->ft(1);
  }
}


template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Particle_Number( const bool option)
{
  // Total number of particles
  if ( ! option )
  {
    if ( Data_Already_Calculated(m_no_of_particles,option) ) return;

    fftw_complex *Psi = m_ft->Getp2In();
    double Particles = 0;

    for (int i=0; i<m_no_of_pts; i++)
      Particles += Psi[i][0]*Psi[i][0]+Psi[i][1]*Psi[i][1];

    Particles *= m_ar;
    m_no_of_particles.back().data = Particles;
  }

  // Number of particles in momentum states
  if ( option )
  {
    if ( Data_Already_Calculated(m_no_of_particles,option) ) return;

    // If empty calculate momentum states
    if ( m_momentum_states.empty() ) Momentum_States_Pos();

    for ( auto field : m_momentum_states )
    {
      fftw_complex *Psi = field->Getp2In();
      double Particles = 0;
      for ( int i=0; i<m_no_of_pts; i++)
        Particles += Psi[i][0]*Psi[i][0]+Psi[i][1]*Psi[i][1];

      Particles *= m_ar;
      m_no_of_particles.back().state_data.push_back( Particles );
    }
  }
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Expval_Position( const bool option )
{
  CPoint<dim> x;

  // Total postion expectation value
  if ( ! option )
  {
    if ( Data_Already_Calculated(m_exp_pos, option) ) return;

    fftw_complex *Psi = m_ft->Getp2In();
    CPoint<dim> tmp;
    for ( int i=0; i<m_no_of_pts; i++)
    {
      x = m_ft->Get_x(i);
      tmp += (x*(Psi[i][0]*Psi[i][0]+Psi[i][1]*Psi[i][1]));
    }
    tmp *= m_ar/Get_Particle_Number();
    m_exp_pos.back().data = tmp;
  }

  // Expectation value for momentum states
  if ( option )
  {
    if ( Data_Already_Calculated(m_exp_pos, option) ) return;
    // If empty calculate momentum states
    if ( m_momentum_states.empty() ) Momentum_States_Pos();

    std::vector<double> N = Get_Particle_Number_States();
    int j = 0;

    for ( auto field : m_momentum_states )
    {
      fftw_complex *Psi = field->Getp2In();
      CPoint<dim> tmp;
      for ( int i=0; i<m_no_of_pts; i++)
      {
        x = field->Get_x(i);
        tmp += (x*(Psi[i][0]*Psi[i][0]+Psi[i][1]*Psi[i][1]));
      }
      tmp *= m_ar / N[j++];
      m_exp_pos.back().state_data.push_back(tmp);
    }
  }
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Expval_Momentum( const bool option )
{
  if ( ! option )
  {
    if ( Data_Already_Calculated(m_exp_mom,option) ) return;

    fftw_complex *Psi = m_ft->Getp2In();
    double N = Get_Particle_Number();

    m_ft->ft(-1);

    CPoint<dim> tmp, k;
    for ( int i=0; i<m_no_of_pts; i++)
    {
      k = m_ft->Get_k(i);
      tmp += (k*(Psi[i][0]*Psi[i][0]+Psi[i][1]*Psi[i][1]));
    }
    tmp *= m_ar_k/N;
    m_exp_mom.back().data = tmp;

    m_ft->ft(1);
  }

  if ( option )
  {
    if ( Data_Already_Calculated(m_exp_mom,option) ) return;
    // If empty calculate momentum states
    if ( m_momentum_states.empty() ) Momentum_States_Pos();

    std::vector<double> N = Get_Particle_Number_States();
    int j = 0;

    for ( auto field : m_momentum_states )
    {
      fftw_complex *Psi = field->Getp2In();

      field->ft(-1);

      CPoint<dim> tmp, k;
      for ( int i=0; i<m_no_of_pts; i++)
      {
        k = field->Get_k(i);
        tmp += (k*(Psi[i][0]*Psi[i][0]+Psi[i][1]*Psi[i][1]));
      }
      tmp *= m_ar_k / N[j++];
      m_exp_mom.back().state_data.push_back(tmp);

      field->ft(1);
    }
  }
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Save_Momentum_States_Pos( const bool )
{
  if ( m_momentum_states.empty() ) Momentum_States_Pos();

  char *header = reinterpret_cast<char *> ( &m_header );
  int state = 1;

  for ( auto field : m_momentum_states )
  {
    //std::string filename;
    //filename = std::to_string( m_header.t ) + "_1_M_" + std::to_string( state++ ) + ".bin";
    char filename[1024];
    sprintf( filename, "%.3f_1_M%d.bin", m_header.t, state++ );
    char *Psi = reinterpret_cast<char *> ( field->Getp2In() );

    std::ofstream file1( filename, std::ofstream::binary | std::ofstream::trunc );
    file1.write( header, sizeof(generic_header) );
    file1.write( Psi, m_no_of_pts*sizeof(fftw_complex) );
    file1.close();
  }
}

//TODO noch anpassen (falscher header)
template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Save_Momentum_States( const bool )
{
  //TODO sollte mit setfix sein
  //m_ft->SetFix(true);
  m_ft->ft(-1);

  char filename[1024];
  sprintf( filename, "%.3f_ana.bin", m_header.t );

  char *header = reinterpret_cast<char *> ( &m_header );
  char *Psi = reinterpret_cast<char *> (m_ft->Getp2In() );

  std::ofstream file1( filename, std::ofstream::binary | std::ofstream::trunc );
  file1.write( header, sizeof(generic_header) );
  file1.write( Psi, m_no_of_pts*sizeof(fftw_complex) );
  file1.close();

  m_ft->ft(1);
  //m_ft->SetFix(false);
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::FWHM( const bool option )
{
  (void)option;
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Phase( const bool option )
{
  double calc_threshold = 0.1;

  if ( ! option )
  {
  }

  if ( option )
  {
    char *header = reinterpret_cast<char *> ( &m_header );

    if ( m_momentum_states.empty() ) Momentum_States_Pos();

    int state = 1;
    int j = 0;

    std::vector<CPoint<dim>> px = Get_Expval_Momentum_States();
    double im, re, re2, im2;
    CPoint<dim> x;
    m_header.bComplex = 0;

    for ( auto field : m_momentum_states )
    {
      double Phase[m_no_of_pts] ;
      fftw_complex *Psi = field->Getp2In();

      for ( int l=0; l<m_no_of_pts; l++ )
      {
        //Subtract exp momentum
        x = field->Get_x(l);
        sincos((-1)*(px[j]*x),&im,&re);
        re2 = Psi[l][0]*re-Psi[l][1]*im;
        im2 = Psi[l][0]*im+Psi[l][1]*re;

        if ( re2*re2+im2*im2 > calc_threshold )
          Phase[l] = std::atan2(im2,re2);
        else
          Phase[l] = 0;
      }

      j++;
      char filename[1024];
      sprintf( filename, "%.3f_1_Phase_%d.bin", m_header.t, state++ );
      char *Phi = reinterpret_cast<char *> ( Phase );

      std::ofstream file1( filename, std::ofstream::binary | std::ofstream::trunc);
      file1.write( header, sizeof(generic_header) );
      file1.write( Phi, m_no_of_pts*sizeof(double) );
      file1.close();
    }
    m_header.bComplex = 1;
  }
}

template <int dim, class T>
void Shared_Ana_Tools <dim,T>::Run_Analysis()
{
  fctptr step_fct=nullptr;

  for ( auto option : m_params->m_analyze )
  {
    if ( option.content == false )
      continue;
    try
    {
      step_fct = this->m_map_fcts.at(option.name);
    }
    catch (const std::out_of_range &oor)
    {
      std::cerr << "Critical Error: Invalid sequence name " << option.name << "\n(" << oor.what() << ")\n";
    }
    (this->*step_fct)(option.separate);
  }
}

