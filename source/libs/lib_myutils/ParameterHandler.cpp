

#include "ParameterHandler.h"
#include "strtk.hpp"
#include "fftw3.h"
#include <cmath>
#include <iostream>
#include <algorithm>

extern double sign( double );
extern double Heaviside( double );
extern double rect( double, double, double );

ParameterHandler::ParameterHandler( const std::string filename )
{
  //Load xml file and get the first node
  if ( !m_xml_doc.load_file(filename.c_str()) )
  {
    std::cerr << "Critical error occured during loading of file " << filename << std::endl;
    throw;
  }

  m_map_freq.insert(std::pair<std::string,int>("none",freq::none));
  m_map_freq.insert(std::pair<std::string,int>("each",freq::each));
  m_map_freq.insert(std::pair<std::string,int>("last",freq::last));
  m_map_freq.insert(std::pair<std::string,int>("packed",freq::packed));

  //Read values from xml
  populate_constants();
  populate_vconstants();
  populate_algorithm();
  populate_simulation();
  populate_sequence();
  populate_analyze();
}

ParameterHandler::~ParameterHandler()
{
  m_map_constants.clear();
  m_map_vconstants.clear();
  m_map_algorithm.clear();
  m_map_simulation.clear();
  m_sequence.clear();
  m_analyze.clear();
}

void ParameterHandler::populate_sequence()
{
  m_sequence.clear();

  std::vector<std::string> vec;
  std::string tmpstr;

  std::string querystr = "/SIMULATION//SEQUENCE//*";
  pugi::xpath_node_set tools = m_xml_doc.select_nodes(querystr.c_str());

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;

    sequence_item item;

    item.name = node.node().name();
    item.content = node.node().child_value();

    item.dt = node.node().attribute("dt").as_double(0.001);
    item.Nk =  node.node().attribute("Nk").as_int(100);;
    item.comp = node.node().attribute("comp").as_int(0);

    tmpstr = node.node().attribute("output_freq").as_string("none");
    item.output_freq = m_map_freq[tmpstr];
    tmpstr = node.node().attribute("rabi_output_freq").as_string("none");
    item.rabi_output_freq = m_map_freq[tmpstr];
    tmpstr = node.node().attribute("pn_freq").as_string("last");
    item.compute_pn_freq = m_map_freq[tmpstr];
    tmpstr = node.node().attribute("custom_freq").as_string("none");
    item.custom_freq = m_map_freq[tmpstr];
    tmpstr = node.node().attribute("analyze").as_string("none");
    item.analyze = m_map_freq[tmpstr];

    item.chirp_mode = node.node().attribute("chirp_mode").as_int(0);
    item.no_of_chirps = node.node().attribute("no_of_chirps").as_int(1);
    item.chirp_min = node.node().attribute("chirp_min").as_double(0);
    item.chirp_max = node.node().attribute("chirp_max").as_double(0);

    vec.clear();
    strtk::parse(item.content,",",vec);

    for ( auto i : vec )
    {
      double val;
      try
      {
        val = stod(i);
      }
      catch ( const std::invalid_argument &ia )
      {
        std::cerr << "Error Parsing xml file: Unable to convert " << i << " to double for element Duration in section Sequence\n";
      }
      item.duration.push_back(val);
    }
    m_sequence.push_back(item);
  }
}

void ParameterHandler::populate_analyze()
{
  m_analyze.clear();

  std::string tmp;

  std::string querystr = "/SIMULATION//ANALYZE//*";
  pugi::xpath_node_set tools = m_xml_doc.select_nodes(querystr.c_str());

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;

    analyze_item item;

    item.name = node.node().name();
    tmp = node.node().child_value();

    std::transform(tmp.begin(), tmp.end(),tmp.begin(), ::toupper);

    if ( tmp == "TRUE" ) item.content = true;
    else item.content = false;
    item.separate = node.node().attribute("separate").as_bool(false);
    m_analyze.push_back(item);
  }
}

void ParameterHandler::populate_constants()
{
  m_map_constants.clear();

  double val;
  std::string tmp, str;

  std::string querystr = "/SIMULATION//CONSTANTS//*";
  pugi::xpath_node_set tools = m_xml_doc.select_nodes(querystr.c_str());

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;

    str = node.node().name();
    tmp = node.node().child_value();

    try
    {
      val = stod(tmp);
    }
    catch ( const std::invalid_argument &ia )
    {
      std::cerr << "Error Parsing xml file: Unable to convert " << tmp << " to double for element <" << str << "> in section CONSTANTS\n";
    }
    m_map_constants.insert ( std::pair<std::string,double>(str,val) );
  }
}

// See populate_constants for more details
void ParameterHandler::populate_vconstants()
{
  m_map_vconstants.clear();

  std::string tmp, str;
  std::vector<std::string> vec;
  std::vector<double> data;

  std::string querystr = "/SIMULATION//VCONSTANTS//*";
  pugi::xpath_node_set tools = m_xml_doc.select_nodes(querystr.c_str());

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;

    vec.clear();
    data.clear();
    str = node.node().name();
    tmp = node.node().child_value();
    strtk::parse(tmp,",",vec);

    for ( auto i : vec )
    {
      try
      {
        double val = stod(i);
        data.push_back(val);
      }
      catch ( const std::invalid_argument &ia )
      {
        std::cerr << "Error Parsing xml file: Unable to convert " << tmp << " to double for element <" << str << "> in section VCONSTANTS\n";
        throw;
      }
    }
    m_map_vconstants.insert ( std::pair<std::string,std::vector<double>>(str,data) );
  }

  /*
  for( auto it1 : m_map_vconstants )
  {
    std::cout << it1.first << "\t";
    for( auto it2 : it1.second )
    {
      std::cout << it2 << "\t";
    }
    std::cout << std::endl;
  }
  */
}

// See populate_constants
void ParameterHandler::populate_algorithm()
{
  m_map_algorithm.clear();

  std::string tmp, str;

  std::string querystr = "/SIMULATION//ALGORITHM//*";
  pugi::xpath_node_set tools = m_xml_doc.select_nodes(querystr.c_str());

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;

    str = node.node().name();
    tmp = node.node().child_value();

    m_map_algorithm.insert ( std::pair<std::string,std::string>(str,tmp) );
  }
}

// See populate_constants
void ParameterHandler::populate_simulation()
{
  m_map_simulation.clear();

  pugi::xpath_node_set tools = m_xml_doc.select_nodes("/SIMULATION/*/text()");

  for (pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    pugi::xpath_node node = *it;
    m_map_simulation[node.parent().name()] = node.parent().child_value();
  }
}

void ParameterHandler::Setup_muParser( mu::Parser &mup )
{
  mup.ClearConst();
  //mup.ClearFun();

  for ( auto i : m_map_constants )
  {
    //std::cout << i.first << ", " << i.second << std::endl;
    mup.DefineConst( i.first.c_str(), i.second );
  }

  mup.DefineFun("Heaviside", Heaviside, false);
  mup.DefineFun("rect", rect, false);
  mup.DefineFun("sign", sign, false);
}

std::string ParameterHandler::Get_simulation( const std::string k )
{
  auto it = m_map_simulation.find(k);
  if ( it == m_map_simulation.end() ) throw std::string( "Error: Could not find the key: " + k + " in section SIMULATION." );
  return (*it).second;
}

double ParameterHandler::Get_Constant( const std::string k )
{
  auto it = m_map_constants.find(k);
  if ( it == m_map_constants.end() ) throw std::string( "Error: Could not find the key: " + k + " in section CONSTANT." );
  return (*it).second;
}

double ParameterHandler::Get_VConstant( const std::string k, const int p )
{
  auto it = m_map_vconstants.find(k);
  if ( it == m_map_vconstants.end() ) throw std::string( "Error: Could not find the key " + k + " in section VCONSTANTS." );
  return (*it).second[p];
}

std::vector<double> ParameterHandler::Get_VConstant( const std::string k)
{
  std::vector<double> retval;
  auto it = m_map_vconstants.find(k);
  if ( it == m_map_vconstants.end() ) throw std::string( "Error: Could not find the key " + k + " in section VCONSTANTS." );
  retval = (*it).second;
  return retval;
}

bool ParameterHandler::Get_MomentumStates( std::vector<double> &vect, const int p)
{
  vect.clear();
  bool found = false;
  std::string tmp = "momentum_state_" + std::to_string(p);
  auto it = m_map_vconstants.find(tmp);

  if ( it != m_map_vconstants.end() )
  {
    vect = (*it).second;
    found = true;
  }
  return found;
}

double ParameterHandler::Get_dt()
{
  double retval=0.001;
  auto it = m_map_algorithm.find("DT");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

double ParameterHandler::Get_t()
{
  double retval=0;
  auto it = m_map_algorithm.find("T");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

double ParameterHandler::Get_epsilon()
{
  double retval=1e-5;
  auto it = m_map_algorithm.find("EPSILON");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

double ParameterHandler::Get_stepsize()
{
  double retval=0.001;
  auto it = m_map_algorithm.find("STEPSIZE");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

double ParameterHandler::Get_xMax()
{
  double retval=10.0;
  auto it = m_map_algorithm.find("XMAX");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

double ParameterHandler::Get_xMin()
{
  double retval=-10.0;
  auto it = m_map_algorithm.find("XMIN");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

double ParameterHandler::Get_yMax()
{
  double retval=10.0;
  auto it = m_map_algorithm.find("YMAX");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

double ParameterHandler::Get_yMin()
{
  double retval=-10.0;
  auto it = m_map_algorithm.find("YMIN");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

double ParameterHandler::Get_zMax()
{
  double retval=10.0;
  auto it = m_map_algorithm.find("ZMAX");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

double ParameterHandler::Get_zMin()
{
  double retval=-10.0;
  auto it = m_map_algorithm.find("ZMIN");
  if ( it != m_map_algorithm.end() ) retval = stod((*it).second);
  return retval;
}

int ParameterHandler::Get_NX()
{
  int retval=256;
  auto it = m_map_algorithm.find("NX");
  if ( it != m_map_algorithm.end() ) retval = stoi((*it).second);
  return retval;
}

int ParameterHandler::Get_NY()
{
  int retval=256;
  auto it = m_map_algorithm.find("NY");
  if ( it != m_map_algorithm.end() ) retval = stoi((*it).second);
  return retval;
}

int ParameterHandler::Get_NZ()
{
  int retval=256;
  auto it = m_map_algorithm.find("NZ");
  if ( it != m_map_algorithm.end() ) retval = stoi((*it).second);
  return retval;
}


int ParameterHandler::Get_NA()
{
  int retval=10;
  auto it = m_map_algorithm.find("NA");
  if ( it != m_map_algorithm.end() ) retval = stoi((*it).second);
  return retval;
}

int ParameterHandler::Get_NK()
{
  int retval=10;
  auto it = m_map_algorithm.find("NK");
  if ( it != m_map_algorithm.end() ) retval = stoi((*it).second);
  return retval;
}

int ParameterHandler::Get_MaxIter()
{
  int retval=10;
  auto it = m_map_algorithm.find("MAXITER");
  if ( it != m_map_algorithm.end() ) retval = stoi((*it).second);
  return retval;
}

void ParameterHandler::Get_Header( generic_header &header, bool bcomplex )
{
  header = {};
  std::string str;
  long long nDims;

  try
  {
    str = Get_simulation("DIM");
  }
  catch (std::string str)
  {
    std::cout << str << std::endl;
  }

  nDims = std::stod(str);

  header.nself = sizeof(generic_header);
  if ( bcomplex )
  {
    header.nDatatyp = sizeof(fftw_complex);
    header.bComplex = 1;
  }
  else
  {
    header.nDatatyp = sizeof(double);
    header.bComplex = 0;
  }
  header.bAtom = 1;
  header.nDims    = nDims;

  switch ( nDims )
  {
  case 1:
    header.nDimX = Get_NX();
    header.nDimY = 1;
    header.nDimZ = 1;
    header.xMax = Get_xMax();
    header.xMin = Get_xMin();
    header.dx = fabs(header.xMax-header.xMin)/double(header.nDimX);
    header.dkx  = 2*M_PI/fabs(header.xMax-header.xMin);
    break;
  case 2:
    header.nDimX = Get_NX();
    header.nDimY = Get_NY();
    header.nDimZ = 1;
    header.xMax = Get_xMax();
    header.xMin = Get_xMin();
    header.yMax = Get_yMax();
    header.yMin = Get_yMin();
    header.dx = fabs(header.xMax-header.xMin)/double(header.nDimX);
    header.dy = fabs(header.yMax-header.yMin)/double(header.nDimY);
    header.dkx  = 2*M_PI/fabs(header.xMax-header.xMin);
    header.dky  = 2*M_PI/fabs(header.yMax-header.yMin);
    break;
  case 3:
    header.nDimX = Get_NX();
    header.nDimY = Get_NY();
    header.nDimZ = Get_NZ();
    header.xMax = Get_xMax();
    header.xMin = Get_xMin();
    header.yMax = Get_yMax();
    header.yMin = Get_yMin();
    header.zMax = Get_zMax();
    header.zMin = Get_zMin();
    header.dx   = fabs(header.xMax-header.xMin)/double(header.nDimX);
    header.dy   = fabs(header.yMax-header.yMin)/double(header.nDimY);
    header.dz   = fabs(header.zMax-header.zMin)/double(header.nDimZ);
    header.dkx  = 2*M_PI/fabs(header.xMax-header.xMin);
    header.dky  = 2*M_PI/fabs(header.yMax-header.yMin);
    header.dkz  = 2*M_PI/fabs(header.zMax-header.zMin);
    break;
  }
  header.dt   = 0.001;
  header.nself_and_data = header.nself + (header.nDimX + header.nDimY)*header.nDatatyp;
}
