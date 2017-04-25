#include <string>
#include <iostream>
#include <fstream>

#include "my_structs.h"
#include "CPoint.h"
#include "cft_1d.h"
#include "cft_2d.h"
#include "cft_3d.h"
#include "shared_ana_tools.h"
#include "ParameterHandler.h"

int main(int argc, char *argv[])
{
  if ( argc > 3 )
  {
    //printf( "No parameter xml file specified.\n" );
    return EXIT_FAILURE;
  }

  //ParameterHandler object from xml
  ParameterHandler params(argv[1]);
  std::string filename;

  if ( argc == 3 )
  {
    filename = argv[2];
  }

  int dim=0;

  try
  {
    std::string tmp = params.Get_simulation("DIM");
    dim = std::stod(tmp);
  }
  catch ( std::string &str )
  {
    std::cout << str << std::endl;
  }

  if (dim == 1)
  {
    if ( argc == 2 )
    {
      Shared_Ana_Tools<1,Fourier::cft_1d> *tool = new Shared_Ana_Tools<1,Fourier::cft_1d>(&params);

      tool->Run_Analysis_in_Directory();
      tool->Write();
      delete tool;
    }
    if ( argc == 3 )
    {
      Shared_Ana_Tools<1,Fourier::cft_1d> *tool = new Shared_Ana_Tools<1,Fourier::cft_1d>(filename,&params);
      tool->Run_Analysis();
      tool->Write();
      delete tool;
    }

  }

  if (dim == 2)
  {
    if ( argc == 2 )
    {
      Shared_Ana_Tools<2,Fourier::cft_2d> *tool = new Shared_Ana_Tools<2,Fourier::cft_2d>(&params);

      tool->Run_Analysis_in_Directory();
      tool->Write();
      delete tool;
    }
    if ( argc == 3 )
    {
      Shared_Ana_Tools<2,Fourier::cft_2d> *tool = new Shared_Ana_Tools<2,Fourier::cft_2d>(filename,&params);
      tool->Run_Analysis();
      tool->Write();
      delete tool;
    }

  }

  if (dim == 3)
  {
    if ( argc == 2 )
    {
      Shared_Ana_Tools<3,Fourier::cft_3d> *tool = new Shared_Ana_Tools<3,Fourier::cft_3d>(&params);

      tool->Run_Analysis_in_Directory();
      tool->Write();
      delete tool;
    }
    if ( argc == 3 )
    {
      Shared_Ana_Tools<3,Fourier::cft_3d> *tool = new Shared_Ana_Tools<3,Fourier::cft_3d>(filename,&params);
      tool->Run_Analysis();
      tool->Write();
      delete tool;
    }

  }
}
