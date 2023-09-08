#include "cantera/thermo.h"
#include <iostream>

//static CommDom* ptr_class = NULL;

extern "C"
{
  void
  cantera_alya_create_( )
  {
    std::unique_ptr<Cantera::ThermoPhase> gas(Cantera::newPhase("h2o2.cti","ohmech"));
    gas->setState_TPX(500.0, 2.0*Cantera::OneAtm, "H2O:1.0, H2:8.0, AR:1.0");
    std::cout << gas->report() << std::endl;
  }

  void
  cantera_alya_delete_()
  {
  }


} // extern "C"
