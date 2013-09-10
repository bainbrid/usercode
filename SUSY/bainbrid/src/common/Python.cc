#include <boost/python.hpp>
#include "RobOps.hh"
#include "RobPlottingOps.hh"
#include "DeadRegionsOps.hh"
#include "ToyMC.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(libbainbrid) {
  
  class_<Operation::RobAlphaT, bases<Operation::_Base> >( "RobAlphaT",  init<float>() );
  
  class_<Operation::RobOps, bases<Operation::_Base> >( "RobOps",  
						       init<const Utils::ParameterSet&>() );
  
  class_<Operation::RobPlottingOps, bases<Operation::_Base> >( "RobPlottingOps",
							       init<const Utils::ParameterSet&>() );
  
  class_<Operation::DeadRegionsOps, bases<Operation::_Base> >( "DeadRegionsOps", 
							       init<const Utils::ParameterSet&>() );
  
  class_<Operation::ToyMC, bases<Operation::_Base> >( "ToyMC", init<const Utils::ParameterSet&>() );
  
}
