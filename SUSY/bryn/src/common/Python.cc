#include <boost/python.hpp>
#include "WeeklyUpdatePlots.hh"
#include "Trigger.hh"
#include "TriggerTurnOns.hh"
#include "EventDump.hh"
#include "EventSelector.hh"
using namespace boost::python;

BOOST_PYTHON_MODULE(libbryn) {

  class_<Operation::WeeklyUpdatePlots, bases<Operation::_Base> >( "WeeklyUpdatePlots",
                    init<const Utils::ParameterSet&>());


  class_<Operation::eventselector_t, bases<Operation::_Base> >( "OP_RunLumiEvSelector",
                    init<const Utils::ParameterSet&>());

  class_<Operation::Trigger, bases<Operation::_Base> >( "PL_TriggerPlots",
                    init<const Utils::ParameterSet&>());

  class_<Operation::TriggerTurnOns, bases<Operation::_Base> >( "PL_TriggerTurnOns",
                    init<const Utils::ParameterSet&>());

  class_<Operation::eventDump, bases<Operation::_Base> >( "EventDump", init<>());

}
