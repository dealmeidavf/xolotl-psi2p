#include "xolotlPerf/papi/PAPIHandlerRegistry.h"
#include "xolotlPerf/papi/PAPITimer.h"
#include "xolotlPerf/standard/EventCounter.h"
#include "xolotlPerf/papi/PAPIHardwareCounter.h"

namespace xolotlPerf
{

std::shared_ptr<ITimer>
PAPIHandlerRegistry::getTimer(std::string name)
{
    // TODO - associate the object we create with the current region.
    std::shared_ptr<ITimer> ret;

    // check if we have already created a timer with this name
    auto iter = allTimers.find(name);
    if( iter != allTimers.end() )
    {
        // We have already created a timer with this name.
        // Return it.
        ret = iter->second;
    }
    else
    {
        // We have not yet created a timer with this name.
        // Build one, and keep track of it.
        ret = std::make_shared<PAPITimer>(name);
        allTimers[name] = ret;
    }
    return ret;
}


std::shared_ptr<IHardwareCounter>
PAPIHandlerRegistry::getHardwareCounter( std::string name,
                        const IHardwareCounter::SpecType& ctrSpec )
{
    // TODO - associate the object we create with the current region
    std::shared_ptr<IHardwareCounter> ret;

    // Check if we have already created a dummy hardware counter set with this name.
    auto iter = allHWCounterSets.find(name);
    if( iter != allHWCounterSets.end() )
    {
        // We have already created a hw counter set with this name.
        // Return it.
        ret = iter->second;
    }
    else
    {
        // We have not yet created a hw counter set with this name.
        // Build one and keep track of it.
        ret = std::make_shared<PAPIHardwareCounter>(name, ctrSpec);
        allHWCounterSets[name] = ret;
    }
    return ret;
}

} // namespace xolotlPerf

