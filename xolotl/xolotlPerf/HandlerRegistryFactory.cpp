#include "XolotlConfigPerf.h"
#include "HandlerRegistryFactory.h"
#include "dummy/DummyHandlerRegistry.h"

#if defined(HAVE_PERFLIB_STD)
#include "standard/StandardHandlerRegistry.h"
#endif // defined(HAVE_PERFLIB_STD)


namespace xolotlPerf
{

static std::shared_ptr<IHandlerRegistry> theHandlerRegistry;

// Create the desired type of handler registry.
bool
initialize( bool useStdRegistry )
{
    bool ret = true;

    if( useStdRegistry )
    {
#if defined(HAVE_PERFLIB_STD)
        // we are to use a standard handler registry
        // (one that collects timings)
        theHandlerRegistry = std::make_shared<StandardHandlerRegistry>();
#else
        // TODO is there another mechanism for writing errors
        // e.g., one that logs error messages?
        std::cerr << "xolotlPerf::initialize: unable to build requested standard handler registry" << std::endl;
        ret = false;
#endif // defined(HAVE_PERFLIB_STD)
    }
    else
    {
        // use a dummy HandlerRegistry for this run
        theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
    }

    return ret;
}




// Provide access to our handler registry.
std::shared_ptr<IHandlerRegistry>
getHandlerRegistry( void )
{
    return theHandlerRegistry;
}


} // end namedpsace xolotlPerf

