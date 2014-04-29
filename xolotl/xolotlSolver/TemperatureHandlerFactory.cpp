#include "TemperatureHandlerFactory.h"
#include "TemperatureHandler.h"
#include "TemperatureProfileHandler.h"
#include <XolotlOptions.h>
#include <fstream>
#include <iostream>

namespace xolotlSolver
{

static std::shared_ptr<ITemperatureHandler> theTemperatureHandler;

// Create the desired type of handler registry.
bool initializeTemperature( bool useConstTempRegistry, bool useTempProfileRegistry,
		xolotlCore::XolotlOptions &options)
{
    bool ret = true;

    if ( useConstTempRegistry && useTempProfileRegistry )
    {
        // A constant temperature value AND a temperature profile cannot both be given.
    	std::cerr << "\nA constant temperature value AND a temperature file cannot both be given.  "
    			"Aborting!" << std::endl;
    	exit( EXIT_FAILURE );
    }
    else if( useConstTempRegistry )
    {
    	auto temp = options.getConstTemperature();
        // we are to use a constant temperature handler
    	std::cout << "\nHandler Temp = " << temp << std::endl;
        theTemperatureHandler = std::make_shared<TemperatureHandler>( temp );
    }
    else if( useTempProfileRegistry )
    {
    	auto tempFileName = options.getTempProfileFilename();
    	std::cout << "\nHandler Temperature file = " << tempFileName << std::endl;
        theTemperatureHandler = std::make_shared<TemperatureProfileHandler>( tempFileName );
        auto theTempHand = std::dynamic_pointer_cast<TemperatureProfileHandler>(theTemperatureHandler);
        theTempHand->initializeTempData(tempFileName.c_str());
    }
    else
    {
    	std::cerr << "Warning: Temperature information has not been given.  Defaulting to constant"
    			" temperature = 1000K " << std::endl;
    	auto temp = options.getConstTemperature();
        // we are to use a constant temperature handler
        theTemperatureHandler = std::make_shared<TemperatureHandler>( temp );
    }

    return ret;
}

// Provide access to our handler registry.
std::shared_ptr<ITemperatureHandler> getTemperatureHandler( xolotlCore::XolotlOptions &options )
{
    if( !theTemperatureHandler )
    {
        // We have not yet been initialized.
        std::cerr << "Warning: xolotlSolver temperature handler requested, but "
        		"library has not been initialized" << std::endl;

        xolotlSolver::initializeTemperature( false, false, options );
    }
    return theTemperatureHandler;
}


} // end namespace xolotlSolver


