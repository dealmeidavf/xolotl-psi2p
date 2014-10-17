#include "TemperatureHandlerFactory.h"
#include "TemperatureHandler.h"
#include "TemperatureProfileHandler.h"
#include <fstream>
#include <iostream>
#include <mpi.h>

namespace xolotlCore {

static std::shared_ptr<ITemperatureHandler> theTemperatureHandler;

// Create the desired type of handler registry.
bool initializeTempHandler(Options &options) {
	// Get the current process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	bool ret = true;

	if (options.useConstTemperatureHandlers()
			&& options.useTemperatureProfileHandlers()) {
		// Only print the error message once when running in parallel
		if (procId == 0) {
			// A constant temperature value AND a temperature profile cannot both be given.
			throw std::string(
					"\nA constant temperature value AND a temperature file cannot both be given.");
		}
	} else if (options.useConstTemperatureHandlers()) {
		auto temp = options.getConstTemperature();
		// we are to use a constant temperature handler
		theTemperatureHandler = std::make_shared<TemperatureHandler>(temp);
	} else if (options.useTemperatureProfileHandlers()) {
		auto tempFileName = options.getTempProfileFilename();
//		std::cout << "\nHandler Temperature file = " << tempFileName << std::endl;
		theTemperatureHandler = std::make_shared<TemperatureProfileHandler>(
				tempFileName);
		theTemperatureHandler->initializeTemperature();
	} else {
		// Only print the error message once when running in parallel
		if (procId == 0) {
			std::cerr
					<< "\nWarning: Temperature information has not been given.  Defaulting to constant"
							" temperature = 1000K " << std::endl;
		}
		auto temp = options.getConstTemperature();
		// we are to use a constant temperature handler
		theTemperatureHandler = std::make_shared<TemperatureHandler>(temp);
	}

	return ret;
}

// Provide access to our handler registry.
std::shared_ptr<ITemperatureHandler> getTemperatureHandler() {
	if (!theTemperatureHandler) {
		// We have not yet been initialized.
		throw std::string("\nxolotlSolver temperature handler requested, but "
				"library has not been initialized.");
	}
	return theTemperatureHandler;
}

} // end namespace xolotlCore
