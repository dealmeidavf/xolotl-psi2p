#ifndef TEMPERATUREHANDLERFACTORY_H
#define TEMPERATUREHANDLERFACTORY_H

#include <memory>
#include <ITemperatureHandler.h>
#include <Options.h>

namespace xolotlFactory {

/**
 * Build the desired type of handler registry.
 *
 * @param options Options for the program
 * @return True iff the handler registry was created successfully.
 */
bool initializeTempHandler(xolotlCore::Options &options);

/**
 * Access the handler registry.
 *
 *  @return The handler registry object.
 */
std::shared_ptr<xolotlCore::ITemperatureHandler> getTemperatureHandler();

} // end namespace xolotlPerf

#endif /* TEMPERATUREHANDLERFACTORY_H */
