#ifndef IOPTIONS_H
#define IOPTIONS_H

// Includes
#include <iostream>
#include <string>
#include <map>
#include "xolotlPerf/xolotlPerf.h"


namespace xolotlCore {

/**
 * IOptions describes the structure needed for the options in xolotl.
 * All private members will be accessed through getters and setters.
 */
class IOptions {

public:

	/**
	 * The destructor
	 */
    virtual ~IOptions() {}

    /**
     * Read the parameters from the given file to set the different
     * xolotl options.
     *
     * @param argc The number of arguments in the argv vector.
     * @param argv Vector of argument strings.
     */
    virtual void readParams(int argc, char* argv[]) = 0;

    /**
     * Show our help message.
     * @param os The output stream upon which to print the help message.
     */
    virtual void showHelp(std::ostream& os) const = 0;

    /**
     * Should the program run after parsing the command line?
     * (Or did parsing the command line suggest the program
     * shouldn't/couldn't run?)
     * @return true is the program should run.
     */
    virtual bool shouldRun() const = 0;

    /**
     * Set the shouldRunFlag.
     * @param flag The value for the shouldRunFlag.
     */
    virtual void setShouldRunFlag(bool flag) = 0;

    /**
     * If program shouldn't run, what should its exit code be?
     * @return the value of the exit code.
     */
    virtual int getExitCode() const = 0;

    /**
     * Set the value for the exit code
     * @param code The value for exit code.
     */
    virtual void setExitCode(int code) = 0;

    /**
     * Get the name of the network file.
     * @return the name of the network file.
     */
    virtual std::string getNetworkFilename() const = 0;

    /**
     * Set the name of the network file.
     * @param name Name for the network file.
     */
    virtual void setNetworkFilename(std::string name) = 0;

    /**
     * Get the Argc for PETSc
     * @return argc.
     */
    virtual int getPetscArgc() const = 0;

    /**
     * Set the Argc for PETSc
     * @param argc The number of options for PETSc.
     */
    virtual void setPetscArgc(int argc) = 0;

    /**
     * Get the Argv for PETSc
     * @return argv.
     */
    virtual char** getPetscArgv() const = 0;

    /**
     * Set the Argv for PETSc
     * @param argv The pointer to the options for PETSc.
     */
    virtual void setPetscArgv(char** argv) = 0;

    /**
     * Should we use const temperature handlers?
     * @return true if xolotl must use a constant temperature.
     */
    virtual bool useConstTemperatureHandlers() const = 0;

    /**
     * Set the constTempFlag.
     * @param flag The value for the constTempFlag.
     */
    virtual void setConstTempFlag(bool flag) = 0;

    /**
     * Obtain the value of the constant temperature to be used.
     * @return The value for the temperature.
     */
    virtual double getConstTemperature() const = 0;

    /**
     * Set the constant temperature.
     * @param temp The value for the constant temperature.
     */
    virtual void setConstTemperature(double temp) = 0;

    /**
     * Should we use temperature profile handlers?
     * @return true if xolotl must use a temperature profile.
     */
    virtual bool useTemperatureProfileHandlers() const = 0;

    /**
     * Set the tempProfileFlag.
     * @param flag The value for the tempProfileFlag.
     */
    virtual void setTempProfileFlag(bool flag) = 0;

    /**
     * Obtain the name of the file containing the temperature profile data.
     * @return The name of the file.
     */
    virtual std::string getTempProfileFilename() const = 0;

    /**
     * Set the name of the profile file to use.
     * @param name The name of the file.
     */
    virtual void setTempProfileFilename(std::string name) = 0;

    /**
     * Should we use the Helium fluence option?
	 * If false, it will not be used.
     * @return true if the Helium fluence option was present in the parameter file,
	 * false if it was not.
     */
    virtual bool useMaxHeliumFluence() const = 0;

    /**
     * Set the heliumFluenceFlag.
     * @param flag The value for the heliumFluenceFlag.
     */
    virtual void setHeliumFluenceFlag(bool flag) = 0;

    /**
     * Obtain the value of the Helium fluence to be used.
     * @return The value of the maximum fluence.
     */
    virtual double getMaxHeliumFluence() const = 0;

    /**
     * Set the value for the maximum fluence to which we want to integrate.
     * @param fluence The value for the maximum fluence.
     */
    virtual void setMaxHeliumFluence(double fluence) = 0;

    /**
     * Should we use the Helium flux option?
	 * If false, it will not be used.
     * @return true if the Helium flux option was present in the parameter file,
	 * false if it was not.
     */
    virtual bool useHeliumFlux() const = 0;

    /**
     * Set the heliumFluxFlag.
     * @param flag The value for the heliumFluenceFlag.
     */
    virtual void setHeliumFluxFlag(bool flag) = 0;

    /**
     * Obtain the value of the Helium flux to be used.
     * @return The value of the flux.
     */
    virtual double getHeliumFlux() const = 0;

    /**
     * Set the value for the flux to use.
     * @param flux The value for the flux.
     */
    virtual void setHeliumFlux(double flux) = 0;

    /**
     * Which type of performance handlers should we use?
     * @return The type of performance handler registry to use.
     */
    virtual xolotlPerf::IHandlerRegistry::RegistryType getPerfHandlerType(void) const = 0;

    /**
     * Set the type of performance handlers to use.
     * @param rtype The type of performance handler registry to use.
     */
    virtual void setPerfHandlerType(xolotlPerf::IHandlerRegistry::RegistryType rtype) = 0;


    /**
     * Should we use the "standard" set of handlers for the visualization?
     * If false, use dummy (stub) handlers.
     * @return true if program should use standard handlers, false if
     * should use dummy handlers.
     */
    virtual bool useVizStandardHandlers() const = 0;

    /**
     * Set the vizStandardHandlersFlag.
     * @param flag The value for the vizStandardHandlersFlag.
     */
    virtual void setVizStandardHandlers(bool flag) = 0;

};//end class IOptions

} /* namespace xolotlCore */

#endif
