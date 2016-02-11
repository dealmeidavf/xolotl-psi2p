#ifndef IOPTIONS_H
#define IOPTIONS_H

// Includes
#include <iostream>
#include <string>
#include <map>
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * IOptions describes the structure needed for the options in Xolotl.
 * All private members will be accessed through getters and setters.
 */
class IOptions {

public:

	/**
	 * The destructor
	 */
	virtual ~IOptions() {
	}

	/**
	 * Read the parameters from the given file to set the different
	 * Xolotl options.
	 *
	 * @param argc The number of arguments in the argv vector
	 * @param argv Vector of argument strings
	 */
	virtual void readParams(int argc, char* argv[]) = 0;

	/**
	 * Show our help message.
	 *
	 * @param os The output stream upon which to print the help message
	 */
	virtual void showHelp(std::ostream& os) const = 0;

	/**
	 * Should the program run after parsing the parameter file?
	 *
	 * @return true is the program should run
	 */
	virtual bool shouldRun() const = 0;

	/**
	 * Set the shouldRunFlag.
	 *
	 * @param flag The value for the shouldRunFlag
	 */
	virtual void setShouldRunFlag(bool flag) = 0;

	/**
	 * If program shouldn't run, what should its exit code be?
	 *
	 * @return the value of the exit code
	 */
	virtual int getExitCode() const = 0;

	/**
	 * Set the value for the exit code.
	 *
	 * @param code The value for exit code
	 */
	virtual void setExitCode(int code) = 0;

	/**
	 * Get the name of the network file.
	 *
	 * @return the name of the network file
	 */
	virtual std::string getNetworkFilename() const = 0;

	/**
	 * Set the name of the network file.
	 *
	 * @param name Name for the network file
	 */
	virtual void setNetworkFilename(const std::string& name) = 0;

	/**
	 * Get the Argc for PETSc.
	 *
	 * @return argc
	 */
	virtual int getPetscArgc() const = 0;

	/**
	 * Set the Argc for PETSc.
	 *
	 * @param argc The number of options for PETSc
	 */
	virtual void setPetscArgc(int argc) = 0;

	/**
	 * Get the Argv for PETSc.
	 *
	 * @return argv
	 */
	virtual char** getPetscArgv() const = 0;

	/**
	 * Set the Argv for PETSc.
	 *
	 * @param argv The pointer to the options for PETSc
	 */
	virtual void setPetscArgv(char** argv) = 0;

	/**
	 * Should we use const temperature handlers?
	 *
	 * @return true if Xolotl must use a constant temperature
	 */
	virtual bool useConstTemperatureHandlers() const = 0;

	/**
	 * Set the constTempFlag.
	 *
	 * @param flag The value for the constTempFlag
	 */
	virtual void setConstTempFlag(bool flag) = 0;

	/**
	 * Obtain the value of the constant temperature to be used.
	 *
	 * @return The value for the temperature
	 */
	virtual double getConstTemperature() const = 0;

	/**
	 * Set the constant temperature.
	 *
	 * @param temp The value for the constant temperature
	 */
	virtual void setConstTemperature(double temp) = 0;

	/**
	 * Should we use temperature profile handlers?
	 *
	 * @return true if Xolotl must use a temperature profile
	 */
	virtual bool useTemperatureProfileHandlers() const = 0;

	/**
	 * Set the tempProfileFlag.
	 *
	 * @param flag The value for the tempProfileFlag
	 */
	virtual void setTempProfileFlag(bool flag) = 0;

	/**
	 * Obtain the name of the file containing the temperature profile data.
	 *
	 * @return The name of the file
	 */
	virtual std::string getTempProfileFilename() const = 0;

	/**
	 * Set the name of the profile file to use.
	 *
	 * @param name The name of the file
	 */
	virtual void setTempProfileFilename(const std::string& name) = 0;

	/**
	 * Should we use the flux amplitude option?
	 * If false, it will not be used.
	 *
	 * @return true if the flux amplitude option was present in
	 * the parameter file, false if it was not
	 */
	virtual bool useFluxAmplitude() const = 0;

	/**
	 * Set the fluxFlag.
	 *
	 * @param flag The value for the fluxFlag
	 */
	virtual void setFluxFlag(bool flag) = 0;

	/**
	 * Obtain the value of the intensity of the flux to be used.
	 *
	 * @return The value of the flux amplitude
	 */
	virtual double getFluxAmplitude() const = 0;

	/**
	 * Set the value for the flux intensity to use.
	 *
	 * @param flux The value for the flux amplitude
	 */
	virtual void setFluxAmplitude(double flux) = 0;

	/**
	 * Should we use a time profile for the helium flux?
	 *
	 * @return True is a time profile file is given for the helium flux
	 */
	virtual bool useFluxTimeProfile() const = 0;

	/**
	 * Set the fluxProfileFlag.
	 *
	 * @param flag The value for the flag
	 */
	virtual void setFluxProfileFlag(bool flag) = 0;

	/**
	 * Obtain the name of the file containing the time profile data for the
	 * helium flux.
	 *
	 * @return The name of the file
	 */
	virtual std::string getFluxProfileName() const = 0;

	/**
	 * Set the name of the time profile file to use.
	 *
	 * @param name The name of the file
	 */
	virtual void setFluxProfileName(const std::string& name) = 0;

	/**
	 * Which type of performance handlers should we use?
	 *
	 * @return The type of performance handler registry to use
	 */
	virtual xolotlPerf::IHandlerRegistry::RegistryType getPerfHandlerType(
			void) const = 0;

	/**
	 * Set the type of performance handlers to use.
	 *
	 * @param rtype The type of performance handler registry to use
	 */
	virtual void setPerfHandlerType(
			xolotlPerf::IHandlerRegistry::RegistryType rtype) = 0;

	/**
	 * Should we use the "standard" set of handlers for the visualization?
	 * If false, use dummy (stub) handlers.
	 *
	 * @return true if program should use standard handlers, false if
	 * should use dummy handlers
	 */
	virtual bool useVizStandardHandlers() const = 0;

	/**
	 * Set the vizStandardHandlersFlag.
	 *
	 * @param flag The value for the vizStandardHandlersFlag
	 */
	virtual void setVizStandardHandlers(bool flag) = 0;

	/**
	 * Obtain the name of the material to be used for the simulation.
	 *
	 * @return The name of the material
	 */
	virtual std::string getMaterial() const = 0;

	/**
	 * Set the name of the material to be used for the simulation.
	 *
	 * @param material The name of the material
	 */
	virtual void setMaterial(const std::string& material) = 0;

	/**
	 * Obtain the value of the concentration for the vacancies.
	 *
	 * @return The concentration value
	 */
	virtual double getInitialVConcentration() const = 0;

	/**
	 * Set the value of the concentration for the vacancies.
	 *
	 * @param conc The value for the concentration
	 */
	virtual void setInitialVConcentration(double conc) = 0;

	/**
	 * Obtain the number of dimensions for the simulation.
	 *
	 * @return The number of dimensions
	 */
	virtual int getDimensionNumber() const = 0;

	/**
	 * Set the number of dimensions for the simulation.
	 *
	 * @param number The number of dimensions
	 */
	virtual void setDimensionNumber(int number) = 0;

	/**
	 * Obtain the vacancy size at which the grouping scheme starts.
	 *
	 * @return The size
	 */
	virtual int getGroupingVMin() const = 0;

	/**
	 * Set the vacancy size at which the grouping scheme starts.
	 *
	 * @param size The size
	 */
	virtual void setGroupingVMin(int size) = 0;

	/**
	 * Obtain the helium width for the grouping scheme.
	 *
	 * @return The width
	 */
	virtual int getGroupingHeWidth() const = 0;

	/**
	 * Set the helium width for the grouping scheme.
	 *
	 * @param width The width
	 */
	virtual void setGroupingHeWidth(int width) = 0;

	/**
	 * Obtain the vacancy width for the grouping scheme.
	 *
	 * @return The width
	 */
	virtual int getGroupingVWidth() const = 0;

	/**
	 * Set the vacancy width for the grouping scheme.
	 *
	 * @param width The width
	 */
	virtual void setGroupingVWidth(int width) = 0;

};
//end class IOptions

} /* namespace xolotlCore */

#endif
