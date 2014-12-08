#ifndef FLUXOPTIONHANDLER_H
#define FLUXOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * FluxOptionHandler handles the flux option.
 */
class FluxOptionHandler : public OptionHandler {
public:

	/**
	 * The default constructor
	 */
    FluxOptionHandler() :
    	OptionHandler("heFlux",
    			"heFlux <value>                    "
    			"This option allows the user to change the Helium flux "
    			"by the factor specified (in nm).\n") {}

	/**
	 * The destructor
	 */
    ~FluxOptionHandler() {}

    /**
     * This method will set the IOptions heliumFluxFlag and heliumFlux
     * to the value given as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The value for the helium flux.
     */
    bool handler(IOptions *opt, std::string arg) {
    	// Set the corresponding flag to true
    	opt->setHeliumFluxFlag(true);

    	// Set the value for the flux
    	double flux = strtod(arg.c_str(), NULL);

    	opt->setHeliumFlux(flux);
    	return true;
    }

};//end class FluxOptionHandler

} /* namespace xolotlCore */

#endif
