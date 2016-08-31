#ifndef FEMATERIALHANDLERFACTORY_H
#define FEMATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <FeFitFluxHandler.h>
#include <DummyAdvectionHandler.h>
#include <DummyTrapMutationHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>
#include <DummyBubbleBurstingHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for iron material.
 */
class FeMaterialFactory : public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	FeMaterialFactory() {}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	FeMaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::FeFitFluxHandler>();
		theAdvectionHandler.push_back(std::make_shared<xolotlCore::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<xolotlCore::DummyTrapMutationHandler>();
		theBubbleBurstingHandler = std::make_shared<xolotlCore::DummyBubbleBurstingHandler>();

		// Switch on the dimension for the diffusion handler
		switch (dim) {
			case 1:
				theDiffusionHandler = std::make_shared<xolotlCore::Diffusion1DHandler>();
				break;
			case 2:
				theDiffusionHandler = std::make_shared<xolotlCore::Diffusion2DHandler>();
				break;
			case 3:
				theDiffusionHandler = std::make_shared<xolotlCore::Diffusion3DHandler>();
				break;
			default:
				// The asked dimension is not good (e.g. -1, 0, 4)
				throw std::string("\nxolotlFactory: Bad dimension for the W100 material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~FeMaterialFactory() {}
};

} // end namespace xolotlFactory

#endif // FEMATERIALHANDLERFACTORY_H
