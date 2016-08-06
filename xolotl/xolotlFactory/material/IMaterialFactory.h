#ifndef IMATERIALHANDLERFACTORY_H
#define IMATERIALHANDLERFACTORY_H

#include <memory>
#include <Options.h>
#include <IFluxHandler.h>
#include <IAdvectionHandler.h>
#include <IDiffusionHandler.h>
#include <IDisplacementHandler.h>

namespace xolotlFactory {

/**
 * Realizations of this interface are responsible for handling the flux and the advection.
 * they are both dependent on the type of material under study.
 */
class IMaterialFactory {
public:

	/**
	 * The destructor
	 */
	~IMaterialFactory() {}

	/**
	 * Initialize the material conditions with the different given options.
	 *
	 * @param options The Xolotl options.
	 */
	virtual void initializeMaterial(xolotlCore::Options &options) = 0;

	/**
	 * Return the flux handler.
	 *
	 * @return The flux handler.
	 */
	virtual std::shared_ptr<xolotlCore::IFluxHandler> getFluxHandler() const = 0;

	/**
	 * Return the advection handler.
	 *
	 * @return The advection handler.
	 */
	virtual std::shared_ptr<xolotlCore::IAdvectionHandler> getAdvectionHandler() const = 0;

	/**
	 * Return the diffusion handler.
	 *
	 * @return The diffusion handler.
	 */
	virtual std::shared_ptr<xolotlCore::IDiffusionHandler> getDiffusionHandler() const = 0;

	/**
	 * Return the displacement handler.
	 *
	 * @return The displacement handler.
	 */
	virtual std::shared_ptr<xolotlCore::IDisplacementHandler> getDisplacementHandler() const = 0;

	/**
	 * Function that create the wanted material factory depending on the given type.
	 *
	 * @param materialType The type of wanted material.
	 * @param dimension The number of dimensions of the problem.
	 * @return The material factory.
	 */
	static std::shared_ptr<IMaterialFactory> createMaterialFactory(const std::string& materialType,
			int dimension);

};

} // end namespace xolotlFactory

#endif // IMATERIALHANDLERFACTORY_H
