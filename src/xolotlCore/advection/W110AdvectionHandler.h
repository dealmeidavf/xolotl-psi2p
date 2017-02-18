#ifndef W110ADVECTIONHANDLER_H
#define W110ADVECTIONHANDLER_H

// Includes
#include "SurfaceAdvectionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class W110AdvectionHandler: public SurfaceAdvectionHandler {

public:

	//! The Constructor
	W110AdvectionHandler() {}

	//! The Destructor
	~W110AdvectionHandler() {}

	/**
	 * This function initialize the list of clusters that will move through advection for a
	 * (110) tungsten material.
	 *
	 * @param network The network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the advecting clusters
	 */
	void initialize(IReactionNetwork *network, int *ofill) {
		// Get all the reactants and their number
		auto reactants = network->getAll();
		int networkSize = reactants->size();
		int dof = network->getDOF();

		// Clear the index and sink strength vectors
		indexVector.clear();
		sinkStrengthVector.clear();

		// Loop on all the reactants
		for (int i = 0; i < networkSize; i++) {
			// Get the i-th cluster
			auto cluster = (PSICluster *) reactants->at(i);
			// Get its diffusion coefficient
			double diffFactor = cluster->getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (xolotlCore::equal(diffFactor, 0.0)) continue;

			// Keep only the helium clusters
			if (cluster->getType() != heType) continue;

			// Get its size
			int heSize = cluster->getSize();

			// Switch on the size to get the sink strength (in eV.nm3)
			double sinkStrength = 0.0;
			switch (heSize) {
				case 1:
					sinkStrength = 0.92e-3;
					break;
				case 2:
					sinkStrength = 1.48e-3;
					break;
				case 3:
					sinkStrength = 6.73e-3;
					break;
				case 4:
					sinkStrength = 6.18e-3;
					break;
				case 5:
					sinkStrength = 33.61e-3;
					break;
				case 6:
					sinkStrength = 37.58e-3;
					break;
				case 7:
					sinkStrength = 41.90e-3;
					break;
		}

			// If the sink strength is still 0.0, this cluster is not advecting
			if (xolotlCore::equal(sinkStrength, 0.0)) continue;

			// Add its index (i) to the vector of indices
			indexVector.push_back(i);

			// Add the sink strength to the vector
			sinkStrengthVector.push_back(sinkStrength);

			// Set the off-diagonal part for the Jacobian to 1
			// Get its id
			int index = cluster->getId() - 1;
			// Set the ofill value to 1 for this cluster
			ofill[index * dof + index] = 1;
		}

		return;
	}

};
//end class W110AdvectionHandler

} /* end namespace xolotlCore */
#endif
