#ifndef IADVECTIONHANDLER_H
#define IADVECTIONHANDLER_H

// Includes
#include <PSICluster.h>
#include <PSIClusterReactionNetwork.h>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for all the physical parts
 * for the advection of mobile helium cluster. The solver call these methods to handle
 * the advection.
 */
class IAdvectionHandler {

public:

	/**
	 * The destructor
	 */
	virtual ~IAdvectionHandler() {}

	/**
	 * This function initialize the list of clusters that will move through advection
	 * and their corresponding sink strength (or driving forces).
	 *
	 * @param network The network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the advecting clusters
	 */
	virtual void initialize(PSIClusterReactionNetwork *network, int *ofill) = 0;

	/**
	 * Compute the flux due to the advection for all the helium clusters,
	 * given the space parameter hx and the position.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param hx The space parameter, here the grid step size in the x direction
	 * @param pos The position on the grid
	 * @param concVector The pointer to the pointer of arrays of concentration at middle,
	 * left, and right grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the advection is computed used to find the next solution
	 */
	virtual void computeAdvection(PSIClusterReactionNetwork *network, double hx,
			std::vector<double> &pos, double **concVector, double *updatedConcOffset) = 0;

	/**
	 * Compute the partial derivatives due to the advection of all the helium clusters given
	 * the space parameter hx and the position.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param hx The space parameter, here the grid step size in the x direction
	 * @param val The pointer to the array that will contain the values of partials
	 * for the advection
	 * @param indices The pointer to the array that will contain the indices of the
	 * advecting cluster in the network
	 * @param pos The position on the grid
	 */
	virtual void computePartialsForAdvection(PSIClusterReactionNetwork *network,
			double hx, double *val, int *indices, std::vector<double> &pos) = 0;

	/**
	 * Get the total number of advecting clusters in the network.
	 *
	 * @return The number of advecting clusters
	 */
	virtual int getNumberOfAdvecting() = 0;

};
//end class IAdvectionHandler

} /* namespace xolotlCore */
#endif
