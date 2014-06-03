#ifndef HECLUSTER_H
#define HECLUSTER_H

// Includes
#include "PSICluster.h"
#include "../../xolotlPerf/xolotlPerf.h"

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of helium.
 */
class HeCluster: public PSICluster {

private:

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size and performance handler registry
	 */
	HeCluster() :
		PSICluster(1)
	{ }

public:

	/**
	 * The constructor. All HeClusters must be initialized with a size.
	 * @param nHe the number of helium atoms in the cluster
	 */
	HeCluster(int nHe, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	~HeCluster();

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of HeCluster.
	 * @return A copy of this HeCluster.
	 */
	virtual std::shared_ptr<Reactant> clone();

protected:

	/**
	 * Computes a row of the reaction connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants alone can form a reaction, the element at the position
	 * of the second reactant is 1, otherwise 0.
	 */
	void createReactionConnectivity();


}; //end class HeCluster

} /* namespace xolotlCore */
#endif
