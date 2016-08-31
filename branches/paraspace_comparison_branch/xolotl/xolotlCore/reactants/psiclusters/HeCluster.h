#ifndef HECLUSTER_H
#define HECLUSTER_H

// Includes
#include "PSICluster.h"
#include <xolotlPerf.h>

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
		PSICluster(1) {}

	/**
	 * This operation "combines" clusters in the sense that it handles all of
	 * the logic and caching required to correctly process the reaction
	 *
	 * (He_c)(V_b) + He_a --> [He_(a+c)][V_(b+1)] + I
	 *
	 * in the case of [He_(a+c)](V_b) not in the network
	 *
	 * This operation fills the reaction connectivity array as well as the
	 * array of combining clusters.
	 *
	 * @param clusters The clusters that can combine with this cluster
	 * (Here it will be HeV clusters and He clusters)
	 * @param productName The name of the product produced in the reaction
	 */
	void combineClusters(std::vector<Reactant *> & clusters,
			const std::string& productName);

public:

	/**
	 * The constructor. All HeClusters must be initialized with a size.
	 *
	 * @param nHe the number of helium atoms in the cluster
	 * @param registry The performance handler registry
	 */
	HeCluster(int nHe, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	~HeCluster() {}

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of HeCluster.
	 *
	 * @return A copy of this HeCluster
	 */
	virtual std::shared_ptr<Reactant> clone();

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux() const;

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getEmissionPartialDerivatives(std::vector<double> & partials) const;

protected:

	/**
	 * Computes a row of the reaction connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants alone can form a reaction, the element at the position
	 * of the second reactant is 1, otherwise 0.
	 */
	void createReactionConnectivity();

	/**
	 * Computes a row of the dissociation connectivity matrix
	 * corresponding to this cluster.
	 *
	 * Connections are made between this cluster and any clusters it affects
	 * in a dissociation reaction.
	 *
	 * The base-class implementation handles dissociation for regular clusters
	 * by processing the reaction.
	 *
	 */
	void createDissociationConnectivity();


}; //end class HeCluster

} /* namespace xolotlCore */
#endif
