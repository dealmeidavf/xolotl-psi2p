// Includes
#include "VCluster.h"
#include <Constants.h>

using namespace xolotlCore;

VCluster::VCluster(int nV) :
		PSICluster(nV) {
	// Set the reactant name appropriately
	name = "Vacancy";
}

VCluster::~VCluster() {
}

void VCluster::createReactionConnectivity() {

	// Local Declarations - Note the reference to the properties map
	std::map<std::string, std::string> props = *(network->properties);
	int numV = size, indexOther;
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int numHeIClusters = std::stoi(props["numHeIClusters"]);
	int numIClusters = std::stoi(props["numIClusters"]);
	std::map<std::string, int> speciesMap;
	int totalSize = 1, firstSize = 0, secondSize = 0;
	int firstIndex = -1, secondIndex = -1;
	std::map<std::string, int> firstSpeciesMap, secondSpeciesMap;
	std::shared_ptr<Reactant> firstReactant, secondReactant;
	std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>reactants = network->reactants;

	/*
	 * This section fills the array of reacting pairs that combine to produce
	 * this cluster. The only reactions that produce V clusters are those V
	 * clusters that are smaller than this.size. Each cluster i combines with
	 * a second cluster of size this.size - i.size.
	 *
	 * Total size starts with a value of one so that clusters of size one are
	 * not considered in this loop.
	 */
	while (totalSize < size) {
		// Increment the base sizes
		++firstSize;
		secondSize = size - firstSize;
		// Update the maps
		firstSpeciesMap["V"] = firstSize;
		secondSpeciesMap["V"] = secondSize;
		// Get the first and second reactants for the reaction
		// first + second = this.
		firstIndex = network->toClusterIndex(firstSpeciesMap);
		firstReactant = reactants->at(firstIndex);
		secondIndex = network->toClusterIndex(secondSpeciesMap);
		secondReactant = reactants->at(secondIndex);
		// Create a ReactingPair with the two reactants
		ReactingPair pair;
		pair.first = std::dynamic_pointer_cast<PSICluster>(firstReactant);
		pair.second = std::dynamic_pointer_cast<PSICluster>(secondReactant);
		// Add the pair to the list
		reactingPairs.push_back(pair);
		// Update the total size. Do not delete this or you'll have an infinite
		// loop!
		totalSize = firstSize + secondSize;
	}

	// Vacancies interact with everything except for vacancies bigger than they
	// would combine with to form vacancies larger than the size limit.

	/* -----  A*He + B*V → (A*He)(B*V) -----
	 * Vacancy clusters can interact with any helium cluster so long as the sum
	 * of the number of helium atoms and vacancies does not produce a cluster
	 * with a size greater than the maximum mixed-species cluster size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numHeOther = 1; numV + numHeOther <= maxMixedClusterSize;
			numHeOther++) {
		speciesMap["He"] = numHeOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(reactants->at(indexOther));
	}

	/* ----- A*V + B*V --> (A+B)*V -----
	 * This cluster should interact with all other clusters of the same type up
	 * to the max size minus the size of this one to produce larger clusters.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numVOther = 1; numV + numVOther <= maxVClusterSize; numVOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["V"] = numVOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(reactants->at(indexOther));
	}

	/* ----- A*I + B*V -----
	 * → (A-B)*I, if A > B
	 * → (B-I)*V, if A < B
	 * → 0, if A = B -----
	 * Vacancies always annihilate interstitials.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numIOther = 1; numIOther <= numIClusters; numIOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["I"] = numIOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(reactants->at(indexOther));
	}

	/* ----- (A*He)(B*V) + C*V → (A*He)[(B+C)*V] -----
	 * Vacancies can interact with a mixed-species cluster so long as the sum of
	 * the number of vacancy atoms and the size of the mixed-species cluster
	 * does not exceed the maximum mixed-species cluster size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (numV == 1) {
		for (int numVOther = 1; numVOther <= maxMixedClusterSize; numVOther++) {
			for (int numHeOther = 1;
					(numHeOther + numVOther + numV) <= maxMixedClusterSize;
					numHeOther++) {
				// Clear the map since we are reusing it
				speciesMap.clear();
				speciesMap["He"] = numHeOther;
				speciesMap["V"] = numVOther;
				int indexOther = network->toClusterIndex(speciesMap);
				reactionConnectivity[indexOther] = 1;
				combiningReactants.push_back(reactants->at(indexOther));
			}
		}
	}

	/* ----- (AHe)*(BI) + (CV) --> (AHe)*(B - C)V -----
	 * Vacancy absorption by HeI under the condition that y - z >= 1
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (numHeIClusters > 0) {
		for (int numIOther = 1; numIOther <= maxMixedClusterSize; numIOther++) {
			for (int numHeOther = 1;
					numIOther + numHeOther <= maxMixedClusterSize;
					numHeOther++) {
				// Clear the map since we are reusing it
				speciesMap.clear();
				bool connects = numIOther - numV >= 1;
				speciesMap["He"] = numHeOther;
				speciesMap["I"] = numIOther;
				int indexOther = network->toClusterIndex(speciesMap);
				reactionConnectivity[indexOther] = (int) connects;
				combiningReactants.push_back(reactants->at(indexOther));
			}
		}
	}

	return;
}

void VCluster::createDissociationConnectivity() {
	// Local Declarations
	int nReactants = network->reactants->size();
	std::map<std::string, int> clusterMap;

	// Vacancy Dissociation
	clusterMap["He"] = 0;
	clusterMap["V"] = size - 1;
	clusterMap["I"] = 0;
	if (size != 1) {
		dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
		clusterMap["V"] = 1;
		dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	}

	// Trap Mutation
	clusterMap["V"] = size + 1;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["I"] = 1;
	clusterMap["V"] = 0;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
}

bool VCluster::isProductReactant(int reactantI, int reactantJ) {

	// Local Declarations, integers for species number for I, J reactants
	int rI_I = 0, rJ_I = 0, rI_He = 0, rJ_He = 0, rI_V = 0, rJ_V = 0;

	// Get the ClusterMap corresponding to
	// the given reactants
	std::map<std::string, int> reactantIMap = network->toClusterMap(reactantI);
	std::map<std::string, int> reactantJMap = network->toClusterMap(reactantJ);

	// Grab the numbers for each species
	// from each Reactant
	rI_I = reactantIMap["I"];
	rJ_I = reactantJMap["I"];
	rI_He = reactantIMap["He"];
	rJ_He = reactantJMap["He"];
	rI_V = reactantIMap["V"];
	rJ_V = reactantJMap["V"];

	// We should have no interstitials, a
	// total of 0 Helium, and a total of
	// size Vacancies
	return ((rI_I + rJ_I) == 0) && ((rI_He + rJ_He) == 0)
			&& ((rI_V + rJ_V) == size);
}

std::map<std::string, int> VCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = 0;
	clusterMap["V"] = size;
	clusterMap["I"] = 0;

	// Return it
	return clusterMap;
}

double VCluster::getReactionRadius() {
	// FIXME Not right...
	return (sqrt(3) / 4) * xolotlCore::latticeConstant;
}
