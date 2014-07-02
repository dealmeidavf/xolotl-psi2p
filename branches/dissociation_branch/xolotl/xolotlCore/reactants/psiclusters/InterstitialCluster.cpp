#include "InterstitialCluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

InterstitialCluster::InterstitialCluster(int nI, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(nI, registry) {

	// Update the composition map
	compositionMap["I"] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "I_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "I";

	// Compute the reaction radius
	double EightPi = 8.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::latticeConstant, 3.0);
	double termOne = 1.15 * (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant;
	double termTwo = pow((3.0 / EightPi) * aCubed * size, (1.0 / 3.0));
	double termThree = pow((3.0 / EightPi) * aCubed, (1.0 / 3.0));
	reactionRadius = termOne + termTwo - termThree;

}

InterstitialCluster::~InterstitialCluster() {
}

std::shared_ptr<Reactant> InterstitialCluster::clone() {
	std::shared_ptr<Reactant> reactant(new InterstitialCluster(*this));
	return reactant;
}

void InterstitialCluster::createReactionConnectivity() {

	// Local Declarations - Note the reference to the properties map
	auto props = network->getProperties();
	int maxIClusterSize = std::stoi(props["maxIClusterSize"]);
	int maxHeIClusterSize = std::stoi(props["maxHeIClusterSize"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int numHeIClusters = std::stoi(props["numHeIClusters"]);
	int firstSize = 0, secondSize = 0, reactantsSize = 0;

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(getId());

	/*
	 * This section fills the array of reacting pairs that combine to produce
	 * this cluster from I clusters that are smaller than this.size. Each
	 * cluster i combines with a second cluster of size this.size - i.size.
	 *
	 * Total size starts with a value of one so that clusters of size one are
	 * not considered in this loop.
	 */
	for (firstSize = 1; firstSize <= (int) size/2; firstSize++) {
		secondSize = size - firstSize;
		// Get the first and second reactants for the reaction
		// first + second = this.
		auto firstReactant = (PSICluster *) network->get("I", firstSize);
		auto secondReactant = (PSICluster *) network->get("I", secondSize);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant) {
			ClusterPair pair(firstReactant,secondReactant);
			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	/* ----- I_a + I_b --> I_(a+b) -----
	 *	Interstitial absorption
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	auto reactants = network->getAll("I");
	combineClusters(reactants,maxIClusterSize,"I");

	/* ----- He_a + I_b --> (He_a)*(I_b)
	 * Interstitials can interact with clusters of He to form HeI clusters.
	 * They cannot cluster with He clusters that are so large that the
	 * combination of the two would produce an HeI cluster above the
	 * maximum size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	reactants = network->getAll("He");
	combineClusters(reactants,maxHeIClusterSize,"HeI");

	/* ----- I_a + V_b -----
	 * --> I_(a-b), if a > b
	 * --> V_(b-a), if a < b
	 * --> 0, if a = b
	 * Interstitial-Vacancy Annihilation
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	reactants = network->getAll("V");
	fillVWithI("V",reactants);
	// Mark the reaction connectivity for the cases where this cluster is
	// produced by the above reaction. This has to be checked for every
	// vacancy.
	reactantsSize = reactants.size();
	for (int i = 0; i < reactantsSize; i++) {
		auto firstReactant = (PSICluster *) reactants[i];
		// Get the interstitial cluster that is bigger than the vacancy
		// and can form this cluster. I only results when it is bigger than V.
		auto secondReactant = (PSICluster *) network->get("I",firstReactant->getSize() + size);
		// Update the connectivity
		if (secondReactant) {
			setReactionConnectivity(firstReactant->getId());
			setReactionConnectivity(secondReactant->getId());
		}
	}

	/* ----- (He_a)(V_b) + (I_c) --> (He_a)[V_(b-c)] -----
	 * Interstitials interact with all mixed-species clusters by
	 * annihilating vacancies.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (numHeVClusters > 0) {
		reactants = network->getAll("HeV");
		replaceInCompound(reactants,"V","I");
	}

	/* ----- (He_a)*(I_b) + I --> (He_a)*[I_(b + 1)] -----
	 * Single interstitial absorption by a HeI cluster under the condition
	 * that (x + y + 1) <= maxSize
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (size == 1 && numHeIClusters > 0) {
		reactants = network->getAll("HeI");
		combineClusters(reactants,maxHeIClusterSize,"HeI");
	}

	return;
}

void InterstitialCluster::createDissociationConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species dissociation
	PSICluster::createDissociationConnectivity();

	// Specific case for the single species cluster
	if (size == 1) {
		// I dissociation of HeI cluster is handled here
		// (He_a)(I_b) --> (He_a)[I_(b-1)] + I_1
		// Get all the HeI clusters of the network
		auto allHeIReactants = network->getAll("HeI");
		for (int i = 0; i < allHeIReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeIReactants.at(i);

			// (He_a)(I_b) is the dissociating one, (He_a)[I_(b-1)] is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			comp[iType] -= 1;
			std::vector<int> compositionVec = { comp[heType], comp[vType],
					comp[iType] };
			auto smallerReactant = network->getCompound("HeI", compositionVec);
			dissociateCluster(allHeIReactants.at(i), smallerReactant);
		}

		// Trap mutation of HeV cluster is handled here
		// (He_a)(V_b) --> He_(a)[V_(b+1)] + I_1
		// Get all the HeV clusters of the network
		auto allHeVReactants = network->getAll("HeV");
		for (int i = 0; i < allHeVReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeVReactants.at(i);

			// (He_a)(V_b) is the dissociating one, (He_a)[V_(b+1)] is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			comp[vType] += 1;
			std::vector<int> compositionVec = { comp[heType], comp[vType],
					comp[iType] };
			auto biggerReactant = network->getCompound("HeV", compositionVec);
			dissociateCluster(allHeVReactants.at(i), biggerReactant);
		}
	}

	return;
}
