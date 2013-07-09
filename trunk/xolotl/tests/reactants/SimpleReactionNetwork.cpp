/*
 * SimpleReactionNetwork.cpp
 *
 *  Created on: Jun 29, 2013
 *      Author: bkj
 */

#include "SimpleReactionNetwork.h"
#include <PSICluster.h>
#include <HeCluster.h>
#include <VCluster.h>
#include <InterstitialCluster.h>
#include <HeVCluster.h>
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace xolotlCore;
using namespace testUtils;

SimpleReactionNetwork::SimpleReactionNetwork() :
	ReactionNetwork() {
	
	// Hard code the size of the largest cluster
	int maxClusterSize = 10;
	int numClusters = maxClusterSize;

	// Add He clusters
	for (int numHe = 1; numHe <= maxClusterSize; numHe++) {
		// Create a He cluster with cluster size numHe
		std::shared_ptr<HeCluster> cluster(new HeCluster(numHe));
		
		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add vacancy clusters
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		// Create a He cluster with cluster size numV
		std::shared_ptr<VCluster> cluster(new VCluster(numV));
		
		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add interstitial clusters
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		// Create a He cluster with cluster size numI
		std::shared_ptr<InterstitialCluster> cluster(
			new InterstitialCluster(numI));
		
		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add HeV clusters
	// I am assuming that HeV clusters can only be created if
	// numV <= numHe.
	
	int numHeVClusters = 0;
	
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = numV; numHe + numV <= maxClusterSize; numHe++) {
			std::shared_ptr<HeVCluster> cluster(new HeVCluster(numHe, numV));
			
			numHeVClusters++;
			reactants->push_back(cluster);
		}
	}
	
	// Setup the properties map
	(*properties)["maxHeClusterSize"] = std::to_string(maxClusterSize);
	(*properties)["maxVClusterSize"] = std::to_string(maxClusterSize);
	(*properties)["maxIClusterSize"] = std::to_string(maxClusterSize);
	(*properties)["maxMixedClusterSize"] = std::to_string(numHeVClusters);
	
	(*properties)["numHeClusters"] = std::to_string(numClusters);
	(*properties)["numVClusters"] = std::to_string(numClusters);
	(*properties)["numIClusters"] = std::to_string(numClusters);
	(*properties)["numHeVClusters"] = std::to_string(numHeVClusters);
}

SimpleReactionNetwork::~SimpleReactionNetwork() {
	// Nothing to do
}

/**
 * This operation creates a SimpleReactionNetwork and makes sure that it is
 * properly registered with the clusters it contains. This operation should
 * always be called instead of constructing a SimpleReactionNetwork manually.
 * @return The reaction network.
 */
std::shared_ptr<xolotlCore::ReactionNetwork> testUtils::getSimpleReactionNetwork() {
	// Create the network
	std::shared_ptr<xolotlCore::ReactionNetwork> network(new SimpleReactionNetwork());
	// Register the reaction network with its clusters
	for (int i = 0; i < network->reactants->size(); i++) {
		network->reactants->at(i)->setReactionNetwork(network);
	}
	return network;
}
