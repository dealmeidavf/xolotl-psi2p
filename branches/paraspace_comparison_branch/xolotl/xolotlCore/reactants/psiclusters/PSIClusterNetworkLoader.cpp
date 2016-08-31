/*
 * PSIClusterNetworkLoader.cpp
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#include "PSIClusterNetworkLoader.h"
#include <TokenizedLineReader.h>
#include <stdio.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <HeCluster.h>
#include <VCluster.h>
#include <InterstitialCluster.h>
#include <HeVCluster.h>
// #include <HeInterstitialCluster.h>
#include "PSIClusterReactionNetwork.h"
#include <xolotlPerf.h>

using namespace xolotlCore;

/**
 * This operation converts a string to a double, taking in to account the fact
 * that the input file may contain keys such as "infinite."
 *
 * @param inString the string to be converted
 * @return the string as a double
 */
static inline double convertStrToDouble(const std::string& inString) {
	return (inString.compare("infinite") == 0) ?
			std::numeric_limits<double>::infinity() :
			strtod(inString.c_str(), NULL);
}

std::shared_ptr<PSICluster> PSIClusterNetworkLoader::createCluster(int numHe,
		int numV, int numI) {
	// Local Declarations
	std::shared_ptr<PSICluster> cluster;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numHe > 0 && numV > 0) {
		// Create a new HeVCluster
		cluster = std::make_shared<HeVCluster>(numHe, numV, handlerRegistry);

		clusters.push_back(cluster);
	}
	else if (numHe > 0 && numI > 0) {
		throw std::string("HeliumInterstitialCluster is not yet implemented.");
		// FIXME! Add code to add it to the list
	}
	else if (numHe > 0) {
		// Create a new HeCluster
		cluster = std::make_shared<HeCluster>(numHe, handlerRegistry);

		clusters.push_back(cluster);
	}
	else if (numV > 0) {
		// Create a new VCluster
		cluster = std::make_shared<VCluster>(numV, handlerRegistry);

		clusters.push_back(cluster);
	}
	else if (numI > 0) {
		// Create a new ICluster
		cluster = std::make_shared<InterstitialCluster>(numI, handlerRegistry);

		// Add it to the ICluster list
		clusters.push_back(cluster);
	}

	return cluster;
}

PSIClusterNetworkLoader::PSIClusterNetworkLoader(
		const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
	handlerRegistry(registry) {
	setInputstream(stream);

	return;
}

void PSIClusterNetworkLoader::setInputstream(
		const std::shared_ptr<std::istream> stream) {
	networkStream = stream;

	return;
}

std::shared_ptr<PSIClusterReactionNetwork> PSIClusterNetworkLoader::load() {
	// Local Declarations
	TokenizedLineReader<std::string> reader;
	std::vector<std::string> loadedLine;
	std::shared_ptr<PSIClusterReactionNetwork> network = std::make_shared<
			PSIClusterReactionNetwork>(handlerRegistry);

	std::string error(
			"PSIClusterNetworkLoader Exception: Insufficient or erroneous data.");
	int numHe = 0, numV = 0, numI = 0;
	double heBindingEnergy = 0.0, vBindingEnergy = 0.0, iBindingEnergy = 0.0,
			migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::shared_ptr<Reactant> > reactants;

	// Load the network if the stream is available
	if (networkStream != NULL) {
		// Load the stream
		reader.setInputStream(networkStream);

		// Loop over each line of the file, which should each be PSIClusters.
		loadedLine = reader.loadLine();
		while (loadedLine.size() > 0) {
			// Check the size of the loaded line
			if (loadedLine.size() < 8)
				// And notify the calling function if the size is insufficient
				throw error;
			// Load the sizes
			if (loadedLine[0][0] != '#') {
				numHe = std::stoi(loadedLine[0]);
				numV = std::stoi(loadedLine[1]);
				numI = std::stoi(loadedLine[2]);
				// Create the cluster
				auto nextCluster = createCluster(numHe, numV, numI);
				// Load the energies
				heBindingEnergy = convertStrToDouble(loadedLine[3]);
				vBindingEnergy = convertStrToDouble(loadedLine[4]);
				iBindingEnergy = convertStrToDouble(loadedLine[5]);
				migrationEnergy = convertStrToDouble(loadedLine[6]);
				diffusionFactor = convertStrToDouble(loadedLine[7]);
				// Set the formation energy
				nextCluster->setBindingEnergy(heBindingEnergy, vBindingEnergy, iBindingEnergy);
				// Set the diffusion factor and migration energy
				nextCluster->setMigrationEnergy(migrationEnergy);
				nextCluster->setDiffusionFactor(diffusionFactor);
				// Add the cluster to the network
				network->add(nextCluster);
				// Add it to the list so that we can set the network later
				reactants.push_back(nextCluster);
			}

			// Load the next line
			loadedLine = reader.loadLine();
		}

		// Set the network for all of the reactants. This MUST be done manually.
		for (auto reactantsIt = reactants.begin();
				reactantsIt != reactants.end(); ++reactantsIt) {
			(*reactantsIt)->setReactionNetwork(network);
		}
	}

	return network;
}
