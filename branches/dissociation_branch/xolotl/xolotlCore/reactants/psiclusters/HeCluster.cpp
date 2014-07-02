// Includes
#include "HeCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

HeCluster::HeCluster(int nHe,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(nHe, registry) {
	// Update the composition map
	compositionMap["He"] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "He";

	// Compute the reaction radius
	double FourPi = 4.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::latticeConstant, 3);
	double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * size,
			(1.0 / 3.0));
	double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
	reactionRadius = 0.3 + termOne - termTwo;

	return;
}

HeCluster::~HeCluster() {
}

std::shared_ptr<Reactant> HeCluster::clone() {
	std::shared_ptr<Reactant> reactant(new HeCluster(*this));
	return reactant;
}

void HeCluster::createReactionConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species reactions
	PSICluster::createReactionConnectivity();

	// Helium-Vacancy clustering
	// He_a + V_b --> (He_a)(V_b)
	// Get all the V clusters from the network
	auto reactants = network->getAll(vType);
	// combineClusters handles V combining with He to form HeV
	combineClusters(reactants, heVType);

	// Helium-Interstitial clustering
	// He_a + I_b --> (He_a)(I_b)
	// Get all the I clusters from the network
	reactants = network->getAll(iType);
	// combineClusters handles I combining with He to form HeI
	combineClusters(reactants, heIType);

	// Helium absorption by HeV clusters
	// He_a + (He_b)(V_c) --> [He_(a+b)](V_c)
	// Get all the HeV clusters from the network
	reactants = network->getAll(heVType);
	// combineClusters handles HeV combining with He to form HeV
	combineClusters(reactants, heVType);

	// Helium absorption by HeI clusters
	// He_a + (He_b)(I_c) --> [He_(a+b)](I_c)
	// Get all the HeI clusters from the network
	reactants = network->getAll(heIType);
	// combineClusters handles HeI combining with He to form HeI
	combineClusters(reactants, heIType);

	return;
}

void HeCluster::createDissociationConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species dissociation
	PSICluster::createDissociationConnectivity();

	// Specific case for the single species cluster
	if (size == 1) {
		// He dissociation of HeV cluster is handled here
		// (He_a)(V_b) --> [He_(a-1)](V_b) + He_1
		// Get all the HeV clusters of the network
		auto allHeVReactants = network->getAll(heVType);
		for (int i = 0; i < allHeVReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeVReactants.at(i);

			// (He_a)(V_b) is the dissociating one, [He_(a-1)](V_b) is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			std::vector<int> compositionVec = { comp[heType] - 1, comp[vType],
					comp[iType] };
			auto smallerReactant = network->getCompound(heVType, compositionVec);
			dissociateCluster(allHeVReactants.at(i), smallerReactant);
		}

		// He dissociation of HeI cluster is handled here
		// (He_a)(I_b) --> [He_(a-1)](I_b) + He_1
		// Get all the HeI clusters of the network
		auto allHeIReactants = network->getAll(heIType);
		for (int i = 0; i < allHeIReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeIReactants.at(i);

			// (He_a)(I_b) is the dissociating one, [He_(a-1)](I_b) is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			std::vector<int> compositionVec = { comp[heType] - 1, comp[vType],
					comp[iType] };
			auto smallerReactant = network->getCompound(heIType, compositionVec);
			dissociateCluster(allHeVReactants.at(i), smallerReactant);
		}
	}

	return;
}
