// Includes
#include "HeCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

HeCluster::HeCluster(int nHe,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(registry) {
	// Set the size
	size = nHe;
	// Update the composition map
	compositionMap[heType] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = heType;

	// Compute the reaction radius
	double FourPi = 4.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::tungstenLatticeConstant, 3);
	double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * size,
			(1.0 / 3.0));
	double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
	reactionRadius = 0.3 + termOne - termTwo;

	return;
}

void HeCluster::combineClusters(std::vector<IReactant *> & clusters,
		const std::string& productName) {
	// Initial declarations
	std::map<std::string, int> secondComposition;

	// Get all the I clusters to loop on them
	auto iClusters = network->getAll(iType);

	// Loop on the potential combining reactants
	for (unsigned int i = 0; i < clusters.size(); i++) {
		// Get the second reactant, its composition and its index
		auto secondCluster = (PSICluster *) clusters[i];
		secondComposition = secondCluster->getComposition();
		// Check that the simple product [He_(a+c)](V_b) doesn't exist
		// b can be 0 so the simple product would be a helium cluster
		PSICluster * simpleProduct;
		if (secondComposition[vType] == 0) {
			simpleProduct = (PSICluster *) network->get(heType,
					size + secondComposition[heType]);
		} else {
			std::vector<int> comp = { size + secondComposition[heType],
					secondComposition[vType], secondComposition[iType] };
			simpleProduct = (PSICluster *) network->getCompound(productName,
					comp);
		}
		if (simpleProduct)
			continue;
		// The simple product doesn't exist so it will go though trap-mutation
		// The reaction is
		// (He_c)(V_b) + He_a --> [He_(a+c)][V_(b+d)] + I_d
		// Loop on the possible I starting by the smallest
		for (auto it = iClusters.begin(); it != iClusters.end(); it++) {
			// Get the size of the I cluster
			int iSize = (*it)->getSize();

			// Create the composition of the potential product
			std::vector<int> comp = { size + secondComposition[heType],
					secondComposition[vType] + 1, secondComposition[iType] };
			auto firstProduct = (PSICluster *) network->getCompound(productName,
					comp);
			// If the first product exists
			if (firstProduct) {
				// This cluster combines with the second reactant
				setReactionConnectivity(secondCluster->getId());
				// Creates the combining cluster
				// The reaction constant will be computed later and is set to 0.0 for now
				CombiningCluster combCluster(secondCluster, 0.0);
				// Push the product onto the list of clusters that combine with this one
				combiningReactants.push_back(combCluster);
				// Break the loop because we want only one reaction
				break;
			}
		}
	}

	return;
}

void HeCluster::createReactionConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species reactions
	PSICluster::createReactionConnectivity();

	// This cluster is always He_a

	// Helium-Vacancy clustering
	// He_a + V_b --> (He_a)(V_b)
	// Get all the V clusters from the network
	auto reactants = network->getAll(vType);
	// combineClusters handles V combining with He to form HeV
	PSICluster::combineClusters(reactants, heVType);

	// Helium-Interstitial clustering
	// He_a + I_b --> (He_a)(I_b)
	// Get all the I clusters from the network
	reactants = network->getAll(iType);
	// combineClusters handles I combining with He to form HeI
	PSICluster::combineClusters(reactants, heIType);

	// Helium absorption by HeV clusters
	// He_a + (He_b)(V_c) --> [He_(a+b)](V_c)
	// Get all the HeV clusters from the network
	reactants = network->getAll(heVType);
	// combineClusters handles HeV combining with He to form HeV
	PSICluster::combineClusters(reactants, heVType);

	// Helium absorption by HeI clusters
	// He_a + (He_b)(I_c) --> [He_(a+b)](I_c)
	// Get all the HeI clusters from the network
	reactants = network->getAll(heIType);
	// combineClusters handles HeI combining with He to form HeI
	PSICluster::combineClusters(reactants, heIType);

	// Helium absorption leading to trap mutation
	// (He_c)(V_b) + He_a --> [He_(a+c)][V_(b+d)] + I_d
	// Get all the HeV clusters from the network
	reactants = network->getAll(heVType);
	// HeCluster::combineClusters handles He combining with HeV to go through trap-mutation
	combineClusters(reactants, heVType);
	// b can be 0 so He clusters can combine with He clusters leading to trap-mutation
	// Get all the He clusters from the network
	reactants = network->getAll(heType);
	combineClusters(reactants, heVType);

	return;
}

void HeCluster::createDissociationConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species dissociation
	PSICluster::createDissociationConnectivity();

	// This cluster is always He_a

	// Vacancy Dissociation
	// (He_a)(V_b) --> He_(a)[V_(b-1)] + V
	// for b == 1
	// Get the HeV cluster
	std::vector<int> compositionVec = { size, 1, 0 };
	auto biggerReactant = (PSICluster *) network->getCompound(heVType,
			compositionVec);
	// Get the single vacancy
	auto singleVacancy = (PSICluster *) network->get(vType, 1);
	// Dissociate
	dissociateCluster(biggerReactant, singleVacancy);

	// Specific case for the single species cluster
	if (size == 1) {
		// He dissociation of HeV cluster is handled here
		// (He_b)(V_c) --> [He_(b-a)](V_c) + He_a
		// for a = 1
		// Get all the HeV clusters of the network
		auto allHeVReactants = network->getAll(heVType);
		for (unsigned int i = 0; i < allHeVReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeVReactants[i];

			// (He_b)(V_c) is the dissociating one, [He_(b-a)](V_c) is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();

			// Skip He_1V_1 because it was counted in the V dissociation
			if (comp[heType] == 1 && comp[vType] == 1)
				continue;

			std::vector<int> compositionVec = { comp[heType] - size,
					comp[vType], 0 };
			auto smallerReactant = (PSICluster *) network->getCompound(heVType,
					compositionVec);
			// Special case for comp[heType] = 1
			if (comp[heType] == 1) {
				smallerReactant = (PSICluster *) network->get(vType,
						comp[vType]);
			}
			dissociateCluster(cluster, smallerReactant);
		}

		// He dissociation of HeI cluster is handled here
		// (He_b)(I_c) --> [He_(b-a)](I_c) + He_a
		// for a = 1
		// Get all the HeI clusters of the network
		auto allHeIReactants = network->getAll(heIType);
		for (unsigned int i = 0; i < allHeIReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeIReactants[i];

			// (He_b)(I_c) is the dissociating one, [He_(b-a)](I_c) is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			std::vector<int> compositionVec =
					{ comp[heType] - 1, 0, comp[iType] };
			auto smallerReactant = (PSICluster *) network->getCompound(heIType,
					compositionVec);
			dissociateCluster(cluster, smallerReactant);
		}
	}

	return;
}
