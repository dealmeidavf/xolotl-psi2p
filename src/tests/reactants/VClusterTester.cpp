/*
 * PSIClusterTester.cpp
 *
 *  Created on: May 6, 2013
 *      Author: Jay Jay Billings
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include <VCluster.h>
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the VCluster.
 */BOOST_AUTO_TEST_SUITE(VCluster_testSuite)

/**
 * This operation checks the ability of the VCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	shared_ptr<ReactionNetwork> network = getSimplePSIReactionNetwork();

	// Prevent dissociation from being added to the connectivity array
	network->disableDissociations();

	// Get the connectivity array from the reactant for a vacancy cluster of size 2.
	auto reactant = (PSICluster *) network->get("V", 2);
	// Check the type name
	BOOST_REQUIRE_EQUAL("V", reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 1, 1, 1, 0, 0,

			// V
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,

			// I
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,

			// HeV
			0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// HeI
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the VCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimplePSIReactionNetwork();

	// Get an V cluster with compostion 0,1,0.
	auto cluster = (PSICluster *) network->get("V", 1);
	// Get one that it combines with (V2)
	auto secondCluster = (PSICluster *) network->get("V", 2);
	// Set the diffusion factor, and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(2.41E+11);
	cluster->setMigrationEnergy(1.66);
	cluster->setTemperature(1000.0);
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(0.0);
	secondCluster->setMigrationEnergy(numeric_limits<double>::infinity());
	secondCluster->setConcentration(0.5);
	secondCluster->setTemperature(1000.0);

	// Compute the rate constants that are needed for the flux
	cluster->computeRateConstants();
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux();
	BOOST_TEST_MESSAGE(
			"InterstitialClusterTester Message: \n" << "Total Flux is " << flux << "\n" << "   -Production Flux: " << cluster->getProductionFlux() << "\n" << "   -Combination Flux: " << cluster->getCombinationFlux() << "\n" << "   -Dissociation Flux: " << cluster->getDissociationFlux() << "\n" << "   -Emission Flux: " << cluster->getEmissionFlux() << "\n");

	BOOST_REQUIRE_CLOSE(444828.3, flux, 0.1);

	return;
}

/**
 * This operation checks the VCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { -2850.42, -3005.08, 0.0, -14316.7, 896815.0,
			257925.0, -1923.81, -2188.27, -2373.78, 356134.7, 377344.4, 224717.0,
			0.0, 0.0, -2054.05 };
	// Get the simple reaction network
	shared_ptr<ReactionNetwork> network = getSimplePSIReactionNetwork(3);

	// Get an V cluster with compostion 0,1,0.
	auto cluster = (PSICluster *) network->get("V", 1);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(2.41E+11);
	cluster->setMigrationEnergy(1.66);
	cluster->setTemperature(1000.0);
	cluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the partial derivatives
	cluster->computeRateConstants();
	// Reinitialize the network for Ids for the partial derivatives
	network->reinitializeNetwork();
	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives();

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 15U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	return;
}

/**
 * This operation checks the reaction radius for VCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the vacancy clsuter
	shared_ptr<VCluster> cluster;

	// The vector of radii to compare with
	double expectedRadii[] = { 0.1372650265, 0.1778340462, 0.2062922619,
			0.2289478080, 0.2480795532 };

	// Check all the values
	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr<VCluster>(new VCluster(i, registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.000001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()

