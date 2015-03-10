#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "W111FitFluxHandler.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W111FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {
	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double) l * 1.25);
	}

	// Create the W111 flux handler
    auto testFitFlux = make_shared<W111FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(grid);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 0.391408, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 0.268688, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 0.139904, 0.01);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
