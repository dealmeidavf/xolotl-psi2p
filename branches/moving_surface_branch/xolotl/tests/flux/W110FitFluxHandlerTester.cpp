#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "W110FitFluxHandler.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W110FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {
	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;
	// Specify the surface position
	int surfacePos = 0;

    auto testFitFlux = make_shared<W110FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(surfacePos, nGridpts, step);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime, surfacePos);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 0.524627, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 0.211160, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 0.064213, 0.01);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
