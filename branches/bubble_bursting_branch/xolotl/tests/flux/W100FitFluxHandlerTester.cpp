#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "W100FitFluxHandler.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W100FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;
	// Specify the surface position
	int surfacePos = 0;

    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step, surfacePos);

	// Create a composition vector
	vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	// Create a vector representing the position of the cluster
	vector<double> x = {1.25, 0.0, 0.0};

	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1, surfacePos);

	BOOST_TEST_MESSAGE( "\nW100FitFluxHandlerTester Message: \n"
						<< "incidentFlux = " << testFlux << " with position "
						<< "(" << x[0] << "," << x[1] << "," << x[2] << ") "
						<< "at time = " << currTime << "\n");
	BOOST_REQUIRE_CLOSE(testFlux, 0.476819, 0.01);
}

BOOST_AUTO_TEST_CASE(checkHeFluence) {

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;
	// Specify the surface position
	int surfacePos = 0;

    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step, surfacePos);

	// Create the composition vector
	vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	// Create a vector representing the position of the cluster
	vector<double> x = {1.25, 0.0, 0.0};

	// Check the flux
	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1, surfacePos);
	BOOST_REQUIRE_CLOSE(testFlux, 0.476819, 0.01);

	// Check that the fluence is 0.0 at the beginning
	BOOST_REQUIRE_EQUAL(testFitFlux->getHeFluence(), 0.0);

	// Increment the helium fluence
	testFitFlux->incrementHeFluence(1.0e-8);
	// Check that the fluence is not 0.0 anymore
	BOOST_REQUIRE_EQUAL(testFitFlux->getHeFluence(), 1.0e-8);
}

BOOST_AUTO_TEST_CASE(checkHeFlux) {

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;
	// Specify the surface position
	int surfacePos = 0;

    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Set the factor to change the Helium flux
    testFitFlux->setHeFlux(2.5);
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step, surfacePos);

    BOOST_REQUIRE_EQUAL(testFitFlux->getHeFlux(), 2.5);

	// Create a composition vector
	vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	// Create a vector representing the position of the cluster
	vector<double> x = {1.25, 0.0, 0.0};

	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1, surfacePos);
	BOOST_REQUIRE_CLOSE(testFlux, 2.5 * 0.476819, 0.01);
}

BOOST_AUTO_TEST_SUITE_END()
