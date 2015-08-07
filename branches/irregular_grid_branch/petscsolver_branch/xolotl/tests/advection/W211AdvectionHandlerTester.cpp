#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <W211AdvectionHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the W211AdvectionHandler.
 */
BOOST_AUTO_TEST_SUITE(W211AdvectionHandler_testSuite)

/**
 * Method checking the initialization and the compute advection methods.
 */
BOOST_AUTO_TEST_CASE(checkAdvection) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();
	// Get its size
	const int size = network->getAll()->size();

	// Create ofill
	int mat[size*size];
	int *ofill = &mat[0];

	// Create the advection handler and initialize it
	W211AdvectionHandler advectionHandler;
	advectionHandler.initialize(network, ofill);

	// Check the total number of advecting clusters
	BOOST_REQUIRE_EQUAL(advectionHandler.getNumberOfAdvecting(), 6);

	// Check the clusters in ofill
	BOOST_REQUIRE_EQUAL(ofill[0], 1);
	BOOST_REQUIRE_EQUAL(ofill[10], 1);
	BOOST_REQUIRE_EQUAL(ofill[20], 1);
	BOOST_REQUIRE_EQUAL(ofill[30], 1);
	BOOST_REQUIRE_EQUAL(ofill[40], 1);
	BOOST_REQUIRE_EQUAL(ofill[50], 1);

	// Set the size parameter in the x direction
	double hx = 1.0;

	// Create the arrays of concentration
	double concentration[3*size];
	double newConcentration[3*size];

	// Initialize their values
	for (int i = 0; i < 3*size; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000K to initialize the diffusion coefficients
	auto reactants = network->getAll();
	for (int i = 0; i < size; i++) {
		auto cluster = (PSICluster *) reactants->at(i);
		cluster->setTemperature(1000.0);
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	double *concOffset = conc + size;
	double *updatedConcOffset = updatedConc + size;

	// Fill the concVector with the pointer to the middle, left, and right grid points
	double **concVector = new double*[3];
	concVector[0] = concOffset; // middle
	concVector[1] = conc; // left
	concVector[2] = conc + 2 * size; // right

	// Set the grid position
	std::vector<double> gridPosition = { hx, 0.0, 0.0 };

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(network, gridPosition,
			concVector, updatedConcOffset, hx, hx);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], -2.02166e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], -3.12551e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2], -5.21336e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], -9.54688e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4], -1.14847e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[5], -6.24245e+09, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[6], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	int indices[nAdvec];
	double val[2*nAdvec];
	// Get the pointer on them for the compute advection method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(network, valPointer,
			indicesPointer, gridPosition, hx, hx);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 0);
	BOOST_REQUIRE_EQUAL(indices[1], 1);
	BOOST_REQUIRE_EQUAL(indices[2], 2);
	BOOST_REQUIRE_EQUAL(indices[3], 3);
	BOOST_REQUIRE_EQUAL(indices[4], 4);
	BOOST_REQUIRE_EQUAL(indices[5], 5);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0],-332783240.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 20798953.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], -403616596.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 25226037.0, 0.01);

	// Get the stencil
	auto stencil = advectionHandler.getStencilForAdvection(gridPosition);

	// Check the value of the stencil
	BOOST_REQUIRE_EQUAL(stencil[0], 1); //x
	BOOST_REQUIRE_EQUAL(stencil[1], 0);
	BOOST_REQUIRE_EQUAL(stencil[2], 0);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
