/**
 * Main.c, currently only able to load clusters
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <Reactant.h>
#include <PSIClusterNetworkLoader.h>
#include <PetscSolver.h>
#include <mpi.h>
#include <MPIUtils.h>
#include <XolotlOptions.h>
#include <HandlerRegistryFactory.h>
#include <HardwareQuantities.h>
#include <IHandlerRegistry.h>
#include <HDF5NetworkLoader.h>

using namespace std;
using std::shared_ptr;

//! This operation prints the start message
void printStartMessage() {
	cout << "Starting Xolotl Plasma-Surface Interactions Simulator" << endl;
	// TODO! Print copyright message
	// TODO! Print date and time
}

std::vector<xolotlPerf::HardwareQuantities> declareHWcounters() {

	// Indicate we want to monitor some important hardware counters.
	std::vector<xolotlPerf::HardwareQuantities> hwq;

	hwq.push_back(xolotlPerf::FP_OPS);
	hwq.push_back(xolotlPerf::L1_CACHE_MISS);

	return hwq;
}

bool initPerf(bool opts, std::vector<xolotlPerf::HardwareQuantities> hwq) {

	bool perfInitOK = xolotlPerf::initialize(opts, hwq);
	if (!perfInitOK) {
		std::cerr
				<< "Unable to initialize requested performance data infrastructure. "
				<< "Aborting" << std::endl;
		return EXIT_FAILURE;
	} else
		return perfInitOK;
}

std::shared_ptr<xolotlSolver::PetscSolver>
setUpSolver( std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry, 
            int argc, char **argv) {

	// Setup the solver
	auto solverInitTimer = handlerRegistry->getTimer("initSolver");
	solverInitTimer->start();
	std::shared_ptr<xolotlSolver::PetscSolver> solver = 
        std::make_shared<xolotlSolver::PetscSolver>(handlerRegistry);
	solver->setCommandLineOptions(argc, argv);
	solver->initialize();
	solverInitTimer->stop();

	return solver;
}

void launchPetscSolver(std::shared_ptr<xolotlSolver::PetscSolver> solver,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry) {

	// Launch the PetscSolver
	auto solverTimer = handlerRegistry->getTimer("solve");
	solverTimer->start();
	// Create object to handle incident flux calculations
	auto fitFluxHandler = std::make_shared<xolotlSolver::FitFluxHandler>();
	// Create object to handle temperature
	auto tempHandler = std::make_shared<xolotlSolver::TemperatureHandler>();
	solver->solve(fitFluxHandler, tempHandler);
	solverTimer->stop();

}

std::shared_ptr<PSIClusterNetworkLoader> setUpNetworkLoader(int rank,
		MPI_Comm comm, std::string networkFilename,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {

	// Create a HDF5NetworkLoader
	std::shared_ptr<HDF5NetworkLoader> networkLoader;
	networkLoader = std::make_shared<HDF5NetworkLoader>(registry);
	// Give the networkFilename to the network loader
	networkLoader->setFilename(networkFilename);

	return networkLoader;
}

//! Main program
int main(int argc, char **argv) {

	// Local Declarations
	int rank;

	// Check the command line arguments.
	// Skip the executable name before parsing.
	argc -= 1;  // one for the executable name
	argv += 1;  // one for the executable name
	XolotlOptions xopts;
	int nOptsUsed = xopts.parseCommandLine(argc, argv);
	if (!xopts.shouldRun()) {
		return xopts.getExitCode();
	}
	argc -= nOptsUsed;
	argv += nOptsUsed;

	// Extract the argument for the file name
	std::string networkFilename = xopts.getNetworkFilename();
	assert(!networkFilename.empty());

	try {
		// Set up our performance data infrastructure.
		// Indicate we want to monitor some important hardware counters.
		auto hwq = declareHWcounters();
		auto perfInitOK = initPerf(xopts.useStandardHandlers(), hwq);

		// Initialize MPI.  We do this instead of leaving it to some
		// other package (e.g., PETSc), because we want to avoid problems
		// with overlapping Timer scopes.
		MPI_Init(&argc, &argv);

		// Access our handler registry to obtain a Timer
		// measuring the runtime of the entire program.
		// NOTE: these long template types could be replaced with 'auto'
		auto handlerRegistry = xolotlPerf::getHandlerRegistry();
		auto totalTimer = handlerRegistry->getTimer("total");
		totalTimer->start();

		// Setup the solver
		auto solver = setUpSolver(handlerRegistry, argc, argv);

		// Load the network
		auto networkLoadTimer = handlerRegistry->getTimer("loadNetwork");
		networkLoadTimer->start();

		// Get the MPI rank
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		// Set up the network loader
		auto networkLoader = setUpNetworkLoader(rank, MPI_COMM_WORLD,
				networkFilename, handlerRegistry);

		// Give the network loader to PETSc as input
		solver->setNetworkLoader(networkLoader);
		networkLoadTimer->stop();

		// Launch the PetscSolver
		launchPetscSolver(solver, handlerRegistry);

		// Finalize our use of the solver.
		auto solverFinalizeTimer = handlerRegistry->getTimer("solverFinalize");
		solverFinalizeTimer->start();
		solver->finalize();
		solverFinalizeTimer->stop();
		totalTimer->stop();

		// Report the performance data about the run we just completed
		// TODO Currently, this call writes EventCounter data to the
		// given stream, but Timer and any hardware counter data is
		// written by the underlying timing library to files, one per process.
		if (rank == 0) {
			handlerRegistry->dump(std::cout);
		}

	} catch (std::string & error) {
		std::cout << error << std::endl;
		std::cout << "Aborting." << std::endl;
		return EXIT_FAILURE;
	}

	// finalize our use of MPI
	MPI_Finalize();

	// Uncomment if GPTL was built with pmpi disabled
	// Output performance data if pmpi is disabled in GPTL
	// Access the handler registry to output performance data
//    auto handlerRegistry = xolotlPerf::getHandlerRegistry();
//    handlerRegistry->dump(rank);

	return EXIT_SUCCESS;
}
