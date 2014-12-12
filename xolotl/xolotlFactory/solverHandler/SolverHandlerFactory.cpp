#include <SolverHandlerFactory.h>
#include <PetscSolver1DHandler.h>
#include <fstream>
#include <iostream>
#include <mpi.h>

namespace xolotlFactory {

static std::shared_ptr<xolotlSolver::ISolverHandler> theSolverHandler;

// Create the desired type of handler registry.
bool initializeDimension(xolotlCore::Options &options) {

	bool ret = true;

	// Get the wanted dimension
	int dim = options.getDimensionNumber();

	// Switch on the dimension
	switch (dim) {
		case 1:
			theSolverHandler = std::make_shared<xolotlSolver::PetscSolver1DHandler>();
			break;
		case 2:
			// To be implemented
			throw std::string("\nxolotlFactory: 2D solver handler is not implemented yet.");
		case 3:
			// To be implemented
			throw std::string("\nxolotlFactory: 3D solver handler is not implemented yet.");
		default:
			// The asked dimension is not good (e.g. -1, 0, 4)
			throw std::string("\nxolotlFactory: Bad dimension for the solver handler.");
	}

	return ret;
}

// Provide access to our handler registry.
std::shared_ptr<xolotlSolver::ISolverHandler> getSolverHandler() {
	if (!theSolverHandler) {
		// We have not yet been initialized.
		throw std::string("\nxolotlFactory: solver requested but "
				"it has not been initialized.");
	}

	return theSolverHandler;
}

} // end namespace xolotlFactory

