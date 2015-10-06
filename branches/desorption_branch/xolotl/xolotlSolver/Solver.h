#ifndef SOLVER_H
#define SOLVER_H

// Includes
#include "ISolver.h"

namespace xolotlSolver {

/**
 * This class and its subclasses realize the ISolver interface to solve the
 * advection-diffusion-reaction problem with currently supported solvers.
 */
class Solver: public ISolver {
protected:

	//! The number command line arguments
	int numCLIArgs;

	//! The command line arguments
	char **CLIArgs;

	//! The network loader that can load the reaction network data.
	PSIClusterNetworkLoader *networkLoader;

	//! The original network created from the network loader.
	PSIClusterReactionNetwork *network;

	//! The original solver handler.
	static ISolverHandler *solverHandler;

public:

	//! The Constructor
	Solver(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	//! The Destructor
	virtual ~Solver();

	/**
	 * This operation transfers the input arguments passed to the program on
	 * startup to the solver. These options are static options specified at
	 * the start of the program whereas the options passed to setOptions() may
	 * change.
	 * @param argc The number of command line arguments
	 * @param argv The array of command line arguments
	 */
	void setCommandLineOptions(int argc, char **argv);

	/**
	 * This operation sets the PSIClusterNetworkLoader that should be used by
	 * the ISolver to load the ReactionNetwork.
	 * @param networkLoader The PSIClusterNetworkLoader that will load the
	 * network.
	 */
	void setNetworkLoader(
			std::shared_ptr<PSIClusterNetworkLoader> networkLoader);

	/**
	 * This operation returns the solver handler for this solver. This
	 * operation is only for use by solver code and is not part of the
	 * ISolver interface.
	 * @return The advection handler for this solver
	 */
	static ISolverHandler *getSolverHandler() {
		return solverHandler;
	}

protected:

    /**
     * The performance handler registry that will be used
     * for this class.
     */
    std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry;

}; //end class Solver

} /* end namespace xolotlSolver */
#endif
