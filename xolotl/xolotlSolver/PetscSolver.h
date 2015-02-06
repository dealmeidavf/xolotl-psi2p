#ifndef PETSCSOLVER_H
#define PETSCSOLVER_H

// Includes
#include "Solver.h"

namespace xolotlSolver {

/**
 * This class realizes the ISolver interface to solve the
 * advection-diffusion-reaction problem with the PETSc solvers from Argonne
 * National Laboratory.
 */
class PetscSolver: public Solver {
private:

	/**
	 * This operation configures the initial conditions of the grid in Xolotl.
	 * @param data The DM (data manager) created by PETSc
	 * @param solutionVector The solution vector that contains the PDE
	 * solution and which needs to be initialized.
	 * @return The error code. 0 if there is no error.
	 */
	PetscErrorCode setupInitialConditions(DM data, Vec solutionVector);

public:

	//! The Constructor
	PetscSolver(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	//! The Destructor
	~PetscSolver();

	/**
	 * This operation sets the run-time options of the solver. The map is a set
	 * of key-value std::string pairs that are interpreted by the solver. These
	 * options may change during execution, but it is up to Solvers to monitor
	 * the map for changes and they may do so at their discretion.
	 * @param options The set of options as key-value pairs with option names
	 * for keys and associated values mapped to those keys. A relevant example
	 * is "startTime" and "0.01" where both are of type std::string.
	 */
	void setOptions(std::map<std::string, std::string> options);

	/**
	 * This operation sets up the mesh that will be used by the solver and
	 * initializes the data on that mesh. This operation will throw an exception
	 * of type std::string if the mesh can not be setup.
	 */
	void setupMesh();

	/**
	 * This operation performs all necessary initialization for the solver
	 * possibly including but not limited to setting up MPI and loading initial
	 * conditions. If the solver can not be initialized, this operation will
	 * throw an exception of type std::string.
	 * @param solverHandler The solver handler
	 */
	void initialize(std::shared_ptr<ISolverHandler> solverHandler);

	/**
	 * This operation directs the Solver to perform the solve. If the solve
	 * fails, it will throw an exception of type std::string.
	 */
	void solve();

	/**
	 * This operation performs all necessary finalization for the solver
	 * including but not limited to cleaning up memory, finalizing MPI and
	 * printing diagnostic information. If the solver can not be finalized,
	 * this operation will throw an exception of type std::string.
	 */
	void finalize();
	
}; //end class PetscSolver

} /* end namespace xolotlSolver */

// Some compilers (e.g., recent versions of Intel) define __func__ 
// to include the namespace or class scope when compiled with the C++11
// support enabled.  Others don't.  Because PETSc's PetscFunctionBeginUser
// does a straight string comparison between what we call the function name
// and what it determines from the compiler, we need a way to provide
// either the scoped name or the unscoped name.
#if defined(__ICC) || defined(__INTEL_COMPILER)
#  define Actual__FUNCT__(sname,fname)  sname "::" fname
#else
#  define Actual__FUNCT__(sname,fname)  fname
#endif /* if it is the Intel compiler */

#endif
