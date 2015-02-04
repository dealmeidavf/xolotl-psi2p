#ifndef PETSCSOLVER2DHANDLER_H
#define PETSCSOLVER2DHANDLER_H

// Includes
#include "PetscSolverHandler.h"

namespace xolotlSolver {

/**
 * This class is a subclass of PetscSolverHandler and implement all the methods needed
 * to solve the ADR equations in 2D using PETSc from Argonne National Laboratory.
 */
class PetscSolver2DHandler: public PetscSolverHandler {
private:
	//! The position of the surface
	std::vector<int> surfacePosition;

	/**
	 * Get the mean position of the surface.
     *
     * @return The mean position of the surface
	 */
	int getMeanSurfacePosition() {
		int size = surfacePosition.size();

		int mean = 0;
		// Compute the mean
		for (int l = 1; l < size - 1; l++) {
			mean += surfacePosition[l];
		}

		return mean / (size - 2);
	}

public:

	//! The Constructor
	PetscSolver2DHandler() {}

	//! The Destructor
	~PetscSolver2DHandler() {}

	/**
	 * Create everything needed before starting to solve.
     * \see ISolverHandler.h
	 */
	void createSolverContext(DM &da, int nx, double hx, int ny,
			double hy, int nz, double hz);

	/**
	 * Initialize the concentration solution vector.
     * \see ISolverHandler.h
	 */
	void initializeConcentration(DM &da, Vec &C) const;

	/**
	 * Compute the new concentrations for the RHS function given an initial
	 * vector of concentrations. Apply the diffusion, advection and all the reactions.
     * \see ISolverHandler.h
	 */
	void updateConcentration(TS &ts, Vec &localC, Vec &F, PetscReal ftime,
			bool &temperatureChanged, bool &hasMoved);

	/**
	 * Compute the off-diagonal part of the Jacobian which is related to cluster's motion.
     * \see ISolverHandler.h
	 */
	void computeOffDiagonalJacobian(TS &ts, Vec &localC, Mat &J) const;

	/**
	 * Compute the diagonal part of the Jacobian which is related to cluster reactions.
     * \see ISolverHandler.h
	 */
	void computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J);

	/**
	 * Get the position of the surface.
     * \see ISolverHandler.h
	 */
	int getSurfacePosition(int j = -1, int k = -1) const {
		return surfacePosition[j];
	}

	/**
	 * Set the position of the surface.
     * \see ISolverHandler.h
	 */
	void setSurfacePosition(int pos, int j = -1, int k = -1) {
		surfacePosition[j] = pos;
	}

}; //end class PetscSolver2DHandler

} /* end namespace xolotlSolver */
#endif
