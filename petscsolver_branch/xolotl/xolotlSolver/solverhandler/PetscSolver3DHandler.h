#ifndef PETSCSOLVER3DHANDLER_H
#define PETSCSOLVER3DHANDLER_H

// Includes
#include "PetscSolverHandler.h"

namespace xolotlSolver {

/**
 * This class is a subclass of PetscSolverHandler and implement all the methods needed
 * to solve the ADR equations in 3D using PETSc from Argonne National Laboratory.
 */
class PetscSolver3DHandler: public PetscSolverHandler {
private:
	//! The position of the surface
	std::vector< std::vector<int> > surfacePosition;

	/**
	 * Get the mean position of the surface.
     *
     * @return The mean position of the surface
	 */
	int getMeanSurfacePosition() {
		int sizeY = surfacePosition.size();
		int sizeZ = surfacePosition.at(0).size();

		int mean = 0;
		// Compute the mean
		for (int j = 0; j < sizeY; j++) {
			for (int k = 0; k < sizeZ; k++) {
				mean += surfacePosition[j][k];
			}
		}

		return mean / (sizeY * sizeZ);
	}
public:

	//! The Constructor
	PetscSolver3DHandler() {}

	//! The Destructor
	~PetscSolver3DHandler() {}

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
	void initializeConcentration(DM &da, Vec &C);

	/**
	 * Compute the new concentrations for the RHS function given an initial
	 * vector of concentrations. Apply the diffusion, advection and all the reactions.
     * \see ISolverHandler.h
	 */
	void updateConcentration(TS &ts, Vec &localC, Vec &F, PetscReal ftime);

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
		return surfacePosition[j][k];
	}

	/**
	 * Set the position of the surface.
     * \see ISolverHandler.h
	 */
	void setSurfacePosition(int pos, int j = -1, int k = -1) {
		surfacePosition[j][k] = pos;
	}

}; //end class PetscSolver3DHandler

} /* end namespace xolotlSolver */
#endif
