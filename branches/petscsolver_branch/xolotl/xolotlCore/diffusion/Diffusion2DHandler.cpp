// Includes
#include "Diffusion2DHandler.h"

namespace xolotlCore {

void Diffusion2DHandler::computeDiffusion(PSIClusterReactionNetwork *network,
		double s, double **concVector, double *updatedConcOffset) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of diffusing cluster
	int nDiff = indexVector.size();

	// Get the number of degrees of freedom which is the size of the network
	int dof = reactants->size();

	// Loop on them
	for (int i = 0; i < nDiff; i++) {
		// Get the diffusing cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;

		// Get the initial concentrations
		double oldConc = concVector[0][index]; // middle
		double oldLeftConc = concVector[1][index]; // left
		double oldRightConc = concVector[2][index]; // right
		double oldBottomConc = concVector[3][index]; // bottom
		double oldTopConc = concVector[4][index]; // top

		// Use a simple midpoint stencil to compute the concentration
		double conc = cluster->getDiffusionCoefficient()
				* (-4.0 * oldConc + oldLeftConc + oldRightConc + oldBottomConc
						+ oldTopConc) * s;

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void Diffusion2DHandler::computePartialsForDiffusion(
		PSIClusterReactionNetwork *network,
		double s, double *val, int *indices) {
	// Get all the reactant
	auto reactants = network->getAll();
	// And the size of the network
	int size = reactants->size();
	// Get the number of diffusing cluster
	int nDiff = indexVector.size();

	// Loop on them
	for (int i = 0; i < nDiff; i++) {
		// Get the diffusing cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster->getDiffusionCoefficient();

		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[i] = index;

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, right, bottom, and top grid point
		val[i * 5] = -4.0 * diffCoeff * s; // middle
		val[(i * 5) + 1] = diffCoeff * s; // left
		val[(i * 5) + 2] = diffCoeff * s; // right
		val[(i * 5) + 3] = diffCoeff * s; // bottom
		val[(i * 5) + 4] = diffCoeff * s; // top
	}

	return;
}

}/* end namespace xolotlCore */
