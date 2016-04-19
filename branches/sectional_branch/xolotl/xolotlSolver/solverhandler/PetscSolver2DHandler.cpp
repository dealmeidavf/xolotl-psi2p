// Includes
#include <PetscSolver2DHandler.h>
#include <HDF5Utils.h>
#include <MathUtils.h>
#include <Constants.h>
#include <SuperCluster.h>

namespace xolotlSolver {

void PetscSolver2DHandler::createSolverContext(DM &da, int nx, double hx, int ny,
		double hy, int, double) {
	PetscErrorCode ierr;

	// Initialize the all reactants pointer
	allReactants = network->getAll();

	// Set the last temperature to 0
	lastTemperature = 0.0;

	// Reinitialize the connectivities in the network after updating the temperature
	// Get the temperature from the temperature handler
	auto temperature = temperatureHandler->getTemperature({0.0, 0.0, 0.0}, 0.0);

	// Set the temperature to compute all the rate constants
	if (!xolotlCore::equal(temperature, lastTemperature)) {
		// Update the temperature for all of the clusters
		int networkSize = network->size();
		for (int i = 0; i < networkSize; i++) {
			// This part will set the temperature in each reactant
			// and recompute the diffusion coefficient
			allReactants->at(i)->setTemperature(temperature);
		}
		for (int i = 0; i < networkSize; i++) {
			// Now that the diffusion coefficients of all the reactants
			// are updated, the reaction and dissociation rates can be
			// recomputed
			auto cluster = (xolotlCore::PSICluster *) allReactants->at(i);
			cluster->computeRateConstants();
		}
		lastTemperature = temperature;
	}

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->size() + 2 * network->getAll("Super").size();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_PERIODIC,
	DMDA_STENCIL_STAR, nx, ny, PETSC_DECIDE, PETSC_DECIDE, dof, 1, NULL, NULL, &da);
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"DMDACreate2d failed.");

	// Set the position of the surface
	for (int j = 0; j < ny; j++) {
		surfacePosition.push_back(0);
		if (movingSurface) surfacePosition[j] = (int) (nx * portion / 100.0);
	}

	// Generate the grid in the x direction
	generateGrid(nx, hx, surfacePosition[0]);

	for (int i = 0; i < grid.size(); i++) {
		std::cout << grid[i] << " ";
	}
	std::cout << std::endl;

	// Initialize the surface of the first advection handler corresponding to the
	// advection toward the surface (or a dummy one if it is deactivated)
	advectionHandlers[0]->setLocation(grid[surfacePosition[0]]);

	// Set the step size
	hY = hy;

	// Set the size of the partial derivatives vectors
	clusterPartials.resize(dof, 0.0);
	reactingPartialsForCluster.resize(dof, 0.0);

	/*  The only spatial coupling in the Jacobian is due to diffusion.
	 *  The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 *  the nonzero coupling between degrees of freedom at one point with degrees
	 *  of freedom on the adjacent point to the left or right. A 1 at i,j in the
	 *  ofill array indicates that the degree of freedom i at a point is coupled
	 *  to degree of freedom j at the adjacent point.
	 *  In this case ofill has only a few diagonal entries since the only spatial
	 *  coupling is regular diffusion.
	 */
	PetscInt *ofill, *dfill;
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &ofill);
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"PetscMalloc (ofill) failed.");
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &dfill);
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"PetscMalloc (dfill) failed.");
	ierr = PetscMemzero(ofill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"PetscMemzero (ofill) failed.");
	ierr = PetscMemzero(dfill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"PetscMemzero (dfill) failed.");

	// Fill ofill, the matrix of "off-diagonal" elements that represents diffusion
	diffusionHandler->initializeOFill(network, ofill);
	// Loop on the advection handlers to account the other "off-diagonal" elements
	for (int i = 0; i < advectionHandlers.size(); i++) {
		advectionHandlers[i]->initialize(network, ofill);
	}

	// Initialize the modified trap-mutation handler and the bubble bursting one here
	// because they add connectivity
	mutationHandler->initialize(surfacePosition[0], network, grid);
	burstingHandler->initialize(surfacePosition[0], network, grid);

	// Get the diagonal fill
	getDiagonalFill(dfill, dof * dof);

	// Load up the block fills
	ierr = DMDASetBlockFills(da, dfill, ofill);
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"DMDASetBlockFills failed.");

	// Free the temporary fill arrays
	ierr = PetscFree(ofill);
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"PetscFree (ofill) failed.");
	ierr = PetscFree(dfill);
	checkPetscError(ierr, "PetscSolver2DHandler::createSolverContext: "
			"PetscFree (dfill) failed.");

	return;
}

void PetscSolver2DHandler::initializeConcentration(DM &da, Vec &C) {
	PetscErrorCode ierr;

	// Pointer for the concentration vector
	PetscScalar ***concentrations;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver2DHandler::initializeConcentration: "
			"DMDAVecGetArrayDOF failed.");

	// Get the local boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr, "PetscSolver2DHandler::initializeConcentration: "
			"DMDAGetCorners failed.");

	// Get the last time step written in the HDF5 file
	int tempTimeStep = -2;
	bool hasConcentrations = xolotlCore::HDF5Utils::hasConcentrationGroup(
			networkName, tempTimeStep);

	// Get the actual surface position if concentrations were stored
	if (hasConcentrations) {
		auto surfaceIndices = xolotlCore::HDF5Utils::readSurface2D(networkName, tempTimeStep);

		// Set the actual surface positions
		for (int i = 0; i < surfaceIndices.size(); i++) {
			surfacePosition[i] = surfaceIndices[i];
		}
	}

	// Get the total size of the grid for the boundary conditions
	PetscInt Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr, "PetscSolver2DHandler::initializeConcentration: "
			"DMDAGetInfo failed.");

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(network, surfacePosition[0], grid);

	// Initialize the grid for the diffusion
	diffusionHandler->initializeDiffusionGrid(advectionHandlers, grid, My, hY);

	// Initialize the grid for the advection
	advectionHandlers[0]->initializeAdvectionGrid(advectionHandlers, grid, My, hY);

	// Pointer for the concentration vector at a specific grid point
	PetscScalar *concOffset;

	// Get all the super clusters
	auto superClusters = network->getAll("Super");
	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->size() + 2 * superClusters.size();

	// Get the single vacancy ID
	auto singleVacancyCluster = network->get(xolotlCore::vType, 1);
	int vacancyIndex = -1;
	if (singleVacancyCluster)
		vacancyIndex = singleVacancyCluster->getId() - 1;

	// Loop on all the grid points
	for (int j = ys; j < ys + ym; j++) {
		for (int i = xs; i < xs + xm; i++) {
			concOffset = concentrations[j][i];

			// Loop on all the clusters to initialize at 0.0
			for (int n = 0; n < dof; n++) {
				concOffset[n] = 0.0;
			}

			// Initialize the vacancy concentration
			if (i > surfacePosition[j] && i < Mx - 1 && vacancyIndex > 0) {
				concOffset[vacancyIndex] = initialVConc;
			}
		}
	}

	// If the concentration must be set from the HDF5 file
	if (hasConcentrations) {
		// Loop on the full grid
		for (int j = 0; j < My; j++) {
			for (int i = 0; i < Mx; i++) {
				// Read the concentrations from the HDF5 file
				auto concVector = xolotlCore::HDF5Utils::readGridPoint(networkName,
						tempTimeStep, i, j);

				// Change the concentration only if we are on the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
					concOffset = concentrations[j][i];
					// Loop on the concVector size
					for (unsigned int l = 0; l < concVector.size(); l++) {
						concOffset[(int) concVector.at(l).at(0)] =
								concVector.at(l).at(1);
					}
				}
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr, "PetscSolver2DHandler::initializeConcentration: "
			"DMDAVecRestoreArrayDOF failed.");

	return;
}

void PetscSolver2DHandler::updateConcentration(TS &ts, Vec &localC, Vec &F,
		PetscReal ftime) {
	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: TSGetDM failed.");

	// Get the total size of the grid in the x direction for the boundary conditions
	int xSize = grid.size();

	// Pointers to the PETSc arrays that start at the beginning (xs, ys) of the
	// local array!
	PetscScalar ***concs, ***updatedConcs;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAVecGetArrayDOFRead (localC) failed.");
	ierr = DMDAVecGetArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAVecGetArrayDOF (F) failed.");

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAGetCorners failed.");

	// Get the total size of the grid
	int Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAGetInfo failed.");

	// The following pointers are set to the first position in the conc or
	// updatedConc arrays that correspond to the beginning of the data for the
	// current grid point. They are accessed just like regular arrays.
	PetscScalar *concOffset, *updatedConcOffset;

	// Set some step size variable
	double sy = 1.0 / (hY * hY);

	// Declarations for variables used in the loop
	double flux;
	int fluxIndex = fluxHandler->getIncidentFluxClusterIndex(), reactantIndex;
	xolotlCore::PSICluster *cluster = NULL;
	double **concVector = new double*[5];
	std::vector<double> gridPosition = { 0.0, 0.0, 0.0 }, incidentFluxVector;
	xolotlCore::PSICluster * bubble;
	xolotlCore::SuperCluster * super;
	int index = 0;
	std::map<std::string, int> comp;
	int heComp = 0;
	double heConc = 0.0;
	double totalHeConc = 0.0;
	// Get all the HeV clusters in the network
	auto bubbles = network->getAll(xolotlCore::heVType);

	// Get all the super clusters
	auto superClusters = network->getAll("Super");
	// Degrees of freedom is the total number of clusters in the network
	const int networkSize = network->size();
	const int dof = networkSize + 2 * superClusters.size();

	// Loop over grid points computing ODE terms for each grid point
	for (int yj = 0; yj < My; yj++) {
		// Set the grid position
		gridPosition[1] = yj * hY;

		// Initialize the helium concentration near the surface
		heConc = 0.0;

		// Loop over grid points
		for (int xi = surfacePosition[yj] + 1; xi < Mx; xi++) {
			// We are only interested in the helium near the surface
			if (grid[xi] - grid[surfacePosition[yj]] > 2.0) continue;

			// Check if we are on the right processor
			if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
				// Get the concentrations at this grid point
				concOffset = concs[yj][xi];
				// Copy data into the PSIClusterReactionNetwork
				network->updateConcentrationsFromArray(concOffset);

				// Loop on the bubbles
				for (int i = 0; i < bubbles.size(); i++) {
					// Get the bubble, its id, and its helium composition
					bubble = (xolotlCore::PSICluster *) bubbles.at(i);
					index = bubble->getId() - 1;
					auto comp = bubble->getComposition();
					heComp = comp[xolotlCore::heType];

					// Sum the helium concentration
					heConc += concOffset[index] * (double) heComp * (grid[xi] - grid[xi-1]);
				}
				// Loop on the super clusters
				for (int i = 0; i < superClusters.size(); i++) {
					super = (xolotlCore::SuperCluster *) superClusters[i];

					// Sum the helium concentration
					heConc += super->getTotalHeliumConcentration() * (grid[xi] - grid[xi-1]);
				}
			}
		}

		// Share the concentration with all the processes
		totalHeConc = 0.0;
		MPI_Allreduce(&heConc, &totalHeConc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		// Set the disappearing rate in the modified TM handler
		mutationHandler->updateDisappearingRate(totalHeConc);

		// Skip if we are not on the right process
		if (yj < ys || yj >= ys + ym) continue;

		// Initialize the flux, advection, and trap-mutation handlers which depend
		// on the surface position at Y
		fluxHandler->initializeFluxHandler(network, surfacePosition[yj], grid);
		advectionHandlers[0]->setLocation(grid[surfacePosition[yj]]);
		mutationHandler->initializeIndex(surfacePosition[yj], network, grid);

		// Get the flux vector which can be different at each Y position
		incidentFluxVector = fluxHandler->getIncidentFluxVec(ftime, surfacePosition[yj]);

		for (int xi = xs; xi < xs + xm; xi++) {
			// Compute the old and new array offsets
			concOffset = concs[yj][xi];
			updatedConcOffset = updatedConcs[yj][xi];

			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi <= surfacePosition[yj] || xi == xSize - 1) {
				for (int i = 0; i < dof; i++) {
					updatedConcOffset[i] = 1.0 * concOffset[i];
				}

				continue;
			}

			// Set the grid position
			gridPosition[0] = grid[xi];

			// Fill the concVector with the pointer to the middle, left, right, bottom, and top grid points
			concVector[0] = concOffset; // middle
			concVector[1] = concs[yj][xi - 1]; // left
			concVector[2] = concs[yj][xi + 1]; // right
			concVector[3] = concs[yj - 1][xi]; // bottom
			concVector[4] = concs[yj + 1][xi]; // top

			// Get the temperature from the temperature handler
			auto temperature = temperatureHandler->getTemperature(gridPosition,
					ftime);

			// Update the network if the temperature changed
			if (!xolotlCore::equal(temperature, lastTemperature)) {
				network->setTemperature(temperature);
				// Update the modified trap-mutation rate that depends on the
				// network reaction rates
				mutationHandler->updateTrapMutationRate(network);
				lastTemperature = temperature;
			}

			// Copy data into the PSIClusterReactionNetwork so that it can
			// compute the fluxes properly. The network is only used to compute the
			// fluxes and hold the state data from the last time step. I'm reusing
			// it because it cuts down on memory significantly (about 400MB per
			// grid point) at the expense of being a little tricky to comprehend.
			network->updateConcentrationsFromArray(concOffset);

			// ----- Account for flux of incoming He of cluster size 1 -----
			updatedConcOffset[fluxIndex] += incidentFluxVector[xi - surfacePosition[yj]];

			// ---- Compute diffusion over the locally owned part of the grid -----
			diffusionHandler->computeDiffusion(network, concVector,
					updatedConcOffset, grid[xi] - grid[xi-1], grid[xi+1] - grid[xi], xi, sy, yj);

			// ---- Compute advection over the locally owned part of the grid -----
			for (int i = 0; i < advectionHandlers.size(); i++) {
				advectionHandlers[i]->computeAdvection(network, gridPosition,
						concVector, updatedConcOffset,
						grid[xi] - grid[xi-1], grid[xi+1] - grid[xi], xi, hY, yj);
			}

			// ----- Compute the modified trap-mutation over the locally owned part of the grid -----
			mutationHandler->computeTrapMutation(network, xi, concOffset,
					updatedConcOffset);

			// ----- Compute all of the new fluxes -----
			for (int i = 0; i < networkSize; i++) {
				cluster = (xolotlCore::PSICluster *) allReactants->at(i);
				// Compute the flux
				flux = cluster->getTotalFlux();
				// Update the concentration of the cluster
				reactantIndex = cluster->getId() - 1;
				updatedConcOffset[reactantIndex] += flux;
			}

			// ---- Moments ----
			for (int i = 0; i < superClusters.size(); i++) {
				cluster = (xolotlCore::PSICluster *) superClusters[i];
				// Compute the helium momentum flux
				flux = cluster->getHeMomentFlux();
				// Update the concentration of the cluster
				reactantIndex = cluster->getHeMomentumId() - 1;
				updatedConcOffset[reactantIndex] += flux;

				// Compute the vacancy momentum flux
				flux = cluster->getVMomentFlux();
				// Update the concentration of the cluster
				reactantIndex = cluster->getVMomentumId() - 1;
				updatedConcOffset[reactantIndex] += flux;
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOFRead (localC) failed.");
	ierr = DMDAVecRestoreArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMDAVecRestoreArrayDOF (F) failed.");
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr, "PetscSolver2DHandler::updateConcentration: "
			"DMRestoreLocalVector failed.");

	// Clear memory
	delete [] concVector;

	return;
}

void PetscSolver2DHandler::computeOffDiagonalJacobian(TS &ts, Vec &localC, Mat &J) const {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver2DHandler::computeOffDiagonalJacobian: "
			"TSGetDM failed.");

	// Get the total size of the grid in the x direction for the boundary conditions
	int xSize = grid.size();

	// Setup some step size variables
	double sy = 1.0 / (hY * hY);

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr, "PetscSolver2DHandler::computeOffDiagonalJacobian: "
			"DMDAGetCorners failed.");

	// Get the total number of diffusing clusters
	const int nDiff = diffusionHandler->getNumberOfDiffusing();

	// Get the total number of advecting clusters
	int nAdvec = 0;
	for (int l = 0; l < advectionHandlers.size(); l++) {
		int n = advectionHandlers[l]->getNumberOfAdvecting();
		if (n > nAdvec) nAdvec = n;
	}

	// Arguments for MatSetValuesStencil called below
	MatStencil row, cols[5];
	PetscScalar diffVals[5 * nDiff];
	PetscInt diffIndices[nDiff];
	PetscScalar advecVals[2 * nAdvec];
	PetscInt advecIndices[nAdvec];
	std::vector<double> gridPosition = { 0.0, 0.0, 0.0 };

	/*
	 Loop over grid points computing Jacobian terms for diffusion and advection
	 at each grid point
	 */
	for (int yj = ys; yj < ys + ym; yj++) {
		// Set the grid position
		gridPosition[1] = yj * hY;

		// Initialize the advection handler which depends
		// on the surface position at Y
		advectionHandlers[0]->setLocation(grid[surfacePosition[yj]]);

		for (int xi = xs; xi < xs + xm; xi++) {
			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi <= surfacePosition[yj] || xi == xSize - 1) continue;

			// Set the grid position
			gridPosition[0] = grid[xi];

			// Get the partial derivatives for the diffusion
			diffusionHandler->computePartialsForDiffusion(network, diffVals, diffIndices,
					grid[xi] - grid[xi-1], grid[xi+1] - grid[xi], xi, sy, yj);

			// Loop on the number of diffusion cluster to set the values in the Jacobian
			for (int i = 0; i < nDiff; i++) {
				// Set grid coordinate and component number for the row
				row.i = xi;
				row.j = yj;
				row.c = diffIndices[i];

				// Set grid coordinates and component numbers for the columns
				// corresponding to the middle, left, right, bottom, and top grid points
				cols[0].i = xi; // middle
				cols[0].j = yj;
				cols[0].c = diffIndices[i];
				cols[1].i = xi - 1; // left
				cols[1].j = yj;
				cols[1].c = diffIndices[i];
				cols[2].i = xi + 1; // right
				cols[2].j = yj;
				cols[2].c = diffIndices[i];
				cols[3].i = xi; // bottom
				cols[3].j = yj - 1;
				cols[3].c = diffIndices[i];
				cols[4].i = xi; // top
				cols[4].j = yj + 1;
				cols[4].c = diffIndices[i];

				ierr = MatSetValuesStencil(J, 1, &row, 5, cols, diffVals + (5 * i), ADD_VALUES);
				checkPetscError(ierr, "PetscSolver2DHandler::computeOffDiagonalJacobian: "
						"MatSetValuesStencil (diffusion) failed.");
			}

			// Get the partial derivatives for the advection
			for (int l = 0; l < advectionHandlers.size(); l++) {
				advectionHandlers[l]->computePartialsForAdvection(network, advecVals,
						advecIndices, gridPosition,
						grid[xi] - grid[xi-1], grid[xi+1] - grid[xi], xi, hY, yj);

				// Get the stencil indices to know where to put the partial derivatives in the Jacobian
				auto advecStencil = advectionHandlers[l]->getStencilForAdvection(gridPosition);

				// Get the number of advecting clusters
				nAdvec = advectionHandlers[l]->getNumberOfAdvecting();

				// Loop on the number of advecting cluster to set the values in the Jacobian
				for (int i = 0; i < nAdvec; i++) {
					// Set grid coordinate and component number for the row
					row.i = xi;
					row.j = yj;
					row.c = advecIndices[i];

					// If we are on the sink, the partial derivatives are not the same
					// Both sides are giving their concentrations to the center
					if (advectionHandlers[l]->isPointOnSink(gridPosition)) {
						cols[0].i = xi - advecStencil[0]; // left?
						cols[0].j = yj - advecStencil[1]; // bottom?
						cols[0].c = advecIndices[i];
						cols[1].i = xi + advecStencil[0]; // right?
						cols[1].j = yj + advecStencil[1]; // top?
						cols[1].c = advecIndices[i];
					}
					else {
						// Set grid coordinates and component numbers for the columns
						// corresponding to the middle and other grid points
						cols[0].i = xi; // middle
						cols[0].j = yj;
						cols[0].c = advecIndices[i];
						cols[1].i = xi + advecStencil[0]; // left or right?
						cols[1].j = yj + advecStencil[1]; // bottom or top?
						cols[1].c = advecIndices[i];
					}

					// Update the matrix
					ierr = MatSetValuesStencil(J, 1, &row, 2, cols, advecVals + (2 * i), ADD_VALUES);
					checkPetscError(ierr, "PetscSolver2DHandler::computeOffDiagonalJacobian: "
							"MatSetValuesStencil (advection) failed.");
				}
			}
		}
	}

	return;
}

void PetscSolver2DHandler::computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"TSGetDM failed.");

	// Get the total size of the grid in the x direction for the boundary conditions
	int xSize = grid.size();

	// Get pointers to vector data
	PetscScalar ***concs;
	ierr = DMDAVecGetArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"DMDAVecGetArrayDOFRead failed.");

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"DMDAGetCorners failed.");

	// Get the total size of the grid
	int Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"DMDAGetInfo failed.");

	// Get all the super clusters
	auto superClusters = network->getAll("Super");
	// Degrees of freedom is the total number of clusters in the network
	const int networkSize = network->size();
	const int dof = networkSize + 2 * superClusters.size();

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset;

	// Arguments for MatSetValuesStencil called below
	MatStencil rowId;
	MatStencil colIds[dof];
	int pdColIdsVectorSize = 0;

	// Store the total number of He clusters in the network for the
	// modified trap-mutation
	int nHelium = network->getAll(xolotlCore::heType).size();

	// Declarations for variables used in the loop
	int reactantIndex;
	xolotlCore::PSICluster * bubble;
	xolotlCore::SuperCluster * super;
	int index = 0;
	std::map<std::string, int> comp;
	int heComp = 0;
	double heConc = 0.0;
	double totalHeConc = 0.0;
	// Get all the HeV clusters in the network
	auto bubbles = network->getAll(xolotlCore::heVType);

	// Loop over the grid points
	for (int yj = 0; yj < My; yj++) {
		// Initialize the helium concentration near the surface
		heConc = 0.0;

		// Loop over grid points
		for (int xi = surfacePosition[yj] + 1; xi < Mx; xi++) {
			// We are only interested in the helium near the surface
			if (grid[xi] - grid[surfacePosition[yj]] > 2.0) continue;

			// Check if we are on the right processor
			if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
				// Get the concentrations at this grid point
				concOffset = concs[yj][xi];
				// Copy data into the PSIClusterReactionNetwork
				network->updateConcentrationsFromArray(concOffset);

				// Loop on the bubbles
				for (int i = 0; i < bubbles.size(); i++) {
					// Get the bubble, its id, and its helium composition
					bubble = (xolotlCore::PSICluster *) bubbles.at(i);
					index = bubble->getId() - 1;
					auto comp = bubble->getComposition();
					heComp = comp[xolotlCore::heType];

					// Sum the helium concentration
					heConc += concOffset[index] * (double) heComp * (grid[xi] - grid[xi-1]);
				}
				// Loop on the super clusters
				for (int i = 0; i < superClusters.size(); i++) {
					super = (xolotlCore::SuperCluster *) superClusters[i];

					// Sum the helium concentration
					heConc += super->getTotalHeliumConcentration() * (grid[xi] - grid[xi-1]);
				}
			}
		}

		// Share the concentration with all the processes
		totalHeConc = 0.0;
		MPI_Allreduce(&heConc, &totalHeConc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		// Set the disappearing rate in the modified TM handler
		mutationHandler->updateDisappearingRate(totalHeConc);

		// Skip if we are not on the right process
		if (yj < ys || yj >= ys + ym) continue;

		// Initialize the trap-mutation handler which depends
		// on the surface position at Y
		mutationHandler->initializeIndex(surfacePosition[yj], network, grid);

		for (int xi = xs; xi < xs + xm; xi++) {
			// Boundary conditions
			// Everything to the left of the surface is empty
			if (xi <= surfacePosition[yj] || xi == xSize - 1) continue;

			// Copy data into the PSIClusterReactionNetwork so that it can
			// compute the new concentrations.
			concOffset = concs[yj][xi];
			network->updateConcentrationsFromArray(concOffset);

			// Update the column in the Jacobian that represents each normal reactant
			for (int i = 0; i < networkSize - superClusters.size(); i++) {
				auto reactant = allReactants->at(i);
				// Get the reactant index
				reactantIndex = reactant->getId() - 1;

				// Set grid coordinate and component number for the row
				rowId.i = xi;
				rowId.j = yj;
				rowId.c = reactantIndex;

				// Get the partial derivatives
				reactant->getPartialDerivatives(clusterPartials);
				// Get the list of column ids from the map
				auto pdColIdsVector = dFillMap.at(reactantIndex);
				// Number of partial derivatives
				pdColIdsVectorSize = pdColIdsVector.size();
				// Loop over the list of column ids
				for (int j = 0; j < pdColIdsVectorSize; j++) {
					// Set grid coordinate and component number for a column in the list
					colIds[j].i = xi;
					colIds[j].j = yj;
					colIds[j].c = pdColIdsVector[j];
					// Get the partial derivative from the array of all of the partials
					reactingPartialsForCluster[j] =
							clusterPartials[pdColIdsVector[j]];
					// Reset the cluster partial value to zero. This is much faster
					// than using memset.
					clusterPartials[pdColIdsVector[j]] = 0.0;
				}
				// Update the matrix
				ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize, colIds,
						reactingPartialsForCluster.data(), ADD_VALUES);
				checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
						"MatSetValuesStencil (reactions) failed.");
			}

			// Update the column in the Jacobian that represents the moment for the super clusters
			for (int i = 0; i < superClusters.size(); i++) {
				auto reactant = (xolotlCore::SuperCluster *) superClusters[i];

				// Get the super cluster index
				reactantIndex = reactant->getId() - 1;
				// Set grid coordinate and component number for the row
				rowId.i = xi;
				rowId.c = reactantIndex;

				// Get the partial derivatives
				reactant->getPartialDerivatives(clusterPartials);
				// Get the list of column ids from the map
				auto pdColIdsVector = dFillMap.at(reactantIndex);
				// Number of partial derivatives
				pdColIdsVectorSize = pdColIdsVector.size();

				// Loop over the list of column ids
				for (int j = 0; j < pdColIdsVectorSize; j++) {
					// Set grid coordinate and component number for a column in the list
					colIds[j].i = xi;
					colIds[j].c = pdColIdsVector[j];
					// Get the partial derivative from the array of all of the partials
					reactingPartialsForCluster[j] =
							clusterPartials[pdColIdsVector[j]];

					// Reset the cluster partial value to zero. This is much faster
					// than using memset.
					clusterPartials[pdColIdsVector[j]] = 0.0;
				}
				// Update the matrix
				ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize,
						colIds, reactingPartialsForCluster.data(), ADD_VALUES);
				checkPetscError(ierr, "PetscSolver1DHandler::computeDiagonalJacobian: MatSetValuesStencil for super cluster failed.");

				// Get the helium momentum index
				reactantIndex = reactant->getHeMomentumId() - 1;
				// Set component number for the row
				rowId.c = reactantIndex;

				// Get the partial derivatives
				reactant->getHeMomentPartialDerivatives(clusterPartials);
				// Get the list of column ids from the map
				pdColIdsVector = dFillMap.at(reactantIndex);
				// Number of partial derivatives
				pdColIdsVectorSize = pdColIdsVector.size();

				// Loop over the list of column ids
				for (int j = 0; j < pdColIdsVectorSize; j++) {
					// Set grid coordinate and component number for a column in the list
					colIds[j].i = xi;
					colIds[j].c = pdColIdsVector[j];
					// Get the partial derivative from the array of all of the partials
					reactingPartialsForCluster[j] =
							clusterPartials[pdColIdsVector[j]];

					// Reset the cluster partial value to zero. This is much faster
					// than using memset.
					clusterPartials[pdColIdsVector[j]] = 0.0;
				}
				// Update the matrix
				ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize,
						colIds, reactingPartialsForCluster.data(), ADD_VALUES);
				checkPetscError(ierr, "PetscSolver1DHandler::computeDiagonalJacobian: MatSetValuesStencil for helium momentum failed.");

				// Get the vacancy momentum index
				reactantIndex = reactant->getVMomentumId() - 1;
				// Set component number for the row
				rowId.c = reactantIndex;

				// Get the partial derivatives
				reactant->getVMomentPartialDerivatives(clusterPartials);
				// Get the list of column ids from the map
				pdColIdsVector = dFillMap.at(reactantIndex);
				// Number of partial derivatives
				pdColIdsVectorSize = pdColIdsVector.size();

				// Loop over the list of column ids
				for (int j = 0; j < pdColIdsVectorSize; j++) {
					// Set grid coordinate and component number for a column in the list
					colIds[j].i = xi;
					colIds[j].c = pdColIdsVector[j];
					// Get the partial derivative from the array of all of the partials
					reactingPartialsForCluster[j] =
							clusterPartials[pdColIdsVector[j]];

					// Reset the cluster partial value to zero. This is much faster
					// than using memset.
					clusterPartials[pdColIdsVector[j]] = 0.0;
				}
				// Update the matrix
				ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize,
						colIds, reactingPartialsForCluster.data(), ADD_VALUES);
				checkPetscError(ierr, "PetscSolver1DHandler::computeDiagonalJacobian: MatSetValuesStencil for vacancy momentum failed.");
			}

			// ----- Take care of the bubble bursting for all the bubbles -----

			// Store the total number of HeV bubbles in the network for the
			// bubble bursting
			int nBubble = burstingHandler->getNBursting(network, xi);

			// Arguments for MatSetValuesStencil called below
			MatStencil row, col;
			PetscScalar burstingVals[nBubble];
			PetscInt burstingIndices[2 * nBubble];

			// Compute the partial derivative from bubble bursting at this grid point
			int nBursting = burstingHandler->computePartialsForBursting(network,
					burstingVals, burstingIndices, xi);

			// Loop on the number of bubbles bursting to set the values
			// in the Jacobian
			for (int i = 0; i < nBursting; i++) {
				// Set grid coordinate and component number for the row and column
				// corresponding to the bursting cluster
				row.i = xi;
				row.c = burstingIndices[2 * i];
				col.i = xi;
				col.c = burstingIndices[(2 * i) + 1];

				ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
						burstingVals + i, ADD_VALUES);
				checkPetscError(ierr, "PetscSolver1DHandler::computeDiagonalJacobian: "
						"MatSetValuesStencil (bursting) failed.");
			}

			// ----- Take care of the modified trap-mutation for all the reactants -----

			// Arguments for MatSetValuesStencil called below
			PetscScalar mutationVals[3 * nHelium];
			PetscInt mutationIndices[3 * nHelium];

			// Compute the partial derivative from modified trap-mutation at this grid point
			int nMutating = mutationHandler->computePartialsForTrapMutation(network,
					mutationVals, mutationIndices, xi);

			// Loop on the number of helium undergoing trap-mutation to set the values
			// in the Jacobian
			for (int i = 0; i < nMutating; i++) {
				// Set grid coordinate and component number for the row and column
				// corresponding to the helium cluster
				row.i = xi;
				row.j = yj;
				row.c = mutationIndices[3 * i];
				col.i = xi;
				col.j = yj;
				col.c = mutationIndices[3 * i];

				ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
						mutationVals + (3 * i), ADD_VALUES);
				checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
						"MatSetValuesStencil (He trap-mutation) failed.");

				// Set component number for the row
				// corresponding to the HeV cluster created through trap-mutation
				row.c = mutationIndices[(3 * i) + 1];

				ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
						mutationVals + (3 * i) + 1, ADD_VALUES);
				checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
						"MatSetValuesStencil (HeV trap-mutation) failed.");

				// Set component number for the row
				// corresponding to the interstitial created through trap-mutation
				row.c = mutationIndices[(3 * i) + 2];

				ierr = MatSetValuesStencil(J, 1, &row, 1, &col,
						mutationVals + (3 * i) + 2, ADD_VALUES);
				checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
						"MatSetValuesStencil (I trap-mutation) failed.");
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOFRead(da, localC, &concs);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"DMDAVecRestoreArrayDOFRead failed.");
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr, "PetscSolver2DHandler::computeDiagonalJacobian: "
			"DMRestoreLocalVector failed.");

	return;
}

} /* end namespace xolotlSolver */