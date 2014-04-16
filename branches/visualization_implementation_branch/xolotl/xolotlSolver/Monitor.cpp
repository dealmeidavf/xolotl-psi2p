// Includes
#include "../xolotlViz/plot/Plot.h"
#include "../xolotlViz/dataprovider/Point.h"
#include "../xolotlViz/plot/ScatterPlot.h"
#include "../xolotlViz/plot/SeriesPlot.h"
#include "../xolotlViz/dataprovider/CvsXDataProvider.h"
#include "../xolotlViz/labelprovider/LabelProvider.h"
#include "PetscSolver.h"
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <vector>
#include <memory>
#include "../xolotlPerf/HandlerRegistryFactory.h"

namespace xolotlSolver {

/* ----- Error Handling Code ----- */

/**
 * This operation checks a Petsc error code and converts it to a bool.
 * @param errorCode The Petsc error code.
 * @return True if everything is OK, false otherwise.
 */
static inline bool checkPetscError(PetscErrorCode errorCode) {
	CHKERRQ(errorCode);
}

//! The pointer to the plot that will be used to visualize the data.
std::shared_ptr<xolotlViz::SeriesPlot> plot;

//! The double that will store the accumulation of helium flux.
double heliumFluence = 0.;

/**
 * This is a monitoring method that will compute the total helium fluence
 */
static PetscErrorCode heliumRetention(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int size = PetscSolver::getNetwork()->size();
	PetscErrorCode ierr;
	PetscInt xs, xm, Mx;

	PetscFunctionBeginUser;

	// Variable to represent the real, or current, time
	PetscReal realTime;
	ierr = TSGetTime(ts, &realTime);

	// Get the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = PetscSolver::getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Setup step size variable
	double hx = 8.0 / (PetscReal) (Mx - 1);

	// Get the helium cluster
	auto heCluster = std::dynamic_pointer_cast<PSICluster>(
			PetscSolver::getNetwork()->get("He", 1));

	// Exit if there is no helium cluster in the network
	if (!heCluster)
		PetscFunctionReturn(0);

	// Get the composition of the cluster
	auto thisComp = heCluster->getComposition();
	std::vector<int> compVec = { thisComp["He"], thisComp["V"], thisComp["I"] };

	// Define the reactant ID
	int reactantIndex = heCluster->getId() - 1;

	// Loop on the grid
	for (int xi = xs; xi < xs + xm; xi++) {
		// Actual position in nm
		double x = xi * hx;

		// Vector representing the position at which the flux will be calculated
		// Currently we are only in 1D
		std::vector<double> gridPosition = { 0, x, 0 };

		// Calculate the incident flux
		auto incidentFlux = fluxHandler->getIncidentFlux(compVec, gridPosition,
				realTime);

		// And add it to the fluence
		heliumFluence += 10000. * incidentFlux * time;
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will create some plots.
 */
static PetscErrorCode monitorSolve(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int size = PetscSolver::getNetwork()->size();
	// The array of cluster names
	std::vector<std::string> names(size);
	PetscErrorCode ierr;
	PetscReal *solutionArray, *gridPointSolution, x, hx;
	Vec localSolution;
	PetscInt xs, xm, Mx;
	int xi, i;

	PetscFunctionBeginUser;

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Fill the array of clusters name because the Id is not the same as
	// reactants->at(i)
	auto reactants = PetscSolver::getNetwork()->getAll();
	std::shared_ptr<PSICluster> cluster;
	for (int i = 0; i < size; i++) {

		// Get the cluster from the list, its id and composition
		cluster = std::dynamic_pointer_cast<PSICluster>(reactants->at(i));
		int id = cluster->getId() - 1;
		auto composition = cluster->getComposition();

		// Create the name
		std::stringstream name;
		name << (cluster->getName()).c_str() << "(" << composition["He"] << ","
				<< composition["V"] << "," << composition["I"] << ") ";

		// Push the header entry on the array
		name >> names[id];
	}

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the local vector, which is capital when running in parallel,
	// and put it into solutionArray
	ierr = DMGetLocalVector(da, &localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMDAVecGetArray(da, localSolution, &solutionArray);
	checkPetscError(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);
	// Setup some step size variables
	hx = 8.0 / (PetscReal) (Mx - 1);

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();
	auto myPointsBis = std::make_shared<std::vector<xolotlViz::Point> >();
	auto myPointsTer = std::make_shared<std::vector<xolotlViz::Point> >();
	auto myPointsQua = std::make_shared<std::vector<xolotlViz::Point> >();
	auto myPointsCin = std::make_shared<std::vector<xolotlViz::Point> >();

	// Choice of the cluster to be plotted
	int iCluster = 7;

	// Print the solution data
	for (xi = xs; xi < xs + xm; xi++) {
		// Dump x
		x = xi * hx;
		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray + size * xi;
		// Update the concentrations in the network to have physics results
		// (non negative)
		PetscSolver::getNetwork()->updateConcentrationsFromArray(
				gridPointSolution);
		// Get the concentrations from the network
		double concentrations[size];
		double * concentration = &concentrations[0];
		PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

		// Create a Point with the concentration[iCluster] as the value
		// and add it to myPoints
		xolotlViz::Point aPoint;
		aPoint.value = concentration[2]; // He
		aPoint.t = time;
		aPoint.x = x;
		myPoints->push_back(aPoint);
		aPoint.value = concentration[11]; // V
		myPointsBis->push_back(aPoint);
		aPoint.value = concentration[12]; // He1V1
		myPointsTer->push_back(aPoint);
		aPoint.value = concentration[13]; // He2V1
		myPointsQua->push_back(aPoint);
		aPoint.value = concentration[29]; // He1V2
		myPointsCin->push_back(aPoint);
	}

	// Get the data provider and give it the points
	plot->getDataProvider(0)->setPoints(myPoints);
	plot->getDataProvider(1)->setPoints(myPointsBis);
	plot->getDataProvider(2)->setPoints(myPointsTer);
	plot->getDataProvider(3)->setPoints(myPointsQua);
	plot->getDataProvider(4)->setPoints(myPointsCin);

	// Change the title of the plot
	std::stringstream title;
	title << "logTS" << timestep << "_" << procId << ".pnm";
	plot->plotLabelProvider->titleLabel = title.str();

	// Render
	plot->render();
	// Restore the array
	VecRestoreArray(solution, &solutionArray);

	PetscFunctionReturn(0);
}

/**
 * This operation sets up a monitor that will call monitorSolve
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetscMonitor(TS ts) {
	PetscErrorCode ierr;
	// Flags to launch the monitors or not
	PetscBool flagPlot, flagRetention;

	PetscFunctionBeginUser;

	// Check the option -plot
	ierr = PetscOptionsHasName(NULL, "-plot", &flagPlot);
	checkPetscError(ierr);

	// Check the option -helium_retention
	ierr = PetscOptionsHasName(NULL, "-helium_retention", &flagRetention);
	checkPetscError(ierr);

	// Don't do anything if no option is set
	if (!flagPlot && !flagRetention)
		PetscFunctionReturn(0);

	// Set the monitor to save plots
	if (flagPlot) {
		// Create a ScatterPlot
		plot = std::make_shared<xolotlViz::SeriesPlot>();

		// Create and set the label provider
		auto labelProvider = std::make_shared<xolotlViz::LabelProvider>();
		labelProvider->axis1Label = "x Position on the Grid";
		labelProvider->axis2Label = "Concentration";

		// Give it to the plot
		plot->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared<xolotlViz::CvsXDataProvider>();
		auto dataProviderBis = std::make_shared<xolotlViz::CvsXDataProvider>();
		auto dataProviderTer = std::make_shared<xolotlViz::CvsXDataProvider>();
		auto dataProviderQua = std::make_shared<xolotlViz::CvsXDataProvider>();
		auto dataProviderCin = std::make_shared<xolotlViz::CvsXDataProvider>();

		// Give it to the plot
		plot->addDataProvider(dataProvider);
		plot->addDataProvider(dataProviderBis);
		plot->addDataProvider(dataProviderTer);
		plot->addDataProvider(dataProviderQua);
		plot->addDataProvider(dataProviderCin);

		// monitorSolve will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSolve, NULL, NULL);
		checkPetscError(ierr);
	}

	// Set the monitor to compute the helium fluence for the retention calculation
	if (flagRetention) {
		// monitorSolve will be called at each timestep
		ierr = TSMonitorSet(ts, heliumRetention, NULL, NULL);
		checkPetscError(ierr);
	}

	PetscFunctionReturn(0);
}

/**
 * This operation computes and prints the helium retention.
 * @param ts The time stepper.
 * @param C The vector of solution.
 */
void computeRetention(TS ts, Vec C) {
	PetscErrorCode ierr;
	// Network size
	const int size = PetscSolver::getNetwork()->size();

	// Get the helium cluster
	auto heCluster = std::dynamic_pointer_cast<PSICluster>(
			PetscSolver::getNetwork()->get("He", 1));

	if (!heCluster) {
		throw std::string(
				"PetscSolver Exception: Cannot compute the retention because there is no helium1 in the network.");
		return;
	}

	// Keep the ID of the helium
	int reactantIndex = heCluster->getId() - 1;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the array of concentration
	PetscReal *solutionArray;
	ierr = DMDAVecGetArray(da, C, &solutionArray);
	checkPetscError(ierr);

	//Get local grid boundaries
	PetscInt xs, xm;
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Store the concentration over the grid
	double heConcentration = 0;

	// Loop on the grid
	for (int xi = xs; xi < xs + xm; xi++) {
		// Get the pointer to the beginning of the solution data for this grid point
		PetscReal *gridPointSolution;
		gridPointSolution = solutionArray + size * xi;
		// Update the concentrations in the network to have physics results
		// (non negative)
		PetscSolver::getNetwork()->updateConcentrationsFromArray(
				gridPointSolution);
		// Get the concentrations from the network
		double concentrations[size];
		double * concentration = &concentrations[0];
		PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

		// Add the current concentration
		heConcentration += concentration[reactantIndex];
	}

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Get the current process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Master process
	if (procId == 0) {
		// Loop on all the other processes
		for (int i = 1; i < worldSize; i++) {
			double otherConcentration = 0.;
			double otherFluence = 0.;
			// Receive the value from the other processes
			MPI_Recv(&otherConcentration, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
			MPI_Recv(&otherFluence, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Add them to the master one
			heConcentration += otherConcentration;
			heliumFluence += otherFluence;
		}

		// Print the result
		std::cout << "Helium retention = "
				<< 100. * heConcentration / heliumFluence << " %" << std::endl;
	}

	else {
		// Send the value of the timer to the master process
		MPI_Send(&heConcentration, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&heliumFluence, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	return;
}

}

/* end namespace xolotlSolver */