#include "FluxHandler.h"
#include <xolotlPerf.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <mpi.h>

namespace xolotlCore {

FluxHandler::FluxHandler() :
		stepSize(0.0),
		heFluence(0.0),
		usingMaxHeFluence(false),
		maxHeFluence(std::numeric_limits<double>::max()),
		heFlux(1.0),
		incidentFluxZero(false),
		useTimeProfile(false),
		normFactor(0.0){

}

void FluxHandler::initializeFluxHandler(int numGridpoints, double step, int surfacePos) {

	// Set the step size
	stepSize = step;

	normFactor = 0.0;
	for (int i = surfacePos + 1; i < numGridpoints - 1; i++) {
		double x = (double) (i - surfacePos) * stepSize;

		normFactor += FitFunction(x) * stepSize;
	}

	// Factor the incident flux will be multiplied by
	double heFluxNormalized = heFlux / normFactor;

	// The first values should always be 0.0 because it is on the left side of the surface
	for (int i = 0; i <= surfacePos; i++) {
		incidentFluxVec.push_back(0.0);
	}

	// End before the last grid point because of the boundary conditions
	for (int i = surfacePos + 1; i < numGridpoints - 1; i++) {
		double x = (double) (i - surfacePos) * stepSize;

		auto incidentFlux = heFluxNormalized * FitFunction(x);

		incidentFluxVec.push_back(incidentFlux);
	}

	// The last value should always be 0.0 because of boundary conditions
	incidentFluxVec.push_back(0.0);

	return;
}

void FluxHandler::recomputeFluxHandler(int surfacePos) {
	// Factor the incident flux will be multiplied by
	double heFluxNormalized = heFlux / normFactor;

	// Get the number of grid points
	int numGridPoints = incidentFluxVec.size();

	// Clear the flux vector
	incidentFluxVec.clear();

	// The first values should always be 0.0 because it is on the left side of the surface
	for (int i = 0; i <= surfacePos; i++) {
		incidentFluxVec.push_back(0.0);
	}

	// Starts a i = 1 because the first value was already put in the vector
	for (int i = surfacePos + 1; i < numGridPoints - 1; i++) {
		double x = (double) (i - surfacePos) * stepSize;

		auto incidentFlux = heFluxNormalized * FitFunction(x);

		incidentFluxVec.push_back(incidentFlux);
	}

	// The last value should always be 0.0 because of boundary conditions
	incidentFluxVec.push_back(0.0);

	return;
}

void FluxHandler::initializeTimeProfile(std::string fileName) {
	// Set use time profile to true
	useTimeProfile = true;

	// Open file dataFile.dat containing the time and amplitude
	std::ifstream inputFile(fileName.c_str());
	std::string line;

	while (getline(inputFile, line)) {
		if (!line.length() || line[0] == '#')
			continue;
		double xamp = 0.0, yamp = 0.0;
		sscanf(line.c_str(), "%lf %lf", &xamp, &yamp);
		time.push_back(xamp);
		amplitude.push_back(yamp);
	}

	return;
}

double FluxHandler::getAmplitude(double currentTime) const {

	double f = 0.0;

	// if x is smaller than or equal to xi[0]
	if (currentTime <= time[0])
		return f = amplitude[0];

	// if x is greater than or equal to xi[n-1]
	if (currentTime >= time[time.size() - 1])
		return f = amplitude[time.size() - 1];

	// loop to determine the interval x falls in, ie x[k] < x < x[k+1]
	for (int k = 0; k < time.size() - 1; k++) {
		if (currentTime < time[k]) continue;
		if (currentTime > time[k + 1]) continue;

		f = amplitude[k]
				+ (amplitude[k + 1] - amplitude[k]) * (currentTime - time[k])
						/ (time[k + 1] - time[k]);
		break;
	}

	return f;
}

double FluxHandler::getIncidentFlux(std::vector<int> compositionVec,
		std::vector<double> position, double currentTime, int surfacePos) {

	// Recompute the flux vector if a time profile is used
	if (useTimeProfile) {
		heFlux = getAmplitude(currentTime);
		recomputeFluxHandler(surfacePos);
	}

	// Get the index number from the position
	int i = position[0] / stepSize;

	// Return the corresponding value
	return incidentFluxVec[i];
}

std::vector<double> FluxHandler::getIncidentFluxVec(double currentTime, int surfacePos) {

	// Recompute the flux vector if a time profile is used
	if (useTimeProfile) {
		heFlux = getAmplitude(currentTime);
		recomputeFluxHandler(surfacePos);
	}

	return incidentFluxVec;
}

void FluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
		std::vector<int> position, double time, double outgoingFlux) {

	return;
}

void FluxHandler::incrementHeFluence(double dt) {

	if (heFluence < maxHeFluence)
	{
		heFluence += heFlux * dt;
	}

	else if (!incidentFluxZero) {
		for (int i = 0; i < incidentFluxVec.size(); i++)
			incidentFluxVec[i] = 0.0;
		incidentFluxZero = true;
	}

	return;
}

double FluxHandler::getHeFluence() const {
	return heFluence;
}

void FluxHandler::setMaxHeFluence(double fluence) {
	usingMaxHeFluence = true;
	maxHeFluence = fluence;

	return;
}

double FluxHandler::getMaxHeFluence() const {
	return maxHeFluence;
}

bool FluxHandler::getUsingMaxHeFluence() {
	return usingMaxHeFluence;
}

void FluxHandler::setHeFlux(double flux) {
	heFlux = flux;
}

double FluxHandler::getHeFlux() const {
	return heFlux;
}

} // end namespace xolotlCore
