#ifndef FLUXHANDLER_H
#define FLUXHANDLER_H

#include "IFluxHandler.h"
#include <vector>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing flux calculations.
 */
class FluxHandler: public IFluxHandler {

protected:

	/**
	 * Vector to hold the incident flux values at each grid
	 * point (x position)
	 */
	std::vector<double> incidentFluxVec;

	/**
	 * Vector to hold the position at each grid
	 * point (x position)
	 */
	std::vector<double> xGrid;

	/**
	 * Size of the surface (dy * dz) on which the flux is integrated
	 * Needed to scale the flux amplitude
	 */
	double elementarySurfaceSize;

	/**
	 * Helium fluence
	 */
	double heFluence;

	/**
	 * The amplitude of the flux
	 */
	double heFlux;

	/**
	 * Are we using a time profile for the amplitude of the helium incoming flux?
	 */
	bool useTimeProfile;

	/**
	 * Value of the fit function integrated on the grid
	 */
	double normFactor;

	/**
	 * Vector to hold the time read from the input
	 * time profile file
	 */
	std::vector<double> time;

	/**
	 * Vector to hold the amplitude read from the input
	 * time profile file
	 */
	std::vector<double> amplitude;

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * It needs to be implemented by the daughter classes.
	 *
	 * @param x The position where to evaluate he fit
	 * @return the evaluated value
	 */
	virtual double FitFunction(double x) {return 0.0;}

	/**
	 * This method returns the value of the helium incident flux amplitude at the
	 * given time when a time profile is used.
	 * @param currentTime The time
	 * @return The value of the helium flux at this time
	 */
	double getAmplitude(double currentTime) const;

	/**
	 * This method recomputes the values of the incident flux vector when
	 * a time profile is given.
	 */
	void recomputeFluxHandler();

public:

	FluxHandler();

	~FluxHandler() {
	}

	/**
	 * Function to calculate and store the incident flux values at each grid point
	 * @param grid The grid on the x axis
	 * @param hy The step size between grid points on the y axis
	 * @param hz The step size between grid points on the z axis
	 */
	virtual void initializeFluxHandler(std::vector<double> grid, double hy = 1.0,
			double hz = 1.0);

	/**
	 * This method reads the values on the time profile file and store them in the
	 * time and amplitude vectors.
	 * @param fileName The name of the file where the values are stored
	 */
	void initializeTimeProfile(std::string fileName);

	/**
	 * This operation returns the incident flux vector
	 * @param currentTime     	 The time
	 * @return incidentFluxVec   The incident flux vector
	 */
	virtual std::vector<double> getIncidentFluxVec(double currentTime);

	/**
	 * Given a specific concentration, position, and time, this operation sets the outgoing
	 * flux to the specified amount.
	 * @param composition  The composition of the cluster
	 * @param position     The position of the cluster
	 * @param time         The time
	 * @return outgoingFlux  The outgoing flux at the given position and time of the cluster with
	 * the specified composition
	 */
	virtual void setOutgoingFlux(std::vector<int> compositionVec,
			std::vector<int> position, double time, double outgoingFlux);

	/**
	 * This operation increments the Helium fluence at the current time step.
	 * @param dt			The length of the time step
	 */
	virtual void incrementHeFluence(double dt);

	/**
	 * This operation returns the Helium fluence
	 * @return	The Helium fluence at current time step
	 */
	virtual double getHeFluence() const;

	/**
	 * This operation sets the factor to change the Helium flux.
	 * @param flux	Helium flux value
	 */
	virtual void setHeFlux(double flux);

	/**
	 * This operation gets the factor that changes the Helium flux.
	 * @return	Helium flux value
	 */
	virtual double getHeFlux() const;

};
//end class FluxHandler

}

#endif
