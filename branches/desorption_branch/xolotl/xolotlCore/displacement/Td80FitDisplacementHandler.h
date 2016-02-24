#ifndef TD80FITDISPLACEMENTHANDLER_H
#define TD80FITDISPLACEMENTHANDLER_H

#include "DisplacementHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IDisplacementHandler interface to calculate the initial vacancy distribution
 * in tungsten material with a threshold energy of 80 eV.
 */
class Td80FitDisplacementHandler: public DisplacementHandler {
private:

	/**
	 * Function that calculates the vacancy at a given position x (in nm).
	 * This function is not normalized.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double VacancyFitFunction(double x) {
		// Value at which the vacancy goes to 0
		double x1 = 7.9;

		if (x > x1) return 0.0;

		double x2 = 2.0 * x/1.000100000000000122e+01-1.0;

		// Compute the fit
		double value = 0.026723
				- 0.047913 * x2
				+ 0.018489 * 1./2. * (3.0 * pow(x2, 2)-1.0)
				+ 0.019202 * 1./2. * (5.0 * pow(x2, 3)-3.0 * x2)
				- 0.030103 * 1./8. * (35.0 * pow(x2, 4)-30.0 * pow(x2, 2) + 3.0)
				+ 0.019085 * 1./8. * (63.0 * pow(x2, 5)-70.0 * pow(x2, 3) + 15.0 * x2)
				- 0.006230 * 1./16. * (231.0 * pow(x2, 6)-315.0 * pow(x2, 4) + 105.0 * (x2, 2) - 5.0);


		return std::max(value,0.0);
	}

public:

	/**
	 * The constructor
	 */
	Td80FitDisplacementHandler() {}

	/**
	 * The Destructor
	 */
	~Td80FitDisplacementHandler() {}

};
//end class Td80FitDisplacementHandler

}

#endif
