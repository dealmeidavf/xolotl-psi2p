/*
 * Constants.h
 *
 *  Created on: May 6, 2013
 *      Author: bkj
 */

#include <memory>

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

namespace xolotlCore {

//! Definitions of fundamental constants used in Xolotl

//! The Boltzmann constant in units of eV K^-1.
static const double kBoltzmann = 8.61733240000000000E-5;

//! Pi, taken from "100000 digits of Pi,"
//! at http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
static const double pi = 3.1415926535897932;

//! Lattice Parameter. Equal to 3.17 Angstroms, taken from Becquart et. al.
//! Journal of Nuclear Materials 403 (2010) 75–88. Given in units here of nm.
static const double latticeConstant = 0.31700000000000000;

/**
 * Statically defined type names. This are used to optimize the lookup
 * calls in the maps so that new strings are not constantly created
 * dynamically. It accounts for a small but significant amount of
 * performance.
 */
static const std::string heType = "He";
static const std::string vType = "V";
static const std::string iType = "I";

} /* end namespace xolotlCore */
#endif /* CONSTANTS_H_ */
