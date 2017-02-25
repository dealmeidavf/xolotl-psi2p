#ifndef IIDENTIFIABLE_H
#define IIDENTIFIABLE_H

#include <string>

namespace xolotlCore 
{

/**
 * An Identifiable is an object that has some identification (e.g., a name
 * and/or an ID).
 */
 class IIdentifiable
 {
  public:

/**
* Obtain the object's given name.
*/
  virtual std::string getName( void ) const = 0;
 };

} // end namespace xolotlCore

#endif // IIDENTIFIABLE_H
