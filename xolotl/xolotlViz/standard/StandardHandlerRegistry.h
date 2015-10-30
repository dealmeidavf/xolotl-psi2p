#ifndef STANDARDHANDLERREGISTRY_H
#define STANDARDHANDLERREGISTRY_H

#include <iostream>
#include <string>
#include <map>
#include "IVizHandlerRegistry.h"


namespace xolotlViz {

/**
 * Factory for creating standard plots using EAVL and MESA libraries.
 * This is used only if the libraries are present and if the user uses
 * the standard registry.
 */
class StandardHandlerRegistry : public IVizHandlerRegistry
{
public:
    /// possible rendering/output types
    enum OutputType
    {
        png,  ///< Raster, with PNG output
        eps   ///< Vector, with EPS output
    };

    /**
     * Construct a StandardHandlerRegistry.
     */
    StandardHandlerRegistry(OutputType otype);

    /**
     * Clean up a StandardHandlerRegistry.
     */
    virtual ~StandardHandlerRegistry();

    /**
     * Obtain a Plot by name.
     *
     * @param name The name of the Plot.
     * @param type The type of plot to return.
     * @return A shared pointer to the newly-created Plot.
     */
    virtual std::shared_ptr<IPlot> getPlot(const std::string& name, PlotType type);

protected:
    OutputType outputType;

};  //end class StandardHandlerRegistry

}//end namespace xolotlViz

#endif
