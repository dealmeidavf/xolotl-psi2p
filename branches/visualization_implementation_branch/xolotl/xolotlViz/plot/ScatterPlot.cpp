// Includes
#include "ScatterPlot.h"
#include "eavl.h"
#include "eavlDataSet.h"
#include "eavlColor.h"
#include "eavlRenderSurfaceOSMesa.h"
#include "eavlScene.h"
#include "eavl1DWindow.h"
#include <iostream>

using namespace xolotlViz;

#define W_WIDTH 1024
#define W_HEIGHT 1024

ScatterPlot::ScatterPlot() {
}

ScatterPlot::~ScatterPlot() {
}

void ScatterPlot::render() {

	// Check if the label provider is set
	if (!plotLabelProvider){
		std::cout << "The LabelProvider is not set!!" << std::endl;
		return;
	}

	// Check if the data provider is set
	if (!plotDataProvider){
		std::cout << "The DataProvider is not set!!" << std::endl;
		return;
	}

	// Get the value that will be plotted on X and Y
	auto xVector = plotDataProvider->getAxis1Vector();
	auto yVector = plotDataProvider->getAxis2Vector();

	// Create the eavlDataSet
    eavlDataSet *data = new eavlDataSet();
    data->SetNumPoints(xVector.size());

    // Give it the xVector
	std::vector< std::vector<double> > coords;
	coords.push_back(xVector);
	std::vector<std::string> coordNames;
	coordNames.push_back("xcoord");
	AddRectilinearMesh(data, coords, coordNames, true, "RectilinearGridCells");

	// Give the yVector to the axisValues
	eavlArray *axisValues = new eavlFloatArray("coords", 1);
	axisValues->SetNumberOfTuples(data->GetNumPoints());
	for (int i = 0; i < yVector.size(); i++){
		axisValues->SetComponentFromDouble(i, 0, yVector.at(i));
	}

	// Add the axisValues to a field of the data set
	eavlField *field = new eavlField(1, axisValues, eavlField::ASSOC_POINTS);
	data->AddField(field);

    // Create an offscreen render surface
    eavlRenderSurface *surface = new eavlRenderSurfaceOSMesa;

    // Pick a background color
    eavlColor bg(1,1,1,1);

    // Create a window
    eavlScene *scene = new eavl1DGLScene();
    eavl1DWindow *window = new eavl1DWindow(bg, surface, scene);
    window->Initialize();
    window->Resize(W_WIDTH,W_HEIGHT);

    // Set up a plot for the data set
    eavlRenderer *plot;
    plot = new eavlCurveRenderer(data, NULL,
                                 eavlColor::magenta,
                                 "",
                                 "coords");
    scene->plots.push_back(plot);

    // Set the view
    scene->ResetView(window);

    // Paint
    window->Paint();

    // Save the final buffer as an image
    char fn[25];
    sprintf(fn, (plotLabelProvider->titleLabel).c_str());
    window->SaveWindowAsPNM(fn);

	return;
}
