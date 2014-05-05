// Includes
#include "Plot.h"

using namespace xolotlViz;

Plot::Plot(std::string name) : xolotlCore::Identifiable(name) {
}

Plot::~Plot() {
}

void Plot::render(std::string fileName) {
	return;
}

void Plot::write(std::string fileName) {
	render(fileName);
	return;
}

void Plot::setPlottingStyle(PlottingStyle style) {
	plotStyle = style;
	return;
}

PlottingStyle Plot::getPlottingStyle() {
	return plotStyle;
}

void Plot::setDataProvider(std::shared_ptr<DataProvider> dataProvider) {
	plotDataProvider = dataProvider;
	return;
}

std::shared_ptr<DataProvider> Plot::getDataProvider() const {
	return plotDataProvider;
}

void Plot::setLabelProvider(std::shared_ptr<LabelProvider> labelProvider) {
	plotLabelProvider = labelProvider;
	return;
}

std::shared_ptr<LabelProvider> Plot::getLabelProvider() const {
	return plotLabelProvider;
}

void Plot::showLegend(bool legendShow) {
	enableLegend = legendShow;
	return;
}

std::string Plot::getLegend() const {
	return " ";
}

void Plot::setLogScale(bool logScale) {
	enableLogScale = logScale;
	return;
}

void Plot::addDataProvider(std::shared_ptr<DataProvider> dataProvider) {
	return;
}

std::shared_ptr<DataProvider> Plot::getDataProvider(int i) const {
	return plotDataProvider;
}

int Plot::getDataProviderNumber() const {
	return 0;
}
