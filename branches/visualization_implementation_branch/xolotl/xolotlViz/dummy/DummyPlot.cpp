// Includes
#include "DummyPlot.h"

using namespace xolotlViz;

DummyPlot::DummyPlot(std::string name) : xolotlCore::Identifiable("unused") {
}

DummyPlot::~DummyPlot() {
}

void DummyPlot::render() {
	return;
}

void DummyPlot::write(std::string fileName) {
	return;
}

void DummyPlot::setPlottingStyle(PlottingStyle style) {
	return;
}

PlottingStyle DummyPlot::getPlottingStyle() {
	return PlottingStyle();
}

void DummyPlot::setDataProvider(std::shared_ptr<DataProvider> dataProvider) {
	plotDataProvider = dataProvider;
	return;
}

std::shared_ptr<DataProvider> DummyPlot::getDataProvider() const {
	return plotDataProvider;
}

void DummyPlot::setLabelProvider(std::shared_ptr<LabelProvider> labelProvider) {
	plotLabelProvider = labelProvider;
	return;
}

std::shared_ptr<LabelProvider> DummyPlot::getLabelProvider() const {
	return plotLabelProvider;
}

void DummyPlot::showLegend(bool legendShow) {
	return;
}

std::string DummyPlot::getLegend() const {
	return " ";
}

void DummyPlot::setLogScale(bool logScale) {
	return;
}

void DummyPlot::addDataProvider(std::shared_ptr<DataProvider> dataProvider){
	plotDataProvider = dataProvider;
}

std::shared_ptr<DataProvider> DummyPlot::getDataProvider(int i) const {
	return plotDataProvider;
}

int DummyPlot::getDataProviderNumber() const {
	return 0;
}
