#include "PSIClusterReactionNetwork.h"
#include "PSICluster.h"
#include <HandlerRegistryFactory.h>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace xolotlCore;
using std::shared_ptr;

void PSIClusterReactionNetwork::setDefaultPropsAndNames() {

	// Shared pointers for the cluster type map
	std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>heVector
	= std::make_shared<std::vector<std::shared_ptr<Reactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>> vVector
	= std::make_shared<std::vector<std::shared_ptr<Reactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>> iVector
	= std::make_shared<std::vector<std::shared_ptr<Reactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>> heVVector
	= std::make_shared<std::vector<std::shared_ptr<Reactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>> heIVector
	= std::make_shared<std::vector<std::shared_ptr<Reactant>>>();

	// Initialize the list of all reactants used by getAll.
	allReactants = std::make_shared<std::vector<std::shared_ptr<Reactant>>>();

	// Initialize default properties
	(*properties)["reactionsEnabled"] = "true";
	(*properties)["dissociationsEnabled"] = "true";
	(*properties)["numHeClusters"] = "0";
	(*properties)["numVClusters"] = "0";
	(*properties)["numIClusters"] = "0";
	(*properties)["numHeVClusters"] = "0";
	(*properties)["numHeIClusters"] = "0";
	(*properties)["maxHeClusterSize"] = "0";
	(*properties)["maxVClusterSize"] = "0";
	(*properties)["maxIClusterSize"] = "0";
	(*properties)["maxHeVClusterSize"] = "0";
	(*properties)["maxHeIClusterSize"] = "0";

	// Initialize the current and last size to 0
	networkSize = 0;
	// Set the reactant names
	names.push_back("He");
	names.push_back("V");
	names.push_back("I");
	// Set the compound reactant names
	compoundNames.push_back("HeV");
	compoundNames.push_back("HeI");

	// Setup the cluster type map
	clusterTypeMap["He"] = heVector;
	clusterTypeMap["V"] = vVector;
	clusterTypeMap["I"] = iVector;
	clusterTypeMap["HeV"] = heVVector;
	clusterTypeMap["HeI"] = heIVector;

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork() :
		ReactionNetwork( )
{

	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork(registry){

	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}


PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		const PSIClusterReactionNetwork &other) :
		ReactionNetwork(other) {

	// The size and ids do not need to be copied. They will be fixed when the
	// reactants are added.

	// Copy the names
	names = other.getNames();
	// Copy the compound names
	compoundNames = other.getCompoundNames();
	// Reset the properties table so that it can be properly updated when the
	// network is filled.
	setDefaultPropsAndNames();
	// Get all of the reactants from the other network and add them to this one
	auto reactants = other.getAll();
	for (int i = 0; i < reactants->size(); i++) {
		add(reactants->at(i)->clone());
	}

}


/**
 * This operation returns a reactant with the given name and size if it
 * exists in the network or null if not.
 * @param name the name of the reactant
 * @param size the size of the reactant
 * @return A shared pointer to the reactant
 */
std::shared_ptr<Reactant> PSIClusterReactionNetwork::get(
		const std::string rName, const int size) const {

	// Local Declarations
	std::map<std::string, int> composition;
	std::shared_ptr<PSICluster> retReactant;

	// Setup the composition map to default values
	composition["He"] = 0;
	composition["V"] = 0;
	composition["I"] = 0;

	// Only pull the reactant if the name and size are valid
	if ((rName == "He" || rName == "V" || rName == "I") && size >= 1) {
		composition[rName] = size;
		// Make sure the reactant is in the map
		if (singleSpeciesMap.count(composition)) {
			retReactant = singleSpeciesMap.at(composition);
		}
	}

	return retReactant;
}

/**
 * This operation returns a compound reactant with the given name and size if it
 * exists in the network or null if not.
 * @param rName the name of the compound reactant
 * @param sizes an array containing the sizes of each piece of the reactant
 * @return A shared pointer to the compound reactant
 */
std::shared_ptr<Reactant> PSIClusterReactionNetwork::getCompound(
		const std::string rName, const std::vector<int> sizes) const {

	// Local Declarations
	std::map<std::string, int> composition;
	std::shared_ptr<PSICluster> retReactant;

	// Setup the composition map to default values
	composition["He"] = 0;
	composition["V"] = 0;
	composition["I"] = 0;

	// Only pull the reactant if the name is valid and there are enough sizes
	// to fill the composition.
	if ((rName == "HeV" || rName == "HeI") && sizes.size() == 3) {
		composition["He"] = sizes[0];
		composition["V"] = sizes[1];
		composition["I"] = sizes[2];
		// Make sure the reactant is in the map
		if (mixedSpeciesMap.count(composition)) {
			retReactant = mixedSpeciesMap.at(composition);
		}
	}

	return retReactant;
}

/**
 * This operation returns all reactants in the network without regard for
 * their composition or whether they are compound reactants. The list may
 * or may not be ordered and the decision is left to implementers.
 * @return The list of all of the reactants in the network
 */
std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>PSIClusterReactionNetwork::getAll() const {

	// Reload the array of all reactants if the size has changed
	if (allReactants->size() != networkSize) {
		// Clear the list
		allReactants->clear();
		// Load the single-species clusters
		for (auto it = singleSpeciesMap.begin(); it != singleSpeciesMap.end(); ++it) {
			allReactants->push_back(it->second);
		}
		// Load the mixed-species clusters
		for (auto it = mixedSpeciesMap.begin(); it != mixedSpeciesMap.end(); ++it) {
			allReactants->push_back(it->second);
		}
	}

//	std::cout << "Num single species = " << singleSpeciesMap.size() << std::endl;
//	std::cout << "Num mixed species = " << mixedSpeciesMap.size() << std::endl;
//	std::cout << "Num species = " << allReactants->size() << std::endl;

	return allReactants;
}

/**
 * This operation returns all reactants in the network with the given name.
 * The list may or may not be ordered and the decision is left to
 * implementers.
 * @param name The reactant or compound reactant name
 * @return The list of all of the reactants in the network or null if the
 * name is invalid.
 */
std::shared_ptr<std::vector<std::shared_ptr<Reactant> > > PSIClusterReactionNetwork::getAll(
		std::string name) const {

	// Local Declarations
	std::shared_ptr<std::vector<std::shared_ptr<Reactant>> > reactants(
			new std::vector<std::shared_ptr<Reactant>>());

	// Only pull the reactants if the name is valid
	if (name == "He" || name == "V" || name == "I" || name == "HeV"
			|| name == "HeI") {
		std::shared_ptr<std::vector<std::shared_ptr<Reactant>> > storedReactants =
				clusterTypeMap.at(name);
		int vecSize = storedReactants->size();
		for (int i = 0; i < vecSize; i++) {
			reactants->push_back(storedReactants->at(i));
		}
	}

	return reactants;
}

/**
 * This operation adds a reactant or a compound reactant to the network.
 * @param reactant The reactant that should be added to the network.
 */
void PSIClusterReactionNetwork::add(std::shared_ptr<Reactant> reactant) {

	// Local Declarations
	int numHe = 0, numV = 0, numI = 0;
	bool isMixed = false;
	std::string numClusterKey, clusterSizeKey, name;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Check the name
		name = reactant->getName();
		if (name == "He" || name == "V" || name == "I" || name == "HeV"
				|| name == "HeI") {
			// Get the composition
			auto composition = reactant->getComposition();
			// Get the species sizes
			numHe = composition.at("He");
			numV = composition.at("V");
			numI = composition.at("I");
			// Determine if the cluster is a compound. If there is more than one
			// type, then the check below will sum to greater than one and we know
			// that we have a mixed cluster.
			isMixed = ((numHe > 0) + (numV > 0) + (numI > 0)) > 1;
			// Only add the element if we don't already have it
			// Add the compound or regular reactant.
			if (isMixed && mixedSpeciesMap.count(composition) == 0) {
				// Put the compound in its map
				mixedSpeciesMap[composition] = std::dynamic_pointer_cast<
						PSICluster>(reactant);
				// Figure out whether we have HeV or HeI and set the keys
				if (numV > 0) {
					numClusterKey = "numHeVClusters";
					clusterSizeKey = "maxHeVClusterSize";
				} else {
					numClusterKey = "numHeIClusters";
					clusterSizeKey = "maxHeIClusterSize";
				}
			} else if (!isMixed && singleSpeciesMap.count(composition) == 0) {
				/// Put the reactant in its map
				singleSpeciesMap[composition] = std::dynamic_pointer_cast<
						PSICluster>(reactant);
				// Figure out whether we have He, V or I and set the keys
				if (numHe > 0) {
					numClusterKey = "numHeClusters";
					clusterSizeKey = "maxHeClusterSize";
				} else if (numV > 0) {
					numClusterKey = "numVClusters";
					clusterSizeKey = "maxVClusterSize";
				} else {
					numClusterKey = "numIClusters";
					clusterSizeKey = "maxIClusterSize";
				}
			} else {
				std::stringstream errStream;
				errStream << "PSIClusterReactionNetwork Message: "
						<< "Duplicate Reactant (He=" << numHe << ",V=" << numV
						<< ",I=" << numI << ") not added!" << std::endl;
				throw errStream.str();
			}
			// Increment the number of total clusters of this type
			int numClusters = std::stoi(properties->at(numClusterKey));
			numClusters++;
			(*properties)[numClusterKey] = std::to_string(
					(long long) numClusters);
			// Increment the max cluster size key
			int maxSize = std::stoi(properties->at(clusterSizeKey));
			int clusterSize = numHe + numV + numI;
			maxSize = std::max(clusterSize, maxSize);
			(*properties)[clusterSizeKey] = std::to_string((long long) maxSize);
			// Update the size
			++networkSize;
			// Set the id for this cluster
			reactant->setId(networkSize);
			// Get the vector for this reactant from the type map
			auto clusters = clusterTypeMap[name];
			clusters->push_back(reactant);
//				std::cout << "Num single species = " << singleSpeciesMap.size()
//						<< std::endl;
//				std::cout << "Num mixed species = " << mixedSpeciesMap.size()
//						<< std::endl;
//				std::cout << "Num species = " << networkSize << std::endl;
//				std::cout << "Real num species = " << (singleSpeciesMap.size() + mixedSpeciesMap.size()) << std::endl;
		}
	}

	return;
}

/**
 * This operation returns the names of the reactants in the network.
 * @return A vector with one each for each of the distinct reactant types
 * in the network.
 */
const std::vector<std::string> & PSIClusterReactionNetwork::getNames() const {

	return names;
}

/**
 * This operation returns the names of the compound reactants in the
 * network.
 * @return A vector with one each for each of the distinct compound
 * reactant types in the network.
 */
const std::vector<std::string> & PSIClusterReactionNetwork::getCompoundNames() const {

	return compoundNames;
}

/**
 * This operation returns a map of the properties of this reaction network.
 * @return The map of properties that has been configured for this
 * ReactionNetwork.
 */
const std::map<std::string, std::string> & PSIClusterReactionNetwork::getProperties() {

	return *properties;
}

/**
 * This operation sets a property with the given key to the specified value
 * for the network. ReactionNetworks may reserve the right to ignore this
 * operation for special key types, most especially those that they manage on
 * their own.
 *
 * Calling this operation will ignore all of the published properties that are
 * configured by default.
 * @param key The key for the property
 * @param value The value to which the key should be set.
 */
void PSIClusterReactionNetwork::setProperty(std::string key,
		std::string value) {

// Check the keys and value before trying to set the property
	if (!key.empty() && !value.empty() && key != "numHeClusters"
			&& key != "numVClusters" && key != "numIClusters"
			&& key != "maxHeClusterSize" && key != "maxVClusterSize"
			&& key != "maxIClusterSize" && key != "maxHeVClusterSize"
			&& key != "maxHeIClusterSize") {
		// Add the property if it made it through that!
		(*properties)[key] = value;
	}

	return;
}

/**
 * This operation returns the size or number of reactants in the network.
 * @return The number of reactants in the network
 */int PSIClusterReactionNetwork::size() {
	return networkSize;
}
