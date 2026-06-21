#include "Checkpoint.H"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include "Variable.H"

Checkpoint::Checkpoint(const char* inDirName, int inMaxFactorSize)
{
    dirName = string(inDirName);
    maxFactorSize = inMaxFactorSize;
}

int
Checkpoint::getIteration()
{
    return iter;
}

bool
Checkpoint::getHasNotConverged()
{
    return hasNotConverged;
}

double
Checkpoint::getInitialScore()
{
    return initialScore;
}

unordered_map<string, double>
Checkpoint::getInitialPLLs()
{
    return initialPLLs;
}

unordered_map<string, int>
Checkpoint::getVariableStatus()
{
    return variableStatus;
}

unordered_map<string, int>
Checkpoint::getGeneModuleIDs()
{
    return geneModuleIDs;
}

vector<pair<string, string>>
Checkpoint::getEdges()
{
    return edges;
}

bool
Checkpoint::loadCheckpointData()
{
    string fileName = dirName + "/checkpoint.txt";
    ifstream inFile(fileName);

    if (!inFile.is_open()) {
        cerr << "Error: could not open checkpoint file: " << fileName << endl;
        return false;
    }

    string label;
    inFile >> label >> iter;
    inFile >> label >> hasNotConverged;
    inFile.close();

    // Increment by one, because we want to start one iteration past the last completed iteration.
	iter++;

    bool success = true;
    success &= loadPLLScore();
    success &= loadModules();
    success &= loadEdges();
    success &= loadLastUpdate();
    return success;
}

bool
Checkpoint::loadPLLScore()
{
    string fileName = dirName + "/pll.txt";
    ifstream inFile(fileName);

	if (!inFile.is_open()) {
		std::cerr << "Error: could not open checkpoint PLL file: " << fileName << std::endl;
		return false;
	}

    double initScore = 0;
    string varName;
    double val;

    while (inFile >> varName >> val) {
        initialPLLs[varName] = val;
        initScore += val;
    }

    initialScore = initScore;

    return true;
}

bool
Checkpoint::loadModules()
{
    string fileName = dirName + "/modules.txt";
	ifstream inFile(fileName);

	if (!inFile.is_open()) {
		std::cerr << "Error: could not open checkpoint module file: " << fileName << std::endl;
		return false;
	}

    string geneName;
    int moduleID;

    // Read clustered modules from modules.txt
    while (inFile >> geneName >> moduleID) {
        geneModuleIDs[geneName] = moduleID;
    }

    return true;
}

bool
Checkpoint::loadLastUpdate()
{
    string fileName = dirName + "/lastUpdate.txt";
    ifstream inFile(fileName);

	if (!inFile.is_open()) {
		std::cerr << "Error: could not open checkpoint last update file: " << fileName << std::endl;
		return false;
	}

    string varName;
    int lastUpdateIter;

    while (inFile >> varName >> lastUpdateIter) {
        variableStatus[varName] = lastUpdateIter;
    }

    return true;
}

bool
Checkpoint::loadEdges()
{
    string fileName = dirName + "/prediction_k" + to_string(maxFactorSize) + ".txt";
    ifstream inFile(fileName);

    if (!inFile.is_open()) {
		std::cerr << "Error: could not open checkpoint graph file: " << fileName << std::endl;
		return false;
	}

    string regName;
    string targetName;
    double weight;

    while (inFile >> regName >> targetName >> weight) {
        edges.push_back(pair<string, string>(regName, targetName));
    }

    return true;
}

void
Checkpoint::writeCheckpointMetadata(int iter, bool notConvergedVal)
{
    string fileName = dirName + "/checkpoint.txt";
    ofstream outFile(fileName);
    outFile << "iter " << iter << endl;
    outFile << "notConverged " << notConvergedVal << endl;
}

void
Checkpoint::writePLLScore(unordered_map<int, double>* currPLL, const unordered_map<int, Variable*>& varSet)
{
	if (currPLL == NULL) {
		return;
	}

    string fileName = dirName + "/pll.txt";
    ofstream outFile(fileName);

	// Force maximum double precision to prevent truncation
    outFile << std::setprecision(std::numeric_limits<double>::max_digits10);

    for(auto iter = currPLL->begin(); iter != currPLL->end(); iter++) {
        auto varIter = varSet.find(iter->first);
        if (varIter == varSet.end()) {
            continue;
        }
        Variable* var = varIter->second;
        outFile << var->getName() << "\t" << iter->second << endl;
	}
}

void
Checkpoint::writeLastUpdate(unordered_map<string, int>& variableStatus)
{
    string fileName = dirName + "/lastUpdate.txt";
    ofstream outFile(fileName);

    for (auto iter = variableStatus.begin(); iter != variableStatus.end(); iter++) {
        outFile << iter->first << "\t" << iter->second << endl;
    }
}

void
Checkpoint::doTar()
{
    char tarCmd[1024];
    sprintf(tarCmd, "tar czf %s.tar.gz %s", dirName.c_str(), dirName.c_str());
    system(tarCmd);
}
