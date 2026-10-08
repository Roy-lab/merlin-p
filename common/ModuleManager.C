#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include "ModuleManager.H"
#include "VariableSet.H"
#include "Variable.H"
#include "gsl/gsl_randist.h"

void
ModuleManager::readModuleMembership(const char* aFName, VariableSet* variableSet)
{
	int largestModuleID = 0;

	ifstream inFile(aFName);
	char buffer[1024];
	while (inFile.good()) {

		inFile.getline(buffer, 1023);

        if (strlen(buffer) <= 0) {
			continue;
		}

        string geneName;
		int moduleID;
		int tokCnt = 0;
        char* tok = strtok(buffer, "\t");

		while (tok != NULL) {
			if (tokCnt == 0) {
				geneName.append(tok);
			} else if (tokCnt == 1) {
				moduleID = atoi(tok);
			}
			tok = strtok(NULL, "\t");
			tokCnt++;
		}

        moduleGeneSet[moduleID].insert(geneName);
		geneModuleID[geneName] = moduleID;

		if (moduleID > largestModuleID) {
			largestModuleID = moduleID;
		}
	}
	inFile.close();

	vector<Variable*>& varSet = variableSet->getVariables();

	for (Variable *variable : varSet) {
		string geneName = variable->getName();
		if(geneModuleID.find(geneName) != geneModuleID.end()) {
			continue;
		}
		largestModuleID++;
		moduleGeneSet[largestModuleID].insert(geneName);
		geneModuleID[geneName] = largestModuleID;
	}

	originalModuleGeneSet = moduleGeneSet;
	originalGeneModuleID = geneModuleID;
}

void
ModuleManager::setDefaultModuleMembership(VariableSet* variableSet)
{
	vector<Variable*>& varSet = variableSet->getVariables();
	int vCnt = varSet.size();

    int moduleCnt = (int) sqrt(vCnt / 2);
	if (moduleCnt > 30)
	{
		moduleCnt = 30;
	}

    // Randomly partition the variables into cluster assignments

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    double step = 1.0 / (double)vCnt;

	vector<int> randIndex;
	map<int, int> usedInit;
	for (int i = 0; i < vCnt; i++)
	{
		double rVal = gsl_ran_flat(rng, 0, 1);
		int rind = (int)(rVal / step);
		while (usedInit.find(rind) != usedInit.end())
		{
			rVal = gsl_ran_flat(rng, 0, 1);
			rind = (int)(rVal / step);
		}
		usedInit[rind] = 0;
		randIndex.push_back(rind);
	}

	int clusterSize = vCnt / moduleCnt;
	for (int e = 0; e < moduleCnt; e++)
	{
		int startInd = e * clusterSize;
		int endInd = (e + 1) * clusterSize;
		if (e == moduleCnt - 1)
		{
			endInd = vCnt;
		}

		for (int i = startInd; i < endInd; i++)
		{
			int dataId = randIndex[i];
			Variable* v = varSet[dataId];
            moduleGeneSet[e].insert(v->getName());
			geneModuleID[v->getName()] = e;
		}
	}

	gsl_rng_free(rng);

	originalModuleGeneSet = moduleGeneSet;
	originalGeneModuleID = geneModuleID;
}

void
ModuleManager::restoreCheckpointModules(const unordered_map<string, int>& checkpointGeneModuleIDs, VariableSet* variableSet)
{
	cout << "Load checkpoint modules..." << endl;

	geneModuleID = checkpointGeneModuleIDs;

	// Clear out the initial module gene sets
	moduleGeneSet.clear();

    // Build the module gene sets
	int largestModuleID = 0;
	for (auto iter = geneModuleID.begin(); iter != geneModuleID.end(); iter++) {
        int moduleID = iter->second;
		moduleGeneSet[moduleID].insert(iter->first);
        if (moduleID > largestModuleID) {
            largestModuleID = moduleID;
        }
	}

	vector<Variable*>& varSet = variableSet->getVariables();

    // Recreate singleton modules (for parentless genes) that were not written to modules.txt
    int genesWithNoNeighborsCount = 0;
    for(auto vIter = varSet.begin(); vIter != varSet.end(); vIter++) {
        string geneName = (*vIter)->getName();
        if(geneModuleID.find(geneName) != geneModuleID.end()) {
            continue;
        }
        largestModuleID++;
        moduleGeneSet[largestModuleID].insert(geneName);
        geneModuleID[geneName] = largestModuleID;
        genesWithNoNeighborsCount++;
    }

	originalModuleGeneSet = moduleGeneSet;
	originalGeneModuleID = geneModuleID;

    cout << "Recovered " << genesWithNoNeighborsCount << " parentless genes not present in modules.txt" << endl;
    cout << "Total modules after recovery: " << moduleGeneSet.size() << endl;
}

void
ModuleManager::resetModuleAssignments()
{
    moduleGeneSet = originalModuleGeneSet;
    geneModuleID = originalGeneModuleID;
}

int
ModuleManager::getModuleID(const string& geneName)
{
    return geneModuleID.at(geneName);
}

const unordered_set<string>&
ModuleManager::getModuleMembers(int moduleID)
{
    return moduleGeneSet.at(moduleID);
}

void
ModuleManager::setModuleMembers(int moduleID, unordered_set<string>& geneSet)
{
    moduleGeneSet[moduleID] = geneSet;

    for (const string& geneName : geneSet)
    {
        geneModuleID[geneName] = moduleID;
    }
}

int
ModuleManager::moduleCount()
{
    return moduleGeneSet.size();
}