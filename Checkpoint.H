#ifndef _CHECKPOINT_
#define _CHECKPOINT_

#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

class Variable;

class Checkpoint
{
    public:
        Checkpoint(const char* dirName, int maxFactorSize);

        bool loadCheckpointData();

        int getIteration();

        bool getHasNotConverged();

        double getInitialScore();

        unordered_map<string, double> getInitialPLLs();

        unordered_map<string, int> getVariableStatus();

        unordered_map<string, int> getGeneModuleIDs();

        vector<pair<string, string>> getEdges();

        void writeCheckpointMetadata(int iter, bool notConvergedVal);

        void writePLLScore(unordered_map<int, double>* currPLL, const vector<Variable*>& varSet);

        void writeLastUpdate(unordered_map<string, int>& variableStatus);

        void doTar();

    private:
        bool loadPLLScore();
        bool loadLastUpdate();
        bool loadModules();
        bool loadEdges();

        int iter;
        bool hasNotConverged;
        double initialScore;
        unordered_map<string, double> initialPLLs;
        unordered_map<string, int> geneModuleIDs;
        vector<pair<string, string>> edges;
        unordered_map<string, int> variableStatus;
        string dirName;
        int maxFactorSize;
};

#endif