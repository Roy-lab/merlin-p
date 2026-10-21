#ifndef _MODULE_MANAGER
#define _MODULE_MANAGER

#include <unordered_map>
#include <unordered_set>
#include <map>

using namespace std;

class VariableSet;

class ModuleManager
{
    public:

        void readModuleMembership(const char* aFName, VariableSet* varSet);

        void setDefaultModuleMembership(VariableSet* varSet);

        void restoreCheckpointModules(const unordered_map<string, int>& checkpointGeneModuleIDs, VariableSet* variableSet);

        void resetModuleAssignments();

        int getModuleID(const string& geneName);

        const unordered_set<string>& getModuleMembers(int moduleID);

        void setModuleMembers(int moduleID, unordered_set<string>& geneSet);

        int moduleCount();

        auto begin() {
            return moduleGeneSet.cbegin();
        }

        auto end() {
            return moduleGeneSet.cend();
        }

    private:
		unordered_map<int, unordered_set<string>> originalModuleGeneSet;
		unordered_map<string, int> originalGeneModuleID;

		unordered_map<int, unordered_set<string>> moduleGeneSet;
		unordered_map<string, int> geneModuleID;

};

#endif