#ifndef _DISTANCE_MANAGER
#define _DISTANCE_MANAGER

#include <unordered_map>
#include <unordered_set>

using namespace std;

class Matrix;
class EvidenceSet;
class FactorGraph;

class DistanceManager
{
    public:
        DistanceManager();

        Matrix* getCorrelationDistances();

        Matrix* getSharedParentDistances();

        void clearFoldData();

        void clearIterationData();

        void initDistances(EvidenceSet* trainSet, int varCount);

        void updateSharedParentDistances(FactorGraph* graph, int varCount);

        virtual void restoreCheckpointSharedParents(unordered_map<int, unordered_map<int, int>> edgeMap);

        virtual void addSharedParents(unordered_map<int, unordered_map<int, int>> edgeMap, int parentID, int childID);

        void addUpdatedThisIteration(int targetID);

    protected:
        // Tracks whether two nodes share at least one parent.
		unordered_map<int, unordered_map<int, int>> sharedParents;

    private:
        // Used for defining modules, the distance metric based on expression.
		Matrix* correlationDistances;

		// Used for defining modules, the distance metric based on shared parents.
		Matrix* sharedParentDistances;

		// Var ids of any nodes with a parent added this iteration.
		unordered_set<int> updatedThisIteration;
};

#endif