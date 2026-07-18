#ifndef _META_LEARNER_
#define _META_LEARNER_

#include <map>
#include <unordered_map>
#include <unordered_set>
#include "CommonTypes.H"
#include "HierarchicalCluster.H"

using namespace std;

class VariableManager;
class EvidenceSource;
class EvidenceSet;
class PotentialManager;
class MetaMove;
class HierarchicalCluster;
class Matrix;
class Checkpoint;

class MetaLearner
{
	public:
		MetaLearner();
		void setShouldLoadCheckpoint(bool shouldLoadCheckpointVal);
		void setMaxFactorSize(int);
		void setBeta1(double aval);
		void setBeta_Motif(double bval);
		void setConvergenceThreshold(double);
		void setRestrictedList(const char*);
		int initEdgePriorMeta(const string& priorName, map<string, map<string, double>*>& graph, unordered_map<int, unordered_map<int, double>*>& edgePrior);
		void doCrossValidation(int);
		void setEvidenceSource(EvidenceSource*);
		void setVariableManager(VariableManager*);
		void setOutputDirName(const char*);
		void setClusteringThreshold(double);
		void setSpecificFold(int);
		void readModuleMembership(const char*);
		void setDefaultModuleMembership();
		int setPriorGraph_All(const char*);

	private:

		void start(int fold, EvidenceSet& trainSet);
		void setupFoldState(Checkpoint& checkpoint, int& iter, bool& notConverged, double& globalScore);
		void setupFoldState(double& globalScore, EvidenceSet& trainSet);
		void restoreCheckpointModules(const unordered_map<string, int>& checkpointGeneModuleIDs);
		void restoreCheckpointGraph(const vector<pair<string, string>>& checkpointEdges);
		void writeFoldProgress(int currFold, int iter, bool notConverged, Checkpoint& checkpoint);
		void initEdgePriorMeta_All();
		double getInitPLLScore(int vId, EvidenceSet& trainSet);
		double getPLLScore();
		double getEdgePrior(int tfID, int targetID);
		void clearFoldSpecData();
		int initEdgeSet();
		bool getNextMove(int vID, int sampleCount, MetaMove& outMove);
		void makeMove(MetaMove& move, int currIteration);
		void dumpAllGraphs(int);
		int setPriorGraph(const char* aFName, map<string, map<string, double>*>& priorGraph);
		void redefineModules(int currFold, EvidenceSet& trainSet);
		double getModuleContribLogistic(string& tgtName, string& tfName);
		void initPhysicalDegree();
		void getEnrichedTFs(map<string, int>& tfSet, map<string, int>* genes, map<string, map<string, double>*>& edgeSet);
		void initDistances(EvidenceSet& trainSet);
		void updateSharedParentDistances();

		FactorGraph* factorGraph;
		PotentialManager* potManager;
		EvidenceSource* evidenceSource;
		VariableManager* varManager;
		
		bool shouldLoadCheckpoint;
		map<string, double> betaMap;
		map<string, unordered_map<int, unordered_map<int, double>*>*> edgePriorMap;
		map<string, map<string, map<string, double>*>*> priorGraphMap;
		int maxFactorSize;
		double convThreshold;
		char restrictedFName[1024];
		unordered_map<int, unordered_map<int, int>> edgeMap;
		map<string,int> restrictedVarList;
		unordered_map<int, double>* currPLL;
		map<int, double> finalScores;
		char outputDirName[1024];
		map<string, double> edgePresenceProb;
		map<int, double> varNeighborhoodPrior;
		double beta1;
		double beta_motif;
		map<int, map<string, int>*> moduleGeneSet;
		unordered_map<string, int> geneModuleID;
		unordered_map<int, unordered_map<string, int>*> moduleIndegree;
		unordered_map<string, int> regulatorModuleOutdegree;
		double clusterThreshold;
		unordered_map<string, int> variableStatus;
		char foldoutDirName[1024];
		int specificFold;
		HierarchicalCluster hc;

		// Used for defining modules, the distance metric based on expression.
		Matrix* correlationDistances;

		// Used for defining modules, the distance metric based on shared parents.
		Matrix* sharedParentDistances;

		// Tracks whether two nodes share at least one parent.
		unordered_map<int, unordered_map<int, int>> sharedParents;

		// Var ids of any nodes with a parent added this iteration.
		unordered_set<int> updatedThisIteration;
};
#endif
