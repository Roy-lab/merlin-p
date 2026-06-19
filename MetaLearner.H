#ifndef _META_LEARNER_
#define _META_LEARNER_

#include <map>
#include "CommonTypes.H"
#include "HierarchicalCluster.H"

using namespace std;

class VariableManager;
class EvidenceManager;
class PotentialManager;
class MetaMove;
class HierarchicalCluster;
class Matrix;

class MetaLearner
{
	public:
		MetaLearner();
		~MetaLearner();
		void setShouldLoadCheckpoint(bool shouldLoadCheckpointVal);
		void setMaxFactorSize_Approx(int);
		void setBeta1(double aval);
		void setBeta_Motif(double bval);
		void setConvergenceThreshold(double);
		void setRestrictedList(const char*);
		int initEdgePriorMeta(const string& priorName, map<string,map<string,double>*>& graph, map<int,INTDBLMAP*>& edgePrior);
		int doCrossValidation(int);
		void setGlobalEvidenceManager(EvidenceManager*);
		void setVariableManager(VariableManager*);
		void setOutputDirName(const char*);
		void setClusteringThreshold(double);
		void setSpecificFold(int);
		int getPredictionError_CrossValid(int);
		void setRandom(bool);
		int readModuleMembership(const char*);
		int setDefaultModuleMembership();
		int setPriorGraph_All(const char*);

	private:

		// Checkpointing
		int writeCheckpointMetadata(int iter, bool notConvergedVal); 
		bool readCheckpointMetadata(int& iter, bool& notConvergedVal); 
		int writePLLScore(); 
		double loadInitPLLScore(); 
		int writeLastUpdate(); 
		int loadLastUpdate(); 
		int doTar(); 
		int readCheckpointModuleMembership(); 
		int populateGraphsFromFile(int maxFactorSizeApprox); 

		int start(int);
		int initEdgePriorMeta_All();
		double getInitPLLScore();
		double getInitPLLScore(int vId);
		double getPLLScore();
		double getEdgePrior(int tfID, int targetID);
		int clearFoldSpecData();
		int initEdgeSet();
		MetaMove* getNextMove(int currK, int vID);
		void makeMove(MetaMove* move, int currIteration);
		void getNewPLLScore(Variable* u, Variable* v, vector<int>& parentIDs, string& edgeKey, double& mbScore, double& scoreImprovement, Potential**);
		int dumpAllGraphs(int,int);
		int setPriorGraph(const char* aFName, map<string,map<string,double>*>& priorGraph);
		int redefineModules();
		double getModuleContribLogistic(string& tgtName, string& tfName);
		int initPhysicalDegree();
		int getEnrichedTFs(map<string,int>& tfSet,map<string,int>* genes,map<string,map<string,double>*>& edgeSet);
		void initCorrelationDistances();

		FactorGraph* factorGraph;
		PotentialManager* potManager;
		EvidenceManager* evidenceManager;
		VariableManager* varManager;

		bool shouldLoadCheckpoint;
		map <string,double> betamap;
		map <string,map<int,INTDBLMAP*>*> edgepriormap;
		map <string,map<string,map<string,double>*>*> priorgraphmap;
		int maxFactorSizeApprox;
		double convThreshold;
		char restrictedFName[1024];
		map<string,int> edgeMap;
		map<string,int> restrictedVarList;
		INTDBLMAP* currPLL;
		map<int,double> finalScores;
		char outputDirName[1024];
		map<string,double> edgePresenceProb;
		map<int,double> varNeighborhoodPrior;
		double beta1;
		double beta_motif;
		bool random;
		map<int,map<string,int>*> moduleGeneSet;
		map<string,int> geneModuleID;
		map<int,map<string,int>*> moduleIndegree;
		map<string,int> regulatorModuleOutdegree;
		double clusterThreshold;
		map<string,int> variableStatus;
		int currFold;
		char foldoutDirName[1024];
		int specificFold;
		HierarchicalCluster hc;
		Matrix* correlationDistances;
};
#endif
