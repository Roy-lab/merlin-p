#ifndef _META_LEARNER_
#define _META_LEARNER_

#include <map>
#include "CommonTypes.H"
#include "HierarchicalCluster.H"

using namespace std;

//This is gamma map for each variable within each datapoint
typedef map<int,INTDBLMAP*> DATA_GAMMA_MAP;
typedef map<int,INTDBLMAP*>::iterator DATA_GAMMA_MAP_ITER;

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
	int setMaxFactorSize_Approx(int);
	int setBeta1(double aval);
	int setBeta_Motif(double bval);
	int setLambda(double bval);
	int setConvergenceThreshold(double);
	int setRestrictedList(const char*);
	int setPreRandomizeSplit();
	int initEdgePriorMeta(const string& priorName, map<string,map<string,double>*>& graph, map<int,INTDBLMAP*>& edgePrior);
	int doCrossValidation(int);
	int setGlobalEvidenceManager(EvidenceManager*);
	int setVariableManager(VariableManager*);
	int setOutputDirName(const char*);
	int setClusteringThreshold(double);
	int setSpecificFold(int);
	int getPredictionError_CrossValid(int);
	int setRandom(bool);
	int readModuleMembership(const char*);
	int setDefaultModuleMembership();
	int setPriorGraph_All(const char*);

	private:
	int start(int);
	int initEdgePriorMeta_All();
	map <string,double> betamap;
	map <string,map<int,INTDBLMAP*>*> edgepriormap;
	map <string,map<string,map<string,double>*>*> priorgraphmap;

	double getInitPLLScore();
	double getInitPLLScore(int vId);
	double getPLLScore();
	double getEdgePrior(int tfID, int targetID);
	int clearFoldSpecData();
	int initEdgeSet();

	MetaMove* getNextMove(int currK, int vID);
	void makeMove(MetaMove* move, int currIteration);
	void getNewPLLScore(Variable* u, Variable* v, vector<int>& parentIDs, string& edgeKey, double& mbScore, double& scoreImprovement, Potential**);

	int dumpAllGraphs(int,int,int);

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

	int maxFactorSizeApprox;
	double convThreshold;
	char restrictedFName[1024];
	map<string,int> edgeMap;
	map<string,int> restrictedVarList;
	INTDBLMAP* currPLL;
	map<int,double> finalScores;
	bool preRandomizeSplit;
	gsl_rng* rnd;
	char outputDirName[1024];
	map<string,double> edgePresenceProb;
	map<int,double> varNeighborhoodPrior;
	double beta1;
	double beta_motif;
	double lambda;
	bool random;
	bool randomData;
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
