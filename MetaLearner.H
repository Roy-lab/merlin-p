#ifndef _META_LEARNER_
#define _META_LEARNER_

#include <map>
#include "CommonTypes.H"
#include "HierarchicalCluster.H"

using namespace std;

//This is gamma map for each variable within each datapoint
typedef map<int,INTDBLMAP*> DATA_GAMMA_MAP;
typedef map<int,INTDBLMAP*>::iterator DATA_GAMMA_MAP_ITER;

class FactorManager;
class VariableManager;
class EvidenceManager;
class PotentialManager;
class MetaMove;
class Graph;
class HierarchicalCluster;

class MetaLearner
{
	public:
	MetaLearner();
	~MetaLearner();
	int setInputFName(const char*);
	int setMaxFactorSize_Approx(int);
	int setPenalty(double);
	int setBeta1(double aval);
	int setBeta_Motif(double bval);
	int setLambda(double bval);
	int setConvergenceThreshold(double);
	int setRestrictedList(const char*);
	int setPreRandomizeSplit();
	int setTargetList(const char*);
	int setNoPrune(bool);
	int initPartitions();
	int initEdgePriorMeta(map<string,map<string,double>*>& graph, map<int,INTDBLMAP*>& edgePrior);
	int doCrossValidation(int);
	int start(int);
	int setGlobalEvidenceManager(EvidenceManager*);
	int setVariableManager(VariableManager*);
	int setOutputDirName(const char*);
	int setClusteringThreshold(double);
	int setSpecificFold(int);
	int getPredictionError_CrossValid(int);
	int setRandom(bool);
	int readModuleMembership(const char*);
	int setDefaultModuleMembership();

	//int setMotifGraph(const char*);
	//int setChIPGraph(const char*);
	//int setBeta_ChIP(double bval);
	int setPriorGraph_All(const char*);
	private:
	
	//map<int,INTDBLMAP*> edgePriors_Motif;
	//map<int,INTDBLMAP*> edgePriors_ChIP;
	//map<string,map<string,double>*> priorGraph_Motif;
	//map<string,map<string,double>*> priorGraph_ChIP;
	//double beta_chip;
	//int initEdgePriorMeta_Motif();
	//int initEdgePriorMeta_ChIP();
	int initEdgePriorMeta_All();
	map <string,double> betamap;
	map <string,map<int,INTDBLMAP*>*> edgepriormap;
	map <string,map<string,map<string,double>*>*> priorgraphmap;
		
	double getInitPLLScore();
	double getPLLScore();
	double getInitPrior();
	double getPriorChange();
	double getEdgePrior(int tfID, int targetID);
	int clearFoldSpecData();
	int initEdgeSet(bool);
	int initCondsetMap_Nopool();
	INTINTMAP* getConditionSet(int cind);

	int collectMoves(int);
	int collectMoves(int,int);
	int sortMoves();
	int makeMoves();
	int attemptMove(MetaMove*, INTINTMAP*);
	int getNewPLLScore(int cid, INTINTMAP& conditionSet, Variable* u, Variable* v, string& edgeKey, double& mbScore, double& scoreImprovement,Potential**);
	double getNewPLLScore_Condition(int csetId, int vId, Potential* newPot);
	double getNewPLLScore_Condition_Tracetrick(int csetId, int vId, int uID, Potential* newPot);

	int genCondSetKey(INTINTMAP& condSet, string& aKey);

	int dumpAllGraphs(int,int,int);
	int getEdgeVars(string&, char*, char*);
	bool checkMBSize(int,int,int);

	int setPriorGraph(const char* aFName, map<string,map<string,double>*>& priorGraph);
	int redefineModules();
	double getModuleContrib(string& tgtName, string& tfName);
	double getModuleContribLogistic(string& tgtName, string& tfName);
	int initPhysicalDegree();
	int getEnrichedTFs(map<string,int>& tfSet,map<string,int>* genes,map<string,map<string,double>*>& edgeSet);
	double getModuleWideScoreImprovement(int cid,INTINTMAP& conditionSet,Variable* u, Variable* v,map<string,int>* moduleMembers);
	
	FactorManager* factorManager;
	FactorGraph* factorGraph;
	PotentialManager* potManager;

	map<int,EvidenceManager*> evMgrSet;

	EvidenceManager* globalEvMgr;

	VariableManager* varManager;
	int maxFactorSizeApprox;
	int maxMBCntGreedy;
	char inputFName[1024];
	double penalty;
	double convThreshold;
	char restrictedFName[1024];
	map<string,INTINTMAP*> edgeConditionMap;
	FactorGraph* fg;
	vector<MetaMove*> moveSet;
	map<string,double> edgeValue;
	map<int,INTINTMAP*> condsetMap;
	map<string,int> condsetKeyIDMap;
	map<int,string> condsetIDKeyMap;
	map<string,int> restrictedVarList;
	INTDBLMAP* currPLL;
	char trueGraphFName[1024];
	map<int,double> finalScores;
	bool preRandomizeSplit;
	gsl_rng* rnd;
	char outputDirName[1024];
	map<string,double> edgePresenceProb;
	map<int,double> varNeighborhoodPrior;
	//These are the edges that got added
	map<string,int> edgeUpdates;
	double beta1;
	double beta_motif;
	double lambda;
	map<string,int> reqdTargetList;
	bool random;
	bool randomData;
	bool noprune;
	map<int,map<string,int>*> moduleGeneSet;
	map<string,int> geneModuleID;
	map<int,map<string,int>*> moduleIndegree;
	map<string,int> regulatorModuleOutdegree;
	map<int,int> idVidMap;
	int globalID;
	double clusterThreshold;
	map<string,int> variableStatus;
	int currFold;
	char foldoutDirName[1024];
	int specificFold;
	HierarchicalCluster hc;
};
#endif
