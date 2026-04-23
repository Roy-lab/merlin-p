#ifndef _FACTOR_MANAGER
#define _FACTOR_MANAGER

#include <vector>
#include <fstream>
using namespace std;
#include "gsl/gsl_randist.h"
#include "CommonTypes.H"
#define SAMPLE_CNT 1000

class Error;
class VariableManager;
class EvidenceManager;
class PotentialManager;
class SlimFactor;
class LatticeStructure;
class Graph;
class FactorGraph;

class FactorManager
{
	public:
		FactorManager();
		~FactorManager();
		int setVariableManager(VariableManager*);
		int setEvidenceManager(EvidenceManager*);
		int setPotentialManager(PotentialManager*);
		int setBaseInstantiation();
		int setBaseInstantiation_Variable();
		int setMaxFactorSize_Approx(int);
		int getMaxFactorSize_Approx();

		int setModelName(const char*);
		int setOutputDir(const char*);
		int setGraph(Graph*);
		int setBeamSize(int);
		int setPenalty(double);
		
		int allocateFactorSpace();
		int allocateFactorSpace_Graph();

		FactorGraph* createInitialFactorGraph();

		int generateSingletonFactors();
		int generateCanonicalFactors();
		int estimateCanonicalParameters(const char*);

		int learnStructure();

		int getFactorIndex(int*, int);
		int getFactorKey(int*, int,string&);

		int getMaximalClusters();
		//Here we populate the factSignMap
		int removeDupFactors();

		double getLikelihood_ChainRule();
		double getLikelihood_ChainRule(FactorGraph*);
		double getLikelihood();

		double getLikelihood_MCMC(FactorGraph*);
		double estimatePartitionFunction(FactorGraph*);
		int randInitSample(INTINTMAP* sample,gsl_rng* r);
		int getNextSample(FactorGraph* fg,INTINTMAP* currSample,INTINTMAP* nextSample,gsl_rng* r);

		bool checkMonotonicity();

		INTINTMAP* getSupersets(int);
		int getSupersets(int,int,INTINTMAP&);
		SlimFactor* getFactorAt(int);
		int deleteMySupersets(int factorID);
		int deleteFactor(SlimFactor*);
		int getActiveSubsetCnt(SlimFactor*,map<int,int>&);
		int updateRefCnt(int,int);
		int clearLattice(INTINTMAP&);

		int readRestrictedVarlist(const char*);
		map<int,Variable*>& getRestrictedVarlist();
	private:
		int makeMarkovBlanketConsistent(SlimFactor*);
		int produceClusters(double,int);
		int produceClusters_NoDup(double,int);

		int getClusterCnt(int,int,int);
		int initFactorSet();
		int populateFactorSet();
		int checkFactorIds();
		int qsort(int*,int,int);
		Error::ErrorCode estimateClusterProperties();

		int addToLattice(SlimFactor*,int**);
		int deleteFromLattice(int);

		double getOverlap(SlimFactor*,SlimFactor*,int&);

		int globalFactorID;	
		int maxFactorSize_Approx;
		VariableManager* vMgr;
		PotentialManager* potMgr;
		EvidenceManager* evMgr;

		//This is used to store the ids of factors which have only one variable
		//This will be used to compute the multi-information
		map<int,SlimFactor*> slimFactorSet;
		map<int,SlimFactor*> canonicalFactorSet;
		map<string,int> factorNameToIDMap;
		map<int,string> factorIDToNameMap;
		//Because we delete the old factors, by our pruning step, but
		//still need the multi-information associated with the factor
		//we store it here.
		map<string,double> delFactors_MI;
		
		int* pairIds;
		int* clusterIds;

		map<int,INTDBLMAP*> markovBlanketSet;
		INTINTMAP defaultInstantiation;
		LatticeStructure lattice;
		map<int,SlimFactor*> goodSlimFactors;

		char outputDir[1024];
		char modelName[1024];

		//Store a map of factor signature and the set of factors that have that signature.
		//Factor signature is defined by the union of variables in the factor and the variables in its Markov
		//blanket. We must not generate clusters for every possible factor and its Markov blanket
		//but per signature. The question we have to address here is the following:
		//Suppose we have two factors F and G with Markov blankets Mf and Mg respectively such that
		//F \Union Mf= G \Union Gf.  Should we pick F \Union Mf or should we pick G \Union Gf.
		//So which of these decompositions estimate the joint F\Union Mf better. We will pick
		//the decomposition that minimizes the interaction between the variables of Mf or in other
		//words the decomposition that lowers the conditional entropy the most. So if H(G|Mg) < H(F|Mf)
		//then we pick decomposition from G or we pick decomposition from F. Unfortunately, this may 
		//increase the interaction between the Markov blanket and the outside world.
		map<string,INTDBLMAP*> signFactGrpMap;
		map<int,string> factSignMap;
		map<string,int> groupFactorRepMap;
		map<string,double> mbSpecific_MI;

		//Only used to pairwise format of factors
		map<string,int> shownPairs;
		Graph* graph;
		double defProb;
		int defInstID;
		INTINTMAP defInstMap;
		int beamSize;
		double mbPenalty;
		map<int,Variable*> restrictedNeighborList;
};
#endif
