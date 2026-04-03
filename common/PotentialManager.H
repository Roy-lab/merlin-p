/* The goal of this class is to take all the MarkovBlankets and the factors and
 *  the data and create all potentials involving valid pairings MarkovBlankets 
 *  and factors. We create all putative potentials upfront because we want to
 *  iterate over the data once and update counts of configurations in one go*/

#ifndef _POTENTIAL_MANAGER
#define _POTENTIAL_MANAGER

#include <vector>
#include <map>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_permutation.h>

using namespace std;

class Potential;
class EvidenceManager;
class SlimFactor;
class Matrix;

//This stores all the potentials
typedef map<int,Potential*> POTFUNC_MAP;
typedef map<int,Potential*>::iterator POTFUNC_MAP_ITER;

class PotentialManager
{
	public:
		PotentialManager();
		~PotentialManager();
		int setEvidenceManager(EvidenceManager*);
		int setRestrictedNeighborSet(map<int,Variable*>&);
		int setLambda(double);
		int init();
		int reset();

		int populatePotential(Potential*);

		double getLikelihood(SlimFactor*,VSET&);
		double getLikelihood(SlimFactor*,VSET&,map<int,int>&);
		int estimateConditionalPotential(SlimFactor*,VSET&,Potential**,STRDBLMAP&);
		int estimateCanonicalPotential(SlimFactor*, VSET&,INTINTMAP&,INTINTMAP&,map<int,SlimFactor*>&);
		int estimateCanonicalPotential_Abbeel(SlimFactor*, VSET&,INTINTMAP&,INTINTMAP&,map<int,SlimFactor*>&);
		int estimateCanonicalPotential_Approximate(SlimFactor*, VSET&,INTINTMAP&,INTINTMAP&,map<int,SlimFactor*>&);
		int resetPotFuncs();
		int estimateCanonicalPotential_Joint(SlimFactor*, VSET&,INTINTMAP&,INTINTMAP&,map<int,SlimFactor*>&);
		Potential* getPotential(int);

		double getSampleLikelihood(map<int,SlimFactor*>&, VSET&, INTINTMAP*);
		int getVariableSample(INTINTMAP&,VSET&,int,SlimFactor*,gsl_rng*);

		int clearJointEntropies();
		int setRandom(bool);

	private:
		int initData();
		int estimateCovariance(int,int);

		POTFUNC_MAP potFuncs;
		//This stores the joint potential for a joint configuration and then can
		//be reused over and over again to generate the conditional
		map<string,Potential*> canonicalPotFunc;
		EvidenceManager* evMgr;
		map<string,double> summaryStat;
		map<string,double> jointEntropies;
		map<int,Variable*> restrictedNeighborSet;
		gsl_matrix* ludecomp;
		gsl_permutation* perm;
		double lambda;
		bool randomData;
		Matrix* data;
		Matrix* meanMat;
		Matrix* covMat;
};
#endif
