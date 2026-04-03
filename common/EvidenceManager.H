#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <map>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

#include "gsl/gsl_randist.h"
using namespace std;

class Error;
class Evidence;

//EMAP stores the map of variable ids and their evidence. This corresponds
//to one line in the datafile
typedef map<int,Evidence*> EMAP;
typedef map<int,Evidence*>::iterator EMAP_ITER;

typedef vector<EMAP*> EVIDENCE_SET;

class VariableManager;
class EvidenceManager
{
	public:
		EvidenceManager();

		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile(const char*);

		//Get the evidence set at the specified id
		EMAP* getEvidenceAt(int);
		EMAP* getRandomEvidenceAt(int);

		int randomizeEvidence(gsl_rng*, VariableManager*);
		int setFoldCnt(int);
		int setPreRandomizeSplit();
		int splitData(int);
		INTINTMAP& getTrainingSet();
		INTINTMAP& getTestSet();
		int populateRandIntegers(gsl_rng*, vector<int>&,int,int);

	private:
		int populateRandIntegers(gsl_rng*, int*,int);
		int populateRandIntegers(gsl_rng* r, int* randInds, INTINTMAP& populateFrom, int size);

		EVIDENCE_SET evidenceSet;
		EVIDENCE_SET randEvidenceSet;
		INTINTMAP trainIndex;
		INTINTMAP testIndex;
		int foldCnt;
		bool preRandomizeSplit;
		int randseed;
};
#endif
