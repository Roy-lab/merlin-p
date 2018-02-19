#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <map>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

#include "Matrix.H"

#include "gsl/gsl_randist.h"
using namespace std;

class Error;
class Evidence;


class VariableManager;
class EvidenceManager
{
	public:
		EvidenceManager();
		~EvidenceManager();
		//This associates the variable information with the variable id
		int setVariableManager(VariableManager*);
		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile(const char*);
		Error::ErrorCode loadEvidenceFromFile_Continuous(const char*);
		int randomizeEvidence(gsl_rng*);
		int getNumberOfEvidences();
		//Get the evidence set at the specified id
		EMAP* getEvidenceAt(int);
		EMAP* getRandomEvidenceAt(int);
		int addEvidence(EMAP*);
		//This is used during EM, updating missing/hidden values using 
		//expected counts. eSetId corresponds to the data entry, vId is the
		//variable whose evidence we are updating and evidData is the set
		//of values and their evidence values (soft evidence)
		int addToEvidence(int eSetId,int vId, INTDBLMAP& evidData);	
		int dumpEvidenceSet(ostream&);
		int dumpSummaryStat(ostream&);
		int getMLSettings(ostream&);

		int setFoldCnt(int);
		int generateValidationSet(const char*,int vSetSize,gsl_rng* r);
		int setPreRandomizeSplit();
		int splitData(int);
		INTINTMAP& getTrainingSet();
		INTINTMAP& getTestSet();
		INTINTMAP& getValidationSet();
		int populateRandIntegers(gsl_rng*, INTINTMAP&,int,int);
		int populateRandIntegers(gsl_rng*, vector<int>&,int,int);
		int standardizeData();
		int partitionData(int,map<int,EvidenceManager*>&,int&,map<int,INTINTMAP*>&);
		int estimateCovariance(int uId, int vId, double& ucov, double& vcov, double& uvcov);
	private:
		Matrix* dataMat;
		int populateEvidence(Evidence**,const char*);
		int populateRandIntegers(gsl_rng*, int*,int);
		int populateRandIntegers(gsl_rng* r, int* randInds, INTINTMAP& populateFrom, int size);
		EVIDENCE_SET evidenceSet;
		EVIDENCE_SET randEvidenceSet;
		VariableManager* vMgr;
		INTINTMAP trainIndex;
		INTINTMAP testIndex;
		INTINTMAP validationIndex;
		int foldCnt;
		bool preRandomizeSplit;
		int randseed;
		int updateDataMat();
};
#endif
