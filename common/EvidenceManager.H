#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <vector>
#include "EvidenceSource.H"

using namespace std;

class Error;

class EvidenceManager : public EvidenceSource
{
	public:

		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile(const char*);

		EvidenceSet getTrainingSet(int foldIndex, int foldCount) override;

		EvidenceSet getTestSet(int foldIndex, int foldCount) override;

	private:
		EvidenceSet getSet(int foldIndex, int foldCount, bool isTestSet);

		vector<vector<double>*> evidenceSet;

};
#endif
