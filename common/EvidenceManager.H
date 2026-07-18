#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <vector>

using namespace std;

class Error;
class EvidenceSet;

class EvidenceManager
{
	public:

		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile(const char*);

		EvidenceSet getTrainingSet(int foldIndex, int foldCount);

		EvidenceSet getTestSet(int foldIndex, int foldCount);

	private:
		EvidenceSet getSet(int foldIndex, int foldCount, bool isTestSet);

		vector<vector<double>*> evidenceSet;

};
#endif
