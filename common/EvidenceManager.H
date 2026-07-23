#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <vector>
#include "EvidenceSource.H"

using namespace std;

class Error;
class EvidenceSet;

class EvidenceManager : public EvidenceSource
{
	public:
		EvidenceManager();
		~EvidenceManager();

		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile(const char*);

		void setupForFold(int foldIndex, int foldCount) override;

		EvidenceSet* getEvidenceSet(SetType type) override;

	private:
		EvidenceSet* createSet(int foldIndex, int foldCount, SetType type);

		vector<vector<double>*> evidenceSet;
		EvidenceSet* trainingSet;
		EvidenceSet* testSet;

};
#endif
