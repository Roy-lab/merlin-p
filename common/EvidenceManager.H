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

		void setExpressionFile(const char*);

		void setupForFold(int foldIndex, int foldCount) override;

		EvidenceSet* getEvidenceSet(SetType type) override;

	private:
		void loadData(int foldIndex, int foldCount);
		void centerData(vector<double>& data, size_t sampleCount);

		vector<double> trainData;
		vector<double> testData;
		size_t trainSampleCount;
		size_t testSampleCount;
		EvidenceSet* trainingSet;
		EvidenceSet* testSet;
		string expressionFile;
};
#endif
