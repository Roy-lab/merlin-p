#ifndef _EVIDENCE_SET
#define _EVIDENCE_SET

#include <vector>
#include <unordered_map>

using namespace std;

class EvidenceSet
{
    public:
        EvidenceSet(vector<double>& dataset, size_t sampleCount);

        int getSampleCount();

        int getVariableCount();

        // Exposed in the header so it can be inlined elsewhere.
        double getEvidenceAt(int varIndex, int sampleIndex)
        {
            return (*dataset)[varIndex * sampleCount + sampleIndex];
        }

    private:
        vector<double>* dataset;
        size_t sampleCount;
        size_t varCount;

};

#endif