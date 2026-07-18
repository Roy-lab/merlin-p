#ifndef _EVIDENCE_SET
#define _EVIDENCE_SET

#include <vector>
#include <unordered_map>

using namespace std;

class EvidenceSet
{
    public:
        EvidenceSet(vector<vector<double>*>& dataset, vector<int> indices);

        vector<double>* getEvidenceAt(int index);

        int getSize();

    private:
        vector<vector<double>*>* dataset;
        vector<int> indices;

};

#endif