#include "EvidenceSet.H"

EvidenceSet::EvidenceSet(vector<vector<double>*>& inDataset, vector<int> inIndices)
{
    dataset = &inDataset;
    indices = inIndices;
}

vector<double>*
EvidenceSet::getEvidenceAt(int index)
{
    int globalIndex = indices[index];
    return (*dataset)[globalIndex];
}

int
EvidenceSet::getSize()
{
    return indices.size();
}