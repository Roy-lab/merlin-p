#include "EvidenceSet.H"

EvidenceSet::EvidenceSet(vector<double>& inDataset, size_t inSampleCount)
{
    dataset = &inDataset;
    sampleCount = inSampleCount;
    varCount = inDataset.size() == 0 ? 0 : inDataset.size() / inSampleCount;
}

int
EvidenceSet::getSampleCount()
{
    return sampleCount;
}

int
EvidenceSet::getVariableCount()
{
    return varCount;
}
