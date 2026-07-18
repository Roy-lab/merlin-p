#ifndef _EVIDENCE_SOURCE
#define _EVIDENCE_SOURCE

class EvidenceSet;

class EvidenceSource
{
    public:
        virtual ~EvidenceSource() = default;

        virtual EvidenceSet getTrainingSet(int foldIndex, int foldCount) = 0;

		virtual EvidenceSet getTestSet(int foldIndex, int foldCount) = 0;
};

#endif