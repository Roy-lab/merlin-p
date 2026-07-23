#ifndef _EVIDENCE_SOURCE
#define _EVIDENCE_SOURCE

class EvidenceSet;

class EvidenceSource
{
    public:
        enum class SetType
        {
            TrainingSet,
            TestSet
        };

        virtual ~EvidenceSource() = default;

        virtual void setupForFold(int foldIndex, int foldCount) = 0;

        virtual EvidenceSet* getEvidenceSet(SetType type) = 0;
};

#endif