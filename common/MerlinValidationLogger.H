#ifndef _MERLIN_VALIDATION_LOGGER
#define _MERLIN_VALIDATION_LOGGER

#include "ValidationLogger.H"

class EvidenceSource;
class PotentialSource;

class MerlinValidationLogger : public ValidationLogger
{
    public:
        MerlinValidationLogger(EvidenceSource* evidenceSource, PotentialSource* potentialSource);

        int logValidationError(int foldID, FactorGraph* factorGraph) override;

    private:
        EvidenceSource* evidenceSource;
        PotentialSource* potentialSource;
};

#endif