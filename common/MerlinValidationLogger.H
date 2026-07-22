#ifndef _MERLIN_VALIDATION_LOGGER
#define _MERLIN_VALIDATION_LOGGER

#include "ValidationLogger.H"

class EvidenceSource;

class MerlinValidationLogger : public ValidationLogger
{
    public:
        MerlinValidationLogger(EvidenceSource* evidenceSource);

        int logValidationError(int foldID, FactorGraph* factorGraph, PotentialSource* potentialSource) override;

    private:
        EvidenceSource* evidenceSource;
};

#endif