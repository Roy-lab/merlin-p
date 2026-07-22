#ifndef _MERLIN_LOGGER
#define _MERLIN_LOGGER

#include "Logger.H"

class EvidenceSource;
class PotentialSource;

class MerlinLogger : public Logger
{
    public:
        MerlinLogger(EvidenceSource* evidenceSource, PotentialSource* potentialSource);

        int logValidationError(int foldID, FactorGraph* factorGraph) override;

        void logVariableMarkovBlankets(int foldID, int maxFactorSize, FactorGraph* factorGraph) override;

    private:
        EvidenceSource* evidenceSource;
        PotentialSource* potentialSource;
};

#endif