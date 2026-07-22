#ifndef _VALIDATION_LOGGER
#define _VALIDATION_LOGGER

#include <vector>

using namespace std;

class VariableSet;
class EvidenceSet;
class FactorGraph;
class PotentialSource;

class ValidationLogger
{
    public:
        void init(const char* outDirName, VariableSet* varSet);

        virtual int logValidationError(int foldID, EvidenceSet* testSet, FactorGraph* factorGraph, PotentialSource* potentialSource);

    protected:

        struct PredictionMetrics {
            double rmse;
            double normRmse;
            double r2;
            double cc;
        };

        PredictionMetrics computePredictionMetrics(const vector<double>& trueVals, const vector<double>& predVals);

        const char* outDirName;
        VariableSet* variableSet;
};

#endif