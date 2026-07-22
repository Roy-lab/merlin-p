#ifndef _LOGGER
#define _LOGGER

#include <vector>

using namespace std;

class VariableSet;
class EvidenceSet;
class FactorGraph;

class Logger
{
    public:
        void setOutDirName(const char* outDirName);

        void setVariableSet(VariableSet* varSet);

        virtual int logValidationError(int foldID, FactorGraph* factorGraph) = 0;

        virtual void logVariableMarkovBlankets(int foldID, int maxFactorSize, FactorGraph* factorGraph) = 0;

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