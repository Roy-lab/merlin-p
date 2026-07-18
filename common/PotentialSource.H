#ifndef _POTENTIAL_SOURCE
#define _POTENTIAL_SOURCE

#include <vector>
#include <unordered_map>

using namespace std;

class EvidenceSet;
class Potential;

class PotentialSource
{
    public:
        virtual ~PotentialSource() = default;

        // Load the provided data, and pre-compute regulator covariances
        virtual void setupForFold(int foldID, int foldCount, vector<int>& regIDs) = 0;

        // Creates a Potential of a single variable.
		virtual Potential* createPotential(int factorID) = 0;

        // Creates a Potential conditioned on a set of parents.
		virtual Potential* createPotential(int factorID, vector<int>& parentIDs) = 0;

        // Computes log likelihood for the factor conditioned upon its existing parents and each candidate parent.
		virtual void computeLLs(int factorID, int sampleSize, vector<int>& existingParents, vector<int>&candidateParents, unordered_map<int, double>&scores) = 0;
};

#endif