#ifndef _POTENTIAL_SOURCE
#define _POTENTIAL_SOURCE

#include <vector>
#include <unordered_map>
#include "EvidenceSource.H"

using namespace std;

class Potential;

class PotentialSource
{
    public:
        virtual ~PotentialSource() = default;

        // Load the provided data, and pre-compute regulator covariances
        virtual void setupForFold(vector<int>& regIDs) = 0;

        // Creates a Potential of a single variable.
        virtual Potential* createPotential(int factorID) = 0;

        // Creates a Potential conditioned on a set of parents.
        virtual Potential* createPotential(int factorID, vector<int>& parentIDs) = 0;

        // Computes log likelihood for the factor conditioned upon its existing parents and each candidate parent.
        virtual void computeLLs(int factorID, vector<int>& existingParents, vector<int>&candidateParents, unordered_map<int, double>&scores) = 0;

        // Evaluates pseudo log likelihood for the factor as univariate gaussian, for the initial edgeless graph.
        virtual double getInitPLLScore(Potential* potential) = 0;

        // Evaluates conditional probability density of the Potential's factor, given its parents, at the sample index.
        virtual double evaluateProbabilityDensity(Potential* potential, int sampleIndex, EvidenceSource::SetType type) = 0;
};

#endif