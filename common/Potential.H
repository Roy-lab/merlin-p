#ifndef _POTENTIAL
#define _POTENTIAL

#include <fstream>
#include <unordered_map>
#include "CommonTypes.H"
#include "Matrix.H"
#include "EvidenceSet.H"

#define PI 3.1472

class Variable;

// Encapsulates the parameters of a conditional Gaussian PDF, and allows
// for evaluating probability density or expectation for a particular sample.
// The conditional Gaussian is represented as a univariate Gaussian with a mean
// evaluated as a linear combination of the conditioning variables.

class Potential
{
	public:
		Potential(int factorID, double variance, unordered_map<int, double>& weights);

		vector<pair<int, double>>& getWeights();

		double getVariance();

		int getFactorID();

		// Returns the expected value of the factor variable, given the values
		// of the conditioning variables.
		double getExpectation(EvidenceSet* evidenceSet, int sampleIndex);

	private:
		int factorID;
		double variance;
		vector<pair<int, double>> weights;
};
#endif
