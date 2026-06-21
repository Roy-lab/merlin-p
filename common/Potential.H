#ifndef _POTENTIAL
#define _POTENTIAL

#include <fstream>
#include <unordered_map>
#include "gsl/gsl_randist.h"
#include "CommonTypes.H"
#include "Matrix.H"

#define PI 3.1472

class Evidence;
class Variable;

// Encapsulates the parameters of a conditional Gaussian PDF, and allows
// for evaluating probability density or expectation for a particular sample.
// The conditional Gaussian is represented as a univariate Gaussian with a mean
// evaluated as a linear combination of the conditioning variables.

class Potential
{
	public:
		Potential(int factorID, double variance, double bias, unordered_map<int, double>& weights);

		vector<pair<int, double>>& getWeights();

		// Returns the expected value of the factor variable, given the values
		// of the conditioning variables.
		double getExpectation(unordered_map<int, Evidence*>*);

		// Evaluates the conditional probability density at the provided factor
		// variable value, given the provided values of the conditioning variables.
		double evaluateProbabilityDensity(unordered_map<int, Evidence*>*);

	private:
		int factorID;
		double variance;
		double bias;
		vector<pair<int, double>> weights;
};
#endif
