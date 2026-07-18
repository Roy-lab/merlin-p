#include <iostream>
#include <math.h>
#include <algorithm>
#include "Potential.H"

#include "gsl/gsl_randist.h"

Potential::Potential(int factorID, double variance, double bias, unordered_map<int, double>& weights)
{
	this->factorID = factorID;
	this->variance = variance;
	this->bias = bias;

	this->weights.reserve(weights.size());

	// Copy over the weights
	for (auto iter = weights.begin(); iter != weights.end(); iter++) {
		this->weights.emplace_back(iter->first, iter->second);
	}

	// Sort the weights by var ID. Keeping these sorted makes it easier to find the intersection of parents between two potentials,
	// which we do when defining modules.
	sort(this->weights.begin(), this->weights.end(), [](const pair<int, double>& a, const pair<int, double>& b) {
		return a.first < b.first;
	});
}

vector<pair<int, double>>&
Potential::getWeights()
{
	return weights;
}

double
Potential::getExpectation(unordered_map<int, double>* evidenceSet)
{
	double mean=0;
	for(auto aIter = weights.begin(); aIter != weights.end(); aIter++)
	{
		if(evidenceSet->find(aIter->first)==evidenceSet->end())
		{
			cerr << "Fatal error! No variable assignment for " << aIter->first << endl;
			exit(-1);
		}
		double aval = (*evidenceSet)[aIter->first];
		mean += aval * aIter->second;
	}
	return mean + bias;
}

double
Potential::evaluateProbabilityDensity(unordered_map<int, double>* evidMap)
{
	if(evidMap->find(factorID)==evidMap->end()) {
		cerr <<"Fatal error! No variable assignment for " << factorID << endl;
		exit(-1);
	}

	double expectation = getExpectation(evidMap);
	double norm = sqrt(2 * PI * variance);
	double x = (*evidMap)[factorID];
	double dev = (x - expectation) * (x - expectation) / (2 * variance);
	double eval = exp(-1.0 * dev);
	double pval = eval / norm;
	return pval;
}
