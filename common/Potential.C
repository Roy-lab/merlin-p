#include <iostream>
#include <math.h>
#include "Evidence.H"
#include "Potential.H"

#include "gsl/gsl_randist.h"

Potential::Potential(int factorID, double variance, double bias, INTDBLMAP& weights)
{
	this->factorID = factorID;
	this->variance = variance;
	this->bias = bias;
	this->weights = weights;
}

INTDBLMAP&
Potential::getWeights()
{
	return weights;
}

double
Potential::getExpectation(map<int,Evidence*>* evidenceSet)
{
	double mean=0;
	for(INTDBLMAP_ITER aIter=weights.begin();aIter!=weights.end();aIter++)
	{
		if(evidenceSet->find(aIter->first)==evidenceSet->end())
		{
			cerr << "Fatal error! No variable assignment for " << aIter->first << endl;
			exit(-1);
		}
		Evidence* evid = (*evidenceSet)[aIter->first];
		double aval = evid->getEvidVal();
		mean += aval * aIter->second;
	}
	return mean + bias;
}

double 
Potential::evaluateProbabilityDensity(map<int,Evidence*>* evidMap)
{
	if(evidMap->find(factorID)==evidMap->end())
	{
		cerr <<"Fatal error! No variable assignment for " << factorID << endl;
		exit(-1);
	}

	double expectation = getExpectation(evidMap);
	double norm = sqrt(2 * PI * variance);
	Evidence* factorEvid = (*evidMap)[factorID];
	double x = factorEvid->getEvidVal();
	double dev = (x - expectation) * (x - expectation) / (2 * variance);
	double eval = exp(-1.0 * dev);
	double pval = eval / norm;
	return pval;
}
