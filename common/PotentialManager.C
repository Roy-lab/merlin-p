#include <iostream>
#include <cstring>
#include <math.h>

#include "CommonTypes.H"
#include "Error.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "PotentialManager.H"

PotentialManager::PotentialManager()
{
	ludecomp=NULL;
	perm=NULL;
	globalCovariances=nullptr;
}

PotentialManager::~PotentialManager()
{
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
	}
	if(globalCovariances!=nullptr)
	{
		delete globalCovariances;
	}
}

int
PotentialManager::init(EvidenceManager* evMgr, bool randomData, vector<int>& varIDs)
{
	if (globalCovariances != nullptr)
	{
		delete globalCovariances;
	}

	INTINTMAP& trainEvidSet = evMgr->getTrainingSet();
	EMAP* evidMap = evMgr->getEvidenceAt(trainEvidSet.begin()->first);
	int varCount = evidMap->size();
	int sampleCount = trainEvidSet.size();

	globalMeans.clear();

	globalCovariances = new Matrix(varCount, varCount);
	globalCovariances->setAllValues(-1);

	// Stores the deviations from the mean for each variable and sample.
	vector<double> deviations(varCount * sampleCount, 0);

	// Copy all the samples into the data matrix
	int sampleIndex = 0;
	for (INTINTMAP_ITER eIter = trainEvidSet.begin(); eIter != trainEvidSet.end(); eIter++)
	{
		EMAP* evidMap=NULL;
		if(randomData)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for (EMAP_ITER vIter = evidMap->begin(); vIter != evidMap->end(); vIter++)
		{
			int vId = vIter->first;
			Evidence* evid = vIter->second;
			double val = evid->getEvidVal();
			deviations[vId * sampleCount + sampleIndex] = val;
		}
		sampleIndex++;
	}

	// Done copying. Now we can go over data and get the means
	for (int i = 0; i < varCount; i++)
	{
		double sampleSum = 0;
		for(int j = 0; j < sampleCount; j++)
		{
			sampleSum += deviations[i * sampleCount + j];
		}
		globalMeans.push_back(sampleSum / sampleCount);
	}

	// Finally, use the means to pre-center the data
	for (int i = 0; i < trainEvidSet.size(); i++)
	{
		for (int j = 0; j < varCount; j++)
		{
			deviations[j * sampleCount + i] -= globalMeans[j];
		}
	}

	int norm = sampleCount - 1;

	// Set covariances along the diagonal.
	for (int i = 0; i < varCount; i++)
	{
		double ssd = 0.001;
		for (int j = 0; j < sampleCount; j++)
		{
			double dev = deviations[i * sampleCount + j];
			ssd += dev * dev;
		}
		globalCovariances->setValue(ssd / norm, i, i);
	}

	// Set covariances between regulators and all other variables.
	for (int i = 0; i < varIDs.size(); i++)
	{
		int regID = varIDs[i];
		for (int j = 0; j < varCount; j++)
		{
			if (regID == j)
			{
				continue;
			}

			double ssd = 0;
			for (int k = 0; k < sampleCount; k++)
			{
				double devI = deviations[regID * sampleCount + k];
				double devJ = deviations[j * sampleCount + k];
				ssd += devI * devJ;
			}

			double covariance = ssd / norm;
			globalCovariances->setValue(covariance, regID, j);
			globalCovariances->setValue(covariance, j, regID);
		}
	}

	ludecomp = gsl_matrix_alloc(MAXFACTORSIZE_ALLOC, MAXFACTORSIZE_ALLOC);
	perm = gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);

	return 0;
}

Potential*
PotentialManager::createPotential(int factorID)
{
	int varCount = globalMeans.size();
	double variance = globalCovariances->getValue(factorID, factorID);
	double bias = globalMeans[factorID];
	INTDBLMAP weights;
	return new Potential(factorID, variance, bias, weights);
}

double
PotentialManager::computeLL(int factorID, vector<int>& parentIDs, int sampleSize, Potential** newPot)
{
	double variance = globalCovariances->getValue(factorID, factorID);
	double bias = globalMeans[factorID];
	INTDBLMAP weights;

	int parentCount = parentIDs.size();

	// Start by collecting a matrix of all the covariances of the conditioning variables,
	// and the marginal variances of the conditioning variables.

	Matrix *parentCovariances = new Matrix(parentCount, parentCount);
	Matrix *parentMarginalVariances = new Matrix(1, parentCount);

	for (int i = 0; i < parentCount; i++)
	{
		int varAID = parentIDs[i];
		double factorCovariance = globalCovariances->getValue(factorID, varAID);
		parentMarginalVariances->setValue(factorCovariance, 0, i);

		for (int j = i; j < parentCount; j++)
		{
			int varBID = parentIDs[j];
			double covariance = globalCovariances->getValue(varAID, varBID);
			parentCovariances->setValue(covariance, i, j);
			parentCovariances->setValue(covariance, j, i);
		}
	}

	// Compute the final values for the variance of the conditional gaussian,
	// plus the regression parameters for the mean of the conditional guassian.

	Matrix* parentCovInverse = parentCovariances->invMatrix(ludecomp, perm);
	Matrix* prod = parentMarginalVariances->multiplyMatrix(parentCovInverse);

	for (int i = 0; i < parentCount; i++)
	{
		int vID = parentIDs[i];
		double aVal = prod->getValue(0, i);
		double bVal = parentMarginalVariances->getValue(0, i);
		double cVal = globalMeans[vID];
		weights[vID] = aVal;
		variance -= aVal * bVal;
		bias -= cVal * aVal;
	}

	delete prod;
	delete parentMarginalVariances;
	delete parentCovariances;
	delete parentCovInverse;

	if(variance < 1e-5)
	{
		variance = 1e-5;
	}

	// If the variance is invalid, then we don't want to attempt adding this edge,
	// so we should just bail out before computing the LL
	if(isnan(variance) || isinf(variance))
	{
		return -1;
	}

	// Now that the conditional Gaussian params are computed, we can create the potential.
	*newPot = new Potential(factorID, variance, bias, weights);

	// Finally, compute the conditional log likelihood.
	return -0.5 * ((sampleSize - 1) + sampleSize * log(2 * PI) + sampleSize * log(variance));
}
