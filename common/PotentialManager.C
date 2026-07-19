#include <iostream>
#include <cstring>
#include <math.h>
#include <unordered_map>

#include "CommonTypes.H"
#include "Error.H"
#include "Potential.H"
#include "EvidenceSource.H"
#include "EvidenceSet.H"
#include "PotentialManager.H"

PotentialManager::PotentialManager(EvidenceSource* source)
{
	globalCovariances = nullptr;
	evidenceSource = source;
}

PotentialManager::~PotentialManager()
{
	if (globalCovariances != nullptr) {
		delete globalCovariances;
	}
}

void PotentialManager::setupForFold(vector<int>& regIDs)
{
	if (globalCovariances != nullptr) {
		delete globalCovariances;
	}

	EvidenceSet* evidenceSet = evidenceSource->getEvidenceSet(EvidenceSource::SetType::TrainingSet);

	vector<double>* evidMap = evidenceSet->getEvidenceAt(0);
	int varCount = evidMap->size();
	int sampleCount = evidenceSet->getSize();

	globalMeans.clear();

	globalCovariances = new Matrix(varCount, varCount);
	globalCovariances->setAllValues(-1);

	// Stores the deviations from the mean for each variable and sample.
	vector<double> deviations(varCount * sampleCount, 0);

	// Copy all the samples into the data matrix
	for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++)
	{
		vector<double>* evidMap = evidenceSet->getEvidenceAt(sampleIndex);
		for (int vID = 0; vID < varCount; vID++)
		{
			double val = (*evidMap)[vID];
			deviations[vID * sampleCount + sampleIndex] = val;
		}
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
	for (int i = 0; i < evidenceSet->getSize(); i++)
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
	for (int i = 0; i < regIDs.size(); i++)
	{
		int regID = regIDs[i];
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
}

Potential* PotentialManager::createPotential(int factorID)
{
	int varCount = globalMeans.size();
	double variance = globalCovariances->getValue(factorID, factorID);
	double bias = globalMeans[factorID];
	unordered_map<int, double> weights;
	return new Potential(factorID, variance, bias, weights);
}

void PotentialManager::computeLLs(int factorID, int sampleSize, vector<int>& existingParents, vector<int>&candidateParents, unordered_map<int, double>&scores) {

	if (existingParents.size() == 0) {
		computeSingleParentLLs(factorID, sampleSize, candidateParents, scores);
		return;
	}

	// In order to compute likelihood, we need the variance of the factor conditioned upon existing parents and the candidate parent.
	// We use the following framing
	// A : Factor
	// B : existing parents
	// C : candidate parent
	// Var(A | B, C) = Var(A | B) - Cov(AC | B) * Var(C | B)^-1 * Cov(CA | B)
	// Var (A | B) = Var(A) - Cov(AB) * Var(B)^-1 * Cov(BA)
	// Var (C | B) = Var(C) - Cov(CB) * Var(B)^-1 * Cov(BC)
	// Cov (AC | B) = Cov(AC) - Cov(AB) * Var(B)^-1 * Cov(BC)

	// This allows us to precompute Var(A | B), Var(B)^-1, and Cov(AB) * Var(B)^-1, because they don't rely on C. Since we only add a single
	// candidate parent at a time, the per candidate work doesn't require any matrix inverts, only scalar and vector multiplication.

	int parentCount = existingParents.size();

	// Var(A)
	double factorVariance = globalCovariances->getValue(factorID, factorID);

	// Cov(AB)
	gsl_vector* existingParentMarginalVariances = gsl_vector_alloc(parentCount);

	// Cov(BB)
	gsl_matrix* existingParentCovariances = gsl_matrix_alloc(parentCount, parentCount);

	for (int i = 0; i < parentCount; i++) {
		int varAID = existingParents[i];
		double marginalCovariance = globalCovariances->getValue(factorID, varAID);
		gsl_vector_set(existingParentMarginalVariances, i, marginalCovariance);

		for (int j = i; j < parentCount; j++) {
			int varBID = existingParents[j];
			double covariance = globalCovariances->getValue(varAID, varBID);
			gsl_matrix_set(existingParentCovariances, i, j, covariance);
			gsl_matrix_set(existingParentCovariances, j, i, covariance);
		}
	}

	gsl_permutation* permutation = gsl_permutation_alloc(parentCount);

	int signum=0;
	gsl_linalg_LU_decomp(existingParentCovariances, permutation, &signum);

	// Var(B)^-1
	gsl_matrix* parentCovInverse = gsl_matrix_alloc(parentCount, parentCount);
	gsl_linalg_LU_invert(existingParentCovariances, permutation, parentCovInverse);

	gsl_matrix_free(existingParentCovariances);
	gsl_permutation_free(permutation);

	// Cov(AB) * Var(B)^-1
	gsl_vector* prod = gsl_vector_alloc(parentCount);
	gsl_blas_dgemv(CblasTrans, 1, parentCovInverse, existingParentMarginalVariances, 0, prod);

	// Cov(AB) * Var(B)^-1 * Cov(BA)
	double dot = 0.0;
	gsl_blas_ddot(prod, existingParentMarginalVariances, &dot);

	// Var(A|B)
	double existingConditionalVariance = factorVariance - dot;

	gsl_vector_free(existingParentMarginalVariances);

	for (int i = 0; i < candidateParents.size(); i++) {
		int candidateID = candidateParents[i];

		// Var(C)
		double candidateVariance = globalCovariances->getValue(candidateID, candidateID);

		// Cov(BC)
		gsl_vector* candidateMarginalVariances = gsl_vector_alloc(parentCount);

		for (int i = 0; i < parentCount; i++) {
			int varAID = existingParents[i];
			double marginalCovariance = globalCovariances->getValue(candidateID, varAID);
			gsl_vector_set(candidateMarginalVariances, i, marginalCovariance);
		}

		// Cov(BC) * Var(B)^-1
		gsl_vector* candidateProd = gsl_vector_alloc(parentCount);
		gsl_blas_dgemv(CblasTrans, 1, parentCovInverse, candidateMarginalVariances, 0, candidateProd);

		// Cov(BC) * Var(B)^-1 * Cov(CB)
		double dot = 0.0;
		gsl_blas_ddot(candidateProd, candidateMarginalVariances, &dot);

		// Var(C|B)
		double candidateConditionalVariance = candidateVariance - dot;

		if (candidateConditionalVariance < 1e-10) {
			continue;
		}

		// Cov(AC)
		double candidateFactorCovariance = globalCovariances->getValue(candidateID, factorID);

		// Cov(AB) * Var(B)^-1 * Cov(BC)
		dot = 0.0;
		gsl_blas_ddot(prod, candidateMarginalVariances, &dot);

		gsl_vector_free(candidateProd);
		gsl_vector_free(candidateMarginalVariances);

		// Cov(AC | B)
		double factorAndCandidateConditionalCov = candidateFactorCovariance - dot;

		// Var(A | B, C)
		double finalVariance = existingConditionalVariance - factorAndCandidateConditionalCov * factorAndCandidateConditionalCov / candidateConditionalVariance;

		if (finalVariance < 1e-5) {
			finalVariance = 1e-5;
		}

		if(isnan(finalVariance) || isinf(finalVariance)) {
			continue;
		}

		scores[candidateID] = -0.5 * ((sampleSize - 1) + sampleSize * log(2 * PI) + sampleSize * log(finalVariance));
	}

	gsl_matrix_free(parentCovInverse);
	gsl_vector_free(prod);
}

// Computes a gaussian log likelihood of factor id conditioned upon a single parent, for each candidate parent.
// The single parent case is broken out into a separate function because it can be done with fast scalar arithmetic.
void PotentialManager::computeSingleParentLLs(int factorID, int sampleSize, vector<int>&candidateParents, unordered_map<int, double>&scores) {

	double factorVariance = globalCovariances->getValue(factorID, factorID);

	for (int i = 0; i < candidateParents.size(); i++) {
		int candidateID = candidateParents[i];

		double candidateVariance = globalCovariances->getValue(candidateID, candidateID);

		if (candidateVariance < 1e-10) {
			continue;
		}

		double factorCandidateCov = globalCovariances->getValue(factorID, candidateID);
		double finalVariance = factorVariance - factorCandidateCov * factorCandidateCov / candidateVariance;

		if (finalVariance < 1e-5) {
			finalVariance = 1e-5;
		}

		if(isnan(finalVariance) || isinf(finalVariance)) {
			continue;
		}

		scores[candidateID] = -0.5 * ((sampleSize - 1) + sampleSize * log(2 * PI) + sampleSize * log(finalVariance));
	}
}

Potential* PotentialManager::createPotential(int factorID, vector<int>& parentIDs)
{
	int parentCount = parentIDs.size();
	double variance = globalCovariances->getValue(factorID, factorID);
	double bias = globalMeans[factorID];

	// Start by collecting a matrix of all the covariances of the conditioning variables,
	// and the marginal variances of the conditioning variables.

	gsl_matrix* parentCovariances = gsl_matrix_alloc(parentCount, parentCount);
	gsl_vector* parentMarginalVariances = gsl_vector_alloc(parentCount);

	for (int i = 0; i < parentCount; i++) {
		int varAID = parentIDs[i];
		double marginalCovariance = globalCovariances->getValue(factorID, varAID);
		gsl_vector_set(parentMarginalVariances, i, marginalCovariance);

		for (int j = i; j < parentCount; j++) {
			int varBID = parentIDs[j];
			double covariance = globalCovariances->getValue(varAID, varBID);
			gsl_matrix_set(parentCovariances, i, j, covariance);
			gsl_matrix_set(parentCovariances, j, i, covariance);
		}
	}

	// Compute the final values for the variance of the conditional gaussian,
	// plus the regression parameters for the mean of the conditional guassian.

	gsl_vector* x = gsl_vector_alloc(parentCount);

	gsl_permutation* permutation = gsl_permutation_alloc(parentCount);

	int signum = 0;
	gsl_linalg_LU_decomp(parentCovariances, permutation, &signum);

	gsl_linalg_LU_solve(parentCovariances, permutation, parentMarginalVariances, x);

	unordered_map<int, double> weights;
	for (int i = 0; i < parentCount; i++) {
		int vID = parentIDs[i];
		double aVal = gsl_vector_get(x, i);
		double bVal = gsl_vector_get(parentMarginalVariances, i);
		double cVal = globalMeans[vID];
		weights[vID] = aVal;
		variance -= aVal * bVal;
		bias -= cVal * aVal;
	}

	gsl_vector_free(x);
	gsl_vector_free(parentMarginalVariances);
	gsl_permutation_free(permutation);
	gsl_matrix_free(parentCovariances);

	return new Potential(factorID, variance, bias, weights);
}

double
PotentialManager::evaluateProbabilityDensity(Potential* potential, int sampleIndex, EvidenceSource::SetType type)
{
	// We can get the evidMap using the sample index
	EvidenceSet* evidenceSet = evidenceSource->getEvidenceSet(type);
	vector<double>* evidence = evidenceSet->getEvidenceAt(sampleIndex);

	int factorID = potential->getFactorID();
	double variance = potential->getVariance();

	double expectation = potential->getExpectation(evidence);
	double norm = sqrt(2 * PI * variance);
	double x = (*evidence)[factorID];
	double dev = (x - expectation) * (x - expectation) / (2 * variance);
	double eval = exp(-1.0 * dev);
	double pval = eval / norm;
	return pval;
}
