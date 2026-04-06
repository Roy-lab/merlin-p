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
	data=NULL;
	meanMat=NULL;
	covMat=NULL;
	ludecomp=NULL;
	perm=NULL;
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
	if(data!=NULL)
	{
		delete data;
	}
	if(meanMat!=NULL)
	{
		delete meanMat;
	}
	if(covMat!=NULL)
	{
		delete covMat;
	}
}

int
PotentialManager::init(EvidenceManager* evMgr, bool randomData)
{
	reset();

	INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
	EMAP* evidMap=evMgr->getEvidenceAt(trainEvidSet.begin()->first);
	int varCnt=evidMap->size();

	// data is the data matrix which will have the variable by sample information
	data=new Matrix(varCnt,trainEvidSet.size());
	meanMat=new Matrix(varCnt,1);
	meanMat->setAllValues(0);
	covMat=new Matrix(varCnt,varCnt);
	covMat->setAllValues(-1);

	// Copy all the samples into the data matrix
	int sampleIndex = 0;
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
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
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double val=evid->getEvidVal();
			data->setValue(val,vId,sampleIndex);
		}
		sampleIndex++;
	}

	// Done copying. Now we can go over the rows of data and get the means
	for(int i=0;i<varCnt;i++)
	{
		double sampleSum=0;
		for(int j=0;j<data->getColCnt();j++)
		{
			sampleSum += data->getValue(i,j);
		}
		double sampleSize=(double) data->getColCnt();
		meanMat->setValue(sampleSum/sampleSize,i,0);
	}

	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}

int
PotentialManager::reset()
{
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
		ludecomp=NULL;
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
		perm=NULL;
	}
	if(data!=NULL)
	{
		delete data;
		data=NULL;
	}
	if(meanMat!=NULL)
	{
		delete meanMat;
		meanMat=NULL;
	}
	if(covMat!=NULL)
	{
		delete covMat;
		covMat=NULL;
	}
	return 0;
}

double
PotentialManager::getCovariance(int uId, int vId)
{
	double cov = covMat->getValue(uId, vId);
	if(cov == -1)
	{
		estimateCovariance(uId, vId);
		cov = covMat->getValue(uId, vId);
	}
	return cov;
}

void
PotentialManager::estimateCovariance(int uId, int vId)
{
	double ssd=0;
	for(int i=0;i<data->getColCnt();i++)
	{
		double vval=data->getValue(vId,i);
		double vmean=meanMat->getValue(vId,0);
		double uval=data->getValue(uId,i);
		double umean=meanMat->getValue(uId,0);
		ssd += (vval-vmean)*(uval-umean);
	}
	if(uId==vId)
	{
		ssd += 0.001;
	}
	double var = ssd/((double)(data->getColCnt()-1));
	covMat->setValue(var,uId,vId);
	covMat->setValue(var,vId,uId);
}

Potential*
PotentialManager::createPotential(int factorID)
{
	double variance = getCovariance(factorID, factorID);
	double bias = meanMat->getValue(factorID, 0);
	INTDBLMAP weights;
	return new Potential(factorID, variance, bias, weights);
}

double
PotentialManager::computeLL(int factorID, vector<int>& parentIDs, int sampleSize, Potential** newPot)
{
	double variance = getCovariance(factorID, factorID);
	double bias = meanMat->getValue(factorID, 0);
	INTDBLMAP weights;

	int parentCount = parentIDs.size();
	int varCount = parentCount + 1;

	// Start by collecting a matrix of the covariances of all the variables in the
	// joint gaussian, a matrix of all the covariances of the conditioning variables,
	// and the marginal variances of the conditioning variables.

	Matrix *covariances = new Matrix(varCount, varCount);
	Matrix *parentCovariances = new Matrix(parentCount, parentCount);
	Matrix *parentMarginalVariances = new Matrix(1, parentCount);

	covariances->setValue(variance, 0, 0);

	for (int i = 0; i < parentCount; i++)
	{
		int varAID = parentIDs[i];
		double factorCovariance = getCovariance(factorID, varAID);
		parentMarginalVariances->setValue(factorCovariance, 0, i);
		covariances->setValue(factorCovariance, 0, i+1);
		covariances->setValue(factorCovariance, i+1, 0);

		for (int j = 0; j < parentCount; j++)
		{
			int varBID = parentIDs[j];
			double covariance = getCovariance(varAID, varBID);
			parentCovariances->setValue(covariance, i, j);
			covariances->setValue(covariance, i+1, j+1);
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
		double cVal = meanMat->getValue(vID, 0);
		weights[vID] = aVal;
		variance -= aVal * bVal;
		bias -= cVal * aVal;
	}

	delete prod;
	delete parentMarginalVariances;

	if(variance < 1e-5)
	{
		variance = 1e-5;
	}

	// If the variance is invalid, then we don't want to attempt adding this edge,
	// so we should just bail out before computing the LL
	if(isnan(variance) || isinf(variance))
	{
		delete covariances;
		delete parentCovariances;
		delete parentCovInverse;
		return -1;
	}

	// Now that the conditional Gaussian params are computed, we can create the potential.
	*newPot = new Potential(factorID, variance, bias, weights);

	// Finally, compute the log likelihood of the conditional Gaussian, by subtracting the
	// LL of the joint Gaussian of the conditioning variables from that of all variables.

	Matrix *covInverse = covariances->invMatrix(ludecomp, perm);
	double determinant = covariances->detMatrix(ludecomp, perm);
	double parentDeterminant = parentCovariances->detMatrix(ludecomp, perm);

	double jointLL = computeJointLL(covariances, covInverse, determinant, sampleSize);
	double jointLLParents = computeJointLL(parentCovariances, parentCovInverse, parentDeterminant, sampleSize);

	delete covariances;
	delete covInverse;
	delete parentCovariances;
	delete parentCovInverse;

	return jointLL - jointLLParents;
}

double
PotentialManager::computeJointLL(Matrix* covariances, Matrix* inverse, double determinant, int sampleSize)
{
	int dim=covariances->getRowCnt();

	// The sum of squared deviations is the covariance matrix * (num samples - 1)
	Matrix sos = Matrix(dim, dim);
	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim; j++)
		{
			sos.setValue(covariances->getValue(i, j) * (sampleSize - 1), i, j);
		}
	}

	// Take trace(inverse covariance * sos)
	Matrix* m = inverse->multiplyMatrix(&sos);
	double trace = 0;
	for(int i = 0; i < dim; i++)
	{
		trace += m->getValue(i, i);
	}
	delete m;

	double constant = sampleSize * log(determinant) + sampleSize * dim * log(2 * PI);
	return -0.5 * (trace + constant);
}
