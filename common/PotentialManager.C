#include <iostream>
#include <cstring>
#include <math.h>

#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

PotentialManager::PotentialManager()
{
	data=NULL;
	meanMat=NULL;
	covMat=NULL;
	ludecomp=NULL;
	perm=NULL;
	lambda=0;
	randomData=false;
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
PotentialManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int 
PotentialManager::setRestrictedNeighborSet(map<int,Variable*>& rSet)
{
	for(map<int,Variable*>::iterator vIter=rSet.begin();vIter!=rSet.end();vIter++)
	{
		restrictedNeighborSet[vIter->first]=vIter->second;
	}
	return 0;
}

int
PotentialManager::setLambda(double lVal)
{
	lambda=lVal;
	return 0;
}

int
PotentialManager::setRandom(bool flag)
{
	randomData=flag;
	return 0;
}

int
PotentialManager::init()
{
	initData();
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

int
PotentialManager::initData()
{
	INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
	EMAP* evidMap=evMgr->getEvidenceAt(trainEvidSet.begin()->first);
	int varCnt=evidMap->size();

	// data is the data matrix which will have the variable by sample information
	if(data==NULL)
	{
		data=new Matrix(varCnt,trainEvidSet.size());
		meanMat=new Matrix(varCnt,1);
		meanMat->setAllValues(0);
		covMat=new Matrix(varCnt,varCnt);
		covMat->setAllValues(-1);
	}

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

	return 0;
}

int
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
	return 0;
}

int
PotentialManager::populatePotential(Potential* aPot)
{
	VSET& potVars=aPot->getAssocVariables();
	for(VSET_ITER vIter=potVars.begin();vIter!=potVars.end(); vIter++)
	{
		double mean=meanMat->getValue(vIter->first,0);
		aPot->updateMean(vIter->first,mean);

		for(VSET_ITER uIter=vIter;uIter!=potVars.end();uIter++)
		{
			double cval=covMat->getValue(vIter->first,uIter->first);
			if(cval==-1)
			{
				estimateCovariance(uIter->first,vIter->first);
				cval=covMat->getValue(vIter->first,uIter->first);
			}
			aPot->updateCovariance(vIter->first,uIter->first,cval);
			aPot->updateCovariance(uIter->first,vIter->first,cval);
		}
	}
	aPot->makeValidJPD(ludecomp,perm);
	return 0;
}

double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet,map<int,int>& visitedVertices )
{
	double dll=0;
	cout <<"Not implemented" << endl;
	return dll;
}

int
PotentialManager::estimateConditionalPotential(SlimFactor* sFactor,VSET& varSet,Potential** pot, STRDBLMAP& counts)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential_Abbeel(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential_Approximate(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int
PotentialManager::resetPotFuncs()
{
	for(map<int,Potential*>::iterator pIter=potFuncs.begin();pIter!=potFuncs.end();pIter++)
	{
		delete pIter->second;
	}
	potFuncs.clear();
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential_Joint(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

Potential*
PotentialManager::getPotential(int fId)
{
	if(potFuncs.find(fId)==potFuncs.end())
	{
		return NULL;
	}
	return potFuncs[fId];
}

double 
PotentialManager::getSampleLikelihood(map<int,SlimFactor*>& factorSet, VSET& varSet, INTINTMAP* sample)
{
	double sampleLL=0;
	cout <<"Not implemented " <<endl;
	return sampleLL;
}

int 
PotentialManager::getVariableSample(INTINTMAP& jointConf,VSET& varSet,int vId,SlimFactor* sFactor, gsl_rng* r)
{
	Potential* pot=NULL;
	if(potFuncs.find(vId)==potFuncs.end())
	{
		pot=new Potential;
		STRDBLMAP counts;
		estimateConditionalPotential(sFactor,varSet,&pot,counts);
		potFuncs[vId]=pot;
	}
	else
	{
		pot=potFuncs[vId];
	}
	//int sample=pot->generateSample(jointConf,vId,r);
	int sample=-1;
	return sample;
}

int
PotentialManager::clearJointEntropies()
{
	jointEntropies.clear();
	return 0;
}
