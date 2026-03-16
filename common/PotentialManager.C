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
			data->setValue(val,vId,eIter->first);
		}
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

//Need to fill in the covariances which are empty
int
PotentialManager::estimateNewCov_EM(map<int,EvidenceManager*> &evMgrSet, Potential* apot, int cid, map<int,map<int,INTDBLMAP*>*>& gammasubset)
{
	int dId=0;
	int covPair=0;
	VSET& varSet=apot->getAssocVariables();
	INTDBLMAP evidCnts;
	INTDBLMAP meanVect;
	map<int,INTDBLMAP*> covMat;
	for(map<int,EvidenceManager*>::iterator evIter=evMgrSet.begin();evIter!=evMgrSet.end();evIter++)
	{
		EvidenceManager* localMgr=evIter->second;
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		//First get the mean and then the variance
		map<int,INTDBLMAP*>* gammaCond=gammasubset[evIter->first];
		for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
		{
			EMAP* evidMap=localMgr->getEvidenceAt(eIter->first);
			INTDBLMAP* gMap=(*gammaCond)[eIter->first];
			double gval=(*gMap)[cid];
			for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end(); vIter++)
			{
				int vId=vIter->first;
				if(evidCnts.find(vId)==evidCnts.end())
				{
					evidCnts[vId]=gval;
				}
				else
				{
					evidCnts[vId]=evidCnts[vId]+gval;
				}
				Evidence* evid=(*evidMap)[vIter->first];
				double val=evid->getEvidVal();
				val=val*gval;
				if(meanVect.find(vId)==meanVect.end())
				{
					meanVect[vId]=val;
				}
				else
				{
					meanVect[vId]=meanVect[vId]+val;
				}
			}
		}
	}
	//Now estimate the mean
	for(INTDBLMAP_ITER idIter=meanVect.begin();idIter!=meanVect.end();idIter++)
	{
		double evidCnt=evidCnts[idIter->first];
		idIter->second=idIter->second/evidCnt;
	}
	//Now the variance
	//Now estimate the variance
	for(map<int,EvidenceManager*>::iterator evIter=evMgrSet.begin();evIter!=evMgrSet.end();evIter++)
	{
		EvidenceManager* localMgr=evIter->second;
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		map<int,INTDBLMAP*>* gammaCond=gammasubset[evIter->first];
		for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[eIter->first];
			double gval=(*gMap)[cid];
			EMAP* evidMap=localMgr->getEvidenceAt(eIter->first);
			for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
			{
				int vId=vIter->first;
				Evidence* evid=(*evidMap)[vIter->first];
				double vval=evid->getEvidVal();
				double vmean=meanVect[vId];
				INTDBLMAP* vcov=NULL;
				if(covMat.find(vId)==covMat.end())
				{
					vcov=new INTDBLMAP;
					covMat[vId]=vcov;
				}
				else
				{
					vcov=covMat[vId];
				}
				for(VSET_ITER uIter=vIter;uIter!=varSet.end();uIter++)
				{
					int uId=uIter->first;
					Evidence* evid1=(*evidMap)[uIter->first];
					double uval=evid1->getEvidVal();
					double umean=meanVect[uId];
					double diffprod=gval*(vval-vmean)*(uval-umean);
					INTDBLMAP* ucov=NULL;
					if(covMat.find(uId)==covMat.end())
					{
						ucov=new INTDBLMAP;
						covMat[uId]=ucov;
					}
					else
					{
						ucov=covMat[uId];
					}
					if(vcov->find(uId)==vcov->end())
					{
						covPair++;
						(*vcov)[uId]=diffprod;
					}
					else
					{
						(*vcov)[uId]=(*vcov)[uId]+diffprod;
					}
					if(uId!=vId)
					{
						if(ucov->find(vId)==ucov->end())
						{
							(*ucov)[vId]=diffprod;
						}
						else
						{
							(*ucov)[vId]=(*ucov)[vId]+diffprod;
						}
					}
				}
			}

		}
	}
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		INTDBLMAP* var=covMat[vIter->first];
		for(VSET_ITER uIter=varSet.begin();uIter!=varSet.end();uIter++)
		{
			if(uIter==vIter)
			{
				double cov=(0.001+(*var)[uIter->first])/evidCnts[vIter->first];
				(*var)[uIter->first]=cov;
			}
			else
			{
				double cov=(*var)[uIter->first]/evidCnts[vIter->first];
				(*var)[uIter->first]=cov;
			}
		}
	}
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end(); vIter++)
	{
		double mean=meanVect[vIter->first];
		INTDBLMAP* covar=covMat[vIter->first];
		apot->updateMean(vIter->first,mean);
		for(VSET_ITER uIter=vIter;uIter!=varSet.end();uIter++)
		{
			if(covar->find(uIter->first)==covar->end())
			{
				cerr <<"No var " << uIter->first << " in covariance of " << vIter->first << endl;
				exit(-1);
			}
			double cval=(*covar)[uIter->first];
			apot->updateCovariance(vIter->first,uIter->first,cval);
			apot->updateCovariance(uIter->first,vIter->first,cval);
		}
	}
	apot->makeValidJPD(ludecomp,perm);
	meanVect.clear();
	for(map<int,INTDBLMAP*>::iterator vIter=covMat.begin();vIter!=covMat.end();vIter++)
	{
		vIter->second->clear();
		delete vIter->second;
	}
	covMat.clear();
	return 0;
}

Error::ErrorCode
PotentialManager::populatePotentialsSlimFactors(map<int,SlimFactor*>& factorSet,VSET& varSet)
{
	//The set of flags to keep status of the potentials that have been calculated
	map<int,bool> doneFlag;
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		doneFlag[fIter->first]=false;
	}
	int popFId=0;
	for(map<int,SlimFactor*>::reverse_iterator rIter=factorSet.rbegin();rIter!=factorSet.rend();rIter++)
	{
		//If we have computed the potential for this flag move one
		if(doneFlag[rIter->first])
		{
			popFId++;
			continue;
		}
		SlimFactor* sFactor=rIter->second;
		if(sFactor->fId==176)
		{
			cout <<"Stop here " << endl;
		}
		//Otherwise create the potential
		Potential* aPotFunc=new Potential;
		for(int j=0;j<sFactor->vCnt;j++)
		{
			Variable* aVar=varSet[sFactor->vIds[j]];
			if(j==sFactor->vCnt-1)
			{
				aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			}
			else
			{
				aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			}
		}
		aPotFunc->potZeroInit();
		populatePotential(aPotFunc);
		aPotFunc->calculateJointEntropy();
		sFactor->jointEntropy=aPotFunc->getJointEntropy();
		if(sFactor->jointEntropy<0)
		{
		//	sFactor->jointEntropy=0;
		//	cout <<"Negative entropy for " << sFactor->fId << endl;
		}
		doneFlag[rIter->first]=true;
		delete aPotFunc;
		if(popFId%100000==0)
		{
			cout <<"Done with " << factorSet.size()-popFId << " factors " << endl;
		}
		popFId++;
	}
	return Error::SUCCESS;
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

int
PotentialManager::estimatePotWeight(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTINTMAP& vIDMatIDMap=apot->getVarMatrixIndexMap();
	map<int,INTDBLMAP*> sumys;
	INTDBLMAP sumxy;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[vId];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			double gammaval=(*gMap)[cid];
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double xpart=(*evidSet)[vId]->getEvidVal()-b;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				if(sumxy.find(vIter->first)==sumxy.end())
				{	
					sumxy[vIter->first]=gammaval*xpart*vval;
				}
				else
				{
					sumxy[vIter->first]=sumxy[vIter->first]+(gammaval*xpart*vval);
				}
				VSET_ITER uIter=vIter;
				for(;uIter!=potVarSet.end();uIter++)
				{
					if(uIter->first==vId)
					{
						continue;
					}
					double uval=(*evidSet)[uIter->first]->getEvidVal();
					double uvval=vval*uval*gammaval;
					uvval=uvval+lambda;
					INTDBLMAP* vrow=NULL;
					INTDBLMAP* urow=NULL;
					if(sumys.find(vIter->first)==sumys.end())
					{
						vrow=new INTDBLMAP;
						sumys[vIter->first]=vrow;
					}
					else
					{
						vrow=sumys[vIter->first];
					}	
					if(vrow->find(uIter->first)==vrow->end())
					{
						(*vrow)[uIter->first]=uvval;
					}
					else
					{
						(*vrow)[uIter->first]=(*vrow)[uIter->first]+uvval;
					}
					if(uIter==vIter)
					{
						continue;
					}
					if(sumys.find(uIter->first)==sumys.end())
					{
						urow=new INTDBLMAP;
						sumys[uIter->first]=urow;
					}
					else
					{
						urow=sumys[uIter->first];
					}
					if(urow->find(vIter->first)==urow->end())
					{
						(*urow)[vIter->first]=uvval;
					}
					else
					{
						(*urow)[vIter->first]=(*urow)[vIter->first]+uvval;
					}
				}
			}
		}
	}
	Matrix* covMatrix=new Matrix(potVarSet.size()-1,potVarSet.size()-1);
	Matrix* xyMatrix=new Matrix(1,potVarSet.size()-1);
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		INTDBLMAP* vrow=vIter->second;
		int rowvid=vIDMatIDMap[vId];
		int rowpos=vIDMatIDMap[vIter->first];
		if(rowpos>rowvid)
		{
			rowpos--;
		}
		double xyval=sumxy[vIter->first];
		xyMatrix->setValue(xyval,0,rowpos);
		for(INTDBLMAP_ITER rIter=vrow->begin();rIter!=vrow->end();rIter++)
		{
			int colvid=vIDMatIDMap[vId];
			int colpos=vIDMatIDMap[rIter->first];
			if(colpos>colvid)
			{
				colpos--;
			}
			covMatrix->setValue(rIter->second,rowpos,colpos);
			covMatrix->setValue(rIter->second,colpos,rowpos);
		}
	}
	Matrix* inv=covMatrix->invMatrix(ludecomp,perm);
	Matrix* wtMatrix=xyMatrix->multiplyMatrix(inv);
	apot->setCondWeight(wtMatrix);
	delete covMatrix;
	delete xyMatrix;
	delete inv;
	delete wtMatrix;
	sumxy.clear();
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		vIter->second->clear();
		delete vIter->second;
	}
	sumys.clear();
	return 0;

}


int
PotentialManager::estimatePotWeight(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTINTMAP& vIDMatIDMap=apot->getVarMatrixIndexMap();
	map<int,INTDBLMAP*> sumys;
	INTDBLMAP sumxy;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double xpart=(*evidSet)[vId]->getEvidVal()-b;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				if(sumxy.find(vIter->first)==sumxy.end())
				{	
					sumxy[vIter->first]=gammaval*xpart*vval;
				}
				else
				{
					sumxy[vIter->first]=sumxy[vIter->first]+(gammaval*xpart*vval);
				}
				VSET_ITER uIter=vIter;
				for(;uIter!=potVarSet.end();uIter++)
				{
					if(uIter->first==vId)
					{
						continue;
					}
					double uval=(*evidSet)[uIter->first]->getEvidVal();
					double uvval=vval*uval*gammaval;
					uvval=uvval+lambda;
					INTDBLMAP* vrow=NULL;
					INTDBLMAP* urow=NULL;
					if(sumys.find(vIter->first)==sumys.end())
					{
						vrow=new INTDBLMAP;
						sumys[vIter->first]=vrow;
					}
					else
					{
						vrow=sumys[vIter->first];
					}	
					if(vrow->find(uIter->first)==vrow->end())
					{
						(*vrow)[uIter->first]=uvval;
					}
					else
					{
						(*vrow)[uIter->first]=(*vrow)[uIter->first]+uvval;
					}
					if(uIter==vIter)
					{
						continue;
					}
					if(sumys.find(uIter->first)==sumys.end())
					{
						urow=new INTDBLMAP;
						sumys[uIter->first]=urow;
					}
					else
					{
						urow=sumys[uIter->first];
					}
					if(urow->find(vIter->first)==urow->end())
					{
						(*urow)[vIter->first]=uvval;
					}
					else
					{
						(*urow)[vIter->first]=(*urow)[vIter->first]+uvval;
					}
				}
			}
		}
	}
	Matrix* covMatrix=new Matrix(potVarSet.size()-1,potVarSet.size()-1);
	Matrix* xyMatrix=new Matrix(1,potVarSet.size()-1);
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		INTDBLMAP* vrow=vIter->second;
		int rowvid=vIDMatIDMap[vId];
		int rowpos=vIDMatIDMap[vIter->first];
		if(rowpos>rowvid)
		{
			rowpos--;
		}
		double xyval=sumxy[vIter->first];
		xyMatrix->setValue(xyval,0,rowpos);
		for(INTDBLMAP_ITER rIter=vrow->begin();rIter!=vrow->end();rIter++)
		{
			int colvid=vIDMatIDMap[vId];
			int colpos=vIDMatIDMap[rIter->first];
			if(colpos>colvid)
			{
				colpos--;
			}
			covMatrix->setValue(rIter->second,rowpos,colpos);
			covMatrix->setValue(rIter->second,colpos,rowpos);
		}
	}
	Matrix* inv=covMatrix->invMatrix(ludecomp,perm);
	Matrix* wtMatrix=xyMatrix->multiplyMatrix(inv);
	apot->setCondWeight(wtMatrix);
	delete covMatrix;
	delete xyMatrix;
	delete inv;
	delete wtMatrix;
	sumxy.clear();
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		vIter->second->clear();
		delete vIter->second;
	}
	sumys.clear();
	return 0;
}


int
PotentialManager::estimatePotWeight_TiedParam(map<int,Potential*>& potSet, int uId, int vId, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas,int componentID)
{
	INTDBLMAP numpart;
	INTDBLMAP denpart;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
			{
				Potential* apot=pIter->second;
				double b=apot->getCondBias();
				VSET& potVarSet=apot->getAssocVariables();
				double xpart=(*evidSet)[uId]->getEvidVal()-b;
				double gammaval=(*gMap)[pIter->first];
				INTDBLMAP& oldWeight=apot->getCondWeight();
				double ypart=0;
				double vval=0;
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==uId)
					{
						continue;
					}
					if(vIter->first==vId)
					{
						vval=(*evidSet)[vIter->first]->getEvidVal();
						continue;
					}
					double aval=(*evidSet)[vIter->first]->getEvidVal();
					double wt=oldWeight[vIter->first];
					ypart=ypart+(aval*wt);
				}
				xpart=gammaval*(xpart-ypart);
				if(numpart.find(pIter->first)==numpart.end())
				{
					numpart[pIter->first]=xpart*vval;
				}
				else
				{
					numpart[pIter->first]=numpart[pIter->first]+(xpart*vval);
				}
				double denval=vval*vval*gammaval;
				if(denpart.find(pIter->first)==denpart.end())
				{
					denpart[pIter->first]=denval;
				}
				else
				{
					denpart[pIter->first]=denpart[pIter->first]+denval;
				}
			}
		}
	}
	double numerator=0;
	double denominator=0;
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* pot=pIter->second;
		double condvar=pot->getCondVariance();
		double npart=numpart[pIter->first];
		npart=npart/condvar;
		numerator=numerator+npart;
		double dpart=denpart[pIter->first];
		dpart=dpart/condvar;
		denominator=denominator+dpart;
	}
	double wt=numerator/denominator;

	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* pot=pIter->second;
		pot->setCondWeightFor(vId,wt);
	}
	numpart.clear();
	denpart.clear();
	return 0;
}


int
PotentialManager::estimatePotWeight_TiedParam(map<int,Potential*>& potSet, int uId, int vId, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	INTDBLMAP numpart;
	INTDBLMAP denpart;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
			{
				Potential* apot=pIter->second;
				double b=apot->getCondBias();
				VSET& potVarSet=apot->getAssocVariables();
				double xpart=(*evidSet)[uId]->getEvidVal()-b;
				double gammaval=(*gMap)[pIter->first];
				INTDBLMAP& oldWeight=apot->getCondWeight();
				double ypart=0;
				double vval=0;
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==uId)
					{
						continue;
					}
					if(vIter->first==vId)
					{
						vval=(*evidSet)[vIter->first]->getEvidVal();
						continue;
					}
					double aval=(*evidSet)[vIter->first]->getEvidVal();
					double wt=oldWeight[vIter->first];
					ypart=ypart+(aval*wt);
				}
				xpart=gammaval*(xpart-ypart);
				if(numpart.find(pIter->first)==numpart.end())
				{
					numpart[pIter->first]=xpart*vval;
				}
				else
				{
					numpart[pIter->first]=numpart[pIter->first]+(xpart*vval);
				}
				double denval=vval*vval*gammaval;
				if(denpart.find(pIter->first)==denpart.end())
				{
					denpart[pIter->first]=denval;
				}
				else
				{
					denpart[pIter->first]=denpart[pIter->first]+denval;
				}
			}
		}
	}
	double numerator=0;
	double denominator=0;
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* pot=pIter->second;
		double condvar=pot->getCondVariance();
		double npart=numpart[pIter->first];
		npart=npart/condvar;
		numerator=numerator+npart;
		double dpart=denpart[pIter->first];
		dpart=dpart/condvar;
		denominator=denominator+dpart;
	}
	double wt=numerator/denominator;

	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* pot=pIter->second;
		pot->setCondWeightFor(vId,wt);
	}
	numpart.clear();
	denpart.clear();
	return 0;
}


int
PotentialManager::estimatePotWeight_TiedParamSet(map<int,Potential*>& potSet, int xId, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	INTDBLMAP sumxy;
	map<int,INTDBLMAP*> sumys;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			INTDBLMAP exclParts;
			for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
			{
				Potential* apot=pIter->second;
				double b=apot->getCondBias();
				VSET& potVarSet=apot->getAssocVariables();
				INTINTMAP sharedMBVars=apot->getSharedMBVars();
				INTDBLMAP& oldWeight=apot->getCondWeight();
				double exclpart=0;
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==xId)
					{
						continue;
					}
					if(sharedMBVars.find(vIter->first)!=sharedMBVars.end())
					{
						continue;
					}
					if(oldWeight.find(vIter->first)==oldWeight.end())
					{
						cerr <<"No weight for  " << vIter->first << " in cond of " << xId << endl;
						exit(-1);
					}
					Evidence* evid=(*evidSet)[vIter->first];
					double aval=oldWeight[vIter->first]*evid->getEvidVal();
					exclpart=exclpart+aval;
				}
				double xpart=(*evidSet)[xId]->getEvidVal()-b-exclpart;
				double c=apot->getCondVariance();
				double gval=(*gMap)[pIter->first];
				xpart=(gval*xpart)/c;
				exclParts[pIter->first]=xpart;
			}
			Potential* apot=potSet.begin()->second;
			INTINTMAP& sharedMBVars=apot->getSharedMBVars();
			for(INTINTMAP_ITER sIter=sharedMBVars.begin();sIter!=sharedMBVars.end();sIter++)
			{
				double aval=(*evidSet)[sIter->first]->getEvidVal();
				double asum=0;
				for(INTDBLMAP_ITER pIter=exclParts.begin();pIter!=exclParts.end();pIter++)
				{
					asum=asum+(aval*pIter->second);
				}
				if(sumxy.find(sIter->first)==sumxy.end())
				{
					sumxy[sIter->first]=asum;
				}
				else
				{
					sumxy[sIter->first]=sumxy[sIter->first]+asum;
				}
				INTINTMAP_ITER qIter=sIter;
				for(;qIter!=sharedMBVars.end();qIter++)
				{
					double bval=(*evidSet)[qIter->first]->getEvidVal();
					double uvval=0;
					for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
					{
						Potential* apot=pIter->second;
						double gammaval=(*gMap)[pIter->first];
						double c=apot->getCondVariance();
						double dval=(gammaval*aval*bval)/c;
						uvval=uvval+dval;
					}
					INTDBLMAP* vrow=NULL;
					INTDBLMAP* urow=NULL;
					if(sumys.find(sIter->first)==sumys.end())
					{
						vrow=new INTDBLMAP;
						sumys[sIter->first]=vrow;
					}
					else
					{
						vrow=sumys[sIter->first];
					}	
					if(vrow->find(qIter->first)==vrow->end())
					{
						(*vrow)[qIter->first]=uvval;
					}
					else
					{
						(*vrow)[qIter->first]=(*vrow)[qIter->first]+uvval;
					}
					if(sIter==qIter)
					{
						continue;
					}
					if(sumys.find(qIter->first)==sumys.end())
					{
						urow=new INTDBLMAP;
						sumys[qIter->first]=urow;
					}
					else
					{
						urow=sumys[qIter->first];
					}
					if(urow->find(qIter->first)==urow->end())
					{
						(*urow)[sIter->first]=uvval;
					}
					else
					{
						(*urow)[sIter->first]=(*urow)[sIter->first]+uvval;
					}

				}
			}
			exclParts.clear();
		}
	}
	Potential* apot=potSet.begin()->second;
	int sharedNeighbors=apot->getSharedMBVars().size();
	Matrix* covMatrix=new Matrix(sharedNeighbors,sharedNeighbors);
	Matrix* xyMatrix=new Matrix(1,sharedNeighbors);
	map<int,int> vIDMatIDMap;
	int rowvid=0;
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		INTDBLMAP* vrow=vIter->second;
		int rowpos=-1;
		if(vIDMatIDMap.find(vIter->first)==vIDMatIDMap.end())
		{
			vIDMatIDMap[vIter->first]=rowvid;
			rowvid++;
		}
		rowpos=vIDMatIDMap[vIter->first];
		double xyval=sumxy[vIter->first];
		xyMatrix->setValue(xyval,0,rowpos);
		for(INTDBLMAP_ITER rIter=vrow->begin();rIter!=vrow->end();rIter++)
		{
			int colpos=-1;
			if(vIDMatIDMap.find(rIter->first)==vIDMatIDMap.end())
			{
				vIDMatIDMap[rIter->first]=rowvid;
				rowvid++;
			}
			colpos=vIDMatIDMap[rIter->first];
			covMatrix->setValue(rIter->second,rowpos,colpos);
			covMatrix->setValue(rIter->second,colpos,rowpos);
		}
	}
	Matrix* inv=covMatrix->invMatrix(ludecomp,perm);
	Matrix* wtMatrix=xyMatrix->multiplyMatrix(inv);
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* lpot=pIter->second;
		INTDBLMAP& potWt=lpot->getCondWeight();
		for(INTINTMAP_ITER vIter=vIDMatIDMap.begin();vIter!=vIDMatIDMap.end();vIter++)
		{
			int pos=vIter->second;
			double wval=wtMatrix->getValue(0,pos);
			lpot->setCondWeightFor(vIter->first,wval);
		}
	}
	vIDMatIDMap.clear();
	delete covMatrix;
	delete xyMatrix;
	delete inv;
	delete wtMatrix;
		
	sumxy.clear();
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		vIter->second->clear();
		delete vIter->second;
	}
	sumys.clear();
	return 0;
}


int
PotentialManager::estimatePotWeight_ExclParam(Potential* pot, int uId, int vId, int cid, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas, int componentID)
{
	double numerator=0;
	double denominator=0;
	double b=pot->getCondBias();
	VSET& potVarSet=pot->getAssocVariables();
	INTDBLMAP& oldWeight=pot->getCondWeight();
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xpart=(*evidSet)[uId]->getEvidVal()-b;
			double ypart=0;
			double vval=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==uId)
				{
					continue;
				}
				if(vIter->first==vId)
				{
					vval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double aval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=oldWeight[vIter->first];
				ypart=ypart+(aval*wt);
			}
			xpart=gammaval*(xpart-ypart);
			numerator=numerator+(xpart*vval);
			double denval=vval*vval*gammaval;
			denominator=denominator+denval;
		}
	}
	double wt=numerator/denominator;
	pot->setCondWeightFor(vId,wt);
	return 0;
}

int
PotentialManager::estimatePotWeight_ExclParam(Potential* pot, int uId, int vId, int cid, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	double numerator=0;
	double denominator=0;
	double b=pot->getCondBias();
	VSET& potVarSet=pot->getAssocVariables();
	INTDBLMAP& oldWeight=pot->getCondWeight();
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xpart=(*evidSet)[uId]->getEvidVal()-b;
			double ypart=0;
			double vval=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==uId)
				{
					continue;
				}
				if(vIter->first==vId)
				{
					vval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double aval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=oldWeight[vIter->first];
				ypart=ypart+(aval*wt);
			}
			xpart=gammaval*(xpart-ypart);
			numerator=numerator+(xpart*vval);
			double denval=vval*vval*gammaval;
			denominator=denominator+denval;
		}
	}
	double wt=numerator/denominator;
	pot->setCondWeightFor(vId,wt);
	return 0;
}


int 
PotentialManager::estimatePotWeight_ExclParamSet(map<int,Potential*>& potSet,int xId, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	map<int,map<int,INTDBLMAP*>*> sumys_set;
	map<int,INTDBLMAP*> sumxy_set;
	Potential* apot=potSet.begin()->second;
	INTINTMAP& sharedMBVars=potSet.begin()->second->getSharedMBVars();
	INTDBLMAP& wt=apot->getCondWeight();
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaSet=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			//First get the shared part
			//The parameters of the shared part should be the same
			double sharedpart=0;
			for(INTINTMAP_ITER vIter=sharedMBVars.begin();vIter!=sharedMBVars.end();vIter++)
			{
				if(wt.find(vIter->first)==wt.end())
				{
					cerr <<"Error! No weight value for " << vIter->first << " in cond of " << xId << endl;
					exit(-1);
				}
				Evidence* evid=(*evidSet)[vIter->first];
				double aval=wt[vIter->first]*evid->getEvidVal();
				sharedpart=sharedpart+aval;
			}
			//Now the specific parts
			for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
			{
				Potential* apot=pIter->second;
				VSET& potVarSet=apot->getAssocVariables();
				double gammaval=(*gMap)[pIter->first];
				double b=apot->getCondBias();
				double xpart=(*evidSet)[xId]->getEvidVal()-b-sharedpart;
				INTDBLMAP* sumxy=NULL;
				map<int,INTDBLMAP*>* sumys=NULL;
				if(sumxy_set.find(pIter->first)==sumxy_set.end())
				{
					sumxy=new INTDBLMAP;
					sumxy_set[pIter->first]=sumxy;
					sumys=new map<int,INTDBLMAP*>;
					sumys_set[pIter->first]=sumys;
				}
				else
				{	
					sumxy=sumxy_set[pIter->first];
					sumys=sumys_set[pIter->first];
				}
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==xId)
					{
						continue;
					}
					if(sharedMBVars.find(vIter->first)!=sharedMBVars.end())
					{
						continue;
					}
					double vval=(*evidSet)[vIter->first]->getEvidVal();
					if(sumxy->find(vIter->first)==sumxy->end())
					{	
						(*sumxy)[vIter->first]=gammaval*xpart*vval;
					}
					else
					{
						(*sumxy)[vIter->first]=(*sumxy)[vIter->first]+(gammaval*xpart*vval);
					}
					VSET_ITER uIter=vIter;
					for(;uIter!=potVarSet.end();uIter++)
					{
						if(uIter->first==xId)
						{
							continue;
						}
						if(sharedMBVars.find(uIter->first)!=sharedMBVars.end())
						{
							continue;
						}
						double uval=(*evidSet)[uIter->first]->getEvidVal();
						double uvval=vval*uval*gammaval;
						INTDBLMAP* vrow=NULL;
						INTDBLMAP* urow=NULL;
						if(sumys->find(vIter->first)==sumys->end())
						{
							vrow=new INTDBLMAP;
							(*sumys)[vIter->first]=vrow;
						}
						else
						{
							vrow=(*sumys)[vIter->first];
						}	
						if(vrow->find(uIter->first)==vrow->end())
						{
							(*vrow)[uIter->first]=uvval;
						}
						else
						{
							(*vrow)[uIter->first]=(*vrow)[uIter->first]+uvval;
						}
						if(uIter==vIter)
						{
							continue;
						}
						if(sumys->find(uIter->first)==sumys->end())
						{
							urow=new INTDBLMAP;
							(*sumys)[uIter->first]=urow;
						}
						else
						{
							urow=(*sumys)[uIter->first];
						}
						if(urow->find(vIter->first)==urow->end())
						{
							(*urow)[vIter->first]=uvval;
						}
						else
						{
							(*urow)[vIter->first]=(*urow)[vIter->first]+uvval;
						}
					}
				}
			}
		}
	}
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* apot=pIter->second;
		VSET& potVarSet=apot->getAssocVariables();
		int unsharedNeighbors=potVarSet.size()-sharedMBVars.size()-1;
		if(unsharedNeighbors==0)
		{
			continue;
		}
		Matrix* covMatrix=new Matrix(unsharedNeighbors,unsharedNeighbors);
		Matrix* xyMatrix=new Matrix(1,unsharedNeighbors);
		map<int,INTDBLMAP*>* sumys=sumys_set[pIter->first];
		INTDBLMAP* sumxy=sumxy_set[pIter->first];
		map<int,int> vIDMatIDMap;
		int rowvid=0;
		for(map<int,INTDBLMAP*>::iterator vIter=sumys->begin();vIter!=sumys->end();vIter++)
		{
			INTDBLMAP* vrow=vIter->second;
			int rowpos=-1;
			if(vIDMatIDMap.find(vIter->first)==vIDMatIDMap.end())
			{
				vIDMatIDMap[vIter->first]=rowvid;
				rowvid++;
			}
			rowpos=vIDMatIDMap[vIter->first];
			double xyval=(*sumxy)[vIter->first];
			xyMatrix->setValue(xyval,0,rowpos);
			for(INTDBLMAP_ITER rIter=vrow->begin();rIter!=vrow->end();rIter++)
			{
				int colpos=-1;
				if(vIDMatIDMap.find(rIter->first)==vIDMatIDMap.end())
				{
					vIDMatIDMap[rIter->first]=rowvid;
					rowvid++;
				}
				colpos=vIDMatIDMap[rIter->first];
				covMatrix->setValue(rIter->second,rowpos,colpos);
				covMatrix->setValue(rIter->second,colpos,rowpos);
			}
		}
		Matrix* inv=covMatrix->invMatrix(ludecomp,perm);
		Matrix* wtMatrix=xyMatrix->multiplyMatrix(inv);
		INTDBLMAP& potWt=apot->getCondWeight();
		for(INTINTMAP_ITER vIter=vIDMatIDMap.begin();vIter!=vIDMatIDMap.end();vIter++)
		{
			int pos=vIter->second;
			double wval=wtMatrix->getValue(0,pos);
			apot->setCondWeightFor(vIter->first,wval);
		}
		vIDMatIDMap.clear();
		delete covMatrix;
		delete xyMatrix;
		delete inv;
		delete wtMatrix;
		
		sumxy->clear();
		delete sumxy;
		for(map<int,INTDBLMAP*>::iterator vIter=sumys->begin();vIter!=sumys->end();vIter++)
		{
			vIter->second->clear();
			delete vIter->second;
		}
		sumys->clear();
		delete sumys;
	}
	sumxy_set.clear();
	sumys_set.clear();
	return 0;
}
/*

int
PotentialManager::estimatePotWeight(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTINTMAP& vIDMatIDMap=apot->getVarMatrixIndexMap();
	map<int,INTDBLMAP*> sumys;
	INTDBLMAP sumxy;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double xpart=(*evidSet)[vId]->getEvidVal()-b;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				if(sumxy.find(vIter->first)==sumxy.end())
				{	
					sumxy[vIter->first]=gammaval*xpart*vval;
				}
				else
				{
					sumxy[vIter->first]=sumxy[vIter->first]+(gammaval*xpart*vval);
				}
				VSET_ITER uIter=vIter;
				for(;uIter!=potVarSet.end();uIter++)
				{
					if(uIter->first==vId)
					{
						continue;
					}
					double uval=(*evidSet)[uIter->first]->getEvidVal();
					double uvval=vval*uval*gammaval;
					uvval=uvval+lambda;
					INTDBLMAP* vrow=NULL;
					INTDBLMAP* urow=NULL;
					if(sumys.find(vIter->first)==sumys.end())
					{
						vrow=new INTDBLMAP;
						sumys[vIter->first]=vrow;
					}
					else
					{
						vrow=sumys[vIter->first];
					}	
					if(vrow->find(uIter->first)==vrow->end())
					{
						(*vrow)[uIter->first]=uvval;
					}
					else
					{
						(*vrow)[uIter->first]=(*vrow)[uIter->first]+uvval;
					}
					if(uIter==vIter)
					{
						continue;
					}
					if(sumys.find(uIter->first)==sumys.end())
					{
						urow=new INTDBLMAP;
						sumys[uIter->first]=urow;
					}
					else
					{
						urow=sumys[uIter->first];
					}
					if(urow->find(vIter->first)==urow->end())
					{
						(*urow)[vIter->first]=uvval;
					}
					else
					{
						(*urow)[vIter->first]=(*urow)[vIter->first]+uvval;
					}
				}
			}
		}
	}
	Matrix* covMatrix=new Matrix(potVarSet.size()-1,potVarSet.size()-1);
	Matrix* xyMatrix=new Matrix(1,potVarSet.size()-1);
	for(map<int,INTDBLMAP*>::iterator vIter=sumys.begin();vIter!=sumys.end();vIter++)
	{
		INTDBLMAP* vrow=vIter->second;
		int rowvid=vIDMatIDMap[vId];
		int rowpos=vIDMatIDMap[vIter->first];
		if(rowpos>rowvid)
		{
			rowpos--;
		}
	}
	return 0;
}*/


int
PotentialManager::estimatePotCondVarBias(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas, int componentID)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTDBLMAP& condWt=apot->getCondWeight();
	double varnum=0;
	double varden=0;
	double biasnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xval=0;
			varden=varden+gammaval;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)) + (lambda*shrinkcorr));
			biasnum=biasnum+(gammaval*(xval-yval));
		}
	}

	double condBias=biasnum/varden;
	b=condBias;
	varnum=0;

	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* gammaCond=gammas[eIter->first];
		map<int,INTDBLMAP*>* gammaSet=(*gammaCond)[componentID];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaSet)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xval=0;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)) + (lambda*shrinkcorr));
		}
	}
	double condCov=varnum/(varden*(1+lambda));
	if(lambda>0)
	{
		condCov=condCov/2;
	}
	apot->setCondVariance(condCov);
	apot->setUnnormCondVariance(varnum);
	apot->setCondBias(condBias);
	return 0;
}


int
PotentialManager::estimatePotCondVarBias(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,INTDBLMAP*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTDBLMAP& condWt=apot->getCondWeight();
	double varnum=0;
	double varden=0;
	double biasnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xval=0;
			varden=varden+gammaval;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)) + (lambda*shrinkcorr));
			biasnum=biasnum+(gammaval*(xval-yval));
		}
	}
	double condBias=biasnum/varden;
	b=condBias;
	varnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,INTDBLMAP*>* gammaCond=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			INTDBLMAP* gMap=(*gammaCond)[tIter->first];
			double gammaval=(*gMap)[cid];
			double xval=0;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)) + (lambda*shrinkcorr));
		}
	}
	double condCov=varnum/(varden*(1+lambda));
	if(lambda>0)
	{
		condCov=condCov/2;
	}
	apot->setCondVariance(condCov);
	apot->setUnnormCondVariance(varnum);
	apot->setCondBias(condBias);
	return 0;
}


int
PotentialManager::estimatePotTiedVar(map<int,Potential*>& potSet,int vId, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas)
{
	double varnum=0;
	double varden=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* dataGammaMap=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			map<int,INTDBLMAP*>* dataGamma=(*dataGammaMap)[tIter->first];
			INTDBLMAP* gMap=(*dataGamma)[vId];
			for(INTDBLMAP_ITER cIter=gMap->begin();cIter!=gMap->end();cIter++)
			{
				Potential* apot=potSet[cIter->first];
				VSET& potVarSet=apot->getAssocVariables();
				INTDBLMAP& condWt=apot->getCondWeight();
				double b=apot->getCondBias();
				double gammaval=(*gMap)[cIter->first];
				double xval=0;
				varden=varden+gammaval;
				EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
				double yval=0;
				double shrinkcorr=0;
				for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
				{
					if(vIter->first==vId)
					{
						xval=(*evidSet)[vIter->first]->getEvidVal();
						continue;
					}
					double vval=(*evidSet)[vIter->first]->getEvidVal();
					double wt=condWt[vIter->first];
					yval=yval+(vval*wt);
				}
				varnum=varnum+(gammaval*((xval-yval-b)*(xval-yval-b)));
			}
		}
	}
	double condCov=varnum/varden;
	for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
	{
		Potential* apot=pIter->second;
		apot->setCondVariance(condCov);
		apot->setUnnormCondVariance(varnum);
	}
	return 0;
}


int
PotentialManager::estimatePotCondBias(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTDBLMAP& condWt=apot->getCondWeight();
	double biasden=0;
	double biasnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* dataGammaMap=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			map<int,INTDBLMAP*>* dataGamma=(*dataGammaMap)[tIter->first];
			INTDBLMAP* gMap=(*dataGamma)[vId];
			double gammaval=(*gMap)[cid];
			double xval=0;
			biasden=biasden+gammaval;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
			}
			biasnum=biasnum+(gammaval*(xval-yval));
		}
	}
	double condBias=biasnum/biasden;
	apot->setCondBias(condBias);
	return 0;
}

int
PotentialManager::estimatePotCondBiasThenVar(Potential* apot,int vId, int cid,double b, map<int,EvidenceManager*>& evSet,map<int,map<int,map<int,INTDBLMAP*>*>*>& gammas)
{
	VSET& potVarSet=apot->getAssocVariables();
	INTDBLMAP& condWt=apot->getCondWeight();
	double den=0;
	double biasnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* dataGammaMap=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			map<int,INTDBLMAP*>* dataGamma=(*dataGammaMap)[tIter->first];
			INTDBLMAP* gMap=(*dataGamma)[vId];
			double gammaval=(*gMap)[cid];
			den=den+gammaval;
			double xval=0;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
			}
			biasnum=biasnum+(gammaval*(xval-yval));
		}
	}


	double condBias=biasnum/den;
	double varnum=0;
	for(map<int,EvidenceManager*>::iterator eIter=evSet.begin();eIter!=evSet.end();eIter++)
	{
		EvidenceManager* localMgr=eIter->second;
		map<int,map<int,INTDBLMAP*>*>* dataGammaMap=gammas[eIter->first];
		INTINTMAP& trainEvidSet=localMgr->getTestSet();
		for(INTINTMAP_ITER tIter=trainEvidSet.begin();tIter!=trainEvidSet.end();tIter++)
		{
			map<int,INTDBLMAP*>* dataGamma=(*dataGammaMap)[tIter->first];
			INTDBLMAP* gMap=(*dataGamma)[vId];
			double gammaval=(*gMap)[cid];
			double xval=0;
			EMAP* evidSet=localMgr->getEvidenceAt(tIter->first);
			double yval=0;
			double shrinkcorr=0;
			for(VSET_ITER vIter=potVarSet.begin();vIter!=potVarSet.end();vIter++)
			{
				if(vIter->first==vId)
				{
					xval=(*evidSet)[vIter->first]->getEvidVal();
					continue;
				}
				double vval=(*evidSet)[vIter->first]->getEvidVal();
				double wt=condWt[vIter->first];
				yval=yval+(vval*wt);
				shrinkcorr=shrinkcorr+(wt*wt);
			}
			varnum=varnum+(gammaval*((xval-yval-condBias)*(xval-yval-condBias)) + (lambda*shrinkcorr));
		}
	}
	double condCov=varnum/(den*(1+lambda));
	if(lambda>0)
	{
		condCov=condCov/2;
	}
	apot->setCondVariance(condCov);
	apot->setUnnormCondVariance(varnum);
	apot->setCondBias(condBias);
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

double
PotentialManager::estimateCanonicalValue(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* condPot)
{
	double pVal=0;
	cout <<"Not implemented " << endl;
	return 0;
}

double
PotentialManager::estimateCanonicalValue_Joint(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* jointPot)
{
	cout <<"Not implemented " << endl;
	return 0;
}
