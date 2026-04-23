#include <iostream>
#include <cstring>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>

#include "Error.H"
#include "Variable.H"
#include "SlimFactor.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "VariableManager.H"

#include "LatticeStructure.H"

#include "PotentialManager.H"
#include "FactorGraph.H"

#include "Vertex.H"
#include "Graph.H"
#include "FactorManager.H"

FactorManager::FactorManager()
{
	globalFactorID=0;
	maxFactorSize_Approx=-1;
	mbPenalty=0;
}

FactorManager::~FactorManager()
{
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		delete fIter->second;
	}
	slimFactorSet.clear();
	lattice.clear();
	factorNameToIDMap.clear();
	factorIDToNameMap.clear();
	mbSpecific_MI.clear();
}

int 
FactorManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

int 
FactorManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int
FactorManager::setPotentialManager(PotentialManager* aPtr)
{
	potMgr=aPtr;
	return 0;
}

int 
FactorManager::setBaseInstantiation()
{
	VSET& variableSet=vMgr->getVariableSet();
	/*for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		Variable* aVar=vIter->second;
		INTVECT& varVals=aVar->getValues();
		defaultInstantiation[aVar->getID()]=varVals[0];
	}*/
	int evidCnt=evMgr->getNumberOfEvidences();
	map<string,int> evidDist;
	map<string,int> evidStrIDMap;
	char confStr[CONSTR_LEN];
	for(int e=0;e<evidCnt;e++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(e);
		string aConfStr;
		
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sprintf(confStr,"-%d=%d-",eIter->first,eIter->second->getHardEvidVal());
			aConfStr.append(confStr);
		}
		if(evidDist.find(aConfStr)==evidDist.end())
		{
			evidDist[aConfStr]=1;
		}
		else
		{
			evidDist[aConfStr]=evidDist[aConfStr]+1;
		}
		evidStrIDMap[aConfStr]=e;
	}
	int maxCnt=0;
	int minCnt=evidCnt;
	int maxevid=-1;
	int minevid=-1;
	for(map<string,int>::iterator aIter=evidDist.begin();aIter!=evidDist.end();aIter++)
	{
		if(aIter->second>maxCnt)
		{
			maxCnt=aIter->second;
			maxevid=evidStrIDMap[aIter->first];
		}
		if(aIter->second<minCnt)
		{
			minCnt=aIter->second;
			minevid=evidStrIDMap[aIter->first];
		}
	}
	//maxevid=0;
	EMAP* evidMap=evMgr->getEvidenceAt(maxevid);
//	EMAP* evidMap=evMgr->getEvidenceAt(minevid);
	string defConfStr;
	for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
	{
		defaultInstantiation[eIter->first]=eIter->second->getHardEvidVal();
		sprintf(confStr,"-%d=%d-",eIter->first,eIter->second->getHardEvidVal());
		defConfStr.append(confStr);
		defInstMap[eIter->first]=eIter->second->getHardEvidVal();
	}
	defProb=(double)evidDist[defConfStr]/(double)evidCnt;
	defInstID=maxevid;
	cout << "Default Inst: ID " << maxevid << " Freq: " << evidDist[defConfStr] << " Def prob: " << defProb << endl;
	cout << "Default Inst: " << defConfStr.c_str() << endl;
	return 0;
}

int 
FactorManager::setBaseInstantiation_Variable()
{
	VSET& variableSet=vMgr->getVariableSet();
	int evidCnt=evMgr->getNumberOfEvidences();
	map<string,int> evidStrIDMap;
	char confStr[CONSTR_LEN];
	INTINTMAP defInstMap;
	string defInstStr;
	defProb=0;
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		map<int,double> vvFreq;
		double total=0;
		for(int e=0;e<evidCnt;e++)
		{
			EMAP* evidMap=evMgr->getEvidenceAt(e);
			EMAP_ITER eIter=evidMap->find(vIter->first);
			int vval=eIter->second->getHardEvidVal();
			if(vvFreq.find(vval)==vvFreq.end())
			{
				vvFreq[vval]=1;
			}
			else
			{
				vvFreq[vval]=vvFreq[vval]+1;
			}
			total=total+1;
		}
		double maxFreq=0;
		int maxVal=0;
		for(map<int,double>::iterator fIter=vvFreq.begin();fIter!=vvFreq.end();fIter++)
		{
			if(fIter->second>maxFreq)
			{
				maxFreq=fIter->second;
				maxVal=fIter->first;
			}
		}
		defInstMap[vIter->first]=maxVal;
		sprintf(confStr,"-%d=%d-",vIter->first,maxVal);
		defInstStr.append(confStr);
		double pval=maxFreq/total;
		//defProb=defProb+log(pval);
		if(pval>defProb)
		{
			defProb=pval;
		}
	}
	//defProb=exp(defProb);
	cout << "Default Inst prob: "<<  defProb << endl;
	cout << "Default Inst: " << defInstStr.c_str() << endl;
	return 0;
}

int 
FactorManager::setMaxFactorSize_Approx(int size)
{
	maxFactorSize_Approx=size;
	return 0;
}

int 
FactorManager::getMaxFactorSize_Approx()
{
	return maxFactorSize_Approx;
}

int
FactorManager::setOutputDir(const char* aPtr)
{
	strcpy(outputDir,aPtr);
	return 0;
}

int
FactorManager::setModelName(const char* aPtr)
{
	strcpy(modelName,aPtr);
	return 0;
}

int
FactorManager::setGraph(Graph* gPtr)
{
	graph=gPtr;
	return 0;
}

int
FactorManager::setBeamSize(int aSize)
{
	beamSize=aSize;
	return 0;
}

int
FactorManager::setPenalty(double p)
{
	mbPenalty=p;
	return 0;
}

int
FactorManager::allocateFactorSpace()
{
	initFactorSet();
	populateFactorSet();
	return 0;
}

int
FactorManager::allocateFactorSpace_Graph()
{
	generateSingletonFactors();
	generateCanonicalFactors();
	return 0;
}

FactorGraph* 
FactorManager::createInitialFactorGraph()
{
	FactorGraph* fg=new FactorGraph;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* factor=fIter->second;
		if(factor->vCnt>1)
		{
			break;
		}
		SlimFactor* newFactor=new SlimFactor;
		newFactor->vCnt=1;
		newFactor->vIds=new int[1];
		newFactor->vIds[0]=factor->vIds[0];
		newFactor->fId=factor->fId;
		newFactor->mutualInfo=factor->mutualInfo;
		newFactor->jointEntropy=factor->jointEntropy;
		newFactor->marginalEntropy=factor->jointEntropy;
		newFactor->mbScore=factor->mbScore;
		newFactor->moveScore=factor->moveScore;
		fg->setFactor(newFactor);
	}
	return fg;
}

int
FactorManager::generateSingletonFactors()
{
	VSET& variableSet=vMgr->getVariableSet();
	//Create the factors and the Markov blanket variables using the neighbours of each variable
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		Vertex* v=graph->getVertex(vIter->second->getName().c_str());
		if(v==NULL)
		{
			cout <<"No vertex corresponding to " << vIter->second->getName() << endl;
			return -1;
		}
		SlimFactor* sFactor=new SlimFactor;
		canonicalFactorSet[globalFactorID]=sFactor;
		sFactor->vIds=new int[1];
		sFactor->vIds[0]=vIter->first;
		sFactor->vCnt=1;
		//sFactor->secondPId=-1;
		sFactor->mutualInfo=0;
		sFactor->jointEntropy=0;
		sFactor->fId=globalFactorID;
		string key;
		getFactorKey(sFactor->vIds,sFactor->vCnt,key);
		factorNameToIDMap[key]=sFactor->fId;
		factorIDToNameMap[sFactor->fId]=key;
		globalFactorID++;

		NINFO_MAP& neighbours=v->getImmediateNeighbours();
		for(NINFO_MAP_ITER nmIter=neighbours.begin();nmIter!=neighbours.end();nmIter++)
		//for(VSET_ITER nmIter=variableSet.begin();nmIter!=variableSet.end();nmIter++)
		{
			int vId=vMgr->getVarID(nmIter->first.c_str());
			//int vId=nmIter->first;
			if(vId==-1)
			{
				cout <<"Did not find variable id in variableManager for " << nmIter->first.c_str() << endl;
				return -1;
			}
			//Take care of the self-loops
			if(vId==sFactor->vIds[0])
			{
				continue;
			}
			sFactor->mergedMB[vId]=0;
		}
	}
	cout <<"Populated " << canonicalFactorSet.size() << " factors and their MB from graph " << endl;
	return 0;
}

int
FactorManager::generateCanonicalFactors()
{
	VSET& variableSet=vMgr->getVariableSet();
	int** subsetSpace=new int* [variableSet.size()];
	for(int i=0;i<variableSet.size();i++)
	{
		subsetSpace[i]=new int[variableSet.size()-1];
	}

	INTINTMAP parentIDs;
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactorSet.begin();aIter!=canonicalFactorSet.end();aIter++)
	{
		parentIDs[aIter->first]=0;
	}
	int fSize=2;
	while(parentIDs.size()>0)
	{
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		INTINTMAP tempParentIDs;
		for(INTINTMAP_ITER pIter=parentIDs.begin();pIter!=parentIDs.end();pIter++)
		{
			SlimFactor* pFactor=canonicalFactorSet[pIter->first];
			for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
			{
				int newVId=vIter->first;
				if(newVId<pFactor->vIds[pFactor->vCnt-1])
				{
					continue;
				}
				//Now add newVId to if newVId is in the Markov blanket of all variables in pFactor
				bool clqChk=true;
				int vid=0;
				while((vid<pFactor->vCnt) && clqChk)
				{
					SlimFactor* qFactor=canonicalFactorSet[pFactor->vIds[vid]];
					if(qFactor->mergedMB.find(newVId)==qFactor->mergedMB.end())
					{
						clqChk=false;
					}
					vid++;
				}
				if(!clqChk)
				{
					continue;
				}
				SlimFactor* sFactor=new SlimFactor;
				sFactor->vIds=new int[pFactor->vCnt+1];
				sFactor->vCnt=pFactor->vCnt+1;
				int fIter=0;
				int dIter=0;
				while((fIter<pFactor->vCnt) && (pFactor->vIds[fIter]<newVId))
				{
					sFactor->vIds[dIter]=pFactor->vIds[fIter];
					fIter++;
					dIter++;
				}
				sFactor->vIds[dIter]=newVId;
				dIter++;
				while(fIter<pFactor->vCnt) 
				{
					sFactor->vIds[dIter]=pFactor->vIds[fIter];
					fIter++;
					dIter++;
				}
				sFactor->fId=globalFactorID;
				string fKey;
				//Get the markov blanket variables 
				for(int i=0;i<sFactor->vCnt;i++)
				{
					SlimFactor* qFactor=canonicalFactorSet[sFactor->vIds[i]];
					for(INTINTMAP_ITER mIter=qFactor->mergedMB.begin();mIter!=qFactor->mergedMB.end();mIter++)
					{
						if(mIter->first!=sFactor->vIds[0])
						{
							sFactor->mergedMB[mIter->first]=mIter->second;
						}
					}
					if(i>0)
					{
						sFactor->mergedMB[qFactor->fId]=0;
					}
				}
				getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
				factorNameToIDMap[fKey]=globalFactorID;
				factorIDToNameMap[globalFactorID]=fKey;
				canonicalFactorSet[globalFactorID]=sFactor;
				globalFactorID++;
				//Here we need to add the super-set and sub-set relationships
				//Specifically, sFactor is a super-set of pFactor
				//pFactor is a sub-set of sFactor
				//sFactor->secondPId=pFactor->fId;
				addToLattice(sFactor,subsetSpace);
				tempParentIDs[sFactor->fId]=0;
			}
		}
		cout <<"Added new " << tempParentIDs.size() << " factors of size "<< fSize << endl;
		fSize++;
		parentIDs.clear();
		for(INTINTMAP_ITER pIter=tempParentIDs.begin();pIter!=tempParentIDs.end();pIter++)
		{
			parentIDs[pIter->first]=0;
		}
		tempParentIDs.clear();
	}
	for(int i=0;i<variableSet.size();i++)
	{
		delete [] subsetSpace[i];
	}
	delete[] subsetSpace;
	cout << "Global factor id " << globalFactorID << endl;
	return 0;
}

//In this function we are going to associate with each canonical parameter, a joint potential table.
//Each entry is a pair corresponding to a configuration of variables of the canonical factor and its value.
//The value is estimated by taking the exponent of the sum over all subsets using the default instantiation.
int
FactorManager::estimateCanonicalParameters(const char* estMethod)
{
	VSET& variableSet=vMgr->getVariableSet();
	char aFName[256];
	sprintf(aFName,"%s/%s.txt",outputDir,estMethod);
	ofstream oFile(aFName);
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactorSet.begin();aIter!=canonicalFactorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		INTINTMAP allSubsets;
		lattice.getAllSubsets(aIter->first,allSubsets);
		allSubsets[sFactor->fId]=0;
		if(strstr(estMethod,"joint")!=NULL)
		{
			potMgr->estimateCanonicalPotential_Joint(sFactor,variableSet,defaultInstantiation,allSubsets,canonicalFactorSet);
		}
		else if(strstr(estMethod,"markovblnkt")!=NULL)
		{
			potMgr->estimateCanonicalPotential(sFactor,variableSet,defaultInstantiation,allSubsets,canonicalFactorSet);
		}
		else if(strstr(estMethod,"approx")!=NULL)
		{
			potMgr->estimateCanonicalPotential_Approximate(sFactor,variableSet,defaultInstantiation,allSubsets,canonicalFactorSet);
		}
		oFile <<"Canonical potential for ID: " << aIter->first <<" ";
		sFactor->canonicalParams->dumpPotential(oFile);
	}
	return 0;
}

int 
FactorManager::learnStructure()
{
	Error::ErrorCode err=estimateClusterProperties();
	if(err!=Error::SUCCESS)
	{
		cout <<Error::getErrorString(err) << endl;
		return 0;
	}
	return 0;
}

//In this function, we go over all good slim factors and prune out the ones that have already been found at higher levels.
int
FactorManager::getMaximalClusters()
{
	//Start from the smallest clusters
	map<int,SlimFactor*>::iterator factorIter=goodSlimFactors.begin();
	while(factorIter!=goodSlimFactors.end())
	{
		int fid=factorIter->first;
		INTINTMAP* supersets=lattice.getSupersets(fid);
		if(supersets!=NULL)
		{
			//If any of the supersets of fid exist in goodSlimFactors then get rid of fid
			//This is because, the dependencies in fid are captured by the superset of fid
			INTINTMAP_ITER iIter=supersets->begin();
			bool foundSuper=false;
			while((iIter!=supersets->end()) && (!foundSuper))
			{
				if(goodSlimFactors.find(iIter->first)!=goodSlimFactors.end())
				{
					foundSuper=true;
				}
				iIter++;
			}
			if(foundSuper)
			{
				goodSlimFactors.erase(factorIter);
			}
		}
		factorIter++;
	}
	return 0;	
}

int
FactorManager::removeDupFactors()
{
	INTINTMAP factorVars;
	string factorVarKey;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		//Need to consider only those factors that may have a Markov blanket 
		if(fIter->second->vCnt==maxFactorSize_Approx)
		{
			break;
		}
		factorVars.clear();
		factorVarKey.clear();
		SlimFactor* sFactor=fIter->second;
		for(int i=0;i<sFactor->vCnt;i++)
		{
			factorVars[sFactor->vIds[i]]=0;
		}
		for(INTINTMAP_ITER mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
		{
			factorVars[mbIter->first]=0;
		}
		char keyPart[32];
		//First append the factor size
		sprintf(keyPart,"%d:",sFactor->vCnt);
		factorVarKey.append(keyPart);
		for(INTINTMAP_ITER vIter=factorVars.begin();vIter!=factorVars.end();vIter++)
		{
			if(vIter!=factorVars.begin())
			{
				sprintf(keyPart,"-%d",vIter->first);
			}
			else
			{
				sprintf(keyPart,"%d",vIter->first);
			}
			factorVarKey.append(keyPart);
		}
		INTDBLMAP* factorSet=NULL;
		if(signFactGrpMap.find(factorVarKey)==signFactGrpMap.end())
		{
			factorSet=new INTDBLMAP;
			signFactGrpMap[factorVarKey]=factorSet;
		}
		else
		{
			factorSet=signFactGrpMap[factorVarKey];
		}
		(*factorSet)[fIter->first]=fIter->second->mbScore;
		factSignMap[fIter->first]=factorVarKey;
	}
	//Now get the non-redundant factors per signature
	char aFName[1024];
	sprintf(aFName,"%s/factsign_k%d.txt",outputDir,maxFactorSize_Approx);
	VSET& variableSet=vMgr->getVariableSet();
	ofstream oFile(aFName);
	oFile <<"GrpSign\tGrpSize\tRepFid\tFactor\tMBlanket" << endl;
	
	for(map<string,INTDBLMAP*>::iterator fSetIter=signFactGrpMap.begin();fSetIter!=signFactGrpMap.end();fSetIter++)
	{
		INTDBLMAP* factorSet=fSetIter->second;
		double minEntropy=10000;
		int minFid=-1;
		SlimFactor* minFactor=NULL;
		for(INTDBLMAP_ITER idIter=factorSet->begin();idIter!=factorSet->end();idIter++)
		{
			double fEntropy=idIter->second;
			if(fEntropy<minEntropy)
			{
				minEntropy=fEntropy;
				minFactor=slimFactorSet[idIter->first];
			}
		}
		if(minFactor==NULL)
		{
			cout <<"No factor for sign " << fSetIter->first.c_str() << endl;
			return -1;
		}
		groupFactorRepMap[fSetIter->first]=minFactor->fId;
		oFile << fSetIter->first.c_str() <<"\t" << fSetIter->second->size() << "\t" << minFactor->fId << "\t";
		minFactor->showFactor(oFile,variableSet,false);
		for(INTINTMAP_ITER mbIter=minFactor->mergedMB.begin();mbIter!=minFactor->mergedMB.end();mbIter++)
		{
			if(mbIter!=minFactor->mergedMB.begin())
			{
				oFile <<"-";
			}
			else
			{
				oFile <<"\t";
			}
			oFile << variableSet[mbIter->first]->getName();
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}

double
FactorManager::getLikelihood_ChainRule()
{
	VSET& variableSet=vMgr->getVariableSet();
	double dataLL=0;
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactorSet.begin();aIter!=canonicalFactorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		double factorLL=potMgr->getLikelihood(sFactor,variableSet);
		dataLL=factorLL+dataLL;
	}

	return dataLL;
}

double 
FactorManager::getLikelihood_ChainRule(FactorGraph* fg)
{
	VSET& variableSet=vMgr->getVariableSet();
	double dataLL=0;
	int factorCnt=fg->getFactorCnt();
	struct timeval begintime;
	gettimeofday(&begintime,NULL);
	
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r,begintime.tv_usec);
	vector<double> lls;
	double meanLL=0;
	for(int j=0;j<50;j++)
	{
		double step=1.0/(double)factorCnt;
		map<int,int> usedInit;
		vector<int> fIds;
		dataLL=0;
		for(int i=0;i<factorCnt;i++)
		{
			double rVal=gsl_ran_flat(r,0,1);
			int rind=(int)(rVal/step);
			while(usedInit.find(rind)!=usedInit.end())
			{
				rVal=gsl_ran_flat(r,0,1);
				rind=(int)(rVal/step);
			}
			usedInit[rind]=0;
			fIds.push_back(rind);
		}
		map<int,int> visitedVertices;
		for(int i=0;i<fIds.size();i++)
		{
			SlimFactor* sFactor=fg->getFactorAt(fIds[i]);
			double factorLL=potMgr->getLikelihood(sFactor,variableSet,visitedVertices);
			dataLL=factorLL+dataLL;
			visitedVertices[sFactor->fId]=0;
		}
		meanLL=meanLL+dataLL;
		lls.push_back(dataLL);
	}
	meanLL=meanLL/((double) lls.size());
	double stdLL=0;
	for(int i=0;i<lls.size();i++)
	{
		double diff=meanLL-lls[i];
		stdLL=stdLL+(diff*diff);
	}
	stdLL=stdLL/(double(lls.size()-1));
	stdLL=sqrt(stdLL);
	cout <<"Mean: " << meanLL <<" Std: " << stdLL << endl;
	dataLL=meanLL;
	gsl_rng_free(r);
	return dataLL;
}

double
FactorManager::getLikelihood()
{
	double dataLL=0;
	int evidCnt=evMgr->getNumberOfEvidences();
	INTDBLMAP potScore;
	VSET& varSet=vMgr->getVariableSet();
	int maxCanSize=0;
	int incorrectLL=0;
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		string aConfStr;
		char confStr[CONSTR_LEN];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sprintf(confStr,"-%d=%d-",eIter->first,eIter->second->getHardEvidVal());
			aConfStr.append(confStr);
		}
		double dataPtll=0;
		for(map<int,SlimFactor*>::iterator fIter=canonicalFactorSet.begin();fIter!=canonicalFactorSet.end();fIter++)
		{
			SlimFactor* sFactor=fIter->second;
			if(sFactor->vCnt==1)
			{
				if(sFactor->mergedMB.size()>maxCanSize)
				{
					maxCanSize=sFactor->mergedMB.size();
				}
			}
			string confKey;
			for(int v=sFactor->vCnt-1;v>=0;v--)
			{
				int vId=sFactor->vIds[v];
				Evidence* evid=(*evidMap)[vId];
				char confStr[CONSTR_LEN];
				if(evid->getType()!=Evidence::HARD)
				{
					cout <<"Handling only hard evidence now " << endl;
					return -1;
				}
				sprintf(confStr,"-%d=%d-",vId,evid->getHardEvidVal());	
				confKey.append(confStr);
			}
			if(sFactor->canonicalParams==NULL)
			{
				cout <<"No canonical params for factor " << endl;
				sFactor->showFactor(cout,varSet);
				return -1;
			}
			double pVal=sFactor->canonicalParams->getJointPotValueForConf(confKey);
			if(pVal==-1)
			{
				cout <<"No value for " << confKey.c_str() << endl;
				return -1;
			}
			if(potScore.find(fIter->first)==potScore.end())
			{
				potScore[fIter->first]=log(pVal);
			}
			else
			{
				potScore[fIter->first]=potScore[fIter->first]+log(pVal);
			}
			dataPtll=dataPtll+log(pVal);
		}
		dataPtll=dataPtll+log(defProb);
		if(dataPtll>0)
		{
			incorrectLL++;
			cout <<"P value " << dataPtll << " for data point " << i << " is greater than 1 " << endl;
			cout << aConfStr.c_str() << endl;
		}
		dataLL=dataLL+dataPtll;
	}
	cout <<"Total Number of data points with log likelihood > 0 " << incorrectLL << endl;
	char aFName[1024];
	sprintf(aFName,"%s/canonical_score_k%d.txt",outputDir,maxCanSize);
	ofstream oFile(aFName);
	for(INTDBLMAP_ITER iIter=potScore.begin();iIter!=potScore.end();iIter++)
	{
		oFile <<iIter->first << "\t" << factorIDToNameMap[iIter->first] <<"\t" <<iIter->second << "\t";
		SlimFactor* sFactor=canonicalFactorSet[iIter->first];
		for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
		{
			if(mIter!=sFactor->mergedMB.begin())
			{
				oFile <<"-";
			}
			oFile << mIter->first;
		}
		oFile << endl;
	}

	oFile.close();

	return dataLL;
}

double 
FactorManager::getLikelihood_MCMC(FactorGraph* fg)
{
	double dataLL=0;
	potMgr->resetPotFuncs();
	double z=estimatePartitionFunction(fg);
	if(z==-1)
	{
		return -1;
	}
	int incorrectLL=0;
	VSET& varSet=vMgr->getVariableSet();
	int evidCnt=evMgr->getNumberOfEvidences();
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		INTINTMAP sample;
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sample[eIter->first]=eIter->second->getHardEvidVal();
		}
		double dataptLL=potMgr->getSampleLikelihood(fg->getAllFactors(),varSet,&sample);
		z=z+exp(dataptLL);
	}
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		INTINTMAP sample;
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sample[eIter->first]=eIter->second->getHardEvidVal();
		}
		double dataptLL=potMgr->getSampleLikelihood(fg->getAllFactors(),varSet,&sample);
		dataptLL=dataptLL-log(z);
		if(dataptLL>0)
		{
			incorrectLL++;
		}
		dataLL=dataptLL+dataLL;
	}
	cout <<"Total Number of data points with log likelihood (MCMC) > 0 " << incorrectLL << endl;

	return dataLL;
}

double 
FactorManager::estimatePartitionFunction(FactorGraph* fg)
{
	double z=0;
	int burnIn=50000;
	int sampleCnt=200000;
	int sampleID=0;
	struct timeval begintime;
	gettimeofday(&begintime,NULL);
	
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r,begintime.tv_usec);

	INTINTMAP* firstSample=new INTINTMAP;
	INTINTMAP* nextSample=new INTINTMAP;
	INTINTMAP* currSample=firstSample;
	randInitSample(currSample,r);
	VSET& varSet=vMgr->getVariableSet();
	while(sampleID <= sampleCnt)
	{
		if(sampleID>burnIn)
		{
			double sampleLikelihood=potMgr->getSampleLikelihood(fg->getAllFactors(),varSet,currSample);
			if(sampleLikelihood==-1)
			{
				return -1;
			}
			z=z+exp(sampleLikelihood);
		}
		sampleID++;
		if(getNextSample(fg,currSample,nextSample,r)==-1)
		{
			return -1;
		}
		currSample=nextSample;
	}
	return z;
}

int
FactorManager::randInitSample(INTINTMAP* sample,gsl_rng* r)
{
	VSET& variableSet=vMgr->getVariableSet();
	
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		INTVECT& values=vIter->second->getValues();
		double step=1.0/(double) values.size();
		int valInd=0;
		double total=step;
		double rVal=gsl_ran_flat(r,0,1);
		while((rVal > total) && (valInd<values.size()))
		{
			total=total+step;
			valInd++;
		}
		if(valInd==values.size())
		{
			cout <<"Could not find sample value " << endl;
			return -1;
		}
		int vval=values[valInd];
		(*sample)[vIter->first]=vval;
	}
	return 0;
}

int 
FactorManager::getNextSample(FactorGraph* fg,INTINTMAP* currSample,INTINTMAP* nextSample,gsl_rng* r)
{
	map<int,SlimFactor*>& factorSet=fg->getAllFactors();
	VSET& varSet=vMgr->getVariableSet();
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		int varSample=potMgr->getVariableSample(*currSample,varSet,fIter->first,fIter->second,r);
		if(varSample==-1)
		{
			return -1;
		}
		(*nextSample)[fIter->first]=varSample;
	}
	return 0;
}

//We want to delete all its supersets  of this factor that will not be needed anymore
int
FactorManager::deleteMySupersets(int fId)
{
	map<int,SlimFactor*>::iterator fIter=slimFactorSet.find(fId);
	if(fIter==slimFactorSet.end())
	{
		cout <<"No factor with id fId " << fId << endl;
		return -1;
	}
	SlimFactor* sFactor=fIter->second;
	INTINTMAP* superSets=lattice.getSupersets(fId);
	if(superSets!=NULL)
	{
		INTINTMAP deletedSupersetIDs;
		for(INTINTMAP_ITER sIter=superSets->begin();sIter!=superSets->end();sIter++)
		{
			int sId=sIter->first;
			map<int,SlimFactor*>::iterator gIter=slimFactorSet.find(sId);
			if(gIter==slimFactorSet.end())
			{
				cout <<"No super factor with id " << sId << " for factor " << fId << endl;
				continue;
			}
			SlimFactor* ssFactor=gIter->second;
			map<int,int> subsetIDs;
			int actSScnt=getActiveSubsetCnt(ssFactor,subsetIDs);
			if((ssFactor->refCnt==0)  && (actSScnt==0))
			//if(ssFactor->refCnt==0)
			{
				if(deleteFactor(ssFactor)==0)
				{
					deletedSupersetIDs[sIter->first]=0;
					for(map<int,int>::iterator aIter=subsetIDs.begin();aIter!=subsetIDs.end();aIter++)
					{
						if(aIter->first==fId)
						{
							continue;
						}
						INTINTMAP* otherSuperSets=lattice.getSupersets(aIter->first);
						if(otherSuperSets==NULL)
						{
							continue;
						}
						INTINTMAP_ITER bIter=otherSuperSets->find(sIter->first);
						if(bIter==otherSuperSets->end())
						{
							continue;
						}
						otherSuperSets->erase(bIter);
					}
				}
			}	
			subsetIDs.clear();
		}
		for(INTINTMAP_ITER dIter=deletedSupersetIDs.begin();dIter!=deletedSupersetIDs.end();dIter++)
		{
			INTINTMAP_ITER sIter=superSets->find(dIter->first);
			if(sIter==superSets->end())
			{
				cerr <<"This factor " << dIter->first << " was a superset of " << sFactor->fId << " but not anymore " << endl;
				exit(-1);
			}
			superSets->erase(sIter);
		}
		deletedSupersetIDs.clear();
	}
	if(sFactor->refCnt==0) 
	{
		if((superSets!=NULL) && (superSets->size()==0))
		{
			deleteFactor(sFactor);
		}
	}
	return 0;
}

int
FactorManager::deleteFactor(SlimFactor* sFactor)
{
	if(sFactor->vCnt==1)
	{
		cout <<"Trying to delete a singleton! factor" << sFactor->fId << endl;
		return -1;
	}
	if(sFactor->fId==9498)
	{
		//cout <<"Stop here " << endl;
	}

	map<int,string>::iterator aIter=factorIDToNameMap.find(sFactor->fId);
	if(aIter==factorIDToNameMap.end())
	{
		cout <<"No factor with id "<< sFactor->fId << endl;
		return -1;
	}
		
	string& factorKey=aIter->second;
	map<string,int>::iterator bIter=factorNameToIDMap.find(factorKey);
	if(bIter==factorNameToIDMap.end())
	{
		cout <<"No factor with name " << factorKey.c_str() << endl;
		return -1;
	}
	factorIDToNameMap.erase(aIter);
	factorNameToIDMap.erase(bIter);
	map<int,SlimFactor*>::iterator fIter=slimFactorSet.find(sFactor->fId);
	slimFactorSet.erase(fIter);
	deleteFromLattice(sFactor->fId);
	delete sFactor;
	return 0;
}

int
FactorManager::getActiveSubsetCnt(SlimFactor* sFactor, map<int,int>& subsetIDs)
{
	int** subsetSpace=new int*[MAXFACTORSIZE_ALLOC];
	for(int i=0;i<MAXFACTORSIZE_ALLOC;i++)
	{
		subsetSpace[i]=new int[MAXFACTORSIZE_ALLOC-1];
	}
	VSET& variableSet=vMgr->getVariableSet();
	sFactor->generateMaximalSubsets(subsetSpace);
	int foundSsets=0;
	for(int sscnt=0;sscnt<sFactor->vCnt;sscnt++)
	{
		string ssKey;
		getFactorKey(subsetSpace[sscnt],sFactor->vCnt-1,ssKey);
		if(factorNameToIDMap.find(ssKey)==factorNameToIDMap.end())
		{
			continue;
		}
		int sId=factorNameToIDMap[ssKey];
		subsetIDs[sId]=0;
		SlimFactor* cFactor=slimFactorSet[sId];
		if(cFactor->refCnt==0)
		{
			continue;
		}
		foundSsets++;
	}
	for(int i=0;i<MAXFACTORSIZE_ALLOC;i++)
	{
		delete[] subsetSpace[i];
	}
	
	delete[] subsetSpace;
	return foundSsets;
}

int
FactorManager::updateRefCnt(int mbId,int status)
{
	map<int,SlimFactor*>::iterator mIter=slimFactorSet.find(mbId);
	if(mIter==slimFactorSet.end())
	{
		return -1;
	}
	SlimFactor* mbFactor=mIter->second;
	mbFactor->refCnt=status;
	return 0;
}

int
FactorManager::clearLattice(INTINTMAP& mbFactors)
{
	map<int,SlimFactor*> delFactors;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		if(fIter->second->vCnt==1)
		{
			continue;
		}
		if(mbFactors.find(fIter->first)!=mbFactors.end())
		{
			continue;
		}
		delFactors[fIter->first]=fIter->second;
	}
	cout <<"Deleting " << delFactors.size() << " factors " << endl;
	for(map<int,SlimFactor*>::iterator dIter=delFactors.begin();dIter!=delFactors.end();dIter++)
	{
		deleteFactor(dIter->second);
	}
	//Reset global factorid
	globalFactorID=mbFactors.rbegin()->first+1;
	delFactors.clear();
	lattice.clear();
	potMgr->clearJointEntropies();
	return 0;
}

int 
FactorManager::readRestrictedVarlist(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	VSET& varSet=vMgr->getVariableSet();
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string varName(buffer);
		int varID=vMgr->getVarID(varName.c_str());
		if(varID==-1)
		{
			continue;
		}
		restrictedNeighborList[varID]=varSet[varID];;
	}
	inFile.close();
	return 0;
}

map<int,Variable*>& 
FactorManager::getRestrictedVarlist()
{
	return restrictedNeighborList;
}

bool 
FactorManager::checkMonotonicity()
{
	bool monotonic=true;
	for(map<int,SlimFactor*>::iterator aIter=slimFactorSet.begin();aIter!=slimFactorSet.end();aIter++)
	{
		int predVCnt=aIter->second->vCnt;
		map<int,SlimFactor*>::iterator bIter=aIter;
		bIter++;
		for(;bIter!=slimFactorSet.end();bIter++)
		{
			int currVCnt=bIter->second->vCnt;
			if(currVCnt<predVCnt) 
			{
				monotonic=false;
				cout <<"Monotonic property violated at prev factor " << aIter->first <<
					" and succ factor " << bIter->first << endl;
			}
			break;
		}
		if(!monotonic)
		{
			break;
		}
	}
	return monotonic;
}

INTINTMAP*
FactorManager::getSupersets(int fId)
{
	return lattice.getSupersets(fId);
}

int
FactorManager::getSupersets(int fId, int level, INTINTMAP& superSetID)
{
	lattice.getSupersets(fId,superSetID,level);
	return 0;
}

SlimFactor*
FactorManager::getFactorAt(int fId)
{
	if(slimFactorSet.find(fId)==slimFactorSet.end())
	{
		return NULL;
	}
	return slimFactorSet[fId];
}

//This function belongs to a bottom up approach of making the Markov blankets consistent.
//We want to make sure that once we are done with the consistency check at one level l, we are guaranteed that every
//factor below l is consistent with the factors at l. 
//For example, consider the factor X and its possible 1-extensions, which can
//be given to us by the superset matrix. Let Y=X\Union Xi be one such extension. To make
//the Markov blankets of X consistent with Y, we need to check that the MBx is a subset of
//Mby \Union Y/X. We can either delete a variable from MBx or add a variable to MBy.
//While adding a variable to MBy will not requires downstream of l, deleting a variable
//from MBx will require checks because we need to make sure that does not make the 
//MBx inconsistent with the Markov blankets of the subsets of X. 
//Using a top down approach also has the same problem if we want to add variables
//to the Markov blanket, because then we have to go up the hierarchy to make sure that
//everything is consistent.
//We chose the bottom up approach and allow only the addition of random variables to the Markov
//blankets. We delete only those R.Vs. from X that are local to the Mx. If we encounter a variable
//that must be deleted which is not local to Mx, that is, was added because a subset of X
//had it in its Markov blanket, we simply add the the variable to the superset of X.
int
FactorManager::makeMarkovBlanketConsistent(SlimFactor* aFactor)
{
	INTINTMAP* superset=lattice.getSupersets(aFactor->fId);
	int delOps=0;
	int addOps=0;
	int lowmi_addOps=0;
	if(superset==NULL)
	{
		return 0;
	}
	INTINTMAP_ITER sIter=superset->begin();
	int* tempVids=new int[maxFactorSize_Approx];
	while(sIter!=superset->end())
	{
		SlimFactor* ssFactor=slimFactorSet[sIter->first];
		
		if( ( (ssFactor->vCnt-aFactor->vCnt)>1) || (ssFactor->vCnt==maxFactorSize_Approx))
		{
			break;
		}
		
		for(INTINTMAP_ITER mbvIter=aFactor->mergedMB.begin();mbvIter!=aFactor->mergedMB.end();mbvIter++)
		{
			int vid=mbvIter->first;
			if((ssFactor->mergedMB.find(vid)==ssFactor->mergedMB.end()) && (!ssFactor->isMemberVariable(vid)))	
			{
				//Now we need to decide whether to add the variable to ssFactor Markov blanket
				//or to delete it from sFactor. Let the variables of ssFactor be Y. Let the variable
				//missing from MBy be V.  Before we add to MBy we need to check that I(Y,V) is a valid factor
				int ssid=0;
				while((ssid<ssFactor->vCnt) && (ssFactor->vIds[ssid]<vid ))
				{
					tempVids[ssid]=ssFactor->vIds[ssid];
					ssid++;
				}
				tempVids[ssid]=vid;
				ssid++;
				while(ssid<ssFactor->vCnt)
				{
					tempVids[ssid]=ssFactor->vIds[ssid];
					ssid++;
				}
				int testfId=getFactorIndex(tempVids,ssFactor->vCnt+1);
				if(slimFactorSet.find(testfId)!=slimFactorSet.end())
				{
					//This just indicates that this variable in ssFactor's Markov blanket has been added
					//to ensure consistency of a Markov blanket of a subset of ssFactor
					ssFactor->mergedMB[vid]=ssFactor->vCnt-1;	
					addOps++;
				}
				else
				{
					if(aFactor->mergedMB[vid]==0)
					{
						INTINTMAP_ITER dIter=aFactor->mergedMB.find(vid);
						aFactor->mergedMB.erase(dIter);
						delOps++;
					}
					else
					{
						ssFactor->mergedMB[vid]=ssFactor->vCnt-1;	
						lowmi_addOps++;
					}
				}
			}
		}
		sIter++;
	}
	delete[] tempVids;
	//if(addOps || delOps || lowmi_addOps)
	if(lowmi_addOps)
	{
		cout<< " Factor " << aFactor->fId <<" addOps " << addOps << " delOps " << delOps << " low mutualinfo addOps " << lowmi_addOps << endl;
	}
	return 0;
}

//Here we generate the clusters using the same technique as the lattice cluster generation algorithm
//but we restrict the growth of clusters only to the Markov blankets.
//First of all we do not need to make checks of non-random mutual information for clusters of size <=maxFactorSize
//because we have already filtered for good clusters in filterClustersWithMI
int
FactorManager::produceClusters(double reqConf,int maxClusterSize)
{
	//Start with all factors of size 2
	/*int currK=2;
	//This plays the role of Fk in the pseudo code
	INTVECT parentIds;
	INTVECT newParentIds;
	INTVECT penultimateMaxIds;
	bool latticeCheckLowerLevel=true;

	//The variables of the new factor
	int* newVids=new int[maxClusterSize];
	//All factors of size 2 which are in slimFactorSet can be blindly added into goodSlimFactors because we have already
	//done the mutual information a.k.a. the support test and confidence does not really make sense here.
	//Since the factor ids grow monotonically with the size of the factors, we can simply iterate over slimfactor set
	//and break out when we reach the first factor whose number of variables is greater than 2.
	map<int,SlimFactor*>::iterator sfIter=slimFactorSet.begin();
	while(sfIter!=slimFactorSet.end())
	{
		if(sfIter->second->vCnt==2)
		{
			goodSlimFactors[sfIter->first]=sfIter->second;
			sfIter->second->confidence=1.0;
			parentIds.push_back(sfIter->first);
		}
		sfIter++;	
		if(sfIter->second->vCnt>2)
		{
			continue;
		}
	}
	if(maxFactorSize_Approx==2)
	{
		map<int,SlimFactor*>::iterator sfIter=slimFactorSet.begin();
		while(sfIter!=slimFactorSet.end())
		{
			if(sfIter->second->vCnt==1)
			{
				//This is needed for growing clusters using Apriori
				penultimateMaxIds.push_back(sfIter->first);
				goodSlimFactors[sfIter->first]=sfIter->second;
				sfIter->second->confidence=1.0;
			}
			sfIter++;
			if(sfIter->second->vCnt>1)
			{
				continue;
			}
		}
	}
	int oldClusterCnt=goodSlimFactors.size();
	int newClusterCnt=1;
	//Use this make map to not reanalyze factors
	INTINTMAP generatedIDs;
	while((newClusterCnt>0) && (currK<=maxClusterSize))//No more factors can be added
	{
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		if(currK<maxFactorSize_Approx)
		{
			//These are factors that have enough mutual information but are not supported
			//at the lower level
			int failedAttempts=0;
			int totalAttempts=0;
			for(int p=0;p<parentIds.size();p++)
			{
				SlimFactor* pFactor=goodSlimFactors[parentIds[p]];
				if(pFactor->fId==160)
				{
					cout <<"Stop here " << endl;
				}
				//Use variables from mergedMB to create a cluster a of size currK+1
				for(INTINTMAP_ITER vIter=pFactor->mergedMB.begin();vIter!=pFactor->mergedMB.end();vIter++)
				{
					int newVId=vIter->first;
					//Do merge-sort of the parent variables and newVId
					int fIter=0;
					int dIter=0;
					while((fIter<pFactor->vCnt) && (pFactor->vIds[fIter]<newVId))
					{
						newVids[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					newVids[dIter]=newVId;
					dIter++;
					while(fIter<pFactor->vCnt) 
					{
						newVids[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					int currFid=getFactorIndex(newVids,pFactor->vCnt+1);
					if(slimFactorSet.find(currFid)==slimFactorSet.end())
					{
						cout << "No factor id  " << currFid << endl;
						return -1;
					}
					SlimFactor* sFactor=slimFactorSet[currFid];
					totalAttempts++;
					//If this check is turned off, we will get the same number of clusters
					//at currK as we did before executing this code.
					if(latticeCheckLowerLevel)
					{
						INTINTMAP* subsets=lattice.getSubsets(sFactor->fId);
						//get all subsets and check for confidence
						//We need to check the last sFactor->vCnt subsets as these
						//will be the ones that correspond to the maximal subsets
						map<int,int>::reverse_iterator rIter=subsets->rbegin();
						int sscnt=0;
						double hitConf=0;
						SlimFactor* aSubset=slimFactorSet[rIter->first];
						while((sscnt<sFactor->vCnt) && ((sFactor->vCnt-aSubset->vCnt) ==1 ))
						{
							int ssId=rIter->first;
							map<int,SlimFactor*>::iterator ssIter=goodSlimFactors.find(ssId);
							if(ssIter!=goodSlimFactors.end())
							{
								hitConf=hitConf+ssIter->second->confidence;
							}
							sscnt++;
							rIter++;
							aSubset=slimFactorSet[rIter->first];
						}
						double conf=hitConf/((double) sFactor->vCnt);
						if(conf>=reqConf)
						{
							goodSlimFactors[sFactor->fId]=sFactor;
							sFactor->confidence=conf;
							newParentIds.push_back(sFactor->fId);
						}
						else
						{
							failedAttempts++;
						}
					}
					else //No need to check confidence at the lower levels where
						//we can compute multi-information correctly
					{
						goodSlimFactors[sFactor->fId]=sFactor;
						newParentIds.push_back(sFactor->fId);
					}
				}//Done with one factor
			}//Done with all factors
			cout <<"Failed to add " << failedAttempts << " out of total " << totalAttempts <<" factors at level " << currK+1 << endl;
			if(currK==(maxFactorSize_Approx-1))
			{
				for(int i=0;i<parentIds.size();i++)
				{
					penultimateMaxIds.push_back(parentIds[i]);
				}
			}
			//Now update parentIds from newParentIds
			parentIds.clear();
			for(int i=0;i<newParentIds.size();i++)
			{
				parentIds.push_back(newParentIds[i]);
			}
			newParentIds.clear();
		}//Upto here we can exactly compute the Markov blankets.
		else
		{
			//The number of variables to add to the factor to create a new factor of size currK+1
			int varsToAdd=currK+1-(maxFactorSize_Approx-1);
			//penultimateMaxIds is made up of all factors of size maxFactorSize-1 that have 
			//Markov blankets associated with them.
			for(int p=0;p<penultimateMaxIds.size();p++)
			{
				SlimFactor* pFactor=goodSlimFactors[penultimateMaxIds[p]];
				//Need to make sure that pFactor's Markov blanket has currK+1-pFactor->vCnt R.V.s
				if(pFactor->mergedMB.size() < varsToAdd)
				{
					continue;
				}
				//Now we need to create subsets of the Markov blanket of size 2.
				
				pFactor->genMBSubsets(varsToAdd);
				int mbssStart=pFactor->mbSubsetStartInd[varsToAdd];
				int mbssEnd=pFactor->mbSubsetStartInd[varsToAdd+1];
				for(int mbsid=mbssStart;mbsid<mbssEnd;mbsid++)
				{
					INTINTMAP* mbset=pFactor->mbSubsets[mbsid];
					//Here we have to create a new factor
					SlimFactor* sFactor=new SlimFactor;
					sFactor->vCnt=currK+1;
					sFactor->vIds=new int[currK+1];
					//Now we need to do a merge sort into the variables of sFactor
					INTINTMAP_ITER mIter=mbset->begin();
					int fIter=0;
					int dIter=0;
					while(dIter!=sFactor->vCnt && mIter!=mbset->end() && fIter!=pFactor->vCnt)
					{
						if(mIter->first<pFactor->vIds[fIter])
						{
							sFactor->vIds[dIter]=mIter->first;
							mIter++;
						}
						else
						{
							sFactor->vIds[dIter]=pFactor->vIds[fIter];
							fIter++;
						}
						dIter++;
					}
					while(mIter!=mbset->end())
					{
						sFactor->vIds[dIter]=mIter->first;
						mIter++;
						dIter++;
					}
					while(fIter!=pFactor->vCnt)
					{
						sFactor->vIds[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					sFactor->fId=getFactorIndex(sFactor->vIds,sFactor->vCnt);
					if(generatedIDs.find(sFactor->fId) !=generatedIDs.end())
					{
						delete sFactor;
						continue;
					}
					generatedIDs[sFactor->fId]=0;
					//Allocate memory for subsets
					int** subsets=new int*[sFactor->vCnt];
					for(int i=0;i<sFactor->vCnt;i++)
					{
						subsets[i]=new int[sFactor->vCnt-1];
					}
					sFactor->generateMaximalSubsets(subsets);
					//Now check the confidence and in the meantime store the subset ids
					//to update the lattice structure
					int* ssIds=new int[sFactor->vCnt];
					double hitConf=0;
					for(int sscnt=0;sscnt<sFactor->vCnt;sscnt++)
					{
						int sId=getFactorIndex(subsets[sscnt],sFactor->vCnt-1);
						ssIds[sscnt]=sId;
						map<int,SlimFactor*>::iterator ssIter=goodSlimFactors.find(sId);
						if(ssIter!=goodSlimFactors.end())
						{
							hitConf=hitConf+ssIter->second->confidence;
						}
					}
					double conf=hitConf/((double)sFactor->vCnt);
					if(conf>=reqConf)
					{
						goodSlimFactors[sFactor->fId]=sFactor;
						sFactor->confidence=conf;
						newParentIds.push_back(sFactor->fId);
						//Update the lattice structure
						for(int i=0;i<sFactor->vCnt;i++)
						{
							if(goodSlimFactors.find(ssIds[i])!=goodSlimFactors.end())
							{
								lattice.addSubset(ssIds[i],sFactor->fId);
								lattice.addSuperset(sFactor->fId,ssIds[i]);
							}
						}
					}
					else
					{
						delete sFactor;
					}
					for(int i=0;i<sFactor->vCnt;i++)
					{
						delete[] subsets[i];
					}
					delete[] subsets;
					delete ssIds;
				}
			}
		}
		currK++;
		newClusterCnt=goodSlimFactors.size()-oldClusterCnt;
		oldClusterCnt=goodSlimFactors.size();
		cout <<"Added " << newClusterCnt<< " new clusters of size " << currK<< endl; 
	}*/
	return 0;
}

int
FactorManager::produceClusters_NoDup(double reqConf,int maxClusterSize)
{
	//Start with all factors of size 2
	/*int currK=2;
	//This plays the role of Fk in the pseudo code
	INTINTMAP parentIds;
	INTINTMAP newParentIds;
	INTINTMAP penultimateMaxIds;
	bool latticeCheckLowerLevel=true;

	//The variables of the new factor
	int* newVids=new int[maxClusterSize];
	//All factors of size 2 which are in slimFactorSet can be blindly added into goodSlimFactors because we have already
	//done the mutual information a.k.a. the support test and confidence does not really make sense here.
	//Since the factor ids grow monotonically with the size of the factors, we can simply iterate over slimfactor set
	//and break out when we reach the first factor whose number of variables is greater than 2.
	map<int,SlimFactor*>::iterator sfIter=slimFactorSet.begin();
	while(sfIter!=slimFactorSet.end())
	{
		if(sfIter->second->vCnt==2)
		{
			//Get the group to which this factor belongs
			string& groupSign=factSignMap[sfIter->first];
			int repId=groupFactorRepMap[groupSign];
			goodSlimFactors[sfIter->first]=sfIter->second;
			sfIter->second->confidence=1.0;
			if(parentIds.find(repId)==parentIds.end())
			{
				parentIds[repId]=0;
			}
		}
		if(sfIter->second->vCnt>2)
		{
			break;
		}
		sfIter++;	
	}
	if(maxFactorSize_Approx==2)
	{
		map<int,SlimFactor*>::iterator sfIter=slimFactorSet.begin();
		while(sfIter!=slimFactorSet.end())
		{
			if(sfIter->second->vCnt==1)
			{
				//This is needed for growing clusters using Apriori
				goodSlimFactors[sfIter->first]=sfIter->second;
				sfIter->second->confidence=1.0;
				string& groupSign=factSignMap[sfIter->first];
				int repId=groupFactorRepMap[groupSign];
				if(penultimateMaxIds.find(repId)==penultimateMaxIds.end())
				{
					penultimateMaxIds[repId]=0;
				}
			}
			sfIter++;
			if(sfIter->second->vCnt>1)
			{
				break;
			}
		}
	}
	int oldClusterCnt=goodSlimFactors.size();
	int newClusterCnt=1;
	//Use this make map to not reanalyze factors
	INTINTMAP generatedIDs;
	while((newClusterCnt>0) && (currK<=maxClusterSize))//No more factors can be added
	{
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		if(currK<maxFactorSize_Approx)
		{
			//These are factors that have enough mutual information but are not supported
			//at the lower level
			int failedAttempts=0;
			int totalAttempts=0;
			for(INTINTMAP_ITER pIter=parentIds.begin();pIter!=parentIds.end();pIter++)
			{
				SlimFactor* pFactor=slimFactorSet[pIter->first];
				//Use variables from mergedMB to create a cluster a of size currK+1
				for(INTINTMAP_ITER vIter=pFactor->mergedMB.begin();vIter!=pFactor->mergedMB.end();vIter++)
				{
					int newVId=vIter->first;
					//Do merge-sort of the parent variables and newVId
					int fIter=0;
					int dIter=0;
					while((fIter<pFactor->vCnt) && (pFactor->vIds[fIter]<newVId))
					{
						newVids[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					newVids[dIter]=newVId;
					dIter++;
					while(fIter<pFactor->vCnt) 
					{
						newVids[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					int currFid=getFactorIndex(newVids,pFactor->vCnt+1);
					if(slimFactorSet.find(currFid)==slimFactorSet.end())
					{
						cout << "Factor " << currFid << " has probably been deleted " <<  endl;
						continue;
					//	return -1;
					}
					//Must select a factor for growth only from the nu dup list
					SlimFactor* sFactor=slimFactorSet[currFid];
					totalAttempts++;
					string& groupSign=factSignMap[sFactor->fId];
					int repId=groupFactorRepMap[groupSign];
					
					//If this check is turned off, we will get the same number of clusters
					//at currK as we did before executing this code.
					INTINTMAP* subsets=lattice.getSubsets(sFactor->fId);
					//get all subsets and check for confidence
					//We need to check the last sFactor->vCnt subsets as these
					//will be the ones that correspond to the maximal subsets
					map<int,int>::reverse_iterator rIter=subsets->rbegin();
					int sscnt=0;
					double hitConf=0;
					//However check for subsets in slimFactorSet
					SlimFactor* aSubset=slimFactorSet[rIter->first];
					while((sscnt<sFactor->vCnt) && ((sFactor->vCnt-aSubset->vCnt) ==1 ))
					{
						int ssId=rIter->first;
						map<int,SlimFactor*>::iterator ssIter=goodSlimFactors.find(ssId);
						if(ssIter!=goodSlimFactors.end())
						{
							hitConf=hitConf+ssIter->second->confidence;
						}
						sscnt++;
						rIter++;
						aSubset=slimFactorSet[rIter->first];
					}
					double conf=hitConf/((double) sFactor->vCnt);
					sFactor->confidence=conf;
					if(latticeCheckLowerLevel)
					{
						if(conf>=reqConf)
						{
							goodSlimFactors[sFactor->fId]=sFactor;
							newParentIds[repId]=0;
						}
						else
						{
							failedAttempts++;
						}
					}
					else //No need to check confidence at the lower levels where
						//we can compute multi-information correctly
					{
						goodSlimFactors[sFactor->fId]=sFactor;
						newParentIds[repId]=0;
					}
				}//Done with one factor
			}//Done with all factors
			cout <<"Failed to add " << failedAttempts << " out of total " << totalAttempts <<" factors at level " << currK+1 << endl;
			if(currK==(maxFactorSize_Approx-1))
			{
				for(INTINTMAP_ITER pIter=parentIds.begin();pIter!=parentIds.end();pIter++)
				{
					penultimateMaxIds[pIter->first]=0;
				}
			}
			//Now update parentIds from newParentIds
			parentIds.clear();
			for(INTINTMAP_ITER pIter=newParentIds.begin();pIter!=newParentIds.end();pIter++)
			{
				parentIds[pIter->first]=0;
			}
			newParentIds.clear();
		}//Upto here we can exactly compute the Markov blankets.
		else
		{
			//The number of variables to add to the factor to create a new factor of size currK+1
			int varsToAdd=currK+1-(maxFactorSize_Approx-1);
			int totalAttempts=0;
			int failedAttempts=0;
			//penultimateMaxIds is made up of all factors of size maxFactorSize-1 that have 
			//Markov blankets associated with them.
			for(INTINTMAP_ITER pIter=penultimateMaxIds.begin();pIter!=penultimateMaxIds.end();pIter++)
			{
				SlimFactor* pFactor=slimFactorSet[pIter->first];
				//Need to make sure that pFactor's Markov blanket has currK+1-pFactor->vCnt R.V.s
				if(pFactor->mergedMB.size() < varsToAdd)
				{
					continue;
				}
				//Now we need to create subsets of the Markov blanket of size 2.
				
				pFactor->genMBSubsets(varsToAdd);
				int mbssStart=pFactor->mbSubsetStartInd[varsToAdd];
				int mbssEnd=pFactor->mbSubsetStartInd[varsToAdd+1];
				for(int mbsid=mbssStart;mbsid<mbssEnd;mbsid++)
				{
					INTINTMAP* mbset=pFactor->mbSubsets[mbsid];
					//Here we have to create a new factor
					SlimFactor* sFactor=new SlimFactor;
					sFactor->vCnt=currK+1;
					sFactor->vIds=new int[currK+1];
					//Now we need to do a merge sort into the variables of sFactor
					INTINTMAP_ITER mIter=mbset->begin();
					int fIter=0;
					int dIter=0;
					while(dIter!=sFactor->vCnt && mIter!=mbset->end() && fIter!=pFactor->vCnt)
					{
						if(mIter->first<pFactor->vIds[fIter])
						{
							sFactor->vIds[dIter]=mIter->first;
							mIter++;
						}
						else
						{
							sFactor->vIds[dIter]=pFactor->vIds[fIter];
							fIter++;
						}
						dIter++;
					}
					while(mIter!=mbset->end())
					{
						sFactor->vIds[dIter]=mIter->first;
						mIter++;
						dIter++;
					}
					while(fIter!=pFactor->vCnt)
					{
						sFactor->vIds[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					sFactor->fId=getFactorIndex(sFactor->vIds,sFactor->vCnt);
					if(generatedIDs.find(sFactor->fId) !=generatedIDs.end())
					{
						delete sFactor;
						continue;
					}
					totalAttempts++;
					generatedIDs[sFactor->fId]=0;
					//Allocate memory for subsets
					int** subsets=new int*[sFactor->vCnt];
					for(int i=0;i<sFactor->vCnt;i++)
					{
						subsets[i]=new int[sFactor->vCnt-1];
					}
					sFactor->generateMaximalSubsets(subsets);
					//Now check the confidence and in the meantime store the subset ids
					//to update the lattice structure
					int* ssIds=new int[sFactor->vCnt];
					double hitConf=0;
					for(int sscnt=0;sscnt<sFactor->vCnt;sscnt++)
					{
						int sId=getFactorIndex(subsets[sscnt],sFactor->vCnt-1);
						ssIds[sscnt]=sId;
						map<int,SlimFactor*>::iterator ssIter=goodSlimFactors.find(sId);
						if(ssIter!=goodSlimFactors.end())
						{
							hitConf=hitConf+ssIter->second->confidence;
						}
					}
					double conf=hitConf/((double)sFactor->vCnt);
					if(conf>=reqConf)
					{
						goodSlimFactors[sFactor->fId]=sFactor;
						sFactor->confidence=conf;
						//Update the lattice structure
						for(int i=0;i<sFactor->vCnt;i++)
						{
							if(goodSlimFactors.find(ssIds[i])!=goodSlimFactors.end())
							{
								lattice.addSubset(ssIds[i],sFactor->fId);
								lattice.addSuperset(sFactor->fId,ssIds[i]);
							}
						}
					}
					else
					{
						failedAttempts++;
						delete sFactor;
					}
					for(int i=0;i<sFactor->vCnt;i++)
					{
						delete[] subsets[i];
					}
					delete[] subsets;
					delete ssIds;
				}
			}
			cout <<"Failed to add " << failedAttempts << " out of total " << totalAttempts <<" factors at level " << currK+1 << endl;
		}
		currK++;
		newClusterCnt=goodSlimFactors.size()-oldClusterCnt;
		oldClusterCnt=goodSlimFactors.size();
		cout <<"Added " << newClusterCnt<< " new clusters of size " << currK<< endl; 
	}*/
	return 0;
}

int 
FactorManager::getFactorIndex(int* vIds, int vCnt)
{
	string aKey;
	getFactorKey(vIds,vCnt,aKey);
	int fId=-1;
	if(factorNameToIDMap.find(aKey)!=factorNameToIDMap.end())
	{
		fId=factorNameToIDMap[aKey];
	}
	return fId;
}

int
FactorManager::getFactorKey(int* vIds, int vCnt, string& key)
{
	for(int v=0;v<vCnt;v++)
	{
		char aBuff[56];
		sprintf(aBuff,"-%d",vIds[v]);
		key.append(aBuff);
	}
	return 0;
}

//vId is the variable id which is used to create the clusters of size cSize, 
//from a total of N variables, starting with vId
int
FactorManager::getClusterCnt(int vId, int cSize, int N)
{
	int clusterCnt=0;
	if(cSize==1)
	{
		clusterCnt=1;
	}
	else if(cSize==2)
	{
		clusterCnt=N-vId-1;
	}
	else if(cSize==3)
	{
		int n=N-cSize+1-vId;
		clusterCnt=(n*(n+1))/2;
	}
	else if(cSize>3)
	{
		//This is the sum of all clusters of size cSize-1 starting
		//from vId+1 to N
		for(int newvId=vId+1;newvId<N;newvId++)
		{
			clusterCnt=clusterCnt+getClusterCnt(newvId,cSize-1,N);
		}
	}
	return clusterCnt;
}

int
FactorManager::initFactorSet()
{	
	int vCnt=vMgr->getVariableSet().size();
	cout <<" Number of factors " << vCnt << endl;
	for(int i=0;i<vCnt;i++)
	{
		SlimFactor* sFactor=new SlimFactor;
		slimFactorSet[i]=sFactor;
		sFactor->vIds=new int[1];
		sFactor->vCnt=1;
		sFactor->mutualInfo=0;
		sFactor->jointEntropy=0;
		sFactor->fId=globalFactorID;
		globalFactorID++;
	}
	cout << "Global factor id " << globalFactorID << endl;
	return 0;
}

//Here we simply write to the memory that we have allocated
int
FactorManager::populateFactorSet()
{
	int currFid=0;
	map<int,Variable*>& variableSet=vMgr->getVariableSet();
	for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		SlimFactor* sFactor=slimFactorSet[currFid];
		currFid++;
		sFactor->vIds[sFactor->vCnt-1]=vIter->first;
		string fKey;
		getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
		factorNameToIDMap[fKey]=sFactor->fId;
		factorIDToNameMap[sFactor->fId]=fKey;
	}
	return 0;
}

int
FactorManager::checkFactorIds()
{
	/*for(int i=0;i<factorCnt;i++)
	{
		SlimFactor* sFactor=slimFactorSet[i];
		int checkFid=getFactorIndex(sFactor->vIds,sFactor->vCnt);
		if(checkFid!=i)
		{
			cout <<"Index not calculated properly for ";
			showFactor(sFactor,cout);
			cout <<"Actual index "<< i << " Calculated index " << checkFid << endl;
		}
	}*/
	return 0;
}

Error::ErrorCode
FactorManager::estimateClusterProperties()
{
	struct timeval begintime;
	struct timeval endtime;
	gettimeofday(&begintime,NULL);
	
	Error::ErrorCode err=potMgr->populatePotentialsSlimFactors(slimFactorSet,vMgr->getVariableSet());

	gettimeofday(&endtime,NULL);
	cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;

	if(err!=Error::SUCCESS)
	{
		return err;
	}
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		sFactor->mbScore=sFactor->jointEntropy;
	}
	return Error::SUCCESS;
}

//This function generates all maximal subsets of this factor
//Then it adds the subset and superset relationships
int
FactorManager::addToLattice(SlimFactor* aFactor,int** subsetSpace)
{
	aFactor->generateMaximalSubsets(subsetSpace);
	for(int i=0;i<aFactor->vCnt;i++)
	{
		int subsetId=getFactorIndex(subsetSpace[i],aFactor->vCnt-1);
		lattice.addSubset(subsetId,aFactor->fId);
		//lattice.addSuperset(aFactor->fId,subsetId);
	}
	return 0;
}

int
FactorManager::deleteFromLattice(int factorId)
{
	lattice.deleteFromLattice(factorId);
	return 0;
}

//totalVars is the number of variables in the union of clusterA and clusterB
double 
FactorManager::getOverlap(SlimFactor* clusterA,SlimFactor* clusterB,int& totalVars)
{
	SlimFactor* smallCluster=clusterB;
	SlimFactor* bigCluster=clusterA;
	if(clusterA->vCnt<clusterB->vCnt)
	{
		smallCluster=clusterA;
		bigCluster=clusterB;
	}

	int commonElements=0;
	for(int i=0;i<smallCluster->vCnt;i++)
	{
		int vInd=0;
		bool found=false;
		while((vInd<bigCluster->vCnt) && (!found))
		{
			if(smallCluster->vIds[i]==bigCluster->vIds[vInd])
			{
				found=true;
			}
			vInd++;
		}
		if(found)
		{
			commonElements++;
		}
	}
	totalVars=smallCluster->vCnt+bigCluster->vCnt-commonElements;
	double overlap=((double) commonElements)/ ((double) smallCluster->vCnt);
	return overlap;
}

//Sort all factors in factorIndSet. These are all the factors of a particular
//size that satisfy some criteria
//This must be done in decreasing order of mutual information
int
FactorManager::qsort(int* factorIndSet,int startind, int endind)
{
	if(endind-startind <=1)
	{
		return 0;
	}
	//pivot is the first element
	int pivot=startind;
	for(int i=startind+1;i<endind;i++)
	{
		//If the element at i is less than element at pivot
		//place element before pivot
		if(slimFactorSet[factorIndSet[i]]->mutualInfo < slimFactorSet[factorIndSet[pivot]]->mutualInfo)
		{
			int temp=factorIndSet[i];
			factorIndSet[i]=factorIndSet[pivot];
			factorIndSet[pivot]=temp;
			pivot++;
		}
	}
	//Now make sure everything after the pivot is greater or equal to the pivot
	qsort(factorIndSet,startind,pivot-1);
	qsort(factorIndSet,pivot+1,endind);
	return 0;
}
