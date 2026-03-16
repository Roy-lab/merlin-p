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

//This function is called when we know the topology and only want to estimate parameters
int 
FactorManager::paramEstimation(const char* estMethod)
{
	setBaseInstantiation();
	VSET& variableSet=vMgr->getVariableSet();
	if(allocateFactorSpace_Graph()==-1)
	{
		return -1;
	}
	estimateCanonicalParameters(estMethod);
	double dll=getLikelihood();
	cout <<"Data likelihood canonical " << dll << endl;
	double dll_exact=getLikelihood_ChainRule();
	cout <<"Data likelihood exact " << dll_exact << endl;
	evaluateMarkovBlanket(-1);
	return 0;	
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

int
FactorManager::generateCanonicalFactors(FactorGraph* fg,map<int,SlimFactor*>& canonicalFactors)
{
	VSET& variableSet=vMgr->getVariableSet();
	int** subsetSpace=new int* [variableSet.size()];
	for(int i=0;i<variableSet.size();i++)
	{
		subsetSpace[i]=new int[variableSet.size()-1];
	}

	INTINTMAP parentIDs;
	map<int,SlimFactor*>& allFactors=fg->getAllFactors();
	for(map<int,SlimFactor*>::iterator aIter=allFactors.begin();aIter!=allFactors.end();aIter++)
	{
		SlimFactor* sFactor=slimFactorSet[aIter->first];
		canonicalFactors[aIter->first]=sFactor;
		sFactor->mergedMB.clear();
		SlimFactor* rFactor=aIter->second;
		for(INTINTMAP_ITER vIter=rFactor->mergedMB.begin();vIter!=rFactor->mergedMB.end();vIter++)
		{
			sFactor->mergedMB[vIter->first]=vIter->second;
		}
		parentIDs[aIter->first]=0;
	}
	cout <<"Before making cliques" << endl;
	for(map<int,SlimFactor*>::iterator sIter=canonicalFactors.begin();sIter!=canonicalFactors.end();sIter++)
	{
		string key;
		SlimFactor* sFactor=sIter->second;
		getFactorKey(sFactor->vIds,sFactor->vCnt,key);
		cout<< "Can. factor: " << key.c_str() <<" MB: ";
		for(INTINTMAP_ITER mvIter=sFactor->mergedMB.begin();mvIter!=sFactor->mergedMB.end();mvIter++)
		{
			cout <<" "<< mvIter->first;
		}
		cout << endl;
	}
	//Now we need to make sure that the MBs of all the variables also form cliques
	int moreChange=1;
	int maxMBSize=maxFactorSize_Approx;
	if(maxMBSize>8)
	{
		maxMBSize=8;
	}
	//Hack to get out of loop
	map<int,int> alreadyDeletedFrom;
	while(moreChange)
	{
		moreChange=0;
		for(map<int,SlimFactor*>::iterator cIter=canonicalFactors.begin();cIter!=canonicalFactors.end();cIter++)
		{
			SlimFactor* cFactor=cIter->second;
			INTINTMAP toErase;
			for(INTINTMAP_ITER uIter=cFactor->mergedMB.begin();uIter!=cFactor->mergedMB.end();uIter++)
			{
				if(toErase.find(uIter->first)!=toErase.end())
				{
					continue;
				}
				SlimFactor* uFactor=canonicalFactors[uIter->first];	
				INTINTMAP_ITER vIter=uIter;
				vIter++;
				for(;vIter!=cFactor->mergedMB.end();vIter++)
				{
					SlimFactor* vFactor=canonicalFactors[vIter->first];
					//Now we need to make sure that u is in v's MB and v is in u's MB
					if(uFactor->mergedMB.find(vIter->first)==uFactor->mergedMB.end())
					{
						if((uFactor->mergedMB.size()<maxMBSize) && (vFactor->mergedMB.size()<maxMBSize))
						{
							if((alreadyDeletedFrom.find(vIter->first)==alreadyDeletedFrom.end()) && 
							   (alreadyDeletedFrom.find(uIter->first)==alreadyDeletedFrom.end()))
						   	{
								uFactor->mergedMB[vIter->first]=vIter->second;
								vFactor->mergedMB[uIter->first]=uIter->second;
							}
						}
						else
						{
							if(uFactor->mergedMB.size()>=maxMBSize)
							{
								toErase[uIter->first]=0;
								alreadyDeletedFrom[uIter->first]=0;
							}
							else
							{
								toErase[vIter->first]=0;
								alreadyDeletedFrom[vIter->first]=0;
							}
							moreChange++;
						}
					}
				}
			}
			for(INTINTMAP_ITER eIter=toErase.begin();eIter!=toErase.end();eIter++)
			{
				//If we are eliminating eIter->first from the MB of dFactor
				//we need to make sure that every variable other than the one
				//we are elminating must also not have an edge to the variable
				//we are eliminating
				SlimFactor* dFactor=canonicalFactors[eIter->first];
				INTINTMAP_ITER dIter=cFactor->mergedMB.find(eIter->first);
				INTINTMAP_ITER fIter=dFactor->mergedMB.find(cFactor->fId);
				if(dIter!=cFactor->mergedMB.end())
				{
					cFactor->mergedMB.erase(dIter);
				}
				if(fIter!=dFactor->mergedMB.end())
				{
					dFactor->mergedMB.erase(fIter);
				}
				for(INTINTMAP_ITER uIter=cFactor->mergedMB.begin();uIter!=cFactor->mergedMB.end();uIter++)
				{
					SlimFactor* eFactor=canonicalFactors[uIter->first];
					INTINTMAP_ITER dIter=eFactor->mergedMB.find(eIter->first);
					INTINTMAP_ITER fIter=dFactor->mergedMB.find(uIter->first);
					if(dIter!=eFactor->mergedMB.end())
					{
						eFactor->mergedMB.erase(dIter);
					}
					if(fIter!=dFactor->mergedMB.end())
					{
						dFactor->mergedMB.erase(fIter);
					}

				}
			}
		}

		cout <<"After making " << moreChange<< " changes " << endl;
		for(map<int,SlimFactor*>::iterator sIter=canonicalFactors.begin();sIter!=canonicalFactors.end();sIter++)
		{
			string key;
			SlimFactor* sFactor=sIter->second;
			getFactorKey(sFactor->vIds,sFactor->vCnt,key);
			cout<< "Can. factor: " << key.c_str() <<" MB: ";
			for(INTINTMAP_ITER mvIter=sFactor->mergedMB.begin();mvIter!=sFactor->mergedMB.end();mvIter++)
			{
				cout <<" "<< mvIter->first;
			}
			cout << endl;
		}
	}
	int fSize=2;
	while(parentIDs.size()>0)
	{
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		INTINTMAP tempParentIDs;
		for(INTINTMAP_ITER pIter=parentIDs.begin();pIter!=parentIDs.end();pIter++)
		{
			SlimFactor* pFactor=canonicalFactors[pIter->first];
			for(INTINTMAP_ITER vIter=pFactor->mergedMB.begin();vIter!=pFactor->mergedMB.end();vIter++)
			{
				int newVId=vIter->first;
				//Now add newVId to if newVId is in the Markov blanket of all variables in pFactor
				bool clqChk=true;
				int vid=0;
				while((vid<pFactor->vCnt) && clqChk)
				{
					SlimFactor* qFactor=canonicalFactors[pFactor->vIds[vid]];
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
				string fKey;
				getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
				if(factorNameToIDMap.find(fKey)!=factorNameToIDMap.end())
				{
					delete sFactor;
					int fId=factorNameToIDMap[fKey];
					sFactor=slimFactorSet[fId];
					sFactor->mergedMB.clear();
				}
				else
				{
					potMgr->populateFactor(slimFactorSet,variableSet,sFactor,false);
					factorNameToIDMap[fKey]=globalFactorID;
					factorIDToNameMap[globalFactorID]=fKey;
					sFactor->fId=globalFactorID;
					slimFactorSet[sFactor->fId]=sFactor;
					globalFactorID++;
				}
				canonicalFactors[sFactor->fId]=sFactor;
				//Here we need to add the super-set and sub-set relationships
				//Specifically, sFactor is a super-set of pFactor
				//Get the markov blanket variables 
				for(int i=0;i<sFactor->vCnt;i++)
				{
					SlimFactor* qFactor=canonicalFactors[sFactor->vIds[i]];
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
	/*for(map<int,SlimFactor*>::iterator sIter=canonicalFactors.begin();sIter!=canonicalFactors.end();sIter++)
	{
		string key;
		SlimFactor* sFactor=sIter->second;
		getFactorKey(sFactor->vIds,sFactor->vCnt,key);
		cout<< "Can. factor: " << key.c_str() <<" MB: ";
		for(INTINTMAP_ITER mvIter=sFactor->mergedMB.begin();mvIter!=sFactor->mergedMB.end();mvIter++)
		{
			cout <<" "<< mvIter->first;
		}
		cout << endl;
	}*/
	
	delete[] subsetSpace;
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
FactorManager::getPseudoLikelihood()
{
	double pseudoLL=0;
	VSET& variableSet=vMgr->getVariableSet();
	for(map<int,SlimFactor*>::iterator aIter=slimFactorSet.begin();aIter!=slimFactorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		if(sFactor->goodMBIDs.size()==0)
		{
			continue;
		}
		double factorLL=potMgr->getPseudoLikelihood(sFactor,variableSet,true);	
		pseudoLL=pseudoLL+factorLL;
	}
	cout <<"Pseudolikelihood: " << pseudoLL << endl;
	return pseudoLL;
}

double
FactorManager::getPseudoLikelihood(FactorGraph* fg,bool train)
{
	double pseudoLL=0;
	VSET& variableSet=vMgr->getVariableSet();
	for(int i=0;i<fg->getFactorCnt();i++)
	{
		SlimFactor* sFactor=fg->getFactorAt(i);
		if(sFactor->goodMBIDs.size()==0)
		{
			continue;
		}
		double factorLL=potMgr->getPseudoLikelihood(sFactor,variableSet,train);	
		pseudoLL=pseudoLL+factorLL;
	}
	//cout <<"Pseudolikelihood: " << pseudoLL << endl;
	return pseudoLL;
}

double
FactorManager::getMVGaussianLikelihood(FactorGraph* fg,bool train)
{
	map<int,SlimFactor*>& factorSet=fg->getAllFactors();
	VSET& variableSet=vMgr->getVariableSet();
	double ll=potMgr->getGaussianLikelihood(factorSet,variableSet,train);
	return ll;
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
FactorManager::getLikelihood(FactorGraph* fg)
{
	double dataLL=0;
	INTDBLMAP potScore;
	VSET& varSet=vMgr->getVariableSet();
	int maxCanSize=0;
	int incorrectLL=0;
	//setBaseInstantiation_Variable();
	setBaseInstantiation();
	map<int,SlimFactor*> canonicalFactors;
	generateCanonicalFactors(fg,canonicalFactors);
	potMgr->resetPotFuncs();
	map<string,int> canFactOrder;
	int maxCanFactSize=0;
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactors.begin();aIter!=canonicalFactors.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		string key;
		getFactorKey(sFactor->vIds,sFactor->vCnt,key);
		if(maxCanFactSize<sFactor->vCnt)
		{
			maxCanFactSize=sFactor->vCnt;
		}
		canFactOrder[key]=aIter->first;
	}
	double den=(double) pow(2.0,(double)maxCanFactSize+2);
	double canThreshold=0.1/den;
	int trivialFactors=0;
	//for(map<string,int>::iterator aIter=canFactOrder.begin();aIter!=canFactOrder.end();aIter++)
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactors.begin();aIter!=canonicalFactors.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		INTINTMAP allSubsets;
		lattice.getAllSubsets(aIter->first,allSubsets);
		allSubsets[sFactor->fId]=0;
		potMgr->estimateCanonicalPotential_Approximate(sFactor,varSet,defaultInstantiation,allSubsets,canonicalFactors);
		//cout <<"Canonical potential for: " << sFactor->fId << endl;
		//sFactor->canonicalParams->dumpPotential(cout);
		sFactor->thresholdToOne(canThreshold);
		if(sFactor->allEntriesInsignificant())
		{
			trivialFactors++;
		}
		//potMgr->estimateCanonicalPotential_Abbeel(sFactor,varSet,defaultInstantiation,allSubsets,canonicalFactors);
		//potMgr->estimateCanonicalPotential(sFactor,varSet,defaultInstantiation,allSubsets,canonicalFactors);
	}
	cout <<"Total number of trivial factors " << trivialFactors << endl;
	int evidCnt=evMgr->getNumberOfEvidences();
	double actualDefProb=0;
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactors.begin();aIter!=canonicalFactors.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		string confKey;
		char confStr[CONSTR_LEN];
		for(int v=sFactor->vCnt-1;v>=0;v--)
		{
			int vId=sFactor->vIds[v];
			int vVal=defInstMap[vId];
			sprintf(confStr,"-%d=%d-",vId,vVal);
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
		if(pVal!=1.0)
		{
			cout <<"Bad pvalue " << pVal << " for potential " << confKey.c_str() << endl;
		}
		actualDefProb=actualDefProb+log(pVal);
	}
	double maxUnnormll=0;
	int maxevid=0;
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		string aConfStr;
		char confStr[CONSTR_LEN];
		int match=0;
		int mismatch=0;

		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sprintf(confStr,"-%d=%d-",eIter->first,eIter->second->getHardEvidVal());
			aConfStr.append(confStr);
			if(eIter->second->getHardEvidVal()==defInstMap[eIter->first])
			{
				match++;
			}
			else
			{
				mismatch++;
			}
		}
		double dataPtll=0;
		for(map<int,SlimFactor*>::iterator aIter=canonicalFactors.begin();aIter!=canonicalFactors.end();aIter++)
		{
			SlimFactor* sFactor=aIter->second;
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
			if(potScore.find(aIter->first)==potScore.end())
			{
				potScore[aIter->first]=log(pVal);
			}
			else
			{
				potScore[aIter->first]=potScore[aIter->first]+log(pVal);
			}
			dataPtll=dataPtll+log(pVal);
		}
		if(dataPtll>maxUnnormll)
		{
			maxUnnormll=dataPtll;
			maxevid=i;
		}
		dataPtll=dataPtll+log(defProb);
		if(dataPtll>0)
		{
			incorrectLL++;
		//	cout << "Pvalue>1\t"<<dataPtll <<"\t" << mismatch <<endl;
		//	cout <<"P value " << dataPtll << " for data point " << i << " is greater than 1 " << endl;
		//	cout << aConfStr.c_str() << endl;
		}
		else
		{
		//	cout << "Pvalue<=1\t"<<dataPtll <<"\t" << mismatch <<endl;
		}
		dataLL=dataLL+dataPtll;
	}
	cout <<"Max unnormalized data ll " << maxUnnormll << " for evidence " << maxevid << endl;
	cout <<"Total Number of data points with log likelihood > 0 " << incorrectLL << endl;
	//Corrected datall would be if the -maxUnnormll=log(defProb). 
	//So subtract for dataLL evidCnt*log(defProb), and then 
	//add evidCnt*(-maxUnnormll)
	if(maxUnnormll>fabs(log(defProb)))
	{
		cout <<"Applying correction" << endl;
		dataLL=dataLL-(evidCnt*log(defProb));
		dataLL=dataLL-(evidCnt*maxUnnormll);
	}

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

int
FactorManager::evaluateMarkovBlanket(double eps_support)
{
	char aFName[1024];
	VSET& varSet=vMgr->getVariableSet();
	potMgr->estimateMarginalEntropies(canonicalFactorSet,varSet,false);
	if(eps_support!=-1)
	{
		sprintf(aFName,"%s/mbscore_score_n%d_x%d_e%.1f.txt",outputDir,beamSize,maxFactorSize_Approx,eps_support);
	}
	else
	{
		sprintf(aFName,"%s/mbscore_score_n%d_x%d.txt",outputDir,beamSize,maxFactorSize_Approx);
	}
	ofstream oFile(aFName);

	oFile<< "Var\tMB\tCondEntr\tMargEntr" << endl;
	for(map<int,SlimFactor*>::iterator fIter=canonicalFactorSet.begin();fIter!=canonicalFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		if(sFactor->mbScore==-1)
		{
			Potential* apot=potMgr->getPotential(sFactor->fId);
			if(apot==NULL)
			{
				cout <<"No potential estimated for factor " << sFactor->fId << endl;
				return -1;
			}
			apot->calculateEntropy();
			sFactor->mbScore=apot->getEntropy();
		}
		oFile <<varSet[sFactor->vIds[0]]->getName() <<"\t";
                for(INTINTMAP_ITER vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end();vIter++)
                {
                        if(vIter!=sFactor->mergedMB.begin())
                        {
                                oFile <<"-";
                        }
                        oFile<< varSet[vIter->first]->getName();
                }
                oFile <<"\t" <<sFactor->mbScore <<"\t" <<sFactor->jointEntropy << endl;
	}
	/*oFile <<"Factor\tFactorStrID\tAllVars\tConditionalEntr\tJointEntr\tMutualInfo" << endl;
	for(map<int,SlimFactor*>::iterator fIter=canonicalFactorSet.begin();fIter!=canonicalFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		potMgr->populateFactor(canonicalFactorSet,varSet,sFactor,false);
		showFactor(sFactor,oFile,false);
		oFile <<"\t";
		for(int i=0;i<sFactor->vCnt;i++)
		{
			oFile <<"-"<<sFactor->vIds[i];
		}
		oFile <<"\t"<< varSet[sFactor->vIds[0]]->getName();
		for(INTINTMAP_ITER vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end();vIter++)
		{
			oFile <<"-"<< varSet[vIter->first]->getName();
		}
		oFile<< "\t" << sFactor->mbScore 
			<< "\t" << sFactor->jointEntropy
			<<"\t" << sFactor->mutualInfo << endl;
	}*/
	oFile.close();
	return 0;
}

int 
FactorManager::showConditionalPotentials(FactorGraph* fg)
{
	char aFName[1024];
	sprintf(aFName,"%s/condpot_k%d",outputDir,maxFactorSize_Approx);
	ofstream oFile(aFName);
	map<int,SlimFactor*>& factorSet=fg->getAllFactors();
	VSET& varSet=vMgr->getVariableSet();
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		Potential* apot=new Potential;
		Variable* currVar=varSet[sFactor->fId];
		apot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
		for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
		{
			Variable* aVar=varSet[mIter->first];
			apot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
		apot->potZeroInit();
		potMgr->populatePotential(apot,false);
		apot->makeValidJPD();
		oFile <<"Potential for " << currVar->getName() << "\t" << currVar->getID() << endl;
		apot->dumpPotential(oFile);
		delete apot;
	}
	oFile.close();
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

double 
FactorManager::getMBScore(SlimFactor* sFactor,int supId)
{
	int* mbVars=new int[MAXFACTORSIZE_ALLOC];
	int mbvarCnt=0;
	SlimFactor* mbFactor=slimFactorSet[supId];
	sFactor->getSetDiff(mbFactor,mbVars,mbvarCnt);
	double mbInfo=getMIFromVars(mbVars,mbvarCnt);
	double proposedScore=sFactor->marginalEntropy-mbFactor->mutualInfo;
	proposedScore=proposedScore+mbInfo;
	double penalty=mbPenalty*log(double(mbvarCnt));
	proposedScore=proposedScore+penalty;
	/*if(proposedScore<0)
	{
		cout<<"Negative entropy for factor " << sFactor->fId << " mbfactor: "<< mbFactor->fId << " Setting to zero" << endl; 
		proposedScore=0;
		delete[] mbVars;
		return -1;
	}*/
	delete[] mbVars;
	return proposedScore;
}

double 
FactorManager::getMBScore(SlimFactor* sFactor)
{
	INTINTMAP fVars;
	fVars[sFactor->fId]=0;
	for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
	{
		fVars[mIter->first]=0;
	}
	double condEntropy=potMgr->getConditionalEntropy(sFactor->fId,fVars,vMgr->getVariableSet());
	//double condEntropy=0;
	double proposedScore=condEntropy;
	double penalty=mbPenalty*log(double(fVars.size()-1));
	proposedScore=proposedScore+penalty;
	fVars.clear();
	return proposedScore;
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

int
FactorManager::makeMBMutuallyConsistent()
{
	VSET& varSet=vMgr->getVariableSet();
	for(map<int,SlimFactor*>::iterator aIter=slimFactorSet.begin();aIter!=slimFactorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		int totalAdd=0;
		int totalDel=0;
		if(sFactor->vCnt>1)
		{
			break;
		}
		if(sFactor->fId==4 || sFactor->fId==5)
		{
			cout << "Stop here " << endl;
		}
		INTINTMAP fVars;
		fVars[sFactor->fId]=0;
		for(INTINTMAP_ITER mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
		{
			fVars[mbIter->first]=0;
		}
		double condEntropy=potMgr->getConditionalEntropy(sFactor->fId,fVars,varSet);
		sFactor->mbScore=condEntropy;
		INTINTMAP deleteMBVarID;
		for(INTINTMAP_ITER mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
		{
			if(mbIter->first<sFactor->fId)
			{
				continue;
			}
			SlimFactor* mbFactor=slimFactorSet[mbIter->first];
			if(mbFactor->mergedMB.find(sFactor->fId)!=mbFactor->mergedMB.end())
			{
				continue;
			}
			INTINTMAP mbVars;
			mbVars[mbFactor->fId]=0;
			for(INTINTMAP_ITER vIter=mbFactor->mergedMB.begin();vIter!=mbFactor->mergedMB.end();vIter++)
			{
				mbVars[vIter->first]=0;
			}
			double mbcondEntropy=potMgr->getConditionalEntropy(mbFactor->fId,mbVars,varSet);
			//Now we have to decide if its better to add sFactor's var to mbFactor or to delete mbFactor's 
			//var from sFactor
			INTINTMAP_ITER vIter=fVars.find(mbIter->first);
			fVars.erase(vIter->first);
			double newcondEntropy=potMgr->getConditionalEntropy(sFactor->fId,fVars,varSet);

			mbVars[sFactor->fId]=0;
			double newmbCondEntropy=potMgr->getConditionalEntropy(mbFactor->fId,mbVars,varSet);

			double gain=mbcondEntropy-newmbCondEntropy;
			double loss=newcondEntropy-condEntropy;
			double gainpenalty=mbPenalty*log(((double) mbFactor->mergedMB.size())/((double)mbFactor->mergedMB.size()+1));
			double losspenalty=mbPenalty*log(((double)sFactor->mergedMB.size())/((double)sFactor->mergedMB.size() + 1));
			gain=gain+gainpenalty;
			loss=loss+losspenalty;
			if(((gain > 0) || (loss > 0)) && (mbFactor->mergedMB.size() <maxFactorSize_Approx ))
			{
				mbFactor->mergedMB[sFactor->fId]=0;
				mbFactor->mbScore=newmbCondEntropy;
				fVars[mbFactor->fId]=0;
				totalAdd++;
			}
			else
			{
				deleteMBVarID[mbIter->first]=0;
				sFactor->mbScore=newcondEntropy;
				totalDel++;
			}
		}
		for(INTINTMAP_ITER dIter=deleteMBVarID.begin();dIter!=deleteMBVarID.end();dIter++)
		{
			INTINTMAP_ITER mIter=sFactor->mergedMB.find(dIter->first);
			sFactor->mergedMB.erase(mIter);
		}
		cout <<"Total additions " << totalAdd << " total deletions " << totalDel << endl;
	}
	return 0;
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
FactorManager::getGoodCandidateMarkovBlankets(SlimFactor* sFactor,INTINTMAP* supersets,int currK, double& minEntropy,double marginalEntropy,INTDBLMAP& goodMBs)
{
	if(supersets==NULL)
	{
		return 0;
	}
	if(supersets->size()==0)
	{
		return 0;
	}
	INTINTMAP_ITER fIter=supersets->begin();
	bool doneFlag=false;
	INTDBLMAP entropyMap;
	vector<int> mbIds;
	double currMinEntropy=minEntropy;
	//This is a temporarary array to store the variables of the Markov blanket whenever need arises
	int* tempMBVars=new int[maxFactorSize_Approx];
	//The supersets are ordered such that the smaller ones are before the larger ones
	//This is because we assign ids to the clusters based on the order in which they are created
	while((!doneFlag) &&(fIter!=supersets->end()))
	{
		int supId=fIter->first;
		SlimFactor* supFactor=slimFactorSet[supId];
		if(supFactor->vCnt>currK)
		{
			doneFlag=true;
			break;
		}//Do stage I
		double proposedEntropy=marginalEntropy-supFactor->mutualInfo;
		double penalty=0;
		if(supFactor->vCnt-sFactor->vCnt>1)
		{
			//Here we need to get the multi-information of the variables only in the Markov blanket
			int mbvarCnt=0;
			sFactor->getSetDiff(supFactor,tempMBVars,mbvarCnt);
			
			double mbInfo=getMIFromVars(tempMBVars,mbvarCnt);
			if(mbInfo<0)
			{
				cout << "Error occured while getting mutual information of sub factors" << endl;
				return -1;
			}
			proposedEntropy=proposedEntropy+mbInfo;
			penalty=mbPenalty*log(double(currK-sFactor->vCnt));
		}
		if(proposedEntropy*(1+penalty)<=currMinEntropy)
		{
			if(proposedEntropy<0)
			{
			//	cout<<"Negative entropy for factor " << sFactor->fId << " mbfactor: "<< supFactor->fId << " setting to 0 "<< endl; 
			//	proposedEntropy=0;
				//return -1;
			}
			currMinEntropy=proposedEntropy;
		}
		entropyMap[supFactor->fId]=proposedEntropy;
		//if(mbIds.size()<beamSize)
		if(mbIds.size()<10)
		{
			mbIds.push_back(supFactor->fId);
			//if(mbIds.size()==beamSize)
			if(mbIds.size()==10)
			{
				for(int i=0;i<mbIds.size();i++)
				{
					for(int j=i+1;j<mbIds.size();j++)
					{
						if(entropyMap[mbIds[i]]>entropyMap[mbIds[j]])
						{
							int tempfid=mbIds[i];
							mbIds[i]=mbIds[j];
							mbIds[j]=tempfid;
						}
					}
				}
			}
		}
		else 
		{
			//Compare with the last element
			int i=mbIds.size()-1;
			bool insertFlag=false;
			while(proposedEntropy<entropyMap[mbIds[i]] && i>=0)
			{
				i--;
				insertFlag=true;
			}

			if(insertFlag)
			{
				i++;
				if(mbIds.size()>1)
				{
					//drop off the last element
					mbIds.pop_back();
					//make a copy of the last element
					mbIds.push_back(mbIds[mbIds.size()-1]);
					for(int j=i;j<mbIds.size()-2;j++)
					{
						mbIds[j]=mbIds[j+1];
					}
				}
				mbIds[i]=supFactor->fId;
			}

		}
		fIter++;
	}
	if(currMinEntropy>=minEntropy)
	{
		return 0;
	}
	double newMinEntropy=currMinEntropy;
	/*for(INTDBLMAP_ITER idIter=entropyMap.begin();idIter!=entropyMap.end();idIter++)
	{
		if(idIter->second<=currMinEntropy)
		{
			goodMBs[idIter->first]=idIter->second;
		}
	}*/
	int i=0; 
	while(i<beamSize) 
	{
		if(i>=mbIds.size())
		{
			break;
		}
		double anentropy=entropyMap[mbIds[i]];
		if((newMinEntropy <= anentropy) && (anentropy<minEntropy))
		{
			newMinEntropy=anentropy;
		}
		if(newMinEntropy>=anentropy)
		{
			goodMBs[mbIds[i]]=anentropy;
		}
		i++;
	}
	for(int m=0;m<mbIds.size();m++)
	{
		SlimFactor* supFactor=slimFactorSet[mbIds[m]];
		for(int v=0;v<supFactor->vCnt;v++)
		{
			int nId=supFactor->vIds[v];
			if(!sFactor->isMemberVariable(nId))
			{
				sFactor->candidateNeighbours[nId]=0;
			}
		}
	}
	minEntropy=newMinEntropy;
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

double
FactorManager::getMIFromVars(int* vIds,int vCnt)
{
	string factorKey;
	getFactorKey(vIds,vCnt,factorKey);
	double mi=-1;
	map<string,int>::iterator aIter=factorNameToIDMap.find(factorKey);
	if(aIter!=factorNameToIDMap.end())
	{
		mi=slimFactorSet[aIter->second]->mutualInfo;
	}
	else if(delFactors_MI.find(factorKey)!=delFactors_MI.end())
	{
		mi=delFactors_MI[factorKey];
	}
	else
	{
		if(mbSpecific_MI.find(factorKey)!=mbSpecific_MI.end())
		{
			mi=mbSpecific_MI[factorKey];
		}
		else
		{
			VSET& varSet=vMgr->getVariableSet();
			Potential* apot=new Potential;
			for(int j=0;j<vCnt;j++)
			{
				Variable* aVar=varSet[vIds[j]];
				if(j==vCnt-1)
				{
					apot->setAssocVariable(aVar,Potential::FACTOR);
				}
				else
				{
					apot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
				}
			}
			apot->potZeroInit();
			potMgr->populatePotential(apot,false);
			apot->calculateJointEntropy();
			mi=-1*apot->getJointEntropy();
			for(int j=0;j<vCnt;j++)
			{
				SlimFactor* sFactor=slimFactorSet[vIds[j]];
				mi=mi+sFactor->jointEntropy;
			}
			mbSpecific_MI[factorKey]=mi;
			delete apot;
		}
	}
	return mi;
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
