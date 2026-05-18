#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"

EvidenceManager::EvidenceManager()
{
	foldCnt=1;
	preRandomizeSplit=false;
	randseed=0;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile(const char* inFName)
{
	ifstream inFile(inFName);
	char* buffer=NULL;
	string buffstr;
	int bufflen=0;
	int lineNo=0;

	// skip the first line (gene headers)
	if(inFile.good())
	{
		getline(inFile,buffstr);
	}

	while(inFile.good())
	{
		getline(inFile,buffstr);

		if(buffstr.length()<=0)
		{
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=buffstr.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());

		//All the evidences for each variable are stored in a map, indexed by the varId
		EMAP* evidMap=new EMAP;
		char* tok=strtok(buffer,"\t");
		//The toks take the form of varid and value

		int vId = 0;
		while(tok!=NULL)
		{
			Evidence* evid = new Evidence;
			evid->assocVariable(vId);
			//double varVal=log(atof(tok));
			double varVal=atof(tok);
			if(isinf(varVal) || isnan(varVal))
			{
				//cout <<"Found nan! " << tok << endl;
				cerr << "Please remove NaNs from the expression data or check the data format. Not a valid number: " << tok << endl;
				exit(-1);
			}
			evid->setEvidVal(varVal);
			(*evidMap)[vId]=evid;
			tok=strtok(NULL,"\t");
			vId++;
		}
		evidenceSet.push_back(evidMap);
		lineNo++;
	}

	inFile.close();

	cout <<"Number of samples read: " << evidenceSet.size() << endl;

	return Error::SUCCESS;
}

//We create a matrix of randomized evidence, where each evidence has some value
//for the random variables. We populate the matrix one random variable at a time.
//We first generate a vector of permuted indices, in the range of 0 to the total
//number of evidences. Then we populate the part of the matrix associated with
//this variable by querying values from the original matrix in the order specified
//by the permuted indices

int
EvidenceManager::randomizeEvidence(gsl_rng* r, VariableManager* vMgr)
{
	//First create all the evidence sets
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=new EMAP;
		randEvidenceSet.push_back(evidMap);
	}
	//Populate variable wise
	VSET& variableSet=vMgr->getVariableSet();
	int* randInds=new int[trainIndex.size()];
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		populateRandIntegers(r,randInds,trainIndex,trainIndex.size());
		int j=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=NULL;
			if(trainIndex.find(i)!=trainIndex.end())
			{
				evidMap=evidenceSet[randInds[j]];
				j++;
			}
			else
			{
				evidMap=evidenceSet[i];
			}
			EMAP* randEvidMap=randEvidenceSet[i];
			Evidence* evid=(*evidMap)[vIter->first];
			(*randEvidMap)[vIter->first]=evid;
		}
	}
	return 0;
}

EMAP*
EvidenceManager::getEvidenceAt(int evId)
{
	if((evId>=evidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return evidenceSet[evId];
}

EMAP*
EvidenceManager::getRandomEvidenceAt(int evId)
{
	if((evId>=randEvidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return randEvidenceSet[evId];
}

int
EvidenceManager::setFoldCnt(int f)
{
	foldCnt=f;
	return 0;
}

int
EvidenceManager::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int
EvidenceManager::splitData(int s)
{
	int testSetSize=evidenceSet.size()/foldCnt;
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1)
	{
		testEndIndex=evidenceSet.size();
	}
	if(foldCnt==1)
	{
		testStartIndex=-1;
		testEndIndex=-1;
	}
	trainIndex.clear();
	testIndex.clear();
	int m=0;
	int* randInds=NULL;
	if(preRandomizeSplit)
	{
		randInds=new int[evidenceSet.size()];
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
		randseed=getpid();
		gsl_rng_set(r,randseed);
		populateRandIntegers(r,randInds,evidenceSet.size());
		gsl_rng_free(r);
		cout <<"Random seed " << randseed << endl;
	}
	for(int i=0;i<evidenceSet.size();i++)
	{
		int eInd=i;
		if(randInds!=NULL)
		{
			eInd=randInds[i];
		}
		if((m>=testStartIndex) && (m<testEndIndex))
		{
			testIndex[eInd]=0;
		}
		else
		{
			trainIndex[eInd]=0;
		}
		m++;
	}
	if(preRandomizeSplit)
	{
		delete[] randInds;
	}
	return 0;
}

INTINTMAP&
EvidenceManager::getTrainingSet()
{
	return trainIndex;
}

INTINTMAP&
EvidenceManager::getTestSet()
{
	return testIndex;
}

int
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds,int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=rind;
	}
	usedInit.clear();
	return 0;
}

int
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds, INTINTMAP& populateFrom, int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	map<int,int> temp;
	INTINTMAP_ITER tIter=populateFrom.begin();
	for(int i=0;i<size;i++)
	{
		int tid=tIter->first;
		temp[i]=tid;
		tIter++;
	}
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=temp[rind];
	}
	usedInit.clear();
	return 0;
}

int
EvidenceManager::populateRandIntegers(gsl_rng* r, vector<int>& randInds,int size, int subsetsize)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<subsetsize;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds.push_back(rind);
	}
	return 0;
}

