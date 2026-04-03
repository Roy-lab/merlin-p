#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>

#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

#include "FactorGraph.H"
#include "MetaMove.H"
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"
#include "HyperGeomPval.H"
#include "Distance.H"
#include "MetaLearner.H"

MetaLearner::MetaLearner()
{
	restrictedFName[0]='\0';
	trueGraphFName[0]='\0';
	preRandomizeSplit=false;
	random=false;
	lambda=0;
	clusterThreshold=0.5;
	specificFold=-1;
	convThreshold=1e-3;
	factorGraph=nullptr;
	currPLL=nullptr;
}

MetaLearner::~MetaLearner()
{
}

int
MetaLearner::setMaxFactorSize_Approx(int aVal)
{
	maxFactorSizeApprox=aVal;
	return 0;
}

int 
MetaLearner::setPenalty(double aVal)
{
	penalty=aVal;
	return 0;
}

int
MetaLearner::setBeta1(double aval)
{
	beta1=aval;
	return 0;
}

int
MetaLearner::initEdgePriorMeta_All()
{
	for(map<string,map<string,map<string,double>*>*>::iterator gIter=priorgraphmap.begin();gIter!=priorgraphmap.end();gIter++)
	{
		map<string,map<string,double>*>* priorgraph = gIter->second;
		map<int,INTDBLMAP*>* edgeprior = new map<int,INTDBLMAP*>();
		edgepriormap[gIter->first] = edgeprior;
		initEdgePriorMeta(*priorgraph,*edgeprior);
	}
	return 0;
}

int
MetaLearner::setPriorGraph_All(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string gname;
		string fname;
		double gbeta;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				gname.append(tok);
			}
			else if(tokCnt==1)
			{
				fname.append(tok);
			}
			else if(tokCnt==2)
			{
				gbeta=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		betamap[gname] = gbeta;
		map<string,map<string,double>*>* priorGraph = new map<string,map<string,double>*>();
		setPriorGraph(fname.c_str(),*priorGraph);
		priorgraphmap[gname] = priorGraph;
	}
	inFile.close();
	return 0;
}

int
MetaLearner::setBeta_Motif(double aval)
{
	beta_motif=aval;
	return 0;
}

int
MetaLearner::setLambda(double l)
{
	lambda=l;
	return 0;
}

int 
MetaLearner::setConvergenceThreshold(double aVal)
{
	convThreshold=aVal;
	return 0;
}

int
MetaLearner::setRestrictedList(const char* aFName)
{
	strcpy(restrictedFName,aFName);
	ifstream inFile(restrictedFName);
	string buffer;
	while(inFile.good())
	{
		getline(inFile,buffer);
		if(buffer.length()<=0)
		{
			continue;
		}
		restrictedVarList[buffer]=0;
	}
	inFile.close();
	return 0;
}


int 
MetaLearner::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int
MetaLearner::setGlobalEvidenceManager(EvidenceManager* anEvMgr)
{
	evidenceManager=anEvMgr;
	return 0;
}

int 
MetaLearner::setVariableManager(VariableManager* aPtr)
{
	varManager=aPtr;
	return 0;
}

int
MetaLearner::setOutputDirName(const char* dirPath)
{
	strcpy(outputDirName,dirPath);
	return 0;
}

int 
MetaLearner::setClusteringThreshold(double aVal)
{
	clusterThreshold=aVal;
	return 0;
}
 

int
MetaLearner::setSpecificFold(int fid)
{
	specificFold=fid;
	return 0;
}

int 
MetaLearner::setPriorGraph(const char* aFName, map<string,map<string,double>*>& priorGraph)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string tfName;
		string tgtName;
		double edgeStrength;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				tfName.append(tok);
			}
			else if(tokCnt==1)
			{
				tgtName.append(tok);
			}
			else if(tokCnt==2)
			{
				edgeStrength=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,double>* tgtSet=NULL;
		if(priorGraph.find(tfName)==priorGraph.end())
		{
			tgtSet=new map<string,double>;
			priorGraph[tfName]=tgtSet;
		}
		else
		{
			tgtSet=priorGraph[tfName];
		}
		(*tgtSet)[tgtName]=edgeStrength;
	}
	inFile.close();
	return 0;
}

int
MetaLearner::setRandom(bool flag)
{
	random=flag;
	return 0;
}

int
MetaLearner::readModuleMembership(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string geneName;
		int moduleID;
		int tokCnt=0;
		char* tok=strtok(buffer,"\t");
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else if(tokCnt==1)
			{
				moduleID=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,int>* geneSet=NULL;
		if(moduleGeneSet.find(moduleID)==moduleGeneSet.end())
		{
			geneSet=new map<string,int>;
			moduleGeneSet[moduleID]=geneSet;
		}
		else
		{
			geneSet=moduleGeneSet[moduleID];
		}
		(*geneSet)[geneName]=0;
		geneModuleID[geneName]=moduleID;
	}
	inFile.close();
	return 0;
}

int
MetaLearner::setDefaultModuleMembership()
{
	VSET& varSet=varManager->getVariableSet();
	int vCnt=varSet.size();
	int moduleCnt=(int) sqrt(vCnt/2);
	if(moduleCnt>30)
	{
		moduleCnt=30;
	}
	map<int,int> matIdvIdMap;
	int mID=0;
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		matIdvIdMap[mID]=vIter->first;
		mID++;
	}
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	//Randomly partition the variables into clusterassignments
	vector<int> randIndex;
	double step=1.0/(double)vCnt;
	map<int,int> usedInit;
	int maxind=0;
	for(int i=0;i<vCnt;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
			if(rind>maxind)
			{
				maxind=rind;
			}
		}
		if(matIdvIdMap.find(rind)==matIdvIdMap.end())
		{
			cout <<"Did not find " << rind << " matrixidmap " << endl;
		}
		int dataind=matIdvIdMap[rind];
		usedInit[rind]=0;
		randIndex.push_back(dataind);
	}
	//For each partition estimate the mean and covariance
	int clusterSize=vCnt/moduleCnt;
	for(int e=0;e<moduleCnt;e++)
	{
		int startInd=e*clusterSize;
		int endInd=(e+1)*clusterSize;
		if(e==moduleCnt-1)
		{
			endInd=clusterSize;
		}
		map<string,int>* geneSet=NULL;
		geneSet=new map<string,int>;
		moduleGeneSet[e]=geneSet;
		for(int i=startInd;i<endInd;i++)
		{
			int dataId=randIndex[i];
			Variable* v=varSet[dataId];
			(*geneSet)[v->getName()]=0;
			geneModuleID[v->getName()]=e;
		}
	}
	randIndex.clear();
	matIdvIdMap.clear();
	usedInit.clear();
	return 0;
}

int
MetaLearner::initEdgePriorMeta(map<string,map<string,double>*>& graph, map<int,INTDBLMAP*>& edgePriors)
{
	VSET& varSet=varManager->getVariableSet();
	for(map<string,int>::iterator rIter=restrictedVarList.begin();rIter!=restrictedVarList.end();rIter++)
	{
		int regId=varManager->getVarID(rIter->first.c_str());
		if(regId==-1)
		{
			continue;
		}
		if(graph.find(rIter->first)==graph.end())
		{
			continue;
		}
		int tfhit=0;
		map<string,double>* tgtSet=graph[rIter->first];
		for(map<string,double>::iterator vIter=tgtSet->begin();vIter!=tgtSet->end();vIter++)
		{
			INTDBLMAP* edgePriorGene=NULL;
			int tgtId=varManager->getVarID(vIter->first.c_str());
			if(tgtId==-1)
			{
				continue;
			}
			if(edgePriors.find(tgtId)==edgePriors.end())
			{
				edgePriorGene=new INTDBLMAP;
				edgePriors[tgtId]=edgePriorGene;
			}
			else
			{
				edgePriorGene=edgePriors[tgtId];
			}
			double ewt=fabs(vIter->second);
			if(edgePriorGene->find(regId)==edgePriorGene->end())
			{
				//(*edgePriorGene)[regId]=vIter->second;
				(*edgePriorGene)[regId]=ewt;
			}
			else
			{
				//(*edgePriorGene)[regId]=(*edgePriorGene)[regId]+vIter->second;
				(*edgePriorGene)[regId]=(*edgePriorGene)[regId]+ewt;
			}
			tfhit++;
		}
		cout <<"TF: "<< rIter->first << " TFhits= " << tfhit << endl;
	}

	return 0;
}

int
MetaLearner::doCrossValidation(int foldCnt)
{
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	rnd=gsl_rng_alloc(gsl_rng_default);

	evidenceManager->setVariableManager(varManager);
	evidenceManager->setFoldCnt(foldCnt);
	evidenceManager->splitData(0);

	potManager = new PotentialManager;
	potManager->setEvidenceManager(evidenceManager);

	//The first key is for the fold number
	//For each fold we have a trained model. For each trained model we have the likelihood on 
	//all the test sets, including the self test.
	int foldBegin=0;
	int foldEnd=foldCnt;
	if(specificFold>-1)
	{
		foldBegin=specificFold;
		foldEnd=specificFold+1;
	}
	for(int f=foldBegin;f<foldEnd;f++)
	{	
		evidenceManager->splitData(f);
		if(random)
		{	
			evidenceManager->randomizeEvidence(r);
		}

		potManager->reset();
		potManager->setRandom(random);
		potManager->init();

		factorGraph = new FactorGraph(varManager);

		char outputDir[1024];
		sprintf(outputDir,"%s/fold%d",outputDirName,f);
		char foldOutputDirCmd[1024];
		sprintf(foldOutputDirCmd,"mkdir %s",outputDir);
		system(foldOutputDirCmd);

		start(f);
		getPredictionError_CrossValid(f);
		clearFoldSpecData();
	}
	gsl_rng_free(r);

	gsl_rng_free(rnd);
	return 0;
}

int
MetaLearner::start(int f)
{
	//Repeat until convergence
	//int currK=1;
	currFold=f;
	sprintf(foldoutDirName,"%s/fold%d",outputDirName,f);
	int maxMBSizeApprox=maxFactorSizeApprox-1;
	int currK=maxMBSizeApprox;
	rnd=gsl_rng_alloc(gsl_rng_default);
	int rseed=getpid();
	gsl_rng_set(rnd,rseed);
	cout <<rseed << endl;
	initEdgePriorMeta_All();
	initEdgeSet();
	initPhysicalDegree();
	int i=0;
	VSET& varSet=varManager->getVariableSet();
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		idVidMap[i]=vIter->first;
		i++;
	}
	if(strlen(trueGraphFName)==0)
	{
		double currGlobalScore=getInitPLLScore();
		double initScore=getInitPrior();
		//currGlobalScore=currGlobalScore+initScore;
		int showid=0;
		int moduleiter=0;
		bool notConvergedTop=true;
		vector<int> randOrder;
		while(moduleiter<1 && notConvergedTop)
		{
			int iter=0;
			bool notConverged=true;
			while(notConverged && iter<50)
			{
				int attemptedMoves=0;
				int subiter=0;
				double scorePremodule=currGlobalScore;
				randOrder.clear();
				evidenceManager->populateRandIntegers(rnd,randOrder,varSet.size(),varSet.size());				
				struct timeval begintime;
				struct timeval endtime;
				struct timezone begintimezone;
				struct timezone endtimezone;
				gettimeofday(&begintime,&begintimezone);
				while(subiter<varSet.size())
				//while(notConverged && subiter<6000)
				{
					int rID=randOrder[subiter];
					if(idVidMap.find(rID)==idVidMap.end())
					{
						cout <<"Variable at  " << rID << " just not found " << endl;
						exit(0);
					}
					//int vID=idVidMap[rID];
					int vID=idVidMap[subiter];
					VSET_ITER vIter=varSet.find(vID);
					if(vIter==varSet.end())
					{
						subiter++;
						continue;
					}
					Variable* v=varSet[vID];
					int lastiter=0;
					if(variableStatus.find(v->getName())!=variableStatus.end())
					{
						lastiter=variableStatus[v->getName()];
						if((iter-lastiter)>=5)
						{
							cout <<"Skipping " << v->getName() << endl;
							subiter++;
							continue;	
						}
					}		
					struct timeval begintime_v;
					struct timeval endtime_v;
					struct timezone begintimezone_v;
					struct timezone endtimezone_v;
					collectMoves(currK,vID);
					if(moveSet.size()==0)
					{
						subiter++;
						continue;
					}
					sortMoves();
					makeMoves();
					double newScore=getPLLScore();
					double diff=newScore-currGlobalScore;
					if(diff<=convThreshold)
					{
					//	notConverged=false;
					}
					//dumpAllGraphs(currK,f,iter);
					currGlobalScore=newScore;
					//cout <<"Current iter " << iter << " Score after beta-theta " << newScore << endl;
					subiter++;
					showid++;
					attemptedMoves++;
					gettimeofday(&endtime_v,&endtimezone_v);
					//printf("Time elapsed for one var %uj secs %d microsec\n",(unsigned int)(endtime_v.tv_sec-begintime_v.tv_sec),(unsigned int)(endtime_v.tv_usec-begintime_v.tv_usec));
				}
				gettimeofday(&endtime,&endtimezone);
				//printf("Time elapsed for all vars %d mins %d secs %d microsec\n", (unsigned int)(endtimezone.tz_minuteswest-begintimezone.tz_minuteswest), (unsigned int)(endtime.tv_sec-begintime.tv_sec,endtime.tv_usec-begintime.tv_usec));
				if((currGlobalScore-scorePremodule)<=convThreshold)
				{
					notConverged=false;
				}
				else
				{
					redefineModules();
				}
				iter++;
				scorePremodule=currGlobalScore;
				dumpAllGraphs(currK,f,iter);
			}
			moduleiter++;
		}
		cout <<"Final Score " << currGlobalScore << endl;
		finalScores[f]=currGlobalScore;
	}
	return 0;
}

double
MetaLearner::getInitPLLScore()
{
	double initScore=0;
	VSET& varSet=varManager->getVariableSet();
	//Initially we just sum up the marginal likelihoods
	currPLL=new INTDBLMAP;
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		if(varNeighborhoodPrior.find(vIter->first)==varNeighborhoodPrior.end())
		{
			continue;
		}
		Variable* var=varSet[vIter->first];
		double newPLL_s=getInitPLLScore(vIter->first);
		double priorScore=varNeighborhoodPrior[vIter->first];
		(*currPLL)[vIter->first]=newPLL_s+priorScore;
		initScore=initScore+(*currPLL)[vIter->first];
	}
	return initScore;
}

double
MetaLearner::getPLLScore()
{
	double gScore=0;
	for(INTDBLMAP_ITER dIter=currPLL->begin();dIter!=currPLL->end();dIter++)
	{
		if(isnan(gScore) || isinf(gScore))
		{
			cout << "Found nan/inf for variable " << dIter->first << endl;
		}
		gScore=gScore+dIter->second;
	}
	return gScore;
}

double 
MetaLearner::getInitPrior()
{
	double graphPrior=0;
	double edgePresence=1/(1+exp(-1*beta1));
	for(map<string,double>::iterator aIter=edgePresenceProb.begin();aIter!=edgePresenceProb.end();aIter++)
	{
		//graphPrior=graphPrior+log(1-edgePresence);
		graphPrior=graphPrior+log(1-aIter->second);
		if(isinf(graphPrior)|| isnan(graphPrior))
		{
			cout <<"Graph prior is "<< graphPrior << " after " << aIter->first << " for " << aIter->second << endl;
		}
	}
	return graphPrior;
}

int
MetaLearner::clearFoldSpecData()
{
	if (factorGraph != nullptr)
	{
		delete factorGraph;
		factorGraph = nullptr;
	}
	edgeMap.clear();
	if (currPLL != nullptr)
	{
		delete currPLL;
		currPLL = nullptr;
	}
	return 0;
}

int
MetaLearner::initEdgeSet()
{
	VSET& varSet=varManager->getVariableSet();
	for(VSET_ITER uIter=varSet.begin();uIter!=varSet.end();uIter++)
	{
		Variable* u=varSet[uIter->first];
		if((restrictedVarList.size()>0) && (restrictedVarList.find(u->getName())==restrictedVarList.end()))
		{
			continue;
		}

		for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
		{
			if(uIter->first==vIter->first)
			{
				continue;
			}
			Variable* v=varSet[vIter->first];
			if(geneModuleID.find(v->getName())==geneModuleID.end())
			{	
				continue;
			}
			string edgeKey;
			//This is going to be a directed graph
			edgeKey.append(u->getName().c_str());
			edgeKey.append("\t");
			edgeKey.append(v->getName().c_str());

			edgeMap[edgeKey]=0;

			double initPrior=getEdgePrior(uIter->first,vIter->first);
			initPrior=1/(1+exp(-1*initPrior));
			if(initPrior<1e-6)
			{
				initPrior=1e-6;
			}
			if(initPrior==1)
			{
				initPrior=1-1e-6;
			}
			edgePresenceProb[edgeKey]=initPrior;
			if(varNeighborhoodPrior.find(vIter->first)==varNeighborhoodPrior.end())
			{
				varNeighborhoodPrior[vIter->first]=log(1-initPrior);
			}
			else
			{
				varNeighborhoodPrior[vIter->first]=varNeighborhoodPrior[vIter->first]+log(1-initPrior);
			}
		}
	}
	cout <<"Restricted varlist size: " << restrictedVarList.size() << endl;
	int n=varSet.size();
	int r=restrictedVarList.size();
	int expEdgeCnt=((r*(r-1))/2) + (r*(n-r)) ;
	cout <<"Inited " << edgeMap.size() << " edges. Expected " << expEdgeCnt << endl;

	//Init the potentials
	for(int f=0;f<factorGraph->getFactorCnt();f++)
	{
		SlimFactor* sFactor=factorGraph->getFactorAt(f);
		sFactor->potFunc=new Potential;
		sFactor->potFunc->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
		sFactor->potFunc->potZeroInit();
		potManager->populatePotential(sFactor->potFunc);
		sFactor->potFunc->initMBCovMean();
	}

	return 0;
}

int
MetaLearner::getPredictionError_CrossValid(int foldid)
{
	VSET& varSet=varManager->getVariableSet();
	char foldoutDirName[1024];
	sprintf(foldoutDirName,"%s/fold%d",outputDirName,foldid);
	INTINTMAP& testSet=evidenceManager->getTestSet();
	map<int,double> varPLL;
	for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
	{
		//for each gc, get the expected value of this datapoint
		EMAP* evidMap=evidenceManager->getEvidenceAt(dIter->first);

		for(map<string,int>::iterator vIter=geneModuleID.begin();vIter!=geneModuleID.end();vIter++)
		{
			int vId=varManager->getVarID(vIter->first.c_str());
			if(vId==-1)
			{
				continue;
			}
			Variable* v=varSet[vId];
			double cll=0;
			SlimFactor* sFactor=factorGraph->getFactorAt(vId);
			Potential* sPot=sFactor->potFunc;
			if(sPot==NULL)
			{
				cout <<"Found null for factor="<< sFactor->fId
					<< " variable=" <<varSet[sFactor->fId]->getName() << endl;
			}
			double pval=sPot->getCondPotValueFor(evidMap);
			if(pval<1e-50)
			{
				pval=1e-50;
			}
			if(isinf(pval) || isnan(pval))
			{
				cout <<"Stop here. Found nan/inf for " << vIter->first << " dtpt "<< dIter->first << endl;
			}
			cll=log(pval);
			if(varPLL.find(vId)==varPLL.end())
			{
				varPLL[vId]=cll;
			}
			else
			{
				varPLL[vId]=varPLL[vId]+cll;
			}
		}
	}
	/*
	for(map<int,double>::iterator pIter=varPLL.begin();pIter!=varPLL.end();pIter++)
	{
		oFile << varSet[pIter->first]->getName() << "\t" << pIter->second << endl;
	}
	pFile << "\tRMSE\tNormRMSE\tCoeff_Det_aka_R^2\tCC"<< endl;
	*/
	Distance d;
	vector<double> truevect;
	vector<double> predvect;
	for(map<string,int>::iterator vIter=geneModuleID.begin();vIter!=geneModuleID.end();vIter++)
	{
		int vId=varManager->getVarID(vIter->first.c_str());
		if(vId==-1)
		{
			continue;
		}
		//pFile <<vIter->first;
		double error=0;
		double norm=0;
		double maxval=-100000;
		double minval=1000000;
		double totalvar=0;
		double truemean=0;
		truevect.clear();
		predvect.clear();

		for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
		{
			EMAP* evidMap=evidenceManager->getEvidenceAt(dIter->first);
			Evidence* evid=(*evidMap)[vId];
			double trueval=evid->getEvidVal();
			truemean=truemean+trueval;
			truevect.push_back(trueval);
		}

		truemean=truemean/((double)testSet.size());

		//First the predicted time course
		SlimFactor* sFactor=factorGraph->getFactorAt(vId);
		Potential* sPot=sFactor->potFunc;
		for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
		{
			EMAP* evidMap=evidenceManager->getEvidenceAt(dIter->first);
			double predval=sPot->predictSample(evidMap);
			Evidence* evid=(*evidMap)[vId];
			double trueval=evid->getEvidVal();
			totalvar=totalvar+((trueval-truemean)*(trueval-truemean));
			//also called residuals
			error=error+((predval-trueval)*(predval-trueval));
			predvect.push_back(predval);
			//norm=norm+(trueval*trueval);
			norm=norm+1;
			if(trueval>maxval)
			{
				maxval=trueval;
			}
			if(trueval<minval)
			{
				minval=trueval;
			}
		}

		//Then the true time course
		/*for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
		{
			EMAP* evidMap=evMgr->getEvidenceAt(dIter->first);
			Evidence* evid=(*evidMap)[vId];
			//pFile <<"\t" <<evid->getEvidVal();
		}*/
		double coeff_det=1-(error/totalvar);
		error=error/norm;
		error=sqrt(error);
		double nmsd=error/(maxval-minval);
		double cc=d.computeCC(truevect,predvect);
		//pFile<< "\t"<<error <<"\t" << nmsd << "\t" << coeff_det <<"\t" << cc <<endl;
	}
	//oFile.close();
	//pFile.close();
	varPLL.clear();
	return 0;
}

int
MetaLearner::collectMoves(int currK,int rind)
{
	for(int i=0;i<moveSet.size();i++)
	{
		delete moveSet[i];
	}
	moveSet.clear();

	int vID=idVidMap[rind];
	VSET& varSet=varManager->getVariableSet();
	VSET_ITER vIter=varSet.find(vID);
	if(vIter==varSet.end())
	{
		return 0;
	}

	Variable* v=vIter->second;
	if(geneModuleID.find(v->getName())==geneModuleID.end())
	{
		return 0;
	}

	map<string,int> testedEdges;
	int moduleID=geneModuleID[v->getName()];
	map<string,int>* moduleMembers=moduleGeneSet[moduleID];
	double bestTargetScore=0;
	double bestScoreImprovement=0;
	Variable* bestu=NULL;
	Potential* bestPot=NULL;
	for(map<string,int>::iterator uIter=restrictedVarList.begin();uIter!=restrictedVarList.end();uIter++)
	{
		int regID=varManager->getVarID(uIter->first.c_str());

		// Ensure we can find the regulator, and that it isnt the same node as the target.
		if(regID==-1 || vIter->first==regID)
		{
			continue;
		}

		Variable* u=varSet[regID];

		string edgeKey;
		edgeKey.append(u->getName().c_str());
		edgeKey.append("\t");
		edgeKey.append(v->getName().c_str());

		if(testedEdges.find(edgeKey)!=testedEdges.end())
		{
			continue;
		}
		testedEdges[edgeKey]=0;

		if(edgeMap.find(edgeKey)==edgeMap.end())
		{
			// If the edge key doesnt exist in the map, then something went wrong during initEdgeSet
			cout <<"No edge " << edgeKey.c_str() << " u " << u->getID() << " v " << v->getID()<< endl;
			exit(0);
		}

		// If the edge already exists, no need to test adding it.
		int edgeValue=edgeMap[edgeKey];
		if(edgeValue==1)
		{
			continue;
		}

		// If v already has the max number of parents, dont test adding another.
		if(!checkMBSize(regID,vIter->first,currK))
		{
			continue;
		}

		double improvement=0;
		double score=0;

		Potential* aPot=NULL;
		getNewPLLScore(u,v,edgeKey,score,improvement,&aPot);

		if(improvement<0)
		{
			continue;
		}

		double moduleWideScoreImprovement=getModuleWideScoreImprovement(u,v,moduleMembers);

		if((bestu==NULL) || (bestScoreImprovement<(improvement+moduleWideScoreImprovement)))
		{
			bestu=u;
			bestTargetScore=score;
			bestScoreImprovement=improvement;
			if(bestPot!=NULL)
			{
				delete bestPot;
			}
			bestPot=aPot;
		}
		else
		{
			delete aPot;
		}
	}

	// We could not find a parent to add to v that would improve the score.
	if((bestu==NULL) || (bestScoreImprovement<=0))
	{
		return 0;
	}

	MetaMove* move=new MetaMove;
	move->setSrcVertex(bestu->getID());
	move->setTargetVertex(v->getID());
	move->setTargetMBScore(bestTargetScore);
	move->setScoreImprovement(bestScoreImprovement);
	move->setDestPot(bestPot);
	moveSet.push_back(move);

	return 0;
}
double
MetaLearner::getModuleWideScoreImprovement(Variable* u, Variable* v,map<string,int>* moduleMembers)
{
	return 0;
	VSET& varSet=varManager->getVariableSet();
	int vID=v->getID();
	int uID=u->getID();
	double moduleScore=0;
	double neighborCnt=0;
	for(map<string,int>::iterator mIter=moduleMembers->begin();mIter!=moduleMembers->end();mIter++)
	{
		int nID=varManager->getVarID(mIter->first.c_str());
		if(nID==vID)
		{
			continue;
		}

		if(nID==uID)
		{
			continue;
		}
		SlimFactor* sFactor=factorGraph->getFactorAt(nID);
		if(sFactor->mergedMB.find(uID)!=sFactor->mergedMB.end())
		{
			continue;
		}
		double scoreImprovement;
		double mbScore;
		Potential* dPot=NULL;
		string edgeKey(u->getName());
		edgeKey.append("\t");
		edgeKey.append(mIter->first.c_str());
		getNewPLLScore(u,varSet[nID],edgeKey,mbScore,scoreImprovement,&dPot);
		if(dPot!=NULL)
		{
			delete dPot;
		}
		moduleScore=moduleScore+scoreImprovement;
		neighborCnt++;
	}
	double score=moduleScore/neighborCnt;
	return score;
}


int
MetaLearner::getNewPLLScore(Variable* u, Variable* v, string& edgeKey, double& mbScore, double& scoreImprovement, Potential** newdPot)
{
	VSET& varSet=varManager->getVariableSet();
	bool dPotDel=true;
	scoreImprovement=0;
	double currPrior=varNeighborhoodPrior[v->getID()];
	double plus=0;
	double minus=0;

	SlimFactor* dFactor=factorGraph->getFactorAt(v->getID());

	// If u is already a parent of v, then we dont need to remove it at the end of this function.
	if(dFactor->mergedMB.find(u->getID())!=dFactor->mergedMB.end())
	{
		dPotDel=false;
	}
	dFactor->mergedMB[u->getID()]=0;

	Potential *dPot=new Potential;
	dPot->setAssocVariable(varSet[dFactor->fId],Potential::FACTOR);
	for(INTINTMAP_ITER mIter=dFactor->mergedMB.begin();mIter!=dFactor->mergedMB.end();mIter++)
	{
		Variable* aVar=varSet[mIter->first];
		dPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		double eprior=getEdgePrior(mIter->first,v->getID());
		double moduleContrib=getModuleContribLogistic((string&)v->getName(),(string&)u->getName());
		double edgeProb=1/(1+exp(-1*(eprior+moduleContrib)));
		double edgeProbOld=1/(1+exp(-1*(eprior)));
		minus=minus+log(1-edgeProbOld);
		plus=plus+log(edgeProb);
	}
	dPot->potZeroInit();
	dPot->setCondBias(dFactor->potFunc->getCondBias());
	dPot->setCondVariance(dFactor->potFunc->getCondVariance());
	dPot->setCondWeight(dFactor->potFunc->getCondWeight());

	currPrior=currPrior+plus-minus;
	double newPLL_d=0;
	*newdPot=dPot;
	potManager->populatePotential(*newdPot);
	(*newdPot)->initMBCovMean();

	if((dPot->getCondVariance()<0) || (isnan(dPot->getCondVariance())) || (isinf(dPot->getCondVariance())))
	{
		scoreImprovement=-1;
	}

	if(scoreImprovement!=-1)
	{
		newPLL_d=getNewPLLScore_Tracetrick(v->getID(),u->getID(),*newdPot);
		newPLL_d=newPLL_d+currPrior;
		double oldPLL_d=(*currPLL)[v->getID()];
		double dImpr=newPLL_d-oldPLL_d;
		if(edgePresenceProb.find(edgeKey)==edgePresenceProb.end())
		{
			cout <<"No edge prior for " << edgeKey.c_str() << endl;
			exit(0);
		}
		
		if(dImpr<=0)
		{
			scoreImprovement=-1;
		}
		else
		{
			if(!dPotDel)
			{	
				cout <<"Trying to add the same regulator " << u->getName() <<" to " << v->getName() << endl;
			}

			scoreImprovement=dImpr;
			mbScore=newPLL_d;
		}
	}
	if (dPotDel)
	{
		SlimFactor* dFactor=factorGraph->getFactorAt(v->getID());
		INTINTMAP_ITER dIter=dFactor->mergedMB.find(u->getID());
		dFactor->mergedMB.erase(dIter);		
	}
	if(scoreImprovement<0)
	{
		delete dPot;
		*newdPot=NULL;
	}
	return 0;
}

double
MetaLearner::getInitPLLScore(int vId)
{
	double pll=0; 
	//Need to fix this to be set automatically
	double paramCnt=0;
	int datasize=0;
	//get parameter prior
	double paramPrior=0;

	INTINTMAP* tSet=&evidenceManager->getTrainingSet();
	datasize=datasize+tSet->size();
	for(INTINTMAP_ITER eIter=tSet->begin();eIter!=tSet->end();eIter++)
	{
		EMAP* evidMap=evidenceManager->getEvidenceAt(eIter->first);

		double cll=0;
		SlimFactor* sFactor=factorGraph->getFactorAt(vId);
		Potential* sPot=sFactor->potFunc;
		double pval=sPot->getCondPotValueFor(evidMap);
		if(isnan(pval))
		{
			cout <<"Pval is nan for datapoint " << eIter->first << endl;
		}
		if(pval<1e-50)
		{
			pval=1e-50;
		}
		cll=cll+pval;
		if(eIter==tSet->begin())
		{
			double vCnt=(double)sPot->getAssocVariables().size();
			paramCnt=paramCnt+(2*vCnt)+((vCnt*(vCnt-1))/2);
			//paramCnt=vCnt-1;
		}

		pll=pll+log(cll);
	}
	//pll=pll-(0.5*paramCnt*log(datasize));
	pll=pll-(lambda*paramCnt*log(datasize));
	return pll;
}

double
MetaLearner::getNewPLLScore_Tracetrick(int vId, int uId, Potential* newPot)
{
	double pll=0; 
	//Need to fix this to be set automatically
	//get parameter prior
	double paramPrior=0;
	Potential* parentPot=new Potential;
	VSET& vars=newPot->getAssocVariables();
	for(VSET_ITER vIter=vars.begin();vIter!=vars.end();vIter++)
	{
		if(vIter->first==vId)
		{
			continue;
		}
		Variable* aVar=vars[vIter->first];
		parentPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
	}
	parentPot->potZeroInit();
	potManager->populatePotential(parentPot);
	//parentPot->initMBCovMean();
	INTINTMAP* tSet=&evidenceManager->getTrainingSet();
	int datasize=tSet->size();
	double jointll1=newPot->computeLL_Tracetrick(datasize);
	double jointll2=parentPot->computeLL_Tracetrick(datasize);
	double vCnt=(double)newPot->getAssocVariables().size();
	double paramCnt=paramCnt+(2*vCnt)+((vCnt*(vCnt-1))/2);
	pll=jointll1-jointll2;
	pll=pll-(lambda*paramCnt*log(datasize));
	delete parentPot;
	return pll;
}

double 
MetaLearner::getEdgePrior(int tfID, int targetID)
{
	INTDBLMAP* regPriors=NULL;
	double prior=beta1;
	double fwt = 0;
	for (map<string,map<int,INTDBLMAP*>*>::iterator pItr=edgepriormap.begin(); pItr!=edgepriormap.end(); pItr++)
	{
		double eweight=0;
		double gbeta = 0;
		map<int,INTDBLMAP*>* edgeprior = pItr->second;
		if(edgeprior->find(targetID)!=edgeprior->end())
		{
			regPriors=(*edgeprior)[targetID];
			if(regPriors->find(tfID)!=regPriors->end())
			{
				eweight=(*regPriors)[tfID];
				gbeta = betamap[pItr->first];
				fwt = fwt + gbeta*eweight;
			}
		}
	}
	prior=beta1+fwt;
	//if(prior<1e-6)
	//{
	//	prior=1e-6;
	//}
	//if(prior==1)
	//{
	//	prior=1-1e-6;
	//}
	return prior;
}

int 
MetaLearner::sortMoves()
{
	for(int m=0;m<moveSet.size();m++)
	{
		for(int n=m+1;n<moveSet.size();n++)
		{
			MetaMove* m1=moveSet[m];
			MetaMove* m2=moveSet[n];
			if(m1->getScoreImprovement()<m2->getScoreImprovement())
			{
				moveSet[m]=m2;
				moveSet[n]=m1;
			}
		}
	}
	return 0;
}

int
MetaLearner::makeMoves()
{
	INTINTMAP* affectedVariables=new INTINTMAP;
	int successMove=0;
	double netScoreDelta=0;
	for(int m=0;m<moveSet.size();m++)
	{
		MetaMove* move=moveSet[m];
		int pool=0;
		if(attemptMove(move,affectedVariables)==0)
		{
			successMove++;
			netScoreDelta=netScoreDelta+move->getScoreImprovement();
		}
		else
		{
			Potential* apot=move->getDestPot();
			delete apot;
		}
	}
	delete affectedVariables;
	//cout <<"Total successful moves " << successMove << " out of total " << moveSet.size() << " with net score improvement " << netScoreDelta<< endl;
	return 0;
}

int
MetaLearner::attemptMove(MetaMove* move, INTINTMAP* affectedVars)
{
	VSET& varSet=varManager->getVariableSet();
	string edgeKey;
	Variable* u=varSet[move->getSrcVertex()];
	Variable* v=varSet[move->getTargetVertex()];
	edgeKey.append(u->getName().c_str());
	edgeKey.append("\t");
	edgeKey.append(v->getName().c_str());

	if(edgeMap.find(edgeKey)==edgeMap.end())
	{
		cout <<"Edge " << edgeKey << " not found " << endl;
		return -1;
	}

	if((affectedVars->find(move->getSrcVertex())!=affectedVars->end()) || (affectedVars->find(move->getTargetVertex())!=affectedVars->end()))
	{
		return -1;
	}
	(*affectedVars)[move->getTargetVertex()]=0;

	SlimFactor* dFactor=factorGraph->getFactorAt(move->getTargetVertex());
	if(dFactor->mergedMB.find(move->getSrcVertex())!=dFactor->mergedMB.end())
	{
		cout <<"Stop !! Trying to add the same edge " <<edgeKey << "   "<< v->getName() << endl;
	}
	dFactor->mergedMB[move->getSrcVertex()]=0;
	delete dFactor->potFunc;
	dFactor->potFunc=move->getDestPot();
	dFactor->updatePartialMeans(dFactor->potFunc->getAllPartialMeans());
	(*currPLL)[dFactor->fId]=move->getTargetMBScore();

	//Get the module and update it's indegree
	int mID=geneModuleID[v->getName()];
	map<string,int>* currIndegree=NULL;
	if(moduleIndegree.find(mID)==moduleIndegree.end())
	{
		currIndegree=new map<string,int>;
		moduleIndegree[mID]=currIndegree;
	}
	else
	{
		currIndegree=moduleIndegree[mID];
	}
	if(currIndegree->find(v->getName())==currIndegree->end())
	{
		//cout <<"Adding new regulator " << u->getName() <<" to module " << mID << endl;
		(*currIndegree)[u->getName()]=1;
	}
	else
	{	
		//cout <<"Updating regulator " << u->getName() <<" to module " << mID << endl;
		(*currIndegree)[u->getName()]=(*currIndegree)[u->getName()]+1;
	}
	if(regulatorModuleOutdegree.find(u->getName())==regulatorModuleOutdegree.end())
	{
		regulatorModuleOutdegree[u->getName()]=1;
	}
	else
	{
		regulatorModuleOutdegree[u->getName()]=regulatorModuleOutdegree[u->getName()]+1;
	}
	//cout << "Made move for " << edgeKey.c_str() << endl;
	edgeMap[edgeKey]=1;
	int curriter=0;
	if(variableStatus.find(v->getName())==variableStatus.end())
	{
		variableStatus[v->getName()]=curriter;
	}
	else
	{
		variableStatus[v->getName()]=variableStatus[v->getName()]+1;
	}
	return 0;
}

int
MetaLearner::dumpAllGraphs(int currK,int foldid,int iter)
{
	VSET& varSet=varManager->getVariableSet();
	char aFName[1024];
	sprintf(aFName,"%s/prediction_k%d.txt",foldoutDirName,currK+1);
	ofstream oFile(aFName);
	factorGraph->dumpVarMB_PairwiseFormat(oFile,varSet);
	oFile.close();
	return 0;
}

bool 
MetaLearner::checkMBSize(int u,int v, int currK)
{
	bool check=true;
	SlimFactor* dFactor=factorGraph->getFactorAt(v);
	if((dFactor->mergedMB.size()>=currK) && (dFactor->mergedMB.find(u)==dFactor->mergedMB.end()))
	{
		check=false;
	}
	return check;
}

int
MetaLearner::initPhysicalDegree()
{
	for(map<int,map<string,int>*>::iterator mIter=moduleGeneSet.begin();mIter!=moduleGeneSet.end();mIter++)
	{
		map<string,int>* indegree=NULL;
		map<string,map<string,int>*> innet;
		map<string,int>* geneSet=mIter->second;
		for(map<string,map<string,map<string,double>*>*>::iterator gpIter=priorgraphmap.begin();gpIter!=priorgraphmap.end();gpIter++)
		{
			map<string,int>enrichedTFs;
			map<string,map<string,double>*>* priorgraph = gpIter->second;
			getEnrichedTFs(enrichedTFs,geneSet,*priorgraph);
			for(map<string,int>::iterator tfIter=enrichedTFs.begin();tfIter!=enrichedTFs.end();tfIter++)
			{
				map<string,double>* motiftgts=(*priorgraph)[tfIter->first];
				map<string,int>* ttgts;
				if (innet.find(tfIter->first) == innet.end())
				{
					ttgts = new map<string,int>();
					innet[tfIter->first] = ttgts;
				}
				else
				{
					ttgts = innet[tfIter->first];
				}
				for(map<string,double>::iterator gIter=motiftgts->begin();gIter!=motiftgts->end();gIter++)
				{
					if(geneSet->find(gIter->first)==geneSet->end())
					{
						continue;
					}
					(*ttgts)[gIter->first] = 0;
				}
			}
		}
		for (map<string,map<string,int>*>::iterator tItr=innet.begin();tItr!=innet.end();tItr++)
		{
			string tf = tItr->first;
			map<string,int>* ttgts = tItr->second;
			if (ttgts->size() == 0)
			{
				continue;
			}
			if (indegree == NULL)
			{
				indegree=new map<string,int>;
			}
			(*indegree)[tf] = ttgts->size();
		}
		
		if(indegree!=NULL)
		{
			moduleIndegree[mIter->first]=indegree;
			cout <<"Physical_indegree for module " << mIter->first << endl;
			for(map<string,int>::iterator dIter=indegree->begin();dIter!=indegree->end();dIter++)
			{
				cout << dIter->first <<"\t" << dIter->second << endl;
				if(regulatorModuleOutdegree.find(dIter->first)==regulatorModuleOutdegree.end())
				{
					regulatorModuleOutdegree[dIter->first]=dIter->second;
				}
				else
				{
					regulatorModuleOutdegree[dIter->first]=regulatorModuleOutdegree[dIter->first]+dIter->second;
				}
			}
		}
	}
	return 0;
}

int
MetaLearner::getEnrichedTFs(map<string,int>& tfSet,map<string,int>* genes,map<string,map<string,double>*>& edgeSet)
{
	VSET& varSet=varManager->getVariableSet();
	int total=varSet.size();
	int k=genes->size();
	HyperGeomPval hgp;
	for(map<string,map<string,double>*>::iterator fIter=edgeSet.begin();fIter!=edgeSet.end();fIter++)
	{
		int vID=varManager->getVarID(fIter->first.c_str());
		if(vID<0)
		{
			continue;
		}
		map<string,double>* tgtSet=fIter->second;
		//int n=tgtSet->size();
		int n=0;
		int hit=0;
		for(map<string,double>::iterator gIter=tgtSet->begin();gIter!=tgtSet->end();gIter++)
		//for(map<string,int>::iterator gIter=genes->begin();gIter!=genes->end();gIter++)
		{
			
			int vID=varManager->getVarID(gIter->first.c_str());
                        if(vID<0)
                        {
                                continue;
                        }
                        n++;
			//if(tgtSet->find(gIter->first)==tgtSet->end())
                	if(genes->find(gIter->first)==genes->end())
			{
				continue;
			}
			hit++;
		}
		double enpval=hgp.getOverRepPval(k,hit,n,total-n);
		if(enpval<0.05 && hit>4)
		//if(hit>0)
		{
			tfSet[fIter->first]=hit;
		}
	}
	return 0;
}

double
MetaLearner::getModuleContribLogistic(string& tgtName, string& tfName)
{
	if((strcmp(tgtName.c_str(),"YOR334W")==0) && strcmp(tfName.c_str(),"YPR133C")==0)
	{
		cout << "Stop here " << endl;
	}

	//return 0;
	double mbeta1=beta1;
	double mbeta2=beta_motif;
	double mprior=1/(1+exp(-(1*mbeta1)));
	if(geneModuleID.find(tgtName)==geneModuleID.end())
	{
		//double improve=log(mprior)-log(1-mprior);
		//return improve;
		return 0;
	}

	int regDegree=0;
	if(regulatorModuleOutdegree.find(tfName)!=regulatorModuleOutdegree.end())
	{
		regDegree=regulatorModuleOutdegree[tfName];
	}
	int moduleID=geneModuleID[tgtName];
	if(moduleIndegree.find(moduleID)==moduleIndegree.end())
	{
		//return log(mprior);
		//double improve=log(mprior)-log(1-mprior);
		//return improve;
		return 0;
	}
	int degree=0;
	double total=0;
	map<string,int>* moddegree=moduleIndegree[moduleID];
	if(moddegree->find(tfName)==moddegree->end())
	{
		//return log(mprior);
		//double improve=log(mprior)-log(1-mprior);
		//return improve;
		//return log(mprior);
		return 0;
	}
	degree=(*moddegree)[tfName];
	double contrib=((double) degree)/((double) regDegree);
	double mprior2=mbeta2*contrib;
	double improve=log(mprior2)-log(mprior);
	//return improve;
	return mprior2;
}

//To redefine the modules we will start with the original set of modules 
//For each original module, find for every gene its pairwise similarity to every other
//gene. merge two nodes that have the greatest pairwise similarity. replace by the merged
//regulatory program. recompute similarity of all nodes to this merged node. repeat with
//finding the next most similar pair of nodes.

int
MetaLearner::redefineModules()
{
	INTINTMAP& tSet=evidenceManager->getTrainingSet();

	map<string,int> genesWithNoNeighbors;

	// Create a node for each member of each module
	for(map<int,map<string,int>*>::iterator gIter=moduleGeneSet.begin();gIter!=moduleGeneSet.end();gIter++)
	{
		map<string,int>* moduleMembers=gIter->second;
		for(map<string,int>::iterator mIter=moduleMembers->begin();mIter!=moduleMembers->end();mIter++)
		{
			int mID=varManager->getVarID(mIter->first.c_str());
			if(mID<0)
			{
				continue;
			}
			SlimFactor* mFactor=factorGraph->getFactorAt(mID);
			INTINTMAP& mbvars1=mFactor->mergedMB;
			INTDBLMAP& regWts=mFactor->potFunc->getCondWeight();

			// If a gene has no neighbors, we dont include it in the clustering algorithm.
			if(mbvars1.size()==0)
			{
				genesWithNoNeighbors[mIter->first]=0;
				continue;
			}

			// Create a node for this gene
			HierarchicalClusterNode* node = hc.getNode(mIter->first);
			if (node == nullptr)
			{
				node = new HierarchicalClusterNode;
				node->nodeName.append(mIter->first);
				hc.addNode(node);

				// Add expression data on the new node
				for(INTINTMAP_ITER eIter=tSet.begin();eIter!=tSet.end();eIter++)
				{
					EMAP* evidMap=evidenceManager->getEvidenceAt(eIter->first);
					Evidence* evid=(*evidMap)[mID];
					double v=evid->getEvidVal();
					node->expr.push_back(v);
				}
			}

			// Add weights for incoming edges onto the node
			for(INTDBLMAP_ITER bIter=regWts.begin();bIter!=regWts.end();bIter++)
			{
				node->attrib[bIter->first]=bIter->second;
			}
		}
	}

	// Perform the new clustering
	map<int,map<string,int>*> newModules;
	hc.cluster(newModules,clusterThreshold);

	// Clear out any data representing the old module assignments
	moduleGeneSet.clear();
	geneModuleID.clear();
	regulatorModuleOutdegree.clear();
	for(map<int,map<string,int>*>::iterator mIter=moduleIndegree.begin();mIter!=moduleIndegree.end();mIter++)
	{
		mIter->second->clear();
		delete mIter->second;
	}
	moduleIndegree.clear();

	char moduleFName[1024];
	sprintf(moduleFName,"%s/fold%d/modules.txt",outputDirName,currFold);
	ofstream modFile(moduleFName);

	// Read in the new module assignments
	int largestModuleID=0;
	VSET& varSet=varManager->getVariableSet();
	for(map<int,map<string,int>*>::iterator mIter=newModules.begin();mIter!=newModules.end();mIter++)
	{
		moduleGeneSet[mIter->first]=mIter->second;
		map<string,int>* geneSet=mIter->second;
		map<string,int>* indegree=new map<string,int>;
		for(map<string,int>::iterator gIter=geneSet->begin();gIter!=geneSet->end();gIter++)
		{
			modFile << gIter->first <<"\t" << mIter->first << endl;
			geneModuleID[gIter->first]=mIter->first;
			int mID=varManager->getVarID(gIter->first.c_str());
			SlimFactor* mFactor=factorGraph->getFactorAt(mID);
			INTINTMAP& mbvars1=mFactor->mergedMB;

			for(INTINTMAP_ITER nIter=mbvars1.begin();nIter!=mbvars1.end();nIter++)
			{
				// Count incoming edges to this module per regulator
				Variable* var=varSet[nIter->first];
				if(indegree->find(var->getName())==indegree->end())
				{
					(*indegree)[var->getName()]=1;
				}
				else
				{
					(*indegree)[var->getName()]=(*indegree)[var->getName()]+1;
				}
				// Count outgoing edges from regulator to any module
				if(regulatorModuleOutdegree.find(var->getName())==regulatorModuleOutdegree.end())
				{
					regulatorModuleOutdegree[var->getName()]=1;
				}	
				else
				{
					regulatorModuleOutdegree[var->getName()]=regulatorModuleOutdegree[var->getName()]+1;
				}
			}
		}
		moduleIndegree[mIter->first]=indegree;
		largestModuleID=mIter->first;
	}
	modFile.close();

	// For any genes with no neighbors, create single gene modules
	for(map<string,int>::iterator gIter=genesWithNoNeighbors.begin();gIter!=genesWithNoNeighbors.end();gIter++)
	{
		largestModuleID++;
		map<string,int>* newmodule=new map<string,int>;
		(*newmodule)[gIter->first]=0;
		moduleGeneSet[largestModuleID]=newmodule;
		geneModuleID[gIter->first]=largestModuleID;
	}
	genesWithNoNeighbors.clear();

	return 0;
}
