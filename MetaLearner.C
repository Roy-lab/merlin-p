#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>
#include <iomanip>
#include <limits>
#include <algorithm>

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
#include "MetaLearner.H"

MetaLearner::MetaLearner()
{
	shouldLoadCheckpoint = false;
	restrictedFName[0]='\0';
	random=false;
	clusterThreshold=0.5;
	specificFold=-1;
	convThreshold=1e-3;
	factorGraph=nullptr;
	currPLL=nullptr;
	correlationDistances=nullptr;
	sharedParentDistances=nullptr;
}

void
MetaLearner::setShouldLoadCheckpoint(bool shouldLoadCheckpointVal)
{
	shouldLoadCheckpoint = shouldLoadCheckpointVal;
}

void
MetaLearner::setMaxFactorSize(int aVal)
{
	maxFactorSize=aVal;
}

void
MetaLearner::setBeta1(double aval)
{
	beta1=aval;
}

void
MetaLearner::initEdgePriorMeta_All()
{
	for(map<string, map<string, map<string, double>*>*>::iterator gIter = priorGraphMap.begin(); gIter != priorGraphMap.end(); gIter++) {
		map<string, map<string, double>*>* priorgraph = gIter->second;
		unordered_map<int, unordered_map<int, double>*>* edgeprior = new unordered_map<int, unordered_map<int, double>*>();
		edgePriorMap[gIter->first] = edgeprior;
		initEdgePriorMeta(gIter->first,*priorgraph,*edgeprior);
	}
}

int
MetaLearner::setPriorGraph_All(const char* aFName)
{
	ifstream inFile(aFName);
    if (!inFile.is_open())
    {
        std::cerr << "Error: Prior config file path incorrect or file cannot be opened: " << aFName << std::endl;
		return Error::DATAFILE_ERR;
    }

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
		betaMap[gname] = gbeta;
		map<string,map<string,double>*>* priorGraph = new map<string,map<string,double>*>();
		int status = setPriorGraph(fname.c_str(), *priorGraph);
		if(status != Error::SUCCESS)
		{
			std::cerr << "Error: failed to load prior graph for " << gname << " from " << fname << std::endl;
			delete priorGraph;
			return status;
		}
		priorGraphMap[gname] = priorGraph;
	}
	inFile.close();
	return Error::SUCCESS;
}

void
MetaLearner::setBeta_Motif(double aval)
{
	beta_motif=aval;
}

void
MetaLearner::setConvergenceThreshold(double aVal)
{
	convThreshold=aVal;
}

void
MetaLearner::setRestrictedList(const char* aFName)
{
	strcpy(restrictedFName,aFName);
	ifstream inFile(restrictedFName);
	string buffer;

	int count = 0; // counter for number of restricted regulators

	while(inFile.good())
	{
		getline(inFile,buffer);
		if(buffer.length()<=0)
		{
			continue;
		}
		restrictedVarList[buffer]=0;
		count++;
	}
	inFile.close();
	std::cout << "Number of regulators read: " << count << std::endl;
}

void
MetaLearner::setGlobalEvidenceManager(EvidenceManager* anEvMgr)
{
	evidenceManager=anEvMgr;
}

void
MetaLearner::setVariableManager(VariableManager* aPtr)
{
	varManager=aPtr;
}

void
MetaLearner::setOutputDirName(const char* dirPath)
{
	strcpy(outputDirName,dirPath);
}

void
MetaLearner::setClusteringThreshold(double aVal)
{
	clusterThreshold=aVal;
}

void
MetaLearner::setSpecificFold(int fid)
{
	specificFold=fid;
}

int
MetaLearner::setPriorGraph(const char* aFName, map<string,map<string,double>*>& priorGraph)
{
	ifstream inFile(aFName);
    if (!inFile.is_open())
    {
        std::cerr << "Error: Prior file path incorrect or file cannot be opened: " << aFName << std::endl;
		return Error::DATAFILE_ERR;
    }

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
	return Error::SUCCESS;
}

void
MetaLearner::setRandom(bool flag)
{
	random=flag;
}

void
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
}

void
MetaLearner::setDefaultModuleMembership()
{
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
	int vCnt = varSet.size();
	int moduleCnt=(int) sqrt(vCnt/2);
	if(moduleCnt>30)
	{
		moduleCnt=30;
	}
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	//Randomly partition the variables into clusterassignments
	vector<int> randIndex;
	double step=1.0/(double)vCnt;
	map<int,int> usedInit;
	for(int i=0;i<vCnt;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randIndex.push_back(rind);
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
	usedInit.clear();
}

int
MetaLearner::initEdgePriorMeta(const string& priorName, map<string,map<string,double>*>& graph, unordered_map<int, unordered_map<int, double>*>& edgePriors)
{
	cout << "Initializing prior: \"" << priorName << "\" " << endl;
	for(map<string,int>::iterator rIter=restrictedVarList.begin();rIter!=restrictedVarList.end();rIter++)
	{
		int regId=varManager->getVarID(rIter->first);
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
			unordered_map<int, double>* edgePriorGene = NULL;
			int tgtId = varManager->getVarID(vIter->first);
			if (tgtId == -1) {
				continue;
			}
			if (edgePriors.find(tgtId) == edgePriors.end()) {
				edgePriorGene = new unordered_map<int, double>;
				edgePriors[tgtId] = edgePriorGene;
			} else {
				edgePriorGene = edgePriors[tgtId];
			}
			double ewt=fabs(vIter->second);
			if (edgePriorGene->find(regId) == edgePriorGene->end()) {
				(*edgePriorGene)[regId] = ewt;
			} else {
				(*edgePriorGene)[regId] += ewt;
			}
			tfhit++;
		}
	}

	return 0;
}

void
MetaLearner::doCrossValidation(int foldCnt)
{
	gsl_rng* r = gsl_rng_alloc(gsl_rng_default);

	evidenceManager->setFoldCnt(foldCnt);
	evidenceManager->splitData(0);

	potManager = new PotentialManager;

	//The first key is for the fold number
	//For each fold we have a trained model. For each trained model we have the likelihood on
	//all the test sets, including the self test.
	int foldBegin = 0;
	int foldEnd = foldCnt;
	if(specificFold > -1) {
		foldBegin = specificFold;
		foldEnd = specificFold + 1;
	}

	for(int f = foldBegin; f < foldEnd; f++) {
		evidenceManager->splitData(f);
		if(random) {
			evidenceManager->randomizeEvidence(r, varManager);
		}

		vector<int> regIDs;
		for (map<string,int>::iterator iter = restrictedVarList.begin(); iter != restrictedVarList.end(); iter++)
		{
			int regID = varManager->getVarID(iter->first);
			regIDs.push_back(regID);
		}

		potManager->init(evidenceManager, random, regIDs);

		factorGraph = new FactorGraph(varManager);

		char outputDir[1024];
		sprintf(outputDir,"%s/fold%d", outputDirName, f);
		char foldOutputDirCmd[1024];
		sprintf(foldOutputDirCmd,"mkdir -p %s",outputDir);
		system(foldOutputDirCmd);

		// Begin identifying regulators/inferring modules for this fold
		start(f);

		getPredictionError_CrossValid(f);
		clearFoldSpecData();
	}
	gsl_rng_free(r);
}

void
MetaLearner::start(int currFold)
{
	sprintf(foldoutDirName, "%s/fold%d", outputDirName, currFold);

	int iter = 0;
	bool notConverged=true;
	bool loadedMetadata = false;
	if(shouldLoadCheckpoint)
	{
		loadedMetadata = readCheckpointMetadata(iter, notConverged);
		if (loadedMetadata)
		{
			iter++;
		}
	}

	initEdgePriorMeta_All(); // populates edgePriorMap
	initEdgeSet(); // populates edgeMap, edgePresenceProb, varNeighborhoodPrior, potentials

	double currGlobalScore = 0;
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();

	if (shouldLoadCheckpoint && loadedMetadata) {

		// overwrites moduleGeneSet, geneModuleID (set by readModuleMembership)
		cout << "Read modules..." << endl;
		readCheckpointModuleMembership();

		// populates moduleIndegree, regulatorModuleOutdegree; updates edgeMap; overwrites potentials
		cout << "Read networks..." << endl;
		populateGraphsFromFile();

		currGlobalScore = loadInitPLLScore(); 

		// populates variableStatus
		loadLastUpdate(); 

	} else {
		// populates moduleIndegree and regulatorModuleOutdegree ONLY IF some initial modules contain >=5 genes. Otherwise they begin empty 
		initPhysicalDegree();

		for (auto vIter = varSet.begin(); vIter != varSet.end(); vIter++) {
			Variable *var = vIter->second;
			variableStatus[var->getName()] = 0;
		}

		currGlobalScore = getInitPLLScore();
	}

	while(notConverged && iter<50) {
		cout << "Beginning regulator identification of iter " << iter << endl;

		int varID = 0;
		double scorePremodule = currGlobalScore;

		while(varID < varSet.size()) {

			Variable* v = varSet[varID];

			// If 5 iterations have passed without finding a score improving parent, then skip.
			int lastiter = variableStatus[v->getName()];
			if((iter - lastiter) >= 5) {
				cout <<"   Skipping gene " << v->getName() << "; no parents added in last 5 iters." << endl;
				varID++;
				continue;
			}

			MetaMove nextMove;
			if (!getNextMove(varID, nextMove)) {
				cout <<"   No move found " << v->getName() << endl;
				varID++;
				continue;
			}

			makeMove(nextMove, iter);

			currGlobalScore = getPLLScore();

			varID++;
		}
		cout << "   Finished identifying regulators with score " << currGlobalScore << endl;

		notConverged = (currGlobalScore - scorePremodule) > convThreshold;

		if(notConverged) {
			cout << "   Network not converged; score improvement of " << (currGlobalScore - scorePremodule) << ". Redefining modules." << endl;
			redefineModules(currFold);
		}

		scorePremodule = currGlobalScore;
		dumpAllGraphs(currFold);

		// Checkpointing
		writeCheckpointMetadata(iter, notConverged);
		writePLLScore();
		writeLastUpdate();
		doTar();

		iter++;
		updatedThisIteration.clear();
	}

	cout << "Final Score " << currGlobalScore << endl;
	finalScores[currFold] = currGlobalScore;
}

double
MetaLearner::getInitPLLScore()
{
	double initScore=0;
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
	//Initially we just sum up the marginal likelihoods
	currPLL=new unordered_map<int, double>;
	for(auto vIter = varSet.begin(); vIter != varSet.end(); vIter++) {
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
	for(auto dIter = currPLL->begin(); dIter != currPLL->end(); dIter++) {
		if(isnan(gScore) || isinf(gScore)) {
			cout << "Found nan/inf for variable " << dIter->first << endl;
		}
		gScore += dIter->second;
	}
	return gScore;
}

void
MetaLearner::clearFoldSpecData()
{
	hc = HierarchicalCluster();
	edgeMap.clear();
	sharedParents.clear();
	if (factorGraph != nullptr) {
		delete factorGraph;
		factorGraph = nullptr;
	}
	if (currPLL != nullptr) {
		delete currPLL;
		currPLL = nullptr;
	}
	if (correlationDistances != nullptr) {
		delete correlationDistances;
		correlationDistances = nullptr;
	}
	if (sharedParentDistances != nullptr) {
		delete sharedParentDistances;
		sharedParentDistances = nullptr;
	}
}

int
MetaLearner::initEdgeSet()
{
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
	for(auto uIter = varSet.begin(); uIter != varSet.end(); uIter++) {
		Variable* u = varSet[uIter->first];
		if((restrictedVarList.size()>0) && (restrictedVarList.find(u->getName())==restrictedVarList.end()))
		{
			continue;
		}

		for(auto vIter = varSet.begin(); vIter != varSet.end(); vIter++) {
			if(uIter->first == vIter->first)
			{
				continue;
			}
			Variable* v=varSet[vIter->first];
			if(geneModuleID.find(v->getName())==geneModuleID.end())
			{
				continue;
			}

			// This is going to be a directed graph. edgeKey looks like "reg_name\tgene_name"
			string edgeKey;
			edgeKey.append(u->getName().c_str());
			edgeKey.append("\t");
			edgeKey.append(v->getName().c_str());

			double initPrior=getEdgePrior(uIter->first,vIter->first);
			initPrior = 1/(1+exp(-1*initPrior));
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
	// cout <<"Restricted varlist size: " << restrictedVarList.size() << endl;
	int n=varSet.size();
	int r=restrictedVarList.size();
	int expEdgeCnt=r*(n-1);

	// Init the potentials
	for(int f=0;f<factorGraph->getFactorCnt();f++)
	{
		SlimFactor* sFactor=factorGraph->getFactorAt(f);
		sFactor->potFunc=potManager->createPotential(sFactor->fId);
	}

	return 0;
}

int
MetaLearner::getPredictionError_CrossValid(int foldid)
{
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
	char foldoutDirName[1024];
	sprintf(foldoutDirName,"%s/fold%d",outputDirName,foldid);
	INTINTMAP& testSet=evidenceManager->getTestSet();
	map<int,double> varPLL;
	for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
	{
		//for each gc, get the expected value of this datapoint
		EMAP* evidMap=evidenceManager->getEvidenceAt(dIter->first);

		for(auto vIter=geneModuleID.begin();vIter!=geneModuleID.end();vIter++)
		{
			int vId=varManager->getVarID(vIter->first);
			if(vId==-1)
			{
				continue;
			}
			Variable* v=varSet[vId];
			SlimFactor* sFactor=factorGraph->getFactorAt(vId);
			Potential* sPot=sFactor->potFunc;
			if(sPot==NULL)
			{
				cout <<"Found null for factor="<< sFactor->fId
					<< " variable=" <<varSet[sFactor->fId]->getName() << endl;
			}
			double pval=sPot->evaluateProbabilityDensity(evidMap);
			if(pval<1e-50)
			{
				pval=1e-50;
			}
			if(isinf(pval) || isnan(pval))
			{
				cout <<"Stop here. Found nan/inf for " << vIter->first << " dtpt "<< dIter->first << endl;
			}
			double cll=log(pval);
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
	vector<double> truevect;
	vector<double> predvect;
	for(auto vIter=geneModuleID.begin();vIter!=geneModuleID.end();vIter++)
	{
		int vId=varManager->getVarID(vIter->first);
		if(vId==-1)
		{
			continue;
		}
		//pFile <<vIter->first;
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
			double predval=sPot->getExpectation(evidMap);
			Evidence* evid=(*evidMap)[vId];
			double trueval=evid->getEvidVal();
			totalvar=totalvar+((trueval-truemean)*(trueval-truemean));
			//also called residuals
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
	}
	//oFile.close();
	//pFile.close();
	varPLL.clear();
	return 0;
}

bool
MetaLearner::getNextMove(int vID, MetaMove& outMove)
{
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
	Variable* v = varSet[vID];
	int maxNumRegs = maxFactorSize - 1;

	if(geneModuleID.find(v->getName()) == geneModuleID.end()) {
		return false;
	}

	// If v already has the max number of parents, dont test adding another.
	SlimFactor* dFactor = factorGraph->getFactorAt(vID);
	if(dFactor->mergedMB.size() >= maxNumRegs) {
		return false;
	}

	// Collect the new set of parents for v
	vector<int> parentIDs;
	for (INTINTMAP_ITER iter = dFactor->mergedMB.begin(); iter != dFactor->mergedMB.end(); iter++) {
		parentIDs.push_back(iter->first);
	}

	double existingParentPlus = 0;
	double existingParentMinus = 0;
	for (vector<int>::iterator iter = parentIDs.begin(); iter != parentIDs.end(); iter++) {
		Variable* parentVar = varSet[*iter];
		double eprior = getEdgePrior(*iter, vID);
		double moduleContrib = getModuleContribLogistic((string&)v->getName(), (string&)parentVar->getName());
		double edgeProb = 1 / (1 + exp(-1 * (eprior + moduleContrib)));
		double edgeProbOld = 1 / (1 + exp(-1 * eprior));
		existingParentMinus += log(1 - edgeProbOld);
		existingParentPlus += log(edgeProb);
	}

	double existingParentPrior = varNeighborhoodPrior[vID] + existingParentPlus - existingParentMinus;

	INTINTMAP* tSet = &evidenceManager->getTrainingSet();
	int datasize = tSet->size();

	int moduleID = geneModuleID[v->getName()];
	map<string, int>* moduleMembers = moduleGeneSet[moduleID];

	// Collect all the candidate parents, and the edge priors for each of them.
	vector<int> candidateParents;
	vector<double> candidatePriors;
	for(map<string,int>::iterator uIter = restrictedVarList.begin(); uIter != restrictedVarList.end(); uIter++) {
		int regID = varManager->getVarID(uIter->first);

		// Ensure we can find the regulator, and that it isnt the same node as the target.
		if(regID == -1 || vID == regID) {
			continue;
		}

		Variable* u = varSet[regID];

		// If the edge already exists, no need to test adding it.
		auto regEdgeIter = edgeMap.find(regID);
		if (regEdgeIter != edgeMap.end()) {
			auto targetEdgeIter = regEdgeIter->second.find(vID);
			if (targetEdgeIter != regEdgeIter->second.end() && targetEdgeIter->second == 1) {
				continue;
			}
		}

		double candidateEdgePrior = getEdgePrior(regID, vID);
		double candidateModuleContrib = getModuleContribLogistic((string&)v->getName(), (string&)u->getName());
		double candidateEdgeProb = 1 / (1 + exp(-1 * (candidateEdgePrior + candidateModuleContrib)));
		double candidateEdgeProbOld = 1 / (1 + exp(-1 * candidateEdgePrior));
		double candidatePlus = log(candidateEdgeProb);
		double candidateMinus = log(1 - candidateEdgeProbOld);
		double candidatePrior = existingParentPrior + candidatePlus - candidateMinus;

		candidateParents.push_back(regID);
		candidatePriors.push_back(candidatePrior);
	}

	// Collect the data likelihood for each candidate parent.
	unordered_map<int, double> candidateScores;
	potManager->computeLLs(vID, datasize, parentIDs, candidateParents, candidateScores);

	double bestScore = 0;
	double bestScoreImprovement = 0;
	Variable* bestRegulator = NULL;

	// Select the best scoring candidate parent.
	for (int i = 0; i < candidateParents.size(); i++) {
		int regID = candidateParents[i];
		double candidatePrior = candidatePriors[i];

		auto scoreIter = candidateScores.find(regID);
		if (scoreIter == candidateScores.end()) {
			continue;
		}

		double condLL = scoreIter->second;
		double score = condLL + candidatePrior;
		double improvement = score - (*currPLL)[vID];

		bool betterMoveExists = (bestRegulator != NULL) && (bestScoreImprovement >= improvement);

		// If there is no score improvement, cleanup aPot and continue.
		if (improvement < 0 || betterMoveExists) {
			continue;
		}

		bestRegulator = varSet[regID];
		bestScore = score;
		bestScoreImprovement = improvement;
	}

	// If we could not find a parent to add to v that would improve the score:
	if(bestRegulator == NULL || bestScoreImprovement <= 0) {
		return false;
	}

	parentIDs.push_back(bestRegulator->getID());

	Potential* potential = potManager->createPotential(vID, parentIDs);

	outMove.setSrcVertex(bestRegulator->getID());
	outMove.setTargetVertex(v->getID());
	outMove.setTargetMBScore(bestScore);
	outMove.setScoreImprovement(bestScoreImprovement);
	outMove.setDestPot(potential);

	return true;
}

double
MetaLearner::getInitPLLScore(int vId)
{
	SlimFactor* sFactor=factorGraph->getFactorAt(vId);
	Potential* sPot=sFactor->potFunc;

	double pll=0;

	INTINTMAP* tSet=&evidenceManager->getTrainingSet();
	for(INTINTMAP_ITER eIter=tSet->begin();eIter!=tSet->end();eIter++)
	{
		EMAP* evidMap=evidenceManager->getEvidenceAt(eIter->first);
		double pval=sPot->evaluateProbabilityDensity(evidMap);
		if(isnan(pval))
		{
			cout <<"Pval is nan for datapoint " << eIter->first << endl;
		}
		if(pval<1e-50)
		{
			pval=1e-50;
		}
		pll += log(pval);
	}

	// The initial graph has no edges, meaning this variable is univariate
	// gaussian, with just 2 params (mean, variance).
	return pll;
}

double
MetaLearner::getEdgePrior(int tfID, int targetID)
{
	double fwt = 0;
	for (map<string, unordered_map<int, unordered_map<int, double>*>*>::iterator pItr = edgePriorMap.begin(); pItr != edgePriorMap.end(); pItr++) {

		unordered_map<int, unordered_map<int, double>*>* edgeprior = pItr->second;
		unordered_map<int, unordered_map<int, double>*>::iterator targetIter = edgeprior->find(targetID);
		if (targetIter == edgeprior->end()) {
			continue;
		}

		unordered_map<int, double>* regPriors = targetIter->second;
		unordered_map<int, double>::iterator regIter = regPriors->find(tfID);
		if (regIter == regPriors->end()) {
			continue;
		}

		double eWeight = regIter->second;
		double gBeta = betaMap[pItr->first];
		fwt += gBeta * eWeight;
	}
	double prior = beta1 + fwt;
	return prior;
}

void
MetaLearner::makeMove(MetaMove& nextMove, int currIteration)
{
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
	Variable* u = varSet[nextMove.getSrcVertex()];
	Variable* v = varSet[nextMove.getTargetVertex()];

	SlimFactor* dFactor = factorGraph->getFactorAt(nextMove.getTargetVertex());

	// Clean up the old potential
	delete dFactor->potFunc;

	// Add the new parent and update the potential
	dFactor->mergedMB[nextMove.getSrcVertex()] = 0;
	dFactor->potFunc = nextMove.getDestPot();

	// Update the current score for this factor
	(*currPLL)[dFactor->fId] = nextMove.getTargetMBScore();

	int mID = geneModuleID[v->getName()];

	// Get or create an indegree map for this module
	unordered_map<string, int>* currIndegree = NULL;
	if(moduleIndegree.find(mID) == moduleIndegree.end()) {
		currIndegree = new unordered_map<string, int>;
		moduleIndegree[mID] = currIndegree;
	} else {
		currIndegree = moduleIndegree[mID];
	}

	// Increment the count of edges from u to the module of v
	if(currIndegree->find(u->getName()) == currIndegree->end()) {
		(*currIndegree)[u->getName()] = 1;
	} else {
		(*currIndegree)[u->getName()] += 1;
	}

	// Increment the count of edges from u
	if(regulatorModuleOutdegree.find(u->getName()) == regulatorModuleOutdegree.end()) {
		regulatorModuleOutdegree[u->getName()] = 1;
	} else {
		regulatorModuleOutdegree[u->getName()] += 1;
	}

	edgeMap[u->getID()][v->getID()] = 1;

	variableStatus[v->getName()] = currIteration;

	updatedThisIteration.push_back(v->getID());

	// Mark all other targets of u as sharing a parent with v.
	unordered_map<int, int> otherTargets = edgeMap[u->getID()];
	for (auto iter = otherTargets.begin(); iter != otherTargets.end(); iter++) {
		int targetID = iter->first;
		sharedParents[v->getID()][targetID] = 1;
		sharedParents[targetID][v->getID()] = 1;
	}
}

void
MetaLearner::dumpAllGraphs(int foldid)
{
	char aFName[1024];
	sprintf(aFName, "%s/prediction_k%d.txt", foldoutDirName, maxFactorSize);
	ofstream oFile(aFName);
	unordered_map<int, Variable*>& varSet=varManager->getVariableSet();
	factorGraph->dumpVarMB(oFile, varSet);
	oFile.close();
}

void
MetaLearner::initPhysicalDegree()
{
	for(map<int, map<string, int>*>::iterator mIter = moduleGeneSet.begin(); mIter != moduleGeneSet.end(); mIter++) {

		int moduleID = mIter->first;
		map<string, int>* moduleGenes = mIter->second;

		// Collect the prior edges from enriched TFs to genes in this module, from all prior graphs.
		unordered_map<string, vector<string>> priorEdges;

		for(map<string, map<string, map<string, double>*>*>::iterator gpIter = priorGraphMap.begin(); gpIter != priorGraphMap.end(); gpIter++) {

			map<string, map<string, double>*>* priorgraph = gpIter->second;

			map<string, int>enrichedTFs;
			getEnrichedTFs(enrichedTFs, moduleGenes, *priorgraph);

			for(map<string, int>::iterator tfIter = enrichedTFs.begin(); tfIter != enrichedTFs.end(); tfIter++) {

				map<string, double>* motiftgts = (*priorgraph)[tfIter->first];
				for(map<string, double>::iterator gIter = motiftgts->begin(); gIter != motiftgts->end(); gIter++) {
					if(moduleGenes->find(gIter->first) == moduleGenes->end()) {
						continue;
					}
					priorEdges[tfIter->first].push_back(gIter->first);
				}
			}
		}

		// Save the count of incoming edges per TF.
		unordered_map<string, int>* indegree = NULL;

		for (auto tItr = priorEdges.begin(); tItr != priorEdges.end(); tItr++) {
			string regID = tItr->first;
			int targetCount = tItr->second.size();
			if (targetCount == 0) {
				continue;
			}
			if (indegree == NULL) {
				indegree = new unordered_map<string, int>;
			}
			(*indegree)[regID] = targetCount;
		}

		if(indegree == NULL) {
			continue;
		}

		moduleIndegree[mIter->first] = indegree;
		cout << "Module " << moduleID << ": " << moduleGenes->size() << " genes, " << indegree->size() << " enriched TFs" << endl;

		// Increment the total count of outgoing edges for each TF

		for(auto dIter = indegree->begin(); dIter != indegree->end(); dIter++) {
			cout << " Enriched TF " << dIter->first << ": " << dIter->second << " target genes in this module across prior networks" <<endl;
			if(regulatorModuleOutdegree.find(dIter->first) == regulatorModuleOutdegree.end()) {
				regulatorModuleOutdegree[dIter->first] = dIter->second;
			} else {
				regulatorModuleOutdegree[dIter->first] = regulatorModuleOutdegree[dIter->first] + dIter->second;
			}
		}
	}
}

void
MetaLearner::getEnrichedTFs(map<string, int>& tfSet, map<string, int>* moduleGenes, map<string, map<string, double>*>& edgeSet)
{
	// A tf is enriched for a module, if it has at least 4 targets within that module, and
	// more module members than expected are targets of the tf.

	// We require 4 in-module targets in order for a tf to be enriched. If there aren't 4 total
	// module genes, then no tf can be enriched.
	if (moduleGenes->size() < 4) {
		return;
	}

	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
	int varCount = varSet.size();
	int moduleVarCount = moduleGenes->size();

	HyperGeomPval hgp;
	for(map<string, map<string, double>*>::iterator fIter = edgeSet.begin(); fIter != edgeSet.end(); fIter++) {
		int uID = varManager->getVarID(fIter->first);
		if(uID < 0) {
			continue;
		}

		map<string,double>* tgtSet = fIter->second;

		int targetsCount = 0;
		int targetsInModuleCount = 0;

		for(map<string, double>::iterator gIter = tgtSet->begin(); gIter != tgtSet->end(); gIter++) {
			int vID = varManager->getVarID(gIter->first);
            if(vID < 0) {
                continue;
            }

            targetsCount++;

            if(moduleGenes->find(gIter->first) != moduleGenes->end()) {
				targetsInModuleCount++;
			}
		}

		double enpval = hgp.getOverRepPval(moduleVarCount, targetsInModuleCount, targetsCount, varCount - targetsCount);
		if (enpval < 0.05 && targetsInModuleCount > 4) {
			tfSet[fIter->first] = targetsInModuleCount;
		}
	}
}

double
MetaLearner::getModuleContribLogistic(string& tgtName, string& tfName)
{
	auto moduleIter = geneModuleID.find(tgtName);
	if(moduleIter == geneModuleID.end()) {
		return 0;
	}

	int moduleID = moduleIter->second;

	auto degreeIter = moduleIndegree.find(moduleID);
	if(degreeIter == moduleIndegree.end()) {
		return 0;
	}

	unordered_map<string,int>* moddegree = degreeIter->second;

	auto regIter = moddegree->find(tfName);
	if(regIter == moddegree->end()) {
		return 0;
	}

	int degree = regIter->second;
	int regDegree = 0;

	auto outDegreeIter = regulatorModuleOutdegree.find(tfName);
	if(outDegreeIter != regulatorModuleOutdegree.end()) {
		regDegree = outDegreeIter->second;
	}

	double contrib= ((double) degree) / ((double) regDegree);
	return beta_motif * contrib;
}

//To redefine the modules we will start with the original set of modules
//For each original module, find for every gene its pairwise similarity to every other
//gene. merge two nodes that have the greatest pairwise similarity. replace by the merged
//regulatory program. recompute similarity of all nodes to this merged node. repeat with
//finding the next most similar pair of nodes.

void
MetaLearner::redefineModules(int currFold)
{
	INTINTMAP& tSet=evidenceManager->getTrainingSet();

	if (correlationDistances == nullptr) {
		initDistances();
	}

	map<string,int> genesWithNoNeighbors;

	// Create a node for each member of each module
	for(map<int, map<string, int>*>::iterator gIter = moduleGeneSet.begin(); gIter != moduleGeneSet.end(); gIter++) {
		map<string, int>* moduleMembers=gIter->second;

		for(map<string, int>::iterator mIter = moduleMembers->begin(); mIter != moduleMembers->end(); mIter++) {
			int mID = varManager->getVarID(mIter->first);
			if(mID < 0) {
				continue;
			}
			SlimFactor* mFactor=factorGraph->getFactorAt(mID);

			// If a gene has no neighbors, we dont include it in the clustering algorithm.
			INTINTMAP& mbvars1 = mFactor->mergedMB;
			if(mbvars1.size() == 0) {
				genesWithNoNeighbors[mIter->first] = 0;
				continue;
			}

			// Create a node for this gene
			HierarchicalClusterNode* node = hc.getNode(mIter->first);
			if (node == nullptr) {
				node = new HierarchicalClusterNode;
				node->nodeName.append(mIter->first);
				node->varID = mID;
				hc.addNode(node);
			}
		}
	}

	updateSharedParentDistances();

	// Perform the new clustering
	map<int,map<string,int>*> newModules;
	hc.cluster(newModules, clusterThreshold, correlationDistances, sharedParentDistances);

	// Clear out any data representing the old module assignments
	moduleGeneSet.clear();
	geneModuleID.clear();
	regulatorModuleOutdegree.clear();
	for(auto mIter=moduleIndegree.begin(); mIter!=moduleIndegree.end(); mIter++)
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
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
	for(map<int,map<string,int>*>::iterator mIter=newModules.begin();mIter!=newModules.end();mIter++)
	{
		moduleGeneSet[mIter->first]=mIter->second;
		map<string,int>* geneSet=mIter->second;
		unordered_map<string, int>* indegree = new unordered_map<string,int>;
		for(map<string,int>::iterator gIter=geneSet->begin();gIter!=geneSet->end();gIter++)
		{
			modFile << gIter->first <<"\t" << mIter->first << endl;
			geneModuleID[gIter->first]=mIter->first;
			int mID=varManager->getVarID(gIter->first);
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

	// For any genes with no neighbors, create single-gene modules
	cout << "   Number of parentless genes: " << genesWithNoNeighbors.size() << endl;
	for(map<string,int>::iterator gIter=genesWithNoNeighbors.begin();gIter!=genesWithNoNeighbors.end();gIter++)
	{
		largestModuleID++;
		map<string,int>* newmodule=new map<string,int>;
		(*newmodule)[gIter->first]=0;
		moduleGeneSet[largestModuleID]=newmodule;
		geneModuleID[gIter->first]=largestModuleID;
	}
	genesWithNoNeighbors.clear();
	cout << "   Finished redefining modules; " << moduleGeneSet.size() << " total modules" << endl;
}

void
MetaLearner::initDistances()
{
	INTINTMAP& samples = evidenceManager->getTrainingSet();
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();

	int varCount = varSet.size();
	int sampleCount = samples.size();

	vector<double> means(varCount, 0);

	for (INTINTMAP_ITER iter = samples.begin(); iter != samples.end(); iter++) {
		EMAP* evidMap = evidenceManager->getEvidenceAt(iter->first);
		for (int i = 0; i < varCount; i++) {
			Evidence* evid=(*evidMap)[i];
			means[i] += evid->getEvidVal();
		}
	}

	for (int i = 0; i < means.size(); i++) {
		means[i] /= sampleCount;
	}

	vector<double> ssd(varCount, 0);
	vector<vector<double>> deviations(varCount, vector<double>(sampleCount, 0));

	int sampleIndex = 0;
	for (INTINTMAP_ITER iter = samples.begin(); iter != samples.end(); iter++) {
		EMAP* evidMap = evidenceManager->getEvidenceAt(iter->first);
		for (int i = 0; i < varCount; i++) {
			double deviation = (*evidMap)[i]->getEvidVal() - means[i];
			deviations[i][sampleIndex] = deviation;
			ssd[i] += deviation * deviation;
		}
		sampleIndex++;
	}

	correlationDistances = new Matrix(varCount, varCount);

	double threshold = sampleCount / 2.0;

	vector<double> invStd(varCount, 0);
	for (int i = 0; i < varCount; i++) {
		invStd[i] = 1.0 / sqrt(ssd[i]);
	}

	for (int i = 0; i < varCount; i++) {
		double* dev_i = deviations[i].data();

		for (int j = i; j < varCount; j++) {
			double* dev_j = deviations[j].data();
			double xy = 0;
			int oppRel = 0;

			for(int k = 0; k < sampleCount; k++) {
				double diff1 = dev_i[k];
				double diff2 = dev_j[k];
				double val = diff1 * diff2;
				xy += val;
				oppRel += (val < 0);
			}

			double cc = abs(xy) * invStd[i] * invStd[j];

			if(oppRel > threshold) {
				cc *= -1;
			}

			cc = 0.5 * (1 - cc);

			correlationDistances->setValue(cc, i, j);
			correlationDistances->setValue(cc, j, i);
		}
	}

	// Initially we have no edges, so no nodes share a parent, and all distances are 1.
	sharedParentDistances = new Matrix(varCount, varCount);
	sharedParentDistances->setAllValues(1);
}

int
MetaLearner::writeCheckpointMetadata(int iter, bool notConvergedVal)
{
    char aFName[1024];
    sprintf(aFName, "%s/checkpoint.txt", foldoutDirName);

    ofstream outFile(aFName);
    outFile << "iter " << iter << endl;
    outFile << "notConverged " << notConvergedVal << endl;
    outFile.close();

    return 0;
}

bool
MetaLearner::readCheckpointMetadata(int& iter, bool& notConvergedVal)
{
    char aFName[1024];
    sprintf(aFName, "%s/checkpoint.txt", foldoutDirName);

    ifstream inFile(aFName);

    if (!inFile.is_open())
    {
        cerr << "Error: could not open checkpoint file: " << aFName << endl;
        return false;
    }

    string label;
    inFile >> label >> iter;
    inFile >> label >> notConvergedVal;
    inFile.close();

    return true;
}

int 
MetaLearner::writePLLScore()
{
	if (currPLL == NULL) {
		return -1;
	}

    char aFName[1024];
    sprintf(aFName, "%s/pll.txt", foldoutDirName);

    ofstream outFile(aFName);

	// Force maximum double precision to prevent truncation
    outFile << std::setprecision(std::numeric_limits<double>::max_digits10);

	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
    for (auto vIter = varSet.begin(); vIter != varSet.end(); vIter++) {
        Variable* var = varSet[vIter->first];
        if (geneModuleID.find(var->getName()) == geneModuleID.end()) {
			continue;
		}
        outFile << var->getName() << "\t" << (*currPLL)[vIter->first] << endl;
    }
    outFile.close();
    return 0;
}

double 
MetaLearner::loadInitPLLScore()
{
    char aFName[1024];
    sprintf(aFName, "%s/pll.txt", foldoutDirName);

    ifstream inFile(aFName);

	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();

	if (currPLL != NULL) {
		delete currPLL;
	}
	currPLL = new unordered_map<int, double>;

    double initScore = 0;
	char buffer[1024];

    while (inFile.good())
    {
        inFile.getline(buffer, 1023);
        if (strlen(buffer) <= 0) 
		{
			continue;
		}
        char* tok = strtok(buffer, "\t");
        int tokCnt = 0;
        Variable* v = NULL;
        double val = 0;
        while (tok != NULL)
        {
            if (tokCnt == 0) 
			{ 
				int vid = varManager->getVarID(tok); 
				v = varSet[vid]; 
			}
            else 
			{ 
				val = atof(tok); 
			}
            tok = strtok(NULL, "\t");
            tokCnt++;
        }
        (*currPLL)[v->getID()] = val;
        initScore += val;
    }
    return initScore;
}

int 
MetaLearner::writeLastUpdate()
{
    char aFName[1024];
    sprintf(aFName, "%s/lastUpdate.txt", foldoutDirName);

    ofstream outFile(aFName);

	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
    for (auto vIter = varSet.begin(); vIter != varSet.end(); vIter++)
    {
        Variable* var = varSet[vIter->first];
        if (geneModuleID.find(var->getName()) == geneModuleID.end())
		{
            continue;
		}
        outFile << var->getName() << "\t" << variableStatus[var->getName()] << endl;
    }
    outFile.close();
    return 0;
}

int 
MetaLearner::loadLastUpdate()
{
    char aFName[1024];
    sprintf(aFName, "%s/lastUpdate.txt", foldoutDirName);

    ifstream inFile(aFName);

	char buffer[1024];
    unordered_map<int, Variable*>& varSet = varManager->getVariableSet();

	while (inFile.good())
    {
        inFile.getline(buffer, 1023);
        if (strlen(buffer) <= 0) 
		{
			continue;
		}
        char* tok = strtok(buffer, "\t");
        int tokCnt = 0;
        Variable* var = NULL;
        int titer = 0;
        while (tok != NULL)
        {
            if (tokCnt == 0) 
			{ 
				int vid = varManager->getVarID(tok); 
				var = varSet[vid]; 
			}
            else 
			{ 
				titer = atoi(tok); 
			}
            tok = strtok(NULL, "\t");
            tokCnt++;
        }
        variableStatus[var->getName()] = titer;
    }
    return 0;
}

int 
MetaLearner::doTar()
{
    char tarCmd[1024];
    sprintf(tarCmd, "tar czf %s.tar.gz %s", foldoutDirName, foldoutDirName);
    system(tarCmd);
    return 0;
}

int
MetaLearner::readCheckpointModuleMembership()
{
	// Clear previous degree distributions
	moduleGeneSet.clear(); 
	geneModuleID.clear();

	char aFName[1024];
    sprintf(aFName, "%s/modules.txt", foldoutDirName);
	ifstream inFile(aFName);
	if (!inFile.is_open())
	{
		std::cerr << "Error: could not open checkpoint module file: " << aFName << std::endl;
		return -1;
	}
	
    char buffer[1024];
    int largestModuleID = -1;

    // Read clustered modules from modules.txt
    while(inFile.good())
    {
        inFile.getline(buffer,1023);
        if(strlen(buffer)<=0)
        {
            continue;
        }

        string geneName;
        int moduleID = -1;

        int tokCnt = 0;
        char* tok = strtok(buffer,"\t");
        while(tok != NULL)
        {
            if(tokCnt == 0)
            {
                geneName.append(tok);
            }
            else if(tokCnt == 1)
            {
                moduleID = atoi(tok);
            }
            tok = strtok(NULL,"\t");
            tokCnt++;
        }
        if(moduleID > largestModuleID)
        {
            largestModuleID = moduleID;
        }

        map<string,int>* geneSet = NULL;
        if(moduleGeneSet.find(moduleID) == moduleGeneSet.end())
        {
            geneSet = new map<string,int>;
            moduleGeneSet[moduleID] = geneSet;
        }
        else
        {
            geneSet = moduleGeneSet[moduleID];
        }
        (*geneSet)[geneName] = 0;
        geneModuleID[geneName] = moduleID;
    }
    inFile.close();

    // Recreate singleton modules (for parentless genes) that were not written to modules.txt
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
    int genesWithNoNeighborsCount = 0;
    for(auto vIter = varSet.begin(); vIter != varSet.end(); vIter++)
    {
        Variable* v = vIter->second;
        string geneName = v->getName();
        if(geneModuleID.find(geneName) != geneModuleID.end())
        {
            continue;
        }
        largestModuleID++;
        map<string,int>* newModule = new map<string,int>;
        (*newModule)[geneName] = 0;
        moduleGeneSet[largestModuleID] = newModule;
        geneModuleID[geneName] = largestModuleID;
        genesWithNoNeighborsCount++;
    }
    cout << "Recovered " << genesWithNoNeighborsCount << " parentless genes not present in modules.txt" << endl;
    cout << "Total modules after recovery: " << moduleGeneSet.size() << endl;
    return 0;
}

int
MetaLearner::populateGraphsFromFile()
{
	// Clear previous degree distributions
	moduleIndegree.clear();
	regulatorModuleOutdegree.clear();

	char aFName[1024];
    sprintf(aFName, "%s/prediction_k%d.txt", foldoutDirName, maxFactorSize);
	ifstream inFile(aFName);
	if (!inFile.is_open())
	{
		std::cerr << "Error: could not open checkpoint graph file: " << aFName << std::endl;
		return -1;
	}
	
	char buffer[1024];
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();

	// Parse edges and track degrees
	while(inFile.good())
	{	
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		Variable* u=NULL; // source variable
		Variable* v=NULL; // target variable
		while(tok!=NULL)
		{
			if(tokCnt==0) // source node
			{
				int vid=varManager->getVarID(tok);
				if (vid < 0)
				{
					break;
				}
				u=varSet[vid];
			}
			else if(tokCnt==1) // target node
			{
				int vid=varManager->getVarID(tok);
				if (vid < 0)
				{
					break;
				}
				v=varSet[vid];
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if (u == NULL || v == NULL)
		{
			continue;
		}

		// Track in-degree counts for the module containing the destination variable 'v'
		int mID=geneModuleID[v->getName()];
		unordered_map<string, int>* currIndegree=NULL;
		if(moduleIndegree.find(mID)==moduleIndegree.end())
		{
			currIndegree = new unordered_map<string, int>;
			moduleIndegree[mID]=currIndegree;
		}
		else
		{
			currIndegree=moduleIndegree[mID];
		}

		// Increment in-degree for source 'u' pointing into module 'mID'
		if(currIndegree->find(u->getName())==currIndegree->end())
		{
			(*currIndegree)[u->getName()]=1;
		}
		else
		{	
			(*currIndegree)[u->getName()]=(*currIndegree)[u->getName()]+1;
		}

		// Track overall out-degree count for the regulator variable 'u'
		if(regulatorModuleOutdegree.find(u->getName())==regulatorModuleOutdegree.end())
		{
			regulatorModuleOutdegree[u->getName()]=1;
		}
		else
		{
			regulatorModuleOutdegree[u->getName()]=regulatorModuleOutdegree[u->getName()]+1;
		}

		// Update the Factor Graph: Add source to target gene's Markov Blanket
		SlimFactor* sFactor=factorGraph->getFactorAt(u->getID());
		SlimFactor* dFactor=factorGraph->getFactorAt(v->getID());
		dFactor->mergedMB[sFactor->fId]=0;

		// Record the unique edge path string in the global edge map
		edgeMap[u->getID()][v->getID()] = 1;
	}

	// Re-initialize potentials: iterate through all factors to rebuild potential distributions using the restored edges
	for(int f=0;f<factorGraph->getFactorCnt();f++)
	{
		SlimFactor* sFactor=factorGraph->getFactorAt(f);

		// Safely clear old potential objects to prevent memory leaks
		if (sFactor->potFunc != NULL)
		{
			delete sFactor->potFunc;
			sFactor->potFunc = NULL;
		}

		if (sFactor->mergedMB.size() == 0)
		{
			sFactor->potFunc = potManager->createPotential(sFactor->fId);
			continue;
		}

		vector<int> parentIDs;
		for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
		{
			parentIDs.push_back(mIter->first);
		}

		Potential* restoredPot = potManager->createPotential(sFactor->fId, parentIDs);

		if (restoredPot == NULL) {
			sFactor->potFunc = potManager->createPotential(sFactor->fId);
		} else {
			sFactor->potFunc = restoredPot;
		}
	}

	inFile.close();
	return 0;
}

void
MetaLearner::updateSharedParentDistances()
{
	unordered_map<int, Variable*>& varSet = varManager->getVariableSet();
	int varCount = varSet.size();

	int updateCount = 0;

	sort(updatedThisIteration.begin(), updatedThisIteration.end());

	vector<bool> willVisit(varCount, false);
	for (int k = 0; k < updatedThisIteration.size(); k++) {
		willVisit[updatedThisIteration[k]] = true;
	}

	vector<double> denoms(varCount, 0);

	for (auto iter = updatedThisIteration.begin(); iter != updatedThisIteration.end(); iter++) {
		int varID = *iter;
		SlimFactor* factorA = factorGraph->getFactorAt(varID);
		vector<pair<int, double>>& weightsA = factorA->potFunc->getWeights();

		double denomA = denoms[varID];
		if (denomA == 0) {
			for (const auto& weight : weightsA) {
				denomA += fabs(weight.second);
			}
			denoms[varID] = denomA;
		}

		unordered_map<int, int>& siblingVarIDs = sharedParents[varID];

		for (auto siblingIter = siblingVarIDs.begin(); siblingIter != siblingVarIDs.end(); siblingIter++) {
			int siblingID = siblingIter->first;

			if (varID == siblingID) {
				continue;
			}

			if (willVisit[siblingID] && siblingID < varID) {
				continue;
			}

			updateCount += 1;

			SlimFactor* factorB = factorGraph->getFactorAt(siblingID);
			vector<pair<int, double>>& weightsB = factorB->potFunc->getWeights();

			double denomB = denoms[siblingID];
			if (denomB == 0) {
				for (const auto& weight : weightsB) {
					denomB += fabs(weight.second);
				}
				denoms[siblingID] = denomB;
			}

			auto itA = weightsA.begin();
			auto itB = weightsB.begin();
			double sharedSign = 0;

			while (itA != weightsA.end() && itB != weightsB.end()) {
				if (itA->first == itB->first) {
					double weight1 = itA->second;
					double weight2 = itB->second;
					if ((weight1 >= 0.0) == (weight2 >= 0.0)) {
						sharedSign += (fabs(weight1) + fabs(weight2)) * 0.5;
					}
					++itA;
					++itB;
				} else if (itA->first < itB->first) {
					++itA;
				} else {
					++itB;
				}
			}

			double distance = 1 - sharedSign / (denomA + denomB - sharedSign);
			sharedParentDistances->setValue(distance, varID, siblingID);
			sharedParentDistances->setValue(distance, siblingID, varID);
		}
	}
}
