#include <math.h>
#include <cstring>
#include "Variable.H"
#include "Error.H"
#include "VariableManager.H"
#include "Potential.H"
#include "SlimFactor.H"
#include "FactorGraph.H"

FactorGraph::FactorGraph(VariableManager* vMgr)
{
	map<int,Variable*>& variableSet=vMgr->getVariableSet();
	int vCnt=variableSet.size();
	// cout << " Number of factors " << vCnt << endl;

	int factorIndex=0;
	for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		SlimFactor* sFactor=new SlimFactor;
		sFactor->vIds=new int[1];
		sFactor->vIds[0]=vIter->first;
		sFactor->vCnt=1;
		sFactor->fId=factorIndex;
		factorSet[sFactor->fId]=sFactor;
		factorIndex++;
	}
}

FactorGraph::~FactorGraph()
{
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		delete fIter->second;
	}
}

int
FactorGraph::getFactorCnt()
{
	return factorSet.size();
}

SlimFactor*
FactorGraph::getFactorAt(int fid)
{
	if(factorSet.find(fid)==factorSet.end())
	{
		return NULL;
	}
	return factorSet[fid];
}

int
FactorGraph::dumpVarMB(ofstream& oFile,VSET& variableSet)
{
	for(map<int,SlimFactor*>::iterator aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		INTDBLMAP& regWts=sFactor->potFunc->getWeights();
		for(INTDBLMAP_ITER mIter=regWts.begin();mIter!=regWts.end();mIter++)
		{
			oFile << variableSet[mIter->first]->getName() << "\t" << variableSet[sFactor->vIds[0]]->getName() << "\t" << mIter->second << endl;
		}
	}
	return 0;
}
