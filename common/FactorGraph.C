#include <math.h>
#include <cstring>
#include "Variable.H"
#include "Error.H"
#include "VariableSet.H"
#include "Potential.H"
#include "SlimFactor.H"
#include "FactorGraph.H"

FactorGraph::FactorGraph(VariableSet* varSet)
{
	vector<Variable*>& variableSet = varSet->getVariables();
	for (int varID = 0; varID < variableSet.size(); varID++) {
		SlimFactor* sFactor = new SlimFactor;
		sFactor->vIds = new int[1];
		sFactor->vIds[0] = varID;
		sFactor->fId = varID;
		factorSet[sFactor->fId] = sFactor;
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
