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
	unordered_map<int, Variable*>& variableSet = vMgr->getVariableSet();
	for(auto vIter = variableSet.begin(); vIter != variableSet.end(); vIter++)
	{
		SlimFactor* sFactor=new SlimFactor;
		sFactor->vIds=new int[1];
		sFactor->vIds[0]=vIter->first;
		sFactor->fId=vIter->first;
		factorSet[sFactor->fId]=sFactor;
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
FactorGraph::dumpVarMB(ofstream& oFile, unordered_map<int, Variable*>& variableSet)
{
	for(map<int,SlimFactor*>::iterator aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		vector<pair<int, double>>& regWts = sFactor->potFunc->getWeights();
		for (const auto& weight : regWts) {
			oFile << variableSet[weight.first]->getName() << "\t" << variableSet[sFactor->vIds[0]]->getName() << "\t" << weight.second << endl;
		}
	}
	return 0;
}
