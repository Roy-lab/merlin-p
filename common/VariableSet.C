#include "VariableSet.H"

VariableSet::VariableSet(unordered_map<string, int> nameMap, vector<Variable*> vars)
{
    varNameIDMap = nameMap;
    variables = vars;
}

int
VariableSet::getVarID(const string& varKey)
{
	if (varNameIDMap.find(varKey) == varNameIDMap.end()) {
		return -1;
	}
	return varNameIDMap[varKey];
}

vector<Variable*>&
VariableSet::getVariables()
{
	return variables;
}
