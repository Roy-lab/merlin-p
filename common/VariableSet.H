#ifndef _VARIABLE_SET
#define _VARIABLE_SET

#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

class Variable;

class VariableSet
{
    public:
        VariableSet(unordered_map<string, int> varNameIDMap, vector<Variable*> variables);

    	int getVarID(const string& varName);

		vector<Variable*>& getVariables();

    private:
        unordered_map<string, int> varNameIDMap;
		vector<Variable*> variables;
};

#endif