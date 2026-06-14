#ifndef _VARIABLE_MANAGER
#define _VARIABLE_MANAGER

#include <string>
#include <unordered_map>

using namespace std;

class Error;
class Variable;
class VariableManager
{
	public:
		Error::ErrorCode readVariables(const char*);
		int getVarID(const string& varName);
		bool isValid(int,int);
		unordered_map<int, Variable*>& getVariableSet();
		Variable* getVariableAt(int);
	private:
		unordered_map<string, int> varNameIDMap;
		unordered_map<int, Variable*> variableSet;
		int globalVarId;
};
#endif
