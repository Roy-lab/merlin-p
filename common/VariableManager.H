#ifndef _VARIABLE_MANAGER
#define _VARIABLE_MANAGER

using namespace std;

class Error;

class VariableManager
{
	public:
		VariableSet* readVariables(const char* fileName, Error::ErrorCode& errorCode);

};
#endif
