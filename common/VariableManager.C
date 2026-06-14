#include <fstream>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

//Reads the schema of the variables

Error::ErrorCode
VariableManager::readVariables(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[400000];

	if(inFile.good())
	{
		inFile.getline(buffer,400000);

		if(strlen(buffer)<=0)
		{
			cout <<"Error: gene expression header is empty" << endl;
			return Error::VARSCHEMA_ERR;
		}

		char* tok=strtok(buffer,"\t");
		int tokCnt=0;

		while(tok!=NULL)
		{
			Variable* var=new Variable;
			var->setID(tokCnt);
			var->setName(tok);
			variableSet[tokCnt]=var;

			string varKey(tok);
			varNameIDMap[varKey]=tokCnt;
			tokCnt++;
			tok=strtok(NULL,"\t");
		}
	}

	inFile.close();

	cout <<"Number of genes read: " << variableSet.size() << endl;

	return Error::SUCCESS;
}

int
VariableManager::getVarID(const string& varKey)
{
	if(varNameIDMap.find(varKey)==varNameIDMap.end()) {
		return -1;
	}
	return varNameIDMap[varKey];
}

bool
VariableManager::isValid(int varID,int varVal)
{
	Variable* rVar=variableSet[varID];
	return rVar->isValidValue(varVal);
}

unordered_map<int, Variable*>&
VariableManager::getVariableSet()
{
	return variableSet;
}

Variable*
VariableManager::getVariableAt(int vId)
{
	if(variableSet.find(vId) == variableSet.end()) {
		cout << "Illegal variable id " << vId << endl;
		return NULL;
	}
	return variableSet[vId];
}
