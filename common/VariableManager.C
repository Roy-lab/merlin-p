#include <fstream>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include "Error.H"
#include "Variable.H"
#include "VariableSet.H"
#include "VariableManager.H"

//Reads the schema of the variables

VariableSet*
VariableManager::readVariables(const char* aFName, Error::ErrorCode& errorCode)
{
	ifstream inFile(aFName);
	char buffer[400000];

	vector<Variable*> variableSet;
	unordered_map<string, int> varNameIDMap;

	if (inFile.good())
	{
		inFile.getline(buffer,400000);

		if (strlen(buffer) <= 0)
		{
			cout << "Error: gene expression header is empty" << endl;
			errorCode = Error::VARSCHEMA_ERR;
			return nullptr;
		}

		char* tok = strtok(buffer, "\t");
		int tokCnt = 0;

		while (tok != NULL)
		{
			Variable* var = new Variable;
			var->setID(tokCnt);
			var->setName(tok);
			variableSet.push_back(var);

			string varKey(tok);
			varNameIDMap[varKey] = tokCnt;

			tokCnt++;
			tok = strtok(NULL, "\t");
		}
	}

	inFile.close();

	cout << "Number of genes read: " << variableSet.size() << endl;

	errorCode = Error::SUCCESS;

	VariableSet* varSet = new VariableSet(varNameIDMap, variableSet);
	return varSet;
}
