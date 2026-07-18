#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include "Error.H"
#include "EvidenceManager.H"
#include "EvidenceSet.H"

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile(const char* inFName)
{
	ifstream inFile(inFName);
	char* buffer = NULL;
	string buffstr;
	int bufflen = 0;

	// skip the first line (gene headers)
	if (inFile.good())
	{
		getline(inFile, buffstr);
	}

	while (inFile.good())
	{
		getline(inFile, buffstr);

		if (buffstr.length() <= 0)
		{
			continue;
		}

		if (bufflen <= buffstr.length())
		{
			if (buffer != NULL)
			{
				delete[] buffer;
			}
			bufflen = buffstr.length() + 1;
			buffer = new char[bufflen];
		}
		strcpy(buffer, buffstr.c_str());

		vector<double>* evidMap = new vector<double>;
		char* tok = strtok(buffer, "\t");

		while (tok != NULL)
		{
			double varVal = atof(tok);
			if (isinf(varVal) || isnan(varVal))
			{
				cerr << "Please remove NaNs from the expression data or check the data format. Not a valid number: " << tok << endl;
				exit(-1);
			}

			evidMap->push_back(varVal);
			tok = strtok(NULL, "\t");
		}
		evidenceSet.push_back(evidMap);
	}

	inFile.close();

	cout <<"Number of samples read: " << evidenceSet.size() << endl;

	return Error::SUCCESS;
}

EvidenceSet
EvidenceManager::getTrainingSet(int foldIndex, int foldCount)
{
	return getSet(foldIndex, foldCount, false);
}

EvidenceSet
EvidenceManager::getTestSet(int foldIndex, int foldCount)
{
	return getSet(foldIndex, foldCount, true);
}

EvidenceSet
EvidenceManager::getSet(int foldIndex, int foldCount, bool isTestSet)
{
	int testSetSize = evidenceSet.size() / foldCount;

	int testStartIndex = foldIndex * testSetSize;
	int testEndIndex = (foldIndex + 1) * testSetSize;

	if (foldIndex == foldCount - 1)
	{
		testEndIndex = evidenceSet.size();
	}

	if (foldCount == 1)
	{
		testStartIndex = -1;
		testEndIndex = -1;
	}

	vector<int> indices;

	for (int i = 0; i < evidenceSet.size(); i++)
	{
		bool isTestIndex = i >= testStartIndex && i < testEndIndex;

		if (isTestIndex == isTestSet)
		{
			indices.push_back(i);
		}
	}

	return EvidenceSet(evidenceSet, indices);
}
