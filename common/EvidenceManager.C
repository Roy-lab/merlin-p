#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include "Error.H"
#include "EvidenceManager.H"
#include "EvidenceSet.H"

EvidenceManager::EvidenceManager()
{
	trainingSet = nullptr;
	testSet = nullptr;
}

EvidenceManager::~EvidenceManager()
{
	if (trainingSet != nullptr) {
		delete trainingSet;
	}
	if (testSet != nullptr) {
		delete testSet;
	}
}

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

	cout << "Number of samples read: " << evidenceSet.size() << endl;

	return Error::SUCCESS;
}

void
EvidenceManager::setupForFold(int foldIndex, int foldCount)
{
	if (trainingSet != nullptr) {
		delete trainingSet;
	}
	if (testSet != nullptr) {
		delete testSet;
	}

	trainingSet = createSet(foldIndex, foldCount, SetType::TrainingSet);
	testSet = createSet(foldIndex, foldCount, SetType::TestSet);
}

EvidenceSet*
EvidenceManager::getEvidenceSet(SetType type)
{
	if (type == SetType::TrainingSet) {
		return trainingSet;
	} else {
		return testSet;
	}
}

EvidenceSet*
EvidenceManager::createSet(int foldIndex, int foldCount, SetType type)
{
	bool isTestSet = type == SetType::TestSet;

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

	EvidenceSet* subset = new EvidenceSet(evidenceSet, indices);
	return subset;
}
