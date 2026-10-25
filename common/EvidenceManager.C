#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include <algorithm>
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

void
EvidenceManager::setExpressionFile(const char* inFName)
{
	expressionFile = inFName;
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

	loadData(foldIndex, foldCount);

	centerData(trainData, trainSampleCount);
	centerData(testData, testSampleCount);

	trainingSet = new EvidenceSet(trainData, trainSampleCount);
	testSet = new EvidenceSet(testData, testSampleCount);
}

void
EvidenceManager::loadData(int foldIndex, int foldCount)
{
	// Clear the previously loaded data
	testData.clear();
	trainData.clear();

	ifstream inFile(expressionFile);
	string line;

	// The first line should contain gene headers
	getline(inFile, line);

	size_t varCount = std::count(line.begin(), line.end(), '\t') + 1;

	// Scan the file once just to count the samples, so that we can prepare properly sized vectors
	size_t sampleCount = 0;
	while (getline(inFile, line)) {
		if (line.empty()) {
			continue;
		}
		sampleCount += 1;
	}

	cout << "Number of samples read: " << sampleCount << endl;

	// Seek back to the start and skip the header again
	inFile.clear();
	inFile.seekg(0);
	getline(inFile, line);

	// Determine the range of the test set samples
	int testStart;
	int testEnd;
	if (foldCount == 1) {
		testStart = -1;
		testEnd = -1;
	} else {
		int foldSize = sampleCount / foldCount;
		testStart = foldIndex * foldSize;
		if (foldIndex == foldCount - 1) {
			testEnd = sampleCount;
		} else {
			testEnd = (foldIndex + 1) * foldSize;
		}
	}

	testSampleCount = testEnd - testStart;
	trainSampleCount = sampleCount - testSampleCount;

	trainData.resize(varCount * trainSampleCount);
	testData.resize(varCount * testSampleCount);

	int sampleIndex = 0;
	int trainSampleIndex = 0;
	int testSampleIndex = 0;

	while (getline(inFile, line)) {
		if (line.empty()) {
			continue;
		}

		bool isTestIndex = sampleIndex >= testStart && sampleIndex < testEnd;

		int varIndex = 0;

		size_t start = 0;
		while (start < line.size()) {
			size_t end = line.find('\t', start);
			string token = line.substr(start, end - start);
			double varVal = stod(token);

			if (isinf(varVal) || isnan(varVal)) {
				cerr << "Please remove NaNs from the expression data or check the data format. Not a valid number: " << token << endl;
				exit(-1);
			}

			if (isTestIndex) {
				testData[varIndex * testSampleCount + testSampleIndex] = varVal;
			} else {
				trainData[varIndex * trainSampleCount + trainSampleIndex] = varVal;
			}

			if (end == string::npos) {
				break;
			}

			start = end + 1;
			varIndex += 1;
		}

		sampleIndex += 1;

		if (isTestIndex) {
			testSampleIndex += 1;
		} else {
			trainSampleIndex += 1;
		}
	}

	inFile.close();
}

void
EvidenceManager::centerData(vector<double>& data, size_t sampleCount)
{
	if (sampleCount == 0) {
		return;
	}

	size_t varCount = data.size() / sampleCount;

	for (size_t i = 0; i < varCount; i++) {
		double sampleSum = 0;
		for(int j = 0; j < sampleCount; j++) {
			sampleSum += data[i * sampleCount + j];
		}
		double mean = sampleSum / sampleCount;
		for (size_t j = 0; j < sampleCount; j++) {
			data[i * sampleCount + j] -= mean;
		}
	}
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
