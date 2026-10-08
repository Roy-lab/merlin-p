#include <math.h>
#include <algorithm>
#include "DistanceManager.H"
#include "Matrix.H"
#include "EvidenceSet.H"
#include "FactorGraph.H"
#include "SlimFactor.H"
#include "Potential.H"

DistanceManager::DistanceManager()
{
    correlationDistances = nullptr;
	sharedParentDistances = nullptr;
}

Matrix*
DistanceManager::getCorrelationDistances()
{
    return correlationDistances;
}

Matrix*
DistanceManager::getSharedParentDistances()
{
    return sharedParentDistances;
}

void
DistanceManager::clearFoldData()
{
    if (correlationDistances != nullptr) {
		delete correlationDistances;
		correlationDistances = nullptr;
	}
	if (sharedParentDistances != nullptr) {
		delete sharedParentDistances;
		sharedParentDistances = nullptr;
	}
    sharedParents.clear();
    clearIterationData();
}

void
DistanceManager::clearIterationData()
{
    updatedThisIteration.clear();
}

void
DistanceManager::initDistances(EvidenceSet* trainSet)
{
    // If distances are already initialized, no need to re-initialize.
	if (correlationDistances != nullptr) {
        return;
    }

	int sampleCount = trainSet->getSampleCount();
	int varCount = trainSet->getVariableCount();

	vector<double> ssd(varCount, 0);

	for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
		for (int varIndex = 0; varIndex < varCount; varIndex++) {
			double deviation = trainSet->getEvidenceAt(varIndex, sampleIndex);
			ssd[varIndex] += deviation * deviation;
		}
	}

	correlationDistances = new Matrix(varCount, varCount);

	double threshold = sampleCount / 2.0;

	vector<double> invStd(varCount, 0);
	for (int i = 0; i < varCount; i++) {
		invStd[i] = 1.0 / sqrt(ssd[i]);
	}

	for (int i = 0; i < varCount; i++) {

		for (int j = i; j < varCount; j++) {
			double xy = 0;
			int oppRel = 0;

			for(int k = 0; k < sampleCount; k++) {
				double diff1 = trainSet->getEvidenceAt(i, k);
				double diff2 = trainSet->getEvidenceAt(j, k);
				double val = diff1 * diff2;
				xy += val;
				oppRel += (val < 0);
			}

			double cc = abs(xy) * invStd[i] * invStd[j];

			if(oppRel > threshold) {
				cc *= -1;
			}

			cc = 0.5 * (1 - cc);

			correlationDistances->setValue(cc, i, j);
			correlationDistances->setValue(cc, j, i);
		}
	}

	// Initially we have no edges, so no nodes share a parent, and all distances are 1.
	sharedParentDistances = new Matrix(varCount, varCount);
	sharedParentDistances->setAllValues(1);
}

void
DistanceManager::updateSharedParentDistances(FactorGraph* factorGraph, int varCount)
{
	vector<int> sortedTargetIDs(updatedThisIteration.begin(), updatedThisIteration.end());
	sort(sortedTargetIDs.begin(), sortedTargetIDs.end());

	vector<bool> willVisit(varCount, false);
	for (int k = 0; k < sortedTargetIDs.size(); k++) {
		willVisit[sortedTargetIDs[k]] = true;
	}

	vector<double> denoms(varCount, 0);

	for (auto iter = sortedTargetIDs.begin(); iter != sortedTargetIDs.end(); iter++) {
		int varID = *iter;
		SlimFactor* factorA = factorGraph->getFactorAt(varID);
		vector<pair<int, double>>& weightsA = factorA->potFunc->getWeights();

		double denomA = denoms[varID];
		if (denomA == 0) {
			for (const auto& weight : weightsA) {
				denomA += fabs(weight.second);
			}
			denoms[varID] = denomA;
		}

		unordered_map<int, int>& siblingVarIDs = sharedParents[varID];

		for (auto siblingIter = siblingVarIDs.begin(); siblingIter != siblingVarIDs.end(); siblingIter++) {
			int siblingID = siblingIter->first;

			if (varID == siblingID) {
				continue;
			}

			if (willVisit[siblingID] && siblingID < varID) {
				continue;
			}

			SlimFactor* factorB = factorGraph->getFactorAt(siblingID);
			vector<pair<int, double>>& weightsB = factorB->potFunc->getWeights();

			double denomB = denoms[siblingID];
			if (denomB == 0) {
				for (const auto& weight : weightsB) {
					denomB += fabs(weight.second);
				}
				denoms[siblingID] = denomB;
			}

			auto itA = weightsA.begin();
			auto itB = weightsB.begin();
			double sharedSign = 0;

			while (itA != weightsA.end() && itB != weightsB.end()) {
				if (itA->first == itB->first) {
					double weight1 = itA->second;
					double weight2 = itB->second;
					if ((weight1 >= 0.0) == (weight2 >= 0.0)) {
						sharedSign += (fabs(weight1) + fabs(weight2)) * 0.5;
					}
					++itA;
					++itB;
				} else if (itA->first < itB->first) {
					++itA;
				} else {
					++itB;
				}
			}

			double distance = 1 - sharedSign / (denomA + denomB - sharedSign);
			sharedParentDistances->setValue(distance, varID, siblingID);
			sharedParentDistances->setValue(distance, siblingID, varID);
		}
	}
}

void
DistanceManager::restoreCheckpointSharedParents(unordered_map<int, unordered_map<int, int>> edgeMap)
{
	for (auto iter = edgeMap.begin(); iter != edgeMap.end(); iter++) {
		unordered_map<int, int>& targets = iter->second;
		for (auto targetIterA = targets.begin(); targetIterA != targets.end(); targetIterA++) {
			for (auto targetIterB = targets.begin(); targetIterB != targets.end(); targetIterB++) {
				sharedParents[targetIterA->first][targetIterB->first] = 1;
			}
		}
	}
}

void 
DistanceManager::addSharedParents(unordered_map<int, unordered_map<int, int>>& edgeMap, int parentID, int childID)
{
	// Mark all other children of parentID as sharing a parent with childID.
	unordered_map<int, int> otherTargets = edgeMap[parentID];
	for (auto iter = otherTargets.begin(); iter != otherTargets.end(); iter++) {
		int siblingID = iter->first;
		sharedParents[childID][siblingID] = 1;
		sharedParents[siblingID][childID] = 1;
	}

	addUpdatedThisIteration(childID);
}

void
DistanceManager::addUpdatedThisIteration(int targetID)
{
    updatedThisIteration.insert(targetID);
}
