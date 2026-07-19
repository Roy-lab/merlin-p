#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Matrix.H"
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"

HierarchicalClusterNode*
HierarchicalCluster::getNode(string nodeName)
{
	if(nodeSet.find(nodeName) == nodeSet.end()) {
		return nullptr;
	}
	return nodeSet[nodeName];
}

void
HierarchicalCluster::addNode(HierarchicalClusterNode* node)
{
	nodeSet[node->nodeName] = node;
}

void
HierarchicalCluster::cluster(map<int,map<string,int>*>& modules, double threshold, Matrix* correlationDistances, Matrix* sharedParentDistances)
{
	//The total number of nodes that can be there in a hierarchical cluster is 2n-1
	int treeNodeCnt = (nodeSet.size() * 2) - 1;

	// Instantiate default distances
	distValues = new double*[treeNodeCnt];
	for (int i = 0; i < treeNodeCnt; i++) {
		distValues[i] = new double[treeNodeCnt];
		for (int j = 0; j < treeNodeCnt; j++) {
			distValues[i][j] = -1000;
		}
	}

	// Load the nodes into a vector for fast iteration.
	vector<HierarchicalClusterNode*> allNodes;
	for(auto iter = nodeSet.begin(); iter != nodeSet.end(); iter++) {
		int nextNodeIndex = allNodes.size();
		allNodes.push_back(iter->second);
		iter->second->id = nextNodeIndex;
	}

	// Populate distances between the leaf nodes.
	vector<Pair> pairs = estimatePairwiseDist(allNodes, correlationDistances, sharedParentDistances, threshold);

	priority_queue<Pair, vector<Pair>, ComparePair> pairQueue(pairs.begin(), pairs.end());

	unordered_map<int, HierarchicalClusterNode*> unmergedNodes;
	for (int i = 0; i < allNodes.size(); i++) {
		unmergedNodes[i] = allNodes[i];
	}

	vector<HierarchicalClusterNode*> internalNodes;
	int nextNodeID = unmergedNodes.size();

	// Merge pairs until threshold distance is reached.
	while(unmergedNodes.size() > 1) {

		Pair nextPair;
		if (!getNextPair(pairQueue, nextPair)) {
			break;
		}

		HierarchicalClusterNode *newNode = createMergeNode(nextPair, unmergedNodes, nextNodeID);
		internalNodes.push_back(newNode);
		addMergeNode(newNode, unmergedNodes, pairs, pairQueue, threshold);

		nextNodeID += 1;
	}

	// Populates modules with the clustering represented by unmergedNodes
	generateModules(unmergedNodes, modules);

	// Reset parent to null on all the cached leaf nodes.
	for(map<string, HierarchicalClusterNode*>::iterator aIter = nodeSet.begin(); aIter != nodeSet.end(); aIter++) {
		aIter->second->parent = nullptr;
	}

	// Clean up all the data
	for (vector<HierarchicalClusterNode*>::iterator it = internalNodes.begin(); it != internalNodes.end(); it++) {
		delete *it;
	}
	for (int i = 0; i < treeNodeCnt; i++) {
		delete [] distValues[i];
	}
	delete [] distValues;
}

vector<HierarchicalCluster::Pair>
HierarchicalCluster::estimatePairwiseDist(vector<HierarchicalClusterNode*>& nodes, Matrix* correlationDistances, Matrix* sharedParentDistances, double threshold)
{
	vector<Pair> pairs;
	for (int i = 0; i < nodes.size(); i++) {

		HierarchicalClusterNode* hcNode1 = nodes[i];

		for(int j = i + 1; j < nodes.size(); j++) {

			HierarchicalClusterNode* hcNode2 = nodes[j];

			double ccdist = correlationDistances->getValue(hcNode1->varID, hcNode2->varID);
			double rdist = sharedParentDistances->getValue(hcNode1->varID, hcNode2->varID);

			double dist = (ccdist + rdist) / 2;
			distValues[i][j] = dist;
			distValues[j][i] = dist;

			if (dist >= threshold) {
				continue;
			}

			Pair p;
			p.node1 = hcNode1;
			p.node2 = hcNode2;
			p.value = dist;
			pairs.push_back(p);
		}
	}
	return pairs;
}

bool
HierarchicalCluster::getNextPair(priority_queue<Pair, vector<Pair>, ComparePair>& pairQueue, Pair& nextPair)
{
	//Keep popping until we reach a pair whose both members have not been visited
	while(!pairQueue.empty()) {
		Pair pair = pairQueue.top();
		pairQueue.pop();

		if (pair.node1->parent == nullptr && pair.node2->parent == nullptr) {
			nextPair = pair;
			return true;
		}
	}

	// There was no unvisited node, so we are done merging.
	return false;
}

HierarchicalClusterNode*
HierarchicalCluster::createMergeNode(Pair& pair, unordered_map<int,HierarchicalClusterNode*>& unmergedNodes, int nextNodeID)
{
	HierarchicalClusterNode* c1 = pair.node1;
	HierarchicalClusterNode* c2 = pair.node2;

	// Create a new node for the merged pair
	HierarchicalClusterNode* c12 = new HierarchicalClusterNode;
	c12->id = nextNodeID;
	c12->left = c1;
	c12->right = c2;
	c12->size = c1->size + c2->size;

	c1->parent = c12;
	c2->parent = c12;

	// Remove the merged nodes from unmergedNodes
	unmergedNodes.erase(c1->id);
	unmergedNodes.erase(c2->id);

	return c12;
}

void
HierarchicalCluster::addMergeNode(HierarchicalClusterNode* node, unordered_map<int, HierarchicalClusterNode*>& unmergedNodes, vector<Pair>& pairs, priority_queue<Pair, vector<Pair>, ComparePair>& pairQueue, double threshold)
{
	HierarchicalClusterNode* c1 = node->left;
	HierarchicalClusterNode* c2 = node->right;

	double* dist_n1 = distValues[c1->id];
	double* dist_n2 = distValues[c2->id];
	double* dist_n12 = distValues[node->id];

	// For every other node, set the distance to the new node and create a new pair.
	for(auto nIter = unmergedNodes.begin(); nIter != unmergedNodes.end(); nIter++) {

		double d1 = dist_n1[nIter->first];
		double d2 = dist_n2[nIter->first];
		double dkm_rdist = ((c1->size * d1) + (c2->size * d2))/((double)(c1->size + c2->size));
		double dist = dkm_rdist;
		dist_n12[nIter->first] = dist;
		double* dist_other = distValues[nIter->first];
		dist_other[node->id] = dist;

		if (dist >= threshold) {
			continue;
		}

		Pair newPair;
		newPair.value = dist;
		newPair.node1 = nIter->second;
		newPair.node2 = node;
		pairQueue.push(newPair);
		pairs.push_back(newPair);
	}

	unmergedNodes[node->id] = node;
}

void
HierarchicalCluster::generateModules(unordered_map<int, HierarchicalClusterNode*>& unmergedNodes, map<int, map<string, int>*>& modules)
{
	int moduleCnt = modules.size();
	for(auto cIter = unmergedNodes.begin(); cIter != unmergedNodes.end(); cIter++) {
		HierarchicalClusterNode* node = cIter->second;
		map<string,int>* moduleMembers = new map<string,int>;
		modules[moduleCnt] = moduleMembers;
		populateMembers(moduleMembers, node);
		moduleCnt += 1;
	}
	cout <<"   Number of non-singleton modules: " << moduleCnt << endl;
}

void
HierarchicalCluster::populateMembers(map<string,int>* members, HierarchicalClusterNode* node)
{
	if(node->left == NULL && node->right == NULL) {
		(*members)[node->nodeName] = 0;
	} else {
		if (node->left != NULL) {
			populateMembers(members, node->left);
		}
		if (node->right != NULL) {
			populateMembers(members, node->right);
		}
	}
}
