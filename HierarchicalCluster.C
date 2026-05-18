#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>

#include "Error.H"
#include "Matrix.H"
#include "Variable.H"
#include "VariableManager.H"
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"

HierarchicalClusterNode*
HierarchicalCluster::getNode(string nodeName)
{
	if(nodeSet.find(nodeName)==nodeSet.end())
	{
		return nullptr;
	}
	return nodeSet[nodeName];
}

void
HierarchicalCluster::addNode(HierarchicalClusterNode* node)
{
	nodeSet[node->nodeName] = node;
}

int
HierarchicalCluster::cluster(map<int,map<string,int>*>& modules, double threshold, Matrix* correlationDistances)
{
	// currNodeSet holds the subset of nodes in the dendrogram that currently have no parent.
	map<int,HierarchicalClusterNode*> currNodeSet;

	// Load leaves into currNodeSet
	for(map<string,HierarchicalClusterNode*>::iterator nIter=nodeSet.begin();nIter!=nodeSet.end();nIter++)
	{
		int nextNodeIndex = currNodeSet.size();
		currNodeSet[nextNodeIndex]=nIter->second;
		nIter->second->id=nextNodeIndex;
	}

	//The total number of nodes that can be there in a hierarchical cluster is 2n-1
	int treenodecnt = (nodeSet.size()*2) - 1;

	// Instantiate default distances
	distvalues = new double*[treenodecnt];
	for (int i = 0; i < treenodecnt; i++)
	{
		distvalues[i] = new double[treenodecnt];
		for (int j = 0; j < treenodecnt; j++)
		{
			distvalues[i][j] = -1000;
		}
	}

	// Populate distances between the leaf nodes.
	vector<Pair*> pairs = estimatePairwiseDist(currNodeSet, correlationDistances, threshold);

	priority_queue<Pair*, vector<Pair*>, ComparePair> pairQueue(pairs.begin(), pairs.end());

	vector<HierarchicalClusterNode*> internalNodes;

	// Merge pairs until threshold distance is reached.
	int nextNodeID = currNodeSet.size();
	while(currNodeSet.size()>1)
	{
		Pair *nextPair = getNextPair(pairQueue);
		if (nextPair == nullptr)
		{
			break;
		}

		HierarchicalClusterNode *newNode = createMergeNode(nextPair, currNodeSet, nextNodeID);
		internalNodes.push_back(newNode);
		addMergeNode(newNode, currNodeSet, pairs, pairQueue, threshold);

		nextNodeID += 1;
	}

	// Populates modules with the clustering represented by currNodeSet
	generateModules(currNodeSet, modules);

	// Reset parent to null on all the cached leaf nodes.
	for(map<string,HierarchicalClusterNode*>::iterator aIter=nodeSet.begin();aIter!=nodeSet.end();aIter++)
	{
		aIter->second->parent = nullptr;
	}

	// Clean up all the data
	for (vector<HierarchicalClusterNode*>::iterator it = internalNodes.begin(); it != internalNodes.end(); it++)
	{
		delete *it;
	}
	for (int i = 0; i < treenodecnt; i++)
	{
		delete [] distvalues[i];
	}
	delete [] distvalues;
	for (int i = 0; i < pairs.size(); i++)
	{
		delete pairs[i];
	}

	return 0;
}

vector<HierarchicalCluster::Pair*>
HierarchicalCluster::estimatePairwiseDist(map<int,HierarchicalClusterNode*>& currNodeSet, Matrix* correlationDistances, double threshold)
{
	vector<double> denoms(currNodeSet.size(), 0);

	for (int i = 0; i < currNodeSet.size(); i++)
	{
		double denom = 0;
		HierarchicalClusterNode* node = currNodeSet[i];
		for(map<int,double>::iterator iter = node->attrib.begin(); iter != node->attrib.end(); iter++)
		{
			denom += fabs(iter->second);
		}
		denoms[i] = denom;
	}

	vector<Pair*> pairs;

	for (int i = 0; i < currNodeSet.size(); i++)
	{
		HierarchicalClusterNode* hcNode1 = currNodeSet[i];
		double den1 = denoms[i];

		for(int j = i + 1; j < currNodeSet.size(); j++)
		{
			HierarchicalClusterNode* hcNode2 = currNodeSet[j];
			double den2 = denoms[j];

			double ccdist = correlationDistances->getValue(hcNode1->varID, hcNode2->varID);

			double sharedSign = 0;
			for(map<int,double>::iterator aIter = hcNode1->attrib.begin(); aIter != hcNode1->attrib.end(); aIter++)
			{
				map<int, double>::iterator bIter = hcNode2->attrib.find(aIter->first);
				if(bIter == hcNode2->attrib.end())
				{
					continue;
				}

				double weight1 = aIter->second;
				double weight2 = bIter->second;
				if(weight1 * weight2 >= 0)
				{
					sharedSign += (fabs(weight1) + fabs(weight2)) * 0.5;
				}
			}

			double rdist = 1 - sharedSign / (den1 + den2 - sharedSign);
			double dist = (ccdist + rdist) / 2;
			distvalues[i][j] = dist;
			distvalues[j][i] = dist;

			if (dist >= threshold)
			{
				continue;
			}

			Pair* p = new Pair;
			p->node1 = hcNode1;
			p->node2 = hcNode2;
			p->value = dist;
			pairs.push_back(p);
		}
	}

	return pairs;
}

HierarchicalCluster::Pair*
HierarchicalCluster::getNextPair(priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue)
{
	//Keep popping until we reach a pair whose both members have not been visited
	while(!pairQueue.empty())
	{
		HierarchicalCluster::Pair* pair = pairQueue.top();

		if (pair->node1->parent == nullptr && pair->node2->parent == nullptr)
		{
			return pair;
		}

		pairQueue.pop();
	}

	// There was no unvisited node, so we are done merging.
	return nullptr;
}

HierarchicalClusterNode*
HierarchicalCluster::createMergeNode(Pair *pair, map<int,HierarchicalClusterNode*>& currNodeSet, int nextNodeID)
{
	HierarchicalClusterNode* c1=pair->node1;
	HierarchicalClusterNode* c2=pair->node2;

	// Create a new node for the merged pair
	HierarchicalClusterNode* c12=new HierarchicalClusterNode;
	c12->id = nextNodeID;
	c12->left=c1;
	c12->right=c2;
	c12->size=c1->size+c2->size;

	c1->parent=c12;
	c2->parent=c12;

	// Remove child nodes from the currNodeSet
	currNodeSet.erase(c1->id);
	currNodeSet.erase(c2->id);

	return c12;
}

void
HierarchicalCluster::addMergeNode(HierarchicalClusterNode* node, map<int,HierarchicalClusterNode*>& currNodeSet, vector<Pair*>& pairs, priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue, double threshold)
{
	HierarchicalClusterNode* c1=node->left;
	HierarchicalClusterNode* c2=node->right;

	double* dist_n1=distvalues[c1->id];
	double* dist_n2=distvalues[c2->id];
	double* dist_n12=distvalues[node->id];

	// For every other node, set the distance to the new node and create a new pair.
	for(map<int,HierarchicalClusterNode*>::iterator nIter=currNodeSet.begin();nIter!=currNodeSet.end();nIter++)
	{
		double d1=dist_n1[nIter->first];
		double d2=dist_n2[nIter->first];
		double dkm_rdist=((c1->size*d1) + (c2->size*d2))/((double)(c1->size + c2->size));
		double dist=dkm_rdist;
		dist_n12[nIter->first]=dist;
		double* dist_other=distvalues[nIter->first];
		dist_other[node->id]=dist;

		if (dist >= threshold)
		{
			continue;
		}

		Pair* newPair=new Pair;
		newPair->value=dist;
		newPair->node1=nIter->second;
		newPair->node2=node;
		pairQueue.push(newPair);
		pairs.push_back(newPair);
	}

	currNodeSet[node->id]=node;
}

int
HierarchicalCluster::generateModules(map<int,HierarchicalClusterNode*>& currNodeSet,map<int,map<string,int>*>& modules)
{
	int moduleCnt=modules.size();
	for(map<int,HierarchicalClusterNode*>::iterator cIter=currNodeSet.begin();cIter!=currNodeSet.end();cIter++)
	{
		HierarchicalClusterNode* node=cIter->second;
		map<string,int>* moduleMembers=new map<string,int>;
		modules[moduleCnt]=moduleMembers;
		populateMembers(moduleMembers,node);
		// cout <<"Module " << moduleCnt << " size: "<< moduleMembers->size() << endl;
		moduleCnt=moduleCnt+1;
	}
	cout <<"   Number of non-singleton modules: " << moduleCnt << endl;
	return 0;
}

int
HierarchicalCluster::populateMembers(map<string,int>* members,HierarchicalClusterNode* node)
{
	if(node->left==NULL && node->right==NULL)
	{
		(*members)[node->nodeName]=0;
	}
	else
	{
		if(node->left!=NULL)
		{
			populateMembers(members,node->left);
		}
		if(node->right!=NULL)
		{
			populateMembers(members,node->right);
		}
	}
	return 0;
}
