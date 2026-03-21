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
#include "Distance.H"
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
HierarchicalCluster::cluster(map<int,map<string,int>*>& modules, double threshold)
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

	priority_queue<Pair*, vector<Pair*>, ComparePair> pairQueue;

	// Populate distances between the leaf nodes.
	estimatePairwiseDist(currNodeSet, pairQueue, treenodecnt);

	vector<HierarchicalClusterNode*> internalNodes;

	// Merge pairs until threshold distance is reached.
	int nextNodeID = currNodeSet.size();
	while(currNodeSet.size()>1)
	{
		Pair *nextPair = getNextPair(pairQueue, threshold);
		if (nextPair == nullptr)
		{
			break;
		}

		HierarchicalClusterNode *newNode = createMergeNode(nextPair, currNodeSet, nextNodeID);
		internalNodes.push_back(newNode);
		addMergeNode(newNode, currNodeSet, pairQueue);

		nextNodeID += 1;
	}

	// Populates modules with the clustering represented by currNodeSet
	generateModules(currNodeSet, modules);

	// Just logs the percent variance explained by the clustering.
	calculatePercentVarianceExplained(modules);

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
	for(int i=0;i<treenodecnt;i++)
	{
		delete [] distvalues[i];
	}
	delete [] distvalues;
	while(!pairQueue.empty())
	{
		Pair* pair = pairQueue.top();
		pairQueue.pop();
		delete pair;
	}

	return 0;
}

int
HierarchicalCluster::estimatePairwiseDist(map<int,HierarchicalClusterNode*>& currNodeSet, priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue, int treenodecnt)
{
	Distance d;

	// Intantiate default distances
	distvalues=new double*[treenodecnt];
	for(int i=0;i<treenodecnt;i++)
	{
		distvalues[i]=new double[treenodecnt];
		for (int j=0;j<treenodecnt;j++)
		{
			distvalues[i][j]=-1000;
		}
	}

	for(int i=0;i<currNodeSet.size();i++)
	{
		for(int j=i+1;j<currNodeSet.size();j++)
		{
			HierarchicalClusterNode* hcNode1=currNodeSet[i];
			HierarchicalClusterNode* hcNode2=currNodeSet[j];

			double ccdist=d.computeCC(hcNode1->expr,hcNode2->expr);
			ccdist=0.5*(1-ccdist);

			double sharedSign=0;
			double den1=0;
			for(map<int,double>::iterator aIter=hcNode1->attrib.begin();aIter!=hcNode1->attrib.end();aIter++)
			{
				den1=den1+fabs(aIter->second);
				if(hcNode2->attrib.find(aIter->first)!=hcNode2->attrib.end())
				{	
					if((aIter->second*hcNode2->attrib[aIter->first])>=0)
					{
						sharedSign=sharedSign+(((fabs(aIter->second)+fabs(hcNode2->attrib[aIter->first])))/2.0);
					}
				}
			}

			double den2=0;
			for(map<int,double>::iterator aIter=hcNode2->attrib.begin();aIter!=hcNode2->attrib.end();aIter++)
			{
				den2=den2+fabs(aIter->second);
			}

			double rdist=1 - (((double)sharedSign)/((double)(den1+den2-sharedSign)));
			double dist=(ccdist+rdist)/2;
			distvalues[i][j]=dist;
			distvalues[j][i]=dist;

			Pair* p=new Pair;
			p->node1=hcNode1;
			p->node2=hcNode2;
			p->value=dist;
			pairQueue.push(p);
		}
	}
	return 0;
}

HierarchicalCluster::Pair*
HierarchicalCluster::getNextPair(priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue, double threshold)
{
	if (pairQueue.empty())
	{
		return nullptr;
	}

	// If the closest pair is further than threshold, then we are done merging.
	Pair* pair = pairQueue.top();
	if (pair->value >= threshold)
	{
		return nullptr;
	}

	//Keep popping until we reach a pair whose both members have not been visited
	while(!pairQueue.empty())
	{
		HierarchicalCluster::Pair* pair = pairQueue.top();

		if (pair->node1->parent == nullptr && pair->node2->parent == nullptr)
		{
			return pair;
		}

		delete pair;
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
HierarchicalCluster::addMergeNode(HierarchicalClusterNode* node, map<int,HierarchicalClusterNode*>& currNodeSet, priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue)
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

		Pair* newPair=new Pair;
		newPair->value=dist;
		newPair->node1=nIter->second;
		newPair->node2=node;
		pairQueue.push(newPair);
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
		cout <<"Module: " << moduleCnt << "\tSize="<< moduleMembers->size() << endl;
		moduleCnt=moduleCnt+1;
	}
	return 0;
}

int
HierarchicalCluster::calculatePercentVarianceExplained(map<int,map<string,int>*>& modules)
{
	double s_total=0;
	double s_err=0;
	map<int,double> globalMean;
	for(map<int,map<string,int>*>::iterator mIter=modules.begin();mIter!=modules.end();mIter++)
	{
		map<int,double> localMean;
		map<string,int>* geneset=mIter->second;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=nodeSet[gIter->first];
			for(int i=0;i<n->expr.size();i++)
			{
				if(localMean.find(i)==localMean.end())
				{
					localMean[i]=n->expr[i];
				}	
				else	
				{
					localMean[i]=localMean[i]+n->expr[i];
				}
			}
		}
		for(map<int,double>::iterator dIter=localMean.begin();dIter!=localMean.end();dIter++)
		{
			if(globalMean.find(dIter->first)==globalMean.end())
			{
				globalMean[dIter->first]=dIter->second;
			}
			else
			{
				globalMean[dIter->first]=globalMean[dIter->first]+dIter->second;
			}
			dIter->second=dIter->second/((double)geneset->size());
		}
		double s_err_m=0;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=nodeSet[gIter->first];
			for(int i=0;i<n->expr.size();i++)
			{
				double diff=n->expr[i]-localMean[i];
				s_err_m=s_err_m+(diff*diff);
			}
		}
		s_err=s_err+s_err_m;
		localMean.clear();
	}
	for(map<int,double>::iterator dIter=globalMean.begin();dIter!=globalMean.end();dIter++)
	{
		dIter->second=dIter->second/((double)nodeSet.size());
	}
	for(map<int,map<string,int>*>::iterator mIter=modules.begin();mIter!=modules.end();mIter++)
	{
		map<string,int>* geneset=mIter->second;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=nodeSet[gIter->first];
			for(int i=0;i<n->expr.size();i++)
			{
				double diff=n->expr[i]-globalMean[i];
				s_total=s_total+(diff*diff);
			}
		}
	}
	double pcv=1.0-(s_err/s_total);
	cout <<"Percent variance explained " << pcv << endl;
	globalMean.clear();
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
