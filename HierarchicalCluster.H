#ifndef _HIERARCHICAL_CLUSTER
#define _HIERARCHICAL_CLUSTER

#include <map>
#include <string>
#include <vector>
#include <queue>

using namespace std;

class HierarchicalClusterNode;

class HierarchicalCluster
{
	public:
		HierarchicalClusterNode* getNode(string nodeName);
		void addNode(HierarchicalClusterNode* node);
		int cluster(map<int,map<string,int>*>& modules, double threshold);

	private:
		struct Pair
		{
			HierarchicalClusterNode* node1;
			HierarchicalClusterNode* node2;
			double value;
		};

		class ComparePair
		{
			public:
			bool operator() (Pair* p1, Pair* p2)
			{
				return p1->value > p2->value;
			}
		};

		int estimatePairwiseDist(map<int,HierarchicalClusterNode*>& currNodeSet, priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue, int treenodecnt);
		void addMergeNode(HierarchicalClusterNode* node, map<int,HierarchicalClusterNode*>& currNodeSet, priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue);
		int calculatePercentVarianceExplained(map<int,map<string,int>*>& modules);

		static HierarchicalClusterNode* createMergeNode(Pair *pair, map<int,HierarchicalClusterNode*>& currNodeSet, int nextNodeID);
		static Pair* getNextPair(priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue, double threshold);
		static int generateModules(map<int,HierarchicalClusterNode*>& currNodeSet,map<int,map<string,int>*>& modules);
		static int populateMembers(map<string,int>* members,HierarchicalClusterNode* node);

		double** distvalues;
		map<string,HierarchicalClusterNode*> nodeSet;
};
#endif
