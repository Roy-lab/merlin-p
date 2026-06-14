#ifndef _HIERARCHICAL_CLUSTER
#define _HIERARCHICAL_CLUSTER

#include <map>
#include <string>
#include <vector>
#include <queue>
#include <unordered_map>

using namespace std;

class HierarchicalClusterNode;
class Matrix;

class HierarchicalCluster
{
	public:
		HierarchicalClusterNode* getNode(string nodeName);

		void addNode(HierarchicalClusterNode* node);

		void cluster(map<int, map<string, int>*>& modules, double threshold, Matrix* correlationDistances, Matrix* sharedParentDistances);

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
			bool operator() (Pair& p1, Pair& p2)
			{
				return p1.value > p2.value;
			}
		};

		vector<Pair> estimatePairwiseDist(vector<HierarchicalClusterNode*>& nodes, Matrix* correlation, Matrix* sharedParentDistances, double threshold);
		void addMergeNode(HierarchicalClusterNode* node, unordered_map<int,HierarchicalClusterNode*>& unmergedNodes, vector<Pair>& pairs, priority_queue<Pair, vector<Pair>, ComparePair>& pairQueue, double threshold);

		static HierarchicalClusterNode* createMergeNode(Pair& pair, unordered_map<int,HierarchicalClusterNode*>& unmergedNodes, int nextNodeID);
		static bool getNextPair(priority_queue<Pair, vector<Pair>, ComparePair>& pairQueue, Pair& nextPair);
		static void generateModules(unordered_map<int, HierarchicalClusterNode*>& unmergedNodes, map<int, map<string, int>*>& modules);
		static void populateMembers(map<string, int>* members, HierarchicalClusterNode* node);

		double** distValues;
		map<string, HierarchicalClusterNode*> nodeSet;
};
#endif
