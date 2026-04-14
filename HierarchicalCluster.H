#ifndef _HIERARCHICAL_CLUSTER
#define _HIERARCHICAL_CLUSTER

#include <map>
#include <string>
#include <vector>
#include <queue>

using namespace std;

class HierarchicalClusterNode;
class Matrix;

class HierarchicalCluster
{
	public:
		HierarchicalClusterNode* getNode(string nodeName);
		void addNode(HierarchicalClusterNode* node);
		int cluster(map<int,map<string,int>*>& modules, double threshold, Matrix* correlationDistances);

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

		vector<Pair*> estimatePairwiseDist(map<int,HierarchicalClusterNode*>& currNodeSet, Matrix* correlation, double threshold);
		void addMergeNode(HierarchicalClusterNode* node, map<int,HierarchicalClusterNode*>& currNodeSet, vector<Pair*>& pairs, priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue, double threshold);

		static HierarchicalClusterNode* createMergeNode(Pair *pair, map<int,HierarchicalClusterNode*>& currNodeSet, int nextNodeID);
		static Pair* getNextPair(priority_queue<Pair*, vector<Pair*>, ComparePair>& pairQueue);
		static int generateModules(map<int,HierarchicalClusterNode*>& currNodeSet,map<int,map<string,int>*>& modules);
		static int populateMembers(map<string,int>* members,HierarchicalClusterNode* node);

		double** distvalues;
		map<string,HierarchicalClusterNode*> nodeSet;
};
#endif
