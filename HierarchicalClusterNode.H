#ifndef _HIERARCHICAL_CLUSTER_NODE
#define _HIERARCHICAL_CLUSTER_NODE

#include <string>

using namespace std;

class HierarchicalClusterNode
{
	public:
		HierarchicalClusterNode();

		HierarchicalClusterNode* left;
		HierarchicalClusterNode* right;
		HierarchicalClusterNode* parent;
		int id;
		int varID;
		string nodeName;
		int size;
};
#endif
