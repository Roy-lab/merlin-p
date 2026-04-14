#ifndef _HIERARCHICAL_CLUSTER_NODE
#define _HIERARCHICAL_CLUSTER_NODE

#include <map>
#include <string>
#include <vector>

using namespace std;

class HierarchicalClusterNode
{
	public:
		HierarchicalClusterNode();

		map<int,double> attrib;
		vector<double> expr;
		HierarchicalClusterNode* left;
		HierarchicalClusterNode* right;
		HierarchicalClusterNode* parent;
		int id;
		int varID;
		string nodeName;
		int size;
};
#endif
