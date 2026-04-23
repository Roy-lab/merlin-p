#include <iostream>
#include "HierarchicalClusterNode.H"

HierarchicalClusterNode::HierarchicalClusterNode()
{
	left=nullptr;
	right=nullptr;
	parent=nullptr;
	size=1;
}
