#ifndef _FACTOR_GRAPH
#define _FACTOR_GRAPH

#include <map>
#include <vector>
#include <fstream>

using namespace std;

class SlimFactor;
class VariableSet;
class Variable;

class FactorGraph
{
	public:
		FactorGraph(VariableSet*);
		~FactorGraph();

		int getFactorCnt();
		SlimFactor* getFactorAt(int);

	private:
		map<int,SlimFactor*> factorSet;
};
#endif
