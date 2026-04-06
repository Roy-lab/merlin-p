#ifndef _FACTOR_GRAPH
#define _FACTOR_GRAPH

#include <map>
#include <fstream>

using namespace std;

class SlimFactor;
class VariableManager;

class FactorGraph
{
	public:
		FactorGraph(VariableManager*);
		~FactorGraph();

		int getFactorCnt();
		SlimFactor* getFactorAt(int);
		int dumpVarMB(ofstream&, VSET&);

	private:
		map<int,SlimFactor*> factorSet;
};
#endif
