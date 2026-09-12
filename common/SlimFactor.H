#ifndef _SLIM_FACTOR
#define _SLIM_FACTOR
#include <iostream>
using namespace std;

#include "CommonTypes.H"

class Potential;

class SlimFactor
{
	public:
		SlimFactor();
		~SlimFactor();

		int fId;
		int* vIds;

		//The merged Markov blanket 
		INTINTMAP mergedMB;

		Potential* potFunc;
};
#endif
