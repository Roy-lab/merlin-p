#include "SlimFactor.H"

SlimFactor::SlimFactor()
{
	vIds = NULL;
	fId = -1;
	potFunc = NULL;
}

SlimFactor::~SlimFactor()
{
	delete[] vIds;
}
