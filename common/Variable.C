#include "Variable.H"

void
Variable::setName(const char* aStr)
{
	name.append(aStr);
}

const string&
Variable::getName()
{
	return name;
}

void
Variable::setID(int aId)
{
	vId=aId;
}

int
Variable::getID()
{
	return vId;
}
