#ifndef _VARIABLE
#define _VARIABLE

#include <string>

using namespace std;

class Variable
{
	public:

		// This is how the variable will be used throughout the system
		void setID(int);
		int getID();

		// This is the name assigned to the variable by the user
		void setName(const char*);
		const string& getName();

	private:
		string name;
		int vId;
};
#endif
