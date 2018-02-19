#include <vector>
#include <map>
#include <string>
using namespace std;

typedef vector<int> INTVECT;
typedef vector<string> STRVECT;

typedef map<string,int> STRINTMAP;
typedef map<string,int>::iterator STRINTMAP_ITER;

typedef map<int,double> INTDBLMAP;
typedef map<int,double>::iterator INTDBLMAP_ITER;

typedef map<string,double> STRDBLMAP;
typedef map<string,double>::iterator STRDBLMAP_ITER;

typedef map<int,double> INTDBLMAP;
typedef map<int,double>::iterator INTDBLMAP_ITER;

typedef map<int,int> INTINTMAP;
typedef map<int,int>::iterator INTINTMAP_ITER;

//EMAP stores the map of variable ids and their evidence. This corresponds
//to one line in the datafile
class Evidence;
typedef vector<Evidence*> EMAP;
typedef vector<Evidence*>::iterator EMAP_ITER;
typedef vector<EMAP*> EVIDENCE_SET;

#define CONSTR_LEN 256
#define MAXFACTORSIZE_ALLOC 200
