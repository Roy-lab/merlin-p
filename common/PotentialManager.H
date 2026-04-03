/* The goal of this class is to take all the MarkovBlankets and the factors and
 *  the data and create all potentials involving valid pairings MarkovBlankets 
 *  and factors. We create all putative potentials upfront because we want to
 *  iterate over the data once and update counts of configurations in one go*/

#ifndef _POTENTIAL_MANAGER
#define _POTENTIAL_MANAGER

using namespace std;

class Potential;
class EvidenceManager;
class Matrix;

class PotentialManager
{
	public:
		PotentialManager();
		~PotentialManager();
		int init(EvidenceManager*, bool);
		int populatePotential(Potential*);

	private:
		int reset();
		int estimateCovariance(int,int);

		gsl_matrix* ludecomp;
		gsl_permutation* perm;
		Matrix* data;
		Matrix* meanMat;
		Matrix* covMat;
};
#endif
