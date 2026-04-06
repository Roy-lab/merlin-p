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

		// When init is called PotentialManager clears its internal data and loads
		// new means and covariances from the provided EvidenceManager.
		int init(EvidenceManager* evMgr, bool random);

		// Creates a Potential of a single variable.
		Potential* createPotential(int factorID);

		// Computes the conditional Gaussian log likelihood of the factor's assignments
		// conditioned upon its parents. It also creates a Potential representing this conditional
		// Gaussian PDF.
		double computeLL(int factorID, vector<int>& parentIDs, int sampleSize, Potential** newPot);

	private:
		int reset();
		double getCovariance(int,int);
		void estimateCovariance(int,int);
		double computeJointLL(Matrix*, Matrix*, double, int);

		gsl_matrix* ludecomp;
		gsl_permutation* perm;
		Matrix* data;
		Matrix* meanMat;
		Matrix* covMat;
};
#endif
