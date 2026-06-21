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
		int init(EvidenceManager* evMgr, bool random, vector<int>& varIDs);

		// Creates a Potential of a single variable.
		Potential* createPotential(int factorID);

		// Creates a Potential conditioned on a set of parents.
		Potential* createPotential(int factorID, vector<int>& parentIDs);

		// Computes a Gaussian log likelihood for the factor conditioned upon its existing parents and each candidate parent.
		void computeLLs(int factorID, int sampleSize, vector<int>& existingParents, vector<int>&candidateParents, unordered_map<int, double>&scores);

	private:
		void computeSingleParentLLs(int factorID, int sampleSize, vector<int>&candidateParents, unordered_map<int, double>&scores);

		vector<double> globalMeans;
		Matrix *globalCovariances;
};
#endif
