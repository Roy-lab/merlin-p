#ifndef _POTENTIAL_MANAGER
#define _POTENTIAL_MANAGER

#include "PotentialSource.H"

using namespace std;

class Potential;
class EvidenceSource;
class Matrix;

class PotentialManager : public PotentialSource
{
	public:
		PotentialManager(EvidenceSource* evidenceSource);
		~PotentialManager();

		void setupForFold(int foldID, int foldCount, vector<int>& regIDs) override;

		Potential* createPotential(int factorID) override;

		Potential* createPotential(int factorID, vector<int>& parentIDs) override;

		// Computes a Gaussian log likelihood for the factor conditioned upon its existing parents and each candidate parent.
		void computeLLs(int factorID, int sampleSize, vector<int>& existingParents, vector<int>&candidateParents, unordered_map<int, double>&scores) override;

	private:
		void computeSingleParentLLs(int factorID, int sampleSize, vector<int>&candidateParents, unordered_map<int, double>&scores);

		vector<double> globalMeans;
		Matrix *globalCovariances;
		EvidenceSource* evidenceSource;
};
#endif
