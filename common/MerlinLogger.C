#include <iostream>
#include <math.h>
#include "MerlinLogger.H"
#include "EvidenceSource.H"
#include "EvidenceSet.H"
#include "FactorGraph.H"
#include "Variable.H"
#include "VariableSet.H"
#include "Potential.H"
#include "SlimFactor.H"
#include "PotentialSource.H"

MerlinLogger::MerlinLogger(EvidenceSource* inEvSource, PotentialSource* inPotentialSource)
{
    evidenceSource = inEvSource;
    potentialSource = inPotentialSource;
}

/**
 * Computes and saves prediction quality metrics for a single cross-validation fold.
 *
 * For each gene/variable in geneModuleID, evaluates the trained model against the
 * held-out test set and writes one output row per gene to: <outputDirName>/fold<foldid>/prediction_stats.txt
 *
 * Columns written:
 *   Gene - gene/variable name
 *   PLL - summed pseudo-log-likelihood over all test samples (density based)
 *   RMSE - root mean squared error between predicted and true values (regression based)
 *   NormRMSE - RMSE normalised by the true-value range [min, max] (regression based)
 *   R2 - coefficient of determination (1 = perfect, 0 = mean baseline) (regression based)
 *   CC - Pearson correlation coefficient between predicted and true values (regression based)
 *
 * Returns 0 on success, -1 if the output file cannot be opened.
 */
int
MerlinLogger::logValidationError(int foldID, FactorGraph* factorGraph)
{
    EvidenceSet* testSet = evidenceSource->getEvidenceSet(EvidenceSource::SetType::TestSet);
    if (testSet->getSize() == 0) {
        return 0;
    }

    vector<Variable*>& varSet = variableSet->getVariables();

    char foldoutDirName[1024];
    sprintf(foldoutDirName, "%s/fold%d", outDirName, foldID);

    char statsFileName[1100];
    sprintf(statsFileName, "%s/prediction_stats.txt", foldoutDirName);

    ofstream statsFile(statsFileName);
    if (!statsFile.is_open())
    {
        cerr << "Error: cannot open output file " << statsFileName << endl;
        return -1;
    }

    // Header row — one line per gene/variable
    statsFile << "Gene\tPLL\tRMSE\tNormRMSE\tR2\tCC" << endl;

    // Log-likelihood
    vector<double> varPLL(varSet.size(), 0);
    for (int sampleIndex = 0; sampleIndex < testSet->getSize(); sampleIndex++)
    {
        for (int varID = 0; varID < varSet.size(); varID++)
        {
            SlimFactor* sFactor = factorGraph->getFactorAt(varID);
            Potential* sPot = sFactor->potFunc;

            double pval = potentialSource->evaluateProbabilityDensity(sPot, sampleIndex, EvidenceSource::SetType::TestSet);
            if (pval < 1e-50) {
                pval = 1e-50;
            }

            varPLL[varID] += log(pval);
        }
    }

    for (int varID = 0; varID < varSet.size(); varID++)
    {
        Variable* variable = varSet[varID];

        // Prediction metrics
        vector<double> truevect;
        vector<double> predvect;

		SlimFactor* sFactor = factorGraph->getFactorAt(varID);
		Potential* sPot = sFactor->potFunc;

		// Collect true values and predicted values for the held-out cells.
        for (int sampleIndex = 0; sampleIndex < testSet->getSize(); sampleIndex++)
        {
            vector<double>* evidence = testSet->getEvidenceAt(sampleIndex);
            double predVal = sPot->getExpectation(evidence);
            predvect.push_back(predVal);
            truevect.push_back((*evidence)[varID]);
        }

        PredictionMetrics met = computePredictionMetrics(truevect, predvect);

        // Write one row: GeneName   PLL   RMSE   NormRMSE   R^2   CC
        statsFile << variable->getName() << "\t"
			<< varPLL[varID] << "\t"
			<< met.rmse << "\t"
			<< met.normRmse << "\t"
			<< met.r2 << "\t"
			<< met.cc << endl;
    }

    statsFile.close();

    cout << "Fold " << foldID << " prediction stats written to " << statsFileName << endl;

    return 0;
}

void
MerlinLogger::logVariableMarkovBlankets(int foldID, int maxFactorSize, FactorGraph* factorGraph)
{
    char foldoutDirName[1024];
    sprintf(foldoutDirName, "%s/fold%d", outDirName, foldID);

	char aFName[1024];
	sprintf(aFName, "%s/prediction_k%d.txt", foldoutDirName, maxFactorSize);

    vector<Variable*>& varSet = variableSet->getVariables();

	ofstream oFile(aFName);

	for (int i = 0; i < factorGraph->getFactorCnt(); i++)
	{
		SlimFactor* sFactor = factorGraph->getFactorAt(i);
		vector<pair<int, double>>& regWts = sFactor->potFunc->getWeights();
		for (const auto& weight : regWts) {
			oFile << varSet[weight.first]->getName() << "\t" << varSet[sFactor->vIds[0]]->getName() << "\t" << weight.second << endl;
		}
	}

	oFile.close();
}