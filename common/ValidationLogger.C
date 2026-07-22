#include <iostream>
#include "ValidationLogger.H"
#include "EvidenceSet.H"
#include "FactorGraph.H"
#include "Variable.H"
#include "VariableSet.H"
#include "Potential.H"
#include "SlimFactor.H"
#include "PotentialSource.H"

void
ValidationLogger::init(const char* dirName, VariableSet* inVarSet)
{
    outDirName = dirName;
    variableSet = inVarSet;
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
ValidationLogger::logValidationError(int foldID, EvidenceSet* testSet, FactorGraph* factorGraph, PotentialSource* potentialSource)
{
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

ValidationLogger::PredictionMetrics
ValidationLogger::computePredictionMetrics(const vector<double>& tv, const vector<double>& pv)
{
    int n = (int)tv.size();
    double tmean = 0, pmean = 0, maxv = -1e9, minv = 1e9;
    for(int i = 0; i < n; i++)
    {
        tmean += tv[i]; pmean += pv[i];
        if(tv[i] > maxv) maxv = tv[i];
        if(tv[i] < minv) minv = tv[i];
    }
    tmean /= n; pmean /= n;

    double ss_res = 0, ss_tot = 0, ss_xy = 0, ss_yy = 0;
    for(int i = 0; i < n; i++)
    {
        double dt = tv[i] - tmean, dp = pv[i] - pmean;
        ss_res += (tv[i] - pv[i]) * (tv[i] - pv[i]);
        ss_tot += dt * dt;
        ss_xy += dt * dp;
        ss_yy += dp * dp;
    }

    PredictionMetrics m;
    m.rmse = sqrt(ss_res / n);
    m.normRmse = (maxv > minv) ? m.rmse / (maxv - minv) : 0.0;
    m.r2 = (ss_tot > 0) ? 1.0 - (ss_res / ss_tot) : 0.0;
    m.cc = (ss_tot > 0 && ss_yy > 0) ? ss_xy / sqrt(ss_tot * ss_yy) : 0.0;
    return m;
}