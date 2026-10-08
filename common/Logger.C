#include <math.h>
#include <iostream>
#include "Logger.H"
#include "ModuleManager.H"
#include "FactorGraph.H"
#include "VariableSet.H"
#include "SlimFactor.H"

void
Logger::setOutDirName(const char* dirName)
{
    outDirName = dirName;
}

void
Logger::setVariableSet(VariableSet* inVarSet)
{
    variableSet = inVarSet;
}

Logger::PredictionMetrics
Logger::computePredictionMetrics(const vector<double>& tv, const vector<double>& pv)
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

void
Logger::logModules(int foldID, ModuleManager& moduleManager, FactorGraph *factorGraph)
{
    char foldoutDirName[1024];
    sprintf(foldoutDirName, "%s/fold%d", outDirName, foldID);

	char moduleFName[1024];
	sprintf(moduleFName, "%s/modules.tsv", foldoutDirName);

	ofstream modFile(moduleFName);

	if(!modFile.is_open()) {
		cerr << "Error: cannot open module output file " << moduleFName << endl;
		return;
	}

    for (auto iter = moduleManager.begin(); iter != moduleManager.end(); iter++)
    {
        int moduleID = iter->first;
        const unordered_set<string>& geneSet = iter->second;

		for (const string& geneName : geneSet)
		{
			int geneID = variableSet->getVarID(geneName);
			if (geneID == -1 || factorGraph->getFactorAt(geneID)->mergedMB.size() == 0)
			{
				continue;
			}
			modFile << geneName << "\t" << moduleID << endl;
		}
	}
}
