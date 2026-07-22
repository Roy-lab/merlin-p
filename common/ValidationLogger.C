#include "ValidationLogger.H"

void
ValidationLogger::setOutDirName(const char* dirName)
{
    outDirName = dirName;
}

void
ValidationLogger::setVariableSet(VariableSet* inVarSet)
{
    variableSet = inVarSet;
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
