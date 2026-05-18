#include <fstream>
#include <iostream>
#include <cstring>
#include <getopt.h>

using namespace std;
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"

#include "FactorGraph.H"
#include "PotentialManager.H"
#include "MetaMove.H"
#include "MetaLearner.H"

#include "Framework.H"

Framework::Framework() = default;
Framework::~Framework() = default;

Error::ErrorCode
Framework::init(int argc, char** argv)
{
	bool dDefault = true;
	bool oDefault = true;
	bool kDefault = true;
	bool cvDefault = true;
	bool regDefault = true;
	bool pDefault = true;
	bool rDefault = true;
	bool moduleDefault = true;
	bool hDefault = true;

	int optret;
	opterr=1;

	while((optret=getopt(argc,argv,"o:k:d:v:l:p:r:c:h:f:q:"))!=-1)
	{
		switch(optret)
		{
			case 'd':
			{
				dDefault=false;
				Error::ErrorCode eCode = varManager.readVariables(optarg);
				if(eCode != Error::SUCCESS)
				{
					cerr << Error::getErrorString(eCode) << endl;
					return eCode;
				}

				eCode = evManager.loadEvidenceFromFile(optarg);
				if(eCode != Error::SUCCESS)
				{
					cerr << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				metaLearner.setGlobalEvidenceManager(&evManager);
				metaLearner.setVariableManager(&varManager);

				break;
			}
			case 'o': // output file - predictions
			{
				oDefault=false;
				metaLearner.setOutputDirName(optarg);
				break;
			}

			case 'k':
			{
				int aSize = atoi(optarg);
				metaLearner.setMaxFactorSize_Approx(aSize);
				kDefault=false;
				break;
			}
			case 'v':
			{
				cvDefault=false;
				cvCnt = atoi(optarg);
				if(cvCnt <= 0)
				{
					cerr << "Cross validation count should be greater than zero.\n";
					return Error::UNKNOWN;
				}
				break;
			}
			case 'l':
			{
				regDefault=false;
				metaLearner.setRestrictedList(optarg);
				break;
			}
			case 'p':
			{
				pDefault=false;
				metaLearner.setBeta1(atof(optarg));
				break;
			}
			case 'r':
			{
				rDefault=false;
				metaLearner.setBeta_Motif(atof(optarg));
				break;
			}
			case 'q':
			{
				metaLearner.setPriorGraph_All(optarg);
				break;
			}
			case 'c':
			{
				moduleDefault=false;
				metaLearner.readModuleMembership(optarg);
				break;
			}
			case 'h':
			{
				hDefault = false;
				metaLearner.setClusteringThreshold(atof(optarg));
				break;
			}
			case 'f':
			{
				metaLearner.setSpecificFold(atoi(optarg));
				break;
			}
			case '?':
			{
				cerr << "Option error " << optopt << endl;
				return Error::UNKNOWN;
			}
			default:
			{
				cerr << "Unhandled option " << optret << endl;
				return Error::UNKNOWN;
			}
		}
	}


	// Validate required arguments
	if(dDefault)
	{
		cerr << "Please specify the name of expression file. (option -d)\n";
		return Error::UNKNOWN;
	}
	if(regDefault)
	{
		cerr << "Please input a file of regulators. (option -l)\n";
		return Error::UNKNOWN;
	}
	if(oDefault)
	{
		cerr << "Please specify the name of output directory. (option -o)\n";
		return Error::UNKNOWN;
	}

	// Apply defaults for unset options
	if(moduleDefault)
	{
		cout << "Setting to default clustering" << endl;
		metaLearner.setDefaultModuleMembership();
	}
	if(cvDefault)
	{
		cvCnt=1;
	}
	if(hDefault)
	{
		metaLearner.setClusteringThreshold(0.6);
	}
	if(pDefault)
	{
		metaLearner.setBeta1(-5);
	}
	if(kDefault)
	{
		metaLearner.setMaxFactorSize_Approx(300);
	}
	if(rDefault)
	{
		metaLearner.setBeta_Motif(4);
	}

	return Error::SUCCESS;
}

int
Framework::start()
{
	metaLearner.doCrossValidation(cvCnt);
	return 0;
}


int
main(int argc, char* argv[])
{
	if(argc<2)
	{
		cout <<"MERLIN-P GRN Inference " <<  endl
			<< "-d gene_expression_file " << endl
			<< "-k maxfactorsize (default size_of_dataset)" << endl
			<< "-v cross_validation_cnt (default 1)" << endl
			<< "-l restricted_regulator_fname" << endl
			<< "-p sparsity_prior (default -5)" << endl
			<< "-r module_prior (default 4)" << endl
			<< "-q prior_config" << endl
			<< "-o outputdirectory" << endl
			<< "-c moduleassignment (default random_partitioning) "<< endl
			<< "-h hierarchical_clustering_threshold (default 0.6)"<< endl
			<< "-f specificfold_torun (default is -1)" << endl ;
		return 0;
	}
	Framework fw;
	if(fw.init(argc,argv)!=Error::SUCCESS)
	{
		return -1;
	}
	fw.start();
	return 0;
}