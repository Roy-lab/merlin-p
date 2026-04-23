LFLAG = -lgsl -lgslcblas 
SRC = Framework.C MetaLearner.C MetaMove.C common/Error.C common/EvidenceManager.C common/Potential.C common/SlimFactor.C common/VariableManager.C common/Evidence.C common/FactorGraph.C common/Matrix.C common/PotentialManager.C common/Variable.C HierarchicalCluster.C HierarchicalClusterNode.C HyperGeomPval.C
LIBPATH=gsl_lib/
INCLPATH1=gsl_incl/
INCLPATH2=common

CC=g++
CFLAGS = -g

merlin: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1) -I $(INCLPATH2)  -L$(LIBPATH) $(LFLAG) $(CFLAGS) -o merlin

clean:
	rm -f merlin *~
