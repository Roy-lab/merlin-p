CC=g++
CFLAGS = -g
LFLAG = -lgsl -lgslcblas 

LIB_SRC = \
	MetaLearner.C \
	MetaMove.C \
	HierarchicalCluster.C \
	HierarchicalClusterNode.C \
	HyperGeomPval.C \
	Checkpoint.C \
	common/Error.C \
	common/EvidenceManager.C \
	common/EvidenceSet.C \
	common/Potential.C \
	common/SlimFactor.C \
	common/VariableManager.C \
	common/VariableSet.C \
	common/FactorGraph.C \
	common/Matrix.C \
	common/PotentialManager.C \
	common/Variable.C \
	common/DistanceManager.C

LIB_OBJ = $(LIB_SRC:.C=.o)

LIBPATH=gsl_lib/
INCLPATH1=gsl_incl/
INCLPATH2=common

merlin: Framework.C libmerlin.a
	$(CC) Framework.C libmerlin.a -I $(INCLPATH1) -I $(INCLPATH2) -L$(LIBPATH) $(LFLAG) $(CFLAGS) -o merlin

libmerlin.a: $(LIB_OBJ)
	rm -f libmerlin.a
	ar rcs libmerlin.a $(LIB_OBJ)
	rm -f $(LIB_OBJ)

%.o: %.C
	$(CC) $(CFLAGS) -I $(INCLPATH1) -I $(INCLPATH2) -c $< -o $@

clean:
	rm -f merlin libmerlin.a $(LIB_OBJ) *~
