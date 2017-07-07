What is MERLIN?
---------------
Modular regulatory network learning with per gene information (MERLIN) is a network inference method that tries to infer a more accurate regulatory network by incorporating a modularity constraint. 

For more information on the method, please refer to the original publication:
> **Integrated Module and Gene-Specific Regulatory Inference Implicates Upstream Signaling Networks**
>
> Roy et al. PLoS Comput Biol, 2013.
>
> http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003252

What is MERLIN+P?
-----------------
MERLIN+Prior (MERLIN+P) is an extension to original MERLIN method that tries to improve the accuracy of inferred network, by incorporating additional prior information.
> **A prior-based integrative framework for functional transcriptional regulatory network inference**
>
> Siahpirani and Roy, Nucl Acids Res, 2016.
>
> http://nar.oxfordjournals.org/content/early/2016/10/28/nar.gkw963.full

What are the inputs of MERLIN?
------------------------------

Required inputs
---------------
The required inputs of MERLIN are a tab delimited expression file, a list of regulators and an output directory.

**Expression file:** This is a tab delimited file containing the expression values of genes of interest in different samples/conditions/time points. Each column correspond to one gene, first row contain the gene name and the rest of rows are the expression values.
For an example, see example/net1_expression.txt

**Regulator list:** This file contains a list of potential regulators (transcription factors, signaling proteins, etc.). The file should contain one regulator per row. The program will infer a network between these regulators and all the genes in the expression file.
For an example, see example/net1_transcription_factors.tsv

**Output directory:** When running MERLIN in a cross validation scheme, this directory will contain one sub-directory per fold, otherwise, there would be only one sub-directory named fold0. Each sub-directory will contain two files:

*inferred network:* The file prediction_k*.txt will contain the inferred network, one edge per line. The first column would be the regulator's name, the second column would be the target gene's name, and the third column would be the regression coefficient.

module assignment: The file modules.txt would contain the final module assignments, one row per gene. The first column contains the gene's name, and the second column would contain the cluster assignment of that gene.

** *Note:* ** If the output directory doesn't exist, the program will run without producing any errors, but will not save the outputs.

Optional inputs
---------------

**Initial cluster assignment:** The file contains the initial cluster assignment for each gene, one gene per row. The first column contains the gene's name and the second column contains the cluster assignment for that gene. For example, see example/net1_net.txt

**Prior networks:** The config file should contains one row per prior network. The first column contains the network's name, the second column contains the path to the prior network, and the third column contains the confidence in that network (the higher the confidence, the more effect that prior will have). The prior file in the second column, should have one edge per row, the first column contains the regulator's name, the second column contains the target gene's name, and the third column contains the confidence in that edge (the higher the confidence, the more we trust in that edge). For an example, see example/net1_config.txt and example/net1_net.txt

How to run
----------
To run, change to MERLIN directory and run:

```
./merlin -d example/net1_expression.txt -l example/net1_transcription_factors.tsv -o out_dir/ -c example/clusterassign.txt -q example/net1_config.txt
```

To get the full list of options, run: ./merlin

How to compile?
---------------

To compile the code, changed to MERLIN directory and run: make

Inferred networks and other resources
-------------------------------------
The expression datasets used in this study, prior networks and the inferred networks were moved to this repository:
https://github.com/Roy-lab/merlin-p_inferred_networks
