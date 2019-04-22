# Inputs and their format for CausalPath

Users need to prepare a minimum of 2 text files to run a CausalPath analysis.

1. Proteomics data
2. Analysis parameters

## Proteomics data

CausalPath reads the proteomics dataset from a tab-delimited text file where the first row contains the column headers. Below are description of those columns.

**ID:** A unique text identifier for each row in the dataset. Ideally, it is better if the ID contains the gene symbols and modification sites if applicable. Those IDs are displayed in the tooltip text of related network parts by visualization software. Make sure that IDs do not contain the characters "|" and "@". Otherwise ChiBE won't be able to display them.

**Symbols:** HGNC symbol of the related gene. If there are more than one gene symbol reltaed to this row, there should be a single space between each symbol.

**Sites:** If the row contains a phosphoprotein measurement, this column should indicate protein sites that are affected. The format for the site is one letter capital aminoacid code, following an integer for the location on the UniProt caconical sequence, such as `Y142`, or `S78`. When there is more than one related site, they should be separated with a pipe (`|`), such like `S151|T153`, to form a site group. If the row is related to more than one gene symbol, there should exists a site group for each symbol, separated with a single space.

**Effect:** If the row is a phosphoprotein measurement, this column can contain the effect of the related phosphorylation on the protein activity. Please use `a` for activating and `i` for inhibiting phosphorylations. If the effect is too complex for a simple classification, then you can leave it blank (preferred) or use a `c` for the complex effect. But `c` will remove this row from possible causes, as CausalPath cannot evaluate complex effects. If this column is left blank for a row, then CausalPath looks it up from its database, which we compiled from PhosphoSitePlus and other resources.

**Value:** The numeric value for the row. There can be more than one value columns, but of course each of them with a unique name. There are many ways to encode values in the data file. They may represent normalized protein reads, or they can be comparison values like fold changes. The nature of the values has to be specified in the parameters file.

## Parameters file

The name of the parameters file have to be `parameters.txt` exactly. Each parameter in this file should be given in a separate line, in the format `parameter-name = parameter-value`. Below are list of possible parameters and their description.

`proteomics-values-file`: Name of the proteomics values file. It should have at least one ID column and one or more columns for experiment values. Platform file and values file can be the same file.

`proteomics-platform-file`: Name of the proteomics platform file. Each row should belong to either a gene's total protein measurement, or a site specific measurement. This file should contain ID, gene symbols, modification sites, and known site effects. Platform file and values file can be the same file.

`id-column`: The name of the ID column in platform and values files.

`symbols-column`: The name of the symbols column in platform file.

`sites-column`: The name of the sites column.

`effect-column`: The name of the effect column.

`value-transformation`: This parameter determines how to use the values in the proteomics file. Options are listed below. When there is only one value column and no transformation is desired, users can select any of the first 3 options, as they have no effect on a single value.

	arithmetic-mean: The arithmetic mean value of the given values is used for significance detection of a single change. There should only be one group of values (marked with value-column), the values have to be distributed around zero, and a threshold value should be provided for significance detection, using the threshold-for-data-significance.

	geometric-mean: The geometric mean value of the given values is used for significance detection of a single change. This is the only case when the geometric mean is used for averaging a group of samples, and it is appropriate if the individual values are formed of some kind of ratios. There should only be one group of values (marked with value-column), the values have to be distributed around zero, and a threshold value should be provided for significance detection, using the threshold-for-data-significance.

	max: The value with maximum absolute is used for the analysis. There should only be one group of values (marked with value-column), the values have to be distributed around zero, and a threshold value should be provided for significance detection, using the threshold-for-data-significance.

	difference-of-means: There should be control and test values, whose difference would be used for significance detection. The threshold for significance (threshold-for-data-significance) should also be provided.

	fold-change-of-mean: There should be control and test values, whose ratio will be converted to fold change and thresholded. The fold change value will be in the range (-inf, -1] + [1, inf). If the data file already contains a fold-change value, then please use the geometric-mean as value transformation. The threshold for significance (threshold-for-data-significance) should also be provided.

	significant-change-of-mean: There should be sufficient amount of control and test values to detect the significance of change with a t-test. Technically there should be more than 3 controls and 3 tests, practically, they should be much more to provide statistical power. The threshold-for-data-significance should be used for a p-value threshold, or alternatively, fdr-threshold-for-data-significance should be used for controlling significance at the false discovery rate level.

	correlation: There should be one group of values (marked with value-column). There must be at least 3 value columns technically, but many more than that practically to have some statistical power for significant correlation. 



`value-column`: Name of a value column. This parameter should be used when there is only one group of experiments to consider in the analysis.

`control-value-column`: Name of a control value column. This parameter should be used when there are control and test value columns in the dataset.

`test-value-column`: Name of a test value column. This parameter should be used when there are control and test value columns in the dataset.

`do-log-transform`: Whether the proteomic values should be log transformed for the analysis. Possible values are 'true' and 'false'. Default is false.

`threshold-for-data-significance`: A threshold value for selecting significant data. Use this parameter only when FDR controlling procedure is already performed outside of CausalPath. This parameter can be set for each different data type separately. The parameter value has to be in the form 'thr-val data-type', such like '1 phosphoprotein' or '2 protein.

`fdr-threshold-for-data-significance`: False discovery rate threshold for data significance. This parameter can be set for each different data type separately. The parameter value has to be in the form 'fdr-val data-type', such like '0.1 phosphoprotein' or '0.05 protein.

`pool-proteomics-for-fdr-adjustment`: Whether to consider proteomic and phosphoproteomic data as a single dataset during FDR adjustment. This is typically the case with RPPA data, and typically not the case with mass spectrometry data. Can be 'true' or 'false'. Default is false.

`correlation-value-threshold`: Option to control correlation with its value. This cannot be used with FDR control, but can be used with p-value control.

`correlation-upper-threshold`: In some types of proteomic data, highest correlations come from errors. A way around is filtering with an upper value.

`pval-threshold-for-correlation`: A p-value threshold for correlation in a correlation-based causality. This parameter should only be used when FDR control is performed outside of CausalPath.

`fdr-threshold-for-correlation`: False discovery rate threshold for the correlations in a correlation-based analysis.

`stdev-threshold-for-data`: This parameter can be set for each different data type separately. The parameter value has to be in the form 'stdev-thr data-type', such like '0.5 phosphoprotein'.

`default-missing-value`: An option to specify a default value for the missing values in the proteomics file.

`minimum-sample-size`: When there are missing values in proteomic file, the comparisons can have different sample sizes for controls and tests. This parameter sets the minimum sample size of the control and test sets.

`calculate-network-significance`: Whether to calculate significances of the properties of the graph. When turned on, a p-value for network size, and also downstream activity enrichment p-values for each gene on the graph are calculated.

`permutations-for-significance`: We will do data randomization to see if the result network is large, or any protein's downstream is enriched. This parameter indicates the number of randomizations we should perform. It should be reasonable high, such as 1000, but not too high.

`fdr-threshold-for-network-significance`: The false discovery rate for network significance calculations for the downstream activity enrichment of genes.

`use-network-significance-for-causal-reasoning`: After calculation of network significances in a non-correlation-based analysis, this option introduces the detected active and inactive proteins as data to be used in the analysis. This applies only to the proteins that already have a changed data on them, and have no previous activity data associated.

`minimum-potential-targets-to-consider-for-downstream-significance`: While calculating downstream significance for each source gene, we may not like to include those genes with just a few qualifying targets to reduce the number of tested hypotheses. These genes may not be significant even all their targets are in the results, and since we use Benjamini-Hochberg procedure to control false discovery rate from multiple hypothesis testing, their presence will hurt the statistical power. Use this parameter to exclude genes with few qualifying targets on the network. Default is 5.

`do-site-matching`: Whether to force site matching in causality analysis. True by default.

`site-match-proximity-threshold`: Phosphorylation relations many times know the target sites. When we observe a change in a site of the target protein which is not targeted by the relation, but the site is very close to a known target site of the relation, this parameter let's us to assume that the relation also applies to those close-by sites.

`site-effect-proximity-threshold`: CausalPath has a database of phosphorylation site effects. When not set, this parameter is 0 by default, which means exact usage of site effects. But sometimes we may see a site with unknown effect is modified, which is very close to another site with a known effect. This parameter let's us to assume those changing sites with unknown effect has the same effect with the neighbor site with known effect. Use responsibly.

`built-in-network-resource-selection`: Determines which network resource to use during the analysis. Multiple network resource should be mentioned together, separated with a space or comma. Possible values are below.

	PC: Pathway Commons v9 for all kinds of relations.

	REACH: Network derived from REACH NLP extraction results for phosphorylation relations.

	PhosphoNetworks: The PhosphoNetworks database for phosphorylations.

	IPTMNet: The IPTMNet database for phosphorylations.

	TRRUST: The TRRUST database for expression relations.

	TFactS: The TFactS database for expression relations.



`relation-filter-type`: Use this parameter to limit the results with a specific type of relation. Possible values are below.

	no-filter: The graph is used with all inferred relations.

	phospho-only: Only phosphorylation relations are desired.

	expression-only: Only expression relations are desired.

	phospho-primary-expression-secondary: All phosphorylation relations are desired. Expression relations are desired only as supplemental, i.e., they have to involve at least one protein that exists in phospho graph.



`gene-focus`: Use this parameter to crop the result network to the neighborhood of certain gene. You should provide gene symbols of these genes in a row separated by a semicolon, such like 'MTOR;RPS6KB1;RPS6'

`mutation-effect-file`: When we have mutations in the analysis, users can provide mutation effects using this parameter, otherwise all mutations are assumed to be inactivating.

`color-saturation-value`: Specifies the value where node colors reach most intense color. Has to be a positive value, and used symmetrically. In the case of value-transformation is significant-change-of-mean, the value is -log(p) with a sign associated to it.

`show-all-genes-with-proteomic-data`: CausalPath generates a result graph, but what about all other significant changes that could not make into the network? CausalPath puts those genes as disconnected nodes in the graph when the analysis is not correlation based. This is true by default but can be turned off by setting to false.

`show-insignificant-data`: Option to make the insignificant protein data on the result graph visible. Seeing these is good for seeing what is being measured, but when they are too much, turning off generates a a better view.

`hide-data-not-part-of-causal-relations`: Limits the data drawn on the result graph to the ones that take part in the identified causal relations.

`data-type-for-expressional-targets`: By default, CausalPath generates explanations only for proteomic changes. But it is possible to explain RNA changes with expressional relations as well, and it is a more direct explanation than total protein measurement. This parameter lets users to control possible data types explainable by expressional relations. Typical values are 'rna' and 'protein'. This parameter can also  be used multiple times to use rna and protein data together.

`generate-data-centric-graph`: An alternative to the gene-centric graph of CausalPath is a data-centric graph where nodes are not genes but the data. This parameter forces to generate this type of result as well. False by default.

`gene-activity`: Use this parameter to assign a specific activity or inactivity to a gene in the analysis. The value has to start with a gene name and one letter code for activity or inactivity, such as 'BRAF a', or 'PTEN i'.

`tf-activity-file`: CausalPath lets users to input results from an inference for transcriptional factor activities, such as PRECEPTS, MARINa or VIPER. For this, the results should be prepared in a file, first column containing TF symbol and the second column whether 'activated' or 'inhibited'. The name of such file should be provided here.

`use-strongest-proteomic-data-per-gene`: When a proteomic experiment outputs too many phosphorylation sites with lots of changes, many proteins have evidences for both activation and inhibition. This produces networks hard to read. A complexity management technique is to turn on this parameter to use only the strongest proteomic feature at the upstream of relations. This is false by default.

`use-missing-proteomic-data-for-test`: Option to use a G-test to check unequal distribution of missing values. If opted, and data is sufficient, the G-test result is combined with t-test result with Fisher's method. But beware. This method assumes that missing values are uniformly distributed to samples. If this is violated, then false positives will appear. If you are not sure, stay away from this option.

`randomized-matrix-directory-for-missing-proteomic-data`: Using randomization is an alternative to using a G-test for interpreting missing data distribution. CausalPath cannot generate those matrices, but it can use pre-generated matrices to compute significances. This operation typically requires a lot of memory.

`missing-value-test-data-sufficiency-threshold`: When we use a G-test in the analysis, we don't want to use it for every proteomic row. Some rows will have insufficient data. To test for sufficiency, we generate an extreme case where missing data is shifted to the smaller group and see if this can provide a p-value small enough. If this extreme shifting cannot make the the p-value small enough (specified with this threshold), then we don't use a G-test for that row.

`custom-resource-directory`: CausalPath downloads some data in the first run and stores in the resource directory. This directory is '.panda' by default. If this needs to be customized, use this parameter.

`tcga-directory`: It is possible to add genomic data from TCGA to CausalPath analysis. This is only useful when the proteomic data have the same sample IDs. Users can load TCGA data into a local directory from Broad Firehose, and provide the directory here. The org.panda.resource.tcga.BroadDownloader in the project https://github.com/PathwayAndDataAnalysis/resource is a utility that can do that.
