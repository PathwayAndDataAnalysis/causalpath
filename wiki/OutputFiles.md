# Possible output files from a CausalPath run and their description 

CausalPath takes a directory name as the program argument, and looks for the `parameters.txt` file in that directory to load all the parameters, including the locations of input data files. The output is generated into the same directory where the `parameters.txt` sits. Below is the list of possible output files that are generated as CasalPath results. They are all tab-delimited text files except the two json files.

**`results.txt`:** Each row provides the prior binary relation and the matching experiment data. Correlation-based and comparison-based analyses use slightly different formats for this file. A correlation-based analysis will provide the correlation value and its p-value, while the comparison-based analysis will provide the change amount and p-values for both source and target data changes corresponding to the row.

**`causative.sif`:** A SIF result file for identified causal relations that can be loaded into and visualized by ChiBE. Make sure the `causative.format` file is also in the same directory.

**`causative.format`:** This file complements `causative.sif` by providing visual features of the result network.

**`causative.json`:** Alternative to loading into ChiBE, users can use the web service to visualize the result networks. Web-service reuires this json file to render the result causal network.

**`causative-data-cetric.sif`:** This file is only generated if the user sets `generate-data-centric-graph = true`. Different from the gene-centric network `causative.sif`, this SIF network uses nodes for the data IDs. Users can visualize this network by loading into ChiBE.

**`causative-data-cetric.format`:** This file complements `causative-data-centric.sif` by providing visual features.

**`conflicting.sif`:** In addition to the causal result network, CausalPath generates a network showing the causal priors with contradicting experiment data. Here contradicting means that the data change is significant within the thresholds, but the change direction is opposite to what we expect by the causal priors. 

**`conflicting.format`:** This file complements `conflicting.sif` by providing visual features of the network.

**`conflicting.json`:** Alternative to loading into ChiBE, users can use the web service to visualize the result networks. Web-service requires this json file to render the result conflicting network.

**`unknown-site-effects.txt`:** One of the major limitations of CausalPath is not knowing functional effects of every measured site. Users have the option to manually curate these effects from the literature and add them into the proteomic data file. This file provides the list of sites with unknown effect, and also there are changes at the downstream that could be explained by the upstream changes only if we knew these site effects. If users manually curate the effect of a site in this file, it can produce an additional edge in the final network with 50% chance.

**`pval-uniformity.txt`:** Whenever the FDR is controlled within the CausalPath analysis, this file is generated to be a graphical guide for the user showing the intensity of the signal in the data. Users can plot this table in a spreadsheet program to see the deviation of the ranked p-values from the x=y line. The more deviation toward x-axis means a signal with higher intensity.

**`significance-pvals.txt`:** This file is generated only if network significance is calculated. The first row provides the overall network significance. The following table provides p-value(s) for each gene downstream. For a correlation-based analysis there is only one p-value per gene, for a comparison-based analysis there are 3 p-values for each gene corresponding the downstream overall size, activation-suggesting downstream and inhibition suggesting downstream. CausalPath uses these values after correcting for multiple hypothesis testing.

**`value-changes.txt`:** This file shows the value changes for each loaded data in a comparison-based analysis. If the original data is used by applying a t-test on the data, then this file will show the t-value and p-value coming from the test, per data row.
