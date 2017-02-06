# Inputs and their format for CausalPath

Users need to prepare a minimum of 2 text files to run a CausalPath analysis.

1. Proteomics data
2. Analysis parameters

## Proteomics data

CausalPath reads the proteomics dataset from a tab-delimited text file where the first row contains the column headers. Below are description of those columns.

**ID:** A unique text identifier for each row in the dataset. Ideally, it is better if the ID contains the gene symbols and modification sites if applicable. Those IDs are displayed in the tooltip text of related network parts by visualization software.

**Symbols:** HGNC symbol of the related gene. If there are more than one gene symbol reltaed to this row, there should be a single space between each symbol.

**Sites:** If the row contains a phosphoprotein measurement, this column should indicate protein sites that are affected. The format for the site is one letter capital aminoacid code, following an integer for the location on the UniProt caconical sequence, such as `Y142`, or `S78`. When there is more than one related site, they should be separated with a pipe (`|`), such like `S151|T153`, to form a site group. If the row is related to more than one gene symbol, there should exists a site group for each symbol, separated with a single space.

**Effect:** If the row is a phosphoprotein measurement, this column can contain the effect of the related phosphorylation on the protein activity. Please use `a` for activating, `i` for inhibiting phosphorylations. If the effect is too complex for a simple classification, then you can use `c` for complex effect. But `c` will remove this row from analysis, as CausalPath cannot evaluate complex effects. If this column is left blank for a row, then CausalPath looks it up from its database, which we compiled from PhosphoSitePlus and other resources.

**Value:** The numeric value for the row. There can be more than one value columns, but of course each of them with a unique name. There are many ways to encode values in the data file. They may represent normalized protein reads, or they can be comparison values like fold changes. The nature of the values has to be specified in the parameters file.

## Parameters file

The name of the parameters file have to be `parameters.txt` exactly. Each parameter in this file should be given in a separate line, in the format `parameter-name = parameter-value`. Below is a sample of initial lines in a parameters file, assuming the data is in a file named `ProteomicsData.txt`.

```
proteomics-values-file = ProteomicsData.txt
id-column = ID
symbols-column = Symbols
sites-column = Sites
effect-column = Effect
```

The names of the columns in the data file can be customized, as long as correctly specified in the parameters file, as above.

The next essential parameter for the parameters file is the indication of how to use the given value(s) in the data file. Don't forget that the causality analysis makes use of value comparisons, and it can not work with a single column of raw measurements. There has to be a comparison in the values, or this software should be able to make a comparison using the given values. We call this a *value transformation*, whose type should be indicated with the parameter `value-transformation`. Below are its possible values. A parameters file has to contain exactly one of the below lines.
```
value-transformation = arithmetic-mean
value-transformation = geometric-mean
value-transformation = difference-of-means
value-transformation = fold-change-of-mean
value-transformation = significant-change-of-mean
value-transformation = correlation
```
**Arithmetic mean:** There is one group of values distributed around zero, and zero means no change. The arithmetic mean of those value columns will be used as the value in the analysis. Users should also supply a threshold value which would mark a significant change. The threshold is supposed to be a positive value, and its mirror negative will be used for downregulations. Number of value columns can be one or more. Below is a sample.
```
value-transformation = arithmetic-mean
value-column = sample1
value-column = sample2
value-column = sample3
threshold-for-data-significance = 2.5
```
**Geometric mean:** There is one group of values, and each value is a fold-change in the range (-inf, -1],[1, inf), where decrease in fold change represented with negative values. The averaging method first maps the negative values into range [0,1], takes the geometric mean, then re-maps the values between 0 and 1 to the range (-inf, -1]. Users should provide a thresold for fold change significance. Below is an example.
```
value-transformation = geometric-mean
value-column = sample1
value-column = sample2
value-column = sample3
threshold-for-data-significance = 2
```
**Difference of means:** There are two groups of values: controls and tests, where each value is a normalized measurement. CausalPath will use the difference of the test mean from the control mean, where the means are arithmetic. Users need to provide a significance threshold for the difference. Below is an example.
```
value-transformation = difference-of-means
test-value-column = sample1
test-value-column = sample2
control-value-column = sample3
control-value-column = sample4
threshold-for-data-significance = 5
```
**Fold change of mean:** There are control and test values, where each value is a normalized measurement. The method will calculate the fold change of the mean of the test values compared to the mean of the control values. Users need to provide a significance threshold for the fold change, which is greater than 1. The threshold will be used symmetrically for the negative fold changes. Below is an example.
```
value-transformation = fold-change-of-mean
test-value-column = sample1
test-value-column = sample2
control-value-column = sample3
control-value-column = sample4
threshold-for-data-significance = 1.5
```
**Correlation:** There is one group of values, technically more than 3 samples, practically many samples. The Pearson correlation will be used for the value. Users need to provide a positive threshold for the correlation value (between 0 and 1), and another threshold for the p-value of the correlation (between 0 and 1). Below is an example.
```
value-transformation = correlation
value-column = sample1
value-column = sample2
value-column = sample3
value-column = sample4
threshold-for-data-significance = 0.01
threshold-for-correlation = 0.5
```


... This page will be completed very soon ...
