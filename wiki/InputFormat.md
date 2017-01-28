# Inputs and their format for CausalPath

Users need to prepare a minimum of 2 text files to run a CausalPath analysis.

1. Proteomics data
2. Analysis parameters

## Proteomics data

CausalPath reads the proteomics dataset from a tab-delimited text file where the first row contains the column headers. Below are description of those columns.

**ID:** A unique text identifier for each row in the dataset. Ideally, it is better if the ID contains the gene symbols and modification sites if applicable. Those IDs are displayed in the tooltip text of related network parts by visualization software.

**Symbols:** HGNC symbol of the related gene. If there are more than one gene symbol reltaed to this row, there should be a single space between each symbol.

**Sites:** If the row contains a phosphoprotein measurement, this column should indicate protein sites that are affected. The format for the site is one letter capital aminoacid code, following an integer for the location on the UniProt caconical sequence, such as `Y142`, or `S78`. When there is more than one related site, they should be separated with a pipe (`|`), such like `S151|T153`, to form a site group. If the row is related to more than one gene symbol, there should exists a site group for each symbol, separated with a single space.

**Effect:** If the row is a phosphoprotein measurement, this column can contain the effect of the related phosphorylation on the protein activity. Please use `a` for activating, `i` for inhibiting phosphorylations (without the quotes). If the effect is too complex for a simple classification, then you can use `c` for complex effect. But `c` will remove this row from analysis, as CausalPath cannot evaluate complex effects. If this column is left blank for a row, then CausalPath looks it up from its database, which we compiled from PhosphoSitePlus and other resources.

**Value:** The numeric value for the row. There can be more than one value columns, but of course each of them with a unique name. There are many ways to encode values in the data file. They may represent normalized protein reads, or they can be comparison values like fold changes. The nature of the values has to be specified in the parameters file.

## Parameters file

The name of the parameters file have to be `parameters.txt` exactly. Each parameter in this file should be given in a separate line, in the format `parameter-name = parameter-value`. Below is a sample of initial lines in a parameters file, assuming the data is in a file named `ProteomicsData.txt`.

```
proteomics-platform-file = ProteomicsData.txt
proteomics-values-file = ProteomicsData.txt
id-column = ID
symbols-column = Symbols
sites-column = Sites
effect-column = Effect
```
Why the data file is repeated? The previous section about data file directs to put gene annotations (first 4 columns) and values in the same file, but this doesn't have to be the case. In the case of RPPA experiments, the first 4 columns are constant for each chip, so it should not have to be repeated for each time. In that case, the annotation can stay in a platform file and values can go in another file, with the condition that both will contain the ID column, with the same header.

The names of the columns in the data file can be customized, as long as correctly specified in the parameters file, as above.



... This page will be completed very soon ...
