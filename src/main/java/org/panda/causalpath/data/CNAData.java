package org.panda.causalpath.data;

/**
 * Series of copy number alterations for a gene.
 */
public class CNAData extends CategoricalData
{
	public CNAData(String id, String symbol)
	{
		super(id, symbol);
	}

	@Override
	public ExperimentData copy()
	{
		CNAData copy = new CNAData(id, getGeneSymbols().iterator().next());
		copy.data = data;
		return copy;
	}
}
