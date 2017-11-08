package org.panda.causalpath.data;

/**
 * Numeric data array for RNA expression.
 */
public class RNAData extends NumericData
{
	public RNAData(String id, String symbol)
	{
		super(id, symbol);
	}

	@Override
	public ExperimentData copy()
	{
		RNAData copy = new RNAData(id, getGeneSymbols().iterator().next());
		copy.vals = vals;
		return copy;
	}

	@Override
	public DataType getType()
	{
		return DataType.RNA;
	}
}
