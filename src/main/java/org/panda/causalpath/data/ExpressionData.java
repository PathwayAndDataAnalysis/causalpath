package org.panda.causalpath.data;

/**
 * Numeric data array for RNA expression.
 */
public class ExpressionData extends NumericData
{
	public ExpressionData(String id, String symbol)
	{
		super(id, symbol);
	}

	@Override
	public ExperimentData copy()
	{
		ExpressionData copy = new ExpressionData(id, getGeneSymbols().iterator().next());
		copy.vals = vals;
		return copy;
	}
}
