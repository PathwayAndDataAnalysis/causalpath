package org.panda.causalpath.data;

import java.util.Set;

/**
 * Base class for numeric data types.
 */
public abstract class NumericData extends ExperimentData
{
	/**
	 * Array for the numerical values.
	 */
	public double[] vals;

	public NumericData(String id, String symbol)
	{
		super(id, symbol);
	}

	public NumericData(String id, Set<String> geneSymbols)
	{
		super(id, geneSymbols);
	}
}
