package org.panda.causalpath.data;

import java.util.Set;

/**
 * Created by babur on 3/24/16.
 */
public abstract class NumericData extends ExperimentData
{
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
