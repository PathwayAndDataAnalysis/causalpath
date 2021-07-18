package org.panda.causalpath.data;

import org.panda.resource.tcga.ProteomicsFileRow;

import java.util.HashSet;

/**
 * Numeric data array for RNA expression.
 */
public class RNAData extends NumericData
{
	public RNAData(String id, String symbol)
	{
		super(id, symbol);
	}

	public RNAData(ProteomicsFileRow row)
	{
		super(row.id, new HashSet<>(row.genes));
		this.vals = row.vals;
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
