package org.panda.causalpath.data;

import org.panda.resource.tcga.ProteomicsFileRow;

import java.util.HashSet;

/**
 * Numeric data array for metabolite abundances. While this class reuses the class structure of other data types, the
 * "gene symbol" field actually stores ChEBI IDs.
 */
public class MetaboliteData extends NumericData
{
	public MetaboliteData(String id, String symbol)
	{
		super(id, symbol);
	}

	public MetaboliteData(ProteomicsFileRow row)
	{
		super(row.id, new HashSet<>(row.genes));
		this.vals = row.vals;
	}

	@Override
	public ExperimentData copy()
	{
		MetaboliteData copy = new MetaboliteData(id, getGeneSymbols().iterator().next());
		copy.vals = vals;
		return copy;
	}

	@Override
	public DataType getType()
	{
		return DataType.METABOLITE;
	}
}
