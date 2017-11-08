package org.panda.causalpath.data;

import org.panda.resource.tcga.ProteomicsFileRow;

import java.util.HashSet;
import java.util.Set;

/**
 * Total protein measurement.
 */
public class ProteinData extends NumericData
{
	public ProteinData(String id, Set<String> geneSymbols)
	{
		super(id, geneSymbols);
	}

	/**
	 * Converts an RPPAData object from the "resource" project.
	 */
	public ProteinData(ProteomicsFileRow row)
	{
		this(row.id, new HashSet<>(row.genes));
		this.vals = row.vals;
	}

	@Override
	public ExperimentData copy()
	{
		ProteinData copy = new ProteinData(id, getGeneSymbols());
		copy.vals = vals;
		return copy;
	}

	@Override
	public DataType getType()
	{
		return DataType.PROTEIN;
	}
}
