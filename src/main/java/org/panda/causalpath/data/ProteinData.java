package org.panda.causalpath.data;

import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.ArrayUtil;

import java.util.HashSet;
import java.util.Set;

/**
 * Total protein measurement.
 */
public class ProteinData extends NumericData
{
	PresenceData pres;

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

	public void initPresenceData(boolean[] consider)
	{
		pres = new PresenceData(id, getGeneSymbols(), getType());
		pres.data = new SingleCategoricalData[vals.length];
		for (int i = 0; i < vals.length; i++)
		{
			pres.data[i] = new Presence(!consider[i] ? ArrayUtil.ABSENT_INT : Double.isNaN(vals[i]) ? 0 : 1);
		}
	}

	@Override
	public ExperimentData copy()
	{
		ProteinData copy = new ProteinData(id, getGeneSymbols());
		copy.vals = vals;
		if (pres != null) copy.pres = (PresenceData) pres.copy();
		return copy;
	}

	@Override
	public DataType getType()
	{
		return DataType.PROTEIN;
	}

	public PresenceData getPresenceData()
	{
		return pres;
	}
}
