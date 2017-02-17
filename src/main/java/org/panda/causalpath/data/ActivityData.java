package org.panda.causalpath.data;

import org.panda.resource.tcga.ProteomicsFileRow;

/**
 * Array of activity values for a gene.
 */
public class ActivityData extends CategoricalData
{
	public ActivityData(String id, String symbol)
	{
		super(id, symbol);
	}

	/**
	 * This constructor is for converting the activity-encoding RPPAData in the resource module to an activity data to
	 * use in this project.
	 */
	public ActivityData(ProteomicsFileRow rppa)
	{
		super(rppa.id, rppa.genes.iterator().next());
		data = new SingleCategoricalData[rppa.vals.length];
		for (int i = 0; i < data.length; i++)
		{
			data[i] = new Activity((int) rppa.vals[i]);
		}
	}

	@Override
	public ExperimentData copy()
	{
		return new ActivityData(id, getGeneSymbols().iterator().next());
	}
}
