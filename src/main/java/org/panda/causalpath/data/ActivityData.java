package org.panda.causalpath.data;

import org.panda.resource.tcga.RPPAData;

/**
 * Created by babur on 3/25/16.
 */
public class ActivityData extends CategoricalData
{
	public ActivityData(String id, String symbol)
	{
		super(id, symbol);
	}

	public ActivityData(RPPAData rppa)
	{
		super(rppa.id, rppa.genes.iterator().next());
		data = new SingleQData[rppa.vals.length];
		for (int i = 0; i < data.length; i++)
		{
			data[i] = new Activity((int) rppa.vals[i]);
		}
	}
}
