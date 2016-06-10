package org.panda.causalpath.data;

/**
 * Created by babur on 3/24/16.
 */
public abstract class CategoricalData extends ExperimentData
{
	public SingleQData[] data;

	public CategoricalData(String id, String symbol)
	{
		super(id, symbol);
	}

	public int[] getCategories()
	{
		int[] c = new int[data.length];

		for (int i = 0; i < c.length; i++)
		{
			c[i] = data[i].getCategory();
		}

		return c;
	}
}
