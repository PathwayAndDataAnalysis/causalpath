package org.panda.causalpath.data;

/**
 * Holds an array of categorical data points.
 */
public abstract class CategoricalData extends ExperimentData
{
	public SingleCategoricalData[] data;

	public CategoricalData(String id, String symbol)
	{
		super(id, symbol);
	}

	/**
	 * Sometimes we need the data array to be less structured for practical purposes. This method represents the
	 * category array as an integer array.
	 */
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
