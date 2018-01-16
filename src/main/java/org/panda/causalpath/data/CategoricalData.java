package org.panda.causalpath.data;

import java.util.Set;
import java.util.stream.IntStream;

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

	public CategoricalData(String id, Set<String> symbols)
	{
		super(id, symbols);
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

	public int getNumberOfCategories()
	{
		return (int) IntStream.of(getCategories()).distinct().count();
	}

	/**
	 * Method getting the category array as parameter for efficiency purposes.
	 */
	public int getNumberOfCategories(int[] cat)
	{
		return (int) IntStream.of(cat).distinct().count();
	}
}
