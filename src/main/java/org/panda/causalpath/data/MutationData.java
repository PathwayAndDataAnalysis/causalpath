package org.panda.causalpath.data;

import org.panda.utility.ArrayUtil;

/**
 * Data array for mutations.
 */
public class MutationData extends CategoricalData
{
	public MutationData(String id, String symbol)
	{
		super(id, symbol);
	}

	/**
	 * Gets the marking of mutated samples.
	 */
	public boolean[] getMutated()
	{
		// Why don't java 8 have boolean streams?

		boolean[] b = new boolean[data.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = data[i].category != 0 && data[i].category != ArrayUtil.ABSENT_INT;
		}
		return b;
	}

	/**
	 * Gets the marking of not-mutated samples.
	 */
	public boolean[] getNotMutated()
	{
		boolean[] b = new boolean[data.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = data[i].category == 0;
		}
		return b;
	}
}
