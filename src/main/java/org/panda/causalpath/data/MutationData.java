package org.panda.causalpath.data;

import org.panda.utility.ArrayUtil;

/**
 * Created by babur on 3/24/16.
 */
public class MutationData extends CategoricalData
{
	public MutationData(String id, String symbol)
	{
		super(id, symbol);
	}

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
