package org.panda.causalpath.data;

import org.panda.utility.ArrayUtil;

/**
 * Created by babur on 3/24/16.
 */
public abstract class SingleQData
{
	int category;

	public static final int ABSENT = ArrayUtil.ABSENT_INT;

	public SingleQData(int category)
	{
		this.category = category;
	}

	/**
	 * Qualitative data has to come in integer categories. For instance 0 and 1 can be used for not-mutated and
	 * mutated, respectively. Or it can be inactivating mutation (-1), no or unknown mutation (0), and activating
	 * mutation (1). Some categories may be constant and some others can be customizable.
	 * @return integer category
	 */
	public int getCategory()
	{
		return category;
	}
}
