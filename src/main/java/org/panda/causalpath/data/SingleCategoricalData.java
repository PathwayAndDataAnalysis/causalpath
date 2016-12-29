package org.panda.causalpath.data;

import org.panda.utility.ArrayUtil;

/**
 * Categorical data is encoded with integer values. For instance 0 and 1 can be used for not-mutated and
 * mutated, respectively. Or it can be inactivating mutation (-1), no or unknown mutation (0), and activating
 * mutation (1). This class holds a single categorical data point.
 */
public abstract class SingleCategoricalData
{
	/**
	 * Integer value for the category.
	 */
	int category;

	/**
	 * Special integer designated to indicate the data is absent, i.e. not measured.
	 */
	public static final int ABSENT = ArrayUtil.ABSENT_INT;

	public SingleCategoricalData(int category)
	{
		this.category = category;
	}

	public int getCategory()
	{
		return category;
	}
}
