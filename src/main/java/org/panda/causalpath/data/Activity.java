package org.panda.causalpath.data;

/**
 *  1: Activated
 *  0: Unknown
 * -1: Inactivated
 */
public class Activity extends SingleCategoricalData
{
	public Activity(int effect)
	{
		super(effect);
	}
}
