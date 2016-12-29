package org.panda.causalpath.data;

/**
 * Data array for methylation.
 */
public class MethylationData extends NumericData
{
	public MethylationData(String id, String symbol)
	{
		super(id, symbol);
	}

	/**
	 * A methylated gene may be under expressed, hence, inactivated.
	 */
	@Override
	public int getEffect()
	{
		return -1;
	}
}
