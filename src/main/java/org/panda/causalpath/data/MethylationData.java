package org.panda.causalpath.data;

/**
 * Created by babur on 4/5/16.
 */
public class MethylationData extends NumericData
{
	public MethylationData(String id, String symbol)
	{
		super(id, symbol);
	}

	@Override
	public int getEffect()
	{
		return -1;
	}
}
