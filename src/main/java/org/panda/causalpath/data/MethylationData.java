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

	@Override
	public ExperimentData copy()
	{
		MethylationData copy = new MethylationData(id, getGeneSymbols().iterator().next());
		copy.vals = vals;
		return copy;
	}

	@Override
	public DataType getType()
	{
		return DataType.DNA_METHYLATION;
	}
}
