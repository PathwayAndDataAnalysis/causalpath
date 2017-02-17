package org.panda.causalpath.data;

import org.panda.utility.ArrayUtil;

/**
 * Data array for mutations.
 */
public class MutationData extends CategoricalData
{
	/**
	 * The effect of the mutation that are represented with this data. If both activating and inhibiting mutations are
	 * considered in an analysis for the same gene, then different MutationData objects should be created for each
	 * effect.
	 */
	private int effect;

	public MutationData(String id, String symbol, int effect)
	{
		super(id, symbol);
		this.effect = effect;
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

	@Override
	public int getEffect()
	{
		return effect;
	}

	@Override
	public boolean isSiteSpecific()
	{
		return true;
	}

	@Override
	public ExperimentData copy()
	{
		MutationData copy = new MutationData(id, getGeneSymbols().iterator().next(), getEffect());
		copy.data = data;
		return copy;
	}
}
