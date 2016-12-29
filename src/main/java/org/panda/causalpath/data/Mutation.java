package org.panda.causalpath.data;

import org.panda.resource.tcga.MutTuple;

import java.util.List;

/**
 * List of mutations of a single gene in a single sample. The integer category of the mutations: 1 is activating, -1 is
 * inactivating, and 0 is unknown effect.
 */
public class Mutation extends SingleCategoricalData
{
	/**
	 * List of detected mutations.
	 */
	protected List<MutTuple> muts;

	public Mutation(int category, List<MutTuple> muts)
	{
		super(category);
		this.muts = muts;
	}

	/**
	 * Constructor for the cases where the effect of mutation is not known.
	 */
	public Mutation(List<MutTuple> muts)
	{
		this(0, muts);
	}

	public List<MutTuple> getDetails()
	{
		return muts;
	}
}
