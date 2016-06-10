package org.panda.causalpath.data;

import org.panda.resource.tcga.MutTuple;

import java.util.List;

/**
 * Created by babur on 3/24/16.
 */
public class Mutation extends SingleQData
{
	protected List<MutTuple> muts;

	public Mutation(int category, List<MutTuple> muts)
	{
		super(category);
		this.muts = muts;
	}

	public Mutation(List<MutTuple> muts)
	{
		this(0, muts);
	}

	public List<MutTuple> getDetails()
	{
		return muts;
	}
}
