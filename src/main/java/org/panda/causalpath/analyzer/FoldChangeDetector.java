package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.utility.ArrayUtil;

/**
 * Change detector based on fold difference.
 */
public class FoldChangeDetector extends ComparisonDetector
{
	public FoldChangeDetector(double threshold, boolean[] control, boolean[] test)
	{
		super(threshold, control, test);
	}

	public double getChangeValue(ExperimentData data)
	{
		if (data instanceof NumericData)
		{
			NumericData nd = (NumericData) data;
			double t = ArrayUtil.mean(nd.vals, test);
			double c = ArrayUtil.mean(nd.vals, control);

			if (t > c) return t / c;

			return -(c / t);
		}
		throw new RuntimeException("Non-numeric data is trying to use a fold-change detector. data = " + data);
	}

	@Override
	public OneDataChangeDetector makeACopy()
	{
		return new FoldChangeDetector(threshold, control, test);
	}
}
