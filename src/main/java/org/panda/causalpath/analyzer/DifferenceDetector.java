package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.causalpath.data.CategoricalData;
import org.panda.utility.ArrayUtil;

/**
 * For change detection based on difference of mean values of two groups.
 */
public class DifferenceDetector extends ComparisonDetector
{
	public DifferenceDetector(double threshold, boolean[] control, boolean[] test)
	{
		super(threshold, control, test);
	}

	public double getChangeValue(ExperimentData data)
	{
		if (data instanceof NumericData)
		{
			NumericData nd = (NumericData) data;
			return ArrayUtil.mean(nd.vals, test) - ArrayUtil.mean(nd.vals, control);
		}
		else if (data instanceof CategoricalData)
		{
			CategoricalData qd = (CategoricalData) data;
			return ArrayUtil.mean(qd.getCategories(), test) - ArrayUtil.mean(qd.getCategories(), control);
		}
		return 0;
	}

	@Override
	public OneDataChangeDetector makeACopy()
	{
		return new DifferenceDetector(threshold, control, test);
	}
}
