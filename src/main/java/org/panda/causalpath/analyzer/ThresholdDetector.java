package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.causalpath.data.CategoricalData;
import org.panda.utility.ArrayUtil;

/**
 * Detects a change if the absolute of mean value is over the threshold.
 *
 * Created by babur on 3/24/16.
 */
public class ThresholdDetector implements OneDataChangeDetector
{
	protected double threshold;

	public ThresholdDetector(double threshold)
	{
		this.threshold = threshold;
	}

	@Override
	public int getChangeSign(ExperimentData data)
	{
		double v = getChangeValue(data);
		if (Math.abs(v) >= threshold) return v > 0 ? 1 : -1;
		else return 0;
	}

	public double getChangeValue(ExperimentData data)
	{
		if (data instanceof NumericData)
		{
			NumericData nd = (NumericData) data;
			return ArrayUtil.mean(nd.vals);
		}
		else if (data instanceof CategoricalData)
		{
			CategoricalData qd = (CategoricalData) data;
			return ArrayUtil.mean(qd.getCategories());
		}
		return 0;
	}
}
