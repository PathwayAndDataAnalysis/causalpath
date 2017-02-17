package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.causalpath.data.CategoricalData;
import org.panda.utility.ArrayUtil;

/**
 * Detects a change if the absolute of mean value is over the threshold.
 */
public class ThresholdDetector implements OneDataChangeDetector
{
	/**
	 * Threshold to use.
	 */
	protected double threshold;

	/**
	 * Whether to use the geometric mean for averaging the values. THis is false by default and only makes sense if the
	 * values are some kind of ratios.
	 */
	protected boolean geometricMean;

	public ThresholdDetector(double threshold)
	{
		this.threshold = threshold;
		this.geometricMean = false;
	}

	public void setGeometricMean(boolean geometricMean)
	{
		this.geometricMean = geometricMean;
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
			return geometricMean ? foldChangeGeometricMean(nd.vals) : ArrayUtil.mean(nd.vals);
		}
		else if (data instanceof CategoricalData)
		{
			CategoricalData qd = (CategoricalData) data;
			return ArrayUtil.mean(qd.getCategories());
		}
		return 0;
	}

	/**
	 * Gets the geometric mean of fold change values that is formatted to span the range (-inf, -1], [1, inf).
	 */
	public double foldChangeGeometricMean(double[] vals)
	{
		double mult = 1;
		int cnt = 0;
		for (double val : vals)
		{
			if (Double.isNaN(val)) continue;
			cnt++;
			mult *= val < 0 ? -1 / val : val;
		}
		return Math.pow(mult, 1D / cnt);
	}
}
