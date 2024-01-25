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
	 * Method to reduce multiple values into one value.
	 */
	AveragingMethod avgMet;

	public ThresholdDetector(double threshold)
	{
		this.threshold = threshold;
	}

	public ThresholdDetector(double threshold, AveragingMethod avgMet)
	{
		this.threshold = threshold;
		this.avgMet = avgMet;
	}

	public void setAveragingMethod(AveragingMethod method)
	{
		this.avgMet = method;
	}

	public void setThreshold(double threshold)
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
			switch (avgMet)
			{
				case FIRST_VALUE: return nd.vals[0];
				case ARITHMETIC_MEAN: return ArrayUtil.mean(nd.vals);
				case FOLD_CHANGE_MEAN: return foldChangeGeometricMean(nd.vals);
				case MAX: return maxOfAbs(nd.vals);
			}
		}
		else if (data instanceof CategoricalData)
		{
			CategoricalData qd = (CategoricalData) data;

			switch (avgMet)
			{
				case FIRST_VALUE: return qd.getCategories()[0];
				case ARITHMETIC_MEAN: return ArrayUtil.mean(qd.getCategories());
				case FOLD_CHANGE_MEAN: throw new RuntimeException("Fold change averaging cannot be applied to categories.");
				case MAX: return maxOfAbs(qd.getCategories());
			}
		}
		return 0;
	}

	@Override
	public OneDataChangeDetector makeACopy()
	{
		return new ThresholdDetector(this.threshold, this.avgMet);
	}

	/**
	 * Gets the geometric mean of fold change values that is formatted to span the range (-inf, -1], [1, inf).
	 */
	public double foldChangeGeometricMean(double[] vals)
	{
		if (vals.length == 1) return vals[0];

		double mult = 1;
		int cnt = 0;
		for (double val : vals)
		{
			if (Double.isNaN(val)) continue;
			cnt++;
			mult *= val < 0 ? -1 / val : val;
		}
		double result = Math.pow(mult, 1D / cnt);
		if (result < 1) result = - 1 / result;
		return result;
	}

	public double maxOfAbs(double[] vals)
	{
		double max = 0;
		double v = 0;

		for (double val : vals)
		{
			double abs = Math.abs(val);
			if (abs > max)
			{
				max = abs;
				v = val;
			}
		}

		return v;
	}

	public double maxOfAbs(int[] vals)
	{
		int max = 0;
		int v = 0;

		for (int val : vals)
		{
			int abs = Math.abs(val);
			if (abs > max)
			{
				max = abs;
				v = val;
			}
		}

		return v;
	}

	public enum AveragingMethod
	{
		ARITHMETIC_MEAN,
		FOLD_CHANGE_MEAN,
		MAX,
		FIRST_VALUE
	}
}
