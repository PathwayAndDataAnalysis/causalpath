package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.causalpath.data.CategoricalData;
import org.panda.utility.ArrayUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Anova;
import org.panda.utility.statistics.ChiSquare;
import org.panda.utility.statistics.Correlation;

import java.util.List;

/**
 * Change detector for experiment data pairs based on their correlation.
 */
public class CorrelationDetector implements TwoDataChangeDetector
{
	/**
	 * Absolute value of the minimum correlation to be detected.
	 */
	protected double correlationThreshold;

	/**
	 * Significance threshold for correlations to be detected.
	 */
	protected double pvalThreshold;

	protected int minimumSampleSize;

	protected double correlationUpperThreshold;

	/**
	 * If no threshold is needed, pass -1 for that threshold.
	 */
	public CorrelationDetector(double correlationThreshold, double pvalThreshold)
	{
		this.correlationThreshold = correlationThreshold;
		this.pvalThreshold = pvalThreshold;
		minimumSampleSize = 3;
		correlationUpperThreshold = -1; // negative value meaning not set
	}

	public void setCorrelationThreshold(double correlationThreshold)
	{
		this.correlationThreshold = correlationThreshold;
	}

	public void setPvalThreshold(double pvalThreshold)
	{
		this.pvalThreshold = pvalThreshold;
	}

	public void setMinimumSampleSize(int minimumSampleSize)
	{
		this.minimumSampleSize = minimumSampleSize;
	}

	public void setCorrelationUpperThreshold(double correlationUpperThreshold)
	{
		this.correlationUpperThreshold = correlationUpperThreshold;
	}

	@Override
	public int getChangeSign(ExperimentData data1, ExperimentData data2)
	{
		// Case 1: Both data are numeric

		if (data1 instanceof NumericData && data2 instanceof NumericData)
		{
			NumericData nd1 = (NumericData) data1;
			NumericData nd2 = (NumericData) data2;

			double[][] v = ArrayUtil.trimNaNs(nd1.vals, nd2.vals);

			if (v[0].length < minimumSampleSize) return 0;

			Tuple corr = Correlation.pearson(v[0], v[1]);

			if (pvalThreshold >= 0)
			{
				if (corr.p > pvalThreshold) return 0;
			}

			// it passes p-val thr at this point, if such a threshold exists.

			if (correlationThreshold > 0 && Math.abs(corr.v) < correlationThreshold) return 0;
			if (correlationUpperThreshold > 0 && Math.abs(corr.v) > correlationUpperThreshold) return 0;
			return (int) Math.signum(corr.v);
		}

		// Case 2: One data is numeric, the other is categorical

		else if (data1 instanceof CategoricalData && data2 instanceof NumericData ||
			data2 instanceof CategoricalData && data1 instanceof NumericData)
		{
			CategoricalData qd = (CategoricalData) (data1 instanceof CategoricalData ? data1 : data2);
			NumericData nd = (NumericData) (data1 == qd ? data2 : data1);

			int[] cat = qd.getCategories();
			double[] val = nd.vals;

			if (cat.length != val.length)
			{
				System.out.println();
			}

			List<double[]> groups = ArrayUtil.separateToCategories(cat, val);
			double pval = Anova.oneWayPval(groups);

			if (pval > pvalThreshold) return 0;

			double[][] v = ArrayUtil.trimNaNs(ArrayUtil.toDouble(cat), val);
			double corr = Correlation.pearsonVal(v[0], v[1]);
			if (correlationThreshold > 0 && Math.abs(corr) < correlationThreshold) return 0;
			return (int) Math.signum(corr);

//			double[][] two = ArrayUtil.separateToBalancedTwo(groups);
//			return ArrayUtil.diffOfMeans(two[0], two[1]) > 0 ? 1 : -1;
		}

		// Case 3: Both data are categorical

		else if (data1 instanceof CategoricalData && data2 instanceof CategoricalData)
		{
			CategoricalData qd1 = (CategoricalData) data1;
			CategoricalData qd2 = (CategoricalData) data2;

			int[] cat1 = qd1.getCategories();
			int[] cat2 = qd2.getCategories();
			double pval = ChiSquare.testDependence(cat1, cat2);
			if (pval > pvalThreshold) return 0;

			double[][] v = ArrayUtil.trimNaNs(ArrayUtil.toDouble(cat1), ArrayUtil.toDouble(cat2));
			double corr = Correlation.pearsonVal(v[0], v[1]);
			if (correlationThreshold > 0 && Math.abs(corr) < correlationThreshold) return 0;

			return (int) Math.signum(corr);
		}

		return 0;
	}
}
