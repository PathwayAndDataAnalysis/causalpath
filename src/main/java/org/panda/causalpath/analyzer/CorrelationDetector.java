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
 * Change detection for experiment data pairs based on their correlation.
 *
 * Created by babur on 3/24/16.
 */
public class CorrelationDetector implements TwoDataChangeDetector
{
	protected double correlationThreshold;
	protected double pvalThreshold;

	/**
	 * If no threshold is needed, pass -1 for that threshold.
	 * @param correlationThreshold
	 * @param pvalThreshold
	 */
	public CorrelationDetector(double correlationThreshold, double pvalThreshold)
	{
		this.correlationThreshold = correlationThreshold;
		this.pvalThreshold = pvalThreshold;
	}

	public void setCorrelationThreshold(double correlationThreshold)
	{
		this.correlationThreshold = correlationThreshold;
	}

	public void setPvalThreshold(double pvalThreshold)
	{
		this.pvalThreshold = pvalThreshold;
	}

	@Override
	public int getChangeSign(ExperimentData data1, ExperimentData data2)
	{
		if (data1 instanceof NumericData && data2 instanceof NumericData)
		{
			NumericData nd1 = (NumericData) data1;
			NumericData nd2 = (NumericData) data2;

			double[][] v = ArrayUtil.trimNaNs(nd1.vals, nd2.vals);
			Tuple corr = Correlation.pearson(v[0], v[1]);

			if (pvalThreshold > 0)
			{
				if (corr.p > pvalThreshold) return 0;
			}

			// it passes p-val thr at this point, if such a threshold exists.

			if (correlationThreshold > 0 && Math.abs(corr.v) < correlationThreshold) return 0;
			return (int) Math.signum(corr.v);
		}
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
