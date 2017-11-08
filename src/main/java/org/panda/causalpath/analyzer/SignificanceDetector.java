package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.causalpath.data.CategoricalData;
import org.panda.causalpath.data.PhosphoProteinData;
import org.panda.utility.ArrayUtil;
import org.panda.utility.Memory;
import org.panda.utility.statistics.ChiSquare;
import org.panda.utility.statistics.TTest;

import java.util.Arrays;

/**
 * Checks the significance of separation of control and test groups.
 */
public class SignificanceDetector extends DifferenceDetector
{
	public SignificanceDetector(double threshold, boolean[] control, boolean[] test)
	{
		super(threshold, control, test);
	}

	@Override
	public int getChangeSign(ExperimentData data)
	{
		double p = getPValue(data);
		return p <= threshold ? getChangeValue(data) > 0 ? 1 : -1 : 0;
	}

	/**
	 * Does a t-test if data is numerical, and does a chi-square test if data is categorical.
	 */
	public double getPValue(ExperimentData data)
	{
		if (data instanceof NumericData)
		{
			double p = TTest.test(ArrayUtil.subset(((NumericData) data).vals, test),
				ArrayUtil.subset(((NumericData) data).vals, control)).p;

//			if (p < 0.01 && data instanceof PhosphoProteinData && !Memory.seenBefore(data.getId()))
//			{
//				System.out.println("\n" + data.id);
//				System.out.println("p = " + p);
//				System.out.println("v = " + getChangeValue(data));
//				reportMissingValueBalance((NumericData) data);
//			}

			return p;
		}
		else if (data instanceof CategoricalData)
		{
			return ChiSquare.testDependence(((CategoricalData) data).getCategories(), control, test);
		}
		throw new RuntimeException("Unhandled kind of experiment data: " + data);
	}

	@Override
	public double getChangeValue(ExperimentData data)
	{
		if (data instanceof NumericData)
		{
			return TTest.test(ArrayUtil.subset(((NumericData) data).vals, test),
				ArrayUtil.subset(((NumericData) data).vals, control)).v;
		}
		return super.getChangeValue(data);
	}

	@Override
	public OneDataChangeDetector makeACopy()
	{
		return new SignificanceDetector(threshold, control, test);
	}


	private void reportMissingValueBalance(NumericData data)
	{
		int cNan = countNan(data.vals, control);
		int tNan = countNan(data.vals, test);

		int cCnt = ArrayUtil.countValue(control, true);
		int tCnt = ArrayUtil.countValue(test, true);

		System.out.println("Control: " + cNan + "/" + cCnt + " -- " + (cNan / (double) (cCnt)));
		System.out.println("Test   : " + tNan + "/" + tCnt + " -- " + (tNan / (double) (tCnt)));
	}

	private int countNan(double[] v, boolean[] use)
	{
		int cnt = 0;
		for (int i = 0; i < v.length; i++)
		{
			if (use[i] && Double.isNaN(v[i])) cnt++;
		}
		return cnt;
	}
}
