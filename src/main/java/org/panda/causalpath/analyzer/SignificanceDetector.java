package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.causalpath.data.CategoricalData;
import org.panda.utility.ArrayUtil;
import org.panda.utility.statistics.ChiSquare;
import org.panda.utility.statistics.TTest;

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
		return p < threshold ? getChangeValue(data) > 0 ? 1 : -1 : 0;
	}

	/**
	 * Does a t-test if data is numerical, and does a chi-square test if data is categorical.
	 */
	private double getPValue(ExperimentData data)
	{
		if (data instanceof NumericData)
		{
			return TTest.getPValOfMeanDifference(ArrayUtil.subset(((NumericData) data).vals, control),
				ArrayUtil.subset(((NumericData) data).vals, test));
		}
		else if (data instanceof CategoricalData)
		{
			return ChiSquare.testDependence(((CategoricalData) data).getCategories(), control, test);
		}
		throw new RuntimeException("Unhandled kind of experiment data: " + data);
	}
}
