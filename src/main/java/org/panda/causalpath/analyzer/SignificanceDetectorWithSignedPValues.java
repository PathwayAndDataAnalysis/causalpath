package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.utility.ArrayUtil;
import org.panda.utility.RandomizedMatrices;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.*;

import java.io.IOException;
import java.util.*;

/**
 * When signed p-values are provided in the data, this class handles them.
 */
public class SignificanceDetectorWithSignedPValues extends SignificanceDetector
{
	public SignificanceDetectorWithSignedPValues(double threshold)
	{
		super(threshold, null, null);
	}

	@Override
	public int getChangeSign(ExperimentData data)
	{
		double p = getPValue(data);
		return p <= threshold ? getChangeValue(data) > 0 ? 1 : -1 : 0;
	}

	public Tuple testDataNaiveOnSingleData(ExperimentData data)
	{
		double signedP = ((NumericData) data).vals[0];

		if (Double.isNaN(signedP)) return Tuple.NaN;

		double sign = Math.signum(signedP);
		double p = Math.abs(signedP);

		return new Tuple(sign * (p == 0D ? 100 : -Math.log(p)), p);
	}

	@Override
	public double getChangeValue(ExperimentData data)
	{
		return testData(data).v;
	}

	@Override
	public OneDataChangeDetector makeACopy()
	{
		SignificanceDetectorWithSignedPValues det = new SignificanceDetectorWithSignedPValues(threshold);
		det.setUseMissingData(useMissingData);
		det.setCategDataSufficiencyThreshold(categDataSufficiencyThreshold);
		det.setMinimumSampleSize(minimumSampleSize);
		det.rmu = rmu;
		return det;
	}
}
