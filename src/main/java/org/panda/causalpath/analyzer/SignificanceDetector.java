package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.utility.ArrayUtil;
import org.panda.utility.RandomizedMatrices;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.*;

import java.io.IOException;
import java.util.*;

/**
 * Checks the significance of separation of control and test groups.
 */
public class SignificanceDetector extends DifferenceDetector
{
	protected int minimumSampleSize;
	protected boolean useMissingData;
	protected double categDataSufficiencyThreshold;

	protected boolean paired;

	RandomMatrixUser rmu;

	public SignificanceDetector(double threshold, boolean[] control, boolean[] test)
	{
		super(threshold, control, test);
		useMissingData = false;
		minimumSampleSize = 0;
		paired = false;
	}

	public void setPaired(boolean paired)
	{
		this.paired = paired;
	}

	public void setUseMissingData(boolean useMissingData)
	{
		this.useMissingData = useMissingData;
	}

	public void setCategDataSufficiencyThreshold(double categDataSufficiencyThreshold)
	{
		this.categDataSufficiencyThreshold = categDataSufficiencyThreshold;
	}

	public void setMinimumSampleSize(int minimumSampleSize)
	{
		this.minimumSampleSize = minimumSampleSize;
	}

	public void setRadomizedMatrices(RandomizedMatrices phosphoRM, RandomizedMatrices totProtRM,
		Set<String> testValueColumn) throws IOException, ClassNotFoundException
	{
		this.rmu = new RandomMatrixUser(phosphoRM, totProtRM, testValueColumn);
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
		return testData(data).p;
	}

	public Tuple testData(ExperimentData data)
	{
		Tuple tup1 = testDataNaive(data);

		if (useMissingData)
		{
			Tuple tup2 = testDataMissing(data);
			tup1 = tup1.getCombined(tup2);
		}

		return tup1;
	}

	/**
	 * Compares controls and tests. Considers independent repeats and merges their p-values using Fisher's combined
	 * probability.
	 *
	 * @param data data to test
	 * @return significance
	 */
	public Tuple testDataNaive(ExperimentData data)
	{
		if (data.hasRepeatData())
		{
			List<Tuple> list = new ArrayList<>();
			Tuple t = testDataNaiveOnSingleData(data);
			if (!t.isNaN()) list.add(t);

			for (ExperimentData repData : data.getRepeatData())
			{
				t = testDataNaiveOnSingleData(repData);
				if (!t.isNaN()) list.add(t);
			}

			if (list.isEmpty()) return new Tuple();
			if (list.size() == 1) return list.get(0);

			double[] p = new double[list.size()];
			int[] s = new int[list.size()];
			double[] v = new double[list.size()];

			for (int i = 0; i < list.size(); i++)
			{
				t = list.get(i);
				p[i] = t.p;
				v[i] = t.v;
				s[i] = Math.signum(v[i]) == Math.signum(v[0]) ? -1 : 1;
			}

			// Get the combined signed p-value
			Tuple sP = StouffersCombinedProbability.combineP2Tailed(p, s);

			double pC = sP.p;
			double vC = Summary.mean(v);
			if (Math.signum(s[0]) != Math.signum(sP.v) && Math.signum(v[0]) == Math.signum(vC))
			{
				vC = -vC;
			}

			return new Tuple(vC, pC);
		}
		else return testDataNaiveOnSingleData(data);
	}

	public Tuple testDataNaiveOnSingleData(ExperimentData data)
	{
		if (data instanceof NumericData)
		{
//			if (true) return new Tuple();

			double[] vals = ((NumericData) data).vals;

			double[] testVals = paired ? ArrayUtil.subset(vals, test, control) : ArrayUtil.subset(vals, test);
			double[] ctrlVals = paired ? ArrayUtil.subset(vals, control, test) : ArrayUtil.subset(vals, control);

			double p = Double.NaN;
			double sign = Double.NaN;

			if (testVals.length >= minimumSampleSize && ctrlVals.length >= minimumSampleSize)
			{
				Tuple ttest = paired ? TTest.testPaired(ctrlVals, testVals) : TTest.test(ctrlVals, testVals);
				p = ttest.p;
				sign = Math.signum(ttest.v);
			}

			return new Tuple(sign * (p == 0 ? 100 : -Math.log(p)), p);
		}
		else if (data instanceof CategoricalData)
		{
			int[] cat = ((CategoricalData) data).getCategories();

			long[][] c = ArrayUtil.convertCategorySubsetsToContingencyTables(cat, control, test);

			if (c.length < 2) return new Tuple();

			double p;

			if (c.length == 2)
			{
				if (c[0].length == 2)
				{
					long[][] e = ArrayUtil.getExtremized(c);
					if (GTest.testDependence(e) > categDataSufficiencyThreshold)
					{
						return new Tuple();
					}
				}

				p = GTest.testDependence(c);
			}
			else
			{
				p = ChiSquare.testDependence(cat, control, test);
			}

			double meanC = ArrayUtil.mean(cat, control);
			double meanT = ArrayUtil.mean(cat, test);
			int sign = (int) Math.signum(meanT - meanC);
			return new Tuple(sign * (p == 0 ? 100 : -Math.log(p)), p);
		}
		throw new RuntimeException("Unhandled kind of experiment data: " + data);
	}

	public Tuple testDataMissing(ExperimentData data)
	{
		if (data instanceof ProteinData)
		{
			data = ((ProteinData) data).getPresenceData();

			if (rmu != null)
			{
				Tuple t = testDataNaive(data);
				if (!t.isNaN())
				{
					t.p = rmu.getPval((PresenceData) data);
				}
				return t;
			}
			else
			{
				return testDataNaive(data);
			}
		}
		return new Tuple();
	}

	@Override
	public double getChangeValue(ExperimentData data)
	{
		return testData(data).v;
	}

	@Override
	public OneDataChangeDetector makeACopy()
	{
		SignificanceDetector det = new SignificanceDetector(threshold, control, test);
		det.setUseMissingData(useMissingData);
		det.setCategDataSufficiencyThreshold(categDataSufficiencyThreshold);
		det.setMinimumSampleSize(minimumSampleSize);
		det.rmu = rmu;
		return det;
	}

	protected class RandomMatrixUser
	{
		RandomizedMatrices totProtRandM;
		RandomizedMatrices phosphoRandM;
		boolean[] totProtTestMark;
		boolean[] phosphoTestMark;

		Map<String, Double> cache;

		RandomMatrixUser(RandomizedMatrices phosphoRM, RandomizedMatrices totProtRM, Set<String> testValueColumn)
			throws IOException, ClassNotFoundException
		{
			phosphoRandM = phosphoRM;
			totProtRandM = totProtRM;
			totProtTestMark = totProtRandM.convertSamplesToBooleanArray(testValueColumn);
			phosphoTestMark = phosphoRandM.convertSamplesToBooleanArray(testValueColumn);
			cache = new HashMap<>();
		}

		double getPval(PresenceData data)
		{
			if (cache.containsKey(data.getId())) return cache.get(data.getId());

			double pval;
			int observed = ArrayUtil.countValue(data.getCategories(), 0, test);
			if (data.getType() == DataType.PHOSPHOPROTEIN)
			{
				pval = phosphoRandM.getPValueForObservedAmount(data.getId(), observed, phosphoTestMark);
			}
			else
			{
				pval = totProtRandM.getPValueForObservedAmount(data.getId(), observed, totProtTestMark);
			}

			cache.put(data.getId(), pval);
			return pval;
		}
	}
}
