package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.utility.ArrayUtil;
import org.panda.utility.BooleanMatrixRandomizer;
import org.panda.utility.RandomizedMatrices;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.*;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;

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

	protected boolean useMissingData;
	protected double categDataSufficiencyThreshold;

	RandomMatrixUser rmu;

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

	public void setUseMissingData(boolean useMissingData)
	{
		this.useMissingData = useMissingData;
	}

	public void setCategDataSufficiencyThreshold(double categDataSufficiencyThreshold)
	{
		this.categDataSufficiencyThreshold = categDataSufficiencyThreshold;
	}

	public void setRandomMatrices(RandomizedMatrices phosphoRandM, RandomizedMatrices totProtRandM,
		List<String> valueColumn) throws IOException, ClassNotFoundException
	{
		this.rmu = new RandomMatrixUser(phosphoRandM, totProtRandM, valueColumn);
	}

	@Override
	public int getChangeSign(ExperimentData data1, ExperimentData data2)
	{
		Tuple corr = calcCorrelation(data1, data2);
		if (Double.isNaN(corr.p)) return 0;

		if (pvalThreshold >= 0 && corr.p > pvalThreshold) return 0;
		if (correlationThreshold > 0 && Math.abs(corr.v) < correlationThreshold) return 0;
		if (correlationUpperThreshold > 0 && Math.abs(corr.v) > correlationUpperThreshold) return 0;
		return (int) Math.signum(corr.v);
	}

	public Tuple calcCorrelation(ExperimentData data1, ExperimentData data2)
	{
		Tuple tup1 = calcCorrelationNaive(data1, data2);

		if (useMissingData)
		{
			Tuple tup2 = calcCorrelationWithMissingData(data1, data2);
			tup1 = tup1.getCombined(tup2);
		}

		return tup1;
	}


	public Tuple calcCorrelationNaive(ExperimentData data1, ExperimentData data2)
	{
		// Case 1: Both data are numeric

		if (data1 instanceof NumericData && data2 instanceof NumericData)
		{
//			if (true) return new Tuple();

			NumericData nd1 = (NumericData) data1;
			NumericData nd2 = (NumericData) data2;

			double[][] v = ArrayUtil.trimNaNs(nd1.vals, nd2.vals);

			if (v[0].length < minimumSampleSize) return new Tuple(Double.NaN, Double.NaN);

			return Correlation.pearson(v[0], v[1]);
		}

		// Case 2: One data is numeric, the other is categorical

		else if (data1 instanceof CategoricalData && data2 instanceof NumericData ||
			data2 instanceof CategoricalData && data1 instanceof NumericData)
		{
			CategoricalData qd = (CategoricalData) (data1 instanceof CategoricalData ? data1 : data2);
			NumericData nd = (NumericData) (data1 == qd ? data2 : data1);

			int[] cat = qd.getCategories();
			int catNum = qd.getNumberOfCategories(cat);
			if (catNum == 1) return new Tuple();

			double[] val = nd.vals;

			List<double[]> groups = ArrayUtil.separateToCategories(cat, val, minimumSampleSize);

			double pval;
			if (catNum == 2)
			{
				pval = TTest.test(groups.get(0), groups.get(1)).p;
			}
			else
			{
				pval = Anova.oneWayPval(groups);
			}

			double[][] v = ArrayUtil.trimNaNs(ArrayUtil.toDouble(cat), val);
			double corr = Correlation.pearsonVal(v[0], v[1]);
			return new Tuple(corr, pval);
		}

		// Case 3: Both data are categorical

		else if (data1 instanceof CategoricalData && data2 instanceof CategoricalData)
		{
			CategoricalData qd1 = (CategoricalData) data1;
			CategoricalData qd2 = (CategoricalData) data2;

			int[] cat1 = qd1.getCategories();
			int[] cat2 = qd2.getCategories();

			long[][] cont = ArrayUtil.convertCategoriesToContingencyTable(cat1, cat2);
			long[][] recon = GTest.reconfigure(cont);

			double pval;

			if (recon != null)
			{
				if (recon[0].length == 2)
				{
					// if an extreme distribution of the matrix cannot produce small enough p-value, skip it
					long[][] ext = ArrayUtil.getExtremized(recon);
					if (GTest.testDependence(ext) > categDataSufficiencyThreshold)
					{
						return new Tuple();
					}
				}
				pval = GTest.testDependence(recon);
			}
			else if (cont.length < 2 || (cont.length > 0 && cont[0].length == 1))
			{
				return new Tuple();
			}
			else
			{
				pval = ChiSquare.testDependence(cont);
			}

			double[][] v = ArrayUtil.trimNaNs(ArrayUtil.toDouble(cat1), ArrayUtil.toDouble(cat2));
			double corr = Correlation.pearsonVal(v[0], v[1]);
			return new Tuple(corr, pval);
		}

		return new Tuple();
	}

	public Tuple calcCorrelationWithMissingData(ExperimentData data1, ExperimentData data2)
	{
		if (data1 instanceof ProteinData || data2 instanceof ProteinData)
		{
			if (data1 instanceof ProteinData) data1 = ((ProteinData) data1).getPresenceData();
			if (data2 instanceof ProteinData) data2 = ((ProteinData) data2).getPresenceData();

			if (rmu != null)
			{
				Tuple t = calcCorrelationNaive(data1, data2);
				if (!t.isNaN())
				{
					t.p = rmu.getPval((PresenceData) data1, (PresenceData) data2);
				}
				return t;
			}
			else
			{
				return calcCorrelationNaive(data1, data2);
			}
		}
		return new Tuple();
	}

	protected class RandomMatrixUser
	{
		RandomizedMatrices totProtRandM;
		RandomizedMatrices phosphoRandM;
		boolean[] totProtUse;
		boolean[] phosphoUse;
		List<String> commonColumns;
		List<String> valueColumn;
		int[] randPhosphoToCPIndexMap;
		int[] randTotProtToCPIndexMap;

		RandomMatrixUser(RandomizedMatrices phosphoRandM, RandomizedMatrices totProtRandM, List<String> valueColumn)
			throws IOException, ClassNotFoundException
		{
			this.totProtRandM = totProtRandM;
			this.phosphoRandM = phosphoRandM;
			totProtUse = totProtRandM.convertSamplesToBooleanArray(valueColumn);
			phosphoUse = phosphoRandM.convertSamplesToBooleanArray(valueColumn);
			commonColumns = totProtRandM.getCommonColumns(phosphoRandM);
			this.valueColumn = valueColumn;
			randPhosphoToCPIndexMap = phosphoRandM.getThisToListIndexMapping(valueColumn);
			randTotProtToCPIndexMap = totProtRandM.getThisToListIndexMapping(valueColumn);
		}

		double getPval(PresenceData data1, PresenceData data2)
		{
			int[] cat1 = data1.getCategories();
			int[] cat2 = data2.getCategories();
			int overlap = 0;
			for (int i = 0; i < cat1.length; i++)
			{
				if (cat1[i] == cat2[i] && cat1[i] != ArrayUtil.ABSENT_INT) overlap++;
			}

			double pval;
			if (data1.getType() == DataType.PHOSPHOPROTEIN)
			{
				if (data2.getType() == DataType.PHOSPHOPROTEIN)
				{
					pval = phosphoRandM.getPValueForOverlapSignificance(data1.getId(), data2.getId(), overlap, phosphoUse);
				} else
				{
					pval = phosphoRandM.getPValueForOverlapSignificance(data1.getId(), totProtRandM, data2.getId(), overlap, commonColumns);
				}
			} else
			{
				if (data2.getType() == DataType.PHOSPHOPROTEIN)
				{
					pval = totProtRandM.getPValueForOverlapSignificance(data1.getId(), phosphoRandM, data2.getId(), overlap, commonColumns);
				} else
				{
					pval = totProtRandM.getPValueForOverlapSignificance(data1.getId(), data2.getId(), overlap, totProtUse);
				}
			}

			return pval;
		}

		double getPval(PresenceData presData, ProteinData protData)
		{
			Tuple obs = calcCorrelationNaive(presData, protData);
			if (obs.isNaN()) return Double.NaN;

			if (presData.getType() == DataType.PHOSPHOPROTEIN)
			{
				return phosphoRandM.getPValueForSignificantT(
					presData.getId(), protData.vals, obs.p, randPhosphoToCPIndexMap, minimumSampleSize);
			}
			else
			{
				return totProtRandM.getPValueForSignificantT(
					presData.getId(), protData.vals, obs.p, randTotProtToCPIndexMap, minimumSampleSize);
			}
		}
	}
}