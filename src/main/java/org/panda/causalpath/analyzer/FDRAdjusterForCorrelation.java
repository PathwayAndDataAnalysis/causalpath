package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.utility.ArrayUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.FDR;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class FDRAdjusterForCorrelation
{
	Set<Set<ExperimentData>> pairs;

	CorrelationDetector cd;

	public FDRAdjusterForCorrelation(Set<Set<ExperimentData>> pairs, CorrelationDetector cd)
	{
		this.pairs = pairs;
		this.cd = cd;
	}

	public void adjustPValueThresholdsForFDR(double fdrForCorrelation)
	{
		Map<String, Double> pvals = new HashMap<>();

		for (Set<ExperimentData> pair : pairs)
		{
			if (allNumeric(pair))
			{
				Iterator<ExperimentData> iter = pair.iterator();
				NumericData nd1 = (NumericData) iter.next();
				NumericData nd2 = (NumericData) iter.next();

				double[][] v = ArrayUtil.trimNaNs(nd1.vals, nd2.vals);
				if (v[0].length >= cd.minimumSampleSize)
				{
					Tuple cor = Correlation.pearson(v[0], v[1]);
					if (!Double.isNaN(cor.p))
					{
						pvals.put(getID(pair), cor.p);
					}
				}
			}
		}

//		System.out.println("Correlation data uniformity:");
//		UniformityChecker.plot(pvals.values().stream().collect(Collectors.toList()));

		double pThr = FDR.getPValueThreshold(pvals, null, fdrForCorrelation);
		System.out.println("Correlation p-value thr = " + pThr);

		cd.setPvalThreshold(pThr);
	}

	private String getID(Set<ExperimentData> pair)
	{
		StringBuilder sb = new StringBuilder();
		pair.stream().sorted((e1, e2) -> e1.getId().compareTo(e2.getId()))
			.forEach(e -> sb.append(e.getId()).append(":"));
		return sb.toString();
	}

	private boolean allNumeric(Set<ExperimentData> pair)
	{
		return !pair.stream().anyMatch(e -> !(e instanceof NumericData));
	}
}
