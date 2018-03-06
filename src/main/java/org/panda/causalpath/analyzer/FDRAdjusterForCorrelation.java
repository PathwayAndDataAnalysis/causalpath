package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.utility.ArrayUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class FDRAdjusterForCorrelation
{
	Set<List<ExperimentData>> pairs;

	CorrelationDetector cd;

	String directory;

	public FDRAdjusterForCorrelation(String directory, Set<List<ExperimentData>> pairs, CorrelationDetector cd)
	{
		this.directory = directory;
		this.pairs = pairs;
		this.cd = cd;
	}

	public void adjustPValueThresholdsForFDR(double fdrForCorrelation) throws IOException
	{
		Map<String, Double> pvals = new HashMap<>();

		for (List<ExperimentData> pair : pairs)
		{
			Iterator<ExperimentData> iter = pair.iterator();
			ExperimentData data1 = iter.next();
			ExperimentData data2 = iter.next();

			Tuple corr = cd.calcCorrelation(data1, data2);
			if (!corr.isNaN())
			{
				pvals.put(getID(pair), corr.p);
			}
		}

		if (directory != null)
		{
			BufferedWriter writer = Files.newBufferedWriter(Paths.get(directory + "/pval-uniformity.txt"));
			UniformityChecker.plot(pvals.values().stream().collect(Collectors.toList()), writer);
			writer.close();
		}

		double pThr = FDR.getPValueThreshold(pvals, null, fdrForCorrelation);
		System.out.println("Correlation p-value thr = " + pThr);

		cd.setPvalThreshold(pThr);

		// debug code for plotting G-test vs t-test 2D histogram----------------
//		Histogram2D h = new Histogram2D(0.005);
//		Histogram hc = new Histogram(0.005);
//		List<Double> list1 = new ArrayList<>();
//		List<Double> list2 = new ArrayList<>();
//		for (Set<ExperimentData> pair : pairs)
//		{
//			Iterator<ExperimentData> iter = pair.iterator();
//			ExperimentData data1 = iter.next();
//			ExperimentData data2 = iter.next();
//
//			Tuple cN = cd.calcCorrelationNaive(data1, data2);
//			Tuple cM = cd.calcCorrelationWithMissingData(data1, data2);
//			if (!cN.isNaN() && !cM.isNaN())// &&
////				cN.p < 0.1 && cM.p < 0.1)
//			{
//				h.count(cN.p, cM.p);
//				list1.add(cN.p);
//				list2.add(cM.p);
//				hc.count(cM.v);
//			}
//		}
//		System.out.println("Pearson correlation " +
//			Correlation.pearson(ArrayUtil.toArray(list1), ArrayUtil.toArray(list2)));
//		System.out.println("Correlation histogram");
//		hc.print();
//		h.plot();
		// end of debug-----------------------------------------------------------------
	}

	private String getID(List<ExperimentData> pair)
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
