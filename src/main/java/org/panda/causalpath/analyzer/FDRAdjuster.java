package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.DataType;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.ProteinData;
import org.panda.causalpath.network.Relation;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Histogram2D;
import org.panda.utility.statistics.UniformityChecker;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class FDRAdjuster
{
	/**
	 * Whether to pool proteomics and phosphoproteomics during FDR control. Typically, this will be false for mass
	 * spectrometry data because they are coming from different experiments, and it will be true for RPPA experiments,
	 * where they are coming from the same experiment.
	 */
	boolean poolProteomics;

	/**
	 * The project directory to generate output.
	 */
	String directory;

	public FDRAdjuster(String directory, boolean poolProteomics)
	{
		this.directory = directory;
		this.poolProteomics = poolProteomics;
	}

	/**
	 * Collects all the data associated with the given relations and sets their p-value threshold to control false
	 * discovery rate at the given level.
	 * @param relations relations in the analysis
	 * @param fdrMap desired FDR threshold for each data type
	 */
	public void adjustPValueThresholdsForRelations(Set<Relation> relations, Map<DataType, Double> fdrMap) throws IOException
	{
		adjustPValueThresholdsOfDatas(relations.stream().map(Relation::getAllData).flatMap(Collection::stream)
			.collect(Collectors.toSet()), fdrMap);
	}

	public void adjustPValueThresholdsOfDatas(Set<ExperimentData> datas, Map<DataType, Double> fdrMap) throws IOException
	{
		Map<DataType, Set<ExperimentData>> dataMap = new HashMap<>();

		// classify data according to type

		for (ExperimentData data : datas)
		{
			OneDataChangeDetector chDet = data.getChDet();

			if (chDet instanceof SignificanceDetector)
			{
				DataType type = data.getType();
				if (!dataMap.containsKey(type)) dataMap.put(type, new HashSet<>());
				dataMap.get(type).add(data);
			}
		}

		// unite proteomics if needed
		if (poolProteomics)
		{
			dataMap.get(DataType.PROTEIN).addAll(dataMap.get(DataType.PHOSPHOPROTEIN));
			dataMap.remove(DataType.PHOSPHOPROTEIN);
		}

		// for each type detect the p-value threshold and set it

		BufferedWriter writer = null;
		if (directory != null)
		{
			writer = Files.newBufferedWriter(Paths.get(directory + "/pval-uniformity.txt"));
			writer.write("Pool proteomics = " + poolProteomics + "\n");
		}

		for (DataType type : dataMap.keySet())
		{
			if (!fdrMap.containsKey(type)) continue;

			Map<ExperimentData, Double> pvalues = new HashMap<>();
			dataMap.get(type).forEach(d ->
			{
				double p = ((SignificanceDetector) d.getChDet()).getPValue(d);
				if (!Double.isNaN(p)) pvalues.put(d, p);
			});

			// Record uniformity

			if (writer != null)
			{
				writer.write("\nData type = " + type + "\n");
				UniformityChecker.plot(pvalues.values().stream().collect(Collectors.toList()), writer);
			}

			double pThr = FDR.getPValueThreshold(pvalues, null, fdrMap.get(type));

//			if (type == DataType.PROTEIN || type == DataType.PHOSPHOPROTEIN)
//			{
//				System.out.println(type + " data uniformity:");
//				UniformityChecker.plot(pvalues.values().stream().collect(Collectors.toList()));
//			}

			System.out.println("type = " + type + "\tpThr = " + pThr);

			dataMap.get(type).stream().map(ExperimentData::getChDet).distinct()
				.forEach(det -> ((ThresholdDetector) det).setThreshold(pThr));
		}

		if (writer != null) writer.close();

//		debug code for plotting t-test versus G-test 2D histogram ----------------------

//		System.out.println("Plotting 2D histograms");
//		for (DataType type : dataMap.keySet())
//		{
//			System.out.println("type = " + type);
//			Histogram2D h = new Histogram2D(0.05);
//			List<Double> list1 = new ArrayList<>();
//			List<Double> list2 = new ArrayList<>();
//			dataMap.get(type).forEach(d ->
//			{
//				double pN = ((SignificanceDetector) d.getChDet()).testDataNaive(d).p;
//				double pM = ((SignificanceDetector) d.getChDet()).testDataMissing(d).p;
//				if (!Double.isNaN(pN) && !Double.isNaN(pM))
//				{
//					h.count(pN, pM);
//					list1.add(pN);
//					list2.add(pM);
//				}
//			});
//			Tuple t = Correlation.pearson(ArrayUtil.toArray(list1), ArrayUtil.toArray(list2));
//			System.out.println("Pearson correlation " + t);
//			h.plot();
//		}

//		debug code ----------------------

	}
}
