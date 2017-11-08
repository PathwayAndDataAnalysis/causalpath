package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.DataType;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.network.Relation;
import org.panda.utility.statistics.FDR;

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

	public FDRAdjuster(boolean poolProteomics)
	{
		this.poolProteomics = poolProteomics;
	}

	public void adjustPValueThresholdsForRelations(Set<Relation> relations, Map<DataType, Double> fdrMap)
	{
		adjustPValueThresholdsOfDatas(relations.stream().map(Relation::getAllData).flatMap(Collection::stream)
			.collect(Collectors.toSet()), fdrMap);
	}

	public void adjustPValueThresholdsOfDatas(Set<ExperimentData> datas, Map<DataType, Double> fdrMap)
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

		for (DataType type : dataMap.keySet())
		{
			if (!fdrMap.containsKey(type)) continue;

			Map<ExperimentData, Double> pvalues = new HashMap<>();
			dataMap.get(type).forEach(d ->
			{
				double p = ((SignificanceDetector) d.getChDet()).getPValue(d);
				if (!Double.isNaN(p)) pvalues.put(d, p);
			});
			double pThr = FDR.getPValueThreshold(pvalues, null, fdrMap.get(type));

//			if (clz == ProteinData.class)
//			{
//				UniformityChecker.plot(new ArrayList<>(pvalues.values()));
//				System.exit(0);
//			}

			System.out.println("type = " + type + "\tpThr = " + pThr);

			dataMap.get(type).stream().map(ExperimentData::getChDet).distinct()
				.forEach(det -> ((ThresholdDetector) det).setThreshold(pThr));
		}
	}
}
