package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.DataType;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.causalpath.network.Relation;
import org.panda.utility.CollectionUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Ozgun Babur
 */
public class RobustnessAnalysis
{
	CausalitySearcher cs;

	Set<Relation> relations;

	Map<DataType, Double> fdrThr;

	double noiseStDev;

	public RobustnessAnalysis(CausalitySearcher cs, Set<Relation> relations, Map<DataType, Double> fdrThr, double noiseStDev)
	{
		this.cs = cs;
		this.relations = relations;
		this.fdrThr = fdrThr;
		this.noiseStDev = noiseStDev;
	}

	public void run(int iterations, String outFile) throws IOException
	{
		Set<Relation> trueRels = getARun(relations);

		int[] trueCnt = new int[iterations];
		int[] falseCnt = new int[iterations];

		for (int i = 0; i < iterations; i++)
		{
			addNoise(relations);

			Set<Relation> newRels = getARun(relations);

			int overlap = CollectionUtil.countOverlap(newRels, trueRels);

			trueCnt[i] = overlap;
			falseCnt[i] = newRels.size() - overlap;
		}

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("Iterations\tTrue relations\tFalse Relations");
		writer.write("\n0\t" + trueRels.size() + "\t0");

		for (int i = 0; i < iterations; i++)
		{
			writer.write("\n" + (i + 1) + "\t" + trueCnt[i] + "\t" + falseCnt[i]);
		}

		writer.close();
	}

	private Set<Relation> getARun(Set<Relation> relations)
	{
		Stream.concat(relations.stream().map(r -> r.sourceData.getDataStream().map(ExperimentData::getChDet)),
			relations.stream().map(r -> r.targetData.getDataStream().map(ExperimentData::getChDet)))
			.filter(det -> det instanceof ThresholdDetector).forEach(det -> ((ThresholdDetector) det).setThreshold(1));

		Set<ExperimentData> data = new HashSet<>();
		cs.setCollectDataUsedForInference(true);
		cs.setCausal(true);
		cs.run(relations);
		data.addAll(cs.getDataUsedForInference());
		cs.setCausal(false);
		cs.run(relations);
		data.addAll(cs.getDataUsedForInference());

		FDRAdjuster adjuster = new FDRAdjuster(false);
		adjuster.adjustPValueThresholdsOfDatas(data, fdrThr);

		cs.setCollectDataUsedForInference(false);
		cs.setCausal(true);
		return cs.run(relations);
	}

	private void addNoise(Set<Relation> relations)
	{
		Set<ExperimentData> datas = getExperimentDatas(relations);

		Random rand = new Random();

		for (ExperimentData data : datas)
		{
			if (data instanceof NumericData)
			{
				NumericData nd = (NumericData) data;

				for (int i = 0; i < nd.vals.length; i++)
				{
					nd.vals[i] += rand.nextGaussian() * noiseStDev;
				}
			}
		}
	}

	private Set<ExperimentData> getExperimentDatas(Set<Relation> relations)
	{
		return Stream.concat(
				relations.stream().map(r -> r.sourceData.getData()).flatMap(Collection::stream),
				relations.stream().map(r -> r.targetData.getData()).flatMap(Collection::stream))
				.collect(Collectors.toSet());
	}
}
