package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.DataType;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.causalpath.network.Relation;
import org.panda.utility.CollectionUtil;
import org.panda.utility.Progress;

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

	boolean correlationBased;
	double fdrThresholdForCorrelation;
	CorrelationDetector corrDet;

	Set<ExperimentData> data;
	Set<Set<ExperimentData>> pairs;

	public RobustnessAnalysis(CausalitySearcher cs, Set<Relation> relations, Map<DataType, Double> fdrThr,
		double noiseStDev, boolean correlationBased, double fdrThresholdForCorrelation, CorrelationDetector corrDet)
	{
		this.cs = cs;
		this.relations = relations;
		this.fdrThr = fdrThr;
		this.noiseStDev = noiseStDev;
		this.correlationBased = correlationBased;
		this.fdrThresholdForCorrelation = fdrThresholdForCorrelation;
		this.corrDet = corrDet;
	}

	public void run(int iterations, String outFile) throws IOException
	{
		cs.setCollectDataUsedForInference(false);
		cs.setCollectDataWithMissingEffect(false);

		Set<Relation> trueRels = getARun(relations);

		int[] trueCnt = new int[iterations];
		int[] falseCnt = new int[iterations];

		Progress p = new Progress(iterations, "Robustness analysis");
		for (int i = 0; i < iterations; i++)
		{
			addNoise(relations);

			Set<Relation> newRels = getARun(relations);

			int overlap = CollectionUtil.countOverlap(newRels, trueRels);

			trueCnt[i] = overlap;
			falseCnt[i] = newRels.size() - overlap;

			p.tick();
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

	private Set<Relation> getARun(Set<Relation> relations) throws IOException
	{
		if (data == null && pairs == null)
		{
			if (!correlationBased)
			{
				Stream.concat(relations.stream().map(r -> r.sourceData.getDataStream().map(ExperimentData::getChDet)),
					relations.stream().map(r -> r.targetData.getDataStream().map(ExperimentData::getChDet)))
					.filter(det -> det instanceof ThresholdDetector).forEach(det -> ((ThresholdDetector) det).setThreshold(1));
			} else
			{
				corrDet.setPvalThreshold(1);
			}

			data = new HashSet<>();
			pairs = new HashSet<>();
			cs.setCollectDataUsedForInference(true);
			cs.setCausal(true);
			cs.run(relations);

			if (correlationBased) pairs.addAll(cs.getPairsUsedForInference());
			else data.addAll(cs.getDataUsedForInference());

			cs.setCausal(false);
			cs.run(relations);

			if (correlationBased) pairs.addAll(cs.getPairsUsedForInference());
			else data.addAll(cs.getDataUsedForInference());

			cs.setCollectDataUsedForInference(false);
			cs.setCausal(true);
		}

		if (correlationBased)
		{
			FDRAdjusterForCorrelation fad = new FDRAdjusterForCorrelation(null, pairs, corrDet);
			fad.adjustPValueThresholdsForFDR(fdrThresholdForCorrelation);
		}
		else
		{
			FDRAdjuster adjuster = new FDRAdjuster(null, false);
			adjuster.adjustPValueThresholdsOfDatas(data, fdrThr);
		}

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
