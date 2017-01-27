package org.panda.causalpath.analyzer;

import org.panda.causalpath.network.GraphFilter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationAndSelectedData;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

/**
 * Calculates the significance of several things in the result network. These are the size of the network overall, the
 * significance of downstream change frequency at the downstream of proteins, and the significance of downstream
 * evidence for activation/inactivation of the protein.
 *
 * @author Ozgun Babur
 */
public class NetworkSignificanceCalculator
{
	Set<Relation> relations;
	RelationTargetCompatibilityChecker rtcc;

	Map<String, Double>[] pvalMaps;

	double graphSizePval;

	boolean causal;

	double significanceThreshold;

	GraphFilter graphFilter;

	public NetworkSignificanceCalculator(Set<Relation> relations)
	{
		this.relations = relations;
		this.rtcc = new RelationTargetCompatibilityChecker();
	}

	public NetworkSignificanceCalculator(Set<Relation> relations, boolean forceSiteMatching, int siteProximityThreshold,
		boolean causal, GraphFilter graphFilter)
	{
		this(relations);
		rtcc.setForceSiteMatching(forceSiteMatching);
		rtcc.setSiteProximityThreshold(siteProximityThreshold);
		this.causal = causal;
//		this.graphFilter = graphFilter;
	}

	public void setSignificanceThreshold(double significanceThreshold)
	{
		this.significanceThreshold = significanceThreshold;
	}

	public void run(int iterations)
	{
		DownstreamCounter dc = new DownstreamCounter(rtcc);
		Map<String, Integer>[] current = dc.run(relations);

		CausalitySearcher cs = new CausalitySearcher(rtcc);
		cs.setCausal(causal);
		Set<RelationAndSelectedData> result = cs.run(relations);
		if (graphFilter != null) result = graphFilter.filter(result);
		int sizeCurrent = (int) result.stream().map(r -> r.relation).distinct().count();

		int sizeCnt = 0;

		Map<String, Integer>[] cnt = new Map[3];
		for (int i = 0; i < 3; i++)
		{
			cnt[i] = new HashMap<>();
		}

		DataLabelShuffler dls = new DataLabelShuffler(relations);
		Set<Relation> rels = dls.getRelations();

		Progress prog = new Progress(iterations, "Calculating significances");

		for (int i = 0; i < iterations; i++)
		{
			dls.shuffle();
			Map<String, Integer>[] run = dc.run(rels);

			for (int j = 0; j < 3; j++)
			{
				for (String gene : run[j].keySet())
				{
					if (!cnt[j].containsKey(gene)) cnt[j].put(gene, 0);
					if (!current[j].containsKey(gene)) current[j].put(gene, 0);

					if (run[j].get(gene) >= current[j].get(gene)) cnt[j].put(gene, cnt[j].get(gene) + 1);
				}
			}

			result = cs.run(rels);
			if (graphFilter != null) result = graphFilter.filter(result);
			if (result.stream().map(r -> r.relation).distinct().count() >= sizeCurrent) sizeCnt++;

			prog.tick();
		}

		graphSizePval = sizeCnt / (double) iterations;

		this.pvalMaps = new Map[3];
		for (int i = 0; i < 3; i++)
		{
			pvalMaps[i] = new HashMap<>();

			for (String gene : current[i].keySet())
			{
				pvalMaps[i].put(gene, current[i].get(gene) == 0 ? 1 : !cnt[i].containsKey(gene) ? 0 :
					cnt[i].get(gene) / (double) iterations);
			}
		}
	}

	public double getOverallGraphSizePval()
	{
		return graphSizePval;
	}

	public Map<String, Double> getDownstreamActivityPvals()
	{
		return pvalMaps[0];
	}
	public Map<String, Double> getActivatoryPvals()
	{
		return pvalMaps[1];
	}
	public Map<String, Double> getInhibitoryPvals()
	{
		return pvalMaps[2];
	}

	public boolean isDownstreamSignificant(String gene)
	{
		return getDownstreamActivityPvals().containsKey(gene) && getDownstreamActivityPvals().get(gene) <= significanceThreshold;
	}
	public boolean isActivatingTargetsSignificant(String gene)
	{
		return getActivatoryPvals().containsKey(gene) && getActivatoryPvals().get(gene) <= significanceThreshold;
	}
	public boolean isInhibitoryTargetsSignificant(String gene)
	{
		return getInhibitoryPvals().containsKey(gene) && getInhibitoryPvals().get(gene) <= significanceThreshold;
	}

	private Double getMinimumPval(String gene)
	{
		double p = 1;
		if (pvalMaps[0].containsKey(gene) && pvalMaps[0].get(gene) < p) p = pvalMaps[0].get(gene);
		if (pvalMaps[1].containsKey(gene) && pvalMaps[1].get(gene) < p) p = pvalMaps[1].get(gene);
		if (pvalMaps[2].containsKey(gene) && pvalMaps[2].get(gene) < p) p = pvalMaps[2].get(gene);
		return p;
	}

	public void writeResults(String filename) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));

		writer.write("Overall graph size pval = " + getOverallGraphSizePval());
		writer.write("\nGene\tDownstream crowded pval\tDownstream suggests activation pval\tDownstream suggests inhibition pval");

		Stream.concat(pvalMaps[0].keySet().stream(), Stream.concat(pvalMaps[1].keySet().stream(), pvalMaps[2].keySet().stream()))
			.distinct().sorted((g1, g2) -> getMinimumPval(g1).compareTo(getMinimumPval(g2))).forEach(gene ->
			FileUtil.lnwrite(ArrayUtil.getString("\t", gene,
				pvalMaps[0].containsKey(gene) ? pvalMaps[0].get(gene) : 1,
				pvalMaps[1].containsKey(gene) ? pvalMaps[1].get(gene) : 1,
				pvalMaps[2].containsKey(gene) ? pvalMaps[2].get(gene) : 1), writer));

		writer.close();
	}
}
