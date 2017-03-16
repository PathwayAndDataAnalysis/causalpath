package org.panda.causalpath.analyzer;

import org.apache.poi.util.IntList;
import org.panda.causalpath.network.GraphFilter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationAndSelectedData;
import org.panda.resource.CancerGeneCensus;
import org.panda.resource.OncoKB;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Calculates the significance of several things in the result network. These are the size of the network overall, the
 * significance of downstream change frequency at the downstream of proteins, and the significance of downstream
 * evidence for activation/inactivation of the protein.
 *
 * @author Ozgun Babur
 */
public class NSCForNonCorr extends NetworkSignificanceCalculator
{
	/**
	 * The result maps that indicate p-values for a node downstream composition.
	 */
	Map<String, Double>[] pvalMaps;

	/**
	 * Constructor with the network.
	 */
	public NSCForNonCorr(Set<Relation> relations)
	{
		super(relations);
	}

	public NSCForNonCorr(Set<Relation> relations, boolean forceSiteMatching, int siteProximityThreshold,
		boolean causal, GraphFilter graphFilter)
	{
		super(relations, forceSiteMatching, siteProximityThreshold, causal, graphFilter);
	}

	/**
	 * Performs a randomization experiment.
	 */
	public void run(int iterations)
	{
		// Get current statistics
		DownstreamCounter dc = new DownstreamCounter(rtcc);
		Map<String, Integer>[] current = dc.run(relations);

		// Get the set of cancer genes in cancer gene databases.
		Set<String> cancerGenes = new HashSet<>(OncoKB.get().getAllSymbols());
		cancerGenes.addAll(CancerGeneCensus.get().getAllSymbols());

		// Get a run with non-randomized data to find current size
		CausalitySearcher cs = new CausalitySearcher(rtcc);
		cs.setCausal(causal);
		Set<RelationAndSelectedData> result = cs.run(relations);
		if (graphFilter != null) result = graphFilter.filter(result);
		long sizeCurrent = result.stream().map(r -> r.relation).distinct().count();
		long currentCGCnt = getGenes(result).stream().filter(cancerGenes::contains).count();

		long sizeCnt = 0;
		long cgCnt = 0;

		// Init counters for the randomizations
		Map<String, Integer>[] cnt = new Map[3];
		for (int i = 0; i < 3; i++)
		{
			cnt[i] = new HashMap<>();
		}

		// Replace relations with a shuffle-safe copy
		DataLabelShuffler dls = new DataLabelShuffler(relations);
		Set<Relation> rels = dls.getRelations();

		Progress prog = new Progress(iterations, "Calculating significances");

		for (int i = 0; i < iterations; i++)
		{
			// Shuffle data labels and count downstream of each gene
			dls.shuffle();
			Map<String, Integer>[] run = dc.run(rels);

			// Count the cases shuffling provided as good results
			for (int j = 0; j < 3; j++)
			{
				for (String gene : run[j].keySet())
				{
					if (!cnt[j].containsKey(gene)) cnt[j].put(gene, 0);
					if (!current[j].containsKey(gene)) current[j].put(gene, 0);

					if (run[j].get(gene) >= current[j].get(gene)) cnt[j].put(gene, cnt[j].get(gene) + 1);
				}
			}

			// Generate a result network for the randomized data
			result = cs.run(rels);
			if (graphFilter != null) result = graphFilter.filter(result);

			// Note if the network is as big
			if (result.stream().map(r -> r.relation).distinct().count() >= sizeCurrent) sizeCnt++;

			long cgCount = getGenes(result).stream().filter(cancerGenes::contains).count();
			if (cgCount >= currentCGCnt) cgCnt++;

			prog.tick();
		}

		// Convert counts to p-values

		graphSizePval = sizeCnt / (double) iterations;
		pvalForCGEnrichment = cgCnt / (double) iterations;

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

	/**
	 * Gets the genes in a result set.
	 */
	private Set<String> getGenes(Set<RelationAndSelectedData> result)
	{
		return result.stream().map(r -> new String[]{r.relation.source, r.relation.target}).flatMap(Arrays::stream)
			.collect(Collectors.toSet());
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
