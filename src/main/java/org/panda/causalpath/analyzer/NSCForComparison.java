package org.panda.causalpath.analyzer;

import org.panda.causalpath.network.Relation;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.FDR;

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
public class NSCForComparison extends NetworkSignificanceCalculator
{
	/**
	 * The result maps that indicate p-values for a node downstream composition.
	 */
	Map<String, Double>[] pvalMaps;


	protected double[] significanceThreshold;

	/**
	 * Constructor with required objects.
	 * @param relations set of relations to process
	 * @param cs the configured causality searcher
	 */
	public NSCForComparison(Set<Relation> relations, CausalitySearcher cs)
	{
		super(relations, cs);
	}

	@Override
	public void setFDRThreshold(double fdrThr)
	{
		significanceThreshold = new double[pvalMaps.length];
		for (int i = 0; i < pvalMaps.length; i++)
		{
			significanceThreshold[i] = FDR.getPValueThreshold(pvalMaps[i], null, fdrThr);
		}
	}

	/**
	 * Performs a randomization experiment.
	 */
	public void run(int iterations)
	{
		// Replace relations with a shuffle-safe copy
		DataLabelShuffler dls = new DataLabelShuffler(relations);
		Set<Relation> rels = dls.getRelations();

		rels = rels.stream().filter(cs::hasConsiderableDownstreamData).collect(Collectors.toSet());

		// Init counter
		DownstreamCounterForComparison dc = new DownstreamCounterForComparison(cs, rels);

		// Get the genes with no sufficient data or their downstream have no such data even to be considered
		Set<String> ignore = dc.getGenesWithNoPotential(rels);

		// Get current statistics
		Map<String, Integer>[] current = dc.run();
		rels.stream().map(r -> r.source).distinct().filter(gene -> !ignore.contains(gene) && !current[0].containsKey(gene))
			.forEach(gene ->
			{
				current[0].put(gene, 0);
				current[1].put(gene, 0);
				current[2].put(gene, 0);
			});

		// Get max possible counts for each source gene
		Map<String, Integer> maxPotential = dc.getGenesPotentialDownstreamMax(rels);

		if (minimumPotentialTargetsToConsider > 1)
		{
			new HashSet<>(current[0].keySet()).stream().filter(gene ->
				!maxPotential.containsKey(gene) || maxPotential.get(gene) < minimumPotentialTargetsToConsider)
				.forEach(gene ->
				{
					current[0].remove(gene);
					current[1].remove(gene);
					current[2].remove(gene);
					ignore.add(gene);
				});
		}

		// Get a run with non-randomized data to find current size
		Set<Relation> result = cs.run(rels);
		long sizeCurrent = result.size();

		long sizeCnt = 0;

		// Init counters for the randomizations
		Map<String, Integer>[] cnt = new Map[3];
		for (int i = 0; i < 3; i++)
		{
			cnt[i] = new HashMap<>();
		}

		Progress prog = new Progress(iterations, "Calculating significances");

		for (int i = 0; i < iterations; i++)
		{
			// Shuffle data labels and count downstream of each gene
			dls.shuffle();
			Map<String, Integer>[] run = dc.run();

			// Count the cases shuffling provided as good results
			for (int j = 0; j < 3; j++)
			{
				for (String gene : run[j].keySet())
				{
					if (ignore.contains(gene))
					{
						continue;
					}

					if (!cnt[j].containsKey(gene)) cnt[j].put(gene, 0);
//					if (!current[j].containsKey(gene)) current[j].put(gene, 0);

					if (run[j].get(gene) >= current[j].get(gene)) cnt[j].put(gene, cnt[j].get(gene) + 1);
				}
			}

			// Generate a result network for the randomized data
			result = cs.run(rels);

			// Note if the network is as big
			if (result.size() >= sizeCurrent) sizeCnt++;

			prog.tick();
		}

		// Convert counts to p-values

		graphSizePval = sizeCnt / (double) iterations;

		this.pvalMaps = new Map[3];
		for (int i = 0; i < 3; i++)
		{
			pvalMaps[i] = new HashMap<>();

			for (String gene : current[i].keySet())
			{
				double pval;

				if (current[i].get(gene) == 0)
				{
					pval = 1;
				}
				else
				{
					int c = !cnt[i].containsKey(gene) ? 0 : cnt[i].get(gene);
					if (c == 0) c++; // we don't want 0 as a p-value because they will always pass multiple hypothesis correction
					pval = c / (double) iterations;
				}

				pvalMaps[i].put(gene, pval);
			}
		}
	}

	/**
	 * Gets the genes in a result set.
	 */
	private Set<String> getGenes(Set<Relation> result)
	{
		return result.stream().map(r -> new String[]{r.source, r.target}).flatMap(Arrays::stream)
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
	public double getDownstreamActivitySignificanceThreshold()
	{
		return significanceThreshold[0];
	}

	public double getActivatorySignificanceThreshold()
	{
		return significanceThreshold[1];
	}

	public double getInhibitorySignificanceThreshold()
	{
		return significanceThreshold[2];
	}

	public boolean isActivatingTargetsSignificant(String gene)
	{
		return getActivatoryPvals().containsKey(gene) &&
			getActivatoryPvals().get(gene) <= getActivatorySignificanceThreshold();
	}

	public boolean isInhibitoryTargetsSignificant(String gene)
	{
		return getInhibitoryPvals().containsKey(gene) &&
			getInhibitoryPvals().get(gene) <= getInhibitorySignificanceThreshold();
	}

	/**
	 * Gets the minimum of 3 pvals for the gene.
	 * @param gene the gene
	 * @return minimum pval
	 */
	private Double getMinimumPval(String gene)
	{
		double p = 1;
		if (pvalMaps[0].containsKey(gene) && pvalMaps[0].get(gene) < p) p = pvalMaps[0].get(gene);
		if (pvalMaps[1].containsKey(gene) && pvalMaps[1].get(gene) < p) p = pvalMaps[1].get(gene);
		if (pvalMaps[2].containsKey(gene) && pvalMaps[2].get(gene) < p) p = pvalMaps[2].get(gene);
		return p;
	}

	/**
	 * Returns significant genes using the current p-value threshold.
	 * @return significant genes
	 */
	public Map<String, Integer> getSignificantGenes()
	{
		Map<String, Integer> sig = new HashMap<>();
		for (String gene : pvalMaps[0].keySet())
		{
			Integer i = null;
			if (getActivatoryPvals().get(gene) <= getActivatorySignificanceThreshold()) i = 1;
			if (getInhibitoryPvals().get(gene) <= getInhibitorySignificanceThreshold())
			{
				i = i == null ? -1 : 0;
			}

			if (i != null) sig.put(gene, i);
		}
		return sig;
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

	public void loadFromFile(String filename) throws IOException
	{
		String line = Files.lines(Paths.get(filename)).findFirst().get();
		this.graphSizePval = Double.valueOf(line.substring(line.lastIndexOf(" ") + 1));

		this.pvalMaps = new Map[3];
		for (int i = 0; i < pvalMaps.length; i++)
		{
			pvalMaps[i] = new HashMap<>();
		}

		Files.lines(Paths.get(filename)).skip(2).map(l -> l.split("\t")).forEach(t ->
		{
			pvalMaps[0].put(t[0], Double.valueOf(t[1]));
			pvalMaps[1].put(t[0], Double.valueOf(t[2]));
			pvalMaps[2].put(t[0], Double.valueOf(t[3]));
		});
	}
}
