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
 * Calculates the significance of several things in the result network. These are the size of the network overall, and
 * the significance of downstream correlation frequency at the downstream of proteins.
 *
 * @author Ozgun Babur
 */
public class NSCForCorrelation extends NetworkSignificanceCalculator
{
	/**
	 * The result maps that indicate p-values for a node downstream composition.
	 */
	Map<String, Double> pvals;

	/**
	 * Constructor with the network.
	 */
	public NSCForCorrelation(Set<Relation> relations)
	{
		super(relations);
	}

	public NSCForCorrelation(Set<Relation> relations, boolean forceSiteMatching, int siteProximityThreshold,
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
		DownstreamCounter dc = new DownstreamCounterForCorrelation(rtcc, causal);
		Map<String, Integer> current = dc.run(relations)[0];

		// Get a run with non-randomized data to find current size
		CausalitySearcher cs = new CausalitySearcher(rtcc);
		cs.setCausal(causal);
		Set<RelationAndSelectedData> result = cs.run(relations);
		if (graphFilter != null) result = graphFilter.filter(result);
		int sizeCurrent = (int) result.stream().map(r -> r.relation).distinct().count();

		int sizeCnt = 0;

		// Init counters for the randomizations
		Map<String, Integer> cnt = new HashMap<>();

		// Replace relations with a shuffle-safe copy
		DataLabelShufflerForCorrelation dls = new DataLabelShufflerForCorrelation(relations);
		Set<Relation> rels = dls.getRelations();

		Progress prog = new Progress(iterations, "Calculating significances");

		for (int i = 0; i < iterations; i++)
		{
			// Shuffle data labels and count downstream of each gene
			dls.shuffle();
			Map<String, Integer> run = dc.run(rels)[0];

			// Count the cases shuffling provided as good results
			for (String gene : run.keySet())
			{
				if (!cnt.containsKey(gene)) cnt.put(gene, 0);
				if (!current.containsKey(gene)) current.put(gene, 0);

				if (run.get(gene) >= current.get(gene)) cnt.put(gene, cnt.get(gene) + 1);
			}

			// Generate a result network for the randomized data
			result = cs.run(rels);
			if (graphFilter != null) result = graphFilter.filter(result);

			// Note if the network is as big
			if (result.stream().map(r -> r.relation).distinct().count() >= sizeCurrent) sizeCnt++;

			prog.tick();
		}

		// Convert counts to p-values

		graphSizePval = sizeCnt / (double) iterations;

		this.pvals = new HashMap<>();

		for (String gene : current.keySet())
		{
			pvals.put(gene, current.get(gene) == 0 ? 1 : !cnt.containsKey(gene) ? 0 :
				cnt.get(gene) / (double) iterations);
		}
	}

	public Map<String, Double> getDownstreamActivityPvals()
	{
		return pvals;
	}

	public void writeResults(String filename) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));

		writer.write("Overall graph size pval = " + getOverallGraphSizePval());
		writer.write("\nGene\tDownstream crowded pval\tDownstream suggests activation pval\tDownstream suggests inhibition pval");

		pvals.keySet().stream().sorted((g1, g2) -> pvals.get(g1).compareTo(pvals.get(g2))).forEach(gene ->
			FileUtil.lnwrite(gene + "\t" + pvals.get(gene), writer));

		writer.close();
	}
}
