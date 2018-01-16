package org.panda.causalpath.analyzer;

import org.panda.causalpath.network.Relation;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.KernelDensityPlot;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

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
	public NSCForCorrelation(Set<Relation> relations, CausalitySearcher cs)
	{
		super(relations, cs);
	}

	@Override
	public void setFDRThreshold(double fdrThr)
	{
		setPvalThreshold(FDR.getPValueThreshold(pvals, null, fdrThr));
	}

	/**
	 * Performs a randomization experiment.
	 */
	public void run(int iterations)
	{
		// Init counter
		DownstreamCounter dc = new DownstreamCounterForCorrelation(cs);

		Set<String> ignore = dc.getGenesWithNoPotential(relations);

		// Get current statistics
		Map<String, Integer> current = dc.run(relations)[0];
		relations.stream().map(r -> r.source).filter(gene -> !ignore.contains(gene) && !current.containsKey(gene))
			.forEach(gene -> current.put(gene, 0));

		// Get max possible counts for each source gene
		Map<String, Integer> maxPotential = dc.getGenesPotentialDownstreamMax(relations);

		if (minimumPotentialTargetsToConsider > 1)
		{
			new HashSet<>(current.keySet()).stream().filter(gene ->
				!maxPotential.containsKey(gene) || maxPotential.get(gene) < minimumPotentialTargetsToConsider)
				.forEach(current::remove);
		}

		// Get a run with non-randomized data to find current size
		Set<Relation> result = cs.run(relations);
		int sizeCurrent = result.size();

		int sizeCnt = 0;

		// Init counters for the randomizations
		Map<String, Integer> cnt = new HashMap<>();

		// Replace relations with a shuffle-safe copy
		DataLabelShufflerForCorrelation dls = new DataLabelShufflerForCorrelation(relations);
		Set<Relation> rels = dls.getRelations();

		Progress prog = new Progress(iterations, "Calculating significances");
		double[] sizes = new double[iterations];

		for (int i = 0; i < iterations; i++)
		{
			System.gc();

			// Shuffle data labels and count downstream of each gene
			dls.shuffle();
			Map<String, Integer> run = dc.run(rels)[0];

			assert !run.keySet().stream().filter(ignore::contains).filter(gene -> run.get(gene) > 0).findAny()
				.isPresent() : "Ignored gene contains non-zero downstream!";

			// Count the cases shuffling provided as good results
			for (String gene : run.keySet())
			{
				if (minimumPotentialTargetsToConsider > 1 && maxPotential.get(gene) < minimumPotentialTargetsToConsider)
				{
					continue;
				}

				assert !ignore.contains(gene);

				if (!cnt.containsKey(gene)) cnt.put(gene, 0);
				if (!current.containsKey(gene)) current.put(gene, 0);

				if (run.get(gene) >= current.get(gene)) cnt.put(gene, cnt.get(gene) + 1);
			}

			// Generate a result network for the randomized data
			result = cs.run(rels);

			// Note if the network is as big
			if (result.size() >= sizeCurrent) sizeCnt++;

			sizes[i] = result.size();

			prog.tick();
		}

		// Plot graph size significance
		Map<String, double[]> map = new HashMap<>();
		map.put(ArrayUtil.mean(sizes) + " (mean of random)", sizes);
		map.put(sizeCurrent + " (actual)", new double[]{sizeCurrent});
		KernelDensityPlot.plot("Graph size distribution", map);

		// Convert counts to p-values

		graphSizePval = sizeCnt / (double) iterations;

		this.pvals = new HashMap<>();

		for (String gene : current.keySet())
		{
			double pval;

			if (current.get(gene) == 0)
			{
				pval = 1;
			}
			else
			{
				int c = !cnt.containsKey(gene) ? 0 : cnt.get(gene);
				if (c == 0) c++; // pval = 0 is non-realistic. bring it to the smallest nonzero value
				pval = c / (double) iterations;
			}

			pvals.put(gene, pval);
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
		writer.write("\nGene\tDownstream crowded pval");

		pvals.keySet().stream().sorted((g1, g2) -> pvals.get(g1).compareTo(pvals.get(g2))).forEach(gene ->
			FileUtil.lnwrite(gene + "\t" + pvals.get(gene), writer));

		writer.close();
	}
}
