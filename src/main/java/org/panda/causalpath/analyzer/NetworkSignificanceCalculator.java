package org.panda.causalpath.analyzer;

import org.panda.causalpath.network.GraphFilter;
import org.panda.causalpath.network.Relation;
import org.panda.utility.statistics.FDR;

import java.io.IOException;
import java.util.Map;
import java.util.Set;

/**
 * Calculates the significance of several things in the result network. These are the size of the network overall, the
 * significance of downstream change frequency at the downstream of proteins, and the significance of downstream
 * evidence for activation/inactivation of the protein.
 *
 * @author Ozgun Babur
 */
public abstract class NetworkSignificanceCalculator
{
	/**
	 * The network.
	 */
	protected Set<Relation> relations;

	/**
	 * Compatibility detector for source-relation-target triplet.
	 */
	protected CausalitySearcher cs;

	/**
	 * P-value for the current graph size.
	 */
	protected double graphSizePval;

	/**
	 * The p-value threshold fora significant change.
	 */
	protected double significanceThreshold;

	protected int minimumPotentialTargetsToConsider;

	/**
	 * Constructor with the network.
	 */
	public NetworkSignificanceCalculator(Set<Relation> relations, CausalitySearcher cs)
	{
		this.relations = relations;
		this.cs = cs;
		cs.initRelationDataMappingMemory();
		this.minimumPotentialTargetsToConsider = 0;
	}

	public void setMinimumPotentialTargetsToConsider(int minimumPotentialTargetsToConsider)
	{
		this.minimumPotentialTargetsToConsider = minimumPotentialTargetsToConsider;
	}

	public void setPvalThreshold(double significanceThreshold)
	{
		this.significanceThreshold = significanceThreshold;
	}

	/**
	 * Sets the FDR threshold for the network significance.
	 *
	 * @param fdrThr the FDR threshold
	 */
	public abstract void setFDRThreshold(double fdrThr);

	/**
	 * Performs a randomization experiment.
	 */
	public abstract void run(int iterations);

	/**
	 * Gets the p-val for the network size.
	 *
	 * @return p-val
	 */
	public double getOverallGraphSizePval()
	{
		return graphSizePval;
	}

	/**
	 * Gets the map of downstream activity size p-vals for each gene.
	 *
	 * @return p-vals map
	 */
	public abstract Map<String, Double> getDownstreamActivityPvals();

	/**
	 * Checks if the amount of downstream activities are significantly large.
	 *
	 * @param gene gene of interest
	 * @return true if significantly large
	 */
	public boolean isDownstreamSignificant(String gene)
	{
		return getDownstreamActivityPvals().containsKey(gene) && getDownstreamActivityPvals().get(gene) <= significanceThreshold;
	}

	/**
	 * Write gene significances to a file.
	 *
	 * @param filename file name
	 * @throws IOException
	 */
	public abstract void writeResults(String filename) throws IOException;

	/**
	 * Reads significances from the given file. This method is mostly for debugging purposes, when we don't want to
	 * spend too much time with calculations but to use results from a previous calculation.
	 * @param filename file name
	 * @throws IOException
	 */
	public abstract void loadFromFile(String filename) throws IOException;
}