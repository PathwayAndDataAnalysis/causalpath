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
public abstract class NetworkSignificanceCalculator
{
	/**
	 * The network.
	 */
	protected Set<Relation> relations;

	/**
	 * Compatibility for source-relation-target triplet, without considering data compatibility.
	 */
	protected RelationTargetCompatibilityChecker rtcc;

	/**
	 * P-value for the current graph size.
	 */
	protected double graphSizePval;

	/**
	 * Doing for either causal or conflicting graph.
	 */
	protected boolean causal;

	/**
	 * The p-value threshold fora significant change.
	 */
	protected double significanceThreshold;

	/**
	 * If set, this filter selects a subset of the detected relations as the result.
	 */
	protected GraphFilter graphFilter;

	/**
	 * Constructor with the network.
	 */
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
		this.graphFilter = graphFilter;
	}

	public void setSignificanceThreshold(double significanceThreshold)
	{
		this.significanceThreshold = significanceThreshold;
	}

	/**
	 * Performs a randomization experiment.
	 */
	public abstract void run(int iterations);

	public double getOverallGraphSizePval()
	{
		return graphSizePval;
	}

	public abstract Map<String, Double> getDownstreamActivityPvals();

	public boolean isDownstreamSignificant(String gene)
	{
		return getDownstreamActivityPvals().containsKey(gene) && getDownstreamActivityPvals().get(gene) <= significanceThreshold;
	}

	public abstract void writeResults(String filename) throws IOException;
}
