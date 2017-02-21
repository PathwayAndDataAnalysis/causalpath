package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.network.Relation;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * This class is for counting significant correlations at the downstream of a protein.
 *
 * @author Ozgun Babur
 */
public class DownstreamCounterForCorrelation extends DownstreamCounter
{
	/**
	 * Are we interested in causal or a conflicting network.
	 */
	boolean causal;

	public DownstreamCounterForCorrelation(RelationTargetCompatibilityChecker rtcc, boolean causal)
	{
		super(rtcc);
		this.causal = causal;
	}

	public Map<String, Integer>[] run(Set<Relation> relations)
	{
		Map<String, Set<String>> total = new HashMap<>();

		relations.stream().filter(this::hasCompatibleCorrelation).forEach(r ->
		{
			if (!total.containsKey(r.source)) total.put(r.source, new HashSet<>());
			total.get(r.source).add(r.target);
		});

		return new Map[]{convertToCounts(total)};
	}

	/**
	 * Checks if any of the source target data pair has a compatible
	 */
	private boolean hasCompatibleCorrelation(Relation r)
	{
		for (ExperimentData sourceData : r.sourceData)
		{
			for (ExperimentData targetData : r.targetData)
			{
				if (rtcc.isCompatible(sourceData, r, targetData) &&
					r.chDet.getChangeSign(sourceData, targetData) * r.getSign() == (causal ? 1 : -1))
				{
					return true;
				}
			}
		}
		return false;
	}
}