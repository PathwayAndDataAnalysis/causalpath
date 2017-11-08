package org.panda.causalpath.analyzer;

import org.panda.causalpath.network.Relation;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * This class is for counting significant correlations at the downstream of a protein.
 *
 * @author Ozgun Babur
 */
public class DownstreamCounterForCorrelation extends DownstreamCounter
{
	public DownstreamCounterForCorrelation(CausalitySearcher cs)
	{
		super(cs);
	}

	public Map<String, Integer>[] run(Set<Relation> relations)
	{
		Map<String, Set<String>> total = new HashMap<>();

		relations.stream().filter(cs::satisfiesCriteria).forEach(r ->
		{
			if (!total.containsKey(r.source)) total.put(r.source, new HashSet<>());
			total.get(r.source).add(r.target);
		});

		return new Map[]{convertToCounts(total)};
	}
}