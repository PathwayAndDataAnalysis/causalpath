package org.panda.causalpath.analyzer;

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
public class DownstreamCounterForCorrelation
{
	CausalitySearcher cs;

	public DownstreamCounterForCorrelation(CausalitySearcher cs)
	{
		this.cs = cs;
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

	protected Map<String, Integer> convertToCounts(Map<String, Set<String>> map)
	{
		return map.keySet().stream().collect(Collectors.toMap(s -> s, s -> map.get(s).size()));
	}

	public Map<String, Integer> getGenesPotentialDownstreamMax(Set<Relation> relations)
	{
		Map<String, Set<String>> map = new HashMap<>();
		relations.stream().filter(cs::hasConsiderableData).forEach(r ->
		{
			if (!map.containsKey(r.source)) map.put(r.source, new HashSet<>());
			map.get(r.source).add(r.target);
		});

		return convertToCounts(map);
	}

	public Set<String> getGenesWithNoPotential(Set<Relation> relations)
	{
		Set<String> consider = relations.stream().filter(cs::hasConsiderableData).map(r -> r.source)
			.collect(Collectors.toSet());

		return relations.stream().map(r -> r.source).filter(gene -> !consider.contains(gene))
			.collect(Collectors.toSet());
	}
}