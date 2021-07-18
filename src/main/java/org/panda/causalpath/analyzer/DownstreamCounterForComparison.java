package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.network.Relation;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * This class is for counting significant changes at the downstream of a protein, as well as how many of those changes
 * suggest that the protein is activated or inactivated.
 *
 * This class does not support correlation-based analysis.
 *
 * @author Ozgun Babur
 */
public class DownstreamCounterForComparison extends DownstreamCounterForCorrelation
{
	Map<Relation, Set<ExperimentData>> rel2Datas;

	public DownstreamCounterForComparison(CausalitySearcher cs, Set<Relation> relations)
	{
		super(cs);

		rel2Datas = new HashMap<>();
		for (Relation rel : relations)
		{
			rel2Datas.put(rel, cs.getExplainableTargetDataWithSiteMatch(rel));
		}
	}

	public Map<String, Integer>[] run()
	{
		Map<String, Set<String>> total = new HashMap<>();
		Map<String, Set<String>> activ = new HashMap<>();
		Map<String, Set<String>> inhib = new HashMap<>();

		Set<Relation> relations = rel2Datas.keySet();

		relations.forEach(r ->
		{
			Set<ExperimentData> td = rel2Datas.get(r);

			if (!td.isEmpty())
			{
				if (!total.containsKey(r.source))
				{
					total.put(r.source, new HashSet<>());
					activ.put(r.source, new HashSet<>());
					inhib.put(r.source, new HashSet<>());
				}
			}

			td.forEach(d ->
			{
				int sign = d.getChangeSign() * r.getSign();

				if (sign != 0)
				{
					total.get(r.source).add(r.target);
				}

				if (sign == 1)
				{
					activ.get(r.source).add(r.target);
				}
				else if (sign == -1)
				{
					inhib.get(r.source).add(r.target);
				}
			});
		});

		return new Map[]{convertToCounts(total), convertToCounts(activ), convertToCounts(inhib)};
	}

	public Map<String, Integer> getGenesPotentialDownstreamMax(Set<Relation> relations)
	{
		Map<String, Set<String>> map = new HashMap<>();
		relations.stream().filter(cs::hasConsiderableDownstreamData).forEach(r ->
		{
			if (!map.containsKey(r.source)) map.put(r.source, new HashSet<>());
			map.get(r.source).add(r.target);
		});

		return convertToCounts(map);
	}

	public Set<String> getGenesWithNoPotential(Set<Relation> relations)
	{
		Set<String> consider = relations.stream().filter(cs::hasConsiderableDownstreamData).map(r -> r.source)
			.collect(Collectors.toSet());

		return relations.stream().map(r -> r.source).filter(gene -> !consider.contains(gene))
			.collect(Collectors.toSet());
	}
}
