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
public class DownstreamCounter
{
	CausalitySearcher cs;

	public DownstreamCounter(CausalitySearcher cs)
	{
		this.cs = cs;
	}

	public Map<String, Integer>[] run(Set<Relation> relations)
	{
		Map<String, Set<String>> total = new HashMap<>();
		Map<String, Set<String>> activ = new HashMap<>();
		Map<String, Set<String>> inhib = new HashMap<>();

		relations.forEach(r ->
		{
			Set<ExperimentData> td = cs.getExplainableTargetData(r);

			if (!td.isEmpty())
			{
				if (!total.containsKey(r.source)) total.put(r.source, new HashSet<>());
				if (!activ.containsKey(r.source)) activ.put(r.source, new HashSet<>());
				if (!inhib.containsKey(r.source)) inhib.put(r.source, new HashSet<>());
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

	protected Map<String, Integer> convertToCounts(Map<String, Set<String>> map)
	{
		return map.keySet().stream().collect(Collectors.toMap(s -> s, s -> map.get(s).size()));
	}
}
