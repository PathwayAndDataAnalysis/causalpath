package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.ProteinData;
import org.panda.causalpath.network.Relation;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

/**
 * This class predicts the accuracy of network propagation by assuming a known proteomics data is unknown and checking
 * the majority vote at the upstream of that node.
 */
public class PropagationAccuracyPredictor
{
	public double run(Set<Relation> causal, Set<Relation> conflicting, CausalitySearcher cs)
	{
		cs.setCausal(true);
		Map<ExperimentData, Set<ExperimentData>> causalUpstr = getTargetToUpstream(causal, cs);
		cs.setCausal(false);
		Map<ExperimentData, Set<ExperimentData>> conflictUpstr = getTargetToUpstream(conflicting, cs);

		// predict accuracy

		int[] successFail = new int[]{0, 0};

		Stream.concat(causalUpstr.keySet().stream(), conflictUpstr.keySet().stream()).distinct().forEach(target ->
		{
			int approve = causalUpstr.containsKey(target) ? causalUpstr.get(target).size() : 0;
			int disapp = conflictUpstr.containsKey(target) ? conflictUpstr.get(target).size() : 0;

			if (approve > 0 && disapp == 0) successFail[0]++;
			else if (approve == 0 && disapp > 0) successFail[1]++;
		});

		return successFail[0] / (double) (successFail[0] + successFail[1]);
	}

	private Map<ExperimentData, Set<ExperimentData>> getTargetToUpstream(Set<Relation> results, CausalitySearcher cs)
	{
		Map<ExperimentData, Set<ExperimentData>> map = new HashMap<>();
		for (Relation rel : results)
		{
			for (ExperimentData target : cs.getExplainableTargetData(rel))
			{
				if (!(target instanceof ProteinData)) continue;

				for (ExperimentData source : cs.getSatisfyingSourceData(rel, target))
				{
					if (!map.containsKey(target)) map.put(target, new HashSet<>());
					map.get(target).add(source);
				}
			}
		}
		return map;
	}
}
