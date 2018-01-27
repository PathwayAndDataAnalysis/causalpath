package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;

import java.util.*;
import java.util.stream.Collectors;

/**
 * This class is for shuffling the rows of the data. It first converts data in simpler format, then provides shuffling
 * multiple times. This class performs the shuffling assuming the analysis is correlation-based.
 *
 * @author Ozgun Babur
 */
public class DataLabelShufflerForCorrelation extends DataLabelShuffler
{
	private static final boolean USE_CACHE = false;

	Map<Object, Map<Object, Integer>> correlationMap;

	public DataLabelShufflerForCorrelation(Set<Relation> relations)
	{
		init(relations);
		correlationMap = new HashMap<>();

		TwoDataChangeDetector det = relations.iterator().next().chDet;
		if (USE_CACHE) det = new ChDet(det);
		this.relations.forEach(r -> r.setChDet(det));
	}

	/**
	 * Makes a copy of the given data that is more efficient for randomization experiments.
	 */
	protected ExperimentData convert(ExperimentData orig, Map<String, ExperimentData> dataMap)
	{
		if (dataMap.containsKey(orig.id)) return dataMap.get(orig.id);

		ExperimentData copy = orig.copy();
		dataMap.put(orig.id, copy);

		return copy;
	}

	private Object getInnerData(ExperimentData data)
	{
		if (data instanceof NumericData) return ((NumericData) data).vals;
		else if (data != null) return ((CategoricalData) data).data;
		return null;
	}

	class ChDet implements TwoDataChangeDetector
	{
		TwoDataChangeDetector originalDet;

		public ChDet(TwoDataChangeDetector originalDet)
		{
			this.originalDet = originalDet;
		}

		@Override
		public int getChangeSign(ExperimentData source, ExperimentData target)
		{
			Object os = getInnerData(source);
			Object ot = getInnerData(target);

			if (!correlationMap.containsKey(os)) correlationMap.put(os, new HashMap<>());
			if (!correlationMap.containsKey(ot)) correlationMap.put(ot, new HashMap<>());

			if (!correlationMap.get(os).containsKey(ot))
			{
				int corr = originalDet.getChangeSign(source, target);
				correlationMap.get(os).put(ot, corr);
				correlationMap.get(ot).put(os, corr);
			}

			return correlationMap.get(os).get(ot);
		}
	}
}
