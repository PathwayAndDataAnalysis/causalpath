package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * This class is for shuffling the rows of the data. It first converts data in simpler format, then provides shuffling
 * multiple times. This class performs the shuffling assuming the analysis is correlation-based.
 *
 * @author Ozgun Babur
 */
public class DataLabelShufflerForCorrelation
{
	Set<Relation> relations;
	Map<Object, Map<Object, Integer>> correlationMap;
	List<double[]> valuesList;
	List<SingleCategoricalData[]> catsList;
	List<NumericData> numDatas;
	List<CategoricalData> catDatas;

	public DataLabelShufflerForCorrelation(Set<Relation> relations)
	{
		initCorrelationMap(relations);
		this.relations = new HashSet<>();

		Map<String, ExperimentData> dataMap = new HashMap<>();

		for (Relation origRel : relations)
		{
			Relation copy = origRel.copy();
			copy.sourceData.addAll(convert(origRel.sourceData, dataMap));
			copy.targetData.addAll(convert(origRel.targetData, dataMap));
			this.relations.add(copy);
		}

		ChDet det = new ChDet();
		this.relations.forEach(r -> r.setChDet(det));

		numDatas = Stream.concat(this.relations.stream().map(r -> r.sourceData).flatMap(Collection::stream),
			this.relations.stream().map(r -> r.targetData).flatMap(Collection::stream))
			.filter(d -> d instanceof NumericData).map(d -> (NumericData) d).distinct()
			.collect(Collectors.toList());

		catDatas = Stream.concat(this.relations.stream().map(r -> r.sourceData).flatMap(Collection::stream),
			this.relations.stream().map(r -> r.targetData).flatMap(Collection::stream))
			.filter(d -> d instanceof CategoricalData).map(d -> (CategoricalData) d).distinct()
			.collect(Collectors.toList());

		valuesList = numDatas.stream().map(d -> d.vals).collect(Collectors.toList());
		catsList = catDatas.stream().map(d -> d.data).collect(Collectors.toList());
	}

	public void shuffle()
	{
		Collections.shuffle(valuesList);
		Collections.shuffle(catsList);
		for (int i = 0; i < numDatas.size(); i++)
		{
			numDatas.get(i).vals = valuesList.get(i);
		}
		for (int i = 0; i < catDatas.size(); i++)
		{
			catDatas.get(i).data = catsList.get(i);
		}
	}

	public Set<Relation> getRelations()
	{
		return relations;
	}

	private Set<ExperimentData> convert(Set<ExperimentData> origDatas, Map<String, ExperimentData> dataMap)
	{
		return origDatas.stream().map(d -> convert(d, dataMap)).filter(d -> d != null).collect(Collectors.toSet());
	}

	/**
	 * Makes a copy of the given data that is more efficient for randomization experiments.
	 */
	private ExperimentData convert(ExperimentData orig, Map<String, ExperimentData> dataMap)
	{
		if (dataMap.containsKey(orig.id)) return dataMap.get(orig.id);

		ExperimentData copy = orig.copy();
		dataMap.put(orig.id, copy);

		return copy;
	}

	private void initCorrelationMap(Set<Relation> relations)
	{
		correlationMap = new HashMap<>();

		// Get experiment data in the relations
		Set<ExperimentData> datas = Stream.concat(relations.stream().map(r -> r.sourceData),
			relations.stream().map(r -> r.targetData)).flatMap(Collection::stream).collect(Collectors.toSet());

		TwoDataChangeDetector det = relations.iterator().next().chDet;

		for (ExperimentData d1 : datas)
		{
			Object o1 = getInnerData(d1);

			for (ExperimentData d2 : datas)
			{
				if (d1.id.compareTo(d2.id) >= 0) continue;

				Object o2 = getInnerData(d2);

				int corr = det.getChangeSign(d1, d2);

				if (!correlationMap.containsKey(o1)) correlationMap.put(o1, new HashMap<>());
				if (!correlationMap.containsKey(o2)) correlationMap.put(o2, new HashMap<>());
				correlationMap.get(o1).put(o2, corr);
				correlationMap.get(o2).put(o1, corr);
			}
		}
	}

	private Object getInnerData(ExperimentData data)
	{
		if (data instanceof NumericData) return ((NumericData) data).vals;
		else if (data != null) return ((CategoricalData) data).data;
		return null;
	}

	class ChDet implements TwoDataChangeDetector
	{
		@Override
		public int getChangeSign(ExperimentData source, ExperimentData target)
		{
			return correlationMap.get(getInnerData(source)).get(getInnerData(target));
		}
	}
}
