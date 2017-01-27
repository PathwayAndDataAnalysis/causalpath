package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * This class is for shuffling the rows of the data. It first converts data in simpler format, then provides shuffling
 * multiple times.
 *
 * @author Ozgun Babur
 */
public class DataLabelShuffler
{
	Set<Relation> relations;
	List<double[]> valuesList;
	List<SingleCategoricalData[]> catsList;
	List<NumericData> numDatas;
	List<CategoricalData> catDatas;


	public DataLabelShuffler(Set<Relation> relations)
	{
		this.relations = new HashSet<>();

		Map<String, ExperimentData> dataMap = new HashMap<>();

		for (Relation origRel : relations)
		{
			Relation copy = origRel.copy();
			copy.sourceData.addAll(convert(origRel.sourceData, dataMap));
			copy.targetData.addAll(convert(origRel.targetData, dataMap));
			this.relations.add(copy);
		}

		CausalityHelper ch = new CausalityHelper();
		this.relations.forEach(r -> r.setChDet(ch));

		numDatas = Stream.concat(this.relations.stream().map(r -> r.sourceData).flatMap(Collection::stream),
			this.relations.stream().map(r -> r.targetData).flatMap(Collection::stream))
			.filter(d -> d instanceof NumericData).map(d -> (NumericData) d).distinct()
			.collect(Collectors.toList());

		catDatas = Stream.concat(this.relations.stream().map(r -> r.sourceData).flatMap(Collection::stream),
			this.relations.stream().map(r -> r.targetData).flatMap(Collection::stream))
			.filter(d -> d instanceof CategoricalData).map(d -> (CategoricalData) d).distinct()
			.collect(Collectors.toList());

		ChDet det = new ChDet();
		numDatas.forEach(d -> d.setChDet(det));
		catDatas.forEach(d -> d.setChDet(det));

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

		if (orig instanceof PhosphoProteinData)
		{
			PhosphoProteinData copy = new PhosphoProteinData(orig.id, orig.getGeneSymbols());
			copy.setSiteMap(((PhosphoProteinData) orig).getSiteMap());
			copy.vals = new double[]{orig.getChangeSign()};
			dataMap.put(copy.id, copy);
			return copy;
		}
		else if (orig instanceof ProteinData)
		{
			ProteinData copy = new ProteinData(orig.id, orig.getGeneSymbols());
			copy.vals = new double[]{orig.getChangeSign()};
			dataMap.put(copy.id, copy);
			return copy;
		}
		else if (orig instanceof MutationData)
		{
			MutationData copy = new MutationData(orig.id, orig.getGeneSymbols().iterator().next(), orig.getEffect());
			copy.data = new SingleCategoricalData[]{new Mutation(orig.getChangeSign(), null)};
			dataMap.put(copy.id, copy);
			return copy;
		}
		else if (orig instanceof CNAData)
		{
			CNAData copy = new CNAData(orig.id, orig.getGeneSymbols().iterator().next());
			copy.data = new SingleCategoricalData[]{new CNA(orig.getChangeSign())};
			dataMap.put(copy.id, copy);
			return copy;
		}
		else if (orig instanceof ExpressionData)
		{
			ExpressionData copy = new ExpressionData(orig.id, orig.getGeneSymbols().iterator().next());
			copy.vals = new double[]{orig.getChangeSign()};
			dataMap.put(copy.id, copy);
			return copy;
		}

		return null;
	}

	class ChDet implements OneDataChangeDetector
	{
		@Override
		public int getChangeSign(ExperimentData data)
		{
			if (data instanceof NumericData)
			{
				return (int) ((NumericData) data).vals[0];
			}
			else
			{
				return ((CategoricalData) data).data[0].getCategory();
			}
		}

		@Override
		public double getChangeValue(ExperimentData data)
		{
			return getChangeSign(data);
		}
	}
}
