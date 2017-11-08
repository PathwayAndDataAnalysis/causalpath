package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;

import java.util.*;
import java.util.stream.Collectors;

/**
 * This class is for shuffling the rows of the data. It first converts data in simpler format, then provides shuffling
 * multiple times. Does not support correlation-based runs.
 *
 * @author Ozgun Babur
 */
public class DataLabelShuffler
{
	Set<Relation> relations;

	Map<DataType, List<NumericData>> numDataLists;
	Map<DataType, List<CategoricalData>> catDataLists;
	Map<DataType, List<double[]>> numDataVals;
	Map<DataType, List<SingleCategoricalData[]>> catDataVals;

	/**
	 * Empty constructor for the extending class.
	 */
	protected DataLabelShuffler()
	{
	}

	public DataLabelShuffler(Set<Relation> relations)
	{
		init(relations);

		CausalityHelper ch = new CausalityHelper();
		this.relations.forEach(r -> r.setChDet(ch));

		ChDet det = new ChDet();
		this.relations.stream().map(Relation::getAllData).flatMap(Collection::stream).distinct()
			.forEach(d -> d.setChDet(det));
	}

	protected void init(Set<Relation> relations)
	{
		this.relations = new HashSet<>();

		Map<String, ExperimentData> dataMap = new HashMap<>();

		for (Relation origRel : relations)
		{
			Relation copy = origRel.copy();
			copy.sourceData = convert(origRel.sourceData, dataMap);
			copy.targetData = convert(origRel.targetData, dataMap);
			this.relations.add(copy);
		}

		numDataLists = new HashMap<>();
		catDataLists = new HashMap<>();
		numDataVals = new HashMap<>();
		catDataVals = new HashMap<>();
		for (DataType type : DataType.values())
		{
			List<ExperimentData> list = this.relations.stream().map(Relation::getAllData).flatMap(Collection::stream)
				.filter(d -> d.getType().equals(type)).distinct().collect(Collectors.toList());

			if (!list.isEmpty())
			{
				if (type.isNumerical())
				{
					numDataLists.put(type, list.stream().map(d -> (NumericData) d).collect(Collectors.toList()));
					numDataVals.put(type, numDataLists.get(type).stream().map(d -> d.vals).collect(Collectors.toList()));
				}
				else
				{
					catDataLists.put(type, list.stream().map(d -> (CategoricalData) d).collect(Collectors.toList()));
					catDataVals.put(type, catDataLists.get(type).stream().map(d -> d.data).collect(Collectors.toList()));
				}
			}
		}
	}

	public void shuffle()
	{
		for (DataType type : numDataLists.keySet())
		{
			List<double[]> vals = numDataVals.get(type);
			Collections.shuffle(vals);

			List<NumericData> list = numDataLists.get(type);
			for (int i = 0; i < list.size(); i++)
			{
				list.get(i).vals = vals.get(i);
			}
		}
		for (DataType type : catDataLists.keySet())
		{
			List<SingleCategoricalData[]> datas = catDataVals.get(type);
			Collections.shuffle(datas);

			List<CategoricalData> list = catDataLists.get(type);
			for (int i = 0; i < list.size(); i++)
			{
				list.get(i).data = datas.get(i);
			}
		}
	}

	public Set<Relation> getRelations()
	{
		return relations;
	}

	private GeneWithData convert(GeneWithData orig, Map<String, ExperimentData> dataMap)
	{
		GeneWithData copy = new GeneWithData(orig.getId());
		for (ExperimentData data : orig.getData())
		{
			copy.add(convert(data, dataMap));
		}
		return copy;
	}

	/**
	 * Makes a copy of the given data that is more efficient for randomization experiments.
	 */
	protected ExperimentData convert(ExperimentData orig, Map<String, ExperimentData> dataMap)
	{
		if (dataMap.containsKey(orig.id)) return dataMap.get(orig.id);

		ExperimentData copy = orig.copy();

		if (copy instanceof NumericData)
		{
			((NumericData) copy).vals = new double[]{orig.getChangeSign()};
		}
		else if (copy instanceof MutationData)
		{
			((MutationData) copy).data = new SingleCategoricalData[]{new Mutation(orig.getChangeSign(), null)};
		}
		else if (copy instanceof CNAData)
		{
			((CNAData) copy).data = new SingleCategoricalData[]{new CNA(orig.getChangeSign())};
		}
		else if (copy instanceof ActivityData)
		{
			((ActivityData) copy).data = new SingleCategoricalData[]{new Activity(orig.getChangeSign())};
		}

		dataMap.put(orig.id, copy);
		return copy;
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

		@Override
		public OneDataChangeDetector makeACopy()
		{
			throw new UnsupportedOperationException("This is not a cloneable detector! Sorry.");
		}
	}
}
