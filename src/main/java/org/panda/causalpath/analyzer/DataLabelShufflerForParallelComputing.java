package org.panda.causalpath.analyzer;

import org.apache.commons.lang.ArrayUtils;
import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;
import org.panda.utility.ArrayUtil;

import java.util.*;
import java.util.stream.Collectors;

/**
 * This class is for shuffling the rows of the data. It first converts data in simpler format, then provides shuffling
 * multiple times. Does not support correlation-based runs.
 *
 * @author Ozgun Babur
 */
public class DataLabelShufflerForParallelComputing implements Cloneable
{
	Map<DataType, Map<ExperimentData, Integer>> changeSignShuffleMap;
	Map<DataType, Integer[]> signArrays;
	Map<ExperimentData, Integer> combinedShuffleMap;

	/**
	 * Empty constructor for the extending class.
	 */
	protected DataLabelShufflerForParallelComputing()
	{
	}

	public DataLabelShufflerForParallelComputing(Set<Relation> relations)
	{
		init(relations);
	}

	protected void init(Set<Relation> relations)
	{
		changeSignShuffleMap = new HashMap<>();
		signArrays = new HashMap<>();
		combinedShuffleMap = new HashMap<>();

		for (DataType type : DataType.values())
		{
			List<ExperimentData> list = relations.stream().map(Relation::getAllData).flatMap(Collection::stream)
				.filter(d -> d.getType().equals(type)).distinct().collect(Collectors.toList());

			if (!list.isEmpty())
			{
				Map<ExperimentData, Integer> signMap = new HashMap<>();
				changeSignShuffleMap.put(type, signMap);
				list.forEach(data -> signMap.put(data, data.getChangeSign()));

				signArrays.put(type, signMap.values().toArray(new Integer[0]));
			}
		}
		updateCombined();
	}

	public void shuffle()
	{
		for (DataType type : changeSignShuffleMap.keySet())
		{
			Map<ExperimentData, Integer> map = changeSignShuffleMap.get(type);
			Set<ExperimentData> datas = map.keySet();

			Integer[] signs = signArrays.get(type);
			Collections.shuffle(Arrays.asList(signs));

			int i = 0;
			for (ExperimentData data : datas)
			{
				map.put(data, signs[i++]);
			}
		}
		updateCombined();
	}

	private void updateCombined()
	{
		changeSignShuffleMap.values().forEach(combinedShuffleMap::putAll);
	}

	public Map<ExperimentData, Integer> getShuffleMap()
	{
		return combinedShuffleMap;
	}

	@Override
	public DataLabelShufflerForParallelComputing clone()
	{
		try
		{
			DataLabelShufflerForParallelComputing clone = (DataLabelShufflerForParallelComputing) super.clone();
			clone.changeSignShuffleMap = new HashMap<>();
			clone.signArrays = new HashMap<>();
			clone.combinedShuffleMap = new HashMap<>(); // leave empty
			for (DataType type : this.changeSignShuffleMap.keySet())
			{
				Map<ExperimentData, Integer> map = this.changeSignShuffleMap.get(type);
				Map<ExperimentData, Integer> mapNew = new HashMap<>(map);
				clone.changeSignShuffleMap.put(type, mapNew);

				Integer[] arrOrig = this.signArrays.get(type);
				Integer[] arrNew = new Integer[arrOrig.length];
				System.arraycopy(arrOrig, 0, arrNew, 0, arrNew.length);
				clone.signArrays.put(type, arrNew);
			}
			return clone;
		}
		catch (CloneNotSupportedException e)
		{
			throw new AssertionError();
		}
	}
}
