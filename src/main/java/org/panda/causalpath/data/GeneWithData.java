package org.panda.causalpath.data;

import org.panda.causalpath.analyzer.OneDataChangeDetector;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A gene with multiple omics data associated.
 *
 * @author Ozgun Babur
 */
public class GeneWithData
{
	/**
	 * ID of the gene, typically the gene symbol
	 */
	String id;

	/**
	 * The associated data, grouped by data type.
	 */
	Map<DataType, Set<ExperimentData>> dataMap;

	/**
	 * Constructor with ID.
	 */
	public GeneWithData(String id)
	{
		this.id = id;
		dataMap = new HashMap<>();
	}

	/**
	 * @return the ID
	 */
	public String getId()
	{
		return id;
	}

	/**
	 * Adds another omics data.
	 * @param data data to add
	 */
	public void add(ExperimentData data)
	{
		DataType type = data.getType();

		if (!dataMap.containsKey(type)) dataMap.put(type, new HashSet<>());
		dataMap.get(type).add(data);
	}

	/**
	 * Adds multiple omic data.
	 * @param col the collection of data
	 */
	public void addAll(Collection<ExperimentData> col)
	{
		if (col == null) return;
		col.forEach(this::add);
	}

	/**
	 * Gets the data that are changed. All data has to be associated with a change detector before calling this method.
	 * @return changed data
	 */
	public Map<ExperimentData, Integer> getChangedData()
	{
		Map<ExperimentData, Integer> map = new HashMap<>();

		for (DataType type : dataMap.keySet())
		{
			map.putAll(getChangedData(type));
		}

		return map;
	}

	/**
	 * Gets changed data of specific type.
	 * @param type the data type
	 * @return changed data of specific type
	 */
	public Map<ExperimentData, Integer> getChangedData(DataType type)
	{
		Map<ExperimentData, Integer> map = new HashMap<>();
		for (ExperimentData data : getData(type))
		{
			int sign = data.getChangeSign();
			if (sign != 0) map.put(data, sign);
		}

		return map;
	}

	/**
	 * Sets the given change detector to all data associated.
	 * @param chDet to set
	 */
	public void setChangeDet(OneDataChangeDetector chDet)
	{
		getDataStream().forEach(d -> d.setChDet(chDet));
	}

	/**
	 * Gets the associated data as stream.
	 * @return data as stream
	 */
	public Stream<ExperimentData> getDataStream()
	{
		return dataMap.values().stream().flatMap(Collection::stream);
	}

	/**
	 * Gets all associated data.
	 * @return associated data
	 */
	public Set<ExperimentData> getData()
	{
		return getDataStream().collect(Collectors.toSet());
	}

	/**
	 * Gets associated data of specific type.
	 * @param type desired data type
	 * @return associated data of the given type
	 */
	public Set<ExperimentData> getData(DataType type)
	{
		if (dataMap.containsKey(type)) return dataMap.get(type);
		return Collections.emptySet();
	}

	/**
	 * Checks if the gene has no data.
	 * @return true if gene has no data
	 */
	public boolean isEmpty()
	{
		return dataMap.isEmpty();
	}

	/**
	 * Tells if this gene is associated with a changed prooteomics or phosphoproteomics data.
	 * @return true if there is any
	 */
	public boolean hasChangedProteomicData()
	{
		return !getChangedData(DataType.PROTEIN).keySet().isEmpty() ||
			!getChangedData(DataType.PHOSPHOPROTEIN).keySet().isEmpty();
	}

	/**
	 * Gets the set of chnaged proteomic data
	 * @return proteomic data with a significant change
	 */
	public Set<ExperimentData> getChangedProteomicData()
	{
		Set<ExperimentData> set = new HashSet<>(getChangedData(DataType.PROTEIN).keySet());
		set.addAll(getChangedData(DataType.PHOSPHOPROTEIN).keySet());
		return set;
	}
}
