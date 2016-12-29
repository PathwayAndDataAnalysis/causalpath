package org.panda.causalpath.resource;

import org.panda.causalpath.analyzer.CausalityHelper;
import org.panda.causalpath.analyzer.OneDataChangeDetector;
import org.panda.causalpath.data.ActivityData;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.PhosphoProteinData;
import org.panda.causalpath.data.ProteinData;
import org.panda.causalpath.network.Relation;
import org.panda.resource.tcga.ProteomicsFileRow;

import java.util.*;

/**
 * Coverts the RPPAData in resource project to appropriate objects for this project and serves them.
 */
public class ProteomicsLoader
{
	/**
	 * Map from genes to related data.
	 */
	Map<String, Set<ExperimentData>> dataMap;

	/**
	 * Set of all data. This collection supposed to hold everything in the dataMap's values.
	 */
	Set<ExperimentData> datas;

	/**
	 * Initializes self using a set of RPPAData.
	 */
	public ProteomicsLoader(Collection<ProteomicsFileRow> rows)
	{
		dataMap = new HashMap<>();
		datas = new HashSet<>();
		rows.stream().distinct().forEach(r ->
		{
			ExperimentData ed = r.isActivity() ? new ActivityData(r) :
				r.isPhospho() ? new PhosphoProteinData(r) : new ProteinData(r);

			for (String sym : ed.getGeneSymbols())
			{
				if (!dataMap.containsKey(sym)) dataMap.put(sym, new HashSet<>());
				dataMap.get(sym).add(ed);
				datas.add(ed);
			}
		});
	}

	/**
	 * Adds the related data to the given relations.
	 */
	public void decorateRelations(Set<Relation> relations)
	{
		CausalityHelper ch = new CausalityHelper();
		for (Relation relation : relations)
		{
			if (dataMap.containsKey(relation.source)) relation.sourceData.addAll(dataMap.get(relation.source));
			if (dataMap.containsKey(relation.target)) relation.targetData.addAll(dataMap.get(relation.target));
			relation.chDet = ch;
		}
	}

	/**
	 * Puts the given change detector to the data that is filtered by the given selector.
	 */
	public void associateChangeDetector(OneDataChangeDetector chDet, DataSelector selector)
	{
		datas.stream().filter(selector::select).forEach(d -> d.setChDet(chDet));
	}

	/**
	 * Function to filter experiment data.
	 */
	public interface DataSelector
	{
		boolean select(ExperimentData data);
	}
}
