package org.panda.causalpath.resource;

import org.panda.causalpath.analyzer.CausalityHelper;
import org.panda.causalpath.analyzer.OneDataChangeDetector;
import org.panda.causalpath.data.ActivityData;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.PhosphoProteinData;
import org.panda.causalpath.data.ProteinData;
import org.panda.causalpath.network.Relation;
import org.panda.resource.tcga.RPPAData;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class RPPALoader
{
	Map<String, Set<ExperimentData>> dataMap;
	Set<ExperimentData> datas;

	public RPPALoader(Collection<RPPAData> rppas)
	{
		dataMap = new HashMap<>();
		datas = new HashSet<>();
		rppas.stream().distinct().forEach(r ->
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

	public void associateChangeDetector(OneDataChangeDetector chDet, DataSelector selector)
	{
		datas.stream().filter(selector::select).forEach(d -> d.setChDet(chDet));
	}

	public interface DataSelector
	{
		boolean select(ExperimentData data);
	}
}
