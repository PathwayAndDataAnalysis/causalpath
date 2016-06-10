package org.panda.causalpath.network;

import org.panda.causalpath.analyzer.TwoDataChangeDetector;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.PhosphoSite;
import org.panda.utility.CollectionUtil;

import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * That is a potential causality edge.
 *
 * Created by babur on 3/24/16.
 */
public class Relation
{
	public String source;
	public String target;
	public RelationType type;

	public Set<ExperimentData> sourceData;
	public Set<ExperimentData> targetData;

	private String mediators;

	public TwoDataChangeDetector chDet;

	// target sites, used when this is phospho relation.
	public Set<PhosphoSite> sites;

	public Relation(String source, String target, RelationType type, String mediators)
	{
		this.source = source;
		this.target = target;
		this.type = type;
		this.mediators = mediators;
		sourceData = new HashSet<>();
		targetData = new HashSet<>();
	}

	public int getSign()
	{
		return type.sign;
	}

	public void pickRelatedExperimentData(Map<String, Set<ExperimentData>> map)
	{
		if (map.containsKey(source)) sourceData.addAll(map.get(source));
		if (map.containsKey(target)) targetData.addAll(map.get(target));
	}

	public String toString()
	{
		return source + "\t" + type.name + "\t" + target + "\t" + mediators +
			(sites == null ? "" : "\t" + CollectionUtil.merge(sites, ";"));
	}

	public Set<String> getTargetWithSites()
	{
		if (sites == null) return Collections.emptySet();
		Set<String> set = new HashSet<>();
		for (PhosphoSite site : sites)
		{
			set.add(target + "_" + site);
		}
		return set;
	}

	@Override
	public int hashCode()
	{
		return source.hashCode() + target.hashCode() + type.hashCode();
	}

	@Override
	public boolean equals(Object obj)
	{
		return obj instanceof Relation && ((Relation) obj).source.equals(source) &&
			((Relation) obj).target.equals(target) && ((Relation) obj).type.equals(type);
	}
}
