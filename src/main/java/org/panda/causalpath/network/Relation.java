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
 * That is a potential causality relation.
 */
public class Relation
{
	/**
	 * Source gene.
	 */
	public String source;

	/**
	 * Target gene.
	 */
	public String target;

	/**
	 * The type of the relation.
	 */
	public RelationType type;

	/**
	 * Set of experiment data associated with source gene.
	 */
	public Set<ExperimentData> sourceData;

	/**
	 * Set of experiment data associated with target gene.
	 */
	public Set<ExperimentData> targetData;

	/**
	 * Pathway Commons IDs of the mediator objects for the relation.
	 */
	private String mediators;

	/**
	 * A change detector that can evaluate pairs of experiment data to see if this relation can explain the causality
	 * between them.
	 */
	public TwoDataChangeDetector chDet;

	/**
	 * Sites for the target. Needed when the relation is a phosphorylation or dephosphorylation.
	 */
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

	/**
	 * Sign of the relation.
	 */
	public int getSign()
	{
		return type.sign;
	}

	/**
	 * Fills the data resources using the given bag.
	 */
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

	public Set<String> getTargetWithSites(int proximityThr)
	{
		if (sites == null) return Collections.emptySet();
		Set<String> set = new HashSet<>();
		for (PhosphoSite site : sites)
		{
			for (int i = 0; i <= proximityThr; i++)
			{
				set.add(target + "_" + (site.getSite() + i));
				set.add(target + "_" + (site.getSite() - i));
			}
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
