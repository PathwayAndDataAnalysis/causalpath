package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationAndSelectedData;
import org.panda.utility.CollectionUtil;

import java.util.*;

/**
 * This class known which experiments data types can be used to detect causality for which
 * relation types.
 *
 * Created by babur on 4/5/16.
 */
public class CausalitySearcher
{
	// phospho relations only
	static Set<Class<? extends ExperimentData>> pSet = Collections.singleton(PhosphoProteinData.class);

	// expression relations only
//	static Set<Class<? extends ExperimentData>> eSet = new HashSet<>(Arrays.asList(ProteinData.class, ExpressionData.class));
	static Set<Class<? extends ExperimentData>> eSet = Collections.singleton(ProteinData.class);

	private boolean forceSiteMatching = true;
	private boolean addInUnknownSigns = false;
	private int siteProximityThreshold = 0;
	private Set<String> genesWithTotalProteinData;
	/**
	 * If that is false, then we are interested in conflicting relations.
	 */
	private boolean causal = true;

	/**
	 * Finds compatible or conflicting relations. The relations have to be associated with experiment data. Both the
	 * experiment data and the relations have to be associated with related change detectors.
	 */
	public Set<RelationAndSelectedData> run(Set<Relation> relations)
	{
		Set<RelationAndSelectedData> selected = new HashSet<>();
		int compatible = causal ? 1 : -1;

		for (Relation relation : relations)
		{
			for (ExperimentData source : relation.sourceData)
			{
				if (skip(source, relation.source)) continue;
//				if (source instanceof ExpressionData) continue;
//				if (!(source instanceof MutationData)) continue;

				for (ExperimentData target : relation.targetData)
				{
					if (skip(target, relation.target)) continue;
					if (relation.type.affectsPhosphoSite && !pSet.contains(target.getClass())) continue;
					if (relation.type.affectsTotalProt && ! eSet.contains(target.getClass())) continue;

					// Check if site-matching constraint holds
					if (forceSiteMatching && relation.type.affectsPhosphoSite && target instanceof PhosphoProteinData &&
						CollectionUtil.intersectionEmpty(relation.getTargetWithSites(siteProximityThreshold),
							((PhosphoProteinData) target).getGenesWithSites())) continue;

					int chgSign = relation.chDet.getChangeSign(source, target);
					int sign = chgSign * source.getEffect() * relation.getSign() * compatible;
					boolean causative = addInUnknownSigns ? sign != -1 : sign == 1;
					if (causative)
					{
						selected.add(new RelationAndSelectedData(relation, source, target));
					}
				}
			}
		}
		return selected;
	}

	private boolean skip(ExperimentData data, String gene)
	{
		return genesWithTotalProteinData != null && genesWithTotalProteinData.contains(gene) &&
		(data instanceof ExpressionData || data instanceof CNAData);
	}

	public void setGenesWithTotalProteinData(Set<String> genesWithTotalProteinData)
	{
		this.genesWithTotalProteinData = genesWithTotalProteinData;
	}

	public void setSiteProximityThreshold(int thr)
	{
		this.siteProximityThreshold = thr;
	}

	public void setCausal(boolean causal)
	{
		this.causal = causal;
	}

	public void setForceSiteMatching(boolean forceSiteMatching)
	{
		this.forceSiteMatching = forceSiteMatching;
	}

	public void setAddInUnknownSigns(boolean add)
	{
		this.addInUnknownSigns = add;
	}
}
