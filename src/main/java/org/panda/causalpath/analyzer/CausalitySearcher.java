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
	private Set<String> genesWithTotalProteinData;

	/**
	 * The relations have to be associated with experiment data. Both the experiment data and the relations have to be
	 * associated with related change detectors.
	 */
	public Set<RelationAndSelectedData> selectCausalRelations(Set<Relation> relations)
	{
		return selectCausalRelations(relations, 1);
	}

	/**
	 * The relations have to be associated with experiment data. Both the experiment data and the relations have to be
	 * associated with related change detectors.
	 */
	public Set<RelationAndSelectedData> selectConflictingRelations(Set<Relation> relations)
	{
		return selectCausalRelations(relations, -1);
	}

	/**
	 * Finds compatible or conflicting relations.
	 * @param compatible 1 if compatible, -1 is conflicting
	 */
	private Set<RelationAndSelectedData> selectCausalRelations(Set<Relation> relations, int compatible)
	{
		Set<RelationAndSelectedData> selected = new HashSet<>();

		for (Relation relation : relations)
		{
			for (ExperimentData source : relation.sourceData)
			{
				if (skip(source, relation.source)) continue;
//				if (source instanceof ExpressionData) continue;
				if (!(source instanceof MutationData)) continue;

				for (ExperimentData target : relation.targetData)
				{
					if (skip(target, relation.target)) continue;
					if (relation.type.affectsPhosphoSite && !pSet.contains(target.getClass())) continue;
					if (relation.type.affectsTotalProt && ! eSet.contains(target.getClass())) continue;

					// Check if site-matching constraint holds
					if (forceSiteMatching && relation.type.affectsPhosphoSite && target instanceof PhosphoProteinData &&
						CollectionUtil.intersectionEmpty(relation.getTargetWithSites(),
							((PhosphoProteinData) target).getGenesWithSites())) continue;

					int chgSign = relation.chDet.getChangeSign(source, target);
					boolean causative = chgSign * source.getEffect() * relation.getSign() * compatible == 1;
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

	public void setForceSiteMatching(boolean forceSiteMatching)
	{
		this.forceSiteMatching = forceSiteMatching;
	}
}
