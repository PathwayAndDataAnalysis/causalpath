package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationAndSelectedData;
import org.panda.utility.CollectionUtil;

import java.util.*;
import java.util.stream.Collectors;

/**
 * This class matches the experiment data with the pathway relations, and detect potential causality.
 */
public class CausalitySearcher
{
	/**
	 * If that is false, then we are interested in conflicting relations.
	 */
	private boolean causal = true;

	/**
	 * When there is a total protein measurement for a gene, we may prefer it to override its RNA expression data. This
	 * set is used to mark those proteins whose RNA expression will be ignored in the analysis.
	 */
	private Set<String> genesWithTotalProteinData;

	private RelationTargetCompatibilityChecker rtcc;

	public CausalitySearcher(RelationTargetCompatibilityChecker rtcc)
	{
		this.rtcc = rtcc;
	}

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
					if (!rtcc.isCompatible(source, relation, target)) continue;

					int chgSign = relation.chDet.getChangeSign(source, target);
					int sign = chgSign * source.getEffect() * relation.getSign();
					if (sign == compatible)
					{
						selected.add(new RelationAndSelectedData(relation, source, target));
					}
				}
			}
		}
		return selected;
	}

	/**
	 * Finds the sites with unknown effect, whose determination of effect can make the result graph larger.
	 */
	public Set<ExperimentData> findDataThatNeedsAnnotation(Set<Relation> relations, Set<RelationAndSelectedData> selected)
	{
		// Find the genes with a downstream in the results.
		Set<String> usedGene = selected.stream().map(r -> r.source.id).collect(Collectors.toSet());

		// The data we want to identify
		Set<ExperimentData> datas = new HashSet<>();

		for (Relation relation : relations)
		{
			if (usedGene.contains(relation.source)) continue;

			for (ExperimentData source : relation.sourceData)
			{
				if (skip(source, relation.source)) continue;

				for (ExperimentData target : relation.targetData)
				{
					if (skip(target, relation.target)) continue;
					if (!rtcc.isCompatible(source, relation, target)) continue;

					int chgSign = relation.chDet.getChangeSign(source, target) * relation.getSign();

					if (chgSign != 0 && source.getEffect() == 0)
					{
						datas.add(source);
					}
				}
			}
		}
		return datas;
	}

	/**
	 * Checks if the experiment data needs to be ignored because there is a total protein data.
	 */
	private boolean skip(ExperimentData data, String gene)
	{
		return genesWithTotalProteinData != null && genesWithTotalProteinData.contains(gene) &&
			(data instanceof ExpressionData || data instanceof CNAData);
	}

	public void setGenesWithTotalProteinData(Set<String> genesWithTotalProteinData)
	{
		this.genesWithTotalProteinData = genesWithTotalProteinData;
	}

	public void setCausal(boolean causal)
	{
		this.causal = causal;
	}
}
