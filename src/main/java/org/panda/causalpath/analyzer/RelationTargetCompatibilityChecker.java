package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;
import org.panda.utility.CollectionUtil;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * This class checks if the pathway relation can have potential to explain the change in the target data. It also checks
 * if the source data is
 *
 * @author Ozgun Babur
 */
public class RelationTargetCompatibilityChecker
{
	/**
	 * Data types to explain activity changes.
	 */
	static Set<Class<? extends ExperimentData>> aSet = new HashSet<>(Arrays.asList(PhosphoProteinData.class, ProteinData.class, MutationData.class));

	/**
	 * Data types to explain phosphorylations.
	 */
	static Set<Class<? extends ExperimentData>> pSet = Collections.singleton(PhosphoProteinData.class);

	/**
	 * Data types to explain total protein changes.
	 */
	static Set<Class<? extends ExperimentData>> eSet = new HashSet<>(Arrays.asList(ProteinData.class, ExpressionData.class));
//	static Set<Class<? extends ExperimentData>> eSet = Collections.singleton(ProteinData.class);

	/**
	 * Parameter to mandate site matching to explain phosphorylations.
	 */
	private boolean forceSiteMatching = true;

	/**
	 * Sometimes a site with unknown effect, but also very close to another site with known effect, is very likely to
	 * have the same effect. This parameter is the site proximity threshold to infer the effect of sites with neighbor
	 * sites.
	 */
	private int siteProximityThreshold = 0;


	public boolean isCompatible(ExperimentData source, Relation rel, ExperimentData target)
	{
		// Source should have potential to explain activity change.
		if (source != null && !aSet.contains(source.getClass())) return false;

		// At least one data has to be proteomics measurement, or source is null
		if (source != null && !(source instanceof ProteinData) && !(target instanceof ProteinData)) return false;

		// The relation type should have potential to explain the change in target

		if (rel.type.affectsTotalProt)
		{
			return eSet.contains(target.getClass());
		}

		if (rel.type.affectsPhosphoSite)
		{
			if (target instanceof PhosphoProteinData)
			{
				return !forceSiteMatching || !CollectionUtil.intersectionEmpty(
					rel.getTargetWithSites(siteProximityThreshold), ((PhosphoProteinData) target).getGenesWithSites());
			}
		}

		return false;
	}

	public void setSiteProximityThreshold(int thr)
	{
		this.siteProximityThreshold = thr;
	}

	public void setForceSiteMatching(boolean forceSiteMatching)
	{
		this.forceSiteMatching = forceSiteMatching;
	}
}
