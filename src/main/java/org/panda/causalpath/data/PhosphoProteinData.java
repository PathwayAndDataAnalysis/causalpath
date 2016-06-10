package org.panda.causalpath.data;

import org.panda.resource.tcga.RPPAData;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by babur on 3/24/16.
 */
public class PhosphoProteinData extends ProteinData
{
	protected Map<String, Set<PhosphoSite>> siteMap;

	public PhosphoProteinData(String id, Set<String> geneSymbols)
	{
		super(id, geneSymbols);
	}

	public PhosphoProteinData(RPPAData rppa)
	{
		super(rppa);

		siteMap = new HashMap<>();

		rppa.sites.keySet().stream().forEach(sym ->
			siteMap.put(sym, rppa.sites.get(sym).stream().map(site ->
				new PhosphoSite(Integer.parseInt(site.substring(1)), site.substring(0, 1),
					rppa.effect == null ? 0 : rppa.effect == RPPAData.SiteEffect.ACTIVATING ? 1 :
					rppa.effect == RPPAData.SiteEffect.INHIBITING ? -1 : 0)).collect(Collectors.toSet())));
	}

	@Override
	public int getEffect()
	{
		boolean activating = siteMap.values().stream().flatMap(Collection::stream).anyMatch(site -> site.effect > 0);
		boolean inhibiting = siteMap.values().stream().flatMap(Collection::stream).anyMatch(site -> site.effect < 0);

		if (activating && !inhibiting) return 1;
		if (!activating && inhibiting) return -1;
		return 0;
	}

	public Set<String> getGenesWithSites()
	{
		if (siteMap == null) return Collections.emptySet();

		Set<String> set = new HashSet<>();
		siteMap.keySet().stream().forEach(gene ->
			set.addAll(siteMap.get(gene).stream().map(site -> gene + "_" + site).collect(Collectors.toList())));

		return set;
	}
}
