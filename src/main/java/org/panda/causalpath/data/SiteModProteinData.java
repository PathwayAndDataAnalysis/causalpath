package org.panda.causalpath.data;

import org.panda.resource.siteeffect.Feature;
import org.panda.resource.tcga.ProteomicsFileRow;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Data specific to protein with a site modified.
 */
public class SiteModProteinData extends ProteinData
{
	/**
	 * Map from gene symbols to corresponding sites on proteins.
	 */
	protected Map<String, Set<ProteinSite>> siteMap;

	protected Feature mod;

	private Set<String> genesWithSites;

	Integer effect;

	public SiteModProteinData(String id, Set<String> geneSymbols, Feature mod)
	{
		super(id, geneSymbols);
		this.mod = mod;
	}

	/**
	 * Constructor to convert an RPPAData object from the "resource" project.
	 */
	public SiteModProteinData(ProteomicsFileRow row)
	{
		super(row);

		siteMap = new HashMap<>();

		try
		{
			row.sites.keySet().stream().forEach(sym ->
				siteMap.put(sym, row.sites.get(sym).stream().map(site ->
					new ProteinSite(Integer.parseInt(site.substring(1)), site.substring(0, 1),
						row.effect == null ? 0 : row.effect == ProteomicsFileRow.SiteEffect.ACTIVATING ? 1 :
							row.effect == ProteomicsFileRow.SiteEffect.INHIBITING ? -1 : 0)).collect(Collectors.toSet())));
		}
		catch (NumberFormatException e)
		{
			e.printStackTrace();
		}

		this.mod = row.mod;
	}

	/**
	 * Effect of this data depends on the effect of the phospho site.
	 */
	@Override
	public int getEffect()
	{
		if (effect == null)
		{
			boolean activating = siteMap.values().stream().flatMap(Collection::stream).anyMatch(site -> site.effect > 0);
			boolean inhibiting = siteMap.values().stream().flatMap(Collection::stream).anyMatch(site -> site.effect < 0);

			if (activating && !inhibiting) effect = 1;
			else if (!activating && inhibiting) effect = -1;
			else effect = 0;
		}
		return effect;
	}

	public Feature getModification()
	{
		return mod;
	}

	/**
	 * Gets string that shows genes with their sites in a set.
	 */
	public Set<String> getGenesWithSites()
	{
		if (this.genesWithSites == null)
		{
			this.genesWithSites = new HashSet<>();

			if (siteMap != null)
			{
				siteMap.keySet().stream().forEach(gene ->
					this.genesWithSites.addAll(
						siteMap.get(gene).stream().map(site -> gene + "_" + site.getSite())
							.collect(Collectors.toList())));
			}
		}

		return this.genesWithSites;
	}

	@Override
	public boolean isSiteSpecific()
	{
		return true;
	}

	@Override
	public ExperimentData copy()
	{
		SiteModProteinData copy = new SiteModProteinData(id, getGeneSymbols(), mod);
		copy.vals = vals;
		copy.setSiteMap(getSiteMap());
		return copy;
	}

	public Map<String, Set<ProteinSite>> getSiteMap()
	{
		return siteMap;
	}

	public void setSiteMap(Map<String, Set<ProteinSite>> siteMap)
	{
		this.siteMap = siteMap;
	}

	@Override
	public DataType getType()
	{
		switch (mod)
		{
			case PHOSPHORYLATION: return DataType.PHOSPHOPROTEIN;
			case ACETYLATION: return DataType.ACETYLPROTEIN;
			case METHYLATION: return DataType.METHYLPROTEIN;
		}
		return null;
	}
}
