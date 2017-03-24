package org.panda.causalpath.resource;

import org.panda.causalpath.data.PhosphoSite;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationType;
import org.panda.resource.PhosphoSitePlus;
import org.panda.resource.network.PhosphoNetworks;
import org.panda.resource.network.SignedPC;
import org.panda.resource.network.SignedREACH;
import org.panda.resource.network.TRRUST;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.PhosphoGraph;

import java.util.*;

/**
 * Provides the signed relations in Pathway Commons as a set of Relation objects.
 */
public class NetworkLoader
{
	/**
	 * Reads 3 of the built-in resource networks.
	 */
	public static Set<Relation> load()
	{
		return load(new HashSet<>(Arrays.asList(ResourceType.PC, ResourceType.PhosphoNetworks, ResourceType.TRRUST)));
	}

	/**
	 * Reads the selected built-in resource networks.
	 */
	public static Set<Relation> load(Set<ResourceType> resourceTypes)
	{
		Set<Relation> relations = new HashSet<>();

		Map<SignedType, Graph> allGraphs = new HashMap<>();

		// Load signed directed graph from Pathway Commons
		if (resourceTypes.contains(ResourceType.PC))
		{
			mergeSecondMapIntoFirst(allGraphs, SignedPC.get().getAllGraphs());
		}

		// Add REACH
		if (resourceTypes.contains(ResourceType.REACH))
		{
			mergeSecondMapIntoFirst(allGraphs, SignedREACH.get().getAllGraphs());
		}

		// Add PhosphoNetworks
		if (resourceTypes.contains(ResourceType.PhosphoNetworks))
		{
			allGraphs.get(SignedType.PHOSPHORYLATES).merge(PhosphoNetworks.get().getGraph());
		}

		// Add TRRUST
		if (resourceTypes.contains(ResourceType.TRRUST))
		{
			allGraphs.get(SignedType.UPREGULATES_EXPRESSION).merge(TRRUST.get().getPositiveGraph());
			allGraphs.get(SignedType.DOWNREGULATES_EXPRESSION).merge(TRRUST.get().getNegativeGraph());
		}

		// Generate relations based on the network

		for (SignedType signedType : allGraphs.keySet())
		{
			Graph graph = allGraphs.get(signedType);

			// take a subset of the network for debugging
//			graph.crop(Arrays.asList("PIK3CA", "AKT1"));

			for (String source : graph.getOneSideSymbols(true))
			{
				for (String target : graph.getDownstream(source))
				{
					RelationType type;
					switch (signedType)
					{
						case PHOSPHORYLATES: type = RelationType.PHOSPHORYLATES; break;
						case DEPHOSPHORYLATES: type = RelationType.DEPHOSPHORYLATES; break;
						case UPREGULATES_EXPRESSION: type = RelationType.UPREGULATES_EXPRESSION; break;
						case DOWNREGULATES_EXPRESSION: type = RelationType.DOWNREGULATES_EXPRESSION; break;
						default: throw new RuntimeException("Is there a new type??");
					}
					Relation rel = new Relation(source, target, type, graph.getMediatorsInString(source, target));

					if (graph instanceof PhosphoGraph)
					{
						PhosphoGraph pGraph = (PhosphoGraph) graph;
						Set<String> sites = pGraph.getSites(source, target);
						if (sites != null && !sites.isEmpty())
						{
							rel.sites = new HashSet<>();

							for (String site : sites)
							{
								Integer effect = PhosphoSitePlus.get().getEffect(target, site);
								rel.sites.add(new PhosphoSite(Integer.parseInt(site.substring(1)),
									site.substring(0, 1), effect == null ? 0 : effect));
							}
						}
					}

					relations.add(rel);
				}
			}
		}
		return relations;
	}

	private static void mergeSecondMapIntoFirst(Map<SignedType, Graph> allGraphs, Map<SignedType, Graph> graphs)
	{
		for (SignedType type : graphs.keySet())
		{
			if (allGraphs.containsKey(type)) allGraphs.get(type).merge(graphs.get(type));
			else allGraphs.put(type, graphs.get(type));
		}
	}

	public enum ResourceType
	{
		PC,
		REACH,
		PhosphoNetworks,
		TRRUST;

		public static Set<ResourceType> getSelectedResources(String s)
		{
			Set<ResourceType> set = new HashSet<>();
			for (String res : s.split(",|;|\\s+|\\|"))
			{
				res = res.trim();
				try
				{
					ResourceType type = valueOf(res);
					set.add(type);
				}
				catch (IllegalArgumentException e)
				{
					throw new RuntimeException("Network resource not recognized: " + res);
				}
			}
			return set;
		}
	}
}
