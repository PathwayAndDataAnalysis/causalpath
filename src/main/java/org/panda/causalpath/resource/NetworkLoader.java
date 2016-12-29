package org.panda.causalpath.resource;

import org.panda.causalpath.data.PhosphoSite;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationType;
import org.panda.resource.PhosphoSitePlus;
import org.panda.resource.network.PhosphoNetworks;
import org.panda.resource.network.SignedPC;
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
	 * Reads Pathway Commons from the SignedPC in the resource project.
	 */
	public static Set<Relation> load()
	{
		Set<Relation> relations = new HashSet<>();

		// Load signed directed graph from Pathway Commons

		SignedPC spc = new SignedPC();
		Map<SignedType, Graph> allGraphs = spc.getAllGraphs();

		// Add PhoshoNetworks and TRRUST

		allGraphs.get(SignedType.PHOSPHORYLATES).merge(PhosphoNetworks.get().getGraph());
		allGraphs.get(SignedType.UPREGULATES_EXPRESSION).merge(TRRUST.get().getPositiveGraph());
		allGraphs.get(SignedType.DOWNREGULATES_EXPRESSION).merge(TRRUST.get().getNegativeGraph());

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
}
