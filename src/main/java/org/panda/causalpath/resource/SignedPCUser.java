package org.panda.causalpath.resource;

import org.panda.causalpath.data.PhosphoSite;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationType;
import org.panda.resource.PhosphoSitePlus;
import org.panda.resource.network.SignedPC;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.PhosphoGraph;

import java.util.*;

/**
 * Created by babur on 4/4/16.
 */
public class SignedPCUser
{
	public static Set<Relation> getSignedPCRelations()
	{
		Set<Relation> relations = new HashSet<>();
		SignedPC spc = new SignedPC();
		Map<SignedType, Graph> allGraphs = spc.getAllGraphs();
		for (SignedType signedType : allGraphs.keySet())
		{
			Graph graph = allGraphs.get(signedType);

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
