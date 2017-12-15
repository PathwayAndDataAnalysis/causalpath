package org.panda.causalpath.resource;

import org.panda.causalpath.data.GeneWithData;
import org.panda.causalpath.data.PhosphoSite;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationType;
import org.panda.resource.network.*;
import org.panda.resource.signednetwork.SignedType;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.utility.CollectionUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.PhosphoGraph;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Provides the signed relations in Pathway Commons as a set of Relation objects.
 */
public class NetworkLoader
{
	/**
	 * Reads 4 of the built-in resource networks.
	 */
	public static Set<Relation> load()
	{
		return load(new HashSet<>(Arrays.asList(ResourceType.PC, ResourceType.PhosphoNetworks, ResourceType.IPTMNet,
			ResourceType.TRRUST, ResourceType.TFactS)));
	}

	/**
	 * Reads the selected built-in resource networks.
	 */
	public static Set<Relation> load(Set<ResourceType> resourceTypes)
	{
		Set<Relation> relations = new HashSet<>();

		Map<SignedType, DirectedGraph> allGraphs = new HashMap<>();

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

		// Add IPTMNet
		if (resourceTypes.contains(ResourceType.IPTMNet))
		{
			addGraph(allGraphs, SignedType.PHOSPHORYLATES, IPTMNet.get().getGraph());
		}

		// Add PhosphoNetworks
		if (resourceTypes.contains(ResourceType.PhosphoNetworks))
		{
			addGraph(allGraphs, SignedType.PHOSPHORYLATES, PhosphoNetworks.get().getGraph());
		}

		// Add TRRUST
		if (resourceTypes.contains(ResourceType.TRRUST))
		{
			addGraph(allGraphs, SignedType.UPREGULATES_EXPRESSION, TRRUST.get().getPositiveGraph());
			addGraph(allGraphs, SignedType.DOWNREGULATES_EXPRESSION, TRRUST.get().getNegativeGraph());
		}

		// Add TFactS
		if (resourceTypes.contains(ResourceType.TFactS))
		{
			addGraph(allGraphs, SignedType.UPREGULATES_EXPRESSION, TFactS.get().getPositiveGraph());
			addGraph(allGraphs, SignedType.DOWNREGULATES_EXPRESSION, TFactS.get().getNegativeGraph());
		}

		// Generate relations based on the network

		for (SignedType signedType : allGraphs.keySet())
		{
			DirectedGraph graph = allGraphs.get(signedType);

			// take a subset of the network for debugging
//			graph.crop(Arrays.asList("HCK", "CD247"));

			SiteEffectCollective sec = new SiteEffectCollective();

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
						case ACTIVATES_GTPASE: type = RelationType.ACTIVATES_GTPASE; break;
						case INHIBITS_GTPASE: type = RelationType.INHIBITS_GTPASE; break;
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
								Integer effect = sec.getEffect(target, site);
								rel.sites.add(new PhosphoSite(Integer.parseInt(site.substring(1)),
									site.substring(0, 1), effect == null ? 0 : effect));
							}
						}
					}

					relations.add(rel);
				}
			}
		}

		// initiate source and target data

		Set<String> genes = relations.stream().map(r -> r.source).collect(Collectors.toSet());
		genes.addAll(relations.stream().map(r -> r.target).collect(Collectors.toSet()));

		Map<String, GeneWithData> map = new HashMap<>();
		genes.forEach(g -> map.put(g, new GeneWithData(g)));

		for (Relation relation : relations)
		{
			relation.sourceData = map.get(relation.source);
			relation.targetData = map.get(relation.target);
		}

		return relations;
	}

	private static void mergeSecondMapIntoFirst(Map<SignedType, DirectedGraph> allGraphs, Map<SignedType, DirectedGraph> graphs)
	{
		for (SignedType type : graphs.keySet())
		{
			if (allGraphs.containsKey(type)) allGraphs.get(type).merge(graphs.get(type));
			else allGraphs.put(type, graphs.get(type));
		}
	}

	private static void addGraph(Map<SignedType, DirectedGraph> allGraphs, SignedType type, DirectedGraph graph)
	{
		if (allGraphs.containsKey(type)) allGraphs.get(type).merge(graph);
		else allGraphs.put(type, graph);
	}

	public enum ResourceType
	{
		PC("Pathway Commons v9 for all kinds of relations."),
		REACH("Network derived from REACH NLP extraction results for phosphorylation relations."),
		PhosphoNetworks("The PhosphoNetworks database for phosphorylations."),
		IPTMNet("The IPTMNet database for phosphorylations."),
		TRRUST("The TRRUST database for expression relations."),
		TFactS("The TFactS database for expression relations.");

		String description;

		ResourceType(String description)
		{
			this.description = description;
		}

		public static Set<ResourceType> getSelectedResources(String s)
		{
			Set<ResourceType> set = new HashSet<>();
			for (String res : s.split(",|;|\\s+|\\|"))
			{
				res = res.trim();

				if (!res.isEmpty())
				{
					try
					{
						ResourceType type = valueOf(res);
						set.add(type);
					} catch (IllegalArgumentException e)
					{
						throw new RuntimeException("Network resource not recognized: " + res);
					}
				}
			}
			return set;
		}

		public static String getUsageInfo()
		{
			StringBuilder sb = new StringBuilder();
			for (ResourceType type : values())
			{
				sb.append("\t").append(type.name()).append(": ").append(type.description).append("\n");
			}
			return sb.toString();
		}
	}


	public static void main(String[] args) throws IOException
	{
		Set<Relation> rels = NetworkLoader.load(new HashSet<>(Arrays.asList(ResourceType.PC, ResourceType.PhosphoNetworks)));
		Map<String, Set<PhosphoSite>> map = new HashMap<>();
		for (Relation rel : rels)
		{
			if (rel.sites != null)
			{
				if (!map.containsKey(rel.target)) map.put(rel.target, new HashSet<>());
				map.get(rel.target).addAll(rel.sites);
			}
		}
		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/babur/Documents/Temp/gene-to-sites.txt"));

		for (String gene : map.keySet())
		{
			Set<PhosphoSite> sites = map.get(gene);
			writer.write("\n" + gene + "\t");
			writer.write(CollectionUtil.merge(sites, " "));
		}

		writer.close();
	}
}
