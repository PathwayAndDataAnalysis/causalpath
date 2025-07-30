package org.panda.causalpath.resource;

import org.panda.causalpath.data.GeneWithData;
import org.panda.causalpath.data.ProteinSite;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationType;
import org.panda.resource.HGNC;
import org.panda.resource.network.*;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.CollectionUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.SiteSpecificGraph;

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
		return load(new HashSet<>(Arrays.asList(ResourceType.PC, ResourceType.PhosphoSitePlus,
			ResourceType.PhosphoNetworks, ResourceType.IPTMNet, ResourceType.TRRUST, ResourceType.TFactS, ResourceType.PCMetabolic)));
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

		// Add PhosphoSitePlus
		if (resourceTypes.contains(ResourceType.PhosphoSitePlus))
		{
			mergeSecondMapIntoFirst(allGraphs, PhosphoSitePlus.get().getAllGraphs());
		}

		// Add REACH
		if (resourceTypes.contains(ResourceType.REACH))
		{
			mergeSecondMapIntoFirst(allGraphs, SignedREACH.get().getAllGraphs());
		}

		// Add IPTMNet
		if (resourceTypes.contains(ResourceType.IPTMNet))
		{
			mergeSecondMapIntoFirst(allGraphs, IPTMNet.get().getAllGraphs());
		}

		// Add PhosphoNetworks
		if (resourceTypes.contains(ResourceType.PhosphoNetworks))
		{
			addGraph(allGraphs, SignedType.PHOSPHORYLATES, PhosphoNetworks.get().getGraph());
		}

		// Add NetworKIN
		if (resourceTypes.contains(ResourceType.NetworKIN))
		{
			addGraph(allGraphs, SignedType.PHOSPHORYLATES, NetworKIN.get().getGraph());
		}

		// Add Rho GEF
		if (resourceTypes.contains(ResourceType.RHOGEF))
		{
			addGraph(allGraphs, SignedType.ACTIVATES_GTPASE, RhoGEF.get().getGraph());
		}

		// Experimental code!!!!!!!!!!!!!!!!!!!
		if (resourceTypes.contains(ResourceType.PCTCGAConsensus))
		{
			DirectedGraph posG = new DirectedGraph("Pos", SignedType.UPREGULATES_EXPRESSION.getTag());
			DirectedGraph negG = new DirectedGraph("Neg", SignedType.DOWNREGULATES_EXPRESSION.getTag());
			String file = "/home/babur/Documents/PC/SignedByTCGAConsensusFiltered.sif";
			posG.load(file, Collections.singleton(posG.getEdgeType()));
			negG.load(file, Collections.singleton(negG.getEdgeType()));
			addGraph(allGraphs, SignedType.UPREGULATES_EXPRESSION, posG);
			addGraph(allGraphs, SignedType.DOWNREGULATES_EXPRESSION, negG);
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

		// Add PC Metabolic
		if (resourceTypes.contains(ResourceType.PCMetabolic))
		{
			addGraph(allGraphs, SignedType.PRODUCES, SignedMetabolic.getGraph(SignedType.PRODUCES));
			addGraph(allGraphs, SignedType.CONSUMES, SignedMetabolic.getGraph(SignedType.CONSUMES));
			addGraph(allGraphs, SignedType.USED_TO_PRODUCE, SignedMetabolic.getGraph(SignedType.USED_TO_PRODUCE));
		}

		// Clean-up conflicts
		cleanUpConflicts(allGraphs);

		// Generate relations based on the network
		relations = addGraphsToRelations(allGraphs, relations);

		return relations;
	}

	public static Set<Relation> addGraphsToRelations(Map<SignedType, DirectedGraph> graphs, Set<Relation> relations)
	{
		for (SignedType signedType : graphs.keySet())
		{
			DirectedGraph graph = graphs.get(signedType);
			if (graph == null)
			{
//				System.out.println("Null graph for type: " + signedType);
				continue;
			}

			// take a subset of the network for debugging
//			graph.crop(Arrays.asList("HCK", "CD247"));

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
						case ACETYLATES: type = RelationType.ACETYLATES; break;
						case DEACETYLATES: type = RelationType.DEACETYLATES; break;
						case DEMETHYLATES: type = RelationType.DEMETHYLATES; break;
						case METHYLATES: type = RelationType.METHYLATES; break;
						case PRODUCES: type = RelationType.PRODUCES; break;
						case CONSUMES: type = RelationType.CONSUMES; break;
						case USED_TO_PRODUCE: type = RelationType.USED_TO_PRODUCE; break;
						default: throw new RuntimeException("Is there a new type??");
					}
					Relation rel = new Relation(source, target, type, graph.getMediatorsInString(source, target));

					if (graph instanceof SiteSpecificGraph)
					{
						SiteSpecificGraph pGraph = (SiteSpecificGraph) graph;
						Set<String> sites = pGraph.getSites(source, target);
						if (sites != null && !sites.isEmpty())
						{
							rel.sites = new HashSet<>();

							for (String site : sites)
							{
								if (!site.isEmpty())
								{
									rel.sites.add(new ProteinSite(Integer.parseInt(site.substring(1)),
										site.substring(0, 1), 0));
								}
							}
						}
					}

					relations.add(rel);
				}
			}
		}

		// clean from non-HGNC

		relations = relations.stream().filter(r -> (r.source.startsWith("CHEBI:") || HGNC.get().getSymbol(r.source) != null) &&
			(r.target.startsWith("CHEBI:") || HGNC.get().getSymbol(r.target) != null)).collect(Collectors.toSet());

		// initiate source and target data

		initMissingSourceTargetData(relations);

		return relations;
	}

	private static void check(Map<SignedType, DirectedGraph> allGraphs)
	{
		DirectedGraph graph = allGraphs.get(SignedType.DEPHOSPHORYLATES);
		if (graph != null)
		{
			Set<String> dwns = graph.getDownstream("ABL1");
			if (dwns.contains("RB1"))
			{
				System.out.println("Here it is!!!");
			}
			else System.out.println("Nope");
		}
		else System.out.println("No graph");
	}

	public static Set<Relation> load(String customFile) throws IOException
	{
		return load(customFile, null);
	}

	public static Set<Relation> load(String customFile, Set<Relation> relations) throws IOException
	{
		if (relations == null) relations = new HashSet<>();

		Files.lines(Paths.get(customFile)).filter(l -> !l.isEmpty() && !l.startsWith("#")).map(Relation::new)
			.forEach(relations::add);

		initMissingSourceTargetData(relations);

		return relations;
	}

	private static void initMissingSourceTargetData(Set<Relation> relations)
	{
		Map<String, GeneWithData> map = new HashMap<>();

		for (Relation relation : relations)
		{
			if (relation.sourceData != null)
			{
				if (!map.containsKey(relation.source)) map.put(relation.source, relation.sourceData);
				if (!map.containsKey(relation.target)) map.put(relation.target, relation.targetData);
			}
		}

		for (Relation relation : relations)
		{
			if (relation.sourceData == null)
			{
				relation.sourceData = map.getOrDefault(relation.source, new GeneWithData(relation.source));
				if (!map.containsKey(relation.source)) map.put(relation.source, relation.sourceData);
				relation.targetData = map.getOrDefault(relation.target, new GeneWithData(relation.target));
				if (!map.containsKey(relation.target)) map.put(relation.target, relation.targetData);
			}
		}
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

	private static void cleanUpConflicts(Map<SignedType, DirectedGraph> graphs)
	{
		cleanUpConflicts((SiteSpecificGraph) graphs.get(SignedType.PHOSPHORYLATES), (SiteSpecificGraph) graphs.get(SignedType.DEPHOSPHORYLATES));
		cleanUpConflicts((SiteSpecificGraph) graphs.get(SignedType.ACETYLATES), (SiteSpecificGraph) graphs.get(SignedType.DEACETYLATES));
		cleanUpConflicts((SiteSpecificGraph) graphs.get(SignedType.METHYLATES), (SiteSpecificGraph) graphs.get(SignedType.DEMETHYLATES));
		cleanUpConflicts(graphs.get(SignedType.UPREGULATES_EXPRESSION), graphs.get(SignedType.DOWNREGULATES_EXPRESSION));
		cleanUpConflicts(graphs.get(SignedType.ACTIVATES_GTPASE), graphs.get(SignedType.INHIBITS_GTPASE));
	}

	private static void cleanUpConflicts(DirectedGraph graph1, DirectedGraph graph2)
	{
		if (graph1 == null || graph2 == null) return;

		Set<String[]> rem1 = new HashSet<>();
		Set<String[]> rem2 = new HashSet<>();

		for (String source : graph1.getOneSideSymbols(true))
		{
			for (String target : graph1.getDownstream(source))
			{
				if (graph2.hasRelation(source, target))
				{
					String[] rel = new String[]{source, target};
					long c1 = paperCnt(graph1.getMediators(source, target));
					long c2 = paperCnt(graph2.getMediators(source, target));

					if (c1 >= c2) rem2.add(rel);
					if (c2 >= c1) rem1.add(rel);
				}
			}
		}

		rem1.forEach(rel -> graph1.removeRelation(rel[0], rel[1]));
		rem2.forEach(rel -> graph2.removeRelation(rel[0], rel[1]));
	}

	private static void cleanUpConflicts(SiteSpecificGraph graph1, SiteSpecificGraph graph2)
	{
		if (graph1 == null || graph2 == null) return;

		for (String source : graph1.getOneSideSymbols(true))
		{
			for (String target : graph1.getDownstream(source))
			{
				if (graph2.hasRelation(source, target))
				{
					// For debugging the relations
//					if (source.equals("IL6ST") && target.equals("STAT3"))
//					{
//						System.out.println();
//					}

					Set<String> s1 = graph1.getSites(source, target);
					Set<String> s2 = graph2.getSites(source, target);

					Set<String> common = CollectionUtil.getIntersection(s1, s2);

					if (!common.isEmpty())
					{
						long c1 = paperCnt(graph1.getMediators(source, target));
						long c2 = paperCnt(graph2.getMediators(source, target));

						for (String site : common)
						{
							if (c1 >= c2) graph2.removeSite(source, target, site);
							if (c2 >= c1) graph1.removeSite(source, target, site);
						}
					}
				}
			}
		}
	}

	private static long paperCnt(Set<String> mediators)
	{
		return mediators.stream().filter(s -> s.startsWith("PMID:")).distinct().count();
	}

	public enum ResourceType
	{
		PC("Pathway Commons v9 for all kinds of relations."),
		PhosphoSitePlus("PhosphoSitePlus database for (de)phosphorylations, (de)acetylations and (de)methylations."),
		REACH("Network derived from REACH NLP extraction results for phosphorylation relations."),
		PhosphoNetworks("The PhosphoNetworks database for phosphorylations."),
		IPTMNet("The IPTMNet database for phosphorylations."),
		RHOGEF("The experimental Rho - GEF relations."),
		PCTCGAConsensus("Unsigned PC relations whose signs are inferred by TCGA studies."),
		TRRUST("The TRRUST database for expression relations."),
		TFactS("The TFactS database for expression relations."),
		NetworKIN("The NetworKIN database for phosphorylation relations."),
		PCMetabolic("Relations involving chemicals in PC"),
		;

		String description;

		ResourceType(String description)
		{
			this.description = description;
		}

		public static Set<ResourceType> getSelectedResources(Set<String> names)
		{
			Set<ResourceType> set = new HashSet<>();
			for (String res : names)
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

		public static Map getValuesAsJson()
		{
			List list = new ArrayList<>();
			for (ResourceType type : values())
			{
				list.add(type.toString());
			}

			Map map = new LinkedHashMap<>();
			map.put("name", "ResourceType");
			map.put("values", list);
			return map;
		}

	}


	public static void main(String[] args) throws IOException
	{
//		writeSitesWithUpstream();
		writeRelations();
	}

	private static void writeRelations() throws IOException
	{
		Set<Relation> rels = NetworkLoader.load();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/ozgunbabur/Documents/causal-priors.txt"));

		for (Relation rel : rels)
		{
			writer.write(rel.toString() + "\n");
		}

		writer.close();
	}

	public static void writeSitesWithUpstream() throws IOException
	{
		Set<Relation> rels = NetworkLoader.load(new HashSet<>(Arrays.asList(ResourceType.PC, ResourceType.PhosphoNetworks)));
		Map<String, Set<ProteinSite>> map = new HashMap<>();
		for (Relation rel : rels)
		{
			if (rel.sites != null)
			{
				if (!map.containsKey(rel.target)) map.put(rel.target, new HashSet<>());
				map.get(rel.target).addAll(rel.sites);
			}
		}
		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/ozgunbabur/Documents/Temp/gene-to-sites.txt"));

		for (String gene : map.keySet())
		{
			Set<ProteinSite> sites = map.get(gene);
			writer.write("\n" + gene + "\t");
			writer.write(CollectionUtil.merge(sites, " "));
		}

		writer.close();
	}

}
