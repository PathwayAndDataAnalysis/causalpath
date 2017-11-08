package org.panda.causalpath.network;

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A graph filter is for selecting a subset of interest from the result graph.
 *
 * @author Ozgun Babur
 */
public class GraphFilter
{
	/**
	 * Relation filter type to use.
	 */
	private RelationFilterType relationFilterType;

	/**
	 * A set of genes to trim the network to their neighborhood.
	 */
	private Set<String> focusGenes;

	/**
	 * Constructor with the type of the relation filter.
	 */
	public GraphFilter(RelationFilterType relationFilterType)
	{
		if (relationFilterType == null)
		{
			System.err.println("Graph filter type is null. No filter will be applied.");
			this.relationFilterType = RelationFilterType.NO_FILTER;
		}
		else this.relationFilterType = relationFilterType;
	}

	/**
	 * Constructor with the name of the relation filter.
	 */
	public GraphFilter(String relationFilterString)
	{
		this(RelationFilterType.get(relationFilterString));
	}

	/**
	 * Constructor with focus genes to crop the graph.
	 */
	public GraphFilter(Set<String> focusGenes)
	{
		this.focusGenes = focusGenes;
		this.relationFilterType = RelationFilterType.NO_FILTER;
	}

	/**
	 * Trims the network either to specific relation types or to the neighborhood of specific genes or both.
	 * @param results the network
	 * @return trimmed network
	 */
	public Set<Relation> filter(Set<Relation> results)
	{
		// Trim with relation type

		switch (relationFilterType)
		{
			case NO_FILTER:
				break;
			case PHOSPHO_ONLY: results = results.stream().filter(r -> r.type.affectsPhosphoSite).collect(Collectors.toSet());
				break;
			case EXPRESSION_ONLY: results = results.stream()
				.filter(r -> r.type.affectsTotalProt).collect(Collectors.toSet());
				break;
			case PHOSPHO_PRIMARY_EXPRESSION_SECONDARY:
			{
				Set<String> genes = getGenesInPhoshoGraph(results);
				results = results.stream().filter(r -> r.type.affectsPhosphoSite ||
					genes.contains(r.source) ||
					genes.contains(r.target))
					.collect(Collectors.toSet());
				break;
			}
			default: throw new RuntimeException("Unhandled relation filter type: " + relationFilterType);
		}

		// Trim with focus genes

		if (focusGenes != null && !focusGenes.isEmpty())
		{
			results = results.stream().filter(r ->
				focusGenes.contains(r.source) ||
				focusGenes.contains(r.target)).collect(Collectors.toSet());
		}

		return results;
	}

	public void setRelationFilterType(RelationFilterType relationFilterType)
	{
		this.relationFilterType = relationFilterType;
	}

	public void setFocusGenes(Set<String> focusGenes)
	{
		this.focusGenes = focusGenes;
	}

	/**
	 * Gets the set of IDs of the genes that has at least one phospho relation coming in or going out.
	 */
	private Set<String> getGenesInPhoshoGraph(Set<Relation> results)
	{
		return results.stream().filter(r -> r.type.affectsPhosphoSite).map(r -> new String[]{r.source, r.target})
			.flatMap(Arrays::stream).collect(Collectors.toSet());
	}

	/**
	 * Type of the relation type filtering.
	 */
	public enum RelationFilterType
	{
		NO_FILTER ("The graph is used with all inferred relations."),
		PHOSPHO_ONLY ("Only phosphorylation relations are desired."),
		EXPRESSION_ONLY ("Only expression relations are desired."),
		PHOSPHO_PRIMARY_EXPRESSION_SECONDARY ("All phosphorylation relations " +
			"are desired. Expression relations are desired only as supplemental, i.e., they have to involve at least " +
			"one protein that exists in phospho graph.");

		String description;

		RelationFilterType(String description)
		{
			this.description = description;
		}

		public String getName()
		{
			return toString().toLowerCase().replaceAll("_", "-");
		}

		public static RelationFilterType get(String s)
		{
			for (RelationFilterType relationFilterType : values())
			{
				if (relationFilterType.getName().equals(s)) return relationFilterType;
			}
			return null;
		}

		public static String getUsageInfo()
		{
			StringBuilder sb = new StringBuilder();
			for (RelationFilterType type : values())
			{
				sb.append("\t").append(type.getName()).append(": ").append(type.description).append("\n");
			}
			return sb.toString();
		}
	}
}
