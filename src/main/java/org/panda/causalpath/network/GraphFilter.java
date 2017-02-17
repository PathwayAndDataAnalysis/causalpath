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

	public Set<RelationAndSelectedData> filter(Set<RelationAndSelectedData> results)
	{
		switch (relationFilterType)
		{
			case NO_FILTER:
			case PHOSPHO_ONLY: results = results.stream().filter(r -> r.relation.type == RelationType.PHOSPHORYLATES ||
				r.relation.type == RelationType.DEPHOSPHORYLATES).collect(Collectors.toSet());
				break;
			case EXPRESSION_ONLY: results = results.stream()
				.filter(r -> r.relation.type == RelationType.UPREGULATES_EXPRESSION ||
				r.relation.type == RelationType.DOWNREGULATES_EXPRESSION).collect(Collectors.toSet());
				break;
			case PHOSPHO_PRIMARY_EXPRESSION_SECONDARY:
			{
				Set<String> genes = getGenesInPhoshoGraph(results);
				results = results.stream().filter(r -> r.relation.type == RelationType.PHOSPHORYLATES ||
					r.relation.type == RelationType.DEPHOSPHORYLATES ||
					genes.contains(r.source.id) ||
					genes.contains(r.target.id))
					.collect(Collectors.toSet());
				break;
			}
			default: throw new RuntimeException("Unhandled relation filter type: " + relationFilterType);
		}

		if (focusGenes != null && !focusGenes.isEmpty())
		{
			results = results.stream().filter(r ->
				focusGenes.contains(r.relation.source) ||
				focusGenes.contains(r.relation.target)).collect(Collectors.toSet());
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
	private Set<String> getGenesInPhoshoGraph(Set<RelationAndSelectedData> results)
	{
		return results.stream().filter(r -> r.relation.type == RelationType.PHOSPHORYLATES ||
			r.relation.type == RelationType.DEPHOSPHORYLATES)
			.map(r -> new String[]{r.source.id, r.target.id}).flatMap(Arrays::stream).collect(Collectors.toSet());
	}


	public enum RelationFilterType
	{
		NO_FILTER ("no-filter", "The graph is used with all inferred relations."),
		PHOSPHO_ONLY ("phospho-only", "Only phosphorylation relations are desired."),
		EXPRESSION_ONLY ("expression-only", "Only expression relations are desired."),
		PHOSPHO_PRIMARY_EXPRESSION_SECONDARY ("phospho-primary-expression-secondary", "All phosphorylation relations " +
			"are desired. Expression relations are desired only as supplemental, i.e., they have to involve at least " +
			"one protein that exists in phospho graph.");

		String name;
		String description;

		RelationFilterType(String name, String description)
		{
			this.name = name;
			this.description = description;
		}

		public static RelationFilterType get(String s)
		{
			for (RelationFilterType relationFilterType : values())
			{
				if (relationFilterType.name.equals(s)) return relationFilterType;
			}
			return null;
		}
	}
}
