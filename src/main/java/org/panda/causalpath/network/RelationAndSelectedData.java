package org.panda.causalpath.network;

import org.panda.causalpath.data.ExperimentData;

/**
 * This is a causative relation. It holds a pathway relation and two experiment data whose dependency can be explained
 * with that relation.
 */
public class RelationAndSelectedData
{
	/**
	 * The edge.
	 */
	public Relation relation;

	/**
	 * Source data.
	 */
	public ExperimentData source;

	/**
	 * Target data.
	 */
	public ExperimentData target;

	public RelationAndSelectedData(Relation relation, ExperimentData source, ExperimentData target)
	{
		this.relation = relation;
		this.source = source;
		this.target = target;
	}

	@Override
	public int hashCode()
	{
		return relation.hashCode() + source.hashCode() + target.hashCode();
	}

	@Override
	public boolean equals(Object obj)
	{
		if (obj instanceof RelationAndSelectedData)
		{
			RelationAndSelectedData r = (RelationAndSelectedData) obj;
			return r.relation.equals(relation) && r.source.equals(source) && r.target.equals(target);
		}
		else return false;
	}
}
