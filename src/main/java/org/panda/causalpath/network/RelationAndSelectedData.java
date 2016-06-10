package org.panda.causalpath.network;

import org.panda.causalpath.data.ExperimentData;

/**
 * Created by babur on 4/18/16.
 */
public class RelationAndSelectedData
{
	public Relation relation;
	public ExperimentData source;
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
