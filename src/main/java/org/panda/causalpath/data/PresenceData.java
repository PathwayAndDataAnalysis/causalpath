package org.panda.causalpath.data;

import java.util.HashSet;
import java.util.Set;

/**
 * Series of copy number alterations for a gene.
 */
public class PresenceData extends CategoricalData
{
	/**
	 * The type of the original data.
	 */
	DataType type;

	public PresenceData(String id, Set<String> symbols, DataType type)
	{
		super(id, symbols);
		this.type = type;
	}

	@Override
	public ExperimentData copy()
	{
		PresenceData copy = new PresenceData(id, new HashSet<>(getGeneSymbols()), type);
		copy.data = data;
		return copy;
	}

	@Override
	public DataType getType()
	{
		return type;
	}
}
