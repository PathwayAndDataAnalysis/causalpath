package org.panda.causalpath.data;

import org.panda.resource.tcga.RPPAData;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by babur on 3/24/16.
 */
public class ProteinData extends NumericData
{
	public ProteinData(String id, Set<String> geneSymbols)
	{
		super(id, geneSymbols);
	}

	public ProteinData(RPPAData rppa)
	{
		this(rppa.id, new HashSet<>(rppa.genes));
		this.vals = rppa.vals;
	}
}
