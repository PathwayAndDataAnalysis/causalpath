package org.panda.causalpath.data;

import org.panda.causalpath.analyzer.OneDataChangeDetector;

import java.util.Collections;
import java.util.Set;

/**
 * Created by babur on 3/24/16.
 */
public abstract class ExperimentData
{
	public String id;
	protected Set<String> geneSymbols;

	public int getEffect()
	{
		return 1;
	}

	protected OneDataChangeDetector chDet;

	public ExperimentData(String id, Set<String> geneSymbols)
	{
		this.id = id;
		this.geneSymbols = geneSymbols;
	}

	public ExperimentData(String id, String geneSymbol)
	{
		this(id, Collections.singleton(geneSymbol));
	}

	public int getChangeSign()
	{
		if (chDet == null) throw new RuntimeException("getChangeSign can be called only after setting the change " +
			"detector.");

		return chDet.getChangeSign(this);
	}

	public double getChangeValue()
	{
		if (chDet == null) throw new RuntimeException("getChangeValue can be called only after setting the change " +
			"detector.");

		return chDet.getChangeValue(this);
	}

	public Set<String> getGeneSymbols()
	{
		return geneSymbols;
	}

	public void setChDet(OneDataChangeDetector chDet)
	{
		this.chDet = chDet;
	}

	public boolean hasChangeDetector()
	{
		return this.chDet != null;
	}

	@Override
	public int hashCode()
	{
		return id.hashCode();
	}

	@Override
	public boolean equals(Object obj)
	{
		return obj instanceof ExperimentData && ((ExperimentData) obj).id.equals(id);
	}
}
