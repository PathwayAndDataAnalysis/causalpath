package org.panda.causalpath.data;

import org.panda.causalpath.analyzer.OneDataChangeDetector;

import java.util.Collections;
import java.util.Set;

/**
 * Base class for a series of experiment data for a gene.
 */
public abstract class ExperimentData
{
	/**
	 * ID of the measurement. This can be an identifier of the antibody in an RPPA experiment, or gene symbol for a
	 * mutation, or an identifier for a mass spectometry peak.
	 */
	public String id;

	/**
	 * List of gene symbols that the measurement applies.
	 */
	protected Set<String> geneSymbols;

	/**
	 * How does a positive value in this experiment data affects the gene's activity? Only negative example is
	 * methylation so far.
	 */
	public int getEffect()
	{
		return 1;
	}

	/**
	 * A plug-in detector for the change in the experiment data. This detector may compare a test and a control group,
	 * or can decide based on a threshold, etc.
	 */
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

	public OneDataChangeDetector getChDet()
	{
		return chDet;
	}

	public String getId()
	{
		return id;
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

	/**
	 * An experiment data is not site specific by default.
	 */
	public boolean isSiteSpecific()
	{
		return false;
	}

	/**
	 * Generates a copy of the object.
	 */
	public abstract ExperimentData copy();

	public abstract DataType getType();

	@Override
	public String toString()
	{
		return id;
	}
}
