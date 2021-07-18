package org.panda.causalpath.data;

/**
 * A site of a protein subject to a modification.
 */
public class ProteinSite
{
	/**
	 * Location on the protein.
	 */
	protected int site;

	/**
	 * The effect on protein activity when it is phosphorylated.
	 */
	protected int effect;

	/**
	 * The one-letter amino acid code of the site.
	 */
	protected String aa;

	public ProteinSite(int site, String aa, int effect)
	{
		this.site = site;
		this.aa = aa;
		this.effect = effect;
	}

	public int getSite()
	{
		return site;
	}

	@Override
	public String toString()
	{
		return aa + site;
	}

	@Override
	public int hashCode()
	{
		return (aa + site).hashCode();
	}

	@Override
	public boolean equals(Object obj)
	{
		return obj instanceof ProteinSite && aa.equals(((ProteinSite) obj).aa) && site == ((ProteinSite) obj).site;
	}
}
