package org.panda.causalpath.data;

/**
 * Created by babur on 3/24/16.
 */
public class PhosphoSite
{
	protected int site;
	protected int effect;

	/**
	 * The one-letter amino acid code of the site.
	 */
	protected String aa;

	public PhosphoSite(int site, String aa, int effect)
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
		return obj instanceof PhosphoSite && aa.equals(((PhosphoSite) obj).aa) && site == ((PhosphoSite) obj).site;
	}
}
