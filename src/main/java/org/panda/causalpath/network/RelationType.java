package org.panda.causalpath.network;

/**
 * Enumeration of the relation types used in the causality framework.
 */
public enum RelationType
{
	UPREGULATES_EXPRESSION(1, false, true, false),
	DOWNREGULATES_EXPRESSION(-1, false, true, false),
	PHOSPHORYLATES(1, true, false, false),
	DEPHOSPHORYLATES(-1, true, false, false),

	// This means source is a GEF that separates GDP from inactive GTPase protein.
	ACTIVATES_GTPASE(1, false, false, true),

	// This means source is a GAP that activates GTP hydrolysis function of the GTPase, which makes GTPase inactive.
	INHIBITS_GTPASE(-1, false, false, true);

	/**
	 * Whether the relation can explain a change in phosphorylation.
	 */
	public boolean affectsPhosphoSite;

	/**
	 * Whether the relation can explain a change in total protein.
	 */
	public boolean affectsTotalProt;

	/**
	 * Whether the relation can chnage GTPase activity.
	 */
	public boolean affectsGTPase;

	/**
	 * Sign of the relation: positive (1) or negative (-1).
	 */
	public int sign;

	RelationType(int sign, boolean affectsPhosphoSite, boolean affectsTotalProt, boolean affectsGTPase)
	{
		this.sign = sign;
		this.affectsPhosphoSite = affectsPhosphoSite;
		this.affectsTotalProt = affectsTotalProt;
		this.affectsGTPase = affectsGTPase;
	}

	public String getName()
	{
		return toString().toLowerCase().replaceAll("_", "-");
	}

	public static RelationType getType(String name)
	{
		try
		{
			return valueOf(name.toUpperCase().replaceAll("-", "_"));
		}
		catch (Exception e)
		{
			return null;
		}
	}
}
