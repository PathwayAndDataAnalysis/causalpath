package org.panda.causalpath.network;

/**
 * Enumeration of the relation types used in the causality framework.
 */
public enum RelationType
{
	UPREGULATES_EXPRESSION("upregulates-expression", 1, false, true),
	DOWNREGULATES_EXPRESSION("downregulates-expression", -1, false, true),
	PHOSPHORYLATES("phosphorylates", 1, true, false),
	DEPHOSPHORYLATES("dephosphorylates", -1, true, false);

	/**
	 * Whether the relation can explain a change in phosphorylation.
	 */
	public boolean affectsPhosphoSite;

	/**
	 * Whether the relation can explain a change in total protein.
	 */
	public boolean affectsTotalProt;

	/**
	 * Sign of the relation: positive (1) or negative (-1).
	 */
	public int sign;

	/**
	 * The name of the relation will be used in the generated SIF graph.
	 */
	public String name;

	RelationType(String name, int sign, boolean affectsPhosphoSite, boolean affectsTotalProt)
	{
		this.name = name;
		this.sign = sign;
		this.affectsPhosphoSite = affectsPhosphoSite;
		this.affectsTotalProt = affectsTotalProt;
	}
}
