package org.panda.causalpath.network;

/**
 * Created by babur on 3/24/16.
 */
public enum RelationType
{
	UPREGULATES_EXPRESSION("upregulates-expression", 1, false, true),
	DOWNREGULATES_EXPRESSION("downregulates-expression", -1, false, true),
	PHOSPHORYLATES("phosphorylates", 1, true, false),
	DEPHOSPHORYLATES("dephosphorylates", -1, true, false);

	public boolean affectsPhosphoSite;
	public boolean affectsTotalProt;
	public int sign;
	public String name;

	RelationType(String name, int sign, boolean affectsPhosphoSite, boolean affectsTotalProt)
	{
		this.name = name;
		this.sign = sign;
		this.affectsPhosphoSite = affectsPhosphoSite;
		this.affectsTotalProt = affectsTotalProt;
	}
}
