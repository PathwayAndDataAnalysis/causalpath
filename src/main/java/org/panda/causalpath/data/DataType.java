package org.panda.causalpath.data;

/**
 * @author Ozgun Babur
 */
public enum DataType
{
	PROTEIN(true),
	PHOSPHOPROTEIN(true),
	RNA(true),
	CNA(false),
	MUTATION(false),
	METHYLATION(true),
	ACTIVITY(false);

	boolean numerical;

	DataType(boolean numerical)
	{
		this.numerical = numerical;
	}

	public boolean isNumerical()
	{
		return numerical;
	}

	public static DataType get(String name)
	{
		return valueOf(name.toUpperCase());
	}
}
