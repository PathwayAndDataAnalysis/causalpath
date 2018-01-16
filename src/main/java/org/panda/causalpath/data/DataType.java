package org.panda.causalpath.data;

import java.util.*;

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

	public String getName()
	{
		return toString().toLowerCase();
	}

	public static DataType get(String name)
	{
		return valueOf(name.toUpperCase());
	}

	public static Map getValuesAsJson()
	{
		List list = new ArrayList<>();
		for (DataType type : values())
		{
			list.add(type.getName());
		}

		Map map = new LinkedHashMap<>();
		map.put("name", "DataType");
		map.put("values", list);
		return map;
	}
}
