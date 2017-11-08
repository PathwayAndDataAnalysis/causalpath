package org.panda.causalpath.data;

import org.panda.resource.tcga.ProteomicsFileRow;

/**
 * Represents the detected or inferred transcriptional factor activity.
 *
 * @author Ozgun Babur
 */
public class ActivityOfTFData extends ActivityData
{
	public ActivityOfTFData(String id, String symbol)
	{
		super(id, symbol);
	}

	public ActivityOfTFData(ProteomicsFileRow rppa)
	{
		super(rppa);
	}
}
