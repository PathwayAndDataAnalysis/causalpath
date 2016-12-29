package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;

/**
 * This two-data change detector checks if the data changed in similar (1) or opposite (-1) directions. 0 means at least
 * one data is not changed. This class helps causality detection when controls and tests are compared.
 */
public class CausalityHelper implements TwoDataChangeDetector
{
	@Override
	public int getChangeSign(ExperimentData data1, ExperimentData data2)
	{
		return data1.getChangeSign() * data2.getChangeSign();
	}
}
