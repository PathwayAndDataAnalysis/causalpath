package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;

/**
 * Created by babur on 3/24/16.
 */
public class CausalityHelper implements TwoDataChangeDetector
{
	@Override
	public int getChangeSign(ExperimentData data1, ExperimentData data2)
	{
		return data1.getChangeSign() * data2.getChangeSign();
	}
}
