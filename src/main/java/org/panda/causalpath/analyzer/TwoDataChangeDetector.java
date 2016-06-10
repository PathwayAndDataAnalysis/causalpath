package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;

/**
 * Detects a change in a pair of experiment data.
 *
 * Created by babur on 3/24/16.
 */
public interface TwoDataChangeDetector
{
	/**
	 * 1: Ends change in the same direction
	 * -1: Ends change in opposite directions
	 * 0: At least one end did not change
	 *
	 * @param data1
	 * @param data2
	 * @return
	 */
	int getChangeSign(ExperimentData data1, ExperimentData data2);
}
