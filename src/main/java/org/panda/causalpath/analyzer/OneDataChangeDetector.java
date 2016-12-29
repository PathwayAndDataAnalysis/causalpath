package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;

/**
 * Detects a change in an experiment data.
 */
public interface OneDataChangeDetector
{
	/**
	 * Returns 1 if the change is positive, -1 if the change is negative, 0 if the change is insignificant (not
	 * detected).
	 */
	int getChangeSign(ExperimentData data);

	/**
	 * Gets the numerical value of the change.
	 */
	double getChangeValue(ExperimentData data);
}
