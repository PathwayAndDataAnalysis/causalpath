package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.ExperimentData;

/**
 * Detects a change in an experiment data.
 *
 * Created by babur on 3/24/16.
 */
public interface OneDataChangeDetector
{
	int getChangeSign(ExperimentData data);
	double getChangeValue(ExperimentData data);
}
