package org.panda.causalpath.analyzer;

/**
 * Base class for change detectors that compare two group of values.
 *
 * Created by babur on 3/24/16.
 */
public abstract class ComparisonDetector extends ThresholdDetector
{
	protected boolean[] control;
	protected boolean[] test;

	public ComparisonDetector(double threshold, boolean[] control, boolean[] test)
	{
		super(threshold);
		this.control = control;
		this.test = test;
	}
}
