package org.panda.causalpath.analyzer;

/**
 * Base class for change detectors that compare two group of values.
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

	public boolean[] getControl()
	{
		return control;
	}

	public boolean[] getTest()
	{
		return test;
	}
}
