package org.panda.causalpath.exceptions;

import org.panda.causalpath.data.ExperimentData;

public class ThresholdNotSetException extends NullPointerException
{
    private ExperimentData data;
    public ThresholdNotSetException() { }

    public ThresholdNotSetException(String s) {
        super(s);
    }

    public ThresholdNotSetException(String s, ExperimentData data) {
        super(s);
        this.data = data;
    }

    public ExperimentData getExperimentData() { return data; }
}
