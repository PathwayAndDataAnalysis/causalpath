package org.panda.causalpath.exceptions;

public class IllegalThresholdDataFormatException extends IllegalArgumentException {
    public IllegalThresholdDataFormatException() { }

    public IllegalThresholdDataFormatException(String s) {
        super(s);
    }

    public IllegalThresholdDataFormatException(String message, Throwable cause) {
        super(message, cause);
    }

    public IllegalThresholdDataFormatException(Throwable cause) {
        super(cause);
    }
}
