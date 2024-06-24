package org.panda.causalpath.exceptions;

public class IllegalFeatureFormatException extends IllegalArgumentException {
    public IllegalFeatureFormatException() { }

    public IllegalFeatureFormatException(String s) {
        super(s);
    }

    public IllegalFeatureFormatException(String message, Throwable cause) {
        super(message, cause);
    }

    public IllegalFeatureFormatException(Throwable cause) {
        super(cause);
    }
}
