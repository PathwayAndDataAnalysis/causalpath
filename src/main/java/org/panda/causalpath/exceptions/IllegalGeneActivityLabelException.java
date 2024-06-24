package org.panda.causalpath.exceptions;

public class IllegalGeneActivityLabelException extends IllegalArgumentException {
    public IllegalGeneActivityLabelException() { }

    public IllegalGeneActivityLabelException(String s) {
        super(s);
    }

    public IllegalGeneActivityLabelException(String message, Throwable cause) {
        super(message, cause);
    }

    public IllegalGeneActivityLabelException(Throwable cause) {
        super(cause);
    }
}
