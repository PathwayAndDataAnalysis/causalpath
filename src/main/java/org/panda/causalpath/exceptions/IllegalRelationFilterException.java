package org.panda.causalpath.exceptions;

public class IllegalRelationFilterException extends IllegalArgumentException {
    public IllegalRelationFilterException() { }

    public IllegalRelationFilterException(String s) {
        super(s);
    }

    public IllegalRelationFilterException(String message, Throwable cause) {
        super(message, cause);
    }

    public IllegalRelationFilterException(Throwable cause) {
        super(cause);
    }
}
