package org.panda.causalpath.exceptions;

public class IllegalDirectoryException extends IllegalArgumentException {
    public IllegalDirectoryException() { }

    public IllegalDirectoryException(String s) {
        super(s);
    }

    public IllegalDirectoryException(String message, Throwable cause) {
        super(message, cause);
    }

    public IllegalDirectoryException(Throwable cause) {
        super(cause);
    }
}
