package org.panda.causalpath.exceptions;

public class FileDoesNotExistException extends IllegalArgumentException {
    public FileDoesNotExistException() { }

    public FileDoesNotExistException(String s) {
        super(s);
    }

    public FileDoesNotExistException(String message, Throwable cause) {
        super(message, cause);
    }

    public FileDoesNotExistException(Throwable cause) {
        super(cause);
    }
}
